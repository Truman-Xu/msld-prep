"""Load, CGenFF-type, and 3D-align a set of ligands; find their shared core.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from typing import Sequence

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdFMCS import AtomCompare, FindMCS

from .datamodel import LigandRecord
from .errors import AlignmentError
from .rtf_records import assign_atom_types, parse_rtf_file

logger = logging.getLogger(__name__)

# RDKit raises these (not a documented, closed set) for bad/unembeddable
# substructure geometry; anything else is unexpected and re-raised.
_EMBEDDING_EXCEPTIONS = (ValueError, RuntimeError)

# Sdf has no native per-atom name field; each ligand's atom names are
# stashed in this same-order sd data field instead (see `write_ligands_sdf`).
_SDF_ATOM_NAMES_PROP = "AtomNames"


def assign_atom_types_from_rtf(mol: Chem.Mol, rtf_path: str) -> None:
    """Stamp CGENFFAtomType/CGENFFPartialCharge onto every atom of `mol`.

    `mol` must have been loaded from a mol2 source (so atoms carry
    `_TriposAtomName`). Matches atoms to the rtf's per-atom records by name.
    Thin wrapper around `rtf_records.parse_rtf_file`/`assign_atom_types`;
    kept here since existing callers import it from `align`.
    """
    rtf_fields = parse_rtf_file(rtf_path)
    assign_atom_types(mol, rtf_fields["atom_types"], rtf_fields["atom_charges"], source_label=rtf_path)


def get_all_atom_types(mols: Sequence[Chem.Mol]) -> dict[str, list[Chem.Atom]]:
    """Group every atom across `mols` by its CGENFFAtomType."""
    all_types: dict[str, list[Chem.Atom]] = {}
    for mol in mols:
        for atom in mol.GetAtoms():
            all_types.setdefault(atom.GetProp("CGENFFAtomType"), []).append(atom)
    return all_types


def assign_isotope_labels(atom_types: dict[str, list[Chem.Atom]]) -> dict[str, int]:
    """Assign each atom an isotope label unique to its CGenFF type (starting at 300).

    Returns the type -> isotope_id mapping used, so callers can log/inspect it.
    """
    type_to_isotope: dict[str, int] = {}
    for isotope_id, (cgenff_type, atoms) in enumerate(atom_types.items(), start=300):
        type_to_isotope[cgenff_type] = isotope_id
        for atom in atoms:
            atom.SetIsotope(isotope_id)
    logger.debug("Assigned %d isotope labels for CGenFF types", len(type_to_isotope))
    return type_to_isotope


def load_ligand(
    mol2_path: str, rtf_path: str, prm_path: str, name: str | None = None, *, is_reference: bool = False
) -> LigandRecord:
    """Load a ligand from mol2 + rtf + prm, typing every atom with its
    CGenFF type/charge. Thin wrapper around `LigandRecord.from_files`,
    kept here since existing callers import it from `align`.
    """
    return LigandRecord.from_files(mol2_path, rtf_path, prm_path, name, is_reference=is_reference)


def write_ligands_sdf(ligands: Sequence[LigandRecord], sdf_path: str) -> None:
    """Write every ligand's current conformer into one combined sdf file
    (one record per ligand, keyed by `_Name` = `ligand.name`), instead of
    one mol2 file per ligand.

    This exists to dodge a real failure mode of the mol2 round trip: our
    mol2 writer (`sybyl_typer.write_mol2_file`) labels every aromatic bond
    "ar" (Sybyl's aromatic bond type), so reloading via
    `Chem.MolFromMol2File` requires RDKit to re-guess a Kekule structure
    from connectivity + that flag alone -- a guess that can fail (a
    kekulization error) for some ring systems even though the original Mol
    was perfectly valid. `Chem.SDWriter` instead writes each bond's already
    -known explicit Kekule order directly, so nothing needs to be re-guessed
    on reload.

    Every atom must already carry a `_TriposAtomName` prop (set when the
    ligand was first loaded/named). Sdf has no native per-atom name field,
    so each record also gets a same-order `AtomNames` sd data field;
    `load_ligands_sdf` restores `_TriposAtomName` from it on reload.
    """
    writer = Chem.SDWriter(sdf_path)
    try:
        for ligand in ligands:
            mol = ligand.mol
            mol.SetProp("_Name", ligand.name)
            atom_names = [atom.GetProp("_TriposAtomName") for atom in mol.GetAtoms()]
            mol.SetProp(_SDF_ATOM_NAMES_PROP, " ".join(atom_names))
            writer.write(mol)
    finally:
        writer.close()
    logger.info("Wrote %d ligand(s) to %s", len(ligands), sdf_path)


def load_ligands_sdf(sdf_path: str) -> dict[str, Chem.Mol]:
    """Read `write_ligands_sdf`'s output back into `{ligand name: Mol}`,
    restoring `_TriposAtomName` on every atom from the `AtomNames` field it
    wrote (sdf preserves atom order exactly as written, so a straight zip
    against that field is correct).
    """
    mols: dict[str, Chem.Mol] = {}
    supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)
    for mol in supplier:
        if mol is None:
            msg = f"Failed to parse a molecule from sdf file `{sdf_path}`."
            logger.error(msg)
            raise AlignmentError(msg)
        atom_names = mol.GetProp(_SDF_ATOM_NAMES_PROP).split()
        for atom, atom_name in zip(mol.GetAtoms(), atom_names):
            atom.SetProp("_TriposAtomName", atom_name)
        mols[mol.GetProp("_Name")] = mol
    return mols


def load_ligand_set_from_sdf(
    sdf_path: str, aligned_dir: str, reference_name: str, analog_names: Sequence[str]
) -> tuple[LigandRecord, list[LigandRecord]]:
    """Reload a full ligand set (reference + analogs) from a combined
    aligned sdf (`write_ligands_sdf`'s output) plus each ligand's
    `.rtf`/`.prm` file in `aligned_dir` -- the step-2 (`build`) counterpart
    to step 1's `pipeline.define_common_core`. See `write_ligands_sdf` for
    why sdf, not mol2.
    """
    mols_by_name = load_ligands_sdf(sdf_path)

    def _load(name: str, *, is_reference: bool) -> LigandRecord:
        if name not in mols_by_name:
            msg = f"Ligand `{name}` not found in sdf file `{sdf_path}`."
            logger.error(msg)
            raise AlignmentError(msg)
        return LigandRecord.from_mol(
            mols_by_name[name],
            os.path.join(aligned_dir, f"{name}.rtf"),
            os.path.join(aligned_dir, f"{name}.prm"),
            name=name,
            is_reference=is_reference,
            mol2_path=sdf_path,
        )

    reference = _load(reference_name, is_reference=True)
    analogs = [_load(name, is_reference=False) for name in analog_names]
    return reference, analogs


def _bond_path_for_atoms(atom_ids: Sequence[int], mol: Chem.Mol) -> list[int]:
    atom_set = set(atom_ids)
    return [
        bond.GetIdx()
        for bond in mol.GetBonds()
        if bond.GetBeginAtomIdx() in atom_set and bond.GetEndAtomIdx() in atom_set
    ]


def align_common_core(mol: Chem.Mol, ref_mol: Chem.Mol, *, complete_rings_only: bool = True) -> Chem.Mol:
    """3D-align `mol` onto `ref_mol` by embedding it tethered to their shared core.

    Mutates `mol` in place (adds/replaces its conformer). Returns the core
    query submol used for the tether, for callers that want to log/inspect it.
    """
    logger.debug(f'Searching for maximum common structure')
    # FindMCS on explicit-H mols is combinatorially explosive (many
    # equivalent H atoms blow up the clique search); the tether/core we
    # need is heavy-atom-only anyway, so search without Hs and match the
    # result back onto the original (still-protonated) ref_mol below.
    mcs = FindMCS([Chem.RemoveHs(mol), Chem.RemoveHs(ref_mol)], completeRingsOnly=complete_rings_only)
    if mcs.queryMol is None:
        # No shared substructure at all (e.g. two genuinely unrelated
        # molecules): FindMCS returns a None queryMol rather than an empty
        # Mol. Raise instead of letting GetSubstructMatch below fail with an
        # opaque Boost.Python.ArgumentError -- align_ligand_set's caller
        # only catches ValueError/RuntimeError as an expected alignment
        # failure to collect and report, not crash on.
        msg = "FindMCS found no shared substructure between the two molecules."
        logger.error(msg)
        raise ValueError(msg)
    substruct_map = ref_mol.GetSubstructMatch(mcs.queryMol)
    bond_path = _bond_path_for_atoms(substruct_map, ref_mol)
    core = Chem.PathToSubmol(ref_mol, bond_path, useQuery=True)
    logger.debug(f'Embedding analog structure with MCS constraint')
    AllChem.ConstrainedEmbed(mol, core, useTethers=True)
    return core


@dataclass
class AlignmentResult:
    aligned: list[LigandRecord] = field(default_factory=list)
    failed: list[tuple[LigandRecord, Exception]] = field(default_factory=list)


def align_ligand_set(ligands: Sequence[LigandRecord], reference: LigandRecord) -> AlignmentResult:
    """Align every ligand in `ligands` onto `reference`, collecting failures instead of raising.

    Failed ligands are logged individually plus in a batch summary, and are
    never included in `.aligned` -- callers (decomposition, group core search)
    should only ever consume `.aligned`.
    """
    result = AlignmentResult()
    for ligand in ligands:
        try:
            align_common_core(ligand.mol, reference.mol)
        except _EMBEDDING_EXCEPTIONS as exc:
            logger.warning("Alignment failed for ligand %s: %r", ligand.name, exc)
            result.failed.append((ligand, exc))
        else:
            result.aligned.append(ligand)

    if result.failed:
        import warnings

        names = [lig.name for lig, _ in result.failed]
        warnings.warn(
            f"{len(result.failed)}/{len(ligands)} ligand(s) excluded from the group core "
            f"due to alignment failure: {names}",
            stacklevel=2,
        )
    return result


def find_shared_core(ligands: Sequence[LigandRecord]) -> Chem.Mol:
    """Find the CGenFF-type-locked maximum common substructure across `ligands`.

    Requires every ligand to already carry isotope labels (via
    `assign_isotope_labels(get_all_atom_types(...))`, run over the exact
    same ligand set, including the reference).
    """
    mols = [ligand.mol for ligand in ligands]
    for ligand in ligands:
        if not any(atom.GetIsotope() for atom in ligand.mol.GetAtoms()):
            msg = (
                f"Ligand {ligand.name} has no isotope-labeled atoms; run "
                "assign_isotope_labels(get_all_atom_types(...)) over the full "
                "ligand set (including the reference) before calling find_shared_core."
            )
            logger.error(msg)
            raise AlignmentError(msg)
    # See align_common_core: explicit Hs make FindMCS combinatorially
    # explosive, so search heavy-atom-only copies; Chem.RemoveHs preserves
    # each remaining atom's isotope label, so CompareIsotopes still works.
    result = FindMCS([Chem.RemoveHs(mol) for mol in mols], atomCompare=AtomCompare.CompareIsotopes, completeRingsOnly=True)
    logger.info(
        "Found shared core across %d ligands: %d atoms, %d bonds",
        len(ligands), result.numAtoms, result.numBonds,
    )
    return result.queryMol


def type_ligands_and_find_core(ligands: Sequence[LigandRecord]) -> Chem.Mol:
    """Isotope-label + find the shared core across `ligands`, assuming every
    one of them is ALREADY 3D-aligned (e.g. each was already run through
    `align_common_core` individually, as `analog_prep.build_analog_from_smiles`
    does). Does not embed/re-align anything.

    Use this instead of `prepare_ligand_set` when the caller has already
    aligned every ligand itself -- `prepare_ligand_set` always re-embeds via
    `align_ligand_set`, which would otherwise redundantly re-run
    `ConstrainedEmbed` a second time. That's wasted work, and worse, RDKit's
    `ConstrainedEmbed` isn't perfectly deterministic: a ligand that aligned
    fine the first time can occasionally fail on a second, independent
    attempt, silently dropping it (or, if it drops all of them, crashing
    `find_shared_core` outright since RDKit's `FindMCS` requires at least 2
    molecules).
    """
    isotope_map = assign_isotope_labels(get_all_atom_types([lig.mol for lig in ligands]))
    logger.debug("Isotope map for shared core search: %s", isotope_map)
    return find_shared_core(ligands)


@dataclass
class PreparedLigandSet:
    reference: LigandRecord
    aligned: list[LigandRecord]
    failed: list[tuple[LigandRecord, Exception]]
    core_query_mol: Chem.Mol


def prepare_ligand_set(reference: LigandRecord, analogs: Sequence[LigandRecord]) -> PreparedLigandSet:
    """Full alignment pipeline: align every analog onto `reference`, type+isotope-label
    the reference and every successfully-aligned analog together, and find their shared core.
    """
    alignment = align_ligand_set(analogs, reference)
    all_ligands = [reference, *alignment.aligned]
    core_query_mol = type_ligands_and_find_core(all_ligands)
    return PreparedLigandSet(
        reference=reference,
        aligned=alignment.aligned,
        failed=alignment.failed,
        core_query_mol=core_query_mol,
    )
