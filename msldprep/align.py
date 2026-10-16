"""Load, CGenFF-type, and 3D-align a set of ligands; find their shared core.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Chem.rdFMCS import AtomCompare, BondCompare, FindMCS

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


def _check_3d(mol: Chem.Mol, source: str) -> None:
    if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D():
        msg = f"{source} has no 3D coordinates; a reference to align onto needs them."
        logger.error(msg)
        raise AlignmentError(msg)


def mol_from_mol2_text(text: str, *, require_3d: bool = True) -> Chem.Mol:
    """Parse one molecule from mol2 text, keeping explicit Hs and sanitizing,
    so a mol2 upload needs no temp file. Atoms carry `_TriposAtomName`.
    Raises AlignmentError if the text is not a single readable molecule, or
    (unless `require_3d` is False) has no 3D coordinates.
    """
    try:
        mol = Chem.MolFromMol2Block(text, sanitize=True, removeHs=False)
    except Exception as exc:  # RDKit's parse errors are not a documented, closed set
        msg = f"Could not parse mol2 text: {exc!r}"
        logger.error(msg)
        raise AlignmentError(msg) from exc
    if mol is None:
        msg = "Could not parse mol2 text."
        logger.error(msg)
        raise AlignmentError(msg)
    if require_3d:
        _check_3d(mol, "The mol2 text")
    return mol


def mol_from_sdf_text(text: str, *, require_3d: bool = True) -> Chem.Mol:
    """Parse exactly one molecule from sdf (or mol block) text, keeping
    explicit Hs and sanitizing, so an sdf upload needs no temp file.
    Raises AlignmentError if the text holds no molecule, more than one, an
    unreadable one, or (unless `require_3d` is False) one without 3D
    coordinates.
    """
    supplier = Chem.SDMolSupplier()
    supplier.SetData(text, sanitize=True, removeHs=False)
    count = len(supplier)
    if count != 1:
        msg = f"Expected one molecule in the sdf text, found {count}."
        logger.error(msg)
        raise AlignmentError(msg)
    mol = supplier[0]
    if mol is None:
        msg = "Could not parse the molecule in the sdf text."
        logger.error(msg)
        raise AlignmentError(msg)
    if require_3d:
        _check_3d(mol, "The sdf text")
    return mol


def _bond_path_for_atoms(atom_ids: Sequence[int], mol: Chem.Mol) -> list[int]:
    atom_set = set(atom_ids)
    return [
        bond.GetIdx()
        for bond in mol.GetBonds()
        if bond.GetBeginAtomIdx() in atom_set and bond.GetEndAtomIdx() in atom_set
    ]


# The seed RDKit's `ConstrainedEmbed` uses when none is given. Passed
# explicitly so the value in use is visible to callers and reproducible.
DEFAULT_EMBED_SEED = 2342


@dataclass(frozen=True)
class McsOptions:
    """The `FindMCS` settings `align_to_reference` exposes. `timeout` is in
    seconds and defaults to 10: on some pairs the search cannot prove that its
    result is the largest and, with RDKit's own limit (None here), would run
    for an hour. After the limit it returns the best structure found so far,
    flagged in `CoreAlignment.mcs_timed_out`.

    Atoms are matched by element only, so charges, hydrogen counts and
    isotope labels never matter. Bond orders do, unless `ignore_bond_orders`
    is set: then any bond matches any bond, which lets a molecule match its
    tautomer or protonation variant even where RDKit perceives the ring
    differently (aromatic in one, alternating single and double bonds in the
    other).
    """

    complete_rings_only: bool = True
    ring_matches_ring_only: bool = False
    timeout: int | None = 10
    ignore_bond_orders: bool = False
    # How many distinct ways the core may sit on a symmetric part of the analog
    # (a ring that can flip) to try; 1 uses the first match only.
    max_orientations: int = 32


@dataclass
class CoreAlignment:
    """What `align_to_reference` found and did.

    `atom_map` maps each core atom's index in the aligned `mol` to the index
    of the matching atom in the reference (heavy atoms only, since the MCS
    search runs without Hs). `rmsd` is the RMSD (angstrom) between those
    atom pairs after the final alignment onto the reference (RDKit's
    `EmbedRMS` measure, also set as that property on the aligned molecule).
    """

    core: Chem.Mol
    atom_map: dict[int, int]
    rmsd: float
    seed: int
    mcs_timed_out: bool = False  # the core is the best found within the time limit, not proven largest
    orientations: int = 1  # distinct ways the core fits the analog, up to `McsOptions.max_orientations`
    orientation: int = 0  # index of the one used


def _embed_onto_reference(mol: Chem.Mol, ref_mol: Chem.Mol, atom_map: dict[int, int], seed: int) -> float:
    """Embed `mol` with the atoms in `atom_map` (index in `mol` to index in
    `ref_mol`) placed on the reference's positions, optimize it with those
    atoms tethered there, and fit it onto them. Returns the RMSD of the
    mapped atoms. Mutates `mol`.

    Adapted from RDKit's `AllChem.ConstrainedEmbed` (BSD-3-Clause, with
    `useTethers=True`): the same steps, but the atom pairs come from the MCS
    that was already found rather than from matching a core query again,
    which would require the charges and bond orders of the reference to
    match exactly.
    """
    ref_conf = ref_mol.GetConformer()
    pairs = list(atom_map.items())
    coord_map = {i: ref_conf.GetAtomPosition(j) for i, j in pairs}
    if AllChem.EmbedMolecule(mol, coordMap=coord_map, randomSeed=seed) < 0:
        raise ValueError("Could not embed molecule.")
    rdMolAlign.AlignMol(mol, ref_mol, atomMap=pairs)
    ff = AllChem.UFFGetMoleculeForceField(mol, confId=0)
    for i, j in pairs:
        point = ref_conf.GetAtomPosition(j)
        extra = ff.AddExtraPoint(point.x, point.y, point.z, fixed=True) - 1
        ff.AddDistanceConstraint(extra, i, 0, 0, 100.0)
    ff.Initialize()
    remaining = 4
    more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
    while more and remaining:
        more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
        remaining -= 1
    rmsd = rdMolAlign.AlignMol(mol, ref_mol, atomMap=pairs)
    mol.SetProp("EmbedRMS", str(rmsd))
    return rmsd


_MAX_MATCHES = 1000  # cap on the substructure matches enumerated for one analog
_OVERLAP_SIGMA = 1.0  # angstrom: width of the overlap score between substituents


def _distinct_orientations(mol: Chem.Mol, query: Chem.Mol, ref_match: Sequence[int], limit: int) -> list[tuple[int, ...]]:
    """The ways `query` fits `mol`, one per distinct placement, the first match
    RDKit finds first. Two matches that only exchange atoms `mol` treats alike
    (the two ortho carbons of an unsubstituted phenyl, say) give the same
    structure and count once; matches that differ where `mol` is not
    symmetric (a ring with a substituent) are kept apart.
    """
    # The query's ring-bond tests need ring information, which a molecule built by
    # RenumberAtoms or edited by hand may lack; computing it again changes nothing else.
    Chem.GetSymmSSSR(mol)
    first = mol.GetSubstructMatch(query)
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    found: list[tuple[int, ...]] = []
    seen: set[frozenset[tuple[int, int]]] = set()
    for match in (first, *mol.GetSubstructMatches(query, uniquify=False, maxMatches=_MAX_MATCHES)):
        signature = frozenset((ranks[a], r) for a, r in zip(match, ref_match))
        if signature in seen:
            continue
        seen.add(signature)
        found.append(tuple(match))
        if len(found) >= limit:
            break
    return found


def _substituent_overlap(mol: Chem.Mol, core_atoms: set[int], ref_mol: Chem.Mol, ref_core_atoms: set[int]) -> float:
    """How well the heavy atoms of `mol` outside the core sit on the heavy atoms of
    the reference outside the core (a sum of Gaussians of their distances; higher
    is better). Zero when either side has none."""
    free = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1 and a.GetIdx() not in core_atoms]
    ref_free = [a.GetIdx() for a in ref_mol.GetAtoms() if a.GetAtomicNum() > 1 and a.GetIdx() not in ref_core_atoms]
    if not free or not ref_free:
        return 0.0
    xyz = mol.GetConformer().GetPositions()[free]
    ref_xyz = ref_mol.GetConformer().GetPositions()[ref_free]
    squared = ((xyz[:, None, :] - ref_xyz[None, :, :]) ** 2).sum(axis=-1)
    return float(np.exp(-squared / (2 * _OVERLAP_SIGMA**2)).sum())


def align_to_reference(
    mol: Chem.Mol, ref_mol: Chem.Mol, *, mcs: McsOptions = McsOptions(), seed: int = DEFAULT_EMBED_SEED
) -> CoreAlignment:
    """3D-align `mol` onto `ref_mol` by embedding it tethered to their shared core.

    The core is their maximum common substructure by element (see
    `McsOptions`), so differing charges or hydrogen counts do not matter.
    Mutates `mol` in place (adds/replaces its conformer). The same inputs
    and `seed` give the same conformer.

    Where the core fits the analog in more than one distinct way (a ring that
    can flip), each way is embedded and the one whose atoms outside the core
    lie closest to the reference's atoms outside the core is kept, so a
    decorated ring turns its substituents toward the reference's. Where that
    makes no difference (the reference has nothing there), the first match
    RDKit returns is used. `orientations` and `orientation` in the result say
    how many ways there were and which one was kept.
    """
    logger.debug(f'Searching for maximum common structure')
    # FindMCS on explicit-H mols is combinatorially explosive (many
    # equivalent H atoms blow up the clique search); the tether/core we
    # need is heavy-atom-only anyway, so search without Hs and match the
    # result back onto the original (still-protonated) molecules below.
    mcs_kwargs = {
        "completeRingsOnly": mcs.complete_rings_only,
        "ringMatchesRingOnly": mcs.ring_matches_ring_only,
    }
    if mcs.ignore_bond_orders:
        mcs_kwargs["bondCompare"] = BondCompare.CompareAny
    if mcs.timeout is not None:
        mcs_kwargs["timeout"] = mcs.timeout
    mcs_result = FindMCS([Chem.RemoveHs(mol), Chem.RemoveHs(ref_mol)], **mcs_kwargs)
    if mcs_result.queryMol is None:
        # No shared substructure at all (e.g. two genuinely unrelated
        # molecules): FindMCS returns a None queryMol rather than an empty
        # Mol. Raise instead of letting GetSubstructMatch below fail with an
        # opaque Boost.Python.ArgumentError -- align_ligand_set's caller
        # only catches ValueError/RuntimeError as an expected alignment
        # failure to collect and report, not crash on.
        msg = "FindMCS found no shared substructure between the two molecules."
        logger.error(msg)
        raise ValueError(msg)
    # The MCS query has heavy atoms only, so it matches the molecules with
    # their hydrogens too, and the indices it returns are indices in them.
    ref_match = ref_mol.GetSubstructMatch(mcs_result.queryMol)
    candidates = _distinct_orientations(mol, mcs_result.queryMol, ref_match, max(1, mcs.max_orientations)) if ref_match else []
    if not ref_match or not candidates or not candidates[0]:
        msg = "The common substructure found by FindMCS could not be matched back onto the molecules."
        logger.error(msg)
        raise ValueError(msg)
    core = Chem.PathToSubmol(ref_mol, _bond_path_for_atoms(ref_match, ref_mol))
    logger.debug(f'Embedding analog structure with MCS constraint')
    ref_core = set(ref_match)
    chosen = 0
    if len(candidates) == 1:
        atom_map = dict(zip(candidates[0], ref_match))
        rmsd = _embed_onto_reference(mol, ref_mol, atom_map, seed)
    else:
        best: tuple[float, Chem.Mol, float] | None = None
        for index, candidate in enumerate(candidates):
            trial = Chem.Mol(mol)
            try:
                trial_rmsd = _embed_onto_reference(trial, ref_mol, dict(zip(candidate, ref_match)), seed)
            except ValueError:
                continue
            score = _substituent_overlap(trial, set(candidate), ref_mol, ref_core)
            # a later orientation must be clearly better to replace the first one
            if best is None or score > best[0] + 1e-3:
                best, chosen = (score, trial, trial_rmsd), index
        if best is None:
            raise ValueError("Could not embed molecule.")
        atom_map = dict(zip(candidates[chosen], ref_match))
        rmsd = best[2]
        mol.RemoveAllConformers()
        mol.AddConformer(Chem.Conformer(best[1].GetConformer()), assignId=True)
        mol.SetProp("EmbedRMS", str(rmsd))
    if mcs_result.canceled:
        logger.warning("The MCS search hit its %s s limit; using the best structure found.", mcs.timeout)
    return CoreAlignment(
        core=core, atom_map=atom_map, rmsd=rmsd, seed=seed, mcs_timed_out=bool(mcs_result.canceled),
        orientations=len(candidates), orientation=chosen,
    )


def align_common_core(mol: Chem.Mol, ref_mol: Chem.Mol, *, complete_rings_only: bool = True) -> Chem.Mol:
    """3D-align `mol` onto `ref_mol` (see `align_to_reference`) and return the
    core query submol used for the tether, for callers that want to
    log/inspect it. Kept for existing callers; new code that needs the atom
    map, RMSD, seed or MCS settings should call `align_to_reference`.
    """
    return align_to_reference(mol, ref_mol, mcs=McsOptions(complete_rings_only=complete_rings_only)).core


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
