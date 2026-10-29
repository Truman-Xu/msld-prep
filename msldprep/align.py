"""Load, CGenFF-type, and 3D-align a set of ligands; find their shared core.
"""

from __future__ import annotations

import itertools
import logging
import os
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
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


def atom_pattern(atom: Chem.Atom) -> tuple[str, tuple[str, ...]]:
    """An atom's CGenFF type and the sorted CGenFF types of all its neighbors,
    hydrogens included: legacy `msld_mcs.py` step (i). Two atoms can share a
    core position only if their patterns are equal (design note (e)), so an
    atom whose surroundings change between ligands, and with them its charge,
    is left out of the core."""
    neighbor_types = sorted(neighbor.GetProp("CGENFFAtomType") for neighbor in atom.GetNeighbors())
    return atom.GetProp("CGENFFAtomType"), tuple(neighbor_types)


def assign_pattern_labels(mols: Sequence[Chem.Mol]) -> dict[tuple[str, tuple[str, ...]], int]:
    """Label every heavy atom of `mols` with an isotope unique to its
    `atom_pattern` (starting at 300), and every hydrogen with 0, so the
    shared-core search can match patterns with `CompareIsotopes` and leave the
    hydrogens out. The labels are only comparable within one call, so label
    the whole ligand set, reference included, together.

    Returns the pattern -> isotope_id mapping used, so callers can log/inspect it.
    """
    pattern_to_isotope: dict[tuple[str, tuple[str, ...]], int] = {}
    for mol in mols:
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                atom.SetIsotope(0)
                continue
            pattern = atom_pattern(atom)
            isotope_id = pattern_to_isotope.setdefault(pattern, 300 + len(pattern_to_isotope))
            atom.SetIsotope(isotope_id)
    logger.debug("Assigned %d isotope labels for atom patterns", len(pattern_to_isotope))
    return pattern_to_isotope


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

    `check_protonation` leaves the search and the pose as they are and only
    reports, in `CoreAlignment.protonation_differences`, the core atoms whose
    protonation differs from the reference's (see `ProtonationDifference`).
    Off by default: in MSLD a different protonation state is usually an
    intended change.
    """

    complete_rings_only: bool = True
    ring_matches_ring_only: bool = False
    timeout: int | None = 10
    ignore_bond_orders: bool = False
    # How many distinct ways the core may sit on a symmetric part of the analog
    # (a ring that can flip) to try; 1 uses the first match only.
    max_orientations: int = 32
    # How many stereoisomers to try when the analog leaves stereocentres unspecified.
    max_stereoisomers: int = 16
    check_protonation: bool = False


@dataclass(frozen=True)
class ProtonationDifference:
    """A core atom of the analog (`atom`) whose formal charge, or number of
    attached atoms counting hydrogens, differs from that of the reference atom
    it matches (`ref_atom`): a protonation or tautomer change. A substituent in
    place of a hydrogen leaves that number as it is, so it is not reported.
    """

    atom: int
    ref_atom: int
    charge: int
    ref_charge: int
    hydrogens: int
    ref_hydrogens: int


def _protonation_differences(mol: Chem.Mol, ref_mol: Chem.Mol, atom_map: dict[int, int]) -> tuple[ProtonationDifference, ...]:
    differences = []
    for atom_idx, ref_idx in sorted(atom_map.items()):
        atom, ref_atom = mol.GetAtomWithIdx(atom_idx), ref_mol.GetAtomWithIdx(ref_idx)
        if atom.GetFormalCharge() != ref_atom.GetFormalCharge() or atom.GetTotalDegree() != ref_atom.GetTotalDegree():
            differences.append(ProtonationDifference(
                atom=atom_idx, ref_atom=ref_idx,
                charge=atom.GetFormalCharge(), ref_charge=ref_atom.GetFormalCharge(),
                hydrogens=atom.GetTotalNumHs(includeNeighbors=True),
                ref_hydrogens=ref_atom.GetTotalNumHs(includeNeighbors=True),
            ))
    return tuple(differences)


@dataclass(frozen=True)
class StereoCenter:
    """What `align_to_reference` did about one tetrahedral stereocentre of the analog.

    `resolved_by` is "analog" when the analog's structure specifies it (kept as
    given); "reference" when it was unspecified, inside the core, and the
    reference has a stereocentre on the atom it maps to, so the configuration
    that fits the reference's geometry was taken; and "energy" otherwise (the
    reference cannot decide it), so the configuration with the lowest
    force-field energy of the aligned pose was taken: nothing the user gave
    decided it, so it deserves a look. "equivalent" marks an unspecified centre
    whose every setting gives the same molecule (nothing was decided). `opposes_reference` is True for a
    specified centre inside the core, on an atom where the reference has a
    stereocentre, whose inverted configuration fits the reference clearly
    better (the analog is probably the opposite configuration there); None
    when that does not apply or could not be checked.
    """

    atom: int
    in_core: bool
    resolved_by: str
    opposes_reference: bool | None = None


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
    stereo: tuple[StereoCenter, ...] = ()
    protonation_differences: tuple[ProtonationDifference, ...] = ()  # filled only with `McsOptions.check_protonation`


def _specified_stereo(mol: Chem.Mol) -> dict[tuple[int, int], tuple]:
    """The specified stereo elements of `mol`, keyed by kind and atom or bond, with
    their descriptor relative to their controlling atoms."""
    return {
        (int(info.type), info.centeredOn): (info.descriptor, tuple(info.controllingAtoms))
        for info in Chem.FindPotentialStereo(mol)
        if info.specified == Chem.StereoSpecified.Specified
    }


def _keeps_stereo(wanted: dict[tuple[int, int], tuple], mol: Chem.Mol) -> bool:
    """Whether the 3D structure of `mol` has every stereo element in `wanted`."""
    placed = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(placed)
    found = _specified_stereo(placed)
    return all(found.get(key) == value for key, value in wanted.items())


def _embed_onto_reference(mol: Chem.Mol, ref_mol: Chem.Mol, atom_map: dict[int, int], seed: int) -> float:
    """Embed `mol` with the atoms in `atom_map` (index in `mol` to index in
    `ref_mol`) placed on the reference's positions, optimize it with those
    atoms held there, and fit it onto them. Returns the RMSD of the mapped
    atoms. Mutates `mol`.

    Adapted from RDKit's `AllChem.ConstrainedEmbed` (BSD-3-Clause, with
    `useTethers=True`): the same steps, but the atom pairs come from the MCS
    that was already found rather than from matching a core query again,
    which would require the charges and bond orders of the reference to
    match exactly.

    The first attempt embeds from random coordinates, then puts the core atoms
    exactly on the reference's positions and keeps them fixed while the rest is
    optimized, so the core has the reference's geometry, ring puckers included.
    A fixed core can force a stereocentre of the analog into the opposite
    configuration when the analog specifies the opposite of the reference; the
    pose is then rejected, and the core is held by a stiff tether instead, and
    if that does not embed, by RDKit's default start and soft tether, which let
    the core give way a little (a higher RMSD, which is how such a conflict is
    recognized). The first pose with the analog's stereochemistry is kept, or
    the last one if none has it.
    """
    ref_conf = ref_mol.GetConformer()
    pairs = list(atom_map.items())
    coord_map = {i: ref_conf.GetAtomPosition(j) for i, j in pairs}
    wanted = _specified_stereo(mol)
    starts: dict[bool, Chem.Mol | None] = {}
    pose = None
    for random_start, force in _EMBED_ATTEMPTS:
        if random_start not in starts:
            start = Chem.Mol(mol)
            embedded = AllChem.EmbedMolecule(start, coordMap=coord_map, randomSeed=seed, useRandomCoords=random_start) >= 0
            starts[random_start] = start if embedded else None
        if starts[random_start] is None:
            continue
        pose = Chem.Mol(starts[random_start])
        rdMolAlign.AlignMol(pose, ref_mol, atomMap=pairs)
        if force is None:
            conf = pose.GetConformer()
            for i, j in pairs:
                conf.SetAtomPosition(i, ref_conf.GetAtomPosition(j))
        ff = AllChem.UFFGetMoleculeForceField(pose, confId=0)
        for i, j in pairs:
            if force is None:
                ff.AddFixedPoint(i)
            else:
                point = ref_conf.GetAtomPosition(j)
                extra = ff.AddExtraPoint(point.x, point.y, point.z, fixed=True) - 1
                ff.AddDistanceConstraint(extra, i, 0, 0, force)
        ff.Initialize()
        remaining = 4
        more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
        while more and remaining:
            more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
            remaining -= 1
        if _keeps_stereo(wanted, pose):
            break
    if pose is None:
        raise ValueError("Could not embed molecule.")
    mol.RemoveAllConformers()
    mol.AddConformer(Chem.Conformer(pose.GetConformer()), assignId=True)
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


_TETHER_FORCE = 1000.0  # kcal/mol/angstrom^2 holding the core on the reference during cleanup
_SOFT_TETHER_FORCE = 100.0  # the softer tether of the last attempt (RDKit's ConstrainedEmbed value)
# How `_embed_onto_reference` tries to place the analog, in order: whether the embedding
# starts from random coordinates, and the force holding the core (None: the core is fixed).
_EMBED_ATTEMPTS = ((True, None), (True, _TETHER_FORCE), (False, _SOFT_TETHER_FORCE))
_RMSD_TOLERANCE = 0.05  # angstrom: fits of the core this close count as equally good
_ENERGY_TIE = 1.0  # kcal/mol: stereoisomer poses this close in energy count as tied (a single pose is noisy)
_OPPOSITE_FITS_BETTER = 0.15  # angstrom: how much better the inverted centre must fit to be flagged


def _reference_stereo_atoms(ref_mol: Chem.Mol) -> set[int]:
    """Indices of the reference's tetrahedral stereocentres, read from its 3D coordinates."""
    reference = Chem.Mol(ref_mol)
    Chem.AssignStereochemistryFrom3D(reference)
    return {
        info.centeredOn
        for info in Chem.FindPotentialStereo(reference)
        if info.type == Chem.StereoType.Atom_Tetrahedral and info.specified == Chem.StereoSpecified.Specified
    }


def _uff_energy(mol: Chem.Mol) -> float:
    """The UFF energy of the conformer, or infinity when UFF lacks parameters for it
    (such a variant then loses every tie to one that has an energy)."""
    if not AllChem.UFFHasAllMoleculeParams(mol):
        return float("inf")
    return float(AllChem.UFFGetMoleculeForceField(mol, confId=0).CalcEnergy())


def _embed_best_stereoisomer(
    mol: Chem.Mol,
    ref_mol: Chem.Mol,
    atom_map: dict[int, int],
    seed: int,
    *,
    given: dict[int, bool],
    ref_stereo: set[int],
    limit: int,
) -> tuple[float, Chem.Mol, Chem.Mol, int]:
    """Embed `mol` onto the reference with its unspecified stereo settled, as
    `(rmsd, aligned copy, the variant it was embedded from, how many different
    variants of the still-unspecified stereo there were)`.

    Unspecified centres inside the core on an atom where the reference has a
    stereocentre are set the way that fits the reference (lowest core RMSD).
    Whatever is still unspecified (other centres, double bonds) is then set
    each way and the pose with the lowest energy kept; poses within
    `_ENERGY_TIE` of the lowest count as tied and the first of them wins, so
    the choice does not wander with the seed when they are alike. Raises
    ValueError if nothing embeds.
    """
    backed = [a for a, specified in given.items() if not specified and atom_map.get(a) in ref_stereo]
    base, fitted = mol, None
    if backed:
        tags = (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        best = None
        for combination in itertools.islice(itertools.product(tags, repeat=len(backed)), limit):
            variant = Chem.Mol(mol)
            for atom, tag in zip(backed, combination):
                variant.GetAtomWithIdx(atom).SetChiralTag(tag)
            trial = Chem.Mol(variant)
            try:
                trial_rmsd = _embed_onto_reference(trial, ref_mol, atom_map, seed)
            except ValueError:
                continue
            if best is None or trial_rmsd < best[0]:
                best = (trial_rmsd, trial, variant)
        if best is None:
            raise ValueError("Could not embed molecule.")
        fitted, base = best, best[2]

    options = StereoEnumerationOptions(onlyUnassigned=True, unique=True, maxIsomers=limit)
    variants = list(EnumerateStereoisomers(base, options=options))
    if len(variants) <= 1 and fitted is not None:
        return (*fitted, 1)
    posed = []
    for variant in variants or [base]:
        trial = Chem.Mol(variant)
        try:
            trial_rmsd = _embed_onto_reference(trial, ref_mol, atom_map, seed)
        except ValueError:
            continue
        posed.append((_uff_energy(trial), trial_rmsd, trial, variant))
    if not posed:
        raise ValueError("Could not embed molecule.")
    lowest = min(entry[0] for entry in posed)
    _, trial_rmsd, trial, variant = next(entry for entry in posed if entry[0] <= lowest + _ENERGY_TIE)
    return trial_rmsd, trial, variant, max(1, len(variants))


def _stereo_report(
    chosen: Chem.Mol, ref_mol: Chem.Mol, atom_map: dict[int, int], ref_stereo: set[int],
    atom: int, specified: bool, rmsd: float, seed: int, free_variants: int,
) -> StereoCenter:
    in_core = atom in atom_map
    backed = in_core and atom_map[atom] in ref_stereo
    if not specified:
        if backed:
            resolved_by = "reference"
        else:
            # with one variant only, every way of setting it is the same molecule: nothing was decided
            resolved_by = "energy" if free_variants > 1 else "equivalent"
        return StereoCenter(atom=atom, in_core=in_core, resolved_by=resolved_by)
    opposes = None
    if backed:
        inverted = Chem.Mol(chosen)
        inverted.GetAtomWithIdx(atom).InvertChirality()
        try:
            opposite_rmsd = _embed_onto_reference(inverted, ref_mol, atom_map, seed)
        except ValueError:
            opposes = False
        else:
            opposes = opposite_rmsd < rmsd - _OPPOSITE_FITS_BETTER
    return StereoCenter(atom=atom, in_core=in_core, resolved_by="analog", opposes_reference=opposes)


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
    # Stereo as the analog gives it, before anything is chosen for it.
    given = {
        info.centeredOn: info.specified == Chem.StereoSpecified.Specified
        for info in Chem.FindPotentialStereo(mol)
        if info.type == Chem.StereoType.Atom_Tetrahedral
    }
    ref_stereo = _reference_stereo_atoms(ref_mol)

    # For each orientation, the pose with its unspecified stereo settled (see
    # `_embed_best_stereoisomer`).
    per_orientation: list[tuple[int, float, float, Chem.Mol, Chem.Mol, int]] = []
    for index, candidate in enumerate(candidates):
        atom_map = dict(zip(candidate, ref_match))
        try:
            trial_rmsd, trial, variant, free_variants = _embed_best_stereoisomer(
                mol, ref_mol, atom_map, seed, given=given, ref_stereo=ref_stereo, limit=max(1, mcs.max_stereoisomers)
            )
        except ValueError:
            continue
        score = _substituent_overlap(trial, set(candidate), ref_mol, ref_core) if len(candidates) > 1 else 0.0
        per_orientation.append((index, trial_rmsd, score, trial, variant, free_variants))
    if not per_orientation:
        raise ValueError("Could not embed molecule.")

    # Among orientations that fit the core equally well, the one whose substituents lie
    # nearest the reference's; a later one must be clearly better than an earlier one.
    fit = min(entry[1] for entry in per_orientation)
    chosen, rmsd, best_score, best_trial, chosen_variant, free_variants = next(e for e in per_orientation if e[1] <= fit + _RMSD_TOLERANCE)
    for entry in per_orientation:
        if entry[1] <= fit + _RMSD_TOLERANCE and entry[2] > best_score + 1e-3:
            chosen, rmsd, best_score, best_trial, chosen_variant, free_variants = entry
    atom_map = dict(zip(candidates[chosen], ref_match))
    mol.RemoveAllConformers()
    mol.AddConformer(Chem.Conformer(best_trial.GetConformer()), assignId=True)
    mol.SetProp("EmbedRMS", str(rmsd))
    Chem.AssignStereochemistryFrom3D(mol)  # the structure now has the stereochemistry it was placed with

    stereo = tuple(
        _stereo_report(chosen_variant, ref_mol, atom_map, ref_stereo, atom, specified, rmsd, seed, free_variants)
        for atom, specified in sorted(given.items())
    )
    if mcs_result.canceled:
        logger.warning("The MCS search hit its %s s limit; using the best structure found.", mcs.timeout)
    protonation = _protonation_differences(mol, ref_mol, atom_map) if mcs.check_protonation else ()
    if protonation:
        logger.warning(
            "Protonation differs from the reference at core atom(s) %s.", ", ".join(str(d.atom) for d in protonation)
        )
    return CoreAlignment(
        core=core, atom_map=atom_map, rmsd=rmsd, seed=seed, mcs_timed_out=bool(mcs_result.canceled),
        orientations=len(candidates), orientation=chosen, stereo=stereo, protonation_differences=protonation,
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
    """Find the maximum common substructure across `ligands` whose atoms have
    the same `atom_pattern` in every ligand (design note (e)). Heavy atoms only:
    a core atom's hydrogens enter through its pattern, and `decompose` adds
    them to the core.

    A ring may be cut, as legacy cuts it, but ring atoms match ring atoms
    only. Any bond matches any bond: the CGenFF types already carry
    aromaticity. Requires every ligand to carry the labels of
    `assign_pattern_labels`, run over the exact same ligand set, including
    the reference.
    """
    mols = [ligand.mol for ligand in ligands]
    for ligand in ligands:
        if not any(atom.GetIsotope() for atom in ligand.mol.GetAtoms()):
            msg = (
                f"Ligand {ligand.name} has no isotope-labeled atoms; run "
                "assign_pattern_labels(...) over the full ligand set (including "
                "the reference) before calling find_shared_core."
            )
            logger.error(msg)
            raise AlignmentError(msg)
    # Heavy-atom copies: RemoveAllHs drops every hydrogen (RemoveHs would keep
    # isotope-labeled ones, #15) and keeps the heavy atoms' labels.
    result = FindMCS(
        [Chem.RemoveAllHs(mol) for mol in mols],
        atomCompare=AtomCompare.CompareIsotopes,
        bondCompare=BondCompare.CompareAny,
        ringMatchesRingOnly=True,
        completeRingsOnly=False,
    )
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
    isotope_map = assign_pattern_labels([lig.mol for lig in ligands])
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
