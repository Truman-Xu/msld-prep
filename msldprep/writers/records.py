"""Pure data-shaping for CHARMM/PDB writers: no file I/O here.

This is where the decision logic that legacy interleaved into its
`fp.write()` calls lives -- e.g. the atom-name translation table (CHARMM
names are capped at 4 characters), and which BOND/IMPR lines are
core-internal, fragment-internal, cross-core-fragment, or genuinely
cross-site (and therefore need to be written out as a commented-out
diagnostic rather than a live record). `writers/pdb.py` and `writers/rtf.py`
turn these records into text; they contain no decisions of their own.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Sequence

from rdkit import Chem

from ..datamodel import CoreDecomposition, CRNResult, Fragment, LigandRecord
from ..naming import allocate_unique_atom_name, new_counters
from ..rtf_records import LonePairRecord

logger = logging.getLogger(__name__)

SEGID = "LIG"


# ---------------------------------------------------------------------------
# Atom name translation (CHARMM names are capped at 4 characters)
# ---------------------------------------------------------------------------


@dataclass
class AtomNameTranslation:
    core: dict[str, str] = field(default_factory=dict)
    fragments: dict[tuple[int, str], dict[str, str]] = field(default_factory=dict)


def compute_atom_name_translation(
    decomposition: CoreDecomposition, ligands: Sequence[LigandRecord]
) -> AtomNameTranslation:
    """Build a globally-unique <=4-char CHARMM name for every core and fragment atom."""
    ligand_by_name = {lig.name: lig for lig in ligands}
    reference = ligand_by_name[decomposition.reference_name]
    ref_name_to_atom = {
        atom.GetProp("_TriposAtomName"): atom
        for atom in reference.mol.GetAtoms()
        if atom.HasProp("_TriposAtomName")
    }

    counters = new_counters()
    translation = AtomNameTranslation()

    for ref_name in decomposition.core_atom_order:
        symbol = ref_name_to_atom[ref_name].GetSymbol()
        translation.core[ref_name] = allocate_unique_atom_name(symbol, counters)
    for lp_by_ligand in decomposition.core_lp_atom_names.values():
        ref_lp_name = lp_by_ligand[decomposition.reference_name]
        translation.core[ref_lp_name] = allocate_unique_atom_name("LP", counters)

    for site in decomposition.sites:
        for ligand_name, fragment in site.fragments.items():
            mol = ligand_by_name[ligand_name].mol
            name_to_atom = {
                atom.GetProp("_TriposAtomName"): atom
                for atom in mol.GetAtoms()
                if atom.HasProp("_TriposAtomName")
            }
            frag_translation: dict[str, str] = {}
            for name in fragment.atom_names:
                frag_translation[name] = allocate_unique_atom_name(name_to_atom[name].GetSymbol(), counters)
            for lp_name in fragment.lp_atom_names:
                frag_translation[lp_name] = allocate_unique_atom_name("LP", counters)
            translation.fragments[(site.index, ligand_name)] = frag_translation

    return translation


def core_translation_for_ligand(decomposition: CoreDecomposition, translation: AtomNameTranslation, ligand_name: str) -> dict[str, str]:
    """This ligand's own local core atom names -> the shared core.rtf translated name."""
    return {
        decomposition.core_atom_names[ligand_name][pos]: translation.core[ref_name]
        for pos, ref_name in enumerate(decomposition.core_atom_order)
    }


def resolve_translated_name(
    raw_name: str,
    *,
    local_core_translation: dict[str, str],
    fragment_translation: dict[str, str],
) -> str | None:
    """Translate a raw atom name that could be either a (this ligand's) core atom or a fragment atom."""
    if raw_name in local_core_translation:
        return local_core_translation[raw_name]
    return fragment_translation.get(raw_name)


# ---------------------------------------------------------------------------
# PDB atom records
# ---------------------------------------------------------------------------


@dataclass
class PdbAtomRecord:
    serial: int
    atom_name: str
    resname: str
    x: float
    y: float
    z: float
    occ: float = 1.0
    b: float = 0.0
    segid: str = SEGID


def _heavy_then_bonded_h_order(mol: Chem.Mol, names: Sequence[str]) -> list[str]:
    """Order `names` as [heavy atom, its bonded H's in `names`, next heavy atom, ...].

    Lone pairs (name prefix "LP") are never included -- they have no fixed
    3D coordinate of their own (ported from msld_crn.py:601-619/717-735).
    """
    name_to_atom = {
        atom.GetProp("_TriposAtomName"): atom
        for atom in mol.GetAtoms()
        if atom.HasProp("_TriposAtomName")
    }
    name_set = set(names)
    ordered: list[str] = []
    added: set[str] = set()
    for name in names:
        if name in added or name.startswith("LP"):
            continue
        atom = name_to_atom[name]
        if atom.GetSymbol() == "H":
            continue
        ordered.append(name)
        added.add(name)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() != "H" or not neighbor.HasProp("_TriposAtomName"):
                continue
            h_name = neighbor.GetProp("_TriposAtomName")
            if h_name in name_set and h_name not in added:
                ordered.append(h_name)
                added.add(h_name)
    return ordered


def compute_pdb_records(mol: Chem.Mol, names: Sequence[str], translation: dict[str, str]) -> list[PdbAtomRecord]:
    """Build PDB ATOM records for `names` (heavy-then-bonded-H order, no lone pairs)."""
    conformer = mol.GetConformer()
    name_to_atom = {
        atom.GetProp("_TriposAtomName"): atom
        for atom in mol.GetAtoms()
        if atom.HasProp("_TriposAtomName")
    }
    ordered_names = _heavy_then_bonded_h_order(mol, names)
    records = []
    for serial, name in enumerate(ordered_names, start=1):
        pos = conformer.GetAtomPosition(name_to_atom[name].GetIdx())
        records.append(
            PdbAtomRecord(
                serial=serial,
                atom_name=translation[name],
                resname=SEGID,
                x=pos.x, y=pos.y, z=pos.z,
            )
        )
    return records


# ---------------------------------------------------------------------------
# RTF atom-block records (core.rtf and site*_pres.rtf ATOM/LONEPAIR lines)
# ---------------------------------------------------------------------------


@dataclass
class AtomBlockEntry:
    translated_name: str
    atom_type: str
    charge: float
    is_lp: bool = False
    lp_record: LonePairRecord | None = None


def compute_core_atom_block(
    decomposition: CoreDecomposition,
    reference: LigandRecord,
    core_charges: dict[str, float],
    translation: AtomNameTranslation,
) -> list[AtomBlockEntry]:
    """Ported from msld_crn.py:833-847: core.rtf's ATOM/LONEPAIR block, reference-sourced."""
    name_to_atom = {
        atom.GetProp("_TriposAtomName"): atom
        for atom in reference.mol.GetAtoms()
        if atom.HasProp("_TriposAtomName")
    }
    entries = []
    for name in _heavy_then_bonded_h_order(reference.mol, decomposition.core_atom_order):
        atom_type = name_to_atom[name].GetProp("CGENFFAtomType")
        entries.append(AtomBlockEntry(translation.core[name], atom_type, core_charges[name]))

    lp_by_host_name = {
        lp.lp_atom: lp for lp in reference.lone_pairs
    }
    for lp_by_ligand in decomposition.core_lp_atom_names.values():
        ref_lp_name = lp_by_ligand[reference.name]
        lp = lp_by_host_name[ref_lp_name]
        atom_type = "LPH"
        entries.append(
            AtomBlockEntry(
                translation.core[ref_lp_name], atom_type, core_charges[ref_lp_name],
                is_lp=True, lp_record=lp,
            )
        )
    return entries


def compute_fragment_atom_block(
    site_index: int,
    ligand: LigandRecord,
    fragment: Fragment,
    fragment_charges: dict[str, float],
    translation: AtomNameTranslation,
) -> list[AtomBlockEntry]:
    """Ported from msld_crn.py:889-903: a fragment's PRES ATOM/LONEPAIR block."""
    name_to_atom = {
        atom.GetProp("_TriposAtomName"): atom
        for atom in ligand.mol.GetAtoms()
        if atom.HasProp("_TriposAtomName")
    }
    frag_translation = translation.fragments[(site_index, ligand.name)]
    entries = []
    for name in _heavy_then_bonded_h_order(ligand.mol, fragment.atom_names):
        atom_type = name_to_atom[name].GetProp("CGENFFAtomType")
        entries.append(AtomBlockEntry(frag_translation[name], atom_type, fragment_charges[name]))

    lp_by_host_name = {lp.lp_atom: lp for lp in ligand.lone_pairs}
    for lp_name in fragment.lp_atom_names:
        lp = lp_by_host_name[lp_name]
        entries.append(
            AtomBlockEntry(
                frag_translation[lp_name], "LPH", fragment_charges[lp_name],
                is_lp=True, lp_record=lp,
            )
        )
    return entries


# ---------------------------------------------------------------------------
# BOND records
# ---------------------------------------------------------------------------


def compute_core_bonds(decomposition: CoreDecomposition, reference: LigandRecord, translation: AtomNameTranslation) -> list[tuple[str, str]]:
    """Ported from msld_crn.py:857-859: core-internal BOND lines, translated names."""
    core_names = set(decomposition.core_atom_order)
    bonds = []
    for bond in reference.mol.GetBonds():
        a = bond.GetBeginAtom().GetProp("_TriposAtomName") if bond.GetBeginAtom().HasProp("_TriposAtomName") else None
        b = bond.GetEndAtom().GetProp("_TriposAtomName") if bond.GetEndAtom().HasProp("_TriposAtomName") else None
        if a in core_names and b in core_names:
            bonds.append((translation.core[a], translation.core[b]))
    return bonds


def compute_fragment_bonds(
    decomposition: CoreDecomposition,
    site_index: int,
    ligand: LigandRecord,
    fragment: Fragment,
    translation: AtomNameTranslation,
) -> list[tuple[str, str]]:
    """Ported from msld_crn.py:912-930: fragment-internal + core-to-fragment BOND lines."""
    local_core = core_translation_for_ligand(decomposition, translation, ligand.name)
    frag_translation = translation.fragments[(site_index, ligand.name)]
    frag_names = set(fragment.atom_names)
    bonds = []
    for bond in ligand.mol.GetBonds():
        begin, end = bond.GetBeginAtom(), bond.GetEndAtom()
        if not (begin.HasProp("_TriposAtomName") and end.HasProp("_TriposAtomName")):
            continue
        a, b = begin.GetProp("_TriposAtomName"), end.GetProp("_TriposAtomName")
        a_in = a in local_core or a in frag_names
        b_in = b in local_core or b in frag_names
        if not (a_in and b_in) or (a in local_core and b in local_core):
            continue
        t_a = resolve_translated_name(a, local_core_translation=local_core, fragment_translation=frag_translation)
        t_b = resolve_translated_name(b, local_core_translation=local_core, fragment_translation=frag_translation)
        if t_a is not None and t_b is not None:
            bonds.append((t_a, t_b))
    return bonds


# ---------------------------------------------------------------------------
# IMPR decisions (the cross-site logic)
# ---------------------------------------------------------------------------


@dataclass
class PresRtfImprDecision:
    atom_names: tuple[str, str, str, str]
    commented_out: bool = False
    reason: str | None = None


def compute_fragment_impr_decisions(
    decomposition: CoreDecomposition,
    site_index: int,
    ligand: LigandRecord,
    fragment: Fragment,
    translation: AtomNameTranslation,
) -> list[PresRtfImprDecision]:
    """Ported from msld_crn.py:931-990: decide whether each IMPR line is
    fragment-local (write live), or spans another site (write commented-out
    with a warning, or drop entirely if it can't be resolved via an anchor).
    """
    local_core = core_translation_for_ligand(decomposition, translation, ligand.name)
    frag_translation = translation.fragments[(site_index, ligand.name)]
    frag_names = set(fragment.atom_names)
    core_names = set(decomposition.core_atom_names[ligand.name])

    # this ligand's site_index -> anchor atom name (for identifying which
    # OTHER site an off-core, off-this-fragment IMPR atom belongs to)
    anchor_to_site: dict[str, int] = {}
    for site in decomposition.sites:
        other_fragment = site.fragments.get(ligand.name)
        if other_fragment and len(other_fragment.anchor_atom_names) == 1:
            anchor_to_site[other_fragment.anchor_atom_names[0]] = site.index

    decisions: list[PresRtfImprDecision] = []
    for impr in ligand.impropers:
        if not any(at in frag_names for at in impr):
            continue

        other_sites: set[int] = set()
        unresolvable = False
        for at in impr:
            if at in core_names or at in frag_names:
                continue
            if at in anchor_to_site:
                other_sites.add(anchor_to_site[at])
            else:
                unresolvable = True

        if unresolvable:
            logger.warning(
                "Cannot print cross-site IMPR %s for site %d fragment %s: involves a "
                "non-anchor atom of another site; dropped.",
                impr, site_index, ligand.name,
            )
            continue
        if len(other_sites) > 1:
            logger.warning(
                "Cannot print IMPR %s for site %d fragment %s: spans three or more sites; skipped.",
                impr, site_index, ligand.name,
            )
            continue

        if not other_sites:
            names = tuple(
                resolve_translated_name(at, local_core_translation=local_core, fragment_translation=frag_translation)
                for at in impr
            )
            decisions.append(PresRtfImprDecision(names))
            continue

        # cross-site: resolve using the LARGEST fragment at the other site (arbitrary
        # but deterministic choice, matching legacy's "for frag2 in range(len(frags[crossite]))"
        # loop, which in practice writes one commented block per other-site fragment --
        # here we do the same, once per other-site fragment).
        other_site = decomposition.sites[next(iter(other_sites))]
        logger.warning(
            "Cross-site IMPR found (site %d fragment %s, spans site %d) -- written as a "
            "comment; manually uncomment as needed.",
            site_index, ligand.name, other_site.index,
        )
        for other_ligand_name, other_fragment in other_site.fragments.items():
            other_frag_translation = translation.fragments[(other_site.index, other_ligand_name)]
            other_anchor = other_fragment.anchor_atom_names[0]
            names = []
            for at in impr:
                if at in core_names:
                    names.append(local_core[at])
                elif at in frag_names:
                    names.append(frag_translation[at])
                else:
                    names.append(other_frag_translation[other_anchor])
            decisions.append(
                PresRtfImprDecision(
                    tuple(names), commented_out=True,
                    reason=f"cross-site: site {site_index} x site {other_site.index}",
                )
            )
    return decisions


# ---------------------------------------------------------------------------
# JSON-serializable file-data bundles (req.: dump site_sub_frag.pdb + pres.rtf
# records to JSON, and a writer that rebuilds the files from that JSON alone)
# ---------------------------------------------------------------------------


@dataclass
class RtfAtomData:
    name: str
    type: str
    charge: float
    is_lp: bool = False
    lp_host1: str | None = None
    lp_host2: str | None = None
    lp_distance: float | None = None


@dataclass
class RtfImprData:
    atoms: list[str]
    commented_out: bool = False


@dataclass
class PdbAtomData:
    serial: int
    name: str
    x: float
    y: float
    z: float


@dataclass
class CoreFileData:
    """Everything needed to write core.rtf + core.pdb, JSON-safe (plain dict fields only)."""

    source_name: str
    format_version: list[int]
    resi_name: str
    net_charge: float
    mass_records: list[list[object]]  # [atom_type, mass, element]
    atoms: list[RtfAtomData]
    bonds: list[list[str]]
    impropers: list[list[str]]
    pdb_atoms: list[PdbAtomData]


@dataclass
class FragmentFileData:
    """Everything needed to write one site{N}_sub{M}_pres.rtf + _frag.pdb, JSON-safe."""

    site_index: int
    frag_number: int
    source_name: str
    format_version: list[int]
    pres_name: str
    net_charge: float
    atoms: list[RtfAtomData]
    bonds: list[list[str]]
    impropers: list[RtfImprData]
    pdb_atoms: list[PdbAtomData]


def _atom_block_to_data(entries: list[AtomBlockEntry]) -> list[RtfAtomData]:
    data = []
    for entry in entries:
        lp = entry.lp_record
        data.append(
            RtfAtomData(
                name=entry.translated_name,
                type=entry.atom_type,
                charge=entry.charge,
                is_lp=entry.is_lp,
                lp_host1=lp.host_atoms[0] if lp else None,
                lp_host2=lp.host_atoms[1] if lp else None,
                lp_distance=lp.distance if lp else None,
            )
        )
    return data


def _pdb_records_to_data(records: list[PdbAtomRecord]) -> list[PdbAtomData]:
    return [PdbAtomData(r.serial, r.atom_name, r.x, r.y, r.z) for r in records]


def compute_core_file_data(
    decomposition: CoreDecomposition,
    ligands: Sequence[LigandRecord],
    crn: CRNResult,
    translation: AtomNameTranslation,
    *,
    resi_name: str = SEGID,
) -> CoreFileData:
    """Compute the full JSON-safe record bundle for core.rtf + core.pdb."""
    from .pdb import compute_pdb_records  # local import: pdb.py doesn't depend on records' bottom half

    ligand_by_name = {lig.name: lig for lig in ligands}
    reference = ligand_by_name[decomposition.reference_name]

    atom_block = compute_core_atom_block(decomposition, reference, crn.core_charges, translation)
    bonds = compute_core_bonds(decomposition, reference, translation)
    core_names = set(decomposition.core_atom_order)
    impropers = [
        [translation.core[at] for at in impr]
        for impr in reference.impropers
        if all(at in core_names for at in impr)
    ]
    mass_records = _collect_core_mass_records(decomposition, ligands)
    pdb_records = compute_pdb_records(reference.mol, decomposition.core_atom_order, translation.core)

    return CoreFileData(
        source_name=reference.name,
        format_version=list(reference.format_version),
        resi_name=resi_name,
        net_charge=round(sum(crn.core_charges.values()), 6),
        mass_records=[[t, m, e] for t, m, e in mass_records],
        atoms=_atom_block_to_data(atom_block),
        bonds=[list(b) for b in bonds],
        impropers=impropers,
        pdb_atoms=_pdb_records_to_data(pdb_records),
    )


def compute_fragment_file_data(
    decomposition: CoreDecomposition,
    site_index: int,
    frag_number: int,
    ligand: LigandRecord,
    crn: CRNResult,
    translation: AtomNameTranslation,
) -> FragmentFileData:
    """Compute the full JSON-safe record bundle for one site{N}_sub{M}_pres.rtf + _frag.pdb."""
    from .pdb import compute_pdb_records

    site = decomposition.sites[site_index]
    fragment = site.fragments[ligand.name]
    fragment_charges = crn.fragment_charges[site_index][ligand.name]

    atom_block = compute_fragment_atom_block(site_index, ligand, fragment, fragment_charges, translation)
    bonds = compute_fragment_bonds(decomposition, site_index, ligand, fragment, translation)
    impr_decisions = compute_fragment_impr_decisions(decomposition, site_index, ligand, fragment, translation)
    frag_translation = translation.fragments[(site_index, ligand.name)]
    pdb_records = compute_pdb_records(ligand.mol, fragment.atom_names, frag_translation)

    return FragmentFileData(
        site_index=site_index,
        frag_number=frag_number,
        source_name=ligand.name,
        format_version=list(ligand.format_version),
        pres_name=f"p{site_index + 1}_{frag_number}",
        net_charge=round(sum(fragment_charges.values()), 6),
        atoms=_atom_block_to_data(atom_block),
        bonds=[list(b) for b in bonds],
        impropers=[
            RtfImprData(list(d.atom_names), d.commented_out) for d in impr_decisions
        ],
        pdb_atoms=_pdb_records_to_data(pdb_records),
    )


def _collect_core_mass_records(
    decomposition: CoreDecomposition, ligands: Sequence[LigandRecord]
) -> list[tuple[str, float, str]]:
    """Ported from writers/rtf.py's `_collect_mass_records` (kept here too so
    `compute_core_file_data` doesn't depend on writers/rtf.py, avoiding a
    circular import -- writers/rtf.py imports this module for its computers).
    """
    ligand_by_name = {lig.name: lig for lig in ligands}
    reference = ligand_by_name[decomposition.reference_name]
    ref_type_by_name = {
        atom.GetProp("_TriposAtomName"): atom.GetProp("CGENFFAtomType")
        for atom in reference.mol.GetAtoms()
        if atom.HasProp("_TriposAtomName")
    }

    seen: set[str] = set()
    records: list[tuple[str, float, str]] = []

    def _add(atom_type: str, source: LigandRecord) -> None:
        if atom_type in seen or atom_type not in source.mass_records:
            return
        seen.add(atom_type)
        mass_record = source.mass_records[atom_type]
        records.append((atom_type, mass_record.mass, mass_record.element))

    for name in decomposition.core_atom_order:
        _add(ref_type_by_name[name], reference)

    for site in decomposition.sites:
        for ligand_name, fragment in site.fragments.items():
            ligand = ligand_by_name[ligand_name]
            type_by_name = {
                atom.GetProp("_TriposAtomName"): atom.GetProp("CGENFFAtomType")
                for atom in ligand.mol.GetAtoms()
                if atom.HasProp("_TriposAtomName")
            }
            for name in fragment.atom_names:
                _add(type_by_name[name], ligand)

    return records
