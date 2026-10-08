"""CHARMM .prm parameter records: per-ligand parsing (via crimm's
PRMParser), merging across a ligand set, and deriving the parameters one
MSLD system's topology actually needs.

Parses into real tuples with typed (float/int) namedtuple values instead
of the opaque strings a hand-rolled parser would produce -- and, unlike a
hand-rolled parser, crimm's PRMParser already handles CHARMM line
continuation. `LigandParameters` (one file's worth) is what
`datamodel.LigandRecord.prm` holds; `MergedPrmLibrary`/`NeededParameters`
(built from many ligands' `LigandParameters` plus a `CoreDecomposition`)
feed `writers/prm.py`'s `full_ligand.prm` output.

`CoreDecomposition`/`LigandRecord` are only used here as type hints (never
constructed or introspected at runtime), imported under `TYPE_CHECKING` so
this module has no runtime dependency on `datamodel.py` -- `datamodel.py`
imports `LigandParameters`/`parse_prm_file` from here at runtime (for
`LigandRecord.prm`/`LigandRecord.from_files`), so a real top-level import
in the other direction would be circular.

Two deliberate fixes relative to legacy (msld_prm.py), per user decision:
1. Conflicting duplicate parameters across .prm files are now actually
   emitted as commented alternative lines in full_ligand.prm (legacy
   detected them, printed "UNcomment as necessary", and then never wrote
   the line).
2. The NONBONDED header reused in full_ligand.prm is now the first one
   encountered in a stable ligand order, not whichever .prm file happened
   to be parsed last (order-dependent in legacy).
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Callable, Iterable, Sequence

import numpy as np
from crimm.IO.PRMParser import (
    angle_par,
    bond_par,
    categorize_lines,
    cmap_par,
    dihe_par,
    impr_par,
    nbfix_par,
    nonbond14_par,
    nonbond_par,
    parse_line_dict,
    ub_par,
)

from .errors import PrmError
from .rtf_records import MassRecord

if TYPE_CHECKING:
    from .datamodel import CoreDecomposition, LigandRecord

logger = logging.getLogger(__name__)

BOND_RADIUS = 3  # angle/dihedral derivation walks this many bonds out


# ---------------------------------------------------------------------------
# Parsing: one .prm file -> structured records (no merge/decision logic)
# ---------------------------------------------------------------------------


@dataclass
class LigandParameters:
    bonds: dict[tuple[str, str], bond_par] = field(default_factory=dict)
    angles: dict[tuple[str, str, str], angle_par] = field(default_factory=dict)
    urey_bradley: dict[tuple[str, str, str], ub_par] = field(default_factory=dict)
    dihedrals: dict[tuple[str, str, str, str], list[dihe_par]] = field(default_factory=dict)
    impropers: dict[tuple[str, str, str, str], impr_par] = field(default_factory=dict)
    cmap: dict[tuple[str, ...], list[cmap_par]] = field(default_factory=dict)
    nonbonded: dict[str, nonbond_par] = field(default_factory=dict)
    nonbonded14: dict[str, nonbond14_par] = field(default_factory=dict)
    nbfix: dict[tuple[str, str], nbfix_par] = field(default_factory=dict)
    nonbonded_header: list[str] = field(default_factory=list)


def _read_nonbonded_header_from_lines(lines: Sequence[str]) -> list[str]:
    """The literal "NONBONDED ..." header line, plus any "-"-continued
    lines immediately after it (e.g. a following "cutnb ..." line) -- not
    exposed by crimm's structured PRMParser output, which only captures
    the per-atom-type entries below it, not the section's own
    nbxmod/cutnb/switching settings. Mirrors rtf_records.py's
    format-version/mass-record scans for data crimm's RTFParser likewise
    doesn't expose.
    """
    header: list[str] = []
    capturing = False
    for raw in lines:
        stripped = raw.strip()
        if not capturing:
            if stripped.upper().startswith("NONBONDED"):
                capturing = True
                header.append(raw.rstrip("\n") + "\n")
                if not stripped.endswith("-"):
                    break
            continue
        header.append(raw.rstrip("\n") + "\n")
        if not stripped.endswith("-"):
            break
    return header


def parse_prm_lines(lines: Iterable[str]) -> LigandParameters:
    """Parse .prm content from an iterable of lines (a file object, or
    `text.splitlines()`). Lines starting with "read " (case-insensitive)
    are dropped first -- `categorize_lines`'s `skip_line` helper only
    special-cases CHARMM's "read rtf" stream-wrapper line, not "read
    param", so a CGenFF unified .str stream's param block would otherwise
    reach `categorize_lines` before any section header has been seen.
    Harmless no-op for an ordinary standalone .prm file, which has no such
    line to begin with.
    """
    all_lines = list(lines)
    filtered = [line for line in all_lines if not line.strip().lower().startswith("read ")]
    parsed = parse_line_dict(categorize_lines(filtered))
    return LigandParameters(
        bonds=parsed["bonds"],
        angles=parsed["angles"],
        urey_bradley=parsed["urey_bradley"],
        dihedrals=parsed["dihedrals"],
        impropers=parsed["improper"],
        cmap=parsed["cmap"],
        nonbonded=parsed["nonbonded"],
        nonbonded14=parsed["nonbonded14"],
        nbfix=parsed["nbfix"],
        nonbonded_header=_read_nonbonded_header_from_lines(filtered),
    )


def parse_prm_text(text: str) -> LigandParameters:
    return parse_prm_lines(text.splitlines())


def parse_prm_file(path: str) -> LigandParameters:
    try:
        with open(path) as fp:
            return parse_prm_lines(fp)
    except OSError as exc:
        msg = f"Could not read parameter file `{path}`: {exc}"
        logger.error(msg)
        raise PrmError(msg) from exc


def load_prm_files(ligand_names: Sequence[str], search_dir: str) -> dict[str, LigandParameters]:
    """Parse `{search_dir}/{name}.prm` once per distinct ligand name."""
    parsed: dict[str, LigandParameters] = {}
    for name in ligand_names:
        if name in parsed:
            continue
        parsed[name] = parse_prm_file(os.path.join(search_dir, f"{name}.prm"))
    return parsed


# ---------------------------------------------------------------------------
# Line formatting -- shared between the duplicate-conflict comments below,
# `writers/prm.py`'s full_ligand.prm output, and `datamodel.LigandRecord`'s
# single-ligand round-trip writer, so all three stay in lockstep.
# ---------------------------------------------------------------------------


def format_bond_line(key: tuple[str, str], value: bond_par) -> str:
    return f"{key[0]} {key[1]} {value.kb:.4f} {value.b0:.4f}"


def format_angle_line(key: tuple[str, str, str], value: angle_par, ub: ub_par | None) -> str:
    line = f"{key[0]} {key[1]} {key[2]} {value.ktheta:.4f} {value.theta0:.4f}"
    if ub is not None:
        line += f" {ub.kub:.4f} {ub.s0:.4f}"
    return line


def format_dihedral_lines(key: tuple[str, str, str, str], terms: list[dihe_par]) -> list[str]:
    return [f"{key[0]} {key[1]} {key[2]} {key[3]} {t.kchi:.4f} {t.n} {t.delta:.4f}" for t in terms]


def format_improper_line(key: tuple[str, str, str, str], value: impr_par) -> str:
    return f"{key[0]} {key[1]} {key[2]} {key[3]} {value.kpsi:.4f} 0 {value.psi0:.4f}"


def format_nonbonded_line(atom_type: str, value: nonbond_par, value14: nonbond14_par | None) -> str:
    if value14 is not None:
        return f"{atom_type} 0.0 {value.epsilon:.4f} {value.rmin_half:.4f} 0.0 {value14.epsilon:.4f} {value14.rmin_half:.4f}"
    return f"{atom_type} 0.0 {value.epsilon:.4f} {value.rmin_half:.4f}"


# ---------------------------------------------------------------------------
# Merging across ligands' .prm files
# ---------------------------------------------------------------------------


@dataclass
class MergedPrmLibrary:
    bonds: dict[tuple[str, str], bond_par] = field(default_factory=dict)
    angles: dict[tuple[str, str, str], angle_par] = field(default_factory=dict)
    urey_bradley: dict[tuple[str, str, str], ub_par] = field(default_factory=dict)
    dihedrals: dict[tuple[str, str, str, str], list[dihe_par]] = field(default_factory=dict)
    impropers: dict[tuple[str, str, str, str], impr_par] = field(default_factory=dict)
    nonbonded: dict[str, nonbond_par] = field(default_factory=dict)
    nonbonded14: dict[str, nonbond14_par] = field(default_factory=dict)
    nonbonded_header: list[str] = field(default_factory=list)
    duplicates: dict[str, list[str]] = field(
        default_factory=lambda: {"BOND": [], "ANGL": [], "DIHE": [], "IMPR": [], "NONB": []}
    )


def _merge_typed_section(
    section_label: str,
    parsed_by_ligand: dict[str, dict],
    ligand_order: Sequence[str],
    format_duplicate: Callable[[object, object], list[str]],
) -> tuple[dict, list[str]]:
    """First-ligand-in-`ligand_order`-wins merge. Keys are compared against
    their reverse too (a bond/angle/dihedral/improper written "backwards"
    in a different ligand's file isn't a conflict); plain string keys
    (nonbonded, keyed by a single atom type) have no reverse, so `key[::-1]`
    is skipped for them.
    """
    merged: dict = {}
    duplicates: list[str] = []
    for ligand_name in ligand_order:
        for key, value in parsed_by_ligand.get(ligand_name, {}).items():
            rev = key[::-1] if isinstance(key, tuple) else key
            existing_key = key if key in merged else (rev if rev in merged else None)
            if existing_key is None:
                merged[key] = value
            elif merged[existing_key] != value:
                logger.warning(
                    "%s: duplicate parameter `%s` with a differing value from ligand %s; "
                    "keeping the first value found, recording this one as a commented alternative.",
                    section_label, key, ligand_name,
                )
                duplicates.extend(format_duplicate(key, value))
    return merged, duplicates


def merge_prm_records(
    prm_by_ligand: dict[str, LigandParameters], *, cgenff: bool, ligand_order: Sequence[str]
) -> MergedPrmLibrary:
    library = MergedPrmLibrary()

    library.bonds, library.duplicates["BOND"] = _merge_typed_section(
        "BOND", {n: p.bonds for n, p in prm_by_ligand.items()}, ligand_order,
        lambda k, v: [f"!{format_bond_line(k, v)}"],
    )
    library.angles, library.duplicates["ANGL"] = _merge_typed_section(
        "ANGLE", {n: p.angles for n, p in prm_by_ligand.items()}, ligand_order,
        lambda k, v: [f"!{format_angle_line(k, v, None)}"],
    )
    library.urey_bradley, _ = _merge_typed_section(
        "UREY-BRADLEY", {n: p.urey_bradley for n, p in prm_by_ligand.items()}, ligand_order,
        lambda k, v: [],  # rides along with its angle line; not separately reported as a conflict
    )
    library.dihedrals, library.duplicates["DIHE"] = _merge_typed_section(
        "DIHEDRAL", {n: p.dihedrals for n, p in prm_by_ligand.items()}, ligand_order,
        lambda k, v: [f"!{line}" for line in format_dihedral_lines(k, v)],
    )
    library.impropers, library.duplicates["IMPR"] = _merge_typed_section(
        "IMPROPER", {n: p.impropers for n, p in prm_by_ligand.items()}, ligand_order,
        lambda k, v: [f"!{format_improper_line(k, v)}"],
    )
    if not cgenff:
        library.nonbonded, library.duplicates["NONB"] = _merge_typed_section(
            "NONBONDED", {n: p.nonbonded for n, p in prm_by_ligand.items()}, ligand_order,
            lambda k, v: [f"!{format_nonbonded_line(k, v, None)}"],
        )
        library.nonbonded14, _ = _merge_typed_section(
            "NONBONDED14", {n: p.nonbonded14 for n, p in prm_by_ligand.items()}, ligand_order,
            lambda k, v: [],
        )
        for name in ligand_order:
            header = prm_by_ligand.get(name, LigandParameters()).nonbonded_header
            if header:
                library.nonbonded_header = header
                break

    return library


# ---------------------------------------------------------------------------
# System topology: which atoms/bonds/impropers exist across core + every fragment
# ---------------------------------------------------------------------------


@dataclass
class LigandTopology:
    atoms: list[tuple[str, str]] = field(default_factory=list)          # (unique_key, atom_type)
    bonds: list[tuple[str, str]] = field(default_factory=list)           # (unique_key, unique_key)
    impropers: list[tuple[str, str, str, str]] = field(default_factory=list)  # atom_types, order preserved
    mass_records: dict[str, MassRecord] = field(default_factory=dict)     # atom_type -> MassRecord


def _core_key(ref_name: str) -> str:
    return f"core:{ref_name}"


def _fragment_key(site_index: int, ligand_name: str, raw_name: str) -> str:
    return f"site{site_index}:{ligand_name}:{raw_name}"


def build_ligand_topology(
    decomposition: "CoreDecomposition",
    ligands: Sequence["LigandRecord"],
) -> LigandTopology:
    """Build the system-wide (core + every site's every fragment) atom/bond/
    improper graph. Real (non-lone-pair) atoms only -- lone pairs never
    appear in a BOND/ANGLE/DIHEDRAL record (CHARMM expresses them via
    LONEPAIR, not classical bonded terms), so they don't affect parameter
    derivation; MASS records are collected separately (`mass_records`)
    since that's independent of connectivity.
    """
    ligand_by_name = {lig.name: lig for lig in ligands}
    reference = ligand_by_name[decomposition.reference_name]
    ref_atom_by_name = {
        atom.GetProp("_TriposAtomName"): atom
        for atom in reference.mol.GetAtoms()
        if atom.HasProp("_TriposAtomName")
    }

    topology = LigandTopology()
    core_key_by_name: dict[str, str] = {}
    for ref_name in decomposition.core_atom_order:
        key = _core_key(ref_name)
        core_key_by_name[ref_name] = key
        atom_type = ref_atom_by_name[ref_name].GetProp("CGENFFAtomType")
        topology.atoms.append((key, atom_type))
        mass = reference.mass_records.get(atom_type)
        if mass is not None:
            topology.mass_records[atom_type] = mass

    for bond in reference.mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if not (a.HasProp("_TriposAtomName") and b.HasProp("_TriposAtomName")):
            continue
        a_name, b_name = a.GetProp("_TriposAtomName"), b.GetProp("_TriposAtomName")
        if a_name in core_key_by_name and b_name in core_key_by_name:
            topology.bonds.append((core_key_by_name[a_name], core_key_by_name[b_name]))

    core_names = set(decomposition.core_atom_order)
    for impr in reference.impropers:
        if all(at in core_names for at in impr):
            types = tuple(ref_atom_by_name[at].GetProp("CGENFFAtomType") for at in impr)
            topology.impropers.append(types)

    for site in decomposition.sites:
        for ligand_name, fragment in site.fragments.items():
            ligand = ligand_by_name[ligand_name]
            atom_by_name = {
                atom.GetProp("_TriposAtomName"): atom
                for atom in ligand.mol.GetAtoms()
                if atom.HasProp("_TriposAtomName")
            }
            local_core_names = decomposition.core_atom_names[ligand_name]
            local_core_key = {
                local_name: _core_key(decomposition.core_atom_order[pos])
                for pos, local_name in enumerate(local_core_names)
            }
            frag_key_by_name: dict[str, str] = {}
            for raw_name in fragment.atom_names:
                key = _fragment_key(site.index, ligand_name, raw_name)
                frag_key_by_name[raw_name] = key
                atom_type = atom_by_name[raw_name].GetProp("CGENFFAtomType")
                topology.atoms.append((key, atom_type))
                mass = ligand.mass_records.get(atom_type)
                if mass is not None:
                    topology.mass_records[atom_type] = mass

            def _key_for(name: str) -> str | None:
                if name in local_core_key:
                    return local_core_key[name]
                return frag_key_by_name.get(name)

            frag_name_set = set(fragment.atom_names)
            for bond in ligand.mol.GetBonds():
                a, b = bond.GetBeginAtom(), bond.GetEndAtom()
                if not (a.HasProp("_TriposAtomName") and b.HasProp("_TriposAtomName")):
                    continue
                a_name, b_name = a.GetProp("_TriposAtomName"), b.GetProp("_TriposAtomName")
                a_in_frag, b_in_frag = a_name in frag_name_set, b_name in frag_name_set
                if not (a_in_frag or b_in_frag):
                    continue
                key_a, key_b = _key_for(a_name), _key_for(b_name)
                if key_a is not None and key_b is not None:
                    topology.bonds.append((key_a, key_b))

            for impr in ligand.impropers:
                if not any(at in frag_name_set for at in impr):
                    continue
                types = []
                ok = True
                for at in impr:
                    if at in local_core_key:
                        types.append(ref_atom_by_name[
                            decomposition.core_atom_order[local_core_names.index(at)]
                        ].GetProp("CGENFFAtomType"))
                    elif at in frag_name_set:
                        types.append(atom_by_name[at].GetProp("CGENFFAtomType"))
                    else:
                        ok = False
                        break
                if ok:
                    topology.impropers.append(tuple(types))

    return topology


# ---------------------------------------------------------------------------
# Derive which merged parameters this topology actually needs
# ---------------------------------------------------------------------------


@dataclass
class NeededParameters:
    masses: list[MassRecord] = field(default_factory=list)
    bonds: dict[tuple[str, str], bond_par] = field(default_factory=dict)
    angles: dict[tuple[str, str, str], angle_par] = field(default_factory=dict)
    urey_bradley: dict[tuple[str, str, str], ub_par] = field(default_factory=dict)
    dihedrals: dict[tuple[str, str, str, str], list[dihe_par]] = field(default_factory=dict)
    impropers: dict[tuple[str, str, str, str], impr_par] = field(default_factory=dict)
    nonbonded: dict[str, nonbond_par] = field(default_factory=dict)
    nonbonded14: dict[str, nonbond14_par] = field(default_factory=dict)
    nonbonded_header: list[str] = field(default_factory=list)
    duplicates: dict[str, list[str]] = field(default_factory=dict)
    order: dict[str, list] = field(default_factory=dict)  # section -> keys, discovery order


def _log_missing(cgenff: bool, message: str) -> None:
    if not cgenff:
        logger.warning(message)
    else:
        logger.info(message)


def derive_needed_parameters(
    topology: LigandTopology, merged: MergedPrmLibrary, *, cgenff: bool
) -> NeededParameters:
    """3-bond-radius connectivity walk deciding which merged BOND/ANGLE/
    DIHEDRAL/IMPROPER/NONBONDED entries this specific system's topology
    actually uses. Missing-parameter messages are logged at WARNING when
    cgenff=False (an explicit force field is expected to be complete) and
    at INFO when cgenff=True (CGenFF/ParamChem streams and the base CHARMM
    CGenFF .prm commonly supply these separately, so a gap here isn't
    necessarily a problem).
    """
    atom_type = dict(topology.atoms)
    adjacency: dict[str, list[str]] = {key: [] for key, _ in topology.atoms}
    for a, b in topology.bonds:
        adjacency.setdefault(a, []).append(b)
        adjacency.setdefault(b, []).append(a)

    needed = NeededParameters(duplicates=merged.duplicates)
    order: dict[str, list] = {"BOND": [], "ANGL": [], "DIHE": [], "IMPR": [], "NONB": []}

    def _lookup(section: dict, key: tuple, target: dict, order_list: list) -> bool:
        rev = key[::-1]
        if key in section and key not in target:
            target[key] = section[key]
            order_list.append(key)
            return True
        if rev in section and rev not in target:
            target[rev] = section[rev]
            order_list.append(rev)
            return True
        return key in target or rev in target

    for at1, _ in topology.atoms:
        for at2 in adjacency.get(at1, []):
            t1, t2 = atom_type[at1], atom_type[at2]
            key = (t1, t2)
            if not _lookup(merged.bonds, key, needed.bonds, order["BOND"]):
                _log_missing(cgenff, f"No BOND parameter found for {at1}-{at2} ({t1} {t2})")

            if t2 in merged.nonbonded and t2 not in needed.nonbonded:
                needed.nonbonded[t2] = merged.nonbonded[t2]
                if t2 in merged.nonbonded14:
                    needed.nonbonded14[t2] = merged.nonbonded14[t2]
                order["NONB"].append(t2)

            for at3 in adjacency.get(at2, []):
                if at3 == at1:
                    continue
                t3 = atom_type[at3]
                key = (t1, t2, t3)
                if not _lookup(merged.angles, key, needed.angles, order["ANGL"]):
                    _log_missing(cgenff, f"No ANGLE parameter found for {at1}-{at2}-{at3} ({t1} {t2} {t3})")
                rev = key[::-1]
                if key in merged.urey_bradley:
                    needed.urey_bradley[key] = merged.urey_bradley[key]
                elif rev in merged.urey_bradley:
                    needed.urey_bradley[rev] = merged.urey_bradley[rev]

                for at4 in adjacency.get(at3, []):
                    if at4 == at2:
                        continue
                    t4 = atom_type[at4]
                    key = (t1, t2, t3, t4)
                    if not _lookup(merged.dihedrals, key, needed.dihedrals, order["DIHE"]):
                        _log_missing(
                            cgenff, f"No DIHEDRAL parameter found for {at1}-{at2}-{at3}-{at4} ({t1} {t2} {t3} {t4})"
                        )

    for impr in topology.impropers:
        key = tuple(impr)
        if not _lookup(merged.impropers, key, needed.impropers, order["IMPR"]):
            _log_missing(cgenff, f"No IMPROPER parameter found for {'-'.join(impr)} ({' '.join(impr)})")

    needed.masses = list(topology.mass_records.values())
    needed.nonbonded_header = merged.nonbonded_header
    needed.order = order
    return needed


# ---------------------------------------------------------------------------
# NBFIX combining rule
# ---------------------------------------------------------------------------


def ab2re(eps1: float, rmin1: float, eps2: float, rmin2: float) -> tuple[float, float]:
    """A/B (attractive/repulsive LJ coefficient) geometric-mean combining rule."""
    if rmin1 == 0 or rmin2 == 0:
        return 0.0, 0.0

    rmin1, rmin2 = rmin1 * 2, rmin2 * 2
    a_ii, a_jj = abs(eps1 * rmin1 ** 12) / 4, abs(eps2 * rmin2 ** 12) / 4
    b_ii, b_jj = abs(eps1 * rmin1 ** 6) / 2, abs(eps2 * rmin2 ** 6) / 2
    a_ij, b_ij = np.sqrt(a_ii * a_jj), np.sqrt(b_ii * b_jj)
    r_ij = (2 * a_ij / b_ij) ** (1 / 6)
    eps_ij = 2 * b_ij / (r_ij ** 6)
    return float(r_ij), float(eps_ij)


def compute_nbfix(
    nonbonded: dict[str, nonbond_par], nonbonded14: dict[str, nonbond14_par] | None = None
) -> dict[tuple[str, str], tuple[float, float, float, float]]:
    """NBFIX combining-rule parameters for every unordered pair of atom types
    (including self-pairs) in `nonbonded`, using each type's regular and
    (if present, else 0.0/0.0) 1-4 Lennard-Jones terms. Pairs with rmin==0
    are omitted.
    """
    nonbonded14 = nonbonded14 or {}
    types = list(nonbonded.keys())
    nbfix: dict[tuple[str, str], tuple[float, float, float, float]] = {}
    for i, type1 in enumerate(types):
        for type2 in types[i:]:
            eps1, rmin1 = nonbonded[type1].epsilon, nonbonded[type1].rmin_half
            eps2, rmin2 = nonbonded[type2].epsilon, nonbonded[type2].rmin_half
            nb14_1, nb14_2 = nonbonded14.get(type1), nonbonded14.get(type2)
            eps14_1, rmin14_1 = (nb14_1.epsilon, nb14_1.rmin_half) if nb14_1 else (0.0, 0.0)
            eps14_2, rmin14_2 = (nb14_2.epsilon, nb14_2.rmin_half) if nb14_2 else (0.0, 0.0)
            reg_r, reg_eps = ab2re(-eps1, rmin1, -eps2, rmin2)
            v14_r, v14_eps = ab2re(-eps14_1, rmin14_1, -eps14_2, rmin14_2)
            if reg_r == 0:
                continue
            nbfix[(type1, type2)] = (-reg_eps, reg_r, -v14_eps, v14_r)
    return nbfix
