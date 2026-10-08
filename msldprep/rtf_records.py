"""Parse the parts of a CHARMM .rtf file RDKit Mol objects can't carry.

Per-atom CGenFF type/charge is stamped directly onto RDKit atoms
(`assign_atom_types`, called from `datamodel.LigandRecord.from_files`/
`analog_prep.mol_to_ligand_record`). Everything else parsed here --
residue net charge, bonds, IMPR (improper dihedral) declarations, and
LONEPAIR pseudo-atom definitions (lone pairs never appear in the .mol2 /
RDKit Mol at all) -- lands as direct fields on `LigandRecord` itself; there
is no separate "rtf aux" wrapper object. `parse_rtf_text`/`parse_rtf_file`
return a plain dict keyed exactly by those `LigandRecord` field names, so
callers just splat it in: `LigandRecord(..., **parse_rtf_file(path))`.

`crimm.IO.RTFParser` already exposes most of this structured (confirmed by
inspection: `topo_dict[...]` has 'total_charge', 'bonds', 'impropers',
'lonepairs' keys); the rtf format-version line and MASS statements are the
one bit of raw text scanning this module still does itself.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Iterable

from crimm.IO.RTFParser import RTFParser
from rdkit import Chem

from .errors import AlignmentError

logger = logging.getLogger(__name__)


@dataclass
class LonePairRecord:
    lp_type: str
    lp_atom: str
    host_atoms: tuple[str, ...]
    distance: float
    angle: float
    dihedral: float
    scale: float


@dataclass
class MassRecord:
    atom_type: str
    mass: float
    element: str


def _read_format_version_from_lines(lines: Iterable[str]) -> tuple[int, int]:
    """The `<major> <minor>` rtf format version line right after the title
    comments (e.g. "36 1"). Not exposed by crimm's structured output, so
    this is the one bit of raw text scanning in this module.
    """
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("*"):
            continue
        parts = stripped.split()
        if len(parts) >= 2 and all(p.lstrip("-").isdigit() for p in parts[:2]):
            return int(parts[0]), int(parts[1])
        break
    logger.debug("Could not find an rtf format version line; defaulting to 36 1.")
    return (36, 1)


def _read_mass_records_from_lines(lines: Iterable[str]) -> dict[str, MassRecord]:
    """MASS statements this .rtf file declares for itself, keyed by atom type.

    Not exposed by crimm's structured output. Many CGenFF-generated
    per-ligand rtf files declare none at all (relying on the master force
    field file for standard types) -- that's expected, not an error.
    """
    records: dict[str, MassRecord] = {}
    for line in lines:
        if not line.startswith("MASS"):
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        _, _, atom_type, mass, element = parts[:5]
        records[atom_type] = MassRecord(atom_type=atom_type, mass=float(mass), element=element)
    return records


def parse_rtf_text(rtf_text: str, *, source_label: str = "<rtf text>") -> dict:
    """Parse a single-residue CHARMM .rtf file's net charge, atom
    types/charges, bonds, IMPR, and LONEPAIR records straight from text --
    no file required (`crimm.IO.RTFParser.RTFParser` natively supports an
    in-memory `rtf_block`).

    Returns a dict keyed exactly by the `LigandRecord` fields it fills:
    `net_charge`, `format_version`, `atom_types`, `atom_charges`, `bonds`,
    `impropers`, `lone_pairs`, `mass_records` -- meant to be splatted
    directly into a `LigandRecord(...)` call.
    """
    rtf = RTFParser(rtf_block=rtf_text)
    if len(rtf.topo_dict) != 1:
        msg = (
            f"Expected exactly one residue in rtf `{source_label}`; "
            f"found {len(rtf.topo_dict)}: {list(rtf.topo_dict.keys())}."
        )
        logger.error(msg)
        raise AlignmentError(msg)

    data = next(iter(rtf.topo_dict.values()))

    atom_type_dict = next(iter(data.get("atoms", {}).values()), {})
    atom_types = {name: rec["atom_type"] for name, rec in atom_type_dict.items()}
    atom_charges = {name: rec["charge"] for name, rec in atom_type_dict.items()}

    bonds = tuple(
        (a, b)
        for order in ("single", "double", "triple", "aromatic")
        for a, b in data.get("bonds", {}).get(order, [])
    )
    impropers = tuple(tuple(impr) for impr in data.get("impropers", []))
    lone_pairs = tuple(
        LonePairRecord(
            lp_type=lp["type"],
            lp_atom=lp["lp_atom"],
            host_atoms=tuple(lp["host_atoms"]),
            distance=lp["distance"],
            angle=lp["angle"],
            dihedral=lp["dihedral"],
            scale=lp["scale"],
        )
        for lp in data.get("lonepairs", [])
    )

    logger.debug(
        "Parsed rtf from %s: net_charge=%s, %d bonds, %d impropers, %d lone pairs",
        source_label, data["total_charge"], len(bonds), len(impropers), len(lone_pairs),
    )

    lines = rtf_text.splitlines()
    return {
        "net_charge": float(data["total_charge"]),
        "format_version": _read_format_version_from_lines(lines),
        "atom_types": atom_types,
        "atom_charges": atom_charges,
        "bonds": bonds,
        "impropers": impropers,
        "lone_pairs": lone_pairs,
        "mass_records": _read_mass_records_from_lines(lines),
    }


def parse_rtf_file(rtf_path: str) -> dict:
    """Parse a single-residue CHARMM .rtf file's net charge, bonds, IMPR, and LONEPAIR records."""
    with open(rtf_path) as fp:
        rtf_text = fp.read()
    return parse_rtf_text(rtf_text, source_label=rtf_path)


def assign_atom_types(
    mol: Chem.Mol, atom_types: dict[str, str], atom_charges: dict[str, float], *, source_label: str = "<rtf>"
) -> None:
    """Stamp CGENFFAtomType/CGENFFPartialCharge onto every atom of `mol`
    from already-parsed `atom_types`/`atom_charges` dicts (e.g. from
    `parse_rtf_file`/`parse_rtf_text`, or a `LigandRecord`'s own fields) --
    no re-parsing.

    `mol` must carry a `_TriposAtomName` prop per atom (either set by
    RDKit's own mol2 reader, or by `sybyl_typer.assign_unique_atom_names`
    for an in-memory-built mol). Matches atoms by name.
    """
    for atom in mol.GetAtoms():
        if not atom.HasProp("_TriposAtomName"):
            msg = (
                "Atom does not have a Tripos atom name; use a mol loaded "
                "from a mol2 source (or named via sybyl_typer.assign_unique_atom_names) "
                "for atom type assignment."
            )
            logger.error(msg)
            raise AlignmentError(msg)
        atom_name = atom.GetProp("_TriposAtomName").upper()
        if atom_name not in atom_types:
            msg = f"Atom name {atom_name} does not exist in rtf `{source_label}`."
            logger.error(msg)
            raise AlignmentError(msg)
        atom.SetProp("CGENFFAtomType", atom_types[atom_name])
        atom.SetDoubleProp("CGENFFPartialCharge", atom_charges[atom_name])
