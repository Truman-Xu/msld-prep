"""Pure formatting: `NeededParameters` -> full_ligand.prm text.

Ported from msld_prm.py:509-545. No decisions here -- `prm_records.py`
already decided what's needed; this only serializes it. Independently
callable, writes `outdir/full_ligand.prm` unchanged in name/format so
downstream consumers (the generated ALF `.inp` script) keep working.
"""

from __future__ import annotations

import logging
import os

from ..prm_records import (
    NeededParameters,
    compute_nbfix,
    format_angle_line,
    format_bond_line,
    format_dihedral_lines,
    format_improper_line,
    format_nonbonded_line,
)

logger = logging.getLogger(__name__)

_DEFAULT_NONBONDED_HEADER = [
    "NONBONDED nbxmod 5 atom cdiel switch vatom vdistance vswitch -\n",
    "cutnb 14.0 ctofnb 12.0 ctonnb 11.5 eps 1.0 e14fac 0.5  geom\n",
]


def write_full_ligand_prm(needed: NeededParameters, outdir: str, *, cgenff: bool) -> None:
    lines = ["* MSLD ligand prm file generated with msldprep (JV,LC)\n* \n\n"]

    if not cgenff:
        lines.append("ATOMS\n")
        for mass in needed.masses:
            lines.append(f"MASS -1 {mass.atom_type} {mass.mass} {mass.element}\n")
        lines.append("\n")

    lines.append("BONDS\n")
    for key in needed.order.get("BOND", []):
        lines.append(format_bond_line(key, needed.bonds[key]) + "\n")
    for line in needed.duplicates.get("BOND", []):
        lines.append(line + "\n")

    lines.append("\nANGLES\n")
    for key in needed.order.get("ANGL", []):
        lines.append(format_angle_line(key, needed.angles[key], needed.urey_bradley.get(key)) + "\n")
    for line in needed.duplicates.get("ANGL", []):
        lines.append(line + "\n")

    lines.append("\nDIHEDRALS\n")
    for key in needed.order.get("DIHE", []):
        for line in format_dihedral_lines(key, needed.dihedrals[key]):
            lines.append(line + "\n")
    for line in needed.duplicates.get("DIHE", []):
        lines.append(line + "\n")

    lines.append("\nIMPROPERS\n")
    for key in needed.order.get("IMPR", []):
        lines.append(format_improper_line(key, needed.impropers[key]) + "\n")
    for line in needed.duplicates.get("IMPR", []):
        lines.append(line + "\n")

    if needed.nonbonded:
        lines.append("\n")
        header = needed.nonbonded_header or _DEFAULT_NONBONDED_HEADER
        lines.extend(header)
        for key in needed.order.get("NONB", []):
            lines.append(format_nonbonded_line(key, needed.nonbonded[key], needed.nonbonded14.get(key)) + "\n")
        for line in needed.duplicates.get("NONB", []):
            lines.append(line + "\n")

    lines.append("\nEND")

    path = os.path.join(outdir, "full_ligand.prm")
    with open(path, "w") as fp:
        fp.writelines(lines)

    if not cgenff:
        _append_nbfix(path, needed.nonbonded, needed.nonbonded14)
    logger.info("Wrote %s", path)


def _append_nbfix(path: str, nonbonded, nonbonded14) -> None:
    """Append an NBFIX section before the closing END, computed from the
    NONBONDED entries just written. Ported from msld_prm.py's `addNBFIX`
    (65-89), operating on the in-memory `nonbonded` dict directly instead
    of re-reading/re-parsing the file it just wrote.
    """
    if not nonbonded:
        return
    nbfix = compute_nbfix(nonbonded, nonbonded14)
    if not nbfix:
        return

    with open(path) as fp:
        lines = fp.readlines()
    with open(path, "w") as fp:
        for i, line in enumerate(lines):
            if "END" not in line and i != len(lines) - 1:
                fp.write(line)
        fp.write("\nNBFIX\n")
        for (type1, type2), (reg_eps, reg_r, v14_eps, v14_r) in nbfix.items():
            fp.write(f"{type1} {type2} {reg_eps:.6f} {reg_r:.6f} {v14_eps:.6f} {v14_r:.6f}\n")
        fp.write("\nEND")
