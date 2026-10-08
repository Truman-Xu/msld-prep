"""CHARMM .rtf writers: core.rtf and site{N}_sub{M}_pres.rtf.

Ported from msld_crn.py:791-992. The BOND/IMPR core-vs-fragment resolution
and cross-site IMPR decision logic lives in `writers/records.py`; this
module is split into "compute a JSON-safe record bundle"
(`compute_core_file_data`/`compute_fragment_file_data`, in records.py) and
"format that bundle as CHARMM text" (`format_core_rtf_text`/
`format_fragment_pres_rtf_text`, here) -- the format half has no dependency
on RDKit Mol objects or the decomposition, so `writers/fragment_json.py` can
reuse it after a JSON round-trip.
"""

from __future__ import annotations

import os
from typing import Sequence

from ..datamodel import CoreDecomposition, CRNResult, LigandRecord
from .records import (
    AtomNameTranslation,
    CoreFileData,
    FragmentFileData,
    RtfAtomData,
    compute_core_file_data,
    compute_fragment_file_data,
)


def _format_atom_line(atom: RtfAtomData) -> str:
    return "ATOM %-4s %-6s %10.6f \n" % (atom.name, atom.type, atom.charge)


def format_core_rtf_text(data: CoreFileData) -> str:
    lines = [
        "* ligand core rtf file generated with msldprep for MSLD (JV,LC)\n",
        f"* (core from {data.source_name})\n* \n",
        f"  {data.format_version[0]} {data.format_version[1]}\n",
    ]
    for atom_type, mass, element in data.mass_records:
        lines.append("MASS -1 %s %s %s\n" % (atom_type, mass, element))
    lines.append("\n")
    lines.append("RESI  %s    %5.3f\n" % (data.resi_name, data.net_charge))
    lines.append("GROUP \n")
    for atom in data.atoms:
        lines.append(_format_atom_line(atom))
    for a, b in data.bonds:
        lines.append("BOND %-4s %-4s\n" % (a, b))
    for impr in data.impropers:
        lines.append("IMPR %-4s %-4s %-4s %-4s\n" % tuple(impr))
    lines.append("PATCH FIRST NONE LAST NONE\n\nEND")
    return "".join(lines)


def format_fragment_pres_rtf_text(data: FragmentFileData) -> str:
    lines = [
        "* fragment patch rtf file generated with msldprep for MSLD (JV,LC)\n",
        f"* (fragment from {data.source_name})\n* \n",
        f"  {data.format_version[0]} {data.format_version[1]}\n\n",
        "PRES %s    %5.3f\n" % (data.pres_name, data.net_charge),
        "GROUP \n",
    ]
    for atom in data.atoms:
        lines.append(_format_atom_line(atom))
    for a, b in data.bonds:
        lines.append("BOND %-4s %-4s\n" % (a, b))
    for impr in data.impropers:
        prefix = "!IMPR" if impr.commented_out else "IMPR"
        lines.append(f"{prefix} " + " ".join("%-4s" % name for name in impr.atoms) + "\n")
    lines.append("\nEND")
    return "".join(lines)


def write_core_rtf_from_data(data: CoreFileData, outdir: str) -> None:
    with open(os.path.join(outdir, "core.rtf"), "w") as fp:
        fp.write(format_core_rtf_text(data))


def write_fragment_pres_rtf_from_data(data: FragmentFileData, outdir: str) -> None:
    path = os.path.join(outdir, f"site{data.site_index + 1}_sub{data.frag_number}_pres.rtf")
    with open(path, "w") as fp:
        fp.write(format_fragment_pres_rtf_text(data))


def write_core_rtf(
    decomposition: CoreDecomposition,
    ligands: Sequence[LigandRecord],
    crn: CRNResult,
    translation: AtomNameTranslation,
    outdir: str,
) -> CoreFileData:
    data = compute_core_file_data(decomposition, ligands, crn, translation)
    write_core_rtf_from_data(data, outdir)
    return data


def write_fragment_pres_rtf(
    decomposition: CoreDecomposition,
    site_index: int,
    frag_number: int,
    ligand: LigandRecord,
    crn: CRNResult,
    translation: AtomNameTranslation,
    outdir: str,
) -> FragmentFileData:
    data = compute_fragment_file_data(decomposition, site_index, frag_number, ligand, crn, translation)
    write_fragment_pres_rtf_from_data(data, outdir)
    return data


def write_all_pres_rtf(
    decomposition: CoreDecomposition,
    ligands: Sequence[LigandRecord],
    crn: CRNResult,
    translation: AtomNameTranslation,
    outdir: str,
) -> dict[tuple[int, str], FragmentFileData]:
    ligand_by_name = {lig.name: lig for lig in ligands}
    results: dict[tuple[int, str], FragmentFileData] = {}
    for site in decomposition.sites:
        for frag_number, ligand_name in enumerate(site.fragments, start=1):
            results[(site.index, ligand_name)] = write_fragment_pres_rtf(
                decomposition, site.index, frag_number, ligand_by_name[ligand_name],
                crn, translation, outdir,
            )
    return results
