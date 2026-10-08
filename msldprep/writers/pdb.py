"""PDB writers: core.pdb, site{N}_sub{M}_frag.pdb, large_lig.pdb.

Ported from msld_crn.py:597-621 (core), :713-737 (fragment), :761-785
(large_lig). The large_lig writer no longer re-reads core.pdb/fragment
pdb files off disk (msld_crn.py:765,779) -- it reuses the same in-memory
record lists write_core_pdb/write_fragment_pdb already built.
"""

from __future__ import annotations

import os
from typing import Sequence

from ..datamodel import CoreDecomposition, LigandRecord
from .records import SEGID, AtomNameTranslation, PdbAtomData, PdbAtomRecord, compute_pdb_records

_PDB_LINE = "ATOM  %5d %-4s %4s%5d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n"


def format_pdb_atom_line(record: PdbAtomRecord) -> str:
    return _PDB_LINE % (
        record.serial, record.atom_name, record.resname, 1,
        record.x, record.y, record.z, record.occ, record.b, record.segid,
    )


def format_pdb_text(atoms: Sequence[PdbAtomData]) -> str:
    """Format a JSON-safe `PdbAtomData` list -- no Mol/translation dependency,
    so `writers/fragment_json.py` can reuse this after a JSON round-trip.
    """
    lines = [
        _PDB_LINE % (a.serial, a.name, SEGID, 1, a.x, a.y, a.z, 1.0, 0.0, SEGID) for a in atoms
    ]
    lines.append("TER\nEND")
    return "".join(lines)


def write_pdb_file(records: Sequence[PdbAtomRecord], path: str) -> None:
    with open(path, "w") as fp:
        for record in records:
            fp.write(format_pdb_atom_line(record))
        fp.write("TER\nEND")


def write_core_pdb(
    decomposition: CoreDecomposition,
    reference: LigandRecord,
    translation: AtomNameTranslation,
    outdir: str,
) -> list[PdbAtomRecord]:
    records = compute_pdb_records(reference.mol, decomposition.core_atom_order, translation.core)
    write_pdb_file(records, os.path.join(outdir, "core.pdb"))
    return records


def write_fragment_pdb(
    site_index: int,
    frag_number: int,
    ligand: LigandRecord,
    fragment_atom_names: Sequence[str],
    translation: AtomNameTranslation,
    outdir: str,
) -> list[PdbAtomRecord]:
    frag_translation = translation.fragments[(site_index, ligand.name)]
    records = compute_pdb_records(ligand.mol, fragment_atom_names, frag_translation)
    path = os.path.join(outdir, f"site{site_index + 1}_sub{frag_number}_frag.pdb")
    write_pdb_file(records, path)
    return records


def write_large_lig_pdb(
    core_records: Sequence[PdbAtomRecord],
    largest_fragment_records_by_site: Sequence[list[PdbAtomRecord]],
    outdir: str,
) -> None:
    with open(os.path.join(outdir, "large_lig.pdb"), "w") as fp:
        for record in core_records:
            fp.write(format_pdb_atom_line(record))
        for records in largest_fragment_records_by_site:
            for record in records:
                fp.write(format_pdb_atom_line(record))
        fp.write("TER\nEND")


def write_all_pdb_files(
    decomposition: CoreDecomposition,
    ligands: Sequence[LigandRecord],
    translation: AtomNameTranslation,
    outdir: str,
    *,
    write_large_lig: bool = True,
) -> None:
    """Convenience: write core.pdb, every site's fragment pdbs, and (optionally) large_lig.pdb."""
    ligand_by_name = {lig.name: lig for lig in ligands}
    reference = ligand_by_name[decomposition.reference_name]
    core_records = write_core_pdb(decomposition, reference, translation, outdir)

    largest_by_site: list[list[PdbAtomRecord]] = []
    for site in decomposition.sites:
        largest_records: list[PdbAtomRecord] = []
        for frag_number, (ligand_name, fragment) in enumerate(site.fragments.items(), start=1):
            records = write_fragment_pdb(
                site.index, frag_number, ligand_by_name[ligand_name],
                fragment.atom_names, translation, outdir,
            )
            if len(records) > len(largest_records):
                largest_records = records
        largest_by_site.append(largest_records)

    if write_large_lig:
        write_large_lig_pdb(core_records, largest_by_site, outdir)
