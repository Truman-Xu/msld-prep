"""JSON round-trip for the per-fragment/core PDB+RTF records.

This is the "save site_sub_frag.pdb + pres.rtf records in dictionaries,
dump as JSON, and provide a writer that reproduces the files from that JSON
alone" checkpoint: `write_records_json` dumps everything
`compute_core_file_data`/`compute_fragment_file_data` (records.py) computed;
`write_files_from_json` rebuilds `core.rtf`/`core.pdb`/every
`site{N}_sub{M}_pres.rtf`/`_frag.pdb` purely from the loaded JSON, with no
dependency on the original RDKit Mol objects, `CoreDecomposition`, or
`CRNResult` still being in memory.
"""

from __future__ import annotations

import json
import logging
import os
from dataclasses import asdict
from typing import Sequence

from ..datamodel import CoreDecomposition, CRNResult, LigandRecord
from . import pdb as pdb_writer
from . import rtf as rtf_writer
from .records import (
    AtomNameTranslation,
    CoreFileData,
    FragmentFileData,
    PdbAtomData,
    RtfAtomData,
    RtfImprData,
    compute_core_file_data,
    compute_fragment_file_data,
)

logger = logging.getLogger(__name__)


def _fragment_key(site_index: int, frag_number: int) -> str:
    return f"site{site_index + 1}_sub{frag_number}"


def compute_all_file_data(
    decomposition: CoreDecomposition,
    ligands: Sequence[LigandRecord],
    crn: CRNResult,
    translation: AtomNameTranslation,
) -> tuple[CoreFileData, dict[str, FragmentFileData]]:
    core_data = compute_core_file_data(decomposition, ligands, crn, translation)
    ligand_by_name = {lig.name: lig for lig in ligands}
    fragment_data: dict[str, FragmentFileData] = {}
    for site in decomposition.sites:
        for frag_number, ligand_name in enumerate(site.fragments, start=1):
            data = compute_fragment_file_data(
                decomposition, site.index, frag_number, ligand_by_name[ligand_name], crn, translation,
            )
            fragment_data[_fragment_key(site.index, frag_number)] = data
    return core_data, fragment_data


def write_records_json(
    core_data: CoreFileData,
    fragment_data: dict[str, FragmentFileData],
    outdir: str,
    filename: str = "records.json",
) -> str:
    payload = {
        "core": asdict(core_data),
        "fragments": {key: asdict(data) for key, data in fragment_data.items()},
    }
    path = os.path.join(outdir, filename)
    with open(path, "w") as fp:
        json.dump(payload, fp, indent=2)
    logger.info("Wrote %s (%d fragment record(s))", path, len(fragment_data))
    return path


def _core_data_from_dict(d: dict) -> CoreFileData:
    return CoreFileData(
        source_name=d["source_name"],
        format_version=d["format_version"],
        resi_name=d["resi_name"],
        net_charge=d["net_charge"],
        mass_records=d["mass_records"],
        atoms=[RtfAtomData(**item) for item in d["atoms"]],
        bonds=d["bonds"],
        impropers=d["impropers"],
        pdb_atoms=[PdbAtomData(**item) for item in d["pdb_atoms"]],
    )


def _fragment_data_from_dict(d: dict) -> FragmentFileData:
    return FragmentFileData(
        site_index=d["site_index"],
        frag_number=d["frag_number"],
        source_name=d["source_name"],
        format_version=d["format_version"],
        pres_name=d["pres_name"],
        net_charge=d["net_charge"],
        atoms=[RtfAtomData(**item) for item in d["atoms"]],
        bonds=d["bonds"],
        impropers=[RtfImprData(**item) for item in d["impropers"]],
        pdb_atoms=[PdbAtomData(**item) for item in d["pdb_atoms"]],
    )


def read_records_json(path: str) -> tuple[CoreFileData, dict[str, FragmentFileData]]:
    with open(path) as fp:
        payload = json.load(fp)
    core_data = _core_data_from_dict(payload["core"])
    fragment_data = {key: _fragment_data_from_dict(value) for key, value in payload["fragments"].items()}
    return core_data, fragment_data


def write_files_from_json(path: str, outdir: str) -> None:
    """Rebuild core.rtf/core.pdb + every site{N}_sub{M}_pres.rtf/_frag.pdb from `path` alone."""
    core_data, fragment_data = read_records_json(path)

    rtf_writer.write_core_rtf_from_data(core_data, outdir)
    with open(os.path.join(outdir, "core.pdb"), "w") as fp:
        fp.write(pdb_writer.format_pdb_text(core_data.pdb_atoms))

    for data in fragment_data.values():
        rtf_writer.write_fragment_pres_rtf_from_data(data, outdir)
        frag_pdb_path = os.path.join(outdir, f"site{data.site_index + 1}_sub{data.frag_number}_frag.pdb")
        with open(frag_pdb_path, "w") as fp:
            fp.write(pdb_writer.format_pdb_text(data.pdb_atoms))

    logger.info("Rebuilt core + %d fragment file(s) from %s", len(fragment_data), path)
