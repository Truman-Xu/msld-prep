import os
import sys

import pytest

from msldprep.writers import charmm_control
from msldprep.writers import pdb as pdb_writer
from msldprep.writers import rtf as rtf_writer

from .conftest import CHEMBL_DIR, LEGACY_REPO_ROOT


@pytest.fixture()
def written_outdir(tmp_path, decomposition, all_ligands, crn_result, translation):
    outdir = str(tmp_path / "crn_out")
    os.makedirs(outdir, exist_ok=True)
    pdb_writer.write_all_pdb_files(decomposition, all_ligands, translation, outdir)
    rtf_writer.write_core_rtf(decomposition, all_ligands, crn_result, translation, outdir)
    rtf_writer.write_all_pres_rtf(decomposition, all_ligands, crn_result, translation, outdir)
    charmm_control.write_lpsites_inp(decomposition, all_ligands, crn_result.core_charges, translation, outdir)
    charmm_control.write_translation_table(decomposition, translation, outdir)
    charmm_control.write_block_counts(decomposition, outdir)
    return outdir


def test_expected_files_exist(written_outdir, decomposition):
    names = set(os.listdir(written_outdir))
    expected = {"core.pdb", "core.rtf", "large_lig.pdb", "lpsites.inp", "translation.txt", "nsubs", "nblocks", "nreps"}
    assert expected <= names
    for site in decomposition.sites:
        for frag_number in range(1, len(site.fragments) + 1):
            assert f"site{site.index + 1}_sub{frag_number}_frag.pdb" in names
            assert f"site{site.index + 1}_sub{frag_number}_pres.rtf" in names


def test_core_rtf_atom_charges_match_crn_result(written_outdir, crn_result, translation):
    text = open(os.path.join(written_outdir, "core.rtf")).read()
    translated_to_ref = {v: k for k, v in translation.core.items()}
    found = 0
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        parts = line.split()
        translated_name, charge = parts[1], float(parts[3])
        ref_name = translated_to_ref[translated_name]
        if ref_name in crn_result.core_charges:
            assert round(charge, 6) == round(crn_result.core_charges[ref_name], 6)
            found += 1
    assert found > 0


def test_downstream_msld_prm_contract(written_outdir, monkeypatch):
    """The real backward-compat check: feed our output into the unmodified,
    out-of-scope msld_prm.MsldPRM and confirm it runs without error -- this
    validates the file-format contract against the actual downstream
    consumer, not a re-implemented mock of it.
    """
    sys.path.insert(0, LEGACY_REPO_ROOT)
    import msld_prm

    monkeypatch.chdir(CHEMBL_DIR)
    msld_prm.MsldPRM(written_outdir, cgenff=True, verbose=False, debug=False)
    assert os.path.exists(os.path.join(written_outdir, "full_ligand.prm"))


def _prm_section_entries(path: str) -> dict[str, int]:
    sections = ("BONDS", "ANGLES", "DIHEDRALS", "IMPROPERS", "NONBONDED")
    counts: dict[str, int] = {}
    current = None
    with open(path) as handle:
        for line in handle:
            words = line.split()
            if not words or words[0].startswith("!"):
                continue
            if words[0] in sections:
                current = words[0]
            elif words[0] == "END":
                current = None
            elif current:
                counts[current] = counts.get(current, 0) + 1
    return counts


@pytest.mark.xfail(
    int(__import__("pandas").__version__.split(".")[0]) >= 3,
    strict=True,
    reason="msld_prm fills its bond matrix by chained assignment, which writes nothing under pandas 3, "
    "so it finds no bonded terms (#17)",
)
def test_downstream_msld_prm_finds_bonded_terms(written_outdir, monkeypatch):
    sys.path.insert(0, LEGACY_REPO_ROOT)
    import msld_prm

    monkeypatch.chdir(CHEMBL_DIR)
    msld_prm.MsldPRM(written_outdir, cgenff=True, verbose=False, debug=False)
    counts = _prm_section_entries(os.path.join(written_outdir, "full_ligand.prm"))
    assert counts.get("BONDS", 0) and counts.get("ANGLES", 0) and counts.get("DIHEDRALS", 0)
