import os

import pytest

from msldprep import align, crn, decompose
from msldprep.writers.records import compute_atom_name_translation

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CHEMBL_DIR = os.path.join(REPO_ROOT, "examples", "2_ChEMBL")

# The original (pre-rewrite) implementation this package replaces -- a
# sibling directory, not a parent of this repo. Only used by tests that
# deliberately check output against its real, unmodified legacy scripts
# (e.g. test_writers.py's msld_prm.py backward-compat check).
LEGACY_REPO_ROOT = os.path.join(os.path.dirname(REPO_ROOT), "msld-py-prep")

REFERENCE_STEM = "11896"
ANALOG_STEMS = ["10749", "10322", "10870"]


def _load(stem: str, *, is_reference: bool = False):
    return align.load_ligand(
        os.path.join(CHEMBL_DIR, f"{stem}.mol2"),
        os.path.join(CHEMBL_DIR, f"{stem}.rtf"),
        os.path.join(CHEMBL_DIR, f"{stem}.prm"),
        name=stem,
        is_reference=is_reference,
    )


@pytest.fixture(scope="session")
def reference_ligand():
    return _load(REFERENCE_STEM, is_reference=True)


@pytest.fixture(scope="session")
def prepared_ligand_set(reference_ligand):
    analogs = [_load(stem) for stem in ANALOG_STEMS]
    return align.prepare_ligand_set(reference_ligand, analogs)


@pytest.fixture(scope="session")
def all_ligands(reference_ligand, prepared_ligand_set):
    return [reference_ligand, *prepared_ligand_set.aligned]


@pytest.fixture(scope="session")
def raw_decomposition(all_ligands, prepared_ligand_set):
    return decompose.decompose_ligand_set(all_ligands, prepared_ligand_set.core_query_mol, REFERENCE_STEM)


@pytest.fixture(scope="session")
def decomposition(raw_decomposition, all_ligands):
    return decompose.apply_manual_overrides(raw_decomposition, all_ligands)


@pytest.fixture(scope="session")
def crn_result(all_ligands, decomposition):
    return crn.compute_crn(all_ligands, decomposition)


@pytest.fixture(scope="session")
def translation(decomposition, all_ligands):
    return compute_atom_name_translation(decomposition, all_ligands)
