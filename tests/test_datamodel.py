import os

import pytest

from msldprep.datamodel import LigandRecord
from msldprep.errors import AlignmentError
from msldprep.prm_records import parse_prm_text
from msldprep.rtf_records import parse_rtf_text

from .conftest import CHEMBL_DIR, REFERENCE_STEM

_RTF_FIELDS = [
    "net_charge", "format_version", "atom_types", "atom_charges", "bonds", "impropers", "lone_pairs", "mass_records",
]


def _load_reference():
    return LigandRecord.from_files(
        os.path.join(CHEMBL_DIR, f"{REFERENCE_STEM}.mol2"),
        os.path.join(CHEMBL_DIR, f"{REFERENCE_STEM}.rtf"),
        os.path.join(CHEMBL_DIR, f"{REFERENCE_STEM}.prm"),
        name=REFERENCE_STEM,
        is_reference=True,
    )


@pytest.mark.xfail(
    strict=True,
    raises=AlignmentError,
    reason="RDKit's mol2 reader rejects a nitro O.co2 (11408) and cannot kekulize a ring NH typed N.ar (15723); "
    "20 of the 89 public ChEMBL files fail (#19)",
)
@pytest.mark.parametrize("stem", ["11408", "15723"])
def test_from_files_reads_every_public_example_form(stem):
    ligand = LigandRecord.from_files(
        os.path.join(CHEMBL_DIR, f"{stem}.mol2"),
        os.path.join(CHEMBL_DIR, f"{stem}.rtf"),
        os.path.join(CHEMBL_DIR, f"{stem}.prm"),
        name=stem,
    )
    assert ligand.mol.GetNumAtoms() > 0


def test_from_files_populates_rtf_and_prm_fields_directly():
    ligand = _load_reference()
    assert ligand.name == REFERENCE_STEM
    assert ligand.is_reference is True
    assert ligand.mol.GetNumAtoms() > 0
    for atom in ligand.mol.GetAtoms():
        assert atom.HasProp("CGENFFAtomType")
        assert atom.HasProp("CGENFFPartialCharge")

    # a real IMPR/LONEPAIR-bearing ligand -- confirms these fields actually
    # got populated, not just left at their empty defaults.
    assert ligand.impropers
    assert ligand.lone_pairs
    assert ligand.atom_types
    assert ligand.prm.bonds
    assert ligand.prm.dihedrals

    assert ligand.mol2_path.endswith(f"{REFERENCE_STEM}.mol2")
    assert ligand.rtf_path.endswith(f"{REFERENCE_STEM}.rtf")
    assert ligand.prm_path.endswith(f"{REFERENCE_STEM}.prm")


def test_to_rtf_text_round_trips_through_the_real_parser():
    ligand = _load_reference()
    reparsed = parse_rtf_text(ligand.to_rtf_text())
    for field in _RTF_FIELDS:
        assert reparsed[field] == getattr(ligand, field), field


def test_write_rtf_writes_a_file_that_round_trips(tmp_path):
    ligand = _load_reference()
    path = str(tmp_path / "roundtrip.rtf")
    ligand.write_rtf(path)
    assert os.path.exists(path)
    with open(path) as fp:
        reparsed = parse_rtf_text(fp.read())
    for field in _RTF_FIELDS:
        assert reparsed[field] == getattr(ligand, field), field


def test_to_prm_text_round_trips_through_the_real_parser():
    ligand = _load_reference()
    reparsed = parse_prm_text(ligand.to_prm_text())
    assert reparsed == ligand.prm


def test_write_prm_writes_a_file_that_round_trips(tmp_path):
    ligand = _load_reference()
    path = str(tmp_path / "roundtrip.prm")
    ligand.write_prm(path)
    assert os.path.exists(path)
    with open(path) as fp:
        reparsed = parse_prm_text(fp.read())
    assert reparsed == ligand.prm
