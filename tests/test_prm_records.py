import os

import pytest
from crimm.IO.PRMParser import angle_par, bond_par, nonbond_par

from msldprep import prm_records
from msldprep.errors import PrmError

from .conftest import CHEMBL_DIR, REFERENCE_STEM


def test_parse_prm_file_real_data():
    path = os.path.join(CHEMBL_DIR, f"{REFERENCE_STEM}.prm")
    parsed = prm_records.parse_prm_file(path)
    assert parsed.bonds
    assert parsed.angles
    assert parsed.dihedrals
    assert parsed.impropers
    # spot-check a known entry (also checked manually against the raw file)
    assert parsed.bonds.get(("CG2R64", "NG311")) == bond_par(390.0, 1.355) or parsed.bonds.get(
        ("NG311", "CG2R64")
    ) == bond_par(390.0, 1.355)


def test_parse_prm_file_missing_raises():
    with pytest.raises(PrmError):
        prm_records.load_prm_files(["does_not_exist_12345"], CHEMBL_DIR)


def test_format_line_helpers_round_trip_field_order():
    assert prm_records.format_bond_line(("A", "B"), bond_par(100.0, 1.5)) == "A B 100.0000 1.5000"
    assert prm_records.format_angle_line(("A", "B", "C"), angle_par(40.0, 109.0), None) == "A B C 40.0000 109.0000"
    from crimm.IO.PRMParser import ub_par

    assert (
        prm_records.format_angle_line(("A", "B", "C"), angle_par(40.0, 109.0), ub_par(5.0, 2.0))
        == "A B C 40.0000 109.0000 5.0000 2.0000"
    )


def test_merge_prm_records_first_ligand_wins_on_conflict():
    a = prm_records.LigandParameters(bonds={("A", "B"): bond_par(100.0, 1.5)})
    b = prm_records.LigandParameters(bonds={("A", "B"): bond_par(200.0, 1.6)})
    merged = prm_records.merge_prm_records({"lig_a": a, "lig_b": b}, cgenff=True, ligand_order=["lig_a", "lig_b"])
    assert merged.bonds[("A", "B")] == bond_par(100.0, 1.5)
    assert merged.duplicates["BOND"] == ["!A B 200.0000 1.6000"]


def test_merge_prm_records_reverse_key_matches():
    a = prm_records.LigandParameters(bonds={("A", "B"): bond_par(100.0, 1.5)})
    b = prm_records.LigandParameters(bonds={("B", "A"): bond_par(100.0, 1.5)})  # same bond, reversed, not a conflict
    merged = prm_records.merge_prm_records({"lig_a": a, "lig_b": b}, cgenff=True, ligand_order=["lig_a", "lig_b"])
    assert merged.bonds == {("A", "B"): bond_par(100.0, 1.5)}
    assert merged.duplicates["BOND"] == []


def test_merge_prm_records_nonbonded_only_when_not_cgenff():
    a = prm_records.LigandParameters(nonbonded={"CG331": nonbond_par(-0.08, 2.04)})
    merged_cgenff = prm_records.merge_prm_records({"a": a}, cgenff=True, ligand_order=["a"])
    assert merged_cgenff.nonbonded == {}
    merged_explicit = prm_records.merge_prm_records({"a": a}, cgenff=False, ligand_order=["a"])
    assert merged_explicit.nonbonded == {"CG331": nonbond_par(-0.08, 2.04)}


def test_build_ligand_topology_and_derive_needed_parameters_real_data(all_ligands, decomposition):
    ligand_names = [lig.name for lig in all_ligands]
    prm_by_ligand = prm_records.load_prm_files(ligand_names, CHEMBL_DIR)
    merged = prm_records.merge_prm_records(prm_by_ligand, cgenff=True, ligand_order=ligand_names)
    topology = prm_records.build_ligand_topology(decomposition, all_ligands)

    assert topology.atoms
    assert topology.bonds
    # every topology atom key should be unique
    keys = [key for key, _ in topology.atoms]
    assert len(keys) == len(set(keys))

    needed = prm_records.derive_needed_parameters(topology, merged, cgenff=True)
    assert needed.bonds
    assert needed.angles
    assert needed.dihedrals
    # every derived bond/angle key must trace back to the merged library
    for key in needed.order["BOND"]:
        assert key in merged.bonds
    for key in needed.order["ANGL"]:
        assert key in merged.angles
    # every dihedral entry is the WHOLE per-atom-tuple term list, not one row per multiplicity
    for key, terms in needed.dihedrals.items():
        assert terms == merged.dihedrals[key]


def test_ab2re_zero_rmin_returns_zero():
    assert prm_records.ab2re(0.1, 0.0, 0.2, 1.5) == (0.0, 0.0)


def test_compute_nbfix_self_consistent():
    nonbonded = {
        "C1": nonbond_par(-0.10, 2.0),
        "C2": nonbond_par(-0.05, 1.8),
    }
    nbfix = prm_records.compute_nbfix(nonbonded)
    assert ("C1", "C1") in nbfix
    assert ("C1", "C2") in nbfix
    for (t1, t2), (eps, r, eps14, r14) in nbfix.items():
        assert eps <= 0  # CHARMM stores eps as non-positive
        assert r >= 0


def test_compute_nbfix_uses_1_4_terms_when_present():
    from crimm.IO.PRMParser import nonbond14_par

    nonbonded = {"C1": nonbond_par(-0.10, 2.0)}
    nonbonded14 = {"C1": nonbond14_par(-0.05, 1.9)}
    with_14 = prm_records.compute_nbfix(nonbonded, nonbonded14)
    without_14 = prm_records.compute_nbfix(nonbonded)
    assert with_14[("C1", "C1")][2:] != without_14[("C1", "C1")][2:]
