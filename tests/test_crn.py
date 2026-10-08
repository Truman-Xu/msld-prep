import pytest

from msldprep import crn
from msldprep.datamodel import CoreDecomposition, Fragment, Site
from msldprep.errors import CRNError


def test_redistribute_charge_proportional_hits_target():
    charges = {"C1": -0.1, "C2": 0.05, "H1": 0.05}
    result = crn.redistribute_charge_proportional(["C1", "C2", "H1"], charges, target_total=0.0)
    assert round(sum(result.values()), 6) == 0.0


def test_redistribute_charge_proportional_skips_lone_pairs():
    charges = {"C1": 0.0, "LP1": 0.05}
    # with only one real atom, the whole correction must land on C1 and
    # LP1 must stay untouched.
    result = crn.redistribute_charge_proportional(["C1", "LP1"], charges, target_total=0.1)
    assert result["LP1"] == 0.05
    assert round(result["C1"], 6) == 0.05


def test_redistribute_charge_proportional_weights_by_magnitude():
    """An atom with a larger existing |charge| should absorb proportionally
    more of the correction than one close to zero -- e.g. C1 (-0.30) should
    take 3x the correction that C2 (-0.10) does."""
    charges = {"C1": -0.30, "C2": -0.10}
    result = crn.redistribute_charge_proportional(["C1", "C2"], charges, target_total=-0.44)
    delta_c1 = result["C1"] - charges["C1"]
    delta_c2 = result["C2"] - charges["C2"]
    assert round(sum(result.values()), 6) == -0.44
    assert delta_c1 < 0 and delta_c2 < 0  # both moved further negative to close the gap
    assert round(delta_c1 / delta_c2, 2) == 3.0


def test_redistribute_charge_proportional_splits_evenly_when_all_zero():
    charges = {"C1": 0.0, "C2": 0.0}
    result = crn.redistribute_charge_proportional(["C1", "C2"], charges, target_total=0.02)
    assert round(result["C1"], 6) == round(result["C2"], 6) == 0.01


def test_redistribute_charge_proportional_raises_when_no_real_atoms():
    with pytest.raises(CRNError):
        crn.redistribute_charge_proportional(["LP1"], {"LP1": 0.05}, target_total=0.1)


def test_redistribute_does_not_mutate_input():
    charges = {"C1": 0.0}
    original = dict(charges)
    crn.redistribute_charge_proportional(["C1"], charges, target_total=0.5)
    assert charges == original


def test_compute_crn_preserves_each_ligands_net_charge(all_ligands, decomposition, crn_result):
    """The core invariant CRN exists to guarantee: after redistribution, each
    ligand's (shared core + its own fragments at every site) charges still
    sum to that ligand's original net charge."""
    for ligand in all_ligands:
        total = sum(crn_result.core_charges.values())
        for site_charges in crn_result.fragment_charges:
            total += sum(site_charges[ligand.name].values())
        assert round(total, 6) == round(ligand.net_charge, 6), (
            f"{ligand.name}: assembled CRN'd charge {total} != net_charge {ligand.net_charge}"
        )


def test_core_lone_pair_charge_included_in_core_sum():
    """Synthetic test for the core-hosted-LP path: no real example in the
    ChEMBL fixture set has an LP on a shared core atom (all sample LONEPAIR
    records are on the varying halogen substituent), so this constructs one
    directly against the datamodel/compute_crn contract instead of going
    through align/decompose.
    """
    decomposition = CoreDecomposition(
        reference_name="A",
        core_query_mol=None,
        core_atom_order=("C1",),
        core_atom_names={"A": ("C1",), "B": ("C1",)},
        core_lp_atom_names={"C1": {"A": "LP1", "B": "LP1"}},
        sites=[Site(index=0, fragments={
            "A": Fragment(ligand_name="A", site_index=0, atom_names=["X1"]),
            "B": Fragment(ligand_name="B", site_index=0, atom_names=["X1"]),
        })],
    )
    ligands = [
        _fake_ligand("A", net_charge=0.0, atom_charges={"C1": -0.05, "LP1": 0.05, "X1": 0.0}),
        _fake_ligand("B", net_charge=0.0, atom_charges={"C1": -0.06, "LP1": 0.05, "X1": 0.01}),
    ]

    result = crn.compute_crn(ligands, decomposition)
    # LP1 must appear in the core charge stats/values, and must be left
    # exactly at its pre-CRN mean -- it's included in the core *sum* (this
    # is the whole point of the fix) but never perturbed by the
    # redistribution, which only ever touches non-lone-pair atoms.
    assert "LP1" in result.core_charges
    assert result.core_charges["LP1"] == 0.05

    # The redistributed C1+LP1 sum must equal ref_charge - site_sum (the
    # actual CRN target), not the raw pre-redistribution core mean --
    # site_sum here is the mean of X1's charges (0.0, 0.01) = 0.005, so the
    # core must net out to 0.0 - 0.005 = -0.005.
    site_avg = round((0.0 + 0.01) / 2, 6)
    expected_core_sum = round(0.0 - site_avg, 6)
    assert round(result.core_charges["C1"] + result.core_charges["LP1"], 6) == expected_core_sum


def _fake_ligand(name, net_charge, atom_charges=None):
    from rdkit import Chem

    from msldprep.datamodel import LigandRecord

    return LigandRecord(name=name, mol=Chem.Mol(), net_charge=net_charge, atom_charges=atom_charges or {})
