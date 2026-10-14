import os

import pytest
from rdkit import Chem

from msldprep.gui import service
from msldprep.gui.service import AlignmentSession, SessionError, parse_analog_table

from .conftest import CHEMBL_DIR, REFERENCE_STEM

# Analogs of the public ChEMBL example series (examples/2_ChEMBL): SMILES of three
# ligands in that directory, written out so the tests do not depend on a reader.
ANALOG_SMILES = {
    "10749": "Oc1ccc2c(Nc3cccc(Br)c3)ncnc2c1",
    "10322": "COc1cc2ncnc(NNc3ccc(F)c(Cl)c3)c2cc1OC",
    "10870": "Clc1ccc2c(NCc3ccc4c(c3)OCO4)ncnc2c1",
}
ANALOG_TEXT = "\n".join(f"{smiles} {name}" for name, smiles in ANALOG_SMILES.items())


def _reference_text() -> str:
    with open(os.path.join(CHEMBL_DIR, f"{REFERENCE_STEM}.mol2")) as handle:
        return handle.read()


@pytest.fixture
def session():
    s = AlignmentSession()
    s.set_reference(REFERENCE_STEM, _reference_text())
    return s


@pytest.fixture
def aligned(session):
    session.set_analogs(ANALOG_TEXT + "\nCCO ethanol\nnonsense bad\n")
    session.align_all()
    return session


def _status(session):
    return {a["name"]: a["status"] for a in session.summary()["analogs"]}


# ---------------------------------------------------------------- the table


def test_parse_analog_table_accepts_either_column_order_and_separators():
    text = "# a comment\n\nCCO ethanol\nmethanol,CO\nCCC\tpropane\n"
    rows = parse_analog_table(text)
    assert [(r.name, r.smiles) for r in rows] == [("ethanol", "CCO"), ("methanol", "CO"), ("propane", "CCC")]
    assert all(r.error is None for r in rows)


def test_parse_analog_table_skips_a_smiles_header_and_names_bare_smiles():
    rows = parse_analog_table("smiles,name\nCCO\n")
    assert [(r.name, r.smiles) for r in rows] == [("analog2", "CCO")]


def test_parse_analog_table_keeps_a_bad_smiles_with_its_error_and_deduplicates_names():
    rows = parse_analog_table("CCO a\nCO a\nnonsense b\n")
    assert [r.name for r in rows] == ["a", "a_2", "b"]
    assert rows[2].error and "nonsense" in rows[2].error


def test_parse_analog_table_makes_names_safe_for_urls():
    assert parse_analog_table("CCO my/ligand 1")[0].name == "my_ligand"


# ------------------------------------------------------------- the session


def test_align_all_aligns_what_it_can_and_reports_the_rest(aligned):
    status = _status(aligned)
    assert [status[n] for n in ANALOG_SMILES] == ["aligned"] * 3
    assert status["ethanol"] == "failed" and status["bad"] == "failed"
    errors = {a["name"]: a["error"] for a in aligned.summary()["analogs"]}
    assert "no shared substructure" in errors["ethanol"]
    assert "nonsense" in errors["bad"]
    assert aligned.summary()["job"] == {"running": False, "done": 5, "total": 5}


def test_aligned_analogs_report_rmsd_core_and_seed(aligned):
    for a in aligned.summary()["analogs"][:3]:
        assert 0 <= a["rmsd"] < 0.1
        assert a["core_atoms"] > 0
        assert a["seed"] == 2342


def test_ligand_view_core_atoms_point_at_matching_atoms_in_the_blocks(aligned):
    view = aligned.ligand_view("10749")
    analog = Chem.MolFromMolBlock(view["molblock"], removeHs=False)
    reference = Chem.MolFromMolBlock(aligned.reference_view()["molblock"], removeHs=False)
    assert len(view["core_atoms"]) == len(view["reference_core_atoms"]) > 0
    # the two lists are each sorted, so they are not pairs; as sets of atoms they agree in elements
    assert sorted(analog.GetAtomWithIdx(i).GetSymbol() for i in view["core_atoms"]) == sorted(
        reference.GetAtomWithIdx(j).GetSymbol() for j in view["reference_core_atoms"]
    )
    assert view["svg"].startswith("<?xml") or "<svg" in view["svg"]


def test_ligand_view_of_a_failed_analog_has_the_error_and_no_structure(aligned):
    view = aligned.ligand_view("ethanol")
    assert view["status"] == "failed" and view["error"]
    assert "molblock" not in view


def test_unknown_analog_is_a_session_error(aligned):
    with pytest.raises(SessionError, match="No analog named"):
        aligned.ligand_view("nope")


def test_realign_with_the_same_seed_reproduces_the_pose_and_a_new_seed_changes_it(aligned):
    first = aligned.ligand_view("10749")["molblock"]
    aligned.realign("10749", 2342)
    assert aligned.ligand_view("10749")["molblock"] == first
    aligned.realign("10749", 7)
    view = aligned.ligand_view("10749")
    assert view["seed"] == 7 and view["molblock"] != first


def test_aligning_again_clears_the_decision_because_the_pose_is_new(aligned):
    aligned.set_decision("10749", "accepted")
    aligned.realign("10749", 3)
    assert aligned.ligand_view("10749")["decision"] == "undecided"


def test_decisions_apply_only_to_aligned_analogs(aligned):
    aligned.set_decision("10749", "accepted")
    aligned.set_decision("10322", "rejected")
    with pytest.raises(SessionError, match="not aligned"):
        aligned.set_decision("ethanol", "accepted")
    with pytest.raises(SessionError, match="Decision must be"):
        aligned.set_decision("10749", "maybe")
    assert [e["name"] for e in aligned.overlay()] == ["10749", "10870"]


def test_export_sdf_lists_reference_then_the_chosen_analogs(aligned):
    aligned.set_decision("10749", "accepted")
    aligned.set_decision("10322", "rejected")

    def names(which):
        supplier = Chem.SDMolSupplier()
        supplier.SetData(aligned.export_sdf(which), removeHs=False)
        return [(m.GetProp("_Name"), m.GetProp("Role")) for m in supplier]

    assert names("accepted") == [(REFERENCE_STEM, "reference"), ("10749", "analog")]
    assert [n for n, _ in names("aligned")] == [REFERENCE_STEM, "10749", "10870"]
    with pytest.raises(SessionError):
        aligned.export_sdf("everything")


def test_exported_analogs_keep_their_aligned_coordinates(aligned):
    supplier = Chem.SDMolSupplier()
    supplier.SetData(aligned.export_sdf("aligned"), removeHs=False)
    exported = {m.GetProp("_Name"): m for m in supplier}
    view = aligned.ligand_view("10870")
    expected = Chem.MolFromMolBlock(view["molblock"], removeHs=False).GetConformer().GetPositions()
    assert (exported["10870"].GetConformer().GetPositions() - expected).max() < 1e-3
    assert float(exported["10870"].GetProp("AlignRMSD")) == pytest.approx(view["rmsd"], abs=1e-4)


def test_changing_the_reference_clears_alignments_but_keeps_unreadable_smiles_failed(aligned):
    aligned.set_reference("again", _reference_text())
    status = _status(aligned)
    assert status["10749"] == "pending" and status["ethanol"] == "pending"
    assert status["bad"] == "failed"


def test_params_reach_the_alignment_and_are_validated(session):
    session.set_params(ring_matches_ring_only=True, timeout=30, seed=11)
    params = session.summary()["params"]
    assert params["ring_matches_ring_only"] is True and params["timeout"] == 30 and params["seed"] == 11
    session.set_analogs("Oc1ccc2c(Nc3cccc(Br)c3)ncnc2c1 a")
    session.align_all()
    assert session.summary()["analogs"][0]["seed"] == 11
    with pytest.raises(SessionError):
        session.set_params(timeout=-1)


def test_alignment_needs_a_reference_and_analogs():
    empty = AlignmentSession()
    with pytest.raises(SessionError, match="reference"):
        empty.align_all()
    empty.set_reference("r", _reference_text())
    with pytest.raises(SessionError, match="analogs"):
        empty.align_all()


def test_a_running_job_blocks_changes(session):
    session.set_analogs(ANALOG_TEXT)
    run = session.begin_alignment()
    assert session.summary()["job"]["running"] is True
    with pytest.raises(SessionError, match="running"):
        session.set_analogs("CCO a")
    with pytest.raises(SessionError, match="running"):
        session.begin_alignment()
    run()
    assert session.summary()["job"]["running"] is False


def test_a_reference_without_3d_coordinates_is_refused():
    flat = Chem.MolFromSmiles("CCO")
    Chem.rdDepictor.Compute2DCoords(flat)
    with pytest.raises(service.AlignmentError, match="3D"):
        AlignmentSession().set_reference("flat", Chem.MolToMolBlock(flat))


def test_reset_forgets_everything(aligned):
    aligned.reset()
    state = aligned.summary()
    assert state["reference"] is None and state["analogs"] == []


def test_importing_the_package_does_not_import_the_gui():
    import subprocess
    import sys

    code = "import sys, msldprep, msldprep.align, msldprep.pipeline; sys.exit(any(m.startswith('msldprep.gui') for m in sys.modules))"
    assert subprocess.run([sys.executable, "-B", "-c", code]).returncode == 0


# ------------------------------------------------- charges on the reference

# The reference of the ChEMBL series written with a protonated ring nitrogen (the
# analogs stay neutral), as a 3D mol block. Without removing the charge the alignment fails;
# see docs/dev/ISSUES_AND_PLANS.md, #9, and the xfail test in tests/test_align.py.
CHARGED_REFERENCE_SMILES = "c1ccc2c(c1)c([nH+]cn2)Nc1cccc(Br)c1"


def _charged_reference_block() -> str:
    mol = Chem.AddHs(Chem.MolFromSmiles(CHARGED_REFERENCE_SMILES))
    Chem.rdDistGeom.EmbedMolecule(mol, randomSeed=1)
    return Chem.MolToMolBlock(mol)


def test_neutralize_removes_charges_and_keeps_the_molecule_3d():
    charged = Chem.MolFromMolBlock(_charged_reference_block(), removeHs=False)
    neutral = service.neutralize(charged)
    assert sum(a.GetFormalCharge() for a in neutral.GetAtoms()) == 0
    assert neutral.GetNumAtoms() == charged.GetNumAtoms() - 1  # the nitrogen gave up its hydrogen
    assert neutral.GetConformer().Is3D()


def test_reference_summary_reports_charge_before_and_after():
    session = AlignmentSession()
    session.set_reference("charged", _charged_reference_block())
    ref = session.summary()["reference"]
    assert (ref["net_charge"], ref["loaded_net_charge"], ref["charged_atoms"]) == (1, 1, 1)
    session.set_neutralize_reference(True)
    ref = session.summary()["reference"]
    assert (ref["net_charge"], ref["loaded_net_charge"], ref["charged_atoms"]) == (0, 1, 0)
    assert session.summary()["neutralize_reference"] is True
    session.set_neutralize_reference(False)
    assert session.summary()["reference"]["net_charge"] == 1


def test_neutralizing_the_reference_lets_neutral_analogs_align():
    session = AlignmentSession()
    session.set_neutralize_reference(True)  # may be set before the reference is loaded
    session.set_reference("charged", _charged_reference_block())
    session.set_analogs(ANALOG_TEXT)
    session.align_all()
    assert set(_status(session).values()) == {"aligned"}


def test_toggling_neutralize_clears_alignments(aligned):
    aligned.set_neutralize_reference(True)
    assert _status(aligned)["10749"] == "pending"


def test_reference_view_lists_its_elements(session):
    elements = session.reference_view()["elements"]
    assert elements == sorted(set(elements)) and {"C", "H", "N"} <= set(elements)


def test_a_failure_caused_by_a_charged_reference_says_so():
    session = AlignmentSession()
    session.set_reference("charged", _charged_reference_block())
    session.set_analogs(ANALOG_TEXT)
    session.align_all()
    errors = [a["error"] for a in session.summary()["analogs"]]
    assert all("try removing its charges" in e for e in errors)
