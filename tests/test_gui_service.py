import io
import os
import zipfile

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


def _sdf_records(text):
    supplier = Chem.SDMolSupplier()
    supplier.SetData(text, removeHs=False)
    return {m.GetProp("_Name"): m for m in supplier}


def _mol2_records(text):
    blocks = ["@<TRIPOS>MOLECULE" + part for part in text.split("@<TRIPOS>MOLECULE")[1:]]
    mols = [Chem.MolFromMol2Block(b, removeHs=False) for b in blocks]
    assert all(m is not None for m in mols)
    return {b.splitlines()[1]: m for b, m in zip(blocks, mols)}


def test_export_sdf_lists_reference_then_the_chosen_analogs(aligned):
    aligned.set_decision("10749", "accepted")
    aligned.set_decision("10322", "rejected")

    def names(which):
        records = _sdf_records(aligned.export("sdf", which).data.decode())
        return [(name, m.GetProp("Role")) for name, m in records.items()]

    assert names("accepted") == [(REFERENCE_STEM, "reference"), ("10749", "analog")]
    assert [n for n, _ in names("aligned")] == [REFERENCE_STEM, "10749", "10870"]
    with pytest.raises(SessionError, match="must be"):
        aligned.export("sdf", "everything")
    with pytest.raises(SessionError, match="format"):
        aligned.export("pdb")


def test_exported_analogs_keep_their_aligned_coordinates(aligned):
    exported = _sdf_records(aligned.export("sdf", "aligned").data.decode())
    view = aligned.ligand_view("10870")
    expected = Chem.MolFromMolBlock(view["molblock"], removeHs=False).GetConformer().GetPositions()
    assert (exported["10870"].GetConformer().GetPositions() - expected).max() < 1e-3
    assert float(exported["10870"].GetProp("AlignRMSD")) == pytest.approx(view["rmsd"], abs=1e-4)


def test_export_mol2_writes_sybyl_typed_structures_in_the_aligned_pose(aligned):
    export = aligned.export("mol2", "aligned")
    assert export.filename == "aligned_aligned.mol2" and export.media_type == "chemical/x-mol2"
    text = export.data.decode()
    assert text.count("@<TRIPOS>MOLECULE") == 4  # reference and three analogs
    records = _mol2_records(text)
    assert list(records) == [REFERENCE_STEM, "10749", "10322", "10870"]
    # the notebook's writer: Sybyl atom types and unique atom names, no charges
    assert "NO_CHARGES" in text and "C.ar" in text
    names = [a.GetProp("_TriposAtomName") for a in records["10749"].GetAtoms()]
    assert len(set(names)) == len(names)
    view = aligned.ligand_view("10749")
    expected = Chem.MolFromMolBlock(view["molblock"], removeHs=False).GetConformer().GetPositions()
    assert (records["10749"].GetConformer().GetPositions() - expected).max() < 1e-3
    # exporting does not change what the session holds
    assert aligned.ligand_view("10749")["molblock"] == view["molblock"]


def test_export_one_structure_per_file_makes_a_zip_named_after_the_structures(aligned):
    for fmt, read in (("sdf", _sdf_records), ("mol2", _mol2_records)):
        export = aligned.export(fmt, "aligned", one_per_file=True)
        assert export.filename == f"aligned_aligned_{fmt}.zip" and export.media_type == "application/zip"
        with zipfile.ZipFile(io.BytesIO(export.data)) as zf:
            assert zf.namelist() == [f"{n}.{fmt}" for n in (REFERENCE_STEM, "10749", "10322", "10870")]
            for member in zf.namelist():
                records = read(zf.read(member).decode())
                assert list(records) == [member.rsplit(".", 1)[0]]  # exactly one structure per file


def test_export_can_leave_out_the_reference(aligned):
    records = _sdf_records(aligned.export("sdf", "aligned", include_reference=False).data.decode())
    assert list(records) == ["10749", "10322", "10870"]
    with zipfile.ZipFile(io.BytesIO(aligned.export("mol2", "aligned", include_reference=False, one_per_file=True).data)) as zf:
        assert REFERENCE_STEM + ".mol2" not in zf.namelist()


def test_a_reference_named_like_an_analog_gets_its_own_file_in_the_zip():
    session = AlignmentSession()
    session.set_reference("10749", _reference_text())
    session.set_analogs(ANALOG_TEXT)
    session.align_all()
    with zipfile.ZipFile(io.BytesIO(session.export("sdf", "aligned", one_per_file=True).data)) as zf:
        assert zf.namelist() == ["10749_reference.sdf", "10749.sdf", "10322.sdf", "10870.sdf"]


def test_export_refuses_when_there_is_nothing_to_write(aligned):
    with pytest.raises(SessionError, match="no analog is marked accepted"):
        aligned.export("sdf", "accepted")
    for name in ANALOG_SMILES:
        aligned.set_decision(name, "rejected")
    with pytest.raises(SessionError, match="no analog is aligned"):
        aligned.export("mol2", "aligned")
    with pytest.raises(SessionError, match="No reference"):
        AlignmentSession().export()


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
    for bad in (-1, 0):
        with pytest.raises(SessionError, match="at least 1 second"):
            session.set_params(timeout=bad)
    assert AlignmentSession().summary()["params"]["timeout"] == 10


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
# analogs stay neutral), as a 3D mol block (#9: a charged reference used to break the match).
CHARGED_REFERENCE_SMILES = "c1ccc2c(c1)c([nH+]cn2)Nc1cccc(Br)c1"


def _charged_reference_block() -> str:
    mol = Chem.AddHs(Chem.MolFromSmiles(CHARGED_REFERENCE_SMILES))
    Chem.rdDistGeom.EmbedMolecule(mol, randomSeed=1)
    return Chem.MolToMolBlock(mol)


def test_the_reference_summary_shows_its_net_charge_and_charged_atoms():
    session = AlignmentSession()
    session.set_reference("charged", _charged_reference_block())
    ref = session.summary()["reference"]
    assert (ref["net_charge"], ref["charged_atoms"]) == (1, 1)


def test_a_charged_reference_aligns_without_removing_its_charges():
    # #9: the charge is no longer part of what an analog has to match
    session = AlignmentSession()
    session.set_reference("charged", _charged_reference_block())
    session.set_analogs(ANALOG_TEXT)
    session.align_all()
    assert set(_status(session).values()) == {"aligned"}
    assert session.summary()["reference"]["net_charge"] == 1


def test_the_gui_ignores_bond_orders_by_default_and_the_setting_can_be_changed(session):
    assert session.summary()["params"]["ignore_bond_orders"] is True
    session.set_params(ignore_bond_orders=False)
    assert session.summary()["params"]["ignore_bond_orders"] is False


def test_tautomers_of_the_reference_align_in_full_with_the_default_settings():
    reference = Chem.AddHs(Chem.MolFromSmiles("O=c1[nH]cnc2ccccc12"))
    Chem.rdDistGeom.EmbedMolecule(reference, randomSeed=3)
    session = AlignmentSession()
    session.set_reference("quinazolinone", Chem.MolToMolBlock(reference))
    session.set_analogs("Oc1ncnc2ccccc12 hydroxy_tautomer")
    session.align_all()
    assert session.summary()["analogs"][0]["core_atoms"] == 11  # all heavy atoms, the C-O bond order differs
    session.set_params(ignore_bond_orders=False)
    session.align_all()
    assert session.summary()["analogs"][0]["core_atoms"] == 10


def test_protonation_changes_are_flagged_only_when_asked_and_the_pose_stays():
    reference = Chem.AddHs(Chem.MolFromSmiles("O=c1cccc[nH]1"))
    Chem.rdDistGeom.EmbedMolecule(reference, randomSeed=3)
    session = AlignmentSession()
    session.set_reference("pyridone", Chem.MolToMolBlock(reference))
    session.set_analogs("Oc1ccccn1 hydroxypyridine")
    assert session.summary()["params"]["check_protonation"] is False
    session.align_all()
    tolerant = session.ligand_view("hydroxypyridine")
    assert tolerant["protonation"] == [] and session.summary()["analogs"][0]["protonation_flag"] is False

    session.set_params(check_protonation=True)
    session.align_all()
    strict = session.ligand_view("hydroxypyridine")
    assert session.summary()["analogs"][0]["protonation_flag"] is True
    assert strict["molblock"] == tolerant["molblock"]
    assert sorted(note["atom"][0] for note in strict["protonation"]) == ["N", "O"]
    assert any("neutral, 0 H; reference N" in note["note"] for note in strict["protonation"])


def test_the_view_says_how_many_orientations_the_core_had_and_which_was_used():
    session = AlignmentSession()
    reference = Chem.AddHs(Chem.MolFromSmiles("COc1ccccc1CC(=O)O"))
    Chem.rdDistGeom.EmbedMolecule(reference, randomSeed=5)
    session.set_reference("ortho", Chem.MolToMolBlock(reference))
    session.set_analogs("Clc1ccccc1CC(=O)O chloro\nc1ccccc1CC(=O)O plain")
    session.align_all()
    chloro, plain = session.ligand_view("chloro"), session.ligand_view("plain")
    assert chloro["orientations"] == 2 and chloro["orientation"] in (0, 1)
    assert (plain["orientations"], plain["orientation"]) == (1, 0)


def test_the_page_is_told_about_stereo_that_needs_a_look():
    session = AlignmentSession()
    reference = Chem.AddHs(Chem.MolFromSmiles("C[C@H](c1ccc(Cl)cc1)C(=O)NC"))
    Chem.rdDistGeom.EmbedMolecule(reference, randomSeed=7)
    session.set_reference("chiral", Chem.MolToMolBlock(reference))
    assert session.summary()["reference"]["stereocenters"] == ["R"]
    session.set_analogs(
        "C[C@H](c1ccc(F)cc1)C(=O)NCC same\nC[C@@H](c1ccc(F)cc1)C(=O)NCC opposite\nCC(c1ccc(F)cc1)C(=O)NCC unspecified"
    )
    session.align_all()
    flags = {a["name"]: a["stereo_flag"] for a in session.summary()["analogs"]}
    assert flags == {"same": False, "opposite": True, "unspecified": False}
    opposite = session.ligand_view("opposite")["stereo"][0]
    assert opposite["flag"] and "opposite" in opposite["note"]
    assert "match the reference" in session.ligand_view("unspecified")["stereo"][0]["note"]


def test_a_reference_without_stereocenters_says_so(session):
    assert session.summary()["reference"]["stereocenters"] == []
