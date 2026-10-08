import os
import shutil

import pytest

from msldprep import align, analog_prep
from msldprep.errors import AlignmentError

from .conftest import CHEMBL_DIR, REFERENCE_STEM

# The CGenFF program is found through MSLDPREP_CGENFF_PATH, else on PATH. Its base parameter
# file is MSLDPREP_CGENFF_BASE_PRM, else par_all36_cgenff.prm in the program's directory.
CGENFF = os.environ.get("MSLDPREP_CGENFF_PATH") or shutil.which("cgenff") or ""
CGENFF_BASE_PRM = os.environ.get("MSLDPREP_CGENFF_BASE_PRM") or os.path.join(
    os.path.dirname(CGENFF), "par_all36_cgenff.prm"
)

TOOLS_AVAILABLE = os.path.isfile(CGENFF) and os.path.isfile(CGENFF_BASE_PRM)
requires_tools = pytest.mark.skipif(
    not TOOLS_AVAILABLE,
    reason="cgenff not found: set MSLDPREP_CGENFF_PATH (and MSLDPREP_CGENFF_BASE_PRM if the "
    "base prm is not next to the program)",
)


@pytest.fixture(scope="module")
def tools():
    return analog_prep.ExternalTools(cgenff=CGENFF, cgenff_base_prm=CGENFF_BASE_PRM)


@pytest.fixture(scope="module")
def fresh_reference():
    return align.load_ligand(
        os.path.join(CHEMBL_DIR, f"{REFERENCE_STEM}.mol2"),
        os.path.join(CHEMBL_DIR, f"{REFERENCE_STEM}.rtf"),
        os.path.join(CHEMBL_DIR, f"{REFERENCE_STEM}.prm"),
        name=REFERENCE_STEM, is_reference=True,
    )


def test_generate_3d_analog_bad_smiles_raises(fresh_reference):
    with pytest.raises(AlignmentError):
        analog_prep.generate_3d_analog("not a smiles!!", fresh_reference, name="bad")


@requires_tools
class TestRealTools:
    def test_mol_to_ligand_record_touches_no_disk(self, tmp_path, fresh_reference, tools, monkeypatch):
        """The streamed builder must never write .mol2/.rtf/.prm files --
        cgenff's default mode is piped via stdin/stdout only."""
        monkeypatch.chdir(tmp_path)
        mol = analog_prep.generate_3d_analog("Fc1ccc(Nc2ncnc3ccccc23)cc1", fresh_reference, name="stream_f")
        before = set(os.listdir(str(tmp_path)))
        ligand = analog_prep.mol_to_ligand_record(mol, "stream_f", tools.cgenff, cgenff_base_prm=tools.cgenff_base_prm)
        after = set(os.listdir(str(tmp_path)))

        assert after == before
        assert ligand.mol2_path is None
        assert ligand.rtf_path is None
        assert ligand.prm_path is None
        assert ligand.atom_types
        assert ligand.prm.bonds or ligand.prm.angles or ligand.prm.dihedrals or ligand.prm.impropers
        for atom in ligand.mol.GetAtoms():
            assert atom.HasProp("CGENFFAtomType")
            assert atom.HasProp("CGENFFPartialCharge")

    def test_build_analog_from_smiles_end_to_end(self, tmp_path, fresh_reference, tools):
        ligand = analog_prep.build_analog_from_smiles(
            "Fc1ccc(Nc2ncnc3ccccc23)cc1", "test_f_analog", fresh_reference, str(tmp_path), tools=tools
        )
        assert ligand.name == "test_f_analog"
        assert ligand.mol.GetNumAtoms() > 0
        for atom in ligand.mol.GetAtoms():
            assert atom.HasProp("CGENFFAtomType")
            assert atom.HasProp("CGENFFPartialCharge")
        for ext in ("mol2", "rtf", "prm"):
            assert os.path.exists(os.path.join(str(tmp_path), f"test_f_analog.{ext}")), ext

    def test_build_analog_atom_names_match_between_mol2_and_rtf(self, tmp_path, fresh_reference, tools):
        """Regression guard for the Cl/Br case-sensitivity bug: every atom
        name in the generated mol2 must exist verbatim in the generated rtf."""
        ligand = analog_prep.build_analog_from_smiles(
            "Clc1ccc(Nc2ncnc3ccccc23)cc1", "test_cl_analog", fresh_reference, str(tmp_path), tools=tools
        )
        rtf_names = set()
        with open(os.path.join(str(tmp_path), "test_cl_analog.rtf")) as fp:
            for line in fp:
                if line.startswith("ATOM"):
                    rtf_names.add(line.split()[1])

        for atom in ligand.mol.GetAtoms():
            assert atom.GetProp("_TriposAtomName") in rtf_names

    def test_build_analogs_from_smiles_collects_failures(self, tmp_path, fresh_reference, tools):
        result = analog_prep.build_analogs_from_smiles(
            {"good": "Fc1ccc(Nc2ncnc3ccccc23)cc1", "bad": "not a smiles!!"},
            fresh_reference, str(tmp_path), tools=tools,
        )
        assert len(result.ligands) == 1
        assert result.ligands[0].name == "good"
        assert len(result.failed) == 1
        assert result.failed[0][0] == "bad"

    def test_run_pipeline_from_smiles_relative_workdir(self, tmp_path, fresh_reference, tools, monkeypatch):
        """Regression guard for a relative-path bug once present in the
        (now-removed) obabel conversion step: run from a cwd where
        workdir/outdir are given as relative paths."""
        from msldprep import pipeline

        monkeypatch.chdir(tmp_path)
        result = pipeline.run_pipeline_from_smiles(
            fresh_reference, {"f_analog": "Fc1ccc(Nc2ncnc3ccccc23)cc1"},
            sysname="relpath_test", outdir="build_rel", workdir="work_rel",
            cgenff=True, tools=tools, toppar_search_dir=CHEMBL_DIR,
            # a single small analog's whole-molecule similarity to the
            # reference is beside the point for a path-handling regression
            # test -- disable grouping so it isn't at the mercy of that
            # unrelated heuristic.
            similarity_threshold=0,
        )
        assert result is not None
        assert "f_analog" in {lig.name for lig in result.ligands}
        assert os.path.exists("build_rel/full_ligand.prm")
