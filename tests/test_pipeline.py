import os

from msldprep import align, io_mcs_text, pipeline

from .conftest import ANALOG_STEMS, CHEMBL_DIR, REFERENCE_STEM

# This file tests pipeline plumbing (file generation, charge conservation,
# path handling), not grouping quality -- the 4-ligand ChEMBL fixture set
# shares a perfectly good MCS core but has a smaller shared core relative to
# its varying substituents than msldprep/grouping.py's default whole-molecule
# similarity threshold tolerates, so every test here disables grouping
# (similarity_threshold=0) to keep testing exactly what it always has:
# one global core across the whole set. Grouping itself is covered by
# tests/test_grouping.py.
_NO_GROUPING = {"similarity_threshold": 0}


def _fresh_ligands(base_dir=CHEMBL_DIR, names=None):
    """pipeline.define_common_core mutates ligand.mol during alignment, so
    each test that calls it needs its own freshly-loaded ligands, not the
    session-scoped fixtures other tests share."""
    names = names if names is not None else ANALOG_STEMS
    reference = align.load_ligand(
        os.path.join(base_dir, f"{REFERENCE_STEM}.mol2"),
        os.path.join(base_dir, f"{REFERENCE_STEM}.rtf"),
        os.path.join(base_dir, f"{REFERENCE_STEM}.prm"),
        name=REFERENCE_STEM, is_reference=True,
    )
    analogs = [
        align.load_ligand(
            os.path.join(base_dir, f"{stem}.mol2"),
            os.path.join(base_dir, f"{stem}.rtf"),
            os.path.join(base_dir, f"{stem}.prm"),
            name=stem,
        )
        for stem in names
    ]
    return reference, analogs


def test_run_pipeline_produces_expected_files(tmp_path):
    reference, analogs = _fresh_ligands()
    outdir = str(tmp_path / "build")
    result = pipeline.run_pipeline(
        reference, analogs, sysname="fxa", outdir=outdir, cgenff=True, toppar_search_dir=CHEMBL_DIR,
        **_NO_GROUPING,
    )
    assert result is not None
    for expected in ("core.rtf", "core.pdb", "full_ligand.prm", "fxa.inp", "config.json", "records.json"):
        assert os.path.exists(os.path.join(outdir, expected)), expected
    # run_pipeline's own artifacts from define_common_core, for later inspection/reuse
    assert os.path.exists(os.path.join(outdir, "MCS_for_MSLD.txt"))
    assert os.path.exists(pipeline.aligned_sdf_path(os.path.join(outdir, "aligned")))


def test_explicit_two_step_workflow_as_separate_process_invocations_would(tmp_path):
    """The replacement for the old save_checkpoint two-run trick: two
    genuinely separate steps, with step 2 reloading ligands from scratch --
    exactly what running as two separate process invocations would do --
    rather than reusing step 1's in-memory (already-aligned) LigandRecords.
    """
    work_dir = str(tmp_path / "work")
    checkpoint = os.path.join(work_dir, "MCS_for_MSLD.txt")
    aligned_dir = os.path.join(work_dir, "aligned")

    reference, analogs = _fresh_ligands()
    core = pipeline.define_common_core(
        reference, analogs, mcs_checkpoint_path=checkpoint, aligned_dir=aligned_dir,
        **_NO_GROUPING,
    )
    assert os.path.exists(checkpoint)
    assert {lig.name for lig in core.ligands} == {REFERENCE_STEM, *ANALOG_STEMS}
    assert os.path.exists(pipeline.aligned_sdf_path(aligned_dir))
    for name in [REFERENCE_STEM, *ANALOG_STEMS]:
        for ext in ("rtf", "prm"):
            assert os.path.exists(os.path.join(aligned_dir, f"{name}.{ext}"))

    outdir = str(tmp_path / "build")
    assert not os.path.exists(os.path.join(outdir, "full_ligand.prm"))

    # step 2, as a fresh process would: reload ligand names from the
    # checkpoint itself, then reload the aligned conformers from the
    # combined sdf file -- no reuse of step 1's in-memory objects, no
    # re-alignment.
    reference_name, analog_names = io_mcs_text.read_ligand_names(checkpoint)
    assert reference_name == REFERENCE_STEM
    assert set(analog_names) == set(ANALOG_STEMS)
    reference2, analogs2 = align.load_ligand_set_from_sdf(
        pipeline.aligned_sdf_path(aligned_dir), aligned_dir, reference_name, analog_names,
    )

    result = pipeline.build_from_core_definition(
        reference2, analogs2,
        sysname="fxa", outdir=outdir, cgenff=True,
        mcs_checkpoint_path=checkpoint, prm_search_dir=aligned_dir, toppar_search_dir=CHEMBL_DIR,
    )
    assert os.path.exists(os.path.join(outdir, "full_ligand.prm"))
    assert result.config.nsites == len(result.decomposition.sites)


def test_run_pipeline_result_charge_conservation(tmp_path):
    """End-to-end sanity check: the CRN invariant (each ligand's assembled
    charge equals its net charge) must hold on the full pipeline's output,
    not just when calling compute_crn directly."""
    reference, analogs = _fresh_ligands()
    outdir = str(tmp_path / "build")
    result = pipeline.run_pipeline(
        reference, analogs, sysname="fxa", outdir=outdir, cgenff=True, toppar_search_dir=CHEMBL_DIR,
        **_NO_GROUPING,
    )
    for ligand in result.ligands:
        total = sum(result.crn_result.core_charges.values())
        for site_charges in result.crn_result.fragment_charges:
            total += sum(site_charges[ligand.name].values())
        assert round(total, 6) == round(ligand.net_charge, 6)
