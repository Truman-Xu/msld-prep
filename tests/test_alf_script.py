from msldprep.config import SystemConfig
from msldprep.writers import alf_script


def test_write_alf_input_script_structure(tmp_path, decomposition, translation):
    config = SystemConfig(
        sysname="fxa", outdir=str(tmp_path), cgenff=True,
        nsubs=[len(site.fragments) for site in decomposition.sites],
    )
    path = alf_script.write_alf_input_script(config, decomposition, translation, str(tmp_path))
    text = open(path).read()

    assert text.startswith("* CHARMM input file for Multi-Site lambda-dynamics")
    assert f"set nsites = {config.nsites} " in text
    assert f"set nblocks = {config.nblocks} " in text
    assert "set fnex = 5.5" in text
    assert "END" in text


def test_baked_vs_runtime_config_values(tmp_path, decomposition, translation):
    """fnex/temp/box/ligseg/resnum/builddir get baked as literal `set` lines;
    nreps/nnodes/ncentral must stay as unresolved @var CHARMM references."""
    config = SystemConfig(
        sysname="fxa", outdir=str(tmp_path), cgenff=True,
        nsubs=[len(site.fragments) for site in decomposition.sites],
        nreps=3, nnodes=5, ncentral=2, temp=310.0,
    )
    path = alf_script.write_alf_input_script(config, decomposition, translation, str(tmp_path))
    text = open(path).read()

    assert "set temp = 310.0" in text
    assert "NREP @nreps" in text
    assert "set nreps" not in text
    assert "set nnodes" not in text
    assert "set ncentral" not in text
    assert "@ncentral" in text  # referenced inside the biasing-potential macro


def test_define_blocks_cover_every_fragment_atom(tmp_path, decomposition, translation):
    config = SystemConfig(
        sysname="fxa", outdir=str(tmp_path), cgenff=True,
        nsubs=[len(site.fragments) for site in decomposition.sites],
    )
    path = alf_script.write_alf_input_script(config, decomposition, translation, str(tmp_path))
    text = open(path).read()

    for site in decomposition.sites:
        for frag_number, (ligand_name, fragment) in enumerate(site.fragments.items(), start=1):
            sname = f"site{site.index + 1}_sub{frag_number}"
            assert f"define {sname} - " in text
            frag_translation = translation.fragments[(site.index, ligand_name)]
            for raw_name in fragment.atom_names:
                assert frag_translation[raw_name] in text


def test_copy_toppar_files_cgenff_false_warns_and_skips(tmp_path, caplog):
    import logging

    with caplog.at_level(logging.WARNING):
        alf_script.copy_toppar_files(str(tmp_path), cgenff=False)
    assert "Manually copy toppar files" in caplog.text
    assert not (tmp_path / "toppar").exists()
