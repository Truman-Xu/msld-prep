import os

from msldprep.config import SystemConfig, read_config_json, write_config_json, write_legacy_config_files


def test_system_config_derived_properties():
    config = SystemConfig(sysname="fxa", outdir="build.fxa", cgenff=True, nsubs=[3, 2, 4])
    assert config.nsites == 3
    assert config.nblocks == 9


def test_config_json_round_trip(tmp_path):
    config = SystemConfig(
        sysname="fxa", outdir=str(tmp_path), cgenff=True, nsubs=[4, 4, 4],
        nreps=2, nnodes=3, ncentral=1, fnex=6.0, box="50.0", temp=300.0,
    )
    write_config_json(config, str(tmp_path))
    reloaded = read_config_json(os.path.join(str(tmp_path), "config.json"))
    assert reloaded == config


def test_write_legacy_config_files(tmp_path):
    config = SystemConfig(sysname="fxa", outdir=str(tmp_path), cgenff=True, nsubs=[4], nnodes=2, ncentral=1)
    write_legacy_config_files(config, str(tmp_path))
    assert open(tmp_path / "name").read() == "fxa"
    assert open(tmp_path / "nnodes").read() == "2"
    assert open(tmp_path / "ncentral").read() == "1"
