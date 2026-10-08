# msld-prep

RDKit-based ligand preparation for multisite lambda dynamics (MSLD) in CHARMM.

Given a reference ligand and a series of analogs that already have CGenFF mol2, rtf and prm
files (or SMILES that CGenFF can parameterize), msldprep:

1. loads and types each ligand, and aligns the analogs in 3D onto the reference;
2. finds the shared core and decomposes each analog into anchors, sites and fragments;
3. renormalizes charges so each fragment is consistent across the series;
4. writes the CHARMM rtf, prm and pdb files and the ALF and `lpsites.inp` inputs.

It replaces the older `msld-py-prep` scripts, and uses crimm's RTF and PRM parsers.
Repository: https://github.com/Truman-Xu/msld-prep

## Install

```bash
conda activate crimm-dev        # the development environment shared with crimm
pip install -e ".[test]" --no-deps
```

The Python package is `msldprep`. The CGenFF program is needed for the SMILES path and is
installed separately; the tool and the tests find it through `MSLDPREP_CGENFF_PATH` (see CONTRIBUTING.md).

## License

GPLv3. See [LICENSE](LICENSE).

## Use

The command line has two steps, with a reviewable checkpoint between them:

```bash
python msld_prep_cli.py from-files --config configs/from_files_example.yaml   # or from-smiles
less work.fxa/MCS_for_MSLD.txt                                                # review or edit
python msld_prep_cli.py build --config configs/build_example.yaml
```

`run_msld_prep.py` and `run_msld_prep_from_smiles.py` run both steps in one call. Annotated
configs are in `configs/`.

## Test

```bash
pytest -q
```
