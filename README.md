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

## Alignment GUI

A browser page for the alignment step only: load a reference with 3D coordinates (mol2 or sdf),
paste analog SMILES, align each analog onto the reference with a constrained embedding, look at
the results in 3D and 2D, accept or reject, and download the structures as sdf or mol2, in one
file or one file per structure. The mol2 files are written the way the notebook does it (Sybyl
types, no charges). No CGenFF typing is involved.

```bash
pip install fastapi uvicorn             # the `[web]` extra, if the environment lacks them
python -m msldprep.gui                  # prints the tunnel command and a URL with an access token
```

On your own computer run the printed `ssh -L ...` command, then open the printed address. The
server listens on localhost only, needs the token, and keeps everything in memory: nothing is
written unless you press the download button. The mol2 files get **new unique atom names** and
no charges, on purpose: CGenFF, which they are written for, only keeps unique short names and
computes its own charges (an option to keep the original names is planned, #10). Try it with `examples/2_ChEMBL/11896.mol2` as the
reference and `examples/gui_demo/analogs.txt` as the analogs. The 3D view needs a browser with
WebGL. 3Dmol.js (BSD-3-Clause) is bundled in `msldprep/gui/static/`, so the page loads nothing
from the internet.

## Test

```bash
pytest -q
```
