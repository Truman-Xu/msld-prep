#! /usr/bin/env python

####
#### Executable script to build MSLD-ready ligand files (msldprep rewrite)
####

"""
Runs the in-memory msldprep pipeline end to end in a single invocation:
align ligands, decompose into a shared core + per-site fragments,
renormalize charges, and write every CHARMM/ALF file.

If you want to review or hand-edit the automatic core/site/fragment
assignment before charges are renormalized and files are written, use
`pipeline.define_common_core` and `pipeline.build_from_core_definition`
directly instead of `run_pipeline` (or the CLI's `from-files`/`build`
subcommands -- see `msld_prep_cli.py`) -- `run_pipeline` (used below) is
the no-review, single-call convenience wrapper over both steps.

All ligand mol2/rtf files, each ligand's own CHARMM .prm parameter file,
and (if cgenff=True) a toppar/ directory with toppar.str must be available
in `base_dir` before running this script.

The `in_frag` list of lists allows you to move core atoms into alchemical
fragments at specific sites. Each nested list corresponds to a single site
attached to the ligand core. The `an_core` list of lists allows you to move
(non-default) anchor atoms at each site back into the core during charge
renormalization -- otherwise every site's anchor atom(s) move into each
ligand's fragment by default (see decompose.apply_manual_overrides).
"""

import logging
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from msldprep import align, pipeline  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

#####################################################################
## (1) Define system and file variables

sysname = "fxa"
base_dir = "../examples/2_ChEMBL"  # directory with {name}.mol2 / {name}.rtf / {name}.prm (+ toppar/)

reference_name = "11896"
analog_names = ["10749", "10322", "10870"]

outdir = f"build.{sysname}"

cgenff = True  # Are CGenFF/ParamChem parameters being used?

in_frag = None  # reflig core atoms to move into each fragment at each site (list of nsite lists)
an_core = None  # anchor atoms at each site to keep in the core (list of nsite lists)

#####################################################################
## (2) Load ligands and run the full pipeline in one call. The common-core
## definition it computes along the way is still written to
## outdir/MCS_for_MSLD.txt and outdir/aligned/ for later inspection/reuse
## (e.g. as input to a separate `build` step -- see this module's docstring).

reference = align.load_ligand(
    os.path.join(base_dir, f"{reference_name}.mol2"),
    os.path.join(base_dir, f"{reference_name}.rtf"),
    os.path.join(base_dir, f"{reference_name}.prm"),
    name=reference_name, is_reference=True,
)
analogs = [
    align.load_ligand(
        os.path.join(base_dir, f"{name}.mol2"),
        os.path.join(base_dir, f"{name}.rtf"),
        os.path.join(base_dir, f"{name}.prm"),
        name=name,
    )
    for name in analog_names
]

result = pipeline.run_pipeline(
    reference, analogs,
    sysname=sysname, outdir=outdir, cgenff=cgenff,
    in_frag=in_frag, an_core=an_core,
    toppar_search_dir=base_dir,
)

#####################################################################
## Finished

print(f"\nMSLD build written to {outdir}/")
print(f"default TOPPAR parameters copied into {outdir}. Check to make sure these work for your system!")
