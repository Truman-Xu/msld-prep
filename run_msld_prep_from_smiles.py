#! /usr/bin/env python

####
#### Build MSLD-ready ligand files starting from a reference ligand plus
#### analogs given only as SMILES strings (msldprep rewrite)
####

"""
For each analog SMILES: build a 3D structure aligned to the reference
ligand's core, write it directly to a mol2 file with a custom Sybyl atom
typer/writer (reproduces OpenBabel's own atom/bond typing, no `obabel` CLI
dependency), and parameterize it with a locally-installed CGenFF binary --
then run the analogs through the same pipeline as run_msld_prep.py (MCS
core-finding, per-site fragment decomposition, charge renormalization, and
every CHARMM/ALF file).

The reference ligand must already have its own .mol2/.rtf/.prm (e.g. an
existing CGenFF-parameterized structure); analogs need only a name and a
SMILES string -- everything else is generated.

Requires a local install of CGenFF (the standalone `cgenff` binary, not the
ParamChem web service), for atom typing/parameterization by analogy. Point to it
with the environment variable MSLDPREP_CGENFF_PATH (the program) and, optionally,
MSLDPREP_CGENFF_BASE_PRM (its par_all36_cgenff.prm; default: the file next to the program).
"""

import logging
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from msldprep import align, analog_prep, pipeline  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

#####################################################################
## (1) Define system and file variables

sysname = "fxa_smiles"

reference_dir = "../examples/2_ChEMBL"  # where {reference_name}.mol2/.rtf/.prm already live
reference_name = "11896"

outdir = f"build.{sysname}"
workdir = f"work.{sysname}"  # generated analog .sdf/.mol2/.rtf/.prm land here

cgenff = True

# name -> SMILES for every analog to generate; edit freely
analog_smiles = {
    "f_analog": "Fc1ccc(Nc2ncnc3ccccc23)cc1",
    "cl_analog": "Clc1ccc(Nc2ncnc3ccccc23)cc1",
    "me_analog": "Cc1ccc(Nc2ncnc3ccccc23)cc1",
}

in_frag = None
an_core = None

cgenff_bin = os.environ.get("MSLDPREP_CGENFF_PATH")
if not cgenff_bin:
    raise SystemExit("Set MSLDPREP_CGENFF_PATH to the standalone cgenff program.")
tools = analog_prep.ExternalTools(
    cgenff=cgenff_bin,
    cgenff_base_prm=os.environ.get("MSLDPREP_CGENFF_BASE_PRM")
    or os.path.join(os.path.dirname(cgenff_bin), "par_all36_cgenff.prm"),
)

#####################################################################
## (2) Load the reference, build every analog from SMILES, and run the
## full pipeline. Analogs that fail to build (bad SMILES, alignment
## failure, or an external-tool error) are logged and excluded; see
## analog_prep.build_analogs_from_smiles if you need that detail
## programmatically instead of via this convenience wrapper.

reference = align.load_ligand(
    os.path.join(reference_dir, f"{reference_name}.mol2"),
    os.path.join(reference_dir, f"{reference_name}.rtf"),
    os.path.join(reference_dir, f"{reference_name}.prm"),
    name=reference_name, is_reference=True,
)

result = pipeline.run_pipeline_from_smiles(
    reference, analog_smiles,
    sysname=sysname, outdir=outdir, workdir=workdir,
    cgenff=cgenff, tools=tools,
    in_frag=in_frag, an_core=an_core,
    toppar_search_dir=reference_dir,
)

#####################################################################
## Finished

print(f"\nMSLD build written to {outdir}/")
print(f"Generated analog files (.mol2/.rtf/.prm) kept in {workdir}/")
print(f"default TOPPAR parameters copied into {outdir}. Check to make sure these work for your system!")
