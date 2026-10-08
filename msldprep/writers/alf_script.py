"""Write the CHARMM/ALF prep input script (`{sysname}.inp`) plus the small
config bookkeeping files ALF's external driver reads directly.

Ported from msld_wrt.py's `writeALF_Files`. The static CHARMM/ALF template
text (header, toppar/core-load boilerplate, the connectivity-deletion
macro, solvent/crystal/LDBI boilerplate, the ALF-3.0+ biasing-potential
macro, the NOE placeholder section) lives in standalone files under
`writers/templates/`, `${placeholder}`-substituted (stdlib `string.Template`,
`$var` syntax -- distinct from CHARMM's own `@var` syntax used throughout
the templates themselves, so there's no collision) with only the
`SystemConfig` fields that legacy baked as literal `set` lines
(`fnex`/`temp`/`box`/`ligseg`/`resnum`/`builddir`/`nsites`/`nblocks`).
`nreps`/`nnodes`/`ncentral` are never substituted into the script -- they
stay as unresolved `@nreps`/`@ncentral` CHARMM variable references for the
external ALF driver to set per-run, exactly as legacy left them.

Data-dependent pieces (the per-fragment rtf-read/patch/build loops, the
MSLD substituent `define` selections, the `ldin` UPOS/UNEG-by-parity
assignment, the exclusion-pair enumeration, and the per-block site-index
list for `msld 0 ...`) remain generated Python code -- ported directly from
`readAtoms`/the corresponding loops in msld_wrt.py, but reading atom names
from `Fragment.atom_names`/`lp_atom_names` (already in memory from Phase 1)
instead of re-reading `{Sname}_pres.rtf` off disk. A few single-line CHARMM
connectives (`rmla`/`msma`) stay as inline string literals -- too small to
warrant their own template file.

`WriteError` (unraised in legacy) becomes `errors.WriterError`; the
`os.system('cp ...')` toppar duplication becomes `shutil.copy`/`copytree`.
"""

from __future__ import annotations

import logging
import os
import shutil
from pathlib import Path
from string import Template
from typing import Sequence

from ..config import SystemConfig
from ..datamodel import CoreDecomposition
from .records import AtomNameTranslation

logger = logging.getLogger(__name__)

_TEMPLATE_DIR = Path(__file__).parent / "templates"


def _load_template(name: str) -> Template:
    return Template((_TEMPLATE_DIR / name).read_text())


def copy_toppar_files(outdir: str, *, cgenff: bool, search_dir: str = ".") -> None:
    """Ported from msld_wrt.py:40-46. `search_dir` defaults to the current
    working directory, matching legacy's cwd-relative `toppar/` lookup.
    """
    if not cgenff:
        logger.warning("Manually copy toppar files into %s!", outdir)
        logger.warning(
            "Scripts assume the NONBONDED line in parameter files is the same across all input ligands."
        )
        return
    toppar_dir = os.path.join(search_dir, "toppar")
    if not os.path.isdir(toppar_dir):
        logger.warning("No toppar/ directory found under %s; skipping toppar copy.", search_dir)
        return
    str_path = os.path.join(toppar_dir, "toppar.str")
    if os.path.exists(str_path):
        shutil.copy(str_path, outdir)
    shutil.copytree(toppar_dir, os.path.join(outdir, "toppar"), dirs_exist_ok=True)


def _snames(decomposition: CoreDecomposition) -> list[str]:
    return [
        f"site{site.index + 1}_sub{frag_number}"
        for site in decomposition.sites
        for frag_number in range(1, len(site.fragments) + 1)
    ]


def _atom_names_by_sname(
    decomposition: CoreDecomposition, translation: AtomNameTranslation
) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for site in decomposition.sites:
        for frag_number, (ligand_name, fragment) in enumerate(site.fragments.items(), start=1):
            frag_translation = translation.fragments[(site.index, ligand_name)]
            result[f"site{site.index + 1}_sub{frag_number}"] = [
                frag_translation[n] for n in [*fragment.atom_names, *fragment.lp_atom_names]
            ]
    return result


def _perturbation_variables_block(config: SystemConfig) -> str:
    lines = ["! perturbation variables\n", f"set nsites = {config.nsites} \n"]
    for i, n in enumerate(config.nsubs, start=1):
        lines.append(f"set nsubs{i} = {n} \n")
    lines.append(f"set nblocks = {config.nblocks} \n\n")
    return "".join(lines)


def _rtf_read_block(snames: Sequence[str]) -> str:
    lines = ["!! Read in the ligand patch rtf files\n"]
    lines += [f"read rtf append card name @builddir/{s}_pres.rtf\n" for s in snames]
    lines.append("\n")
    return "".join(lines)


def _patch_build_block(decomposition: CoreDecomposition) -> str:
    lines = ["!! Read in the ligand fragment pdb files\n", "ic generate\n\n"]
    for site in decomposition.sites:
        for frag_number in range(1, len(site.fragments) + 1):
            lines.append(f"patch p{site.index + 1}_{frag_number} @ligseg @resnum setup\n")
            lines.append(
                f"read coor pdb resid name @builddir/site{site.index + 1}_sub{frag_number}_frag.pdb\n"
            )
            lines.append("ic param\n")
            lines.append("ic build\n\n")
    return "".join(lines)


def _define_block(snames: Sequence[str], atom_names_by_sname: dict[str, list[str]]) -> str:
    lines = ["!! Define MSLD substituent selections\n"]
    for sname in snames:
        lines.append(f"define {sname} - \n   select ( - \n")
        for at in atom_names_by_sname[sname]:
            lines.append(f"   atom @ligseg @resnum {at} .or. -\n")
        lines.append("   none ) end\n\n")
    return "".join(lines)


def _call_block(snames: Sequence[str]) -> str:
    lines = [f"   Call {i + 2} sele {sname} show end\n" for i, sname in enumerate(snames)]
    lines.append("\n")
    return "".join(lines)


def _ldin_block(decomposition: CoreDecomposition) -> str:
    lines = []
    block = 1
    for site in decomposition.sites:
        nsub_here = len(site.fragments)
        for sub in range(1, nsub_here + 1):
            block += 1
            frac = 1.0 / nsub_here
            if sub == 1:
                rex = "NONE"
            elif sub % 2:
                rex = "UPOS 7.0"
            else:
                rex = "UNEG 7.0"
            lines.append(
                f"   ldin {block:<2d} {frac:.4f} 0.0  5.0  @lams{site.index + 1}s{sub} 5.0 !RX! {rex}\n"
            )
    lines.append("\n")
    return "".join(lines)


def _exclusion_block(decomposition: CoreDecomposition) -> str:
    lines: list[str] = []
    numexcl, excl = 1, 0
    lines.append(f"   set excl{numexcl} =")
    block0 = 1
    for site in decomposition.sites:
        nsub_here = len(site.fragments)
        for block1 in range(block0 + 1, block0 + nsub_here + 1):
            for block2 in range(block1 + 1, block0 + nsub_here + 1):
                lines.append(f" {block1} {block2}")
                excl += 1
                if excl == 11:
                    excl = 0
                    numexcl += 1
                    lines.append(f"\n   set excl{numexcl} =")
        block0 += nsub_here
    lines.append("\n   excl")
    lines += [f" @excl{i + 1}" for i in range(numexcl)]
    lines.append("\n\n")
    return "".join(lines)


def _msld_line_block(decomposition: CoreDecomposition) -> str:
    lines = ["   msld 0"]
    for site in decomposition.sites:
        lines += [f"  {site.index + 1}"] * len(site.fragments)
    lines.append("  fnex @fnex\n")
    lines.append("   msma\n\n")
    return "".join(lines)


def _noe_site_lines(decomposition: CoreDecomposition) -> str:
    lines = []
    for site in decomposition.sites:
        lines.append(f"!! ! Site {site.index + 1}\n")
        if len(site.fragments) >= 2:
            lines.append("!!    assign sele atom LIG 1 Atom1 end sele atom LIG 1 Atom2 end -\n")
            lines.append("!!    kmin 1.0 rmin 0.0 kmax 1.0 rmax 0.2 fmax 2.0 rswitch 2.0 sexp 1.0\n")
    return "".join(lines)


def write_alf_input_script(
    config: SystemConfig, decomposition: CoreDecomposition, translation: AtomNameTranslation, outdir: str
) -> str:
    """Assemble and write `outdir/{sysname}.inp`. Returns the path written."""
    snames = _snames(decomposition)
    atom_names_by_sname = _atom_names_by_sname(decomposition, translation)

    parts = [
        _load_template("header.charmm.tmpl").safe_substitute(
            fnex=config.fnex, sysname=config.sysname, builddir=config.builddir,
            box=config.box, temp=config.temp, ligseg=config.ligseg, resnum=config.resnum,
        ),
        _perturbation_variables_block(config),
        _load_template("toppar_and_core_load.charmm.tmpl").safe_substitute(),
        _rtf_read_block(snames),
        _patch_build_block(decomposition),
        "!! Read in LonePair sites (if applicable)\nstream @builddir/lpsites.inp\n\n",
        _define_block(snames, atom_names_by_sname),
        _load_template("system_assembly.charmm.tmpl").safe_substitute(),
        _load_template("block_header.charmm.tmpl").safe_substitute(),
        _call_block(snames),
        _load_template("block_settings.charmm.tmpl").safe_substitute(),
        _ldin_block(decomposition),
        _exclusion_block(decomposition),
        "   !!rmla bond thet dihe impr \n   rmla bond thet impr \n",
        _msld_line_block(decomposition),
        _load_template("biasing_potential.charmm.tmpl").safe_substitute(),
        _load_template("noe_placeholder.charmm.tmpl").safe_substitute(
            noe_site_lines=_noe_site_lines(decomposition)
        ),
    ]

    path = os.path.join(outdir, f"{config.sysname}.inp")
    with open(path, "w") as fp:
        fp.write("".join(parts))
    logger.info("Wrote %s", path)
    return path


def write_alf_files(
    config: SystemConfig,
    decomposition: CoreDecomposition,
    translation: AtomNameTranslation,
    outdir: str,
    *,
    toppar_search_dir: str = ".",
) -> None:
    """Convenience: copy_toppar_files + write_legacy_config_files + write_alf_input_script."""
    from ..config import write_legacy_config_files

    copy_toppar_files(outdir, cgenff=config.cgenff, search_dir=toppar_search_dir)
    write_legacy_config_files(config, outdir)
    write_alf_input_script(config, decomposition, translation, outdir)
