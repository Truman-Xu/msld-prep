"""Orchestration: the full msldprep pipeline from loaded ligands to a
complete MSLD build directory, split into two explicit steps:

1. `define_common_core` -- align every analog onto the reference, find and
   decompose the shared core, and export the result for review: a
   human-editable checkpoint text file (`io_mcs_text.export_mcs_text`) plus
   a self-contained ligand directory (`{aligned_dir}/aligned_ligands.sdf` +
   `{name}.rtf`/`.prm` per ligand) holding each ligand's ALIGNED conformer.
2. `build_from_core_definition` -- import that (possibly hand-edited)
   checkpoint, apply manual core/fragment overrides, renormalize charges,
   and write every CHARMM/ALF file.

These replace an earlier single `run_pipeline` that decided which of the
two to do based on whether a checkpoint file already existed on disk (a
legacy `msld_py_prep.py`-style "run the same command twice" trick) -- that
implicit, file-existence-driven behavior is gone; callers now call exactly
the step(s) they want, explicitly. `run_pipeline`/`run_pipeline_from_smiles`
remain as thin, in-process convenience wrappers over the two steps for
callers who don't need to review anything in between (most of the test
suite, casual scripts).

Why step 2 needs the *aligned* conformers, not the originals: RDKit's
`ConstrainedEmbed` (used to tether each analog's 3D embedding to the
reference's core) isn't perfectly deterministic, so re-running alignment
in a second, separate process could produce a conformer that's subtly
different from -- or occasionally fails where the first run succeeded --
what was actually reviewed in the checkpoint. Loading step 2's ligands from
`define_common_core`'s own `{aligned_dir}` output sidesteps that entirely:
the exact conformer that was decomposed and checkpointed is what gets
written to the final PDB files, whether step 2 runs in the same process or
a completely separate one. Those aligned conformers are written and
reloaded as a single combined sdf file rather than one mol2 file per
ligand, since RDKit's mol2 reader occasionally fails to re-kekulize an
aromatic ring reloaded that way (see `align.write_ligands_sdf`).

Used by both the test suite and `rewrite/msld_prep_cli.py`.
"""

from __future__ import annotations

import logging
import os
import shutil
from dataclasses import dataclass

from . import align, crn, decompose, grouping, io_mcs_text
from . import analog_prep
from . import prm_records as prm_mod
from .config import SystemConfig, write_config_json
from .datamodel import CoreDecomposition, CRNResult, LigandRecord
from .errors import DecompositionError
from .writers import alf_script, charmm_control, fragment_json, pdb
from .writers import prm as prm_writer
from .writers import rtf as rtf_writer
from .writers.records import compute_atom_name_translation

logger = logging.getLogger(__name__)


@dataclass
class PipelineResult:
    ligands: list[LigandRecord]
    decomposition: CoreDecomposition
    crn_result: CRNResult
    config: SystemConfig


@dataclass
class CoreDefinition:
    ligands: list[LigandRecord]  # aligned; mol2_path now points at aligned_dir's combined sdf
    decompositions: list[CoreDecomposition]  # one per core group (grouping.group_ligands_by_similarity)
    mcs_checkpoint_path: str
    aligned_dir: str

    @property
    def decomposition(self) -> CoreDecomposition:
        """Convenience for the common single-group case; raises if grouping
        actually split this ligand set into more than one core group -- use
        `.decompositions` (or the `build` CLI command / a `core_id=`-scoped
        `build_from_core_definition` call per group) instead.
        """
        if len(self.decompositions) != 1:
            core_ids = [d.core_id for d in self.decompositions]
            raise DecompositionError(
                f"This ligand set was split into {len(self.decompositions)} core groups "
                f"{core_ids}, not one -- use `.decompositions` and build each group "
                "separately (e.g. the `build` CLI command, which does this automatically)."
            )
        return self.decompositions[0]


ALIGNED_SDF_FILENAME = "aligned_ligands.sdf"


def aligned_sdf_path(aligned_dir: str) -> str:
    """Path to the combined sdf file `_write_aligned_ligand_dir` writes
    inside `aligned_dir` -- what a step-2 (`build`) caller should pass to
    `align.load_ligand_set_from_sdf` to reload these ligands later.
    """
    return os.path.join(aligned_dir, ALIGNED_SDF_FILENAME)


def _write_aligned_ligand_dir(ligands: list[LigandRecord], aligned_dir: str) -> None:
    """Write every ligand's current (aligned) conformer into one combined
    `aligned_sdf_path(aligned_dir)` sdf file (see `align.write_ligands_sdf`
    for why sdf rather than one mol2 file per ligand), and make sure its
    .rtf/.prm are there too -- copied from wherever they already live, or
    (for a ligand with no rtf_path/prm_path at all, e.g. built purely
    in-memory via `analog_prep.mol_to_ligand_record` directly) materialized
    fresh via `LigandRecord.write_rtf`/`write_prm`. Either way, `aligned_dir`
    alone ends up with everything a later, separate `build_from_core_definition`
    call needs to reload these exact ligands. Mutates each `ligand.mol2_path`
    to the combined sdf path, and `rtf_path`/`prm_path` to point at `aligned_dir`.
    """
    os.makedirs(aligned_dir, exist_ok=True)
    for ligand in ligands:
        rtf_path = os.path.join(aligned_dir, f"{ligand.name}.rtf")
        if ligand.rtf_path and os.path.abspath(ligand.rtf_path) != os.path.abspath(rtf_path):
            shutil.copy(ligand.rtf_path, rtf_path)
        elif not ligand.rtf_path:
            ligand.write_rtf(rtf_path)
        ligand.rtf_path = rtf_path

        prm_path = os.path.join(aligned_dir, f"{ligand.name}.prm")
        if ligand.prm_path and os.path.abspath(ligand.prm_path) != os.path.abspath(prm_path):
            shutil.copy(ligand.prm_path, prm_path)
        elif not ligand.prm_path:
            ligand.write_prm(prm_path)
        ligand.prm_path = prm_path

    sdf_path = aligned_sdf_path(aligned_dir)
    align.write_ligands_sdf(ligands, sdf_path)
    for ligand in ligands:
        ligand.mol2_path = sdf_path


def define_common_core(
    reference: LigandRecord,
    analogs: list[LigandRecord],
    *,
    mcs_checkpoint_path: str,
    aligned_dir: str,
    skip_pairwise_alignment: bool = False,
    similarity_threshold: float = grouping.DEFAULT_SIMILARITY_THRESHOLD,
) -> CoreDefinition:
    """Step 1 of 2: align `analogs` onto `reference`, split them into one or
    more core groups, find/decompose each group's own shared core, and
    export the result for `build_from_core_definition` to consume (see
    module docstring for why this includes the aligned conformers, not just
    the checkpoint text).

    `reference`/`analogs` should be freshly loaded (via `align.load_ligand`)
    -- `align_common_core` mutates each ligand's `mol` (embeds a conformer)
    in place.

    Set `skip_pairwise_alignment=True` if every one of `analogs` has
    already been individually aligned to `reference` (e.g. by
    `analog_prep.build_analog_from_smiles`) -- this skips the redundant
    second `align_common_core` embed `align.align_ligand_set` would
    otherwise do. `analogs` is then used as-is, with no alignment-failure
    filtering.

    `similarity_threshold` is `grouping.group_ligands_by_similarity`'s
    Tanimoto/Morgan cutoff for splitting `analogs` into per-group shared
    cores instead of forcing one global MCS across all of them (see that
    function's docstring for why); pass 0 to disable grouping and always
    compute the single global core, matching pre-grouping behavior.
    """
    if skip_pairwise_alignment:
        aligned_analogs = analogs
    else:
        aligned_analogs = align.align_ligand_set(analogs, reference).aligned
    all_ligands = [reference, *aligned_analogs]

    groups, excluded = grouping.group_ligands_by_similarity(
        reference, aligned_analogs, threshold=similarity_threshold
    )
    lone_groups = [g for g in groups if not g.analogs]
    if lone_groups:
        logger.warning(
            "%d core group(s) ended up with only their own reference ligand (no shared "
            "core to find against a single ligand) and were dropped: %s",
            len(lone_groups), [g.reference.name for g in lone_groups],
        )
        excluded = [*excluded, *(g.reference for g in lone_groups)]
    groups = [g for g in groups if g.analogs]
    if not groups:
        msg = (
            f"No usable core group survived similarity grouping (threshold {similarity_threshold}) "
            f"for reference {reference.name!r} -- every analog ended up without a partner. "
            "Lower --similarity-threshold (or pass 0 to disable grouping) and try again."
        )
        logger.error(msg)
        raise DecompositionError(msg)

    decompositions = []
    for group in groups:
        group_ligands = [group.reference, *group.analogs]
        group_core_query_mol = align.type_ligands_and_find_core(group_ligands)
        decompositions.append(decompose.decompose_ligand_set(
            group_ligands, group_core_query_mol, group.reference.name, core_id=group.core_id,
        ))

    checkpoint_dir = os.path.dirname(mcs_checkpoint_path)
    if checkpoint_dir:
        os.makedirs(checkpoint_dir, exist_ok=True)
    io_mcs_text.export_mcs_text(
        decompositions, mcs_checkpoint_path,
        excluded=[lig.name for lig in excluded], similarity_threshold=similarity_threshold,
    )
    _write_aligned_ligand_dir(all_ligands, aligned_dir)

    logger.info(
        "Exported %d core group definition(s) to %s and aligned ligand files to %s/ -- "
        "review/edit the checkpoint if needed, then run build_from_core_definition "
        "(or the `build` CLI command) against %s/.",
        len(decompositions), mcs_checkpoint_path, aligned_dir, aligned_dir,
    )
    return CoreDefinition(
        ligands=all_ligands, decompositions=decompositions,
        mcs_checkpoint_path=mcs_checkpoint_path, aligned_dir=aligned_dir,
    )


def build_from_core_definition(
    reference: LigandRecord,
    analogs: list[LigandRecord],
    *,
    sysname: str,
    outdir: str,
    cgenff: bool,
    mcs_checkpoint_path: str,
    core_id: int | None = None,
    in_frag: list[list[str]] | None = None,
    an_core: list[list[str]] | None = None,
    prm_search_dir: str = ".",
    toppar_search_dir: str = ".",
) -> PipelineResult:
    """Step 2 of 2: import the (possibly hand-edited) common-core
    definition from `mcs_checkpoint_path`, apply manual overrides,
    renormalize charges, and write every CHARMM/ALF file.

    Does no alignment of its own: `reference`/`analogs` must carry the
    conformer `define_common_core` computed (i.e. typically loaded via
    `align.load_ligand_set_from_sdf` from the `aligned_sdf_path(aligned_dir)`
    file it wrote), not freshly re-embedded -- see the module docstring for why.

    `reference`/`analogs` must be exactly the ligand set for one core group
    -- if `mcs_checkpoint_path` holds more than one (see
    `io_mcs_text.read_core_ids`), pass `core_id` to pick which, and build
    each group separately (the `build` CLI command does this automatically).
    """
    os.makedirs(outdir, exist_ok=True)
    all_ligands = [reference, *analogs]
    decomposition = io_mcs_text.import_mcs_text(mcs_checkpoint_path, all_ligands, core_id=core_id)
    decomposition = decompose.apply_manual_overrides(decomposition, all_ligands, in_frag, an_core)
    crn_result = crn.compute_crn(all_ligands, decomposition)

    translation = compute_atom_name_translation(decomposition, all_ligands)
    pdb.write_all_pdb_files(decomposition, all_ligands, translation, outdir)
    rtf_writer.write_core_rtf(decomposition, all_ligands, crn_result, translation, outdir)
    rtf_writer.write_all_pres_rtf(decomposition, all_ligands, crn_result, translation, outdir)
    charmm_control.write_lpsites_inp(decomposition, all_ligands, crn_result.core_charges, translation, outdir)
    charmm_control.write_translation_table(decomposition, translation, outdir)
    charmm_control.write_block_counts(decomposition, outdir)

    core_data, fragment_data = fragment_json.compute_all_file_data(
        decomposition, all_ligands, crn_result, translation
    )
    fragment_json.write_records_json(core_data, fragment_data, outdir)

    ligand_names_in_order = [lig.name for lig in all_ligands]
    prm_by_ligand = prm_mod.load_prm_files(ligand_names_in_order, prm_search_dir)
    merged = prm_mod.merge_prm_records(prm_by_ligand, cgenff=cgenff, ligand_order=ligand_names_in_order)
    topology = prm_mod.build_ligand_topology(decomposition, all_ligands)
    needed = prm_mod.derive_needed_parameters(topology, merged, cgenff=cgenff)
    prm_writer.write_full_ligand_prm(needed, outdir, cgenff=cgenff)

    config = SystemConfig(
        sysname=sysname, outdir=outdir, cgenff=cgenff,
        nsubs=[len(site.fragments) for site in decomposition.sites],
    )
    write_config_json(config, outdir)
    alf_script.write_alf_files(
        config, decomposition, translation, outdir, toppar_search_dir=toppar_search_dir
    )

    logger.info("Pipeline complete: %s", outdir)
    return PipelineResult(
        ligands=all_ligands, decomposition=decomposition, crn_result=crn_result, config=config
    )


def define_common_core_from_smiles(
    reference: LigandRecord,
    analog_smiles: dict[str, str],
    *,
    mcs_checkpoint_path: str,
    aligned_dir: str,
    workdir: str,
    tools: analog_prep.ExternalTools = analog_prep.ExternalTools(),
    similarity_threshold: float = grouping.DEFAULT_SIMILARITY_THRESHOLD,
) -> CoreDefinition:
    """Step 1 of 2, from SMILES: build every analog in `analog_smiles`
    (name -> SMILES) via a streamed local CGenFF call
    (`analog_prep.build_analogs_from_smiles`, materializing each analog's
    aligned mol2/rtf/prm into `workdir`), then `define_common_core` on the
    result (`skip_pairwise_alignment=True` -- each analog was already
    aligned to `reference` while it was being built).

    Analogs that fail to build (bad SMILES, alignment failure, or an
    external-tool error) are logged and excluded rather than raising --
    see `analog_prep.build_analogs_from_smiles` for details. See
    `define_common_core` for `similarity_threshold`.
    """
    os.makedirs(workdir, exist_ok=True)
    built = analog_prep.build_analogs_from_smiles(analog_smiles, reference, workdir, tools=tools)
    return define_common_core(
        reference, built.ligands, mcs_checkpoint_path=mcs_checkpoint_path, aligned_dir=aligned_dir,
        skip_pairwise_alignment=True, similarity_threshold=similarity_threshold,
    )


def run_pipeline(
    reference: LigandRecord,
    analogs: list[LigandRecord],
    *,
    sysname: str,
    outdir: str,
    cgenff: bool,
    in_frag: list[list[str]] | None = None,
    an_core: list[list[str]] | None = None,
    prm_search_dir: str | None = None,
    toppar_search_dir: str = ".",
    skip_pairwise_alignment: bool = False,
    similarity_threshold: float = grouping.DEFAULT_SIMILARITY_THRESHOLD,
) -> PipelineResult:
    """Convenience: `define_common_core` then immediately
    `build_from_core_definition`, in one process, with no manual review in
    between -- for callers who just want a complete build in a single call.

    Only works for a ligand set that ends up as exactly one core group
    (the common case, and the only case that needs no review): a set that
    similarity grouping actually splits into several groups raises, since
    picking apart which group's fragments mean what -- and whether the
    auto-picked local reference per group makes sense -- is exactly the
    kind of thing the explicit `define_common_core` + `build_from_core_definition`
    two-step workflow (or the `build` CLI command, which builds every group
    automatically) exists for.

    The checkpoint and aligned-ligand directory `define_common_core` always
    produces still land in `{outdir}/MCS_for_MSLD.txt` and `{outdir}/aligned/`
    respectively (in case you want to inspect or reuse them later); this
    function reuses the already-aligned `LigandRecord`s directly rather
    than reloading from those files, since it never leaves the process.

    `prm_search_dir` defaults to `{outdir}/aligned` (where `define_common_core`
    just copied every ligand's .prm) -- pass it explicitly only if your
    ligands' .prm files should be sourced from somewhere else instead.
    """
    mcs_checkpoint_path = os.path.join(outdir, "MCS_for_MSLD.txt")
    aligned_dir = os.path.join(outdir, "aligned")
    core = define_common_core(
        reference, analogs, mcs_checkpoint_path=mcs_checkpoint_path, aligned_dir=aligned_dir,
        skip_pairwise_alignment=skip_pairwise_alignment, similarity_threshold=similarity_threshold,
    )
    decomposition = core.decomposition  # raises if grouping produced more than one core group
    group_names = set(decomposition.core_atom_names)
    core_reference = next(lig for lig in core.ligands if lig.name == decomposition.reference_name)
    core_analogs = [
        lig for lig in core.ligands if lig.name in group_names and lig.name != decomposition.reference_name
    ]
    return build_from_core_definition(
        core_reference, core_analogs,
        sysname=sysname, outdir=outdir, cgenff=cgenff,
        mcs_checkpoint_path=mcs_checkpoint_path, core_id=decomposition.core_id,
        in_frag=in_frag, an_core=an_core,
        prm_search_dir=prm_search_dir or aligned_dir, toppar_search_dir=toppar_search_dir,
    )


def run_pipeline_from_smiles(
    reference: LigandRecord,
    analog_smiles: dict[str, str],
    *,
    sysname: str,
    outdir: str,
    workdir: str,
    cgenff: bool = True,
    tools: analog_prep.ExternalTools = analog_prep.ExternalTools(),
    in_frag: list[list[str]] | None = None,
    an_core: list[list[str]] | None = None,
    toppar_search_dir: str = ".",
    similarity_threshold: float = grouping.DEFAULT_SIMILARITY_THRESHOLD,
) -> PipelineResult:
    """Convenience: `define_common_core_from_smiles` then immediately
    `build_from_core_definition`, in one process, with no manual review in
    between. See `run_pipeline` for why this only works for a ligand set
    that ends up as exactly one core group.

    `reference`'s own `.rtf`/`.prm` are expected to already exist (it's a
    normal, already-parameterized `LigandRecord`, e.g. from `align.load_ligand`).
    For visibility into exactly which analogs failed to build, call
    `analog_prep.build_analogs_from_smiles` yourself and pass its `.ligands`
    to `run_pipeline` directly instead of using this wrapper.
    """
    mcs_checkpoint_path = os.path.join(outdir, "MCS_for_MSLD.txt")
    aligned_dir = os.path.join(outdir, "aligned")
    core = define_common_core_from_smiles(
        reference, analog_smiles, mcs_checkpoint_path=mcs_checkpoint_path, aligned_dir=aligned_dir,
        workdir=workdir, tools=tools, similarity_threshold=similarity_threshold,
    )
    decomposition = core.decomposition  # raises if grouping produced more than one core group
    group_names = set(decomposition.core_atom_names)
    core_reference = next(lig for lig in core.ligands if lig.name == decomposition.reference_name)
    core_analogs = [
        lig for lig in core.ligands if lig.name in group_names and lig.name != decomposition.reference_name
    ]
    return build_from_core_definition(
        core_reference, core_analogs,
        sysname=sysname, outdir=outdir, cgenff=cgenff,
        mcs_checkpoint_path=mcs_checkpoint_path, core_id=decomposition.core_id,
        in_frag=in_frag, an_core=an_core,
        prm_search_dir=aligned_dir, toppar_search_dir=toppar_search_dir,
    )
