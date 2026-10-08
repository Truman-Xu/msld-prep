"""Build a fully-parameterized LigandRecord for a SMILES-only analog, or
for any RDKit Mol.

Bridges "I have a SMILES string" (or "I have an RDKit Mol") to "I have an
aligned, CGenFF-typed LigandRecord ready for the rest of the pipeline":

1. Build a 3D conformer tethered to a reference ligand's core
   (`align.align_common_core`, SMILES path only).
2. Type it with Sybyl atom/bond types and unique atom names entirely
   in-memory (`sybyl_typer.annotate_sybyl_types`/`assign_unique_atom_names`)
   and format it to mol2 text (`sybyl_typer.format_mol2_text`) -- no
   `obabel` CLI dependency, no intermediate SDF or mol2 file.
3. Parameterize by piping that mol2 text directly into a locally-installed
   CGenFF binary's stdin and capturing its unified `.str` stream from
   stdout (`mol_to_ligand_record`) -- cgenff's default (non `-x`) mode
   already streams to stdout/reads from stdin (confirmed: byte-identical
   output to file-based invocation), so nothing cgenff produces ever
   touches disk either; the stream is split into its rtf/param blocks and
   parsed straight into a `LigandRecord` (`rtf_records.parse_rtf_text` /
   `prm_records.parse_prm_text`).
4. (`build_analog_from_smiles` only) Materialize `.mol2`/`.rtf`/`.prm`
   files into `workdir` afterward, via `sybyl_typer.write_mol2_file` and
   `LigandRecord.write_rtf`/`write_prm` -- so the on-disk artifact contract
   existing callers (`pipeline.run_pipeline_from_smiles`'s `prm_search_dir`
   convention, user inspection) rely on is preserved, even though nothing
   was written to disk *during* CGenFF parameterization itself.
"""

from __future__ import annotations

import logging
import os
import subprocess
import warnings
from dataclasses import dataclass, field

from rdkit import Chem

from . import prm_records, rtf_records
from .align import align_common_core
from .datamodel import LigandRecord
from .errors import AlignmentError, ExternalToolError, MsldPrepError
from .sybyl_typer import annotate_sybyl_types, assign_unique_atom_names, format_mol2_text, write_mol2_file

logger = logging.getLogger(__name__)


@dataclass
class ExternalTools:
    """Paths to the external tools this module shells out to."""

    cgenff: str = "cgenff"
    cgenff_base_prm: str | None = None  # passed as `-p charmm <file>`; omit to skip


def generate_3d_analog(smiles: str, reference: LigandRecord, *, name: str) -> Chem.Mol:
    """Build a 3D-embedded RDKit Mol for `smiles`, tethered to `reference`'s core."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        msg = f"Could not parse SMILES for analog {name}: {smiles!r}"
        logger.error(msg)
        raise AlignmentError(msg)
    mol = Chem.AddHs(mol)
    try:
        logger.info(f'Aligning common core for {name}')
        align_common_core(mol, reference.mol)
    except (ValueError, RuntimeError) as exc:
        msg = f"Could not align analog {name} ({smiles}) to the reference: {exc!r}"
        logger.error(msg)
        raise AlignmentError(msg) from exc
    return mol


def _run_cgenff_streamed(cgenff: str, mol2_text: str, *, cgenff_base_prm: str | None = None) -> str:
    """Run cgenff on `mol2_text` piped via stdin, returning its captured
    stdout (cgenff's default mode: a unified .str stream). No file I/O --
    confirmed cgenff's default stdin/stdout behavior produces output
    byte-identical to its file-based invocation for the same molecule.
    """
    cmd = [cgenff]
    if cgenff_base_prm:
        cmd += ["-p", "charmm", cgenff_base_prm]
    try:
        result = subprocess.run(cmd, input=mol2_text, capture_output=True, text=True, timeout=600)
    except FileNotFoundError as exc:
        msg = f"cgenff not found: `{cgenff}`"
        logger.error(msg)
        raise ExternalToolError(msg) from exc
    except subprocess.TimeoutExpired as exc:
        msg = f"cgenff timed out: {' '.join(cmd)}"
        logger.error(msg)
        raise ExternalToolError(msg) from exc
    if result.returncode != 0:
        msg = f"cgenff failed (exit {result.returncode}): {' '.join(cmd)}\n{result.stderr}"
        logger.error(msg)
        raise ExternalToolError(msg)
    return result.stdout


def _split_cgenff_stream(text: str) -> tuple[str, str]:
    """Split a CGenFF unified .str stream into its rtf block (from "read
    rtf card append" onward) and param block (from "read param card ...
    append" onward). Each is independently parseable as-is: `RTFParser`
    already special-cases a "read rtf" line, and `prm_records.parse_prm_lines`
    strips any "read " wrapper line itself; both parsers stop cleanly at
    their own embedded "END" regardless of what text (the other block,
    the trailing "RETURN") follows it.
    """
    try:
        rtf_start = text.index("read rtf card append")
        param_start = text.index("read param", rtf_start)
    except ValueError as exc:
        msg = "cgenff output did not look like a CGenFF .str stream (missing rtf/param markers)"
        logger.error(msg)
        raise ExternalToolError(msg) from exc
    return text[rtf_start:], text[param_start:]


def mol_to_ligand_record(
    mol: Chem.Mol, name: str, cgenff: str, *, cgenff_base_prm: str | None = None, is_reference: bool = False
) -> LigandRecord:
    """Type `mol` (one conformer, explicit Hs) via a streamed CGenFF call
    and return a complete `LigandRecord` -- no .mol2/.rtf/.prm files ever
    touch disk. See module docstring for the full streaming path.
    """
    annotate_sybyl_types(mol)
    assign_unique_atom_names(mol)
    mol2_text = format_mol2_text(mol, mol_name=name)

    stream = _run_cgenff_streamed(cgenff, mol2_text, cgenff_base_prm=cgenff_base_prm)
    rtf_block, param_block = _split_cgenff_stream(stream)

    source_label = f"<cgenff stream for {name}>"
    rtf_fields = rtf_records.parse_rtf_text(rtf_block, source_label=source_label)
    rtf_records.assign_atom_types(mol, rtf_fields["atom_types"], rtf_fields["atom_charges"], source_label=source_label)
    prm = prm_records.parse_prm_text(param_block)

    ligand = LigandRecord(name=name, mol=mol, prm=prm, is_reference=is_reference, **rtf_fields)
    logger.info("Typed %s via streamed CGenFF (net_charge=%s)", name, ligand.net_charge)
    return ligand


def build_analog_from_smiles(
    smiles: str,
    name: str,
    reference: LigandRecord,
    workdir: str,
    *,
    tools: ExternalTools = ExternalTools(),
) -> LigandRecord:
    """SMILES -> aligned 3D -> streamed CGenFF (`mol_to_ligand_record`,
    nothing touches disk) -> LigandRecord, then materialize `.mol2`/`.rtf`/
    `.prm` into `workdir` from that LigandRecord for downstream consumers
    that still expect real files (e.g. `pipeline.run_pipeline_from_smiles`'s
    `prm_search_dir` convention) and for user inspection.
    """
    os.makedirs(workdir, exist_ok=True)
    mol = generate_3d_analog(smiles, reference, name=name)
    ligand = mol_to_ligand_record(mol, name, tools.cgenff, cgenff_base_prm=tools.cgenff_base_prm)

    ligand.mol2_path = os.path.join(workdir, f"{name}.mol2")
    ligand.rtf_path = os.path.join(workdir, f"{name}.rtf")
    ligand.prm_path = os.path.join(workdir, f"{name}.prm")
    write_mol2_file(ligand.mol, ligand.mol2_path, mol_name=name)
    ligand.write_rtf(ligand.rtf_path)
    ligand.write_prm(ligand.prm_path)

    return ligand


@dataclass
class AnalogBuildResult:
    ligands: list[LigandRecord] = field(default_factory=list)
    failed: list[tuple[str, Exception]] = field(default_factory=list)


def build_analogs_from_smiles(
    smiles_by_name: dict[str, str],
    reference: LigandRecord,
    workdir: str,
    *,
    tools: ExternalTools = ExternalTools(),
) -> AnalogBuildResult:
    """Build every analog in `smiles_by_name`, collecting failures instead of raising."""
    result = AnalogBuildResult()
    for name, smiles in smiles_by_name.items():
        try:
            result.ligands.append(build_analog_from_smiles(smiles, name, reference, workdir, tools=tools))
        except MsldPrepError as exc:
            logger.warning("Failed to build analog %s from SMILES %r: %r", name, smiles, exc)
            result.failed.append((name, exc))

    if result.failed:
        names = [n for n, _ in result.failed]
        warnings.warn(
            f"{len(result.failed)}/{len(smiles_by_name)} analog(s) failed to build from SMILES: {names}",
            stacklevel=2,
        )
    return result
