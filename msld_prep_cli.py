#! /usr/bin/env python

####
#### Command-line tool for the msldprep pipeline
####

"""
Three subcommands, mirroring `pipeline.py`'s explicit two-step split (no
implicit "run the same command twice" workflow -- each command does
exactly what its name says, every time):

  from-files    step 1: define the common core from ligands that already
                have mol2/rtf/prm files (pipeline.define_common_core).
                Analogs too structurally different to share one core with
                the reference are automatically split into their own core
                group(s) by Tanimoto/Morgan similarity -- see
                --similarity-threshold -- each recorded as its own "CORE
                ID" section in the checkpoint.
  from-smiles   step 1: define the common core, building a reference
                ligand's analogs from SMILES via a streamed local CGenFF
                call first (pipeline.define_common_core_from_smiles). Same
                --similarity-threshold grouping as from-files.
  build         step 2: read a common-core definition written by either of
                the above (optionally hand-edited in between), apply
                overrides, renormalize charges, and write every CHARMM/ALF
                file (pipeline.build_from_core_definition) -- once per CORE
                ID group found in the checkpoint (each into its own
                outdir/core<N>/ when there's more than one; pass --core-id
                to build just one). Ligand names are read from the
                checkpoint itself -- no need to repeat --reference/--analogs.

Every option can be given on the command line, in a YAML config file via
`--config`, or both -- config-file values are used as defaults, and any
flag actually given on the command line overrides them. See
`configs/from_files_example.yaml` / `configs/from_smiles_example.yaml` /
`configs/build_example.yaml` for annotated examples.

Examples:
    # Review the automatic core/site/fragment assignment before building:
    python msld_prep_cli.py from-files --config configs/from_files_example.yaml
    less work.fxa/MCS_for_MSLD.txt        # edit if needed
    python msld_prep_cli.py build --config configs/build_example.yaml

    # Or skip straight to build.py with the same --work-dir/--checkpoint --
    # nothing forces a pause in between.

    python msld_prep_cli.py from-smiles --config configs/from_smiles_example.yaml
    python msld_prep_cli.py build --config configs/build_example.yaml \\
        --aligned-dir work.fxa_smiles/aligned --checkpoint work.fxa_smiles/MCS_for_MSLD.txt
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys

import yaml
from rdkit import Chem

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from msldprep import align, analog_prep, grouping, io_mcs_text, pipeline  # noqa: E402
from msldprep.datamodel import LigandRecord  # noqa: E402
from msldprep.errors import MsldPrepError  # noqa: E402

logger = logging.getLogger(__name__)


def _json_list(text: str) -> list:
    try:
        value = json.loads(text)
    except json.JSONDecodeError as exc:
        raise argparse.ArgumentTypeError(f"expected a JSON list, e.g. '[[\"N1\"],[]]', got {text!r}: {exc}") from exc
    if not isinstance(value, list):
        raise argparse.ArgumentTypeError(f"expected a JSON list, got {text!r}")
    return value


def _name_value_pairs_to_dict(pairs: list[str]) -> dict[str, str]:
    """['f_analog=Fc1ccc...', ...] -> {'f_analog': 'Fc1ccc...', ...}.

    Splits on the *first* '=' only, since SMILES strings routinely contain
    '=' themselves (double bonds) but analog names should not.
    """
    result: dict[str, str] = {}
    for pair in pairs:
        if "=" not in pair:
            raise argparse.ArgumentTypeError(f"expected NAME=SMILES, got {pair!r}")
        name, smiles = pair.split("=", 1)
        result[name] = smiles
    return result


def _parse_smiles_file(path: str) -> dict[str, str]:
    """Parse a text file of analog SMILES: one entry per line, `SMILES
    [NAME]` (whitespace-separated; anything after NAME is ignored, so a
    trailing comment/extra column is harmless). NAME is auto-generated as
    "analog<line_number>" (with a warning) if omitted -- every analog still
    needs a name to use as its file stem/ligand identifier. Blank lines
    and lines starting with '#' are skipped.
    """
    result: dict[str, str] = {}
    with open(path) as fp:
        for lineno, raw_line in enumerate(fp, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(maxsplit=1)
            smiles = parts[0]
            if len(parts) > 1:
                name = parts[1].split(maxsplit=1)[0]
            else:
                name = f"analog{lineno}"
                logger.warning("%s line %d: no name given for SMILES %r; using %r.", path, lineno, smiles, name)
            if name in result:
                logger.warning("%s line %d: duplicate analog name %r; overwriting the earlier entry.", path, lineno, name)
            result[name] = smiles
    return result


def _load_reference_from_file(
    path: str, name: str | None, cgenff: str, cgenff_base_prm: str | None
) -> LigandRecord:
    """Load a reference ligand directly from a .mol2 or .sdf/.mol file
    (must already carry a real 3D conformer -- it's the alignment target
    every analog gets tethered to) and type it via the same streamed local
    CGenFF call analogs use (`analog_prep.mol_to_ligand_record`), rather
    than requiring pre-existing .rtf/.prm files.
    """
    ext = os.path.splitext(path)[1].lower()
    if ext == ".mol2":
        mol = Chem.MolFromMol2File(path, removeHs=False, sanitize=True)
    elif ext in (".sdf", ".sd", ".mol"):
        mol = Chem.MolFromMolFile(path, removeHs=False, sanitize=True)
        if mol is not None:
            mol = Chem.AddHs(mol, addCoords=True)
    else:
        raise SystemExit(f"--reference-file must be a .mol2 or .sdf/.mol file, got {path!r}")
    if mol is None:
        raise SystemExit(f"Could not parse a molecule with a 3D conformer from {path}")

    ligand_name = name or os.path.splitext(os.path.basename(path))[0]
    return analog_prep.mol_to_ligand_record(
        mol, ligand_name, cgenff, cgenff_base_prm=cgenff_base_prm, is_reference=True
    )


def _add_config_and_log_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--config", metavar="PATH", help="YAML config file; CLI flags override its values")
    parser.add_argument("--log-level", default=None, help="Python logging level (default: INFO)")


def build_parser() -> tuple[argparse.ArgumentParser, dict[str, argparse.ArgumentParser]]:
    parser = argparse.ArgumentParser(
        description="Build MSLD-ready ligand files (msldprep pipeline).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- step 1: from-files -------------------------------------------------
    from_files = subparsers.add_parser(
        "from-files", help="Step 1: define the common core from ligands that already have mol2/rtf/prm files"
    )
    _add_config_and_log_args(from_files)
    from_files.add_argument("--sysname", help="System name, used only to compute the default --work-dir")
    from_files.add_argument("--base-dir", help="Directory with {name}.mol2/{name}.rtf/{name}.prm")
    from_files.add_argument("--reference", help="Reference ligand name (its mol2/rtf/prm stem in --base-dir)")
    from_files.add_argument(
        "--analogs", nargs="+", metavar="NAME",
        help="Analog ligand names (each one's mol2/rtf/prm stem in --base-dir). Pass a single "
             "'*' (quote it -- an unquoted * is expanded by your shell before Python ever sees "
             "it) to include every {stem}.mol2 in --base-dir that has a matching {stem}.str "
             "(CGenFF's own combined-topology-and-parameter stream file), except --reference "
             "itself; a {stem}.mol2 with no matching .str logs a warning and is skipped.",
    )
    from_files.add_argument(
        "--work-dir", help="Directory for this step's output. Default: work.<sysname>"
    )
    from_files.add_argument(
        "--checkpoint", help="Path for the exported common-core checkpoint. Default: <work-dir>/MCS_for_MSLD.txt"
    )
    from_files.add_argument(
        "--aligned-dir",
        help="Directory for each ligand's aligned sdf/rtf/prm, which `build` reads from later. "
             "Default: <work-dir>/aligned",
    )
    from_files.add_argument(
        "--similarity-threshold", type=float,
        help="Tanimoto/Morgan similarity cutoff for splitting analogs into per-group shared "
             f"cores instead of forcing one global core across all of them. Default: {grouping.DEFAULT_SIMILARITY_THRESHOLD}. "
             "Pass 0 to disable grouping and always compute a single global core.",
    )

    # --- step 1: from-smiles -------------------------------------------------
    from_smiles = subparsers.add_parser(
        "from-smiles",
        help="Step 1: define the common core, building a reference ligand's analogs from SMILES first",
    )
    _add_config_and_log_args(from_smiles)
    from_smiles.add_argument("--sysname", help="System name, used only to compute the default --work-dir")
    from_smiles.add_argument(
        "--reference-dir",
        help="Directory with the reference ligand's already-parameterized mol2/rtf/prm "
             "(mutually exclusive with --reference-file)",
    )
    from_smiles.add_argument(
        "--reference-file",
        help="Path to the reference ligand as a single .mol2 or .sdf/.mol file (must already "
             "carry a real 3D conformer -- every analog gets aligned to it) -- typed via the "
             "same streamed local CGenFF call as the analogs, instead of requiring pre-existing "
             ".rtf/.prm (mutually exclusive with --reference-dir)",
    )
    from_smiles.add_argument(
        "--reference",
        help="Reference ligand name: its mol2/rtf/prm stem in --reference-dir, or (with "
             "--reference-file) just an identifier -- defaults to --reference-file's basename",
    )
    from_smiles.add_argument(
        "--analog-smiles", nargs="+", metavar="NAME=SMILES",
        help="One or more NAME=SMILES pairs, e.g. f_analog=Fc1ccc(Nc2ncnc3ccccc23)cc1. "
             "Prefer the config file for many analogs -- it is a plain YAML mapping there. "
             "Merged with --analog-smiles-file if both are given; overrides it on a name clash.",
    )
    from_smiles.add_argument(
        "--analog-smiles-file",
        help="Path to a text file of analog SMILES, one per line: 'SMILES [NAME]' "
             "(whitespace-separated; '#' comments and blank lines are skipped; NAME is "
             "auto-generated with a warning if omitted). Merged with --analog-smiles.",
    )
    from_smiles.add_argument(
        "--work-dir",
        help="Directory generated analog mol2/rtf/prm files are written into before alignment. "
             "Default: work.<sysname>",
    )
    from_smiles.add_argument(
        "--checkpoint", help="Path for the exported common-core checkpoint. Default: <work-dir>/MCS_for_MSLD.txt"
    )
    from_smiles.add_argument(
        "--aligned-dir",
        help="Directory for each ligand's aligned sdf/rtf/prm, which `build` reads from later. "
             "Default: <work-dir>/aligned",
    )
    from_smiles.add_argument("--cgenff-bin", help="Path to the standalone cgenff binary (default: $MSLDPREP_CGENFF_PATH, else `cgenff` on PATH)")
    from_smiles.add_argument("--cgenff-base-prm", help="Path to CGenFF's base .prm file, passed as `-p charmm` (default: $MSLDPREP_CGENFF_BASE_PRM, else omitted)")
    from_smiles.add_argument(
        "--similarity-threshold", type=float,
        help="Tanimoto/Morgan similarity cutoff for splitting analogs into per-group shared "
             f"cores instead of forcing one global core across all of them. Default: {grouping.DEFAULT_SIMILARITY_THRESHOLD}. "
             "Pass 0 to disable grouping and always compute a single global core.",
    )

    # --- step 2: build ---------------------------------------------------------
    build = subparsers.add_parser(
        "build", help="Step 2: read a common-core definition and write every CHARMM/ALF file"
    )
    _add_config_and_log_args(build)
    build.add_argument("--sysname", help="System name (used for the default --outdir and the .inp filename)")
    build.add_argument(
        "--checkpoint", help="Path to the common-core checkpoint written by from-files/from-smiles (possibly hand-edited)"
    )
    build.add_argument(
        "--aligned-dir",
        help="Directory with every ligand's aligned sdf/rtf/prm, as written by from-files/from-smiles "
             "(ligand names are read from --checkpoint, not repeated here)",
    )
    build.add_argument(
        "--outdir",
        help="Output directory for the MSLD build. Default: build.<sysname>. If --checkpoint holds "
             "more than one CORE ID (see --core-id), every group is built automatically, each into "
             "its own outdir/core<N>/ subdirectory.",
    )
    build.add_argument(
        "--core-id", type=int,
        help="Build only this CORE ID from --checkpoint instead of every group found there "
             "(see from-files/from-smiles's --similarity-threshold for how a checkpoint can end "
             "up with more than one).",
    )
    cgenff_group = build.add_mutually_exclusive_group()
    cgenff_group.add_argument(
        "--cgenff", dest="cgenff", action="store_true", default=None,
        help="CGenFF/ParamChem parameters are being used (default)",
    )
    cgenff_group.add_argument(
        "--no-cgenff", dest="cgenff", action="store_false", help="An explicit (non-CGenFF) force field is being used"
    )
    build.add_argument(
        "--in-frag", type=_json_list, metavar="JSON",
        help='Reflig core atoms to move into each fragment at each site: a JSON list of '
             'per-site lists of atom names, e.g. \'[["N1"], []]\'. Prefer the config file '
             'for this option -- it is a plain YAML nested list there.',
    )
    build.add_argument(
        "--an-core", type=_json_list, metavar="JSON",
        help="Anchor atoms to keep in the core at each site during charge renormalization "
             "(otherwise every site's anchor atom(s) move into the fragment by default): "
             "same JSON-list-of-lists shape as --in-frag.",
    )
    build.add_argument(
        "--prm-search-dir", help="Directory to search for each ligand's .prm file. Default: --aligned-dir"
    )
    build.add_argument(
        "--toppar-search-dir",
        help="Directory to copy toppar/ files from when --cgenff is set (typically wherever "
             "--base-dir/--reference-dir pointed at in step 1). Default: current directory.",
    )

    return parser, {"from-files": from_files, "from-smiles": from_smiles, "build": build}


def _load_config(path: str) -> dict:
    with open(path) as fp:
        config = yaml.safe_load(fp)
    if config is None:
        return {}
    if not isinstance(config, dict):
        raise SystemExit(f"Config file {path} must contain a YAML mapping of option -> value, not a {type(config).__name__}.")
    return config


def _apply_config_defaults(args: argparse.Namespace, config: dict, parser: argparse.ArgumentParser) -> None:
    """Fill in any option left at None (i.e. not given on the command line)
    from `config`. Options actually passed on the command line always win.
    """
    # argparse has no public API to list a parser's option dests; `_actions`
    # is the conventional (if underscore-prefixed) way to introspect this.
    valid_dests = {action.dest for action in parser._actions}
    for key, value in config.items():
        if key not in valid_dests:
            logger.warning("Config key %r is not a recognized option for this subcommand; ignored.", key)
            continue
        if getattr(args, key, None) is None:
            setattr(args, key, value)


def _require(args: argparse.Namespace, fields: list[str]) -> None:
    missing = [f for f in fields if getattr(args, f, None) is None]
    if missing:
        flags = ", ".join(f"--{f.replace('_', '-')}" for f in missing)
        raise SystemExit(f"Missing required option(s): {flags} (set on the command line or in --config)")


def _discover_wildcard_analogs(base_dir: str, reference: str) -> list[str]:
    """Every ligand stem in `base_dir` with both a .mol2 and a matching
    .str file (CGenFF's own combined-topology-and-parameter stream
    output), excluding `reference`. Logs a warning for every .mol2 found
    with no matching .str.
    """
    stems = sorted({
        os.path.splitext(fname)[0] for fname in os.listdir(base_dir) if fname.endswith(".mol2")
    })
    analogs = []
    for stem in stems:
        if stem == reference:
            continue
        if os.path.exists(os.path.join(base_dir, f"{stem}.str")):
            analogs.append(stem)
        else:
            logger.warning("%s.mol2 found in %s with no matching %s.str; skipping.", stem, base_dir, stem)
    return analogs


def _is_wildcard(analogs) -> bool:
    return analogs == "*" or analogs == ["*"]


def _load_analogs(base_dir: str, names: list[str], *, tolerate_failures: bool) -> list:
    """Load each of `names` from `base_dir`. With `tolerate_failures=False`
    (an explicit --analogs list), any load failure raises immediately --
    you asked for exactly these ligands. With `tolerate_failures=True`
    (wildcard-discovered analogs), a ligand that fails to load (a mol2 with
    a matching .str but no matching .rtf/.prm after all, a malformed mol2,
    etc.) is logged as a warning and skipped instead of aborting the whole
    run over one bad ligand -- matching the same "warn and keep going"
    spirit as the .str-existence check `_discover_wildcard_analogs` already
    applies.
    """
    analogs = []
    failed = []
    for name in names:
        try:
            analogs.append(align.load_ligand(
                os.path.join(base_dir, f"{name}.mol2"),
                os.path.join(base_dir, f"{name}.rtf"),
                os.path.join(base_dir, f"{name}.prm"),
                name=name,
            ))
        except (MsldPrepError, OSError) as exc:
            if not tolerate_failures:
                raise
            logger.warning("Could not load analog %s from %s; skipping. (%s)", name, base_dir, exc)
            failed.append(name)
    if failed:
        logger.warning("%d/%d wildcard-discovered analog(s) failed to load and were skipped: %s", len(failed), len(names), failed)
    return analogs


def _run_from_files(args: argparse.Namespace) -> int:
    _require(args, ["sysname", "base_dir", "reference", "analogs"])
    work_dir = args.work_dir or f"work.{args.sysname}"
    checkpoint = args.checkpoint or os.path.join(work_dir, "MCS_for_MSLD.txt")
    aligned_dir = args.aligned_dir or os.path.join(work_dir, "aligned")

    analog_names = args.analogs
    wildcard = _is_wildcard(analog_names)
    if wildcard:
        analog_names = _discover_wildcard_analogs(args.base_dir, args.reference)
        logger.info("Discovered %d analog(s) via wildcard in %s: %s", len(analog_names), args.base_dir, analog_names)

    reference = align.load_ligand(
        os.path.join(args.base_dir, f"{args.reference}.mol2"),
        os.path.join(args.base_dir, f"{args.reference}.rtf"),
        os.path.join(args.base_dir, f"{args.reference}.prm"),
        name=args.reference, is_reference=True,
    )
    analogs = _load_analogs(args.base_dir, analog_names, tolerate_failures=wildcard)
    similarity_threshold = (
        args.similarity_threshold if args.similarity_threshold is not None else grouping.DEFAULT_SIMILARITY_THRESHOLD
    )

    pipeline.define_common_core(
        reference, analogs, mcs_checkpoint_path=checkpoint, aligned_dir=aligned_dir,
        similarity_threshold=similarity_threshold,
    )

    print(f"\nCommon-core definition written to {checkpoint}")
    print(f"Aligned ligand files written to {aligned_dir}/")
    print("Review/edit the checkpoint if needed, then run the `build` command.")
    return 0


def _run_from_smiles(args: argparse.Namespace) -> int:
    _require(args, ["sysname"])
    if not args.reference_file and not (args.reference_dir and args.reference):
        raise SystemExit(
            "Missing required option(s): --reference-file, or both --reference-dir and --reference "
            "(set on the command line or in --config)"
        )
    if not args.analog_smiles and not args.analog_smiles_file:
        raise SystemExit(
            "Missing required option(s): --analog-smiles and/or --analog-smiles-file "
            "(set on the command line or in --config)"
        )

    work_dir = args.work_dir or f"work.{args.sysname}"
    checkpoint = args.checkpoint or os.path.join(work_dir, "MCS_for_MSLD.txt")
    aligned_dir = args.aligned_dir or os.path.join(work_dir, "aligned")
    tools = analog_prep.ExternalTools(
        cgenff=args.cgenff_bin or os.environ.get("MSLDPREP_CGENFF_PATH") or "cgenff",
        cgenff_base_prm=args.cgenff_base_prm or os.environ.get("MSLDPREP_CGENFF_BASE_PRM"),
    )

    analog_smiles = _parse_smiles_file(args.analog_smiles_file) if args.analog_smiles_file else {}
    if args.analog_smiles:  # from the CLI: ["NAME=SMILES", ...] -> {NAME: SMILES}; overrides the file on a name clash
        inline = _name_value_pairs_to_dict(args.analog_smiles) if isinstance(args.analog_smiles, list) else args.analog_smiles
        for name in set(analog_smiles) & set(inline):
            logger.warning("Analog name %r given both in --analog-smiles-file and --analog-smiles; using the latter.", name)
        analog_smiles.update(inline)

    if args.reference_file:
        if args.reference_dir:
            logger.warning("Both --reference-file and --reference-dir given; using --reference-file.")
        reference = _load_reference_from_file(args.reference_file, args.reference, tools.cgenff, tools.cgenff_base_prm)
    else:
        reference = align.load_ligand(
            os.path.join(args.reference_dir, f"{args.reference}.mol2"),
            os.path.join(args.reference_dir, f"{args.reference}.rtf"),
            os.path.join(args.reference_dir, f"{args.reference}.prm"),
            name=args.reference, is_reference=True,
        )

    similarity_threshold = (
        args.similarity_threshold if args.similarity_threshold is not None else grouping.DEFAULT_SIMILARITY_THRESHOLD
    )
    pipeline.define_common_core_from_smiles(
        reference, analog_smiles,
        mcs_checkpoint_path=checkpoint, aligned_dir=aligned_dir, workdir=work_dir, tools=tools,
        similarity_threshold=similarity_threshold,
    )

    print(f"\nCommon-core definition written to {checkpoint}")
    print(f"Aligned ligand files written to {aligned_dir}/")
    print("Review/edit the checkpoint if needed, then run the `build` command.")
    return 0


def _run_build(args: argparse.Namespace) -> int:
    _require(args, ["sysname", "checkpoint", "aligned_dir"])
    outdir = args.outdir or f"build.{args.sysname}"
    prm_search_dir = args.prm_search_dir or args.aligned_dir
    toppar_search_dir = args.toppar_search_dir or "."
    cgenff = True if args.cgenff is None else args.cgenff

    available_core_ids = io_mcs_text.read_core_ids(args.checkpoint)
    if args.core_id is not None:
        if args.core_id not in available_core_ids:
            raise SystemExit(
                f"--core-id {args.core_id} not found in {args.checkpoint} "
                f"(available: {available_core_ids})"
            )
        core_ids = [args.core_id]
    else:
        core_ids = available_core_ids

    multi = len(core_ids) > 1
    for core_id in core_ids:
        reference_name, analog_names = io_mcs_text.read_ligand_names(args.checkpoint, core_id=core_id)
        reference, analogs = align.load_ligand_set_from_sdf(
            pipeline.aligned_sdf_path(args.aligned_dir), args.aligned_dir, reference_name, analog_names,
        )
        group_outdir = os.path.join(outdir, f"core{core_id}") if multi else outdir

        pipeline.build_from_core_definition(
            reference, analogs,
            sysname=args.sysname, outdir=group_outdir, cgenff=cgenff,
            mcs_checkpoint_path=args.checkpoint, core_id=core_id,
            in_frag=args.in_frag, an_core=args.an_core,
            prm_search_dir=prm_search_dir, toppar_search_dir=toppar_search_dir,
        )
        print(f"\nMSLD build for CORE ID {core_id} (reference {reference_name}) written to {group_outdir}/")

    return 0


def main(argv: list[str] | None = None) -> int:
    parser, subparsers_by_name = build_parser()
    args = parser.parse_args(argv)

    if args.config:
        _apply_config_defaults(args, _load_config(args.config), subparsers_by_name[args.command])
    if args.log_level is None:
        args.log_level = "INFO"

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(levelname)s %(name)s: %(message)s",
    )

    if args.command == "from-files":
        return _run_from_files(args)
    if args.command == "from-smiles":
        return _run_from_smiles(args)
    if args.command == "build":
        return _run_build(args)
    raise AssertionError(args.command)  # argparse's `required=True` on the subparsers makes this unreachable


if __name__ == "__main__":
    sys.exit(main())
