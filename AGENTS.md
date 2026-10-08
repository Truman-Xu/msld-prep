# Instructions for AI agents

msld-prep (Python package `msldprep`) prepares small-molecule ligand series for multisite
lambda dynamics (MSLD) in CHARMM: it loads CGenFF-parameterized ligands, aligns them in 3D,
finds the shared core and the per-site fragments, renormalizes charges, and writes the CHARMM
and ALF input files. It is an RDKit-based, in-memory rewrite of the older `msld-py-prep`
scripts, and it builds on crimm's RTF and PRM parsers.

**Before doing anything, read [`docs/dev/AGENT_HANDOFF.md`](docs/dev/AGENT_HANDOFF.md)**
sections 0 to 5: protocol, the owner's standing instructions, environment facts, traps,
current state and claims. Standards are in
[`docs/dev/BEST_PRACTICES.md`](docs/dev/BEST_PRACTICES.md); developer rules, which bind you
too, are in [`CONTRIBUTING.md`](CONTRIBUTING.md). When you finish, log what you did in the
handoff. There is no `ROADMAP.md` or `NOTES.md` yet; do not create them until the owner asks.

The rules that matter most, in case you read nothing else:

- Do not commit or push unless the owner asks. That permission does not carry over to a later
  session.
- Sign commits, never files: a `Co-Authored-By:` trailer naming your exact model; no author
  names, banners or session notes in any file.
- `/home/ziqiaoxu/msld-prep-tools` is a symlink to this directory and is published to other
  cluster users. Never remove or repoint it.
- Do not edit anything in `/home/ziqiaoxu/crimm`, `/home/ziqiaoxu/crimm-ml` or
  `/home/ziqiaoxu/msld-py-prep`. They are separate projects. This project only reads them. If
  a change is needed there, say so in your handoff entry and let the owner decide.
- Default test run: `pytest -p no:cacheprovider -q` (about 15 s). CGenFF-backed tests skip
  unless `MSLDPREP_CGENFF_PATH` is set, and OpenBabel-backed ones unless `MSLDPREP_OBABEL_PATH`
  is (or `obabel` is on `PATH`); a skipped test is not a pass. Python and environment
  details are in the handoff, section 3. The env to use, `crimm-dev`, is not on `PATH`.
- Library code logs or raises; it does not `print`. Raise a subclass of `MsldPrepError`
  (`msldprep/errors.py`). Do not use `assert` to validate input.
- RDKit, numpy and crimm are core dependencies of the pipeline. Anything heavier, or only
  needed by one writer or tool, is imported inside the function that uses it.
- `msldprep/sybyl_typing_data.py` is generated. Change `scripts/generate_sybyl_typing_data.py`
  or its inputs in `data/`, then regenerate; never edit the output by hand.
- Output file formats (rtf, prm, pdb, the `lpsites.inp` and ALF scripts, `MCS_for_MSLD.txt`)
  are contracts with CHARMM and downstream scripts. A change to what is written is its own
  change, with a test, and is called out in the handoff.
- A bug you find but do not fix gets an `xfail(strict=True)` test with the reason and a
  GitHub issue. Remove the marker in the change that fixes it.
- One concern per change. A fix does not ride along with a rename.
- A `TODO` in code carries an issue number.
- Run output (`work.*/`, `build.*/`) is git-ignored; never commit it.
- Do not run CGenFF-backed or long jobs on the shared login node without telling the owner.
- Personal notes of the owner live in `CLAUDE.local.md`, which is not in the repository.
