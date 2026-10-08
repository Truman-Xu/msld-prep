# Good practices for msld-prep

Standards for keeping msldprep correct and maintainable. They are targets, not a description
of the current state; [AGENT_HANDOFF.md](AGENT_HANDOFF.md) has the state. Where something is a
suggestion to decide on and not a clear win, it says so.

## The short list

1. Every bug fix comes with a test that fails without it.
2. Scientific output (charges, core and site assignment, parameters) is checked against
   reference files, ideally ones produced by the legacy tools or by CHARMM.
3. Files the pipeline writes are contracts; changing them is a separate, described change.
4. Library code logs or raises; it does not `print`. Errors derive from `MsldPrepError`.
5. Run output (`work.*/`, `build.*/`) is never committed. The one generated module the package
   needs (`sybyl_typing_data.py`) is regenerated, not hand-edited.
6. A `TODO` in the code points to an issue.
7. No direct pushes to `main` once a remote exists; changes arrive by pull request.

## Code

- Type-hint public functions. State **units** for physical quantities: Å, e (elementary
  charge), degrees. Say in the docstring whether a function modifies its argument.
- Never use `assert` for input or state validation. Raise the specific `MsldPrepError`
  subclass that fits (`AlignmentError`, `DecompositionError`, `CRNError`, `WriterError`,
  `PrmError`, `ExternalToolError`).
- Never catch bare `Exception` without re-raising or logging what was caught. Catch the narrow
  set RDKit raises, as `align._EMBEDDING_EXCEPTIONS` does.
- Log through `logging.getLogger(__name__)`. The package ships a `NullHandler` only; handlers
  belong to the application.
- No module-level mutable state. Pass configuration in (`SystemConfig`), do not read it from
  globals.
- Do not write to the current working directory by default. Take an explicit directory.
- Give RDKit properties that carry meaning between stages a named constant (for example
  `_SDF_ATOM_NAMES_PROP`, the `CGENFFAtomType` property), not repeated string literals.
- Choose a linter and formatter (suggestion: ruff, run through pre-commit) before the code
  grows. Reformat the codebase once, in one commit with no other changes.

## Science-specific code

- Charge renormalization must conserve the total charge of each ligand to a stated tolerance.
  Test that, per ligand and per fragment. Write the tolerance down.
- Anything that embeds or aligns in 3D (RDKit `ConstrainedEmbed`) is not fully deterministic.
  Pass a seed where RDKit accepts one, and test properties of the result (core RMSD below a
  bound, atom correspondence) rather than exact coordinates.
- Atom names are part of the output contract. Test name translation between the reference and
  each analog (`compute_atom_name_translation`), including collisions.
- Keep CGenFF atom-type handling in one place. Types travel through RDKit as isotope labels
  and properties; a second encoding would break the grouping silently.
- Compare floats with tolerances chosen for the quantity.
- Cite the method next to the implementation: MSLD, ALF, the charge renormalization scheme,
  CGenFF. Record the CGenFF and force-field versions used in the output where the format has
  room for a comment.

## Testing

- **Layers.** Fast offline tests on the small ChEMBL fixtures run on every change. Tests that
  need CGenFF, the legacy `msld-py-prep` checkout or OpenBabel skip themselves with a reason
  when the tool is absent, and say so; a skipped test is not a pass.
- **Fixtures.** Keep a few small ligand series in the repository covering the hard cases:
  charged groups, multiple sites, an analog that falls in a separate core group, chirality,
  aromatic ring reload, lone pairs.
- **Golden files.** Check written rtf, prm, pdb and inp files against stored references.
- **Invariants.** Cheap properties that hold for any input: total charge is conserved,
  every fragment atom has exactly one owner, write-then-read returns the same records.
- **Regression.** Every bug fix adds the test that would have caught it.
- **Notebooks.** `align.ipynb` is committed, with its outputs, at the owner's request
  (2026-10-08). It is 8.9 MB; do not add further executed notebooks without asking, and
  consider stripping outputs from it later (suggestion: `nbstripout`).

## Scientific validity and reproducibility

- Record provenance in output: msldprep version, CGenFF version, the config used.
- When a bug could have produced wrong parameters or charges, say so plainly in the changelog
  and describe how users can check systems they already built. In scientific software this
  matters more than the fix.

## Packaging and release

- One version source, in `pyproject.toml`. Add `msldprep.__version__` when a release is cut.
- Keep a `CHANGELOG.md` in the Keep a Changelog layout once there is a release.
- State the supported Python and RDKit range. RDKit output (aromaticity perception, mol2
  reading) changes between releases; test the oldest and newest you support.
- Check that data files needed at run time are inside the package; `data/` is a build input
  for `sybyl_typing_data.py`, not a run-time dependency.
- Licence: GPLv3 (owner decision, 2026-10-08). Check that it is compatible with the licence of
  the code carried over from `msld-py-prep`.

## Project management

- Track work as GitHub issues with acceptance criteria. Branches and pull requests link to them.
- Short-lived branches, deleted on merge. Squash-merge with a message that explains why.
- Commit and push at the end of each session so nothing lives only on the cluster. Only when
  the owner asks.
- `CONTRIBUTING.md`, issue forms, a pull-request template and `CODEOWNERS` exist. The remote
  is https://github.com/Truman-Xu/msld-prep (empty until the first push).

## Working with AI coding tools

- Generated code gets the same review as any other, with extra attention to numerical
  details: atom order, index base, units, sign conventions.
- Ask for tests in the same change, and read the tests first.
- Keep work orders for agents as GitHub issues, not loose Markdown files at the repo root.
- Do not commit agent scratch files or local tool settings; add them to `.gitignore`, or for
  files personal to one developer, to `.git/info/exclude`.
- The root files for agents are `AGENTS.md` (the instructions) and `CLAUDE.md` (one line,
  `@AGENTS.md`).
