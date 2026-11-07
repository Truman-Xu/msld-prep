# AGENT_HANDOFF

Shared memory for AI agents working on msld-prep (Python package `msldprep`). Any model, any tool. Written for agents
first; humans read it second.

```yaml
schema: agent-handoff/1
repo: github.com/Truman-Xu/msld-prep   # public; first push 2026-10-08
default_branch: main
last_updated: 2026-11-07
log: open issue labelled agent-log   # gh issue list --label agent-log --state open
```

## 0. Protocol

Read this section fully. Skim the rest for what your task touches.

This file holds what changes rarely: protocol, orientation, the owner's standing
instructions, the environment and gotchas. What changes every session (STATE, CLAIMS and
LOG) is in the open GitHub issue labelled `agent-log` (section 5), so that recording a
session needs no commit.

### On arrival

1. Read sections 1 to 4. They are short and change rarely.
2. Read the `agent-log` issue: its description (CLAIMS and STATE) and its newest three
   comments (LOG). `N=$(gh issue list --label agent-log --state open --json number -q
   '.[0].number')`, then `gh issue view $N` and `gh issue view $N --comments`.
3. Run `git status --short` and `git log --oneline -5`. If they disagree with STATE, trust
   git and correct STATE.
4. Check CLAIMS before editing a file. If another agent holds it, do not edit it.

### On departure

1. Post one LOG entry as a comment on the `agent-log` issue (template in the issue
   description): `gh issue comment $N --body-file <file>`.
2. Update STATE so it is true now. Re-read the description immediately before, edit it in
   place, and write it back with `gh issue edit $N --body-file <file>`. Do not append.
3. Release your CLAIMS in the same edit.
4. Move anything you learned that will still be true next month into FACTS or GOTCHAS of
   this file. That is a tracked change; commit it with your work.

Without network access or `gh`, write the entry and the new STATE to local files and post
them in your next session that has both.

When the issue grows long (about 50 comments) or a milestone closes, open the next one
(label `agent-log`) with the current STATE, link it from the old one, then close the old one
and move the pin.

### Writing rules

- **Tag every claim with how you know it.** Use exactly one of:
  - `[ran]` you executed it and saw the result this session
  - `[read]` you read the code or file; did not execute
  - `[told]` the owner said so
  - `[inferred]` your conclusion; say from what
- **Date everything** as `YYYY-MM-DD`. Never write "today", "recently" or "currently".
- **Cite locations** as `path:line` or `path::symbol`. Prefer the symbol; line numbers rot.
- **State what you did not check.** A missing check written down is worth more than a
  confident guess.
- **One fact, one place.** Link to another document instead of copying from it.
- **Do not edit or delete another agent's entry.** If it is wrong, say so in your own
  comment: `SUPERSEDED YYYY-MM-DD by <agent>: <why>`, with a link to the comment.
- **Identify yourself** as `<model>/<tool>/<short task name>`, for example
  `claude-sonnet/claude-code/repo-setup`. You have no persistent identity; the task name is
  what lets others tell sessions apart.
- **The log is public**, like this file. Comments appear under the GitHub account the agent
  runs as, usually a person's, so the `<model>/<tool>/<task>` heading of your comment is
  what tells sessions apart. Nothing about a particular machine or account in the log;
  that goes in the owner's `CLAUDE.local.md`.
- Keep entries terse. No preamble, no summaries of summaries.

### What goes where

| Kind of information | Place |
| --- | --- |
| How the code is built, conventions, architecture | `docs/dev/NOTES.md` (**not created yet**; the owner will ask for it) |
| What work is planned, its order, its status | `docs/dev/ROADMAP.md` (**not created yet**; the owner will ask for it) |
| Active issues and planned features (not a roadmap; no order yet) | `docs/dev/ISSUES_AND_PLANS.md` |
| Design decisions for the re-design, with the owner's verdicts and reasons | `docs/dev/DESIGN_NOTES.md`, kept identical in content to `design-notes.md` in the owner's vault `/home/ziqiaoxu/obsidian/msld-prep/`; edit both in the same session |
| Standards the work should meet | `docs/dev/BEST_PRACTICES.md` |
| Reference papers, as PDFs | `papers/` (git-ignored; never commit, the repository is public). If a paper can't be read online, tell the owner, who will download it `[told 2026-10-30]` |
| Anything about crimm, which msldprep depends on | `/home/ziqiaoxu/crimm/docs/dev/AGENT_HANDOFF.md` (read-only for you) |
| Who is doing what now, current state, what was tried, what failed, decisions | the `agent-log` issue (section 5) |
| Questions an agent could not settle | an issue labelled `question` |
| Protocol, standing instructions, environment facts, traps | this file |

Until NOTES and ROADMAP exist, put structure facts under Orientation or Gotchas here, and open
questions in issues labelled `question`. Do not create the two files on your own.

## 1. Orientation

- msld-prep (package `msldprep`) prepares ligand series for multisite lambda dynamics (MSLD) in CHARMM. It is an
  in-memory, RDKit-based rewrite of the older scripts in `/home/ziqiaoxu/msld-py-prep`.
  `[read msldprep/__init__.py, tests/conftest.py]`
- GPLv3. `[told 2026-10-08]` The repository is github.com/Truman-Xu/msld-prep. `[told 2026-10-08]`
- Owner: Truman Xu (`Truman-Xu`). `[told 2026-10-08]` for the name; git `user.name` is
  `Truman-Xu` `[ran 2026-10-08]`.
- Package map `[read 2026-10-08]`, in pipeline order:
  - `align.py`: load a ligand, type it from its rtf, align analogs onto the reference with
    RDKit, find the maximum common substructure.
  - `analog_prep.py`, `sybyl_typer.py`, `sybyl_typing_data.py`: build a parameterized ligand
    from SMILES or an RDKit `Mol`, using in-memory Sybyl typing and a streamed local CGenFF
    call. `sybyl_typing_data.py` is generated by `scripts/generate_sybyl_typing_data.py`.
  - `decompose.py`, `grouping.py`: core, anchor, site and fragment decomposition; split of
    dissimilar analogs into separate core groups.
  - `crn.py`: charge renormalization.
  - `rtf_records.py`, `prm_records.py`, `datamodel.py`: parsing and records.
  - `io_mcs_text.py`: the reviewable `MCS_for_MSLD.txt` checkpoint.
  - `writers/`: CHARMM and ALF output (rtf, prm, pdb, `lpsites.inp`, control and fragment
    json).
  - `gui/`: optional alignment web GUI (`python -m msldprep.gui`): `service.py` (in-memory
    `AlignmentSession`, no web types), `app.py` (FastAPI, token-protected), `static/` (page and
    bundled 3Dmol.js 2.4.2). Nothing in the pipeline imports it. `[read 2026-10-14]`
  - `pipeline.py`: the two explicit steps `define_common_core` and
    `build_from_core_definition`, plus thin `run_pipeline*` wrappers.
  - Entry points at the repo root: `msld_prep_cli.py` (`from-files`, `from-smiles`, `build`),
    `run_msld_prep.py`, `run_msld_prep_from_smiles.py`; configs in `configs/`.
- **Dependency on crimm.** The dependency is the latest PyPI release of crimm for now; the
  version is pinned after the crimm overhaul is finished. `[told 2026-10-08]`
  This project is listed in crimm's handoff (Orientation, "What goes where") with the names it
  imports. `msldprep/rtf_records.py` imports `crimm.IO.RTFParser.RTFParser`;
  `msldprep/prm_records.py` imports from `crimm.IO.PRMParser`; `tests/test_prm_records.py`
  imports `angle_par`, `bond_par`, `nonbond_par`, `ub_par`, `nonbond14_par` from it.
  `[read, grep 2026-10-08]` crimm is a separate project. This repo only reads it; if crimm
  has to change, say so in your entry and let the owner decide. The owner plans to log this
  project in crimm's handoff later `[told 2026-10-08]`; do not do it unprompted.
- Other sibling projects: `/home/ziqiaoxu/crimm-ml` (not used here), `/home/ziqiaoxu/msld-py-prep`
  (the legacy implementation; one test imports its `msld_prm.py`; read-only for you).

## 2. Standing instructions from the owner

Follow these unless the owner says otherwise in your own session.

| Date | Instruction | Tag |
| --- | --- | --- |
| 2026-11-07 | Agent activity (state, claims, session log) is not kept in a tracked file. It goes in the open `agent-log` issue; open questions are issues labelled `question` (section 5). | `[told]` |
| 2026-11-02 | Commit and push the day's work (this permission covers that push only). The owner will ask someone else to follow the GUI README (#8 item 4) and bring back their feedback. | `[told]` |
| 2026-11-02 | **#16:** add hydrogens to a mol2 reference that lacks them, with a warning (done, `e938bda`). **#17:** pin pandas below 3 for the legacy tests (done, `fb7d4b2`). **#19:** read what RDKit refuses through OpenBabel, the owner's idea (done, `a3be5d6`). **#18** first (done, `b6f60bf`). **(p):** one source file per run, the owner's direction. | `[told]` |
| 2026-10-30 | Update the handoff and agent memory, file issues for (j) to (n) (done, #20 to #24), commit and push (this permission covers that push only). | `[told]` |
| 2026-10-30 | **Second flag level in (j) adopted:** above 0.01 e or abs(dipole change) above 0.5 D, as well as the flag above 0.001 e; provisional, to be re-checked in lesson 6; neither stops a run by default. Not implemented (#20). | `[told]` |
| 2026-10-30 | **(n) confirmed:** the exact legacy rule is the default, single-point anchors out included. **(o) verdict:** keep the connected core for now, revisit in lesson 6. Recorded in `docs/dev/DESIGN_NOTES.md`; not implemented. | `[told]` |
| 2026-10-30 | crimm's hotfix for `RDKConverter` is pushed: update the docs here (done). It is `1816fbb` on crimm's branch `fix/59-rdkit-plain-residue`, Refs crimm #59, replacing `8075b75`. **Pushed is not published:** no PR, not on crimm's `master`, not on PyPI (checked the same day). Keep the `PYTHONPATH` override until it is. | `[told]`, `[ran]` |
| 2026-10-30 | **Design verdicts (j) to (n):** (j) report Q-RMSD per ligand, the largest atom change and the spread at each core position, in e; flag above 0.001 e; no stop by default, optional config limit; the percent check goes. (k) cut-ring part closed by (e); the rest is a report, the core is not optimized automatically. (l) deduplicate, average per substituent, identity by fragment graph with attachment points. (m) spread evenly in one exact pass. (n) keep the legacy behavior, anchors in: a site attached at two or more points keeps them in the core, a single-point anchor still moves out (`msld_crn.py:250-254`); there must be an easy option or function to move them out. Details in `docs/dev/DESIGN_NOTES.md`. **Not implemented yet;** issues #20 (j), #21 (k), #22 (l), #23 (m), #24 (n). | `[told]` |
| 2026-10-29 | **Design verdicts (g), (h):** (g) adopted and extended: flag shrinking cores and one-substituent sites, and design an algorithm that proposes several common cores for a diverse pool, each member covered by its core above a percent-of-atoms cutoff (today's `grouping.group_ligands_by_similarity` uses fingerprint similarity instead); details decided the same day: coverage on heavy atoms, cutoff 0.7 by default or a requested number of cores, analogs may be in several proposals with an option for unique ones (each analog in its best-covering core), search by growing; (h) adopted: no `input()`/`quit()` in library code, the GUI may prompt, the script path follows (b) (stop with a report and a typed error, decision as code or config). Not implemented yet. The notes (e) to (i) were committed and pushed at the owner's request (that permission covered that push only). | `[told]` |
| 2026-10-29 | **Design verdicts (e), (f), (i):** (e) adopted: the legacy neighbor-type rule defines core atoms, implemented on `FindMCS` with pattern labels (prototype `scripts/pattern_mcs.py` in the owner's vault); (f) adopted: atoms more than 1.1 Å off in an input pose stay in the core and the pose is flagged; (i) settled by (e): a core heavy atom's hydrogens are core. Details in `docs/dev/DESIGN_NOTES.md`. **Implemented 2026-10-29** (`fd6cbaa` (e), (i), #15; `b003554` (f)). **Sites are legacy's:** one connected non-core piece per site, so separate substituents on neighboring core atoms are separate sites (confirmed by the owner the same day). | `[told]` |
| 2026-10-28 | **#1 is a flag, not a filter:** with `check_protonation` on, the alignment is unchanged and the differing core atoms are reported (not excluded from the core, not a refused analog). Done, `248b646`. | `[told]` |
| 2026-10-28 | **#10 as written:** only names an atom already has are kept (a mol2 reference); analog core atoms do not inherit the reference's names. Done, `25d4299`. **#15 stays parked** until the re-design settles the core-hydrogen rule. | `[told]` |
| 2026-10-28 | Make the commits (#1, #10, #16 test) and push; later, commit and push the #17 test and the lesson-3 docs (each permission covered that push only). Then lesson 3. | `[told]` |
| 2026-10-21 | **Design verdicts:** (a) settings come from a config file; (b) the script path is a session object run top to bottom, corrections as code or config, stops only for decisions, optional checkpoint, one implementation behind script, GUI and CLI: **adopted as the design**; (d) atom renames live in the data, not in a `translation.txt`. Details in `docs/dev/DESIGN_NOTES.md`. | `[told]` |
| 2026-10-21 | **Keep the design notes in two places:** `docs/dev/DESIGN_NOTES.md` and the vault's `design-notes.md`, same content. | `[told]` |
| 2026-10-21 | Exclude `zoe_msld_prep.ipynb` from git (done, `.git/info/exclude`). The owner's PI plans to publish it in a workshop for the legacy code on GitHub; until then treat it as unpublished (G12). | `[told]` |
| 2026-10-21 | Commit and push the session's work (this permission covers that push only). | `[told]` |
| 2026-10-21 | **#14 (faster core search) is postponed** until the whole pipeline is re-designed and re-implemented: no optimization work before then. | `[told]` |
| 2026-10-21 | **#1 means tolerance**: the structure search must treat different protonation states as the same substructure (done in alignment). Next session: an option, **off by default**, to enforce matching protonation states, because a different state is usually an intended alchemical change. | `[told]` |
| 2026-10-21 | **#2 is limited to the structure alignment code**; the owner is satisfied with it there (closed). Chirality through the whole pipeline is tested after the re-design. | `[told]` |
| 2026-10-21 | Write the GUI README (#8); done as `msldprep/gui/README.md`. **#10 (keep original atom names, as an option when saving mol2) is for the next session.** | `[told]` |
| 2026-10-21 | Running the full suite on the login node is fine. `origin` switched to SSH at the owner's request. | `[told]` |
| 2026-10-21 | **Learning plan started (Q13).** Lessons and design notes live in the owner's Obsidian vault, `/home/ziqiaoxu/obsidian/msld-prep/` (`learning-plan.md`, `lesson-NN-*.md`, `design-notes.md`), not in this repository. Lessons may run the legacy code on public examples. Design verdicts so far: GUI is the main checkpoint, but a script path is also wanted, without the legacy "run twice" flow; combine machine-read output files; keep human-read files to those needing judgement; `translation.txt` should not be needed. | `[told]` |
| 2026-10-08 | This is a separate repository, set up and configured the way crimm is: `AGENTS.md`, `CLAUDE.md`, `CONTRIBUTING.md`, handoff and best-practice docs, issue forms, `CODEOWNERS`, `environment-dev.yml`, `pyproject.toml`. | `[told]` |
| 2026-10-08 | Do not touch anything in the crimm project. | `[told]` |
| 2026-10-08 | Do not write `ROADMAP.md` or `NOTES.md` until the owner asks. | `[told]` |
| 2026-10-08 | The owner will later have an agent log this project in crimm's handoff docs so other agents know it exists. **Done 2026-10-08** by claude-sonnet/claude-code/msld-prep-intro (entry in crimm's `docs/dev/AGENT_HANDOFF.md`, written at the owner's request; crimm's agents may now name this project). | `[told]` |
| 2026-10-08 | SUPERSEDED in part 2026-10-30: the fix is now `1816fbb`, pushed on `fix/59-rdkit-plain-residue`, and the branch `RDKConverter-Het` is deleted; the reminder for when it is published still stands. Old text: The `RDKConverter-Het` hotfix in crimm (`8075b75`) will not be published anytime soon. Crimm's agents are to remind the owner that this project needs it when it is published. | `[told]` |
| 2026-10-08 | Commit the handoff edit, push all commits, turn `TODO.md` into GitHub issues, update the dev docs. This permission covers that push only; it does not carry over to later sessions. | `[told]` |
| 2026-10-08 | Licence is GPLv3. | `[told]` |
| 2026-10-08 | The project is named `msld-prep`; its remote is `https://github.com/Truman-Xu/msld-prep.git`. The importable package stays `msldprep` `[inferred: a hyphen cannot be in an import name]`. | `[told]` |
| 2026-10-08 | `build.*/` and `work.*/` are never committed and are in `.gitignore`. | `[told]` |
| 2026-10-08 | **Do not publish the structures of the analogs.** The lab's analog series is unpublished. Never commit, push, or paste into docs, tests, issues, PRs or commit messages any file or output that contains those structures: `align.ipynb`, `analogues.txt`, `sil.*`, `cellz*.yaml`, or anything built from them. The owner may name more. | `[told]` |
| 2026-10-08 | `align.ipynb` is not committed. It was removed from the local commits, and the history rewritten, before anything was pushed. The first commit is a snapshot. | `[told]` |
| 2026-10-08 | `examples/` is committed. | `[told]` |
| 2026-10-08 | The directory is `/home/ziqiaoxu/msld-prep`, and `/home/ziqiaoxu/msld-prep-tools` stays as a symbolic link to it, because the old path was published to other users of the cluster and is being tested by the lab. Never remove or repoint that link. | `[told]` |
| 2026-10-08 | Keep `crimm-dev` (which carries the local editable crimm, not the PyPI release) as the working env for now. | `[told]` |
| 2026-10-08 | Install paths of CGenFF and OpenBabel come from environment variables, never from source or example configs. | `[told]` |
| 2026-10-08 | Use the latest PyPI release of crimm for now; pin a version after the crimm overhaul. | `[told]` |
| 2026-10-08 | The working env is `crimm-dev` if it satisfies every dependency of this project. It does (see section 3). | `[told]` |
| 2026-10-12 | Project goal: a streamlined pipeline with only the inputs and outputs MSLD needs, in-memory objects for intermediate stages (fewer saved files), and a robust, user-friendly web GUI with visualization tailored to this process. | `[told]` |
| 2026-10-12 | Working item 1: a web GUI for **only the alignment step** that `align.ipynb` runs. Runs on the cluster, reached by SSH tunnel. Stack: FastAPI backend plus small JS front end with 3Dmol.js; the later full app reuses the API. Input: reference with 3D coordinates as **mol2 or sdf** (no CRD, so no crimm hotfix) plus analog SMILES. Output: aligned sdf and one Sybyl-typed mol2 per analog. No CGenFF call in the GUI. Plan: `/home/ziqiaoxu/.claude/plans/thank-you-for-the-valiant-sonnet.md` (milestones M0 to M5); not started. | `[told]` |
| 2026-10-16 | **Scope of the web GUI: structure alignment with MCS only, for now.** Also: (a) for a truly symmetric structure any match is fine; for a decorated ring the analog's substituents are oriented by the reference's substituents (done, #11), and **when protein pocket structures are brought in, the ring flip is to be decided by the pocket's vacant space: raise this then**; (b) stereochemistry: annotate the reference's stereocentres and map them onto the analogs; an unspecified centre outside the reference-backed core is chosen by energy and flagged; a specified configuration opposite to the reference is kept and flagged (done in alignment, #2); (c) a faster MCS for large molecules, starting from the owner's junction-tree idea (#14); **updated at the day's end: the owner decides on Monday 2026-10-19 whether to work on it, so ask then rather than start it**; (d) the Uncharger option is removed (done, #12); (e) `test_structures/` holds the owner's test structures and is git-ignored, never committed; (f) #5 to #7 closed as built. | `[told]` |
| 2026-10-16 | **The web GUI uses American English** in all text it shows (page, notes and messages from `msldprep/gui/`): "stereocenter", "center", "color", "behavior". Applied the same day. | `[told]` |
| 2026-10-16 | Commit and push all of the day's work (this permission covers that push only). | `[told]` |
| 2026-10-16 | **Push to GitHub over SSH** from now on: `git push git@github.com:Truman-Xu/msld-prep.git main`, or `origin` once it points there (`git remote set-url origin git@github.com:Truman-Xu/msld-prep.git`, offered, not done). Pushing still needs the owner's say-so each session. See section 3. | `[told]` |
| 2026-10-14 | **Reminder for the next session (DONE 2026-10-16, #12): remove the GUI option "Remove charges from the reference"** (RDKit `Uncharger`). It was a workaround for #9; since that fix charges do not matter for the match or the alignment, so it serves no purpose. Tracked as #12 (F-13). Say so at the start of the next session; the owner's own words were "set up a reminder for me". Also saved in agent memory. | `[told]` |
| 2026-10-14 | The owner reports that different protonation states of the same ring (`c(=O)[nH]c` against `c(=O)nc`) make RDKit's substructure match treat structures as different, hurting alignment quality, and wants it resolved first: relax the match on protonation states, or log the states, match, align and put the protons back. Local files for this: `sil.mol2` (reference) and `panel_a.smi` (analogs); both belong to the unpublished series (G12). | `[told]` |
| 2026-10-14 | The mol2 export's renaming of every atom and its missing charges are **deliberate**, because CGenFF only accepts unique short names and computes its own charges. Document it (done: README, `ISSUES_AND_PLANS.md` \"Deliberate behaviors\"), and plan an option to keep the original names (#10, F-11). Rewrite #4 and #8 to match the work (done). Push the commits (done). | `[told]` |
| 2026-10-14 | **Every commit that finishes an issue ends its message with `Closes #N`**, on its own line above the `Co-Authored-By` trailer; a commit that only works toward an issue says `Part of #N`. Written into `CONTRIBUTING.md`. Apply it from the next commit on. | `[told]` |
| 2026-10-14 | Next GUI requests: save the aligned structures as mol2 (the notebook's Sybyl typing and mol2 writer), and an option "one structure per file" for individual sdf or mol2 files; file a GitHub issue for the unfixed charged-reference bug (done: #9); make a separate commit for each issue worked on. The owner reports that all visualization features work. Pushed the four earlier commits at the owner's request; the commits made after that request are local. | `[told]` |
| 2026-10-14 | Three requests for the GUI and the docs: (1) a GUI option to remove charges from the reference with `Uncharger`, as the notebook does; (2) a toggle to color the reference by element, gray when off; (3) a dev document of active issues and planned features, **not a roadmap** (the roadmap waits until the owner understands the original pipeline), whose first issue is matching the common structure regardless of atom charge (a charged nitrogen gives "molecule doesn't match the core"). Done as `docs/dev/ISSUES_AND_PLANS.md`; `ROADMAP.md` is still not created. | `[told]` |
| 2026-10-14 | **Item 1 scope narrowed, and a demo is due soon.** The alignment GUI needs only the constrained embed of each analog onto the reference: no atom typing (no Sybyl or CGenFF) and no common core across the series. So no per-analog mol2 export; the download is an sdf (reference first, then analogs). The learning plan (item 2, Q13) is postponed by the owner until after the demo. | `[told]` |
| 2026-10-12 | Working item 2: the owner is not satisfied with the current rewrite design but has not run MSLD and does not yet understand the original implementation. The next deliverable is a **learning plan** for the legacy `msld-py-prep` (not a software roadmap): detailed explanation of each module, charge renormalization, and the purpose of each output file. Design decisions and old-versus-new tests with real MSLD runs come after, with colleagues helping on CHARMM. See Q13. | `[told]` |

Working agreements carried over from crimm (not owner instructions for this project; change
them by logging why) `[inferred: same owner, same setup]`:

- Commit only when asked; never push unless asked. Permission for one session does not carry
  over to the next.
- Sign commits, never files.
- One concern per change. A fix does not ride along with a rename.
- A bug you find but do not fix gets a test marked `xfail(strict=True)` with the reason.

## 3. Environment facts

Machine: the owner's cluster (Rocky Linux 8). Hostname `gollum` `[ran 2026-10-08]`. Paths are
specific to it.

| Fact | Tag |
| --- | --- |
| `python3` on `PATH` is 3.6.8 and is not usable. Use `/home/ziqiaoxu/.conda/envs/crimm-dev/bin/python` (Python 3.12.15, RDKit 2026.09.1, NumPy 2.5.3, pytest 9.1.1, OpenBabel 3.2.1, scipy 1.18.1, PyYAML 6.0.3). | `[ran 2026-10-08]` |
| The env was `crimm-dev-new` when this repo was set up; by the end of the session it had been swapped into `crimm-dev` (`crimm-dev-new/bin/python` no longer exists, same package versions). Re-check with `ls ~/.conda/envs` before relying on either name. | `[ran 2026-10-08]` |
| In that env `crimm` is an editable install that resolves to `/home/ziqiaoxu/crimm/crimm/__init__.py`. | `[ran 2026-10-08]` |
| `torch210` also has crimm, RDKit 2025.09.4 and PyYAML, but no OpenBabel. The suite was not run there. | `[ran 2026-10-08]` imports only |
| `crimm-dev` has every import the package and tests use: rdkit, numpy, yaml, pytest, pandas, openbabel, scipy, crimm, `crimm.IO.RTFParser`, `crimm.IO.PRMParser`. So no `environment-dev.yml` exists in this repo; I deleted the unbuilt one I had written. | `[ran 2026-10-08]` |
| crimm 2024.7a2 is installed editable from `/home/ziqiaoxu/crimm`, with its `.pth` and finder in the **user site** (`~/.local/lib/python3.12/site-packages`), not in `crimm-dev` itself. Any Python 3.12 of the owner that sees the user site gets the local checkout, not the PyPI release. The local `pyproject.toml` says 2026.2.1. | `[ran 2026-10-08]` pip show |
| **2026-10-30:** the worktree `/home/ziqiaoxu/crimm-hotfix` is on branch `fix/59-rdkit-plain-residue` at `1816fbb`, which is what this project imports; this project's suite and crimm's `tests/test_rdkit_converter.py` pass on it. `msldprep` does not call the converter, so only `align.ipynb` needs the fix. Earlier text, branch and commit now out of date: since 2026-10-12 the owner's sessions in this project put crimm's git worktree `/home/ziqiaoxu/crimm-hotfix` (branch `RDKConverter-Het`, `8075b75`) first on `PYTHONPATH`, through the `env` block of `.claude/settings.local.json` (local, git-excluded). It shadows the editable install. Check with `python -B -c "import crimm; print(crimm.__file__)"`. The worktree exists because the editable `/home/ziqiaoxu/crimm` has that branch switched off. | `[ran 2026-10-12]`, `[told 2026-10-12]` |
| Since 2026-10-12 `crimm-dev` also has the packages for the GUI: fastapi 0.143.0, uvicorn 0.54.0, starlette 1.7.0, pydantic 2.14.1, python-multipart 0.0.32, httpx 0.28.1 (httpx was already there; `TestClient` needs it). Installed with `pip install` at the owner's request, nothing else changed. `pip check` reports `nglview 4.0 has requirement numpy<2.3, but you have numpy 2.5.3`; that predates this install (numpy was not touched). Offline suite after the install: 78 passed, 58 skipped. | `[ran 2026-10-12]` |
| The suite also passes, 136 passed, against the latest PyPI crimm (2026.2.2, wheel installed with `--no-deps` into a scratch directory and put first on `PYTHONPATH`). | `[ran 2026-10-08]` |
| Tools are found through environment variables, else `PATH`: `MSLDPREP_CGENFF_PATH` (program), `MSLDPREP_CGENFF_BASE_PRM` (base prm; default `par_all36_cgenff.prm` beside the program; the root script and tests use that default, the CLI omits `-p` if unset), `MSLDPREP_OBABEL_PATH` (obabel). With none found, 5 CGenFF tests and 53 OpenBabel tests skip. Installs that work here: `/export/apps/RockyOS8/cgenff/src/2026_May/silcsbio.2026.1/cgenff/cgenff` and `/export/apps/RockyOS8/openbabel/3.1.1/bin/obabel`; crimm's handoff names a different CGenFF (`.../src/2026/silcsbio.2026.1-alpha/...`). `obabel` is not on `PATH` in a non-activated shell. | `[ran 2026-10-08]` |
| `tests/test_writers.py::test_downstream_msld_prm_contract` imports `msld_prm` from `/home/ziqiaoxu/msld-py-prep` and runs it. | `[read]`, `[ran]` passes |
| SUPERSEDED 2026-10-16 by the owner's rule to push over SSH: `ssh -T git@github.com` authenticates as `Truman-Xu` with `~/.ssh/id_rsa`, and `git push git@github.com:Truman-Xu/msld-prep.git main` worked (from the agent's shell only with the sandbox disabled, since it needs the network). Old text: a plain `git push` fails here: the remote is https, there is no credential helper and no display for askpass. Worked, without changing any git config: `git -c credential.helper= -c "credential.helper=!/home/ziqiaoxu/.conda/envs/crimm-dev/bin/gh auth git-credential" push origin main`, with `GIT_TERMINAL_PROMPT=0`. Check `git push --dry-run` first (G14). | `[ran 2026-10-12]`, `[ran 2026-10-16]` |
| Branch `main`, tracking `origin/main`, **public**, issues enabled, pushed 2026-10-08. `origin` is `git@github.com:Truman-Xu/msld-prep.git` since 2026-10-21 (was https); on 2026-10-21 the push worked from the agent's shell without disabling the sandbox. `gh` is at `/home/ziqiaoxu/.conda/envs/crimm-dev/bin/gh`, logged in as `Truman-Xu` (scopes `repo`, `read:org`, `admin:public_key`, `gist`). `gh` is not on `PATH`. | `[ran 2026-10-08]` |

Commands that work `[ran 2026-10-08]`:

```bash
PY=/home/ziqiaoxu/.conda/envs/crimm-dev/bin/python
export MSLDPREP_CGENFF_PATH=/export/apps/RockyOS8/cgenff/src/2026_May/silcsbio.2026.1/cgenff/cgenff
export MSLDPREP_OBABEL_PATH=/export/apps/RockyOS8/openbabel/3.1.1/bin/obabel
$PY -B -m pytest -p no:cacheprovider -q        # 136 passed, 0 skipped, about 11 s
env -u MSLDPREP_CGENFF_PATH -u MSLDPREP_OBABEL_PATH $PY -B -m pytest -p no:cacheprovider -q   # 78 passed, 58 skipped
```

## 4. Gotchas

Things that cost a previous agent time or would produce wrong results silently.

| ID | Trap | Tag |
| --- | --- | --- |
| G1 | Running anything that imports crimm without `python -B` rewrites `.pyc` files inside `/home/ziqiaoxu/crimm`. They are gitignored and harmless, but they are writes into a project you must not touch. Use `-B` or `PYTHONDONTWRITEBYTECODE=1`. | `[ran 2026-10-08: 41 files refreshed by a run without -B]` |
| G2 | Because crimm is an editable install from `/home/ziqiaoxu/crimm`, whichever crimm branch is checked out there decides what msldprep sees. It was `cleanup/phase0-tests-and-dev-docs` on 2026-10-08. A test failure here may come from that, not from this repo. | `[ran 2026-10-08]` |
| G3 | `ConstrainedEmbed` is not fully deterministic, which is why `build` reloads the aligned conformers from `{aligned_dir}/aligned_ligands.sdf` instead of re-aligning. Do not "simplify" step 2 to re-run alignment. | `[read pipeline.py docstring]` |
| G3a | SUPERSEDED in part 2026-10-12 by claude-sonnet/claude-code/m0-align-api, for G3's claim only: RDKit's `ConstrainedEmbed` defaults to `randomseed=2342`, and the same analog aligned twice onto the same reference in one process gave identical coordinates (3 ChEMBL analogs, `np.allclose` and exact equality in the new tests). So the observed unreliability is not run-to-run noise for identical input. What was *not* checked: whether a different atom order, RDKit version, platform or reference conformer changes the embed, which may be what G3 saw. Keep reloading the sdf in `build`; do not rely on this to re-align. `align_to_reference(..., seed=)` now exposes the seed. | `[ran 2026-10-12]` |
| G16 | The GUI runs one in-memory `AlignmentSession` per server, so two people on one server share and overwrite each other's work. The access token is the only protection, and the login node is shared: do not leave a server running with a guessable token (I tested with `demo123` and stopped it). The default port is a free one picked at start, because 8765 was already taken by another user's server on 2026-10-14. | `[ran 2026-10-14]` |
| G17 | An agent cannot see the 3D view here: headless Firefox 140 on this node has no WebGL even with software-GL preferences (`webgl: False`), so `draw()` in `static/app.js` was **not** run in a browser. What was run in a real browser (Firefox through Marionette on port 2828, a throwaway profile): the whole page flow (load, align, select, accept, align again, download link) with no JavaScript errors, and 3Dmol's per-atom selection by `serial` (29 atoms, 17 core, 0 mismatches). To drive a headless browser from here, speak the Marionette protocol directly; there is no selenium or playwright. `pkill -f <pattern>` matches the agent's own shell command and kills it; use `pgrep`, then `kill <pid>`. | `[ran 2026-10-14]` |
| G17a | SUPERSEDES G17's "no WebGL" for this node, `[ran 2026-10-14]`: WebGL works in Firefox 140 with Mesa's software renderer (`llvmpipe`) on a virtual display. Nothing was installed: `/usr/bin/Xvnc` and Mesa's `swrast_dri.so` are system packages. Recipe: `Xvnc :77 -geometry 1500x1000 -depth 24 -rfbport -1 -nolisten tcp -SecurityTypes None -rfbunixpath <scratch>/vnc.sock &`, then `DISPLAY=:77 LIBGL_ALWAYS_SOFTWARE=1 firefox --marionette --no-remote --profile <throwaway profile> about:blank &`, with these in the profile's `user.js`: `webgl.force-enabled` true, `webgl.disabled` false, `layers.acceleration.disabled` true, `webgl.forbid-software` false, `webgl.disable-angle` false. Headless mode (`MOZ_HEADLESS`) did not get WebGL; software-WebRender prefs crashed Firefox. Drive it with `scratchpad`-style Marionette code (G17); `TakeScreenshot` then shows the 3D view. **Stop Firefox (Marionette listens on 127.0.0.1:2828) and Xvnc when done**: other users of the login node can reach loopback ports and the X socket. Kill by pid (`pgrep -x firefox`, `pgrep -x Xvnc`), not `pkill -f`. | `[ran 2026-10-14]` |
| G18 | To check how Safari's engine behaves, use WebKitGTK 2.54, which is installed with Python bindings for the system Python: `/usr/bin/python3` with `gi.require_version('WebKit2', '4.0')`, a `WebKit2.WebView` in a `Gtk.Window` on the Xvnc display (G17a), `run_javascript` to read the page (it cannot return a Promise: store the result in `window` and read it in a second call), and `get_snapshot` for a picture (WebGL works). `document.elementFromPoint` can miss an overlap that Safari has: WebKitGTK drew the hidden file inputs 1 px wide, so measure the boxes of the elements as well. Stop Xvnc afterwards. | `[ran 2026-10-16]` |
| G15 | **No longer true for alignment since 2026-10-14** (`align_to_reference` matches by element and ignores isotope labels); still true for anything that matches a query built from labeled atoms. Original text: `align.prepare_ligand_set` / `assign_isotope_labels` set isotope labels on the **reference** `Mol` too. A labeled reference no longer matches an unlabeled analog, so `align_common_core` / `align_to_reference` then raise `ValueError: molecule doesn't match the core`. Align before labeling, or use a freshly loaded reference (the new tests do). The session-scoped `reference_ligand` fixture in `tests/conftest.py` is labeled once `prepared_ligand_set` has run. The GUI must use a copy and never label. | `[ran 2026-10-12]` |
| G19 | `Chem.RemoveHs` keeps hydrogens that carry an isotope label (RDKit's default `RemoveHsParameters`). Once `assign_isotope_labels` has run, every "heavy-atom copy" made with `RemoveHs` still has all its hydrogens: `align.find_shared_core` searched with them (#15, **fixed 2026-10-29**: hydrogens are labeled 0 and the copies are made with `Chem.RemoveAllHs`). Also, `FindMCS` query atoms built with `CompareIsotopes` report atomic number 0, so count atoms instead of testing elements. | `[ran 2026-10-21]` |
| G20 | **Fixed 2026-11-02 (#18), the file still carries none:** `import_mcs_text` folds lone pairs in from each ligand's rtf by host atom, so `build` has them again; anything that reads `MCS_for_MSLD.txt` without the ligands still sees none. Original text: `MCS_for_MSLD.txt` does not carry lone pairs: `export_mcs_text` writes none and `import_mcs_text` reads them only from SITE lines, so anything built after re-reading the checkpoint (`build_from_core_definition`, the `build` CLI command) has no lone pairs (#18, strict xfail). A decomposition used directly in memory still has them. Also: the checkpoint's CORE lines may carry `# TIE` and `# FAR` comments (since `9d00ab2`); the importer ignores them, legacy `msld_crn.py` would not. | `[ran 2026-10-29]` |
| G21 | **Fixed 2026-11-02 (#19), when `obabel` is found:** `mol2_reader` reads them through OpenBabel; without it they are refused with the cause, so a test environment without `obabel` still cannot load these 20. Original text: 20 of the 89 mol2 files in `examples/2_ChEMBL` do not load through `LigandRecord.from_files` (#19, strict xfail), 8 of them in its `mol_list.txt`, so any ChEMBL figure over the full list is over fewer ligands unless they are recovered. To recover them outside the library: `obabel -imol2 X.mol2 -osdf`, read with `sanitize=False`, set +1 on a 4-valent neutral N bonded to `O-`, then `SanitizeMol`; atom order is kept, so names can be set from the mol2. | `[ran 2026-10-30]` |
| G22 | `/export/apps/RockyOS8/openbabel/3.1.1/bin/obabel` run directly fails on `libcoordgen.so.3`. `tests/test_sybyl_typer.py::_obabel_env` builds `LD_LIBRARY_PATH` for it, so the suite passes. For a direct call use `/home/ziqiaoxu/.conda/envs/crimm-dev/bin/obabel` (3.2.1). | `[ran 2026-10-30]` |
| G4 | Aligned ligands are stored as one sdf, not per-ligand mol2, because RDKit's mol2 reader sometimes fails to re-kekulize an aromatic ring. Atom names ride in the `AtomNames` data field. | `[read align.py]` |
| G5 | CGenFF atom types are carried on RDKit atoms as the `CGENFFAtomType` property. **Since 2026-10-29** the isotope labels (from 300) stand for an atom's *pattern*, its type plus its neighbors' types (`align.atom_pattern`, `assign_pattern_labels`), on heavy atoms only; hydrogens get 0. They are written into `aligned_ligands.sdf` with the structures. Do not use isotopes for anything else in this package. | `[read]`, `[ran 2026-10-29]` |
| G6 | `msldprep/sybyl_typing_data.py` is generated from `data/atomtyp.txt` and `data/types.txt`. Edit the generator or inputs, then rerun it. | `[read]` |
| G7 | Docstrings and configs still say `rewrite/msld_prep_cli.py` and `rewrite/scripts/...`, from when this code lived in a `rewrite/` directory. Four matches in `.py` and `.yaml` files. The real paths have no `rewrite/`. | `[ran grep 2026-10-08]` |
| G8 | Untracked does not mean disposable. `build.*/` and `work.*/` are git-ignored but still hold the owner's runs; do not delete them. `TODO.md` was untracked; its two items are issues #1 and #2, and the owner deleted the file on 2026-10-14. | `[inferred]` |
| G9 | Other sessions may edit `docs/dev/*.md` concurrently. Re-read a file immediately before editing it, and make targeted edits, never whole-file rewrites. | `[inferred: same hazard recorded in crimm's handoff]` |
| G10 | **2026-11-02 (#17):** the `test` extra pins `pandas<3`, and both legacy contract tests skip under pandas 3, so in `crimm-dev` (pandas 3.0.6) they skip rather than pass vacuously; under pandas 2.3.3 they pass. Original text: `test_downstream_msld_prm_contract` puts `/home/ziqiaoxu/msld-py-prep` on `sys.path` and emits pandas `ChainedAssignmentError` warnings from that legacy file. The warnings are in legacy code, not in this repo. **They are not harmless (2026-10-28):** under pandas 3 the chained assignments write nothing, so legacy `msld_prm` finds no bonds, angles or dihedrals and the test passes without testing anything (#17, strict xfail `test_downstream_msld_prm_finds_bonded_terms`). Legacy `msld_crn.py` does not run at all under pandas 3 (`KeyError: 0`, line 326). To run legacy code, use `/home/ziqiaoxu/.conda/envs/autogrow/bin/python` (Python 3.10, pandas 1.5.0, numpy 1.23; no RDKit check made) for `msld_crn`/`msld_prm`; `msld_mcs` runs in `crimm-dev`. | `[ran 2026-10-08]`, `[ran 2026-10-28]` |
| G11 | `/home/ziqiaoxu/msld-prep-tools` is a symlink to `/home/ziqiaoxu/msld-prep` (published to cluster users; see standing instructions). Use the new path in docs and new scripts; do not break the link. | `[told 2026-10-08]`, `[ran]` |
| G12 | The owner's untracked files `align.ipynb`, `analogues.txt`, `sil.mol2`, `sil.prm`, `sil.rtf`, `cellz.yaml` and `cellz_build.yaml` belong to the unpublished analog series (the notebook holds its structures, `analogues.txt` lists them, `cellz*.yaml` point to them). They are listed in `.git/info/exclude`, so `git add .` skips them. That file is local; a fresh clone has no such protection. `cellz.yaml` also holds a CGenFF install path. | `[read]` heads of the files, `[ran]` exclude |
| G13 | **Fix available 2026-10-30:** crimm `1816fbb` replaces `8075b75` (Q12); the notebook has not been re-run with it. Original text: `align.ipynb` reads its inputs from absolute paths in another lab member's directory, and was run only as far as cell 9. There it fails: `crimm/Adaptors/RDKitConverter.py::heterogen_to_rdkit` reads `het_res.pdbx_description`, which a `Residue` read from a CRD file lacks. The owner's hotfix is commit `8075b75` on crimm's local branch `RDKConverter-Het`, not the checked-out branch (`cleanup/phase0-tests-and-dev-docs`). Cells after 9 were not run. | `[ran 2026-10-08]` failure; `[inferred]` cause |
| G14 | Rewriting history (`git filter-branch`, then `gc --prune=now`) was done on 2026-10-08 while nothing was pushed. **The repository is public and was pushed the same day, so from now on history rewrites are not a clean undo**: a published file stays in forks and caches. Check `git ls-files`, `git log --stat` and `git push --dry-run` before every push. | `[ran]` |

## 5. State, claims, open questions and log

These change every session, so they are not kept in this file. Since 2026-11-07 `[told]`:

- **State, claims and the log** are in the open issue labelled `agent-log`, pinned on the
  issue list (first one: #27). Its description holds Claims and State, rewritten in place;
  each comment is one session's log entry. Find it with
  `gh issue list --label agent-log --state open`.
- **Open questions** are issues labelled `question`. Former questions still cited elsewhere:
  Q12 is #25, Q13 is #26 (cited as "handoff Q13" in `DESIGN_NOTES.md` and
  `ISSUES_AND_PLANS.md`). Q3 was dropped: `TODO.md` no longer exists. The others were
  answered.
- **The earlier log**, 2026-10-08 to 2026-11-07, with State, Claims and the open questions as
  they stood, is sections 5 to 8 of this file at commit 2df1ac9:
  `git show 2df1ac9:docs/dev/AGENT_HANDOFF.md`.
