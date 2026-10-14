# Issues and plans

Last updated 2026-10-14. **This is not a roadmap.** A roadmap needs an order and reasons for the
order, and the order depends on how the original `msld-py-prep` pipeline works and what the new
pipeline's inputs and outputs should be. The owner is still learning that (see
[AGENT_HANDOFF.md](AGENT_HANDOFF.md), Q13), so this file does only two jobs: it lists the
**active issues** and the **planned features**, and it gives the templates for adding more. The
roadmap is written later, when the owner has decided the design; a `ROADMAP.md` is not created
until then.

Where things live:

| Kind of information | Where |
| --- | --- |
| Status and discussion of one issue | the GitHub issue, `https://github.com/Truman-Xu/msld-prep/issues/<N>` |
| What is open, and the evidence for it | this file |
| Who did what, traps, decisions | [AGENT_HANDOFF.md](AGENT_HANDOFF.md) |
| Standards the work should meet | [BEST_PRACTICES.md](BEST_PRACTICES.md) |

## How to use this file

- **A bug or a missing piece of behaviour** goes in "Active issues" the day it is found, even
  before it is understood. Give it a GitHub issue (the issue form "Task"), a strict `xfail`
  test that states the behaviour that is wanted, and a row below. An item without an issue
  number has a provisional ID, `P-<n>`, until the issue exists; then the row takes the issue
  number and the `xfail` reason is updated.
- **An idea for a new capability** goes in "Planned features". It stays there, unordered, until
  the roadmap exists.
- **Evidence beats opinion.** Say what was run, on which input, and what came out. Say what was
  not checked. Use only the public examples (`examples/2_ChEMBL`) in a report: the lab's
  unpublished analog series must never appear in an issue, a test or this file
  (handoff section 2, G12).
- **Closing** an item: move it to "Done" with the date and the commit, and remove the `xfail`
  marker in the same change.

Template for an issue:

```markdown
### <ID> Short title (#N)

- **Status:** open | investigating | fix planned | fixed (date, commit)
- **Found:** date, by whom, how
- **What happens:** the symptom, with the exact message
- **Evidence:** inputs, commands, output; what was run and what was only read
- **Cause:** known, suspected (say why), or unknown
- **Wanted:** the behaviour that should hold, as a checkable sentence
- **Workaround:** if there is one
- **Test:** the `xfail` test, or why none is possible
- **Not checked:** what is still open
```

Template for a feature:

```markdown
| ID | Feature | Why (who needs it) | Depends on | Issue | Notes |
```

## Active issues

| ID | Title | Kind | Status | Test |
| --- | --- | --- | --- | --- |
| #9 | Alignment fails with "molecule doesn't match the core" when the **reference** has a charged atom | bug | open, cause found, workaround in the GUI | `tests/test_align.py::test_align_to_reference_matches_atoms_whatever_their_charge` (strict `xfail`) |
| #1 | Reference and analog protonation states can mismatch | correctness | open, not reproduced | none yet |
| #2 | Respect chirality of molecules through the pipeline | correctness | open, not reproduced | none yet |

### #9 Matching the common structure regardless of charge

- **Status:** open. GitHub issue #9, created 2026-10-14.
- **Found:** 2026-10-14 by the owner, using the alignment GUI: a reference with a charged
  nitrogen gave "molecule doesn't match the core". Reproduced and diagnosed the same day with
  public data (below).
- **What happens:** `msldprep.align.align_to_reference` (and so `align_common_core`,
  `analog_prep.generate_3d_analog` and the GUI) raises `ValueError: molecule doesn't match the
  core` for every analog when the reference carries a formal charge that the analog does not
  have. The message points at nothing: the molecules do share a large structure.
- **Evidence** (`[ran 2026-10-14]`, RDKit 2026.09.1, reference and analog taken from the
  public ChEMBL example series):
  - reference `c1ccc2c(c1)c([nH+]cn2)Nc1cccc(Br)c1` (ring nitrogen protonated), analog
    `Oc1ccc2c(Nc3cccc(Br)c3)ncnc2c1` (neutral): `ValueError: molecule doesn't match the core`;
  - reference `c1ccc2c(c1)c(ncn2)[NH2+]c1cccc(Br)c1` (exocyclic nitrogen protonated), same
    analog: same error;
  - the neutral reference `c1ccc2c(c1)c(ncn2)Nc1cccc(Br)c1` with the same analog works (core
    of 18 atoms);
  - the **reverse** (neutral reference, charged analog `Oc1ccc2c(Nc3cccc(Br)c3)[nH+]cnc2c1`)
    also works, so only a charged reference breaks the match.
- **Cause** (`[ran]`, diagnosed step by step on the first case): `FindMCS` compares elements
  only, so it finds the full 18-atom common structure and the search succeeds. The next step,
  `Chem.PathToSubmol(ref_mol, bond_path, useQuery=True)`, turns each reference atom into a
  query that requires its atomic number **and, for a charged atom, its formal charge**
  (`AtomAnd(AtomAtomicNum 7, AtomFormalCharge 1)`). `ConstrainedEmbed` then calls
  `mol.GetSubstructMatch(core)` on the neutral analog, finds nothing, and raises. The code is in
  `msldprep/align.py::align_to_reference`; the raise is inside RDKit's `ConstrainedEmbed`.
  Whether other reference properties (explicit hydrogen counts, isotope labels, see G15) can
  break the match the same way was **not checked**.
- **Wanted:** the common structure is found and used by element and connectivity, whatever the
  charges of the atoms in it. If a reference and an analog differ in net charge, say so in a
  clear message (this overlaps #1) instead of reporting a failed match.
- **Workaround:** the GUI option "Remove charges from the reference" (RDKit's `Uncharger`, as
  in the notebook) neutralizes the reference first; with it, the first two cases above align
  (`[ran]`). The GUI also adds a hint to the error when the reference is charged. Removing the
  charge changes the molecule, so use it only when the analogs really are neutral.
- **Directions for a fix, not decided:** (a) build the core query from element-only atoms
  instead of `PathToSubmol(useQuery=True)`; (b) use the atom pairs `FindMCS` already gives us
  and embed with them directly (`EmbedMolecule(coordMap=...)`) instead of letting
  `ConstrainedEmbed` match again; (c) neutralize both sides before matching. Each changes what
  is aligned in the charged case, so the owner decides after seeing what a charged reference
  should mean for MSLD (#1).
- **Test:** `tests/test_align.py::test_align_to_reference_matches_atoms_whatever_their_charge`
  (two cases, `xfail(strict=True)`, reason names this entry), plus
  `test_a_charged_analog_aligns_onto_a_neutral_reference`, which passes today and must keep
  passing.
- **Not checked:** the effect on the rest of the pipeline (typing, the shared core search,
  charge renormalization) when charged atoms are present at all.

## Planned features

Unordered. Nothing here is scheduled.

| ID | Feature | Why | Depends on | Issue | Notes |
| --- | --- | --- | --- | --- | --- |
| F-1 | Alignment GUI, first slice | demo to the lab; check alignments by eye | none | #3 to #8 | Built 2026-10-14, see "Done". Open: 3D view not yet seen in a browser by an agent; issues #4 to #8 describe a wider plan than what was built |
| F-2 | Show the protein pocket in the GUI for clash checks | the workshop notebook asks for it | F-1 | none | load a PDB for display only |
| F-3 | Export Sybyl-typed mol2 per analog from the GUI | what the notebook hands to CGenFF | F-1 | #7 | dropped from the first slice on 2026-10-14 (no atom typing wanted yet) |
| F-4 | Several users or sessions on one GUI server | today one session per server (G16) | F-1 | none | needs session ids in the API |
| F-5 | Take already-parameterized sets (mol2 + rtf + prm) in the GUI | reuse existing CGenFF output | F-1 | none | |
| F-6 | Carry the GUI's aligned set into `define_common_core` | one tool from alignment to the build | F-1, design decisions | none | blocked on the learning plan |
| F-7 | Fewer saved files: keep intermediate stages in memory | owner's goal (handoff section 2, 2026-10-12) | learning plan (Q13): which files CHARMM and ALF really read | none | |
| F-8 | A web GUI for the whole preparation, not only alignment | owner's goal | F-6, F-7 | none | |
| F-9 | Compare the old and new pipelines on real MSLD runs | decide the design with evidence | learning plan, CHARMM basics from colleagues | none | |
| F-10 | Learning plan for the original `msld-py-prep` | prerequisite of the roadmap | none | none | postponed until after the GUI demo; seed in handoff Q13 |

## Open design questions

From handoff Q13; **none verified**. Answer them while learning the original pipeline; each
answer may create an issue or a feature above.

1. The original merges identical fragments per site (for example `NSUBS 19 24` for 31 ligands);
   msld-prep appears not to (`nsubs` 34 34 34 for 34 ligands). Intended or a regression?
2. Does msld-prep's proportional charge redistribution meet the same constraints as the
   original's cyclic 1e-6 loop, and does the difference matter?
3. The original's `.str` splitting has no counterpart on the `from-files` path; a config
   comment about `.str` may be stale.
4. Which output files do CHARMM or ALF really read, and which are only review artefacts?
5. `bookending_corrections/` and `Lg_Solvate.sh` of the original have no counterpart here.

## Roadmap

Not written. It waits for the learning plan (F-10) and for the owner's decisions on the
pipeline's inputs and outputs. When it exists it will take the items above, put them in an
order with reasons, and link each to its issue.

## Done

| Date | What | Where |
| --- | --- | --- |
| 2026-10-12 | `align_to_reference` (atom map, RMSD, seed, MCS options) and the text loaders | #3, `e65a2e4`, `3d9c979` |
| 2026-10-14 | Alignment GUI first slice: load a reference and analog SMILES, constrained embed, 3D and 2D views, accept or reject, align again with a seed, sdf download | F-1; `msldprep/gui/`, uncommitted when written |
| 2026-10-14 | GUI option to remove charges from the reference with `Uncharger` (the notebook's step) | `msldprep/gui/service.py::neutralize` |
| 2026-10-14 | GUI toggle to color the reference by element; gray when off | `msldprep/gui/static/app.js` |
