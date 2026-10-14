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
| #9 | Alignment fails with "molecule doesn't match the core" when the **reference** has a charged atom | bug | **fixed 2026-10-14**, closes with the commit that says `Closes #9` | `tests/test_align.py::test_align_to_reference_matches_atoms_whatever_their_charge` (the `xfail` marker is removed) |
| P-1 | Which symmetric match is embedded is arbitrary, so a ring can come out flipped | bug (suspected) | open, seen once, not judged | none yet |
| #1 | Reference and analog protonation states can mismatch | correctness | open, not reproduced | none yet |
| #2 | Respect chirality of molecules through the pipeline | correctness | open, not reproduced | none yet |

### #9 Matching the common structure regardless of charge

- **Status:** fixed 2026-10-14 (see **Fix**); GitHub issue #9, created the same day.
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
- **Related, found 2026-10-14 on the owner's real data:** a reference file whose rings carry
  formal charges (net +3) is perceived by RDKit as non-aromatic in the fused ring system, while the
  analogs' rings are aromatic. Strict bond comparison then mismatched 10 ring bonds and
  `completeRingsOnly` dropped both rings (core 20 of 32 atoms). `McsOptions.ignore_bond_orders`
  fixes it (core 31 to 33 atoms). Removing charges does not (`Uncharger` leaves the N+ of an
  iminium-like ring nitrogen). Cleaning such a reference (re-assigning aromatic bond orders) is not
  planned; it matters for #1 (protonation states).
- **Workaround:** the GUI option "Remove charges from the reference" (RDKit's `Uncharger`, as
  in the notebook) neutralizes the reference first; with it, the first two cases above align
  (`[ran]`). The GUI also adds a hint to the error when the reference is charged. Removing the
  charge changes the molecule, so use it only when the analogs really are neutral.
- **Fix** (2026-10-14, `[ran]`): `align_to_reference` no longer rebuilds a core query from the
  reference and matches it again. It matches the query that `FindMCS` returned (element-only
  atoms) back onto both molecules and embeds with those atom pairs, in the helper
  `_embed_onto_reference`, adapted from RDKit's `ConstrainedEmbed`. Charges, hydrogen counts and
  isotope labels (so also handoff G15) no longer matter. On 33 demo analogs the new code gave the
  same atom map and coordinates within 0.03 angstrom for 32, and a different but equally good
  symmetric map for one (see P-1).
- **Test:** `tests/test_align.py::test_align_to_reference_matches_atoms_whatever_their_charge`
  (two cases; it was a strict `xfail` and became a plain test with the fix), plus
  `test_a_charged_analog_aligns_onto_a_neutral_reference` and
  `test_the_core_is_matched_by_element_so_an_isotope_labeled_reference_still_aligns`.
- **Not checked:** the effect on the rest of the pipeline (typing, the shared core search,
  charge renormalization) when charged atoms are present at all.

### P-1 The symmetric match that is embedded is arbitrary

- **Status:** open, provisional ID (no GitHub issue yet).
- **Found:** 2026-10-14 while checking the fix for #9.
- **What happens:** when the common structure matches a molecule in more than one way (a
  symmetric ring, or two equivalent substituents), `GetSubstructMatch` returns the first match
  and that one is embedded. Which one comes first is an accident of the query. Two ways of
  building the query for the same MCS gave different atom maps for one of 34 demo analogs,
  and the two poses differ by up to 4 angstrom in places, with the same core RMSD (0.033).
- **Cause:** suspected, `[read]`: all matches have equal core RMSD for a symmetric core, so
  nothing in the code prefers one.
- **Wanted:** among the matches, embed the one that puts the atoms outside the core where the
  reference has atoms (or the one with the least clash), or let the user flip it in the GUI.
- **Not checked:** which of the two poses is closer to the reference's non-core atoms; whether
  the owner's flipped-ring cases (workshop notebook, "check the ring is not flipped") are this
  problem; how often it happens on the real series.

## Planned features

Unordered. Nothing here is scheduled.

| ID | Feature | Why | Depends on | Issue | Notes |
| --- | --- | --- | --- | --- | --- |
| F-1 | Alignment GUI, first slice | demo to the lab; check alignments by eye | none | #3 to #8 | Built 2026-10-14, see "Done"; the 3D view was checked in Firefox with software WebGL on 2026-10-14. Issues #4 and #8 were rewritten on 2026-10-14 to match what was built; #5 to #7 still describe the first plan and are partly superseded |
| F-2 | Show the protein pocket in the GUI for clash checks | the workshop notebook asks for it | F-1 | none | load a PDB for display only |
| F-3 | Export Sybyl-typed mol2 per analog from the GUI | what the notebook hands to CGenFF | F-1 | #7 | **Done 2026-10-14**, see "Done" |
| F-4 | Several users or sessions on one GUI server | today one session per server (G16) | F-1 | none | needs session ids in the API |
| F-5 | Take already-parameterized sets (mol2 + rtf + prm) in the GUI | reuse existing CGenFF output | F-1 | none | |
| F-6 | Carry the GUI's aligned set into `define_common_core` | one tool from alignment to the build | F-1, design decisions | none | blocked on the learning plan |
| F-7 | Fewer saved files: keep intermediate stages in memory | owner's goal (handoff section 2, 2026-10-12) | learning plan (Q13): which files CHARMM and ALF really read | none | |
| F-8 | A web GUI for the whole preparation, not only alignment | owner's goal | F-6, F-7 | none | |
| F-9 | Compare the old and new pipelines on real MSLD runs | decide the design with evidence | learning plan, CHARMM basics from colleagues | none | |
| F-11 | Option to keep the original atom names in the mol2 export | the owner's choice to rename is deliberate (CGenFF), but some uses need the original names | F-3 | #10 | default stays "new unique names"; refuse kept names that are duplicated or longer than four characters |
| F-12 | Manual tuning of an aligned analog in the GUI: rotate a bond to fix a small misalignment | the constrained embed can leave a free part in a slightly wrong pose | F-1 | none yet | 3Dmol.js has no bond-rotation tool (checked 2026-10-14: no torsion or dihedral code in the bundled 2.4.2), only atom picking (`setClickable`, `setHoverable`) and `setCoordinates`. Likely design: pick a bond in the page, send it with an angle to the server, rotate the atoms on the side that holds no matched-core atom with RDKit (`rdMolTransforms`), return the new coordinates, show clashes. The core must stay fixed. Not decided |
| F-10 | Learning plan for the original `msld-py-prep` | prerequisite of the roadmap | none | none | postponed until after the GUI demo; seed in handoff Q13 |

## Deliberate behaviors

Things that look like gaps but are choices. Do not "fix" them; change them only by an item above.

- **mol2 export renames every atom and writes no charges** (owner, 2026-10-14). The GUI's mol2
  files are written for CGenFF. CGenFF keeps only unique names of at most four characters
  (it silently renumbers others) and computes its own charges, so the writer
  (`sybyl_typer.assign_unique_atom_names`, `format_mol2_text`) gives every atom a fresh unique
  name, the reference included, and writes `NO_CHARGES`. The reference's names from its own file
  are therefore not kept. An option to keep them is F-11 (#10). The sdf export keeps the
  structure as it is and adds each analog's RMSD and seed.
- **A realigned analog loses its accept or reject mark**, because the pose is new.
- **The reference is used with hydrogens**: added with coordinates if the file has none, and
  re-added after removing its charges.

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
| 2026-10-14 | Alignment matches by element and ignores charges and hydrogen counts; option to ignore bond orders so tautomers and protonation variants match in full (on by default in the GUI) | #9; `msldprep/align.py::align_to_reference`, `McsOptions.ignore_bond_orders` |
| 2026-10-14 | Collapsible 2D drawing of the reference in panel 1; with an analog selected, the atoms it matched are marked | `msldprep/gui/static/app.js::renderReference2D`, `ligand_view()["reference_svg"]` |
| 2026-10-14 | GUI download as sdf or mol2 (the notebook's Sybyl typing and mol2 writer), one file or one file per structure (zip), with or without the reference | F-3, #7; `msldprep/gui/service.py::AlignmentSession.export` |
