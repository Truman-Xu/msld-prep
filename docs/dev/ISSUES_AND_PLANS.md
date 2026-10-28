# Issues and plans

Last updated 2026-10-21. **This is not a roadmap.** A roadmap needs an order and reasons for the
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
| Design decisions for the re-design | [DESIGN_NOTES.md](DESIGN_NOTES.md) |
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
| #9 | Alignment fails with "molecule doesn't match the core" when the **reference** has a charged atom | bug | **fixed 2026-10-14**, closed by `e670117` | `tests/test_align.py::test_align_to_reference_matches_atoms_whatever_their_charge` |
| #11 | Which symmetric match is embedded is arbitrary, so a ring can come out flipped | bug | **fixed 2026-10-16**, closed by `56f697c` | `tests/test_align.py::test_a_decorated_ring_turns_its_substituent_toward_the_reference_substituent` |
| #13 | The core search kept RDKit's 3600 s limit, so one alignment could take an hour | bug | **fixed 2026-10-16**, closed by `9002a9e` (10 s by default) | `tests/test_align.py::test_the_core_search_has_a_ten_second_limit_by_default` |
| #15 | The shared-core search keeps every hydrogen: `Chem.RemoveHs` keeps isotope-labeled atoms, so `align.find_shared_core` runs `FindMCS` with all hydrogens, against its own comment | bug | open, found 2026-10-21; may explain part of the slowness of #14. **Not a one-line fix:** the hydrogens in the query are what puts core hydrogens in the core today (and hydrogens of differing type make hydrogen-only sites); the fix needs an explicit hydrogen rule, see the comment on #15 and `DESIGN_NOTES.md` (i) | `tests/test_align.py::test_the_shared_core_search_leaves_hydrogens_out` (strict `xfail`) |
| #1 | Option to require matching protonation states between reference and analog | enhancement | **done 2026-10-28**: `McsOptions.check_protonation` (off by default) reports core atoms whose formal charge or total degree differs, in `CoreAlignment.protonation_differences`; the match and pose are unchanged (owner's choice 2026-10-28: flag, not exclude or refuse). GUI setting and per-analog lines. Tolerance stays the default (#9, `ignore_bond_orders`) | `tests/test_align.py::test_the_protonation_check_flags_core_atoms_and_leaves_the_pose_as_it_is`, `test_the_protonation_check_does_not_flag_a_substituent_in_place_of_a_hydrogen`, `tests/test_gui_service.py::test_protonation_changes_are_flagged_only_when_asked_and_the_pose_stays` |
| #2 | Respect chirality of molecules in the structure alignment | correctness | **closed 2026-10-21**, scope limited to alignment by the owner, who is satisfied with it; the whole pipeline is tested after the re-design. Before: alignment settles and reports stereo (`b51a361`) and puts the core exactly on the reference, falling back when that would invert a centre (`f43e56d`); decomposition and the CHARMM writers were not examined, and the GUI is limited to alignment | `tests/test_align.py` (stereochemistry section, `test_the_core_keeps_the_reference_geometry_ring_puckers_included`) |
| #16 | GUI: a mol2 reference without hydrogens is used without them: RDKit's mol2 reader marks atoms `NoImplicit`, so `AddHs` in `set_reference` adds none | bug | open, found 2026-10-28 while working on #10; fix or refuse is the owner's call | `tests/test_gui_service.py::test_a_mol2_reference_without_hydrogens_gets_them` (strict `xfail`) |
| #17 | The legacy `msld_prm` contract test passes without testing anything under pandas 3: chained assignment writes nothing, so legacy finds no bonded terms | bug (test) | open, found 2026-10-28 in lesson 3; legacy is read-only here, so the fix is a pandas < 3 interpreter for the test, an upstream fix, or retiring the test: the owner's call | `tests/test_writers.py::test_downstream_msld_prm_finds_bonded_terms` (strict `xfail` when pandas >= 3) |
| #14 | Speed up the core search for large molecules (junction tree and other ideas) | enhancement | **postponed 2026-10-21** by the owner until the whole pipeline is re-designed and re-implemented | none yet |

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
- **Workaround, used until the fix and removed with #12:** a GUI option "Remove charges from the
  reference" (RDKit's `Uncharger`). It is not needed any more.
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

### #11 The symmetric match that is embedded is arbitrary (fixed)

- **Status:** fixed 2026-10-16 (`56f697c`, closes #11). When the core fits the analog in more than
  one distinct way (a ring that can flip), each way is embedded and the one whose atoms outside the
  core lie nearest the reference's atoms outside the core is kept. On the 33 public demo analogs
  the first match had the meta-substituted aniline ring flipped for 28 (the analog's halogen
  3 to 6 angstrom from the reference's iodine; now 0.2 to 1.2). Where the reference has nothing to
  overlap, the first match is kept (a true symmetric structure: any match is fine, owner,
  2026-10-16). Not decided by the pocket yet: see F-2.
- **Not checked:** the owner's own series beyond the three analogs of `panel_a`; cores with more
  than 32 distinct orientations (the search stops at `McsOptions.max_orientations`).

### #2 Stereochemistry (closed 2026-10-21 for alignment)

- **Status:** the alignment step is done 2026-10-16 (`b51a361`, `a44e480`). The owner limited the
  issue to alignment on 2026-10-21 and closed it; the other stages (`decompose.py`, the writers) are
  tested after the re-design.
- **What alignment does now:** a specified centre is kept as given. An unspecified centre inside the
  core, on an atom where the reference has a stereocentre (read from its 3D coordinates), is set the
  way that fits the reference. Any other unspecified centre or double bond is set each way and the
  pose with the lowest UFF energy is kept, flagged in the GUI as a guess. A specified centre whose
  inverted configuration fits the reference clearly better is kept and flagged as opposing the
  reference. Owner's choices, 2026-10-16.
- **Core geometry** (`f43e56d`, 2026-10-16): the owner saw the reference's piperidine chairs come
  out inverted or twisted in the analogs. That is a conformation, not a configuration (the
  reference has no stereocentre there), so chirality tags could not fix it; the cause was RDKit's
  default embedding with a coordinate map, which started the core about 1 angstrom off, and the
  soft tether, which could not flip a ring back. The core atoms are now put exactly on the
  reference's positions and kept fixed while the rest is optimized. When that would invert a centre
  the analog specifies (the opposite configuration to the reference), the pose is rejected and the
  core is held by a stiff tether, then by RDKit's default start and soft tether, so the analog keeps
  its configuration and is flagged. On the owner's 22 analogs: core RMSD 0.000 and no ring torsion
  off (before: 15 had a ring more than 30 degrees off). Checked by the owner in the GUI.
- **Evidence** (`[ran]`): in a toy test the unspecified centre came out R or S with the seed and now
  always matches the reference; in the owner's `stereoisomer` set the specified cis and trans ring
  isomers survived in every run (they already did), and the unspecified centres are settled the
  same way (not at one seed in five, where UFF energies differ by about 8 kcal/mol from conformer
  noise).
- **Not checked:** stereo through the written files (sdf and mol2 carry the 3D coordinates; the
  mol2 writer has no chirality tags), through `decompose.py` and the CHARMM writers; double-bond
  stereo is enumerated but not reported; the energy rule uses one pose per isomer, so a close call
  is noise.

## Planned features

Unordered. Nothing here is scheduled.

| ID | Feature | Why | Depends on | Issue | Notes |
| --- | --- | --- | --- | --- | --- |
| F-1 | Alignment GUI, first slice | demo to the lab; check alignments by eye | none | #3 to #8 | Built 2026-10-14, see "Done"; the 3D view was checked in Firefox with software WebGL on 2026-10-14. Issues #4 and #8 were rewritten on 2026-10-14 to match what was built; #5 to #7 still describe the first plan and are partly superseded |
| F-2 | Show the protein pocket in the GUI for clash checks | the workshop notebook asks for it | F-1 | none | load a PDB for display only. **When the pocket is added, revisit how a flipped ring is decided**: today it is chosen by the overlap of the analog's substituents with the reference ligand's (#11); the owner wants the pocket's free space to decide it then (2026-10-16) |
| F-3 | Export Sybyl-typed mol2 per analog from the GUI | what the notebook hands to CGenFF | F-1 | #7 | **Done 2026-10-14**, see "Done" |
| F-4 | Several users or sessions on one GUI server | today one session per server (G16) | F-1 | none | needs session ids in the API |
| F-5 | Take already-parameterized sets (mol2 + rtf + prm) in the GUI | reuse existing CGenFF output | F-1 | none | |
| F-6 | Carry the GUI's aligned set into `define_common_core` | one tool from alignment to the build | F-1, design decisions | none | blocked on the learning plan |
| F-7 | Fewer saved files: keep intermediate stages in memory | owner's goal (handoff section 2, 2026-10-12) | learning plan (Q13): which files CHARMM and ALF really read | none | |
| F-8 | A web GUI for the whole preparation, not only alignment | owner's goal | F-6, F-7 | none | |
| F-9 | Compare the old and new pipelines on real MSLD runs | decide the design with evidence | learning plan, CHARMM basics from colleagues | none | |
| F-11 | Option to keep the original atom names in the mol2 export | the owner's choice to rename is deliberate (CGenFF), but some uses need the original names | F-3 | #10 | **Done 2026-10-28**, see "Done". Default stays "new unique names"; kept names that are duplicated (ignoring case) or longer than four characters are refused with the atoms named. Only a reference read from mol2 has names to keep; analogs from SMILES are named as before (owner, 2026-10-28: as the issue says, no inheriting of reference names by the core) |
| F-12 | Manual tuning of an aligned analog in the GUI: rotate a bond to fix a small misalignment | the constrained embed can leave a free part in a slightly wrong pose | F-1 | none yet | 3Dmol.js has no bond-rotation tool (checked 2026-10-14: no torsion or dihedral code in the bundled 2.4.2), only atom picking (`setClickable`, `setHoverable`) and `setCoordinates`. Likely design: pick a bond in the page, send it with an angle to the server, rotate the atoms on the side that holds no matched-core atom with RDKit (`rdMolTransforms`), return the new coordinates, show clashes. The core must stay fixed. Not decided |
| F-14 | **Speed up the core search for large molecules; postponed 2026-10-21 until after the re-design** (fix #15 first, then profile) (junction tree comparison as in Jin et al., `seedSmarts`, scaffold substructure first, caching) | aligning large molecules is slow; the MCS search is the likely cause (profile first: the orientation and stereo search embed several candidates too) | none | #14 | risks of the junction-tree idea are listed in #14: the tree is not unique, labels must not carry bond orders, atoms inside matched rings still need mapping, results may differ from `FindMCS` |
| F-10 | Learning plan for the original `msld-py-prep` | prerequisite of the roadmap | none | none | **started 2026-10-21**; lessons live in the owner's Obsidian vault (`/home/ziqiaoxu/obsidian/msld-prep/`); the design notes are there and in [DESIGN_NOTES.md](DESIGN_NOTES.md). Lessons 1 to 3 done (lesson 3 on 2026-10-28: CRN) |

## Deliberate behaviors

Things that look like gaps but are choices. Do not "fix" them; change them only by an item above.

- **mol2 export renames every atom and writes no charges** (owner, 2026-10-14). The GUI's mol2
  files are written for CGenFF. CGenFF keeps only unique names of at most four characters
  (it silently renumbers others) and computes its own charges, so the writer
  (`sybyl_typer.assign_unique_atom_names`, `format_mol2_text`) gives every atom a fresh unique
  name, the reference included, and writes `NO_CHARGES`. The reference's names from its own file
  are therefore not kept unless "Keep the reference's atom names" is chosen (F-11, #10). The sdf export keeps the
  structure as it is and adds each analog's RMSD and seed.
- **A realigned analog loses its accept or reject mark**, because the pose is new.
- **The reference is used with hydrogens**, added with coordinates if the file has none. Not true
  for a mol2 reference without hydrogens: none are added (#16, a bug).
- **A symmetric core is placed by the reference's substituents** (#11): among the equivalent ways
  the core fits the analog, the one whose substituents lie nearest the reference's is embedded.
  Pocket space is meant to take over this decision later (F-2).
- **Unspecified stereo is settled, not left to chance** (#2): by the reference where it has a
  stereocentre on the mapped atom, otherwise by energy, and the page says which ("check stereo").
  A specified configuration is never overruled, even against the reference.
- **The core search stops after 10 seconds by default** (#13) and uses the best structure found,
  noted in the details.
- **The core sits exactly on the reference** (#2, `f43e56d`), so the RMSD of the core is 0.000 for
  almost every analog. A non-zero RMSD means a fallback was needed, in practice a stereocentre whose
  configuration conflicts with the reference's.
- **The GUI's text is in American English** (owner, 2026-10-16): "stereocenter", "color".

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
| 2026-10-14 | GUI option to remove charges from the reference with `Uncharger` (the notebook's step); **removed 2026-10-16** (#12, `d013379`) because charges no longer matter for the match | `msldprep/gui/service.py::neutralize` |
| 2026-10-14 | GUI toggle to color the reference by element; gray when off | `msldprep/gui/static/app.js` |
| 2026-10-14 | Alignment matches by element and ignores charges and hydrogen counts; option to ignore bond orders so tautomers and protonation variants match in full (on by default in the GUI) | #9; `msldprep/align.py::align_to_reference`, `McsOptions.ignore_bond_orders` |
| 2026-10-14 | Collapsible 2D drawing of the reference in panel 1; with an analog selected, the atoms it matched are marked | `msldprep/gui/static/app.js::renderReference2D`, `ligand_view()["reference_svg"]` |
| 2026-10-14 | GUI download as sdf or mol2 (the notebook's Sybyl typing and mol2 writer), one file or one file per structure (zip), with or without the reference | F-3, #7; `msldprep/gui/service.py::AlignmentSession.export` |
| 2026-10-16 | Core search limited to 10 seconds by default; the details say when the limit was hit | #13, `9002a9e` |
| 2026-10-16 | Orientation of a symmetric core chosen by the overlap of substituents with the reference's | #11, `56f697c`, `0b6af57` |
| 2026-10-16 | Stereo settled during alignment (by the reference, else by energy) and reported in the GUI; the reference's stereocentres shown | #2 (part), `b51a361`, `a44e480` |
| 2026-10-16 | Issues #5, #6, #7 closed as built; GUI scope limited to structure alignment with MCS (owner) | GitHub |
| 2026-10-16 | Core placed exactly on the reference so rings keep its puckers, with a fallback that keeps an opposite configuration | #2 (part), `f43e56d` |
| 2026-10-28 | Opt-in report of protonation changes in the core (off by default), in alignment and the GUI | #1; `McsOptions.check_protonation` |
| 2026-10-28 | mol2 export option to keep the reference's atom names (off by default) | F-11, #10; `sybyl_typer.assign_unique_atom_names(keep_existing=)`, `AlignmentSession.export(keep_names=)` |
| 2026-10-16 | `python -m msldprep.gui --reload` for development; page files served with `no-cache`, so edits show on a refresh | `decffa7` |
| 2026-10-16 | A second click on the selected analog deselects it (reference alone in the view) | `3a38fac` |
| 2026-10-16 | Hidden file inputs kept inside their buttons; they had covered the "Settings" and "2D structure" toggles in Safari when the left panel scrolled | `af93fd5` |
| 2026-10-16 | GUI text in American English ("stereocenter") | `f1f9edd` |
| 2026-10-21 | GUI README: walk-through, contents of the downloads, troubleshooting | #8 (part), `msldprep/gui/README.md`, uncommitted when written |
