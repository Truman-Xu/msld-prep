# Design notes for the msld-prep re-design

Design decisions for the re-design of msld-prep, collected while the owner learns the legacy
`msld-py-prep` (handoff Q13). Each note says what was observed, the owner's verdict with its date,
and the reason. "Open" means no verdict yet.

**Two copies, kept identical in content:** this file, and `design-notes.md` in the owner's Obsidian
vault (`/home/ziqiaoxu/obsidian/msld-prep/`), where the lessons it refers to also live. Change both
in the same session. In the vault copy the lesson names are links.

## From lesson 1

### (a) Settings are data, not edits to the program

**Observed:** to configure a legacy run you open `msld_py_prep.py` and change Python variables in
it (`sysname`, `molfile`, `cgenff`, `inFrag`, `AnCore`, the cutoff). So each project needs its own
edited copy of the tool. A build doesn't record which settings made it, apart from that copy.
Nothing checks the values before the run, and a typo is a Python error in the middle of the run.

**Proposal:** settings come from outside the program, as function arguments or a config file, and
are checked before anything runs. The settings used are written into the build. msld-prep already
goes part of the way (`configs/*.yaml`, `msld_prep_cli.py`).

**Verdict 2026-10-21: adopted.** Settings come from a config file. Reason: a project's settings
become data that can be checked, kept with the build, and reused, and the tool itself stays
unedited.

### (b) A person checks the core, in the GUI or in a script

**Verdict 2026-10-21:** agreed. The GUI is the main checkpoint. There must also be a way to do the
whole preparation from a script. Running the same script twice, as the legacy code does, breaks
the flow and is counter-intuitive.

**Verdict 2026-10-21: the five proposals below are adopted as the design** of the script path.
Reason: one object model serves the script, the GUI and the command line, review happens in
memory, and a rerun is reproducible.

**The script path:**

1. **One session object, run top to bottom.** Each step returns an object you can look at in a
   notebook, so the review happens in memory instead of between two runs:

   ```python
   from msldprep import Series

   series = Series.from_files("1.mol2", ["2.mol2", "3.mol2"])   # or Series.from_smiles(...)
   series.align()
   core = series.find_core()
   core.show()                     # 2D grid with the core highlighted, 3D overlay
   core.report()                   # sites, nsubs, flags that need a decision
   series.renormalize().report()   # charge changes per atom and fragment
   series.write("build.faah")
   ```

2. **Corrections are written as code or config, not as edits to an output file.** For example
   `core.move_to_fragment("1", ["C013"])` or an `overrides:` list in the config, matched by atom
   name. Rerunning from the start then gives the same build, and the corrections are recorded.
3. **The run stops only when a decision is needed.** It stops for a one-substituent site, a large
   charge change, a poor core fit, or an ambiguous atom, with a report that says what to decide.
   A `strict=False` lets it continue with warnings.
4. **The checkpoint file becomes optional.** `core.save("core.json")` and `Core.load(...)`, for
   stopping and coming back or handing over to the GUI. It is not the only way to pass data
   between steps.
5. **One implementation behind three faces.** The script API, the GUI (a view on the same
   `Series`) and the command line (`msldprep build config.yaml`) all call the same objects.

### (c) "Common core" needs one explicit definition

**Verdict 2026-10-21:** agreed. Lesson 2 shows the two pipelines define it differently; see (e).

### (d) Fewer saved files

**Verdict 2026-10-21:** agreed. Combine the machine-read files as much as possible. Keep
human-readable files to the minimum that really needs human judgement. `translation.txt` should
not be needed if the data and code are done right.

**Constraint, agreed 2026-10-21** (this is what "implemented right" means): every fragment is patched into the same residue as the core, so atom
names must be unique across the core and all fragments. Renaming can't always be avoided, since
two ligands both have a `C009`. What can be avoided is a separate file of renames for people: the
mapping can live in the data and be shown on demand (GUI, `core.report()`).

## From lesson 2

### (e) Which rule defines a core atom?

Legacy: same CGenFF type, same neighbor types, and (for ties) the same position, in every ligand.
msld-prep: same CGenFF type in a whole-ring graph MCS. The legacy rule moves a ring carbon whose
substituent changes into the fragment (its charge differs a lot, −0.115 vs +0.063). msld-prep keeps
it in the core and leaves the difference to CRN. Lesson 3 should measure what each costs in charge
changes.

**Verdict 2026-10-29: adopted, the legacy rule.** An atom is a core atom only if its CGenFF type and
the types of all its neighbors (hydrogens included) are the same in every ligand; an atom whose
neighbor has a different type is dropped. Reason: an atom whose surroundings change also changes its
charge, so the rule keeps such atoms out of the core. In lesson 3 the legacy core needed almost no
charge correction (2 of 31 positions differed).

**How:** on top of RDKit's `FindMCS`, not as a port of the legacy algorithm. Each heavy atom is
labeled with its pattern (own type plus the sorted types of its neighbors), and `FindMCS` matches the
labels on heavy-atom copies (`CompareIsotopes`, any bond matches any bond, `ringMatchesRingOnly=True`,
`completeRingsOnly=False` so a ring can be cut as legacy cuts it). Where the core fits a ligand in
several ways, the placement that agrees best in 3D over all ligands is taken. Prototype
`scripts/pattern_mcs.py` in the owner's vault (2026-10-29, on public FAAH copies): ligands 13, 31, 1
give legacy's heavy-atom core exactly in 0.01 s; all 32 give legacy's core plus the 5 positions
legacy dropped by distance, in 0.02 s. One difference from legacy: `FindMCS` keeps one connected
core, so an atom that changes in the middle of the molecule cuts off the smaller side, where legacy
could keep both.

**Implemented 2026-10-29** (`fd6cbaa`), with (i). On the ChEMBL test set the core and the sites equal
legacy's. Sites are now legacy's too: a connected piece of non-core atoms, so a cut ring is one
site, and separate substituents on neighboring core atoms are separate sites (confirmed by the
owner 2026-10-29: the correct way to handle a cut ring).

### (f) Positions: a check, not a definition

The legacy core depends on the input conformers: a rotated phenyl in 17 of 32 FAAH poses makes a
one-substituent site. Proposal: decide the core by chemistry, and use positions to flag poor fits
for the person to look at.

**Verdict 2026-10-29: adopted.** Atoms that match by (e) stay in the core even when they lie more
than 1.1 Å from the reference in an input pose; the pose is flagged for the person to look at.
Positions are used only to choose among equivalent placements. Reason: on all 32 FAAH ligands the
legacy distance cut removes 5 positions that match by chemistry (one ortho carbon, up to 2.34 Å off
in 6 ligands; four CH of a phenyl that is rotated in many input poses, up to 2.23 Å), and
msld-prep's alignment puts the core exactly on the reference anyway.

**Implemented 2026-10-29** (`b003554`): `decompose.choose_core_placements` and `far_core_atoms`;
the flags are in `CoreDecomposition.far_core_atoms` and the log, not yet in the checkpoint file.

### (g) Flag what needs judgement

Adding ligands can shrink the core. A site with one substituent, or a core that shrinks when one
ligand is added, should be reported with the ligand that caused it.

**Verdict 2026-10-29: adopted, and extended.** These cases are flagged. What is really wanted is an
algorithm that **proposes several common cores for a pool of diverse analogs**, each with its own
member analogs, under a **cutoff on the percent of each member's atoms that the proposed core
covers**. An analog joins a core only if that core covers at least the cutoff of it; an analog no
core covers well enough gets a core of its own or is reported as left out. Reason: one core forced
over a diverse pool shrinks to whatever survives everywhere, and the coverage cutoff states directly
what makes a core worth sharing.

msld-prep today (`grouping.group_ligands_by_similarity`) splits the pool **before** any core is
found, by whole-molecule similarity (Morgan radius-2 fingerprints, Tanimoto ≥ 0.7): the reference's
group first, then groups seeded greedily among the rest. Fingerprint similarity only stands in for
what the owner wants; the cutoff should be measured on the core itself. Fingerprints may still serve
as a fast first guess.

**Decided 2026-10-29:**

- **Coverage counts heavy atoms only:** the core's heavy atoms divided by the analog's heavy atoms.
- **Cutoff 0.7 by default.** Alternatively the person gives the **number of common cores** wanted,
  instead of a cutoff.
- **An analog may appear in several proposals.** An option makes the proposals **unique**: each
  analog stays only in the proposal whose core covers it best (highest percent overlap).
- **Search by growing:** a proposal grows by adding analogs while every member stays at or above the
  cutoff; when no analog can be added, the next proposal starts.
- The cores follow the pattern rule of (e).

**Left for the implementation:** how a growing proposal picks its seed and the next analog to add;
how growing meets a requested number of cores (which proposals are kept, and what happens to
analogs none of them covers); how the unique option breaks a tie in coverage.

### (h) No keyboard prompts or `quit()` in library code

Legacy stops and asks with `input()` for an ambiguous atom, or exits. In the re-design an ambiguity
is part of the report (g), decided in the GUI or by an override (b2).

**Verdict 2026-10-29: adopted.** No `input()`, `quit()` or `sys.exit()` in library code. The GUI may
prompt for a decision. The script path follows (b): the run goes on until a decision is really
needed, then stops with a report of what to decide (b3), and the decision is given as code or config
and recorded (b2). A `strict=False` run continues with warnings. Library code raises a typed error
for that stop, never ends the process. Reason: a prompt or an exit inside the library breaks the
GUI, scripts and notebooks alike, and leaves no record of what was decided.

### (i) Which hydrogens belong to the core? (raised by #15)

Today msld-prep decides it by accident: the shared-core search keeps every hydrogen (#15), so a
hydrogen is core only when its CGenFF type matches in every ligand. A hydrogen whose type differs
(for example HGR61 against HGR62 on the same core carbon) becomes an **alchemical site of its own**:
on the ChEMBL test set, two of the three sites contain nothing but hydrogens. Searching heavy atoms
only, with nothing else changed, makes it worse: every core hydrogen becomes a fragment atom.

Proposal: search heavy atoms only, then apply an explicit rule. The hydrogens of a core heavy atom
are core when that atom carries the same number of hydrogens, of the same type, in every ligand.
Otherwise the atom is an anchor and its hydrogens join the fragment attached there. No
hydrogen-only sites. Whether a type mismatch should instead push the heavy atom out of the core,
as the legacy neighbor-type rule does, belongs to (e).

**Settled 2026-10-29 as a consequence of (e).** The search runs on heavy atoms only, and a heavy
atom's pattern includes its hydrogens' types. So a core heavy atom carries the same number of
hydrogens, of the same types, in every ligand, and **its hydrogens are core**. A heavy atom whose
hydrogens differ in number or type does not match, so it and its hydrogens go to the fragment. No
hydrogen-only sites. This also gives the fix of #15 its direction. **Implemented 2026-10-29**
(`fd6cbaa`, closes #15).

## From lesson 3

### (j) Measure charge changes per atom, not as a percent of a fragment's total

Both pipelines flag a fragment when its total changes by more than 5% of its original total. For
a near-neutral fragment that means nothing: legacy prints "−44.43% ** CHECK" for 0.0013 e, and
msld-prep warns "−833%" for 0.008 e. Proposal: report the largest change on any atom and the RMS
change per ligand, in e, with a threshold the owner sets (the paper's figure of about 0.001 e RMS
is a candidate, to be re-checked).

**Verdict 2026-10-30: adopted (#20).** The percent check goes. Every run reports, in e: **Q-RMSD per ligand**
(RMS of new minus CGenFF charge over the ligand's atoms), **the largest change on any atom**, and **the
charge spread at each core position**. A ligand above **0.001 e** is flagged. Nothing stops by default;
an optional config limit can stop a run. Reason: the percent check looks at the wrong thing. It raises
false alarms (FAAH 1, 13, 31 with anchors in the core: 2 of 3 fragments flagged, for at most 0.0013 e,
with every ligand under 0.001 e) and misses real changes (ChEMBL 34 ligands with anchors out: 9 of 68
fragments flagged, while all 34 ligands are above 0.001 e and one core atom moves 0.055 e). It does not
see what core averaging does to single atoms.

**Why a flag and not a stop.** The paper (`papers/nihms-1851979.pdf`, read in full 2026-10-30) says that
above 0.001 e "large deviations in both Δμ and ΔΔG can occur", while most ΔΔG stay under 0.5 kcal/mol
and some outliers exceed 1.0. Its Table 1 gives mean Q-RMSD 0.002 e for FAAH, 0.003 e for ChEMBL, 0.022 e
for 12fused (its worst series, 0.45 kcal/mol) and 0.007 e over all series. So 0.001 e is below what the
paper treats as routine; a stop there would halt most real series. The paper recommends bookending
whenever CRN is used. It also reports dipole changes (Δμ), which this verdict leaves out.

**Checked against the SI 2026-10-30** (`papers/NIHMS1851979-supplement-Supplementary_Information.pdf`,
Tables S1 to S6, 191 molecules). **Q-RMSD covers all of a ligand's atoms, lone pairs included.** The
charges in legacy's own `examples/1_FAAH/build.name`, against each ligand's rtf, reproduce all 32 values of
Table S6 to six decimals; leaving out lone pairs or hydrogens does not. This is also a reference to test
against in lesson 6.

What the flag does on those 191 molecules, with hydration ΔΔG between original and renormalized charges:

| Q-RMSD (e) | Molecules | Median abs(ΔΔG) | Above 0.5 kcal/mol | Above 1.0 kcal/mol |
| --- | --- | --- | --- | --- |
| up to 0.001 | 54 | 0.009 | 0 | 0 |
| 0.001 to 0.002 | 60 | 0.091 | 1 | 0 |
| 0.002 to 0.005 | 13 | 0.120 | 1 | 0 |
| 0.005 to 0.01 | 17 | 0.094 | 2 | 1 |
| 0.01 to 0.02 | 29 | 0.286 | 9 | 1 |
| above 0.02 | 18 | 1.183 | 14 | 13 |

So 0.001 e is a sound line for "nothing to do": no molecule under it exceeds 0.173 kcal/mol, and it misses
none of the 27 above 0.5. But it fires on 137 of 191 (72%), and only 20% of those exceed 0.5.

**Verdict 2026-10-30: add a second, stronger flag (#20).** 0.001 e stays as the line below which nothing needs
doing. A ligand gets the stronger flag when **Q-RMSD is above 0.01 e or abs(Δμ) is above 0.5 D**, so the run
also computes the dipole change, from the coordinates. Neither flag stops a run by default. Above 0.01 e alone
fires on 47 and catches 23 of the 27. Adding abs(Δμ) above 0.5 D fires on 52 and catches 25 of the 27, including
all 15 above 1.0. **The thresholds are provisional; re-check them in lesson 6.** 41 of the 47 above 0.01 e are
from one series (12fused), outside it only 2 of 6 exceed 0.5, and these are hydration results, not binding.

Measurements in (j) to (n) are from 2026-10-30 at `61030bc`, in memory (so #18 does not affect them), on
input poses, with CGenFF charges as the baseline; no CHARMM run.

### (k) Don't cut a charge-coupled ring into separate sites

On FAAH 13, 31, 1, legacy makes one site of the substituted phenyl. msld-prep keeps the whole ring
in the core, then moves the substituted ring carbons out as anchors, which gives **two sites**.
CGenFF shifts charge across the ring (ligand 13: +0.045 e in one half, −0.047 e in the other), and
each half has to be forced to its own site average. Result: largest atom change 0.016 e against
0.002 e, RMS 0.0034 e against 0.00024 e (14×). Proposal: a site should not split atoms whose
charges move together, and the CRN cost should be part of choosing the core. This belongs with
(e) and (i).

**Correction 2026-10-30:** the numbers above predate rule (e). With (e), FAAH 1, 13, 31 is one site, so
the split is gone. With anchors in the core msld-prep now gives legacy's figures (RMS 0.00023 e, largest
0.002 e).

**Verdict 2026-10-30: the cut-ring part is closed by (e); the rest is a report, not an optimizer (#21).** The
core is not chosen by CRN cost automatically. The report names the core atoms with the widest charge
spread and what moving each one would do, and the person decides through an override (b2). Reason: most
of the change now comes from averaging core charges (ChEMBL 34 ligands: RMS 0.0083 e on core atoms,
0.0012 e on fragment atoms), because the pattern rule looks only one bond out and an atom can pass it with
a varying charge. Moving that one atom can help a lot: on ChEMBL 34 ligands, moving `N011` out takes
ligands above 0.001 e from 34 to 6 or 9 (depending on the site, whose bonding was not checked) and mean
Q-RMSD from 0.0023 e to 0.0014 or 0.0016 e. The paper describes the same case (12fused: "sufficient atom
type overlap" but "poorer electrostatic overlap") and the same remedy (`InFrag`).

### (l) The site average: per substituent or per ligand?

Legacy averages over the unique fragments at a site (each substituent counts once). msld-prep
averages over ligands, so a substituent shared by ten ligands counts ten times, and adding a
duplicate ligand changes every substituent's charges. Since MSLD simulates each unique substituent
once, per-substituent weighting matches the simulation. Tied to whether msld-prep deduplicates
fragments at all (lesson 6).

**Verdict 2026-10-30: deduplicate, and average per substituent (#22).** Each distinct substituent at a site
appears once and counts once, as in legacy and the paper ("fragments at each site are filtered to remove
possible duplicates"; the target is the average "of all substituents at a single perturbation site").
Reason: the charges then depend on which substituents exist, not on which combinations happen to be
listed, and MSLD samples every combination anyway. Today msld-prep writes one block per ligand: ChEMBL 34
ligands get `NSUBS 34 34` where only 14 and 27 substituents are distinct, so 68 blocks instead of 41. The
weighting itself costs almost nothing there: a core correction of 0.0018 e in total, mean Q-RMSD 0.00236 e
against 0.00234 e. This changes what is written (`NSUBS`, patches, checkpoint), and the link from each
ligand to its substituent at each site must be kept for reading results.

**The identity test must compare the fragment graph with its attachment points:** atom types, bonds, and
which core atom each attached atom bonds to. Two simpler tests failed on FAAH: types with distance to the
core merged para-CF₃ with meta-CF₃, and so did a graph without attachment points, because a cut ring is
the same chain either way. The wrong merge put 0.087 e on one atom. On ChEMBL 34 ligands, the substituents
the graph test calls identical have identical charges.

### (m) How to spread a correction over a fragment's atoms

Legacy spreads it evenly; msld-prep in proportion to |charge|, which loads highly charged atoms
(a CF₃ carbon takes 3.4× the even share). Both meet the constraints, and the amount to spread is
set by the decomposition, so this is the smaller question. To decide: even, proportional, or a
least-squares fit that keeps each atom closest to its CGenFF charge.

**Verdict 2026-10-30: spread evenly (#23).** One exact pass, not legacy's ±1e-6 loop. Lone pairs still do not
move. Reason: for a fixed total, plain least squares gives equal shares, so "even" and "least squares" are
the same option, and proportional is the weighted version with weights 1/|q|. Even was as good or slightly
better everywhere (mean Q-RMSD on ChEMBL 34 ligands 0.00234 e against 0.00238 e; on 4 ligands 0.00492 e
against 0.00507 e), and it does not load highly charged atoms. Fitting the site totals by least squares as
well gave no gain (0.00235 e against 0.00236 e), so that is dropped. This changes the charges written.

## From lesson 4

### (n) A site attached at two points: anchors in the core or in the fragments?

On FAAH 1, 13, 31 the fragment is attached to the core at two atoms (`C008`, `C013`; legacy's anchor
`DUM`). Legacy moves an anchor into the fragments only when it is a single atom
(`msld_crn.py:253`), so these stay in the core (core 31 atoms, ligand 1's fragment 8). msld-prep
moves every anchor out (`apply_manual_overrides`; core 28, fragment 11). Both can be overridden
(`AnCore` / `an_core`). Moving anchors out gives each substituent its own charges at the attachment;
keeping them makes the shared core bigger. Here their charges barely differ between ligands
(lesson 3: `C008` 0.008/0.005/0.008, `C013` −0.121/−0.120/−0.120). On all 32 FAAH ligands they differ
more (spread 0.021 e and 0.011 e).

**Verdict 2026-10-30: keep the legacy behavior, anchors in (#24).** A site attached at two or more points keeps
its anchors in the core. As in legacy (`msld_crn.py:250-254`), a single-point anchor still moves into the
fragments. Reason: the owner prefers fragments smaller to begin with. **There must be an easy option or
function to move the anchors out.** Today's `apply_manual_overrides` moves every anchor out by default,
so this default changes. **Confirmed by the owner 2026-10-30:** this exact rule, single-point anchors out
included, is the default.

The cost is known and accepted, and it is why the easy option is wanted. Mean Q-RMSD, per substituent and
spread evenly:

| Series | Legacy rule | Anchors out | Fragment atoms, legacy rule / out |
| --- | --- | --- | --- |
| FAAH, 32 ligands | 0.00259 e, 28 of 32 above 0.001 e | 0, no charge changes | 8-20 / 11-23 |
| ChEMBL, 4 ligands | 0.00611 e, 4 of 4 | 0.00492 e, 4 of 4 | 6-20 / 9-22 |
| ChEMBL, 34 ligands | 0.00285 e, 34 of 34 | 0.00236 e, 34 of 34 | 8-22 / 10-24 |

With anchors out, the FAAH core charges are identical across all 32 ligands and every substituent keeps its
CGenFF charges. The legacy-rule figures are close to the paper's Table 1 (FAAH 0.002 e, ChEMBL 0.003 e).
What the larger fragments do to sampling was not measured.

## From the review of lesson 3

### (o) A core in two pieces: keep the second piece?

**Observed 2026-10-30** on the 31 ligands of `examples/2_ChEMBL/mol_list.txt`. Legacy's core is
`C017`, `N009`, `C010`, `H023`, `N011` with `NSUBS 19 24`: two sites, so 456 combinations from 43
substituents. msld-prep's core is the same minus `C017`, with `NSUBS 31`: one site, so 31 combinations from
31 substituents. All five atoms pass rule (e) in all 31 ligands, but `C017` is not bonded to the other four.
`FindMCS` keeps one connected core, the difference (e) already notes, and with `C017` gone the rest of the
molecule is one connected piece, so one site. Legacy's lone core atom is what splits it into two sites.

This costs msld-prep the multisite gain on this series. Whether a lone core atom is good is less clear.
Legacy's own log flags it: "Large RMSD for CORE atom C005 - RMSD > 0.8 (1.532)", and `C005` is ligand 9538's
name for this position, assuming the log uses the first ligand's names. So the atom that gives the second
site fits poorly in 3D; its charge fit and its effect on sampling are not known. Legacy's figures were read
from its `MCS_for_MSLD.txt` and log, not re-run. Eight of these ligands need #19's workaround to load.

**Verdict 2026-10-30: keep the connected core for now, and revisit in lesson 6.** Reason: the multisite gain
is large, but there is no evidence yet on what a lone, poorly fitting core atom does to charges or sampling.
Lesson 6 runs old against new on the full series, which can supply it.
