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
hydrogen-only sites. This also gives the fix of #15 its direction.

## From lesson 3

### (j) Measure charge changes per atom, not as a percent of a fragment's total

Both pipelines flag a fragment when its total changes by more than 5% of its original total. For
a near-neutral fragment that means nothing: legacy prints "−44.43% ** CHECK" for 0.0013 e, and
msld-prep warns "−833%" for 0.008 e. Proposal: report the largest change on any atom and the RMS
change per ligand, in e, with a threshold the owner sets (the paper's figure of about 0.001 e RMS
is a candidate, to be re-checked). **Open.**

### (k) Don't cut a charge-coupled ring into separate sites

On FAAH 13, 31, 1, legacy makes one site of the substituted phenyl. msld-prep keeps the whole ring
in the core, then moves the substituted ring carbons out as anchors, which gives **two sites**.
CGenFF shifts charge across the ring (ligand 13: +0.045 e in one half, −0.047 e in the other), and
each half has to be forced to its own site average. Result: largest atom change 0.016 e against
0.002 e, RMS 0.0034 e against 0.00024 e (14×). Proposal: a site should not split atoms whose
charges move together, and the CRN cost should be part of choosing the core. This belongs with
(e) and (i). **Open.**

### (l) The site average: per substituent or per ligand?

Legacy averages over the unique fragments at a site (each substituent counts once). msld-prep
averages over ligands, so a substituent shared by ten ligands counts ten times, and adding a
duplicate ligand changes every substituent's charges. Since MSLD simulates each unique substituent
once, per-substituent weighting matches the simulation. Tied to whether msld-prep deduplicates
fragments at all (lesson 6). **Open.**

### (m) How to spread a correction over a fragment's atoms

Legacy spreads it evenly; msld-prep in proportion to |charge|, which loads highly charged atoms
(a CF₃ carbon takes 3.4× the even share). Both meet the constraints, and the amount to spread is
set by the decomposition, so this is the smaller question. To decide: even, proportional, or a
least-squares fit that keeps each atom closest to its CGenFF charge. **Open.**
