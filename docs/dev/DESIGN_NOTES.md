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
changes. **Open.**

### (f) Positions: a check, not a definition

The legacy core depends on the input conformers: a rotated phenyl in 17 of 32 FAAH poses makes a
one-substituent site. Proposal: decide the core by chemistry, and use positions to flag poor fits
for the person to look at. **Open.**

### (g) Flag what needs judgement

Adding ligands can shrink the core. A site with one substituent, or a core that shrinks when one
ligand is added, should be reported with the ligand that caused it. **Open.**

### (h) No keyboard prompts or `quit()` in library code

Legacy stops and asks with `input()` for an ambiguous atom, or exits. In the re-design an ambiguity
is part of the report (g), decided in the GUI or by an override (b2). **Open.**

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
as the legacy neighbor-type rule does, belongs to (e). **Open.**

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
