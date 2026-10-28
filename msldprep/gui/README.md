# Alignment GUI

A browser page for one step of the preparation: placing each analog onto a reference ligand in 3D.
You load a reference with 3D coordinates, paste the analogs as SMILES, and the page aligns each
analog by its common structure with the reference. You look at every result, accept or reject it,
and download the structures.

What the page does **not** do: atom typing with CGenFF, charges, or the common core of the whole
series. Those come later in the pipeline. The page only decides where the atoms of each analog sit.

## Start it

The server runs on the cluster and you reach it from your own computer through an SSH tunnel.

```bash
pip install fastapi uvicorn      # only if the environment lacks them
python -m msldprep.gui
```

It prints two things: an `ssh -L ...` command and an address with an access token, like this:

```text
Alignment GUI is starting. From your own computer, open a tunnel:
  ssh -L 40123:127.0.0.1:40123 you@gollum
then open this address in your browser:
  http://127.0.0.1:40123/?token=Qx...
```

Run the `ssh` command in a terminal **on your own computer** and leave it open. Then open the
address in your browser. Ctrl-C in the server's terminal stops the server.

Options:

| Option | Default | Use |
| --- | --- | --- |
| `--port N` | a free port | pick the port yourself |
| `--token T` | a random token | a fixed token, for example to keep a bookmark working across restarts. Pick a long one: other users of the login node can reach the port |
| `--host H` | `127.0.0.1` | listen on another interface. It prints a warning: anyone who can reach that port and has the token can use the page |
| `--reload` | off | for development: restarts when a Python file of `msldprep` changes, with the same port and token |

Changes to the page's own files (`msldprep/gui/static/`) show when you refresh the browser.

## Walk-through with the demo

The demo uses public data: the ChEMBL series of `examples/2_ChEMBL`.

1. **Reference.** Choose `examples/2_ChEMBL/11896.mol2`. Below the button appear its atom count
   (with hydrogens), its net charge, any charged atoms and its stereocenters. "2D structure" opens
   a drawing of it.
2. **Analogs.** Open `examples/gui_demo/analogs.txt`, or paste lines of the form `SMILES name`.
   A table with a name column also works. Press **Use these analogs**. The demo has 35 lines.
3. **Align.** Press **Align all**. The list fills in as each analog finishes. With the demo, 33
   are aligned and two fail on purpose, to show how a failure looks:
   - `ethanol`: "FindMCS found no shared substructure between the two molecules." Ethanol has no
     ring, and by default only whole rings count as common structure.
   - `bad_smiles`: "Could not parse SMILES 'nonsense'."
4. **Look.** Click an analog, or move through the list with the arrow keys. Click the selected
   analog again to deselect it, and the view shows the reference alone. In the 3D view the reference
   is gray, the analog's matched core is orange and the rest of the analog is teal. The right panel
   shows:

   | Field | Meaning |
   | --- | --- |
   | RMSD of core | distance between the analog's core atoms and the reference's. The core is placed exactly on the reference, so this is normally 0.000 Å |
   | Core atoms | how many heavy atoms were matched |
   | Seed | the random seed of the embedding |
   | Orientation | only for a symmetric core: which of the equivalent matches was used. It is chosen by how well the analog's substituents overlap the reference's |
   | Stereo | one line per stereocenter of the analog; a line in red needs your attention (below) |
   | Protonation | only with "Flag protonation changes in the core" on: one red line per core atom whose formal charge or number of hydrogens differs from the reference atom it matches, for example `N7: neutral, 0 H; reference N7: neutral, 1 H`. The list shows "check protonation" for that analog |
   | Core search | appears when the search hit its time limit, so the core may be smaller than it could be |

   Below the table, a 2D drawing shows the matched core in blue.
5. **Review.** Press **Accept** (key `A`) or **Reject** (key `R`). If the part of the analog outside
   the core sits badly, for example a side chain clashing with where the protein would be, enter
   another seed and press **Align again**. The core stays where it is. Only the free part moves.
6. **Download.** Choose which structures (accepted only, or all aligned that were not rejected), the
   format, whether to include the reference, whether to get one file per structure and, for mol2,
   whether to keep the reference's atom names (below).

The toolbar above the 3D view has **Overlay all aligned** (every aligned analog as thin lines, to see
the whole series at once), **Hydrogens**, **Color reference by element** (off shows it in gray, so
the analog stands out) and **Reset view**.

### The stereo lines

| Line | What happened | Action |
| --- | --- | --- |
| "as given" | the SMILES specified the center, and it was kept | none |
| "not specified, set to match the reference" | the center is in the core and the reference has it, so the reference's configuration was used | none |
| "not specified (...), chosen by energy; a guess, specify it to be sure" | no reference to follow: the configuration with the lower energy was used | specify it in the SMILES if it matters |
| "as given, but the opposite configuration fits the reference better" | the SMILES says the opposite of the reference. Kept as given | check that this is intended |

### Settings (panel 3)

| Setting | Default | Effect |
| --- | --- | --- |
| Match whole rings only | on | a ring is either fully in the common structure or not at all |
| Ring atoms match ring atoms only | off | a ring atom never matches a chain atom |
| Ignore bond orders | on | atoms match by element and connectivity only, so tautomers and protonation variants still match. Leave it on for a reference read from a protonated structure: its rings may be drawn as alternating single and double bonds instead of aromatic |
| Flag protonation changes in the core | off | the alignment is the same either way. With it on, a core atom whose formal charge, or number of attached atoms counting hydrogens, differs from the reference's is reported (panel 4). A substituent in place of a hydrogen is not reported. Leave it off when the series is meant to change protonation states, as MSLD series often are |
| Core search time limit | 10 s | after this, the best common structure found so far is used |
| Seed | 2342 | the seed of the first alignment. The same input and seed give the same result |

Charges and hydrogen counts never affect the match; "Flag protonation changes in the core" only reports them.

## What the downloads contain

- **sdf**: one record per structure, the reference first if included. The record's title is the
  structure's name. Each record has a `Role` field (`reference` or `analog`), and each analog has
  `AlignRMSD` and `AlignSeed`. Atom names are not stored.
- **mol2**: Sybyl atom and bond types (the same typing as the notebook), **new unique atom names
  and no charges**. Both are deliberate: these files are written for CGenFF, which needs unique
  short atom names and computes its own charges. With **Keep the reference's atom names** checked,
  a reference read from mol2 keeps the names in its file. Analogs, which come from SMILES, have no
  names of their own and get new ones as usual. If a kept name is longer than four characters or
  used twice (ignoring case), the download is refused with the atoms named, because CGenFF would
  renumber them. An sdf reference has no atom names, so the option changes nothing for it.
- **One structure per file**: a zip with one `<name>.sdf` or `<name>.mol2` per structure. If an
  analog has the same name as the reference, the reference's file gets `_reference` added.

The file names are `aligned_accepted.sdf`, `aligned_aligned_mol2.zip` and so on, after your choices.

## Troubleshooting

| Problem | Cause and fix |
| --- | --- |
| The browser cannot connect | The tunnel is not open, or it was opened on the cluster instead of your own computer. Run the printed `ssh -L` command on your computer and keep it open. The port in the address must be the port in the tunnel |
| "Missing or wrong access token." | The address lost its `?token=...` part. Copy the whole address printed by the server. After a restart without `--token`, the token is new |
| "Address already in use" at start | Another server uses the port you gave with `--port`. Leave `--port` out and a free one is picked |
| The 3D view stays empty | The browser has no WebGL (some remote desktops and locked-down browsers turn it off). Use another browser or computer |
| Every analog fails with a charged reference | Fixed (#9): charges no longer affect the match. If you still see it, report it with the reference file |
| The core is much smaller than expected | Check that "Ignore bond orders" is on, and look for "Core search hit its time limit" in the details. Raise the limit for large molecules |
| "An alignment is running; wait for it to finish." | One alignment job at a time. Wait for the list to finish |
| Everything is gone | The server keeps everything in memory. A restart, `--reload` or **Start over** forgets the reference, the analogs and the results. Download before you stop |
| Someone else's structures appear | One server holds one session. Two people on the same server share it and overwrite each other's work. Each person should start their own server |

## Notes

- Nothing is written on the server. Files exist only when you press Download.
- 3Dmol.js 2.4.2 (BSD-3-Clause) is bundled in `static/`, so the page loads nothing from the
  internet.
- The page is built on `msldprep.gui.service.AlignmentSession`, which has no web code and can be
  used from Python directly.
