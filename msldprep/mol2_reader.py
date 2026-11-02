"""Read mol2 into RDKit, from a file or from text, with OpenBabel as the
fallback for what RDKit's own mol2 reader refuses.

RDKit refuses two patterns that CGenFF-ready mol2 files do contain: a nitro
group typed `N.pl3` with an `O.co2` (rejected while the atoms are read), and
a ring NH typed `N.ar` (the ring then cannot be kekulized). Legacy reads
mol2 with its own parser, so these ligands work there.

The fallback pipes the mol2 text through `obabel` to sdf and reads that, so
nothing is written to disk. OpenBabel
writes a nitro group as `N(=O)[O-]` with no charge on the nitrogen, so a
4-valent neutral N bonded to an O- gets +1 before sanitizing. The result is
used only if it has the mol2's atoms, in the same order, elements and
positions; the atom names are then copied from the mol2. Atom types and
charges come from the rtf afterwards, as for any other ligand.

`obabel` is found through the `MSLDPREP_OBABEL_PATH` environment variable,
else on `PATH`. Without it, mol2 RDKit refuses is refused with its cause.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess

import numpy as np
from rdkit import Chem

from .errors import AlignmentError

logger = logging.getLogger(__name__)

OBABEL_PATH_ENV = "MSLDPREP_OBABEL_PATH"
_COORD_TOLERANCE = 1e-3  # angstrom; mol2 and sdf both carry 4 decimals


def find_obabel() -> str | None:
    """The `obabel` program from `MSLDPREP_OBABEL_PATH`, else `PATH`, or None."""
    path = os.environ.get(OBABEL_PATH_ENV) or shutil.which("obabel")
    return path if path and os.path.isfile(path) else None


def _obabel_env(obabel: str) -> dict[str, str]:
    # A cluster install may need its own lib directories to load (libcoordgen, for one).
    env = dict(os.environ)
    prefix = os.path.dirname(os.path.dirname(os.path.abspath(obabel)))
    lib_dirs = [os.path.join(prefix, "lib"), os.path.join(prefix, "lib64")]
    existing = env.get("LD_LIBRARY_PATH", "")
    env["LD_LIBRARY_PATH"] = ":".join([*lib_dirs, existing]) if existing else ":".join(lib_dirs)
    return env


def _mol2_atoms(text: str) -> tuple[list[str], list[str], np.ndarray]:
    """Names, elements and coordinates of the ATOM section, in file order."""
    names, elements, coords = [], [], []
    in_atoms = False
    for line in text.splitlines():
        if line.startswith("@<TRIPOS>"):
            in_atoms = line.startswith("@<TRIPOS>ATOM")
            continue
        fields = line.split()
        if in_atoms and len(fields) >= 6:
            names.append(fields[1])
            coords.append([float(x) for x in fields[2:5]])
            elements.append(fields[5].split(".")[0])
    return names, elements, np.array(coords)


def _rdkit_failure(text: str) -> str:
    """Why RDKit refuses `text`, in words a user can act on."""
    mol = Chem.MolFromMol2Block(text, removeHs=False, sanitize=False)
    if mol is None:
        return (
            "RDKit rejects its atom types while reading them, for example an O.co2 "
            "that is not part of a carboxylate, as in a nitro group"
        )
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:  # RDKit's sanitize errors are not a documented, closed set
        return f"RDKit cannot sanitize it: {str(exc).splitlines()[0]}"
    return "RDKit could not read it"


def _charge_nitro_nitrogens(mol: Chem.Mol) -> int:
    """Give +1 to each 4-valent neutral N bonded to an O-; return how many."""
    charged = 0
    for atom in mol.GetAtoms():
        valence = sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
        if (
            atom.GetSymbol() == "N"
            and atom.GetFormalCharge() == 0
            and valence == 4
            and any(n.GetSymbol() == "O" and n.GetFormalCharge() == -1 for n in atom.GetNeighbors())
        ):
            atom.SetFormalCharge(1)
            charged += 1
    return charged


def _read_through_obabel(text: str, obabel: str) -> Chem.Mol:
    try:
        result = subprocess.run(
            [obabel, "-imol2", "-osdf"],
            input=text, capture_output=True, text=True, timeout=60, env=_obabel_env(obabel),
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise AlignmentError(f"obabel could not convert it: {exc!r}") from exc
    if result.returncode != 0 or not result.stdout.strip():
        raise AlignmentError(f"obabel could not convert it: {result.stderr.strip()}")

    mol = Chem.MolFromMolBlock(result.stdout, sanitize=False, removeHs=False)
    if mol is None:
        raise AlignmentError("RDKit cannot read obabel's sdf either")
    _charge_nitro_nitrogens(mol)
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:  # RDKit's sanitize errors are not a documented, closed set
        raise AlignmentError(f"obabel's structure does not sanitize: {str(exc).splitlines()[0]}") from exc

    names, elements, coords = _mol2_atoms(text)
    if mol.GetNumAtoms() != len(names):
        raise AlignmentError(f"obabel's structure has {mol.GetNumAtoms()} atoms, the mol2 {len(names)}")
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    if [s.upper() for s in symbols] != [e.upper() for e in elements]:
        raise AlignmentError("obabel's structure has its atoms in a different order")
    if np.abs(mol.GetConformer().GetPositions() - coords).max() > _COORD_TOLERANCE:
        raise AlignmentError("obabel's structure has moved atoms")
    for atom, name in zip(mol.GetAtoms(), names):
        atom.SetProp("_TriposAtomName", name)
    return mol


def read_mol2_text(text: str, *, source: str = "mol2 text") -> Chem.Mol:
    """Read one molecule from mol2 `text`, with explicit hydrogens, sanitized,
    atoms named by `_TriposAtomName`.

    Falls back to OpenBabel (see the module docstring) when RDKit refuses it,
    with a warning. Raises AlignmentError naming the cause when neither can
    read it; `source` names the input in messages.
    """
    try:
        mol = Chem.MolFromMol2Block(text, removeHs=False, sanitize=True)
    except Exception:  # RDKit's parse errors are not a documented, closed set
        mol = None
    if mol is not None:
        return mol

    if "@<TRIPOS>ATOM" not in text:
        msg = f"Failed to load {source}: it has no @<TRIPOS>ATOM section."
        logger.error(msg)
        raise AlignmentError(msg)
    cause = _rdkit_failure(text)
    obabel = find_obabel()
    if obabel is None:
        msg = (
            f"Failed to load {source}: {cause}. OpenBabel can read such files: "
            f"set {OBABEL_PATH_ENV} to the obabel program, or put obabel on PATH."
        )
        logger.error(msg)
        raise AlignmentError(msg)
    try:
        mol = _read_through_obabel(text, obabel)
    except AlignmentError as exc:
        msg = f"Failed to load {source}: {cause}, and the OpenBabel fallback failed: {exc}"
        logger.error(msg)
        raise AlignmentError(msg) from exc
    logger.warning("Read %s through OpenBabel, because %s.", source, cause)
    return mol


def read_mol2_file(path: str) -> Chem.Mol:
    """`read_mol2_text` on the contents of the file at `path`."""
    with open(path) as fp:
        return read_mol2_text(fp.read(), source=f"mol2 file `{path}`")
