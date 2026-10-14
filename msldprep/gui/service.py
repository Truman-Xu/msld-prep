"""State and logic behind the alignment GUI: one in-memory `AlignmentSession`
per server, no web types. The FastAPI layer in `app.py` is a thin shell over
this module, so everything here is testable without a server or a browser.

Scope is the notebook's alignment step only: a reference with 3D coordinates,
analogs given as SMILES, each embedded onto the reference's shared substructure
(`align.align_to_reference`). No CGenFF typing and no shared core across the
whole series.
"""

from __future__ import annotations

import io
import logging
import re
import threading
from dataclasses import dataclass, field, replace
from typing import Callable

from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Draw import rdMolDraw2D

from ..align import (
    DEFAULT_EMBED_SEED,
    CoreAlignment,
    McsOptions,
    align_to_reference,
    mol_from_mol2_text,
    mol_from_sdf_text,
)
from ..errors import AlignmentError, MsldPrepError

logger = logging.getLogger(__name__)

# RDKit raises these (not a closed set) when a molecule cannot be embedded.
_EMBEDDING_EXCEPTIONS = (ValueError, RuntimeError)

PENDING, ALIGNED, FAILED = "pending", "aligned", "failed"
UNDECIDED, ACCEPTED, REJECTED = "undecided", "accepted", "rejected"
DECISIONS = (UNDECIDED, ACCEPTED, REJECTED)


class SessionError(MsldPrepError):
    """Raised for a request the session cannot honor (unknown ligand, busy, bad input)."""


def _net_charge(mol: Chem.Mol) -> int:
    return sum(atom.GetFormalCharge() for atom in mol.GetAtoms())


def neutralize(mol: Chem.Mol) -> Chem.Mol:
    """`mol` with its charges removed by RDKit's `Uncharger`, as the notebook
    does for its reference. Hydrogens are taken off for the neutralization (a
    protonated nitrogen must lose its hydrogen atom, not just its charge) and
    put back with coordinates afterwards, so a 3D molecule stays 3D.
    """
    try:
        plain = Chem.RemoveHs(mol)
        uncharged = rdMolStandardize.Uncharger().uncharge(plain)
        return Chem.AddHs(uncharged, addCoords=True)
    except Exception as exc:  # RDKit's errors here are not a documented, closed set
        raise AlignmentError(f"Could not remove the charges from the reference: {exc!r}") from exc


@dataclass
class AnalogRow:
    """One parsed line of the analog table."""

    name: str
    smiles: str
    error: str | None = None


def _smiles_parses(text: str) -> bool:
    with rdBase.BlockLogs():  # a bad line is reported to the user, not to stderr
        return Chem.MolFromSmiles(text) is not None


def _clean_name(name: str) -> str:
    # Names end up in URLs and file names: keep word characters, dot, plus, minus.
    return re.sub(r"[^\w.+\-]", "_", name) or "analog"


def parse_analog_table(text: str) -> list[AnalogRow]:
    """Parse pasted or uploaded analog text into rows.

    One analog per line, fields separated by commas, tabs or spaces. Either
    `SMILES NAME` or `NAME SMILES` is accepted: the first field is taken as the
    SMILES unless it does not parse and the second does. A single field is a
    SMILES and gets the name `analog<line number>`. Blank lines and lines
    starting with `#` are skipped, as is a first line naming a `smiles`
    column. A line whose SMILES does not parse is kept, with its error, so the
    user sees it. Duplicate names get a numeric suffix.
    """
    rows: list[AnalogRow] = []
    seen: dict[str, int] = {}
    for lineno, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = [f for f in re.split(r"[,\t ]+", line) if f]
        if lineno == 1 and any(f.lower() == "smiles" for f in fields):
            continue
        if len(fields) == 1:
            smiles, name = fields[0], f"analog{lineno}"
        else:
            smiles, name = fields[0], fields[1]
            if not _smiles_parses(smiles) and _smiles_parses(fields[1]):
                smiles, name = fields[1], fields[0]
        name = _clean_name(name)
        count = seen.get(name, 0) + 1
        seen[name] = count
        if count > 1:
            name = f"{name}_{count}"
        error = None if _smiles_parses(smiles) else f"Could not parse SMILES {smiles!r}."
        rows.append(AnalogRow(name=name, smiles=smiles, error=error))
    return rows


@dataclass
class AnalogEntry:
    name: str
    smiles: str
    status: str = PENDING
    decision: str = UNDECIDED
    error: str | None = None
    mol: Chem.Mol | None = None
    alignment: CoreAlignment | None = None
    unparseable: bool = False  # the SMILES itself is bad, whatever the reference


@dataclass
class AlignParams:
    mcs: McsOptions = field(default_factory=McsOptions)
    seed: int = DEFAULT_EMBED_SEED


@dataclass
class JobState:
    running: bool = False
    done: int = 0
    total: int = 0


def _depiction_svg(mol: Chem.Mol, highlight: set[int], size: tuple[int, int] = (340, 260)) -> str:
    """Flat 2D drawing of `mol` without hydrogens, with the atoms (and the
    bonds between them) whose indices in `mol` are in `highlight` marked.
    """
    flat = Chem.Mol(mol)
    for atom in flat.GetAtoms():
        atom.SetIntProp("_orig", atom.GetIdx())
    flat = Chem.RemoveHs(flat)
    flat.RemoveAllConformers()
    rdDepictor.Compute2DCoords(flat)
    atoms = [a.GetIdx() for a in flat.GetAtoms() if a.GetIntProp("_orig") in highlight]
    chosen = set(atoms)
    bonds = [
        b.GetIdx() for b in flat.GetBonds() if b.GetBeginAtomIdx() in chosen and b.GetEndAtomIdx() in chosen
    ]
    drawer = rdMolDraw2D.MolDraw2DSVG(*size)
    drawer.drawOptions().setHighlightColour((0.55, 0.78, 0.98, 1.0))
    drawer.DrawMolecule(
        flat,
        highlightAtoms=atoms,
        highlightBonds=bonds,
    )
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


class AlignmentSession:
    """The reference, the analogs and their alignments. Safe to call from
    several threads: state is guarded by a lock and the slow embedding runs
    outside it.
    """

    def __init__(self) -> None:
        self._lock = threading.RLock()
        self.reset()

    # -- state ---------------------------------------------------------------

    def reset(self) -> None:
        with self._lock:
            self.reference: Chem.Mol | None = None  # as used: neutralized if asked, with hydrogens
            self.reference_loaded: Chem.Mol | None = None  # as read from the file
            self.reference_name: str | None = None
            self.neutralize_reference = False
            self.analogs: dict[str, AnalogEntry] = {}
            self.params = AlignParams()
            self.job = JobState()

    def _require_idle(self) -> None:
        if self.job.running:
            raise SessionError("An alignment is running; wait for it to finish.")

    def _entry(self, name: str) -> AnalogEntry:
        try:
            return self.analogs[name]
        except KeyError:
            raise SessionError(f"No analog named {name!r}.") from None

    # -- inputs --------------------------------------------------------------

    def set_reference(self, name: str, text: str) -> None:
        """Load the reference from mol2 or sdf text (told apart by content).
        Hydrogens are added if the file has none. Clears earlier alignments,
        which no longer refer to this reference.
        """
        if "@<TRIPOS>MOLECULE" in text:
            mol = mol_from_mol2_text(text)
        else:
            mol = mol_from_sdf_text(text)
        with self._lock:
            self._require_idle()
            self.reference_loaded = mol
            self.reference_name = _clean_name(name)
            self._rebuild_reference()

    def set_neutralize_reference(self, neutralize: bool) -> None:
        """Remove (or keep) the charges of the reference, as the notebook does with
        `Uncharger`. Works before a reference is loaded too; it then applies to
        the next one. Earlier alignments no longer fit the changed reference and
        are cleared."""
        with self._lock:
            self._require_idle()
            self.neutralize_reference = bool(neutralize)
            if self.reference_loaded is not None:
                self._rebuild_reference()

    def _rebuild_reference(self) -> None:
        loaded = self.reference_loaded
        self.reference = neutralize(loaded) if self.neutralize_reference else Chem.AddHs(loaded, addCoords=True)
        for entry in self.analogs.values():
            self._clear_alignment(entry)

    def set_analogs(self, text: str) -> list[AnalogRow]:
        rows = parse_analog_table(text)
        if not rows:
            raise SessionError("No analogs found: give one SMILES per line, optionally followed by a name.")
        with self._lock:
            self._require_idle()
            self.analogs = {}
            for row in rows:
                entry = AnalogEntry(name=row.name, smiles=row.smiles)
                if row.error:
                    entry.status, entry.error, entry.unparseable = FAILED, row.error, True
                self.analogs[row.name] = entry
        return rows

    def set_params(
        self,
        *,
        complete_rings_only: bool | None = None,
        ring_matches_ring_only: bool | None = None,
        timeout: int | None = None,
        seed: int | None = None,
    ) -> None:
        """Change settings used by the next alignment; None leaves one as is.
        A `timeout` of 0 returns RDKit's default."""
        with self._lock:
            self._require_idle()
            mcs = self.params.mcs
            if complete_rings_only is not None:
                mcs = replace(mcs, complete_rings_only=complete_rings_only)
            if ring_matches_ring_only is not None:
                mcs = replace(mcs, ring_matches_ring_only=ring_matches_ring_only)
            if timeout is not None:
                if timeout < 0:
                    raise SessionError("The MCS timeout cannot be negative.")
                mcs = replace(mcs, timeout=timeout or None)
            self.params = AlignParams(mcs=mcs, seed=self.params.seed if seed is None else seed)

    # -- alignment -----------------------------------------------------------

    @staticmethod
    def _clear_alignment(entry: AnalogEntry) -> None:
        entry.mol, entry.alignment = None, None
        entry.decision = UNDECIDED
        if not entry.unparseable:
            entry.status, entry.error = PENDING, None

    def _align_entry(self, name: str, reference: Chem.Mol, params: AlignParams, seed: int) -> None:
        """Embed one analog; record the result or the reason it failed."""
        with self._lock:
            smiles = self._entry(name).smiles
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmiles(smiles)
        alignment, error = None, None
        if mol is None:
            error = f"Could not parse SMILES {smiles!r}."
        else:
            mol = Chem.AddHs(mol)
            try:
                alignment = align_to_reference(mol, reference, mcs=params.mcs, seed=seed)
            except _EMBEDDING_EXCEPTIONS as exc:
                error = str(exc) or exc.__class__.__name__
                if "match the core" in error and any(a.GetFormalCharge() for a in reference.GetAtoms()):
                    # #9: the reference's charges are part of the core the analog must match
                    error += " The reference has charged atoms: try removing its charges."
        with self._lock:
            entry = self._entry(name)
            entry.decision = UNDECIDED
            if alignment is None:
                entry.status, entry.error, entry.mol, entry.alignment = FAILED, error, None, None
                logger.info("Alignment failed for %s: %s", name, error)
            else:
                entry.status, entry.error, entry.mol, entry.alignment = ALIGNED, None, mol, alignment

    def begin_alignment(self) -> Callable[[], None]:
        """Check that an alignment of every analog can start, mark the job as
        running and return the function that does the work. The web layer
        calls this in the request and runs the returned function in a thread,
        so a client polling `job` right after the request never sees an idle
        session before the work starts."""
        with self._lock:
            self._require_idle()
            if self.reference is None:
                raise SessionError("Load a reference first.")
            if not self.analogs:
                raise SessionError("Add analogs first.")
            reference, params = self.reference, self.params
            names = list(self.analogs)
            self.job = JobState(running=True, done=0, total=len(names))

        def run() -> None:
            try:
                for name in names:
                    self._align_entry(name, reference, params, params.seed)
                    with self._lock:
                        self.job.done += 1
            finally:
                with self._lock:
                    self.job.running = False

        return run

    def align_all(self) -> None:
        """Align every analog with the current settings, blocking until done."""
        self.begin_alignment()()

    def realign(self, name: str, seed: int) -> None:
        """Embed one analog again, with `seed`; the current settings apply."""
        with self._lock:
            self._require_idle()
            if self.reference is None:
                raise SessionError("Load a reference first.")
            self._entry(name)
            reference, params = self.reference, self.params
        self._align_entry(name, reference, params, seed)

    def set_decision(self, name: str, decision: str) -> None:
        if decision not in DECISIONS:
            raise SessionError(f"Decision must be one of {DECISIONS}.")
        with self._lock:
            entry = self._entry(name)
            if entry.status != ALIGNED and decision != UNDECIDED:
                raise SessionError(f"{name} is not aligned, so it cannot be {decision}.")
            entry.decision = decision

    # -- views ---------------------------------------------------------------

    def summary(self) -> dict:
        with self._lock:
            mcs = self.params.mcs
            return {
                "reference": None
                if self.reference is None
                else {
                    "name": self.reference_name,
                    "atoms": self.reference.GetNumAtoms(),
                    "net_charge": _net_charge(self.reference),
                    "loaded_net_charge": _net_charge(self.reference_loaded),
                    "charged_atoms": sum(1 for a in self.reference.GetAtoms() if a.GetFormalCharge()),
                },
                "neutralize_reference": self.neutralize_reference,
                "params": {
                    "complete_rings_only": mcs.complete_rings_only,
                    "ring_matches_ring_only": mcs.ring_matches_ring_only,
                    "timeout": mcs.timeout or 0,
                    "seed": self.params.seed,
                },
                "job": {"running": self.job.running, "done": self.job.done, "total": self.job.total},
                "analogs": [
                    {
                        "name": e.name,
                        "smiles": e.smiles,
                        "status": e.status,
                        "decision": e.decision,
                        "error": e.error,
                        "rmsd": e.alignment.rmsd if e.alignment else None,
                        "seed": e.alignment.seed if e.alignment else None,
                        "core_atoms": len(e.alignment.atom_map) if e.alignment else None,
                    }
                    for e in self.analogs.values()
                ],
            }

    def reference_view(self) -> dict:
        with self._lock:
            if self.reference is None:
                raise SessionError("No reference loaded.")
            return {
                "name": self.reference_name,
                "molblock": Chem.MolToMolBlock(self.reference),
                "svg": _depiction_svg(self.reference, set()),
                "elements": sorted({a.GetSymbol() for a in self.reference.GetAtoms()}),
            }

    def ligand_view(self, name: str) -> dict:
        """Everything the page needs to show one analog: its 3D block, the
        reference atoms and analog atoms that form the matched core (indices
        into the blocks as sent), and a 2D drawing with the core marked."""
        with self._lock:
            entry = self._entry(name)
            view = {
                "name": entry.name,
                "smiles": entry.smiles,
                "status": entry.status,
                "decision": entry.decision,
                "error": entry.error,
            }
            if entry.status != ALIGNED:
                return view
            alignment = entry.alignment
            view.update(
                molblock=Chem.MolToMolBlock(entry.mol),
                core_atoms=sorted(alignment.atom_map),
                reference_core_atoms=sorted(alignment.atom_map.values()),
                rmsd=alignment.rmsd,
                seed=alignment.seed,
                svg=_depiction_svg(entry.mol, set(alignment.atom_map)),
            )
            return view

    def overlay(self) -> list[dict]:
        """3D blocks of every aligned analog that was not rejected."""
        with self._lock:
            return [
                {"name": e.name, "molblock": Chem.MolToMolBlock(e.mol)}
                for e in self.analogs.values()
                if e.status == ALIGNED and e.decision != REJECTED
            ]

    # -- output --------------------------------------------------------------

    def export_sdf(self, which: str = "aligned") -> str:
        """The reference followed by the analogs, as one sdf text: `which` is
        "accepted" (only those marked accepted) or "aligned" (every aligned
        analog not marked rejected). Each record carries its RMSD and seed."""
        if which not in ("aligned", "accepted"):
            raise SessionError('`which` must be "aligned" or "accepted".')
        with self._lock:
            if self.reference is None:
                raise SessionError("No reference loaded.")
            buffer = io.StringIO()
            writer = Chem.SDWriter(buffer)
            try:
                reference = Chem.Mol(self.reference)
                reference.SetProp("_Name", self.reference_name or "reference")
                reference.SetProp("Role", "reference")
                writer.write(reference)
                for entry in self.analogs.values():
                    if entry.status != ALIGNED or entry.decision == REJECTED:
                        continue
                    if which == "accepted" and entry.decision != ACCEPTED:
                        continue
                    mol = Chem.Mol(entry.mol)
                    mol.SetProp("_Name", entry.name)
                    mol.SetProp("Role", "analog")
                    mol.SetProp("AlignRMSD", f"{entry.alignment.rmsd:.4f}")
                    mol.SetProp("AlignSeed", str(entry.alignment.seed))
                    writer.write(mol)
            finally:
                writer.close()
            return buffer.getvalue()
