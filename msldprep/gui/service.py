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
import zipfile
from dataclasses import dataclass, field, replace
from typing import Callable

from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor
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
from ..sybyl_typer import SybylTypingError, annotate_sybyl_types, assign_unique_atom_names, format_mol2_text

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


def _stereocenters(mol: Chem.Mol) -> list[str]:
    """The CIP labels of a 3D molecule's stereocentres, read from its coordinates, like "R" or "S"."""
    copy = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(copy)
    return [label for _, label in Chem.FindMolChiralCenters(copy, includeUnassigned=False, useLegacyImplementation=False)]


def _stereo_notes(mol: Chem.Mol, alignment: CoreAlignment) -> list[dict]:
    """One entry per stereocentre of an aligned analog worth showing, with a sentence
    for the page and `flag` set where the user should look."""
    notes = []
    for center in alignment.stereo:
        name = f"{mol.GetAtomWithIdx(center.atom).GetSymbol()}{center.atom + 1}"
        if center.resolved_by == "equivalent":
            continue
        if center.resolved_by == "reference":
            text, flag = f"{name}: not specified, set to match the reference", False
        elif center.resolved_by == "energy":
            reason = "the reference has no stereocentre there" if center.in_core else "outside the matched core"
            text, flag = f"{name}: not specified ({reason}), chosen by energy; a guess, specify it to be sure", True
        elif center.opposes_reference:
            text, flag = f"{name}: as given, but the opposite configuration fits the reference better", True
        else:
            text, flag = f"{name}: as given", False
        notes.append({"atom": name, "resolved_by": center.resolved_by, "in_core": center.in_core,
                      "opposes_reference": center.opposes_reference, "note": text, "flag": flag})
    return notes


EXPORT_FORMATS = ("sdf", "mol2")


@dataclass
class Export:
    """A file ready to send: its suggested name, media type and bytes."""

    filename: str
    media_type: str
    data: bytes


def _structure_text(mol: Chem.Mol, name: str, props: dict[str, str], fmt: str) -> str:
    """One structure as sdf or mol2 text. `mol` is a copy the caller owns."""
    if fmt == "sdf":
        mol.SetProp("_Name", name)
        for key, value in props.items():
            mol.SetProp(key, value)
        buffer = io.StringIO()
        writer = Chem.SDWriter(buffer)
        try:
            writer.write(mol)
        finally:
            writer.close()
        return buffer.getvalue()
    try:
        annotate_sybyl_types(mol)
        assign_unique_atom_names(mol)
        return format_mol2_text(mol, mol_name=name)
    except SybylTypingError as exc:
        raise SessionError(f"Could not write {name} as mol2: {exc}") from exc


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
    # The GUI matches across tautomers and protonation states by default: bond orders
    # differ between them, and a strict match then drops the atoms around the change.
    mcs: McsOptions = field(default_factory=lambda: McsOptions(ignore_bond_orders=True))
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
            self.reference: Chem.Mol | None = None  # with hydrogens
            self.reference_name: str | None = None
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
            self.reference = Chem.AddHs(mol, addCoords=True)
            self.reference_name = _clean_name(name)
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
        ignore_bond_orders: bool | None = None,
        timeout: int | None = None,
        seed: int | None = None,
    ) -> None:
        """Change settings used by the next alignment; None leaves one as is.
        `timeout` is the core search's time limit in seconds, at least 1."""
        with self._lock:
            self._require_idle()
            mcs = self.params.mcs
            if complete_rings_only is not None:
                mcs = replace(mcs, complete_rings_only=complete_rings_only)
            if ring_matches_ring_only is not None:
                mcs = replace(mcs, ring_matches_ring_only=ring_matches_ring_only)
            if ignore_bond_orders is not None:
                mcs = replace(mcs, ignore_bond_orders=ignore_bond_orders)
            if timeout is not None:
                if timeout < 1:
                    raise SessionError("The core search time limit must be at least 1 second.")
                mcs = replace(mcs, timeout=timeout)
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
                    "charged_atoms": sum(1 for a in self.reference.GetAtoms() if a.GetFormalCharge()),
                    "stereocenters": _stereocenters(self.reference),
                },
                "params": {
                    "complete_rings_only": mcs.complete_rings_only,
                    "ring_matches_ring_only": mcs.ring_matches_ring_only,
                    "ignore_bond_orders": mcs.ignore_bond_orders,
                    "timeout": mcs.timeout,
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
                        "mcs_timed_out": bool(e.alignment and e.alignment.mcs_timed_out),
                        "stereo_flag": bool(e.alignment and any(n["flag"] for n in _stereo_notes(e.mol, e.alignment))),
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
                mcs_timed_out=alignment.mcs_timed_out,
                orientations=alignment.orientations,
                orientation=alignment.orientation,
                stereo=_stereo_notes(entry.mol, alignment),
                svg=_depiction_svg(entry.mol, set(alignment.atom_map)),
                reference_svg=_depiction_svg(self.reference, set(alignment.atom_map.values())),
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

    def export(
        self, fmt: str = "sdf", which: str = "aligned", *, include_reference: bool = True, one_per_file: bool = False
    ) -> Export:
        """The structures chosen by `which`, written as `fmt` ("sdf" or "mol2").

        `which` is "accepted" (only analogs marked accepted) or "aligned" (every
        aligned analog not marked rejected). With `include_reference` the
        reference comes first. One file holds them all, or with `one_per_file`
        a zip holds one file per structure, named after it.

        sdf records carry the analog's RMSD and seed. mol2 is written the way
        the notebook does (`sybyl_typer`): Sybyl atom and bond types, freshly
        generated unique atom names, no charges.
        """
        if fmt not in EXPORT_FORMATS:
            raise SessionError(f"The format must be one of {EXPORT_FORMATS}.")
        if which not in ("aligned", "accepted"):
            raise SessionError('`which` must be "aligned" or "accepted".')
        with self._lock:
            if self.reference is None:
                raise SessionError("No reference loaded.")
            analogs = [
                (e.name, Chem.Mol(e.mol), {"Role": "analog", "AlignRMSD": f"{e.alignment.rmsd:.4f}", "AlignSeed": str(e.alignment.seed)})
                for e in self.analogs.values()
                if e.status == ALIGNED and e.decision != REJECTED and (which == "aligned" or e.decision == ACCEPTED)
            ]
            if not analogs:
                raise SessionError(
                    "Nothing to export: no analog is aligned and not rejected."
                    if which == "aligned"
                    else "Nothing to export: no analog is marked accepted."
                )
            items = []
            if include_reference:
                reference_name = self.reference_name or "reference"
                if one_per_file and reference_name in {name for name, _, _ in analogs}:
                    reference_name += "_reference"  # keep its file apart from an analog of the same name
                items.append((reference_name, Chem.Mol(self.reference), {"Role": "reference"}))
            items += analogs

        texts = [(name, _structure_text(mol, name, props, fmt)) for name, mol, props in items]
        if one_per_file:
            archive = io.BytesIO()
            with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as zf:
                for name, text in texts:
                    zf.writestr(f"{name}.{fmt}", text)
            return Export(f"aligned_{which}_{fmt}.zip", "application/zip", archive.getvalue())
        data = "".join(text for _, text in texts).encode()
        media = "chemical/x-mdl-sdfile" if fmt == "sdf" else "chemical/x-mol2"
        return Export(f"aligned_{which}.{fmt}", media, data)
