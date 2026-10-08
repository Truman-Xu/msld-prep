"""Shared data structures that flow between msldprep pipeline stages.

Plain dataclasses only -- no pandas. RDKit `Mol`/`Atom` objects are the
source of truth for connectivity/geometry/CGenFF typing; these dataclasses
capture the additional structure (site/fragment/anchor assignment, charges)
that RDKit has no native concept of.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

from rdkit import Chem

from .errors import AlignmentError
from .prm_records import (
    LigandParameters,
    format_angle_line,
    format_bond_line,
    format_dihedral_lines,
    format_improper_line,
    format_nonbonded_line,
    parse_prm_file,
)
from .rtf_records import LonePairRecord, MassRecord, assign_atom_types, parse_rtf_file

logger = logging.getLogger(__name__)


def _format_lonepair_line(lp) -> str:
    parts = ["LONEPAIR", lp.lp_type, lp.lp_atom, *lp.host_atoms, "DIST", f"{lp.distance:.4f}"]
    if lp.angle:
        parts += ["ANGL", f"{lp.angle:.4f}"]
    if lp.dihedral:
        parts += ["DIHE", f"{lp.dihedral:.4f}"]
    parts += ["SCAL", f"{lp.scale}"]
    return " ".join(parts) + "\n"


@dataclass
class LigandRecord:
    """One ligand: identity, its fully CGenFF-typed/isotope-labeled RDKit
    Mol, and the CHARMM topology/parameter data RDKit has no native concept
    of -- net charge, IMPR/LONEPAIR declarations, and bonded/nonbonded
    force-field terms (`prm`) -- as direct fields, not a nested "aux"
    object.

    Prefer `LigandRecord.from_files` (mol2+rtf+prm on disk) or
    `analog_prep.mol_to_ligand_record` (an in-memory RDKit Mol, typed via a
    streamed CGenFF call, nothing touches disk) over constructing this
    directly.
    """

    name: str
    mol: Chem.Mol
    net_charge: float
    prm: LigandParameters = field(default_factory=LigandParameters)
    format_version: tuple[int, int] = (36, 1)
    atom_types: dict[str, str] = field(default_factory=dict)
    atom_charges: dict[str, float] = field(default_factory=dict)
    bonds: tuple[tuple[str, str], ...] = field(default_factory=tuple)
    impropers: tuple[tuple[str, str, str, str], ...] = field(default_factory=tuple)
    lone_pairs: tuple[LonePairRecord, ...] = field(default_factory=tuple)
    mass_records: dict[str, MassRecord] = field(default_factory=dict)
    is_reference: bool = False
    mol2_path: str | None = None
    rtf_path: str | None = None
    prm_path: str | None = None

    @classmethod
    def from_files(
        cls, mol2_path: str, rtf_path: str, prm_path: str, name: str | None = None, *, is_reference: bool = False
    ) -> "LigandRecord":
        """Load a ligand from mol2 + rtf + prm, typing every atom with its
        CGenFF type/charge and attaching the parsed rtf/prm records.
        """
        mol = Chem.MolFromMol2File(mol2_path, removeHs=False, sanitize=True)
        if mol is None:
            msg = f"Failed to load mol2 file `{mol2_path}` into RDKit."
            logger.error(msg)
            raise AlignmentError(msg)

        return cls.from_mol(
            mol, rtf_path, prm_path, name=name or mol2_path, is_reference=is_reference, mol2_path=mol2_path
        )

    @classmethod
    def from_mol(
        cls,
        mol: Chem.Mol,
        rtf_path: str,
        prm_path: str,
        name: str,
        *,
        is_reference: bool = False,
        mol2_path: str | None = None,
    ) -> "LigandRecord":
        """Type an already-loaded `mol` against rtf + prm, attaching the
        parsed rtf/prm records. `mol` must carry a `_TriposAtomName` prop on
        every atom already (e.g. from RDKit's own mol2 reader, or restored
        by `align.load_ligands_sdf`) -- shared tail of `from_files` (mol2
        source) and `align.load_ligand_set_from_sdf` (sdf source).
        """
        rtf_fields = parse_rtf_file(rtf_path)
        assign_atom_types(mol, rtf_fields["atom_types"], rtf_fields["atom_charges"], source_label=rtf_path)
        prm = parse_prm_file(prm_path)

        return cls(
            name=name,
            mol=mol,
            prm=prm,
            is_reference=is_reference,
            mol2_path=mol2_path,
            rtf_path=rtf_path,
            prm_path=prm_path,
            **rtf_fields,
        )

    def to_rtf_text(self) -> str:
        """Format a standalone single-residue CHARMM .rtf for this ligand,
        matching the structure CGenFF's own `-x` split-output mode produces
        (RESI/[MASS]/ATOM/BOND/IMPR/LONEPAIR/END). A functional round-trip,
        not a byte-identical one: per-atom penalty score comments aren't
        preserved (discarded during parsing, and not meaningful once
        re-derived data no longer matches the original CGenFF run).
        """
        lines = [
            "* ligand rtf file generated with msldprep\n",
            f"* (residue {self.name})\n* \n",
            f"  {self.format_version[0]} {self.format_version[1]}\n\n",
        ]
        for mass in self.mass_records.values():
            lines.append(f"MASS -1 {mass.atom_type} {mass.mass} {mass.element}\n")
        if self.mass_records:
            lines.append("\n")

        lines.append(f"RESI {self.name}    {self.net_charge:.3f}\n")
        lines.append("GROUP \n")
        for atom_name, atom_type in self.atom_types.items():
            lines.append(f"ATOM {atom_name:<4} {atom_type:<6} {self.atom_charges[atom_name]:10.6f}\n")

        lines.append("\n")
        for a, b in self.bonds:
            lines.append(f"BOND {a:<4} {b:<4}\n")
        for impr in self.impropers:
            lines.append("IMPR " + " ".join(f"{a:<4}" for a in impr) + "\n")

        if self.lone_pairs:
            lines.append("\n")
            for lp in self.lone_pairs:
                lines.append(_format_lonepair_line(lp))

        lines.append("\nEND")
        return "".join(lines)

    def write_rtf(self, path: str) -> None:
        with open(path, "w") as fp:
            fp.write(self.to_rtf_text())
        logger.info("Wrote %s", path)

    def to_prm_text(self) -> str:
        """Format a standalone CHARMM .prm for this ligand from `self.prm`,
        matching the structure CGenFF's own `-x` split-output mode produces
        (BONDS/ANGLES/DIHEDRALS/IMPROPERS/[NONBONDED]/[NBFIX]/END). A
        functional round-trip, not a byte-identical one: per-term penalty
        score / analogy-source comments aren't preserved.
        """
        p = self.prm
        lines = ["* ligand prm file generated with msldprep\n*\n\n"]

        lines.append("BONDS\n")
        for key, bp in p.bonds.items():
            lines.append(format_bond_line(key, bp) + "\n")

        lines.append("\nANGLES\n")
        for key, ap in p.angles.items():
            lines.append(format_angle_line(key, ap, p.urey_bradley.get(key)) + "\n")

        lines.append("\nDIHEDRALS\n")
        for key, terms in p.dihedrals.items():
            for line in format_dihedral_lines(key, terms):
                lines.append(line + "\n")

        lines.append("\nIMPROPERS\n")
        for key, ip in p.impropers.items():
            lines.append(format_improper_line(key, ip) + "\n")

        if p.nonbonded:
            lines.append("\nNONBONDED\n")
            for atom_type, nb in p.nonbonded.items():
                lines.append(format_nonbonded_line(atom_type, nb, p.nonbonded14.get(atom_type)) + "\n")

        if p.nbfix:
            lines.append("\nNBFIX\n")
            for (a, b), nf in p.nbfix.items():
                lines.append(f"{a} {b} {nf.emin:.4f} {nf.rmin:.4f}\n")

        lines.append("\nEND")
        return "".join(lines)

    def write_prm(self, path: str) -> None:
        with open(path, "w") as fp:
            fp.write(self.to_prm_text())
        logger.info("Wrote %s", path)


@dataclass
class Fragment:
    """One ligand's substituent (R-group) atoms at a given site."""

    ligand_name: str
    site_index: int
    atom_indices: list[int] = field(default_factory=list)
    atom_names: list[str] = field(default_factory=list)
    lp_atom_names: list[str] = field(default_factory=list)
    anchor_atom_names: list[str] = field(default_factory=list)


@dataclass
class Site:
    """A single alchemical site: one Fragment per ligand that has one."""

    index: int
    fragments: dict[str, Fragment] = field(default_factory=dict)


@dataclass
class CoreDecomposition:
    """Full core/anchor/site/fragment decomposition across a ligand set."""

    reference_name: str
    core_query_mol: Chem.Mol
    core_atom_order: tuple[str, ...]
    core_atom_names: dict[str, tuple[str, ...]]
    # host core atom name (in core_atom_order) -> {ligand_name: this ligand's LP atom name}
    core_lp_atom_names: dict[str, dict[str, str]] = field(default_factory=dict)
    sites: list[Site] = field(default_factory=list)
    # names of ligands whose core match was an ambiguous tie (identical anchor
    # pattern across >1 candidate match), resolved by picking the first one
    tied_core_matches: list[str] = field(default_factory=list)
    # which core group this decomposition belongs to, when a ligand set was
    # split by grouping.group_ligands_by_similarity (see io_mcs_text's "CORE
    # ID" sections); 1 for the single-group case, which is also the default
    # for every decomposition built without grouping at all.
    core_id: int = 1


@dataclass
class CRNResult:
    """Charge-renormalized charges, plus the bookkeeping needed to write PRES charge totals."""

    core_charges: dict[str, float]
    core_charge_stats: dict[str, tuple[float, float]]
    fragment_charges: list[dict[str, dict[str, float]]]
    site_average: list[float]
    site_integer_charge: list[int]
    charge_changed: dict[str, bool]
