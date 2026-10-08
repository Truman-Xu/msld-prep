"""Pure-Python/RDKit reimplementation of OpenBabel's Sybyl atom/bond typer.

Ported directly from OpenBabel's own source, so it reproduces what
`obabel ... -O out.mol2` would write, without shelling out to the `obabel`
binary:

- `OBAtomTyper::AssignHyb`/`AssignTypes` (src/typer.cpp) -- the two-pass
  SMARTS rule engine (INTHYB then EXTTYP, both "last match matched wins",
  driven off data/atomtyp.txt) that assigns each atom an internal type
  name (e.g. "N3", "Car", "O.co2").
- `OBTypeTable::Translate` (src/data.cpp) -- the internal-type -> Sybyl
  (SYB) lookup, driven off data/types.txt, first-matching-row-wins.
- `OBBond::IsAmide`/`IsCarbonyl` and the mol2 writer's bond-label choice
  (src/bond.cpp, src/formats/mol2format.cpp) -- the "1"/"2"/"3"/"ar"/"am"
  bond-type string written to a mol2 BOND record.

The rule tables themselves (`sybyl_typing_data.py`) are OpenBabel's own
data/atomtyp.txt and data/types.txt (confirmed byte-identical to
`openbabel/data/{atomtyp,types}.txt` upstream), pre-parsed once by
`scripts/generate_sybyl_typing_data.py` into a plain Python module -- no
text parsing happens at runtime, only SMARTS compilation (RDKit has no
serializable query-mol format) and rule application.

Known, deliberate divergence: OpenBabel computes its own aromaticity
perception from scratch (typer.cpp's `OBAromaticTyper`, a bespoke
SSSR/Hueckel-electron-counting walk) before running the INTHYB/EXTTYP
passes. Reimplementing that model too is out of scope here; RDKit's own
aromaticity perception (`Atom.GetIsAromatic()`/`Bond.GetIsAromatic()`) is
used instead. The two models agree on ordinary drug-like aromatic rings
(validated against real `obabel` output); exotic ring systems could in
principle type differently.

Also out of scope: OpenBabel's bond-order-PERCEPTION-from-3D-geometry
(src/bondtyper.cpp) -- that solves a different problem (guessing bond
orders for formats that don't carry them, e.g. plain XYZ/PDB). Every Mol
this module handles already carries explicit, correct bond orders from
RDKit's own SMILES parsing, so it isn't needed.
"""

from __future__ import annotations

import logging

from rdkit import Chem

from .errors import MsldPrepError
from .naming import allocate_unique_atom_name, new_counters
from .sybyl_typing_data import ATN_TO_INT, EXTTYP_RULES, INT_TO_SYB, INTHYB_RULES

logger = logging.getLogger(__name__)


class SybylTypingError(MsldPrepError):
    """Raised when Sybyl atom/bond typing or mol2 writing cannot be completed."""


# The one EXTTYP rule whose `^N` hybridization query sits on an atom
# INSIDE a recursive SMARTS sub-pattern (the central carbon), not the
# outer/typed atom (the oxygen). Recursive-SMARTS atoms aren't exposed in
# GetSubstructMatches' match tuple, so this can't go through the generic
# strip-^N-and-postfilter-match[0] path used for every other rule; it's
# reproduced directly in `_match_carboxylate_oxygen` instead.
_CARBOXYLATE_O_SMARTS = "[$([#8D1]~[#6^2]~[#8D1])]"


_compiled_inthyb_rules: list[tuple[Chem.Mol, int]] | None = None
_compiled_exttyp_rules: list[tuple[str, Chem.Mol, int | None, str]] | None = None


def _load_rules() -> tuple[
    list[tuple[Chem.Mol, int]], list[tuple[str, Chem.Mol, int | None, str]], dict[str, str], dict[int, str]
]:
    """Compile `sybyl_typing_data`'s SMARTS strings into RDKit queries, once
    per process (RDKit query mols aren't picklable/cacheable across
    processes, so this compilation step can't itself be pushed into the
    generated data module -- only the SMARTS text can be).
    """
    global _compiled_inthyb_rules, _compiled_exttyp_rules
    if _compiled_inthyb_rules is None:
        inthyb: list[tuple[Chem.Mol, int]] = []
        for smarts, value in INTHYB_RULES:
            query = Chem.MolFromSmarts(smarts)
            if query is None:
                logger.warning("Could not parse INTHYB SMARTS %r; skipping", smarts)
                continue
            inthyb.append((query, value))

        exttyp: list[tuple[str, Chem.Mol, int | None, str]] = []
        for smarts, cleaned, required_hyb, value in EXTTYP_RULES:
            query = Chem.MolFromSmarts(cleaned)
            if query is None:
                logger.warning("Could not parse EXTTYP SMARTS %r; skipping", smarts)
                continue
            exttyp.append((smarts, query, required_hyb, value))

        _compiled_inthyb_rules, _compiled_exttyp_rules = inthyb, exttyp
        logger.debug(
            "Compiled %d INTHYB rules, %d EXTTYP rules, %d INT->SYB entries",
            len(inthyb), len(exttyp), len(INT_TO_SYB),
        )
    assert _compiled_exttyp_rules is not None
    return _compiled_inthyb_rules, _compiled_exttyp_rules, INT_TO_SYB, ATN_TO_INT


def _assign_hybridization(mol: Chem.Mol, inthyb_rules: list[tuple[Chem.Mol, int]]) -> dict[int, int]:
    """Port of `OBAtomTyper::AssignHyb` (typer.cpp:180-227): apply every
    INTHYB rule in file order (last match wins), then fall back to a
    degree-based guess for any atom no rule touched.
    """
    hyb: dict[int, int] = {}
    for query, value in inthyb_rules:
        for match in mol.GetSubstructMatches(query):
            hyb[match[0]] = value

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in hyb:
            continue
        degree = atom.GetDegree()  # explicit Hs assumed, matching OB's GetExplicitDegree()
        if degree <= 2:
            hyb[idx] = 1
        elif degree == 3:
            hyb[idx] = 2
        elif degree == 4:
            hyb[idx] = 3
        else:
            hyb[idx] = degree
    return hyb


def _match_carboxylate_oxygen(mol: Chem.Mol, hyb: dict[int, int]) -> set[int]:
    """EXTTYP `[$([#8D1]~[#6^2]~[#8D1])]` -> O.co2: a degree-1 oxygen bonded
    to a hyb=2 carbon that is itself bonded to another degree-1 oxygen.
    """
    matches: set[int] = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.GetDegree() != 1:
            continue
        carbon = atom.GetNeighbors()[0]
        if carbon.GetAtomicNum() != 6 or hyb.get(carbon.GetIdx()) != 2:
            continue
        other_terminal_oxygens = [
            nbr
            for nbr in carbon.GetNeighbors()
            if nbr.GetIdx() != atom.GetIdx() and nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 1
        ]
        if other_terminal_oxygens:
            matches.add(atom.GetIdx())
    return matches


def _assign_ext_types(
    mol: Chem.Mol,
    exttyp_rules: list[tuple[str, Chem.Mol, int | None, str]],
    hyb: dict[int, int],
) -> dict[int, str]:
    """Port of `OBAtomTyper::AssignTypes` (typer.cpp:135-178): apply every
    EXTTYP rule in file order (last match wins, matched atom = match[0]),
    then the guanidinium C2->C+ special case.
    """
    ext_type: dict[int, str] = {}
    for smarts_text, query, required_hyb, value in exttyp_rules:
        if smarts_text == _CARBOXYLATE_O_SMARTS:
            for idx in _match_carboxylate_oxygen(mol, hyb):
                ext_type[idx] = value
            continue
        for match in mol.GetSubstructMatches(query):
            target = match[0]
            if required_hyb is not None and hyb.get(target) != required_hyb:
                continue
            ext_type[target] = value

    # Guanidinium special case (typer.cpp:157-177, "Fixes PR#1800964"): a
    # "C2"-typed carbon with exactly 3 neighbors typed Npl/N2/Ng+ is the
    # symmetric central carbon of a guanidinium group, retyped cationic.
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if ext_type.get(idx) != "C2":
            continue
        guanidine_n = sum(
            1 for nbr in atom.GetNeighbors() if ext_type.get(nbr.GetIdx()) in ("Npl", "N2", "Ng+")
        )
        if guanidine_n == 3:
            ext_type[idx] = "C+"
    return ext_type


def assign_internal_types(mol: Chem.Mol) -> dict[int, str]:
    """OpenBabel's internal atom type name per atom index (e.g. "N3",
    "Car", "O.co2") -- the value `atom->GetType()` would return, before
    translation to Sybyl. Exposed separately from `assign_sybyl_types`
    because the mol2 bond-type rule needs to test for the untranslated
    "O.co2" internal type specifically (see `sybyl_bond_type`).
    """
    inthyb_rules, exttyp_rules, _, atn_to_int = _load_rules()
    hyb = _assign_hybridization(mol, inthyb_rules)
    ext_type = _assign_ext_types(mol, exttyp_rules, hyb)

    internal: dict[int, str] = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        int_type = ext_type.get(idx)
        if int_type is None:
            # No EXTTYP rule matched (typical for H, halogens, metals):
            # port of OBAtom::GetType()'s ATN->INT fallback (atom.cpp:475-488).
            int_type = atn_to_int.get(atom.GetAtomicNum())
            if int_type is None:
                msg = f"No Sybyl typing fallback for atomic number {atom.GetAtomicNum()} (atom idx {idx})"
                logger.error(msg)
                raise SybylTypingError(msg)
        internal[idx] = int_type
    return internal


def assign_sybyl_types(mol: Chem.Mol) -> dict[int, str]:
    """Sybyl (SYB) atom type per atom index, reproducing what OpenBabel's
    own `obabel ... -O out.mol2` would write to the ATOM section's type
    column (see module docstring for the one known divergence: aromaticity
    perception uses RDKit's model, not OpenBabel's).
    """
    _, _, int_to_syb, _ = _load_rules()
    internal = assign_internal_types(mol)
    sybyl: dict[int, str] = {}
    for idx, int_type in internal.items():
        syb = int_to_syb.get(int_type)
        if syb is None:
            msg = f"No SYB column entry for internal type {int_type!r} (atom idx {idx})"
            logger.error(msg)
            raise SybylTypingError(msg)
        sybyl[idx] = syb
    return sybyl


def _is_carbonyl_bond(bond: Chem.Bond) -> bool:
    """Port of `OBBond::IsCarbonyl` (bond.cpp:454-464)."""
    if bond.GetBondType() != Chem.BondType.DOUBLE:
        return False
    return {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()} == {6, 8}


def _is_amide_bond(bond: Chem.Bond) -> bool:
    """Port of `OBBond::IsAmide` (bond.cpp:192-222)."""
    begin, end = bond.GetBeginAtom(), bond.GetEndAtom()
    if begin.GetAtomicNum() == 6 and end.GetAtomicNum() == 7:
        carbon, nitrogen = begin, end
    elif begin.GetAtomicNum() == 7 and end.GetAtomicNum() == 6:
        carbon, nitrogen = end, begin
    else:
        return False
    if bond.GetBondType() != Chem.BondType.SINGLE:
        return False
    if nitrogen.GetDegree() != 3:  # must be a degree-3 nitrogen
        return False
    return any(_is_carbonyl_bond(b) for b in carbon.GetBonds())


def sybyl_bond_type(bond: Chem.Bond, internal_types: dict[int, str]) -> str:
    """Sybyl mol2 BOND-section type label ("1"/"2"/"3"/"ar"/"am"), a port
    of mol2format.cpp:785-794. `internal_types` must be the untranslated
    internal type names from `assign_internal_types` -- the "O.co2" check
    there runs on `atom->GetType()`, i.e. before Sybyl translation (which,
    for O.co2 specifically, happens to be an identity mapping anyway, but
    the source checks the internal value so this does too).
    """
    i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    if bond.GetIsAromatic() or internal_types.get(i) == "O.co2" or internal_types.get(j) == "O.co2":
        return "ar"
    if _is_amide_bond(bond):
        return "am"
    return str(int(round(bond.GetBondTypeAsDouble())))


ATOM_TYPE_PROP = "SybylType"
BOND_TYPE_PROP = "SybylBondType"


def annotate_sybyl_types(mol: Chem.Mol) -> None:
    """Compute Sybyl atom types and mol2 bond-type labels for every atom
    and bond in `mol`, the same way OpenBabel's own `obabel ... -O
    out.mol2` would type them (see module docstring for the one known
    divergence: aromaticity perception uses RDKit's model, not
    OpenBabel's). Mutates `mol` in place: sets a string prop `ATOM_TYPE_PROP`
    ("SybylType") on every atom and `BOND_TYPE_PROP` ("SybylBondType") on
    every bond -- `write_mol2_file` reads these back rather than
    recomputing them, so typing and writing can be called (and tested)
    independently.
    """
    internal_types = assign_internal_types(mol)
    _, _, int_to_syb, _ = _load_rules()
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        int_type = internal_types[idx]
        syb = int_to_syb.get(int_type)
        if syb is None:
            msg = f"No SYB column entry for internal type {int_type!r} (atom idx {idx})"
            logger.error(msg)
            raise SybylTypingError(msg)
        atom.SetProp(ATOM_TYPE_PROP, syb)
    for bond in mol.GetBonds():
        bond.SetProp(BOND_TYPE_PROP, sybyl_bond_type(bond, internal_types))


def assign_unique_atom_names(mol: Chem.Mol) -> dict[int, str]:
    """Allocate a fresh globally-unique <=4-char name per atom (via
    `naming.allocate_unique_atom_name`) and set it as that atom's
    `"_TriposAtomName"` prop -- the same prop RDKit's own mol2 reader sets,
    so downstream code that expects to find names there
    (`rtf_records.assign_atom_types`, etc.) works identically
    whether `mol` was loaded from a real mol2 file or built and named
    in-memory. CGenFF (the downstream consumer) preserves whatever unique
    names it's given but silently renumbers plain/duplicate element-symbol
    names on its own, which is what this exists to avoid.

    Returns the RDKit atom index -> assigned name mapping.
    """
    counters = new_counters()
    names: dict[int, str] = {}
    for atom in mol.GetAtoms():
        name = allocate_unique_atom_name(atom.GetSymbol(), counters)
        atom.SetProp("_TriposAtomName", name)
        names[atom.GetIdx()] = name
    return names


def format_mol2_text(mol: Chem.Mol, *, mol_name: str = "LIG") -> str:
    """Format `mol` (already typed via `annotate_sybyl_types` and named via
    `assign_unique_atom_names`, with one conformer of 3D coords and
    explicit Hs) as Tripos mol2 text -- no file I/O, so callers that need
    the mol2 purely in-memory (e.g. to pipe into a subprocess's stdin)
    don't need a throwaway file just to get it. Charges are written as 0.0
    with charge_type NO_CHARGES -- CGenFF computes its own charges from a
    bond-increment model and does not read this column.
    """
    if mol.GetNumConformers() == 0:
        msg = f"Cannot format a mol2 block for `{mol_name}`: molecule has no conformer."
        logger.error(msg)
        raise SybylTypingError(msg)
    if mol.GetNumAtoms() > 0:
        first = mol.GetAtomWithIdx(0)
        if not first.HasProp(ATOM_TYPE_PROP):
            msg = f"`{mol_name}` has not been typed yet -- call annotate_sybyl_types(mol) first."
            logger.error(msg)
            raise SybylTypingError(msg)
        if not first.HasProp("_TriposAtomName"):
            msg = f"`{mol_name}` has not been named yet -- call assign_unique_atom_names(mol) first."
            logger.error(msg)
            raise SybylTypingError(msg)

    conf = mol.GetConformer()
    lines = [
        "@<TRIPOS>MOLECULE",
        mol_name,
        f"{mol.GetNumAtoms():>5}{mol.GetNumBonds():>6}     0     0     0",
        "SMALL",
        "NO_CHARGES",
        "",
        "@<TRIPOS>ATOM",
    ]
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        lines.append(
            f"{idx + 1:>7} {atom.GetProp('_TriposAtomName'):<8} {pos.x:>9.4f} {pos.y:>9.4f} {pos.z:>9.4f} "
            f"{atom.GetProp(ATOM_TYPE_PROP):<6} {1:>4} {mol_name:<8} {0.0:>9.4f}"
        )
    lines.append("@<TRIPOS>BOND")
    for bond_id, bond in enumerate(mol.GetBonds(), start=1):
        label = bond.GetProp(BOND_TYPE_PROP)
        lines.append(f"{bond_id:>6} {bond.GetBeginAtomIdx() + 1:>5} {bond.GetEndAtomIdx() + 1:>5}   {label}")

    return "\n".join(lines) + "\n"


def write_mol2_file(mol: Chem.Mol, path: str, *, mol_name: str = "MOL") -> dict[int, str]:
    """Write `mol` (already typed by `annotate_sybyl_types`, with one
    conformer of 3D coords and explicit Hs) as a Tripos mol2 file, without
    shelling out to `obabel`. Names every atom via `assign_unique_atom_names`
    first, unless that's already been done (idempotent).

    Returns the RDKit atom index -> assigned mol2 name mapping.
    """
    if mol.GetNumAtoms() > 0 and not mol.GetAtomWithIdx(0).HasProp("_TriposAtomName"):
        names = assign_unique_atom_names(mol)
    else:
        names = {atom.GetIdx(): atom.GetProp("_TriposAtomName") for atom in mol.GetAtoms()}

    with open(path, "w") as fp:
        fp.write(format_mol2_text(mol, mol_name=mol_name))

    logger.info("Wrote mol2 file %s (%d atoms, %d bonds)", path, mol.GetNumAtoms(), mol.GetNumBonds())
    return names
