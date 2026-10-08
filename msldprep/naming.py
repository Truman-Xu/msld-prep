"""Shared unique-atom-name allocation, used both by the CHARMM writers
(translating raw atom names into <=4-char CHARMM-safe names for the shared
core + fragments, `writers/records.py`) and by `analog_prep.py` (assigning
unique Tripos atom names to a freshly-generated analog before it's handed
to CGenFF, since CGenFF preserves whatever unique names it's given but
otherwise renumbers plain, non-unique element-symbol names on its own).

Each element symbol (and the literal "LP") draws from its own counter, so
a system's naming capacity scales with how many distinct symbols it uses,
not with its total atom count. `compute_atom_name_translation` shares one
`counters` store across an entire multi-ligand system's core plus every
ligand's per-site fragment atoms -- a large congeneric series can easily
carry more than 999 atoms in total even with far fewer than 999 of any
single element, so a single counter shared across every symbol (an earlier
version of this module, ported from msld_crn.py:570-596/689-711's `pnum`
scheme) ran out well before any individual element actually needed to.
"""

from __future__ import annotations

import logging

from .errors import MsldPrepError

logger = logging.getLogger(__name__)


class NamingError(MsldPrepError):
    """Raised when a unique <=4-char atom name can no longer be allocated."""


def new_counters() -> dict[str, list[int]]:
    """A fresh, empty per-symbol counter store -- `allocate_unique_atom_name`
    lazily creates each symbol's own counter state the first time it sees it.
    """
    return {}


def allocate_unique_atom_name(symbol: str, counters: dict[str, list[int]]) -> str:
    """Allocate the next unique <=4-char atom name for `symbol`.

    `counters` is mutated in place so repeated calls with the same dict
    keep producing globally-unique names across however many atoms share
    it (e.g. a whole msld system's core + every fragment, or a single
    standalone molecule) -- but each symbol draws from its own counter
    state, so exhausting one element's names doesn't affect any other
    symbol's remaining budget. Single-character element symbols (C, N, O,
    H, ...) get `{symbol}{n:03d}` (max 999 per symbol); two-character
    symbols (Cl, Br, ...) and the literal "LP" get a letter+digit suffix
    cycling through up to 260 combinations (also per symbol).
    Always uppercases `symbol` first: CHARMM atom names are conventionally
    all-uppercase, and CGenFF's own output re-uppercases two-letter element
    symbols (e.g. "Cl" -> "CL") regardless of what case it's given in the
    input mol2 -- generating a mixed-case name here would silently stop
    matching what CGenFF echoes back in its .rtf output for that atom.
    """
    symbol = symbol.upper()
    name = "LP" if symbol.startswith("LP") else symbol
    if len(name) == 1:
        state = counters.setdefault(name, [1])
        n = state[0]
        if n > 999:
            msg = f"Too many `{name}` atoms (max 999)."
            logger.error(msg)
            raise NamingError(msg)
        state[0] += 1
        return f"{name}{n:03d}"
    if len(name) == 2:
        state = counters.setdefault(name, [65, 48])  # [letter_ord, digit_ord]
        if state[0] > 90:
            state[0] = 65
            state[1] += 1
        if state[1] > 57:
            msg = f"Too many `{name}` atoms (max 260)."
            logger.error(msg)
            raise NamingError(msg)
        newat = f"{name}{chr(state[1])}{chr(state[0])}"
        state[0] += 1
        return newat
    msg = f"Cannot allocate a <=4-char atom name for symbol `{symbol}`."
    logger.error(msg)
    raise NamingError(msg)
