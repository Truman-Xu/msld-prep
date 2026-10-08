#! /usr/bin/env python
"""Regenerate msldprep/sybyl_typing_data.py from OpenBabel's own
data/atomtyp.txt and data/types.txt.

`msldprep/sybyl_typer.py` used to parse those two text files itself on
every process start (comment-stripping, tokenizing, a `^N` hybridization-
query regex extraction, first-occurrence-wins table scans). None of that
depends on which molecule is being typed, so it's pulled out into this
one-time build step instead: run this script whenever `data/atomtyp.txt`/
`data/types.txt` change (e.g. a newer OpenBabel release), and it rewrites
`msldprep/sybyl_typing_data.py` as a plain, directly-importable Python
module -- no file I/O or text parsing left in the runtime typing path.

Usage: python rewrite/scripts/generate_sybyl_typing_data.py
"""

from __future__ import annotations

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
ATOMTYP_PATH = DATA_DIR / "atomtyp.txt"
TYPES_PATH = DATA_DIR / "types.txt"
OUTPUT_PATH = REPO_ROOT / "msldprep" / "sybyl_typing_data.py"

_HYB_QUERY_RE = re.compile(r"\^(\d+)")
# Only a standalone '#...' (preceded by whitespace or start-of-line) is a
# comment -- SMARTS atomic-number primitives like "[#7^3]" also use '#',
# always immediately after a non-whitespace character ('[').
_COMMENT_RE = re.compile(r"(?<!\S)#.*$")


def _strip_comment(line: str) -> str:
    return _COMMENT_RE.sub("", line).strip()


def _extract_outer_hyb_query(smarts: str) -> tuple[str, int | None]:
    match = _HYB_QUERY_RE.search(smarts)
    if match is None:
        return smarts, None
    return _HYB_QUERY_RE.sub("", smarts, count=1), int(match.group(1))


def parse_atomtyp(path: Path) -> tuple[list[tuple[str, int]], list[tuple[str, str, int | None, str]]]:
    """Returns (inthyb_rules, exttyp_rules), both in file order:
    inthyb_rules: [(smarts, hyb_value), ...]
    exttyp_rules: [(original_smarts, cleaned_smarts_for_rdkit, required_hyb, ext_type_name), ...]
    """
    inthyb: list[tuple[str, int]] = []
    exttyp: list[tuple[str, str, int | None, str]] = []
    for raw_line in path.read_text().splitlines():
        line = _strip_comment(raw_line)
        if not line:
            continue
        tokens = line.split()
        if len(tokens) < 3:
            continue
        category, smarts, value = tokens[0], tokens[1], tokens[2]
        if category == "INTHYB":
            inthyb.append((smarts, int(value)))
        elif category == "EXTTYP":
            cleaned, required_hyb = _extract_outer_hyb_query(smarts)
            exttyp.append((smarts, cleaned, required_hyb, value))
    return inthyb, exttyp


def parse_types_table(path: Path) -> tuple[dict[str, str], dict[int, str]]:
    """First-occurrence-wins INT->SYB and ATN->INT lookups, mirroring
    OpenBabel's `OBTypeTable::Translate` (a linear first-match scan).
    """
    lines = [line for line in path.read_text().splitlines() if line.strip() and not line.startswith("#")]
    header = lines[0].split()
    int_col = header.index("INT")
    atn_col = header.index("ATN")
    syb_col = header.index("SYB")

    int_to_syb: dict[str, str] = {}
    atn_to_int: dict[int, str] = {}
    for line in lines[1:]:
        fields = line.split()
        if len(fields) != len(header):
            continue
        int_name, atn, syb = fields[int_col], int(fields[atn_col]), fields[syb_col]
        int_to_syb.setdefault(int_name, syb)
        atn_to_int.setdefault(atn, int_name)
    return int_to_syb, atn_to_int


def _format_module(
    inthyb_rules: list[tuple[str, int]],
    exttyp_rules: list[tuple[str, str, int | None, str]],
    int_to_syb: dict[str, str],
    atn_to_int: dict[int, str],
) -> str:
    lines = [
        '"""Sybyl atom-typing rule tables, generated from OpenBabel\'s own',
        "data/atomtyp.txt and data/types.txt.",
        "",
        "DO NOT HAND-EDIT. Regenerate with:",
        "    python rewrite/scripts/generate_sybyl_typing_data.py",
        'after updating data/atomtyp.txt / data/types.txt. See sybyl_typer.py',
        "for how these are used.",
        '"""',
        "",
        "from __future__ import annotations",
        "",
        "# (smarts, hyb_value), in atomtyp.txt file order -- last SMARTS matched",
        "# for a given atom wins.",
        "INTHYB_RULES: list[tuple[str, int]] = [",
    ]
    for smarts, value in inthyb_rules:
        lines.append(f"    ({smarts!r}, {value!r}),")
    lines.append("]")
    lines.append("")
    lines.append(
        "# (original_smarts, cleaned_smarts_for_rdkit, required_hyb, sybyl_internal_type),"
    )
    lines.append(
        "# in atomtyp.txt file order -- last EXTTYP rule matched for a given atom"
    )
    lines.append(
        "# wins. `cleaned_smarts_for_rdkit` has any `^N` hybridization query"
    )
    lines.append(
        "# stripped out (RDKit's SMARTS engine has no equivalent primitive);"
    )
    lines.append(
        "# `required_hyb` is that stripped value, re-checked in Python instead."
    )
    lines.append("EXTTYP_RULES: list[tuple[str, str, int | None, str]] = [")
    for smarts, cleaned, required_hyb, value in exttyp_rules:
        lines.append(f"    ({smarts!r}, {cleaned!r}, {required_hyb!r}, {value!r}),")
    lines.append("]")
    lines.append("")
    lines.append("# OpenBabel internal atom-type name -> Sybyl (SYB) type.")
    lines.append("INT_TO_SYB: dict[str, str] = {")
    for int_name, syb in int_to_syb.items():
        lines.append(f"    {int_name!r}: {syb!r},")
    lines.append("}")
    lines.append("")
    lines.append(
        "# Atomic number -> OpenBabel internal atom-type name, used as a fallback"
    )
    lines.append(
        "# for atoms no EXTTYP rule matches (H, halogens, metals, ...)."
    )
    lines.append("ATN_TO_INT: dict[int, str] = {")
    for atn, int_name in atn_to_int.items():
        lines.append(f"    {atn!r}: {int_name!r},")
    lines.append("}")
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    inthyb_rules, exttyp_rules = parse_atomtyp(ATOMTYP_PATH)
    int_to_syb, atn_to_int = parse_types_table(TYPES_PATH)
    module_text = _format_module(inthyb_rules, exttyp_rules, int_to_syb, atn_to_int)
    OUTPUT_PATH.write_text(module_text)
    print(
        f"Wrote {OUTPUT_PATH} "
        f"({len(inthyb_rules)} INTHYB rules, {len(exttyp_rules)} EXTTYP rules, "
        f"{len(int_to_syb)} INT->SYB entries, {len(atn_to_int)} ATN->INT entries)"
    )


if __name__ == "__main__":
    main()
