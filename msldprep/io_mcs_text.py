"""Optional manual-review checkpoint: export/import a CoreDecomposition as a
human-editable text file, format-compatible with legacy's `MCS_for_MSLD.txt`
(`writeMCS`, msld_mcs_rdecomp.py:707-741).

The whole msldprep pipeline is in-memory end to end and never needs this --
call these only if you want to hand-inspect or hand-edit the core/site/
fragment assignment (legacy's two-run workflow: generate, review, edit,
re-run) between decomposition and charge renormalization.

Legacy's `ANCHOR ATOMS` section is lossy: any anchor with more than one atom
is collapsed to the literal string "DUM" (confirmed against real output,
`examples/2_ChEMBL/MCS_for_MSLD.txt`), so it can't be losslessly read back.
`import_mcs_text` treats that section as advisory only and re-derives the
authoritative anchor atom name(s) structurally, from each ligand's actual
bond graph -- this also makes it robust to hand-edited CORE/SITE sections,
which is the entire point of the manual-review workflow.

Multiple core groups (`grouping.group_ligands_by_similarity` splitting one
ligand set into several independently-cored subsets) share one checkpoint
file: each group's NSUBS/REFLIG/CORE/ANCHOR ATOMS/SITE blocks are preceded
by a "CORE ID n" line. A checkpoint with only one group omits "CORE ID"
entirely and is byte-identical to the pre-grouping single-core format.
"""

from __future__ import annotations

import logging
from typing import Sequence

from rdkit import Chem

from .datamodel import CoreDecomposition, Fragment, LigandRecord, Site
from .errors import DecompositionError

logger = logging.getLogger(__name__)

_DUM = "DUM"


def export_mcs_text(
    decompositions: Sequence[CoreDecomposition],
    path: str,
    *,
    excluded: Sequence[str] = (),
    similarity_threshold: float | None = None,
) -> None:
    """Write `decompositions` to `path` in legacy `MCS_for_MSLD.txt` format.

    A single decomposition (the common, ungrouped case) is written exactly
    as before, with no "CORE ID" line. Multiple decompositions -- one per
    core group, as produced by `grouping.group_ligands_by_similarity` --
    each get their own "CORE ID n" section instead, preceded by a header
    comment recording the similarity threshold and any excluded ligand
    names (ligands with no partner within threshold of anything, dropped
    from every group; see `grouping.group_ligands_by_similarity`).
    """
    decompositions = sorted(decompositions, key=lambda d: d.core_id)
    multi_group = len(decompositions) > 1
    total_mols = len({name for d in decompositions for name in d.core_atom_names}) + len(excluded)

    lines = [
        "# Maximum Common Substructure Search for MSLD (msldprep)\n",
        f"# {total_mols} molecules processed\n\n",
    ]
    if multi_group or excluded:
        threshold_str = f"{similarity_threshold:.2f}" if similarity_threshold is not None else "?"
        lines.append(
            f"# split into {len(decompositions)} core group(s) by Tanimoto/Morgan similarity "
            f"(threshold {threshold_str})\n"
        )
        if excluded:
            lines.append(f"# excluded (no partner within threshold): {' '.join(excluded)}\n")
        lines.append("\n")

    for decomposition in decompositions:
        if multi_group:
            lines.append(f"CORE ID {decomposition.core_id}\n")
        lines.extend(_format_decomposition_block(decomposition))

    lines.append("END")

    with open(path, "w") as fp:
        fp.writelines(lines)


def _format_decomposition_block(decomposition: CoreDecomposition) -> list[str]:
    """NSUBS/REFLIG/CORE/ANCHOR ATOMS/SITE-n-FRAGMENTS lines for one decomposition (no END)."""
    mol_names = list(decomposition.core_atom_names.keys())
    nsubs = [len(site.fragments) for site in decomposition.sites]

    lines = [
        f"NSUBS {' '.join(str(n) for n in nsubs)}\n\n",
        f"REFLIG {decomposition.reference_name}\n\n",
        "CORE\n",
    ]
    tied = set(decomposition.tied_core_matches)
    for mol_name in mol_names:
        names = decomposition.core_atom_names[mol_name]
        line = f"{mol_name} {' '.join(names)}"
        if mol_name in tied:
            line += "  # TIE: ambiguous core match (identical anchor pattern); first match picked automatically"
        lines.append(line + "\n")

    lines.append("\nANCHOR ATOMS\n")
    for mol_name in mol_names:
        anchor_strs = []
        for site in decomposition.sites:
            fragment = site.fragments.get(mol_name)
            anchors = fragment.anchor_atom_names if fragment else []
            anchor_strs.append(anchors[0] if len(anchors) == 1 else _DUM)
        lines.append(f"{mol_name} {' '.join(anchor_strs)}\n")
    lines.append("\n")

    for site in decomposition.sites:
        lines.append(f"SITE {site.index + 1} FRAGMENTS\n")
        for mol_name in mol_names:
            fragment = site.fragments.get(mol_name)
            names = fragment.atom_names if fragment else []
            lines.append(f"{mol_name} {' '.join(names)}\n")
        lines.append("\n")

    return lines


def _split_into_groups(lines: list[str]) -> dict[int, list[str]]:
    """Split raw checkpoint lines on "CORE ID n" marker lines.

    No markers at all (the single-group, legacy-compatible format) -> the
    whole file is one implicit group, `{1: lines}`. Content before the
    first marker (the header comments) is discarded either way -- it's
    never part of a group's own NSUBS/REFLIG/CORE/... blocks.
    """
    markers: list[tuple[int, int]] = []
    for i, line in enumerate(lines):
        parts = line.strip().split()
        if parts[:2] != ["CORE", "ID"]:
            continue
        if len(parts) != 3 or not parts[2].isdigit():
            msg = f"Malformed `CORE ID` line (expected `CORE ID <n>`): {line!r}"
            logger.error(msg)
            raise DecompositionError(msg)
        markers.append((i, int(parts[2])))

    if not markers:
        return {1: lines}

    boundaries = [i for i, _ in markers] + [len(lines)]
    return {
        core_id: lines[start + 1 : end]
        for (start, core_id), end in zip(markers, boundaries[1:])
    }


def _parse_one_group(path: str, core_id: int, lines: list[str]) -> dict:
    mols: list[str] = []
    core: dict[str, list[str]] = {}
    sites_raw: list[dict[str, list[str]]] = []
    reflig = None
    nsubs: list[int] = []

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("NSUB"):
            nsubs = [int(x) for x in line.split()[1:]]
        elif line.startswith("REFL"):
            reflig = line.split()[1]
        elif line.startswith("CORE"):
            i += 1
            while i < len(lines) and lines[i].strip() and not lines[i].startswith("ANCH"):
                parts = lines[i].split("#", 1)[0].split()
                if parts:
                    mols.append(parts[0])
                    core[parts[0]] = parts[1:]
                i += 1
            continue
        elif line.startswith("ANCH"):
            i += 1
            while i < len(lines) and lines[i].strip() and not lines[i].startswith("SITE"):
                i += 1
            continue
        elif line.startswith("SITE"):
            site_data: dict[str, list[str]] = {}
            i += 1
            while i < len(lines) and lines[i].strip() and not lines[i].startswith(("SITE", "END")):
                parts = lines[i].split()
                site_data[parts[0]] = parts[1:]
                i += 1
            sites_raw.append(site_data)
            continue
        elif line.startswith("END"):
            break
        i += 1

    if reflig is None or reflig not in mols:
        msg = f"`{path}` (CORE ID {core_id}) has no valid REFLIG line matching a molecule in its CORE section."
        logger.error(msg)
        raise DecompositionError(msg)
    if len(sites_raw) != len(nsubs):
        msg = (
            f"`{path}` (CORE ID {core_id}): NSUBS declares {len(nsubs)} sites but "
            f"{len(sites_raw)} SITE blocks were found."
        )
        logger.error(msg)
        raise DecompositionError(msg)

    return {"mols": mols, "core": core, "sites": sites_raw, "reflig": reflig}


def _parse_sections(path: str) -> dict[int, dict]:
    with open(path) as fp:
        lines = fp.readlines()
    return {core_id: _parse_one_group(path, core_id, group_lines) for core_id, group_lines in _split_into_groups(lines).items()}


def _resolve_group(parsed: dict[int, dict], core_id: int | None, path: str) -> tuple[int, dict]:
    if core_id is not None:
        if core_id not in parsed:
            msg = f"`{path}` has no CORE ID {core_id} (available: {sorted(parsed)})."
            logger.error(msg)
            raise DecompositionError(msg)
        return core_id, parsed[core_id]
    if len(parsed) == 1:
        (only_id, only_group), = parsed.items()
        return only_id, only_group
    msg = f"`{path}` has {len(parsed)} core group(s) ({sorted(parsed)}); pass core_id= to pick one."
    logger.error(msg)
    raise DecompositionError(msg)


def read_core_ids(path: str) -> list[int]:
    """Every CORE ID present in `path`, sorted -- `[1]` for a single-group
    (legacy-format) checkpoint with no "CORE ID" markers at all. Lets a
    caller (e.g. the CLI's `build` command) discover how many core groups
    there are before deciding which one(s) to build.
    """
    with open(path) as fp:
        lines = fp.readlines()
    return sorted(_split_into_groups(lines))


def read_ligand_names(path: str, *, core_id: int | None = None) -> tuple[str, list[str]]:
    """(reference_name, analog_names) declared in one core group of an MCS
    checkpoint file -- lets a caller (e.g. the CLI's `build` command) know
    which ligands to load before calling `import_mcs_text`, without having
    to be told again on the command line.

    `core_id` selects which group when the checkpoint has more than one
    (see `read_core_ids`); optional (and ignored) when it has just the one.
    """
    parsed = _parse_sections(path)
    _, group = _resolve_group(parsed, core_id, path)
    reflig = group["reflig"]
    return reflig, [name for name in group["mols"] if name != reflig]


def import_mcs_text(path: str, ligands: Sequence[LigandRecord], *, core_id: int | None = None) -> CoreDecomposition:
    """Read a (possibly hand-edited) MCS text file back into a CoreDecomposition.

    Re-derives anchor atom names structurally from each ligand's Mol rather
    than trusting the file's (lossy) ANCHOR ATOMS section. `core_id` picks
    which core group to read when the checkpoint has more than one (see
    `read_core_ids`); optional (and ignored) when it has just the one.
    """
    parsed_all = _parse_sections(path)
    resolved_core_id, parsed = _resolve_group(parsed_all, core_id, path)
    ligand_by_name = {lig.name: lig for lig in ligands}

    missing = [name for name in parsed["mols"] if name not in ligand_by_name]
    if missing:
        msg = f"`{path}` (CORE ID {resolved_core_id}) references ligand(s) not present in `ligands`: {missing}."
        logger.error(msg)
        raise DecompositionError(msg)

    reference = ligand_by_name[parsed["reflig"]]
    core_atom_order = tuple(parsed["core"][parsed["reflig"]])
    core_atom_names = {name: tuple(parsed["core"][name]) for name in parsed["mols"]}

    for mol_name, names in core_atom_names.items():
        if len(names) != len(core_atom_order):
            msg = (
                f"`{path}` (CORE ID {resolved_core_id}): ligand {mol_name} has {len(names)} CORE atoms, "
                f"expected {len(core_atom_order)} (matching REFLIG {parsed['reflig']})."
            )
            logger.error(msg)
            raise DecompositionError(msg)
        ligand = ligand_by_name[mol_name]
        atom_names = {
            atom.GetProp("_TriposAtomName")
            for atom in ligand.mol.GetAtoms()
            if atom.HasProp("_TriposAtomName")
        }
        bad = [name for name in names if name not in atom_names and not name.startswith("LP")]
        if bad:
            msg = f"`{path}`: ligand {mol_name}'s CORE atom name(s) {bad} not found in its structure."
            logger.error(msg)
            raise DecompositionError(msg)

    sites = []
    for site_index, site_data in enumerate(parsed["sites"]):
        fragments: dict[str, Fragment] = {}
        for mol_name in parsed["mols"]:
            ligand = ligand_by_name[mol_name]
            names = site_data.get(mol_name, [])
            real_names = [n for n in names if not n.startswith("LP")]
            lp_names = [n for n in names if n.startswith("LP")]

            name_to_atom = {
                atom.GetProp("_TriposAtomName"): atom
                for atom in ligand.mol.GetAtoms()
                if atom.HasProp("_TriposAtomName")
            }
            bad = [n for n in real_names if n not in name_to_atom]
            if bad:
                msg = (
                    f"`{path}`: ligand {mol_name}'s SITE {site_index + 1} atom name(s) "
                    f"{bad} not found in its structure."
                )
                logger.error(msg)
                raise DecompositionError(msg)

            anchors = _derive_anchor_names(ligand.mol, set(core_atom_names[mol_name]), real_names)
            fragments[mol_name] = Fragment(
                ligand_name=mol_name,
                site_index=site_index,
                atom_indices=[name_to_atom[n].GetIdx() for n in real_names],
                atom_names=real_names,
                lp_atom_names=lp_names,
                anchor_atom_names=anchors,
            )
        sites.append(Site(index=site_index, fragments=fragments))

    core_query_mol = _rebuild_core_query_mol(reference.mol, core_atom_order)

    return CoreDecomposition(
        reference_name=parsed["reflig"],
        core_query_mol=core_query_mol,
        core_atom_order=core_atom_order,
        core_atom_names=core_atom_names,
        core_lp_atom_names={},
        sites=sites,
        core_id=resolved_core_id,
    )


def _derive_anchor_names(mol: Chem.Mol, core_names: set[str], fragment_names: list[str]) -> list[str]:
    name_to_atom = {
        atom.GetProp("_TriposAtomName"): atom
        for atom in mol.GetAtoms()
        if atom.HasProp("_TriposAtomName")
    }
    anchors: list[str] = []
    for name in fragment_names:
        atom = name_to_atom.get(name)
        if atom is None:
            continue
        for neighbor in atom.GetNeighbors():
            if not neighbor.HasProp("_TriposAtomName"):
                continue
            neighbor_name = neighbor.GetProp("_TriposAtomName")
            if neighbor_name in core_names and neighbor_name not in anchors:
                anchors.append(neighbor_name)
    return anchors


def _rebuild_core_query_mol(reference_mol: Chem.Mol, core_atom_order: Sequence[str]) -> Chem.Mol:
    """Reconstruct a core query submol from the reference ligand's core atoms.

    The text format has no SMARTS/query-mol field to round-trip, so this
    rebuilds the induced subgraph over `core_atom_order` from the reference
    ligand's own structure -- semantically equivalent for any downstream use
    (e.g. re-matching this decomposition against a new ligand later).
    """
    name_to_idx = {
        atom.GetProp("_TriposAtomName"): atom.GetIdx()
        for atom in reference_mol.GetAtoms()
        if atom.HasProp("_TriposAtomName")
    }
    core_indices = {name_to_idx[name] for name in core_atom_order if name in name_to_idx}
    bond_path = [
        bond.GetIdx()
        for bond in reference_mol.GetBonds()
        if bond.GetBeginAtomIdx() in core_indices and bond.GetEndAtomIdx() in core_indices
    ]
    return Chem.PathToSubmol(reference_mol, bond_path, useQuery=True)
