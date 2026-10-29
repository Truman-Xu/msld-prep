"""Anchor/site/fragment decomposition off a shared core.

`align.find_shared_core` finds the heavy-atom maximum common substructure
across a ligand set, atoms matched by CGenFF type and neighbor types
(design note (e)), but stops there -- it doesn't say which atoms attach
where. This module does: for each ligand, match the shared core, add each
core heavy atom's hydrogens to it (design note (i)), find the core atoms
with a substituent hanging off them ("anchors"), and give each connected
non-core piece to the site of the anchors it is attached to.

Design choice (see plan): this skips RDKit's `RGroupDecomposition` and goes
straight to a single isotope-locked substructure match + shortest-path BFS
assignment (ported from msld_mcs_rdecomp.py's `r_group_decomposition`,
which legacy code only used as a *fallback* when RGD failed -- exactly the
symmetric-core situation our isotope-locked core tends to still hit). Site
identity across ligands is established once, structurally, from the shared
core_query_mol's own bond graph -- not by hoping each ligand's independent
anchor-detection order happens to line up.
"""

from __future__ import annotations

import itertools
import logging
import warnings
from typing import Sequence

from rdkit import Chem

from .datamodel import CoreDecomposition, Fragment, LigandRecord, Site
from .errors import DecompositionError

logger = logging.getLogger(__name__)


def match_core(
    mol: Chem.Mol,
    core_query: Chem.Mol,
    *,
    ligand_name: str,
    reference_pattern: tuple[bool, ...] | None = None,
) -> tuple[tuple[int, ...], bool]:
    """Return the atom-index tuple (in query-atom order) where `mol` matches `core_query`,
    plus whether that match was an ambiguous tie resolved by picking the first candidate.

    If more than one match exists (a genuinely symmetric core), picks the
    match whose off-core-neighbor pattern is closest to `reference_pattern`
    and warns; if still tied, warns and picks the first candidate instead of
    raising -- callers should surface `tied=True` to the user (e.g. via an
    annotated `io_mcs_text` export) since the choice is arbitrary.
    """
    matches = mol.GetSubstructMatches(core_query, uniquify=True)
    if not matches:
        msg = f"Ligand {ligand_name} does not contain the shared core substructure."
        logger.error(msg)
        raise DecompositionError(msg)
    if len(matches) == 1 or reference_pattern is None:
        return matches[0], False

    def off_core_pattern(match: tuple[int, ...]) -> tuple[bool, ...]:
        match_set = set(match)
        return tuple(
            any(nbr.GetIdx() not in match_set for nbr in mol.GetAtomWithIdx(idx).GetNeighbors())
            for idx in match
        )

    scored = sorted(
        matches,
        key=lambda m: sum(a != b for a, b in zip(off_core_pattern(m), reference_pattern)),
    )
    best, second = scored[0], scored[1]
    if off_core_pattern(best) == off_core_pattern(second):
        warnings.warn(
            f"Ligand {ligand_name} has {len(matches)} ambiguous core matches with "
            "identical anchor patterns; picked the first one arbitrarily. Review "
            "this ligand's CORE line in the exported MCS_for_MSLD.txt (flagged there "
            "with a comment).",
            stacklevel=2,
        )
        return best, True
    warnings.warn(
        f"Ligand {ligand_name} has {len(matches)} core matches; picked the one whose "
        "anchor pattern best matches the reference ligand.",
        stacklevel=2,
    )
    return best, False


def _anchor_query_positions(mol: Chem.Mol, heavy_indices: tuple[int, ...], core_set: set[int]) -> set[int]:
    """Query positions of the core heavy atoms that have a neighbor outside the core."""
    positions = set()
    for pos, atom_idx in enumerate(heavy_indices):
        neighbors = mol.GetAtomWithIdx(atom_idx).GetNeighbors()
        if any(nbr.GetIdx() not in core_set for nbr in neighbors):
            positions.add(pos)
    return positions


def _hydrogens_by_type(mol: Chem.Mol, atom_idx: int) -> dict[str, list[int]]:
    hydrogens: dict[str, list[int]] = {}
    for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
        if nbr.GetAtomicNum() == 1:
            hydrogens.setdefault(nbr.GetProp("CGENFFAtomType"), []).append(nbr.GetIdx())
    return hydrogens


def _with_core_hydrogens(
    mol: Chem.Mol, heavy_indices: tuple[int, ...], ref_mol: Chem.Mol, ref_heavy_indices: tuple[int, ...], *, ligand_name: str
) -> tuple[int, ...]:
    """The core of `mol`: each matched heavy atom followed by its hydrogens, in
    the reference's order for that position (design note (i)).

    The pattern rule guarantees a core heavy atom the same number of hydrogens,
    of the same types, as the reference's; hydrogens of one type are paired
    with the reference's by 3D distance (the ligands are aligned onto it).
    """
    ref_conf, conf = ref_mol.GetConformer(), mol.GetConformer()
    core: list[int] = []
    for atom_idx, ref_idx in zip(heavy_indices, ref_heavy_indices):
        core.append(atom_idx)
        own, ref_own = _hydrogens_by_type(mol, atom_idx), _hydrogens_by_type(ref_mol, ref_idx)
        if {t: len(h) for t, h in own.items()} != {t: len(h) for t, h in ref_own.items()}:
            msg = (
                f"Ligand {ligand_name}'s core atom {mol.GetAtomWithIdx(atom_idx).GetProp('_TriposAtomName')} "
                "has different hydrogens from the reference's atom at that core position."
            )
            logger.error(msg)
            raise DecompositionError(msg)
        ref_order = sorted(h for hs in ref_own.values() for h in hs)
        paired: dict[int, int] = {}
        for h_type, ref_hs in ref_own.items():
            best = min(
                itertools.permutations(own[h_type]),
                key=lambda perm: sum(
                    (conf.GetAtomPosition(h) - ref_conf.GetAtomPosition(r)).LengthSq() for h, r in zip(perm, ref_hs)
                ),
            )
            paired.update(zip(ref_hs, best))
        core.extend(paired[r] for r in ref_order)
    return tuple(core)


def _fragment_components(mol: Chem.Mol, core_set: set[int]) -> list[list[int]]:
    """The connected pieces of `mol` left when the core is taken out."""
    components: list[list[int]] = []
    seen: set[int] = set()
    for atom in mol.GetAtoms():
        start = atom.GetIdx()
        if start in core_set or start in seen:
            continue
        component, stack = [], [start]
        seen.add(start)
        while stack:
            idx = stack.pop()
            component.append(idx)
            for nbr in mol.GetAtomWithIdx(idx).GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx not in core_set and nbr_idx not in seen:
                    seen.add(nbr_idx)
                    stack.append(nbr_idx)
        components.append(sorted(component))
    return components


def _attached_positions(mol: Chem.Mol, component: list[int], heavy_indices: tuple[int, ...]) -> set[int]:
    """Query positions of the core heavy atoms bonded to `component`."""
    members = set(component)
    return {
        pos for pos, atom_idx in enumerate(heavy_indices)
        if any(nbr.GetIdx() in members for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors())
    }


def _group_anchor_positions(anchor_positions: set[int], joined: list[set[int]]) -> list[list[int]]:
    """Group anchor query-positions into sites, as legacy `msld_mcs.py` does:
    anchors share a site when one piece of non-core atoms is attached to all of
    them in some ligand (`joined`), so a ring cut by the core is one fragment.
    Separate substituents on neighboring core atoms are separate sites.
    Grouping is transitive and uses only query positions, so it is the same
    for every ligand by construction.
    """
    parent = {pos: pos for pos in anchor_positions}

    def find(pos: int) -> int:
        while parent[pos] != pos:
            parent[pos] = parent[parent[pos]]
            pos = parent[pos]
        return pos

    def union(a: int, b: int) -> None:
        parent[find(a)] = find(b)

    for positions in joined:
        first, *rest = sorted(positions)
        for pos in rest:
            union(first, pos)

    groups: dict[int, list[int]] = {}
    for pos in sorted(anchor_positions):
        groups.setdefault(find(pos), []).append(pos)
    return sorted(groups.values(), key=lambda group: group[0])


def assign_fragment_atoms(
    mol: Chem.Mol, core_local_indices: set[int], heavy_indices: tuple[int, ...], anchor_groups_query: list[list[int]]
) -> list[list[int]]:
    """Assign each non-core piece of `mol` to the site of the anchors it is attached to.

    Within a site, atoms are ordered by their shortest-path distance to the
    site's first anchor.
    """
    site_of_position = {pos: site_idx for site_idx, group in enumerate(anchor_groups_query) for pos in group}
    site_atoms: list[list[int]] = [[] for _ in anchor_groups_query]
    for component in _fragment_components(mol, core_local_indices):
        sites = {site_of_position[pos] for pos in _attached_positions(mol, component, heavy_indices)}
        if len(sites) != 1:
            msg = f"A fragment of {len(component)} atoms is attached to sites {sorted(sites)}, not to exactly one."
            logger.error(msg)
            raise DecompositionError(msg)
        site_atoms[sites.pop()].extend(component)

    for site_idx, atoms in enumerate(site_atoms):
        representative = heavy_indices[anchor_groups_query[site_idx][0]]
        site_atoms[site_idx] = sorted(
            atoms, key=lambda idx: (len(Chem.rdmolops.GetShortestPath(mol, idx, representative)), idx)
        )
    return site_atoms


def _fold_in_lone_pairs(
    ligand: LigandRecord,
    core_atom_order: tuple[str, ...],
    ligand_core_names: tuple[str, ...],
    fragments_by_site: list[Fragment],
) -> dict[str, str]:
    """Fold this ligand's LONEPAIR records into the decomposition.

    Fragment-hosted LPs are appended to the owning Fragment in place.
    Core-hosted LPs are returned as {reference_core_name: this_ligand's_lp_name},
    keyed by the *reference* ligand's name for that core position (so LP
    identity across ligands can be correlated positionally, the same way
    core_atom_names is correlated to core_atom_order).
    """
    host_to_ref_name = dict(zip(ligand_core_names, core_atom_order))
    name_to_fragment = {
        name: frag for frag in fragments_by_site for name in frag.atom_names
    }
    core_lp_names: dict[str, str] = {}
    for lp in ligand.lone_pairs:
        host = lp.host_atoms[0]
        if host in host_to_ref_name:
            core_lp_names[host_to_ref_name[host]] = lp.lp_atom
        elif host in name_to_fragment:
            name_to_fragment[host].lp_atom_names.append(lp.lp_atom)
        else:
            logger.warning(
                "Lone pair %s (ligand %s) hosted on %s, which is neither a core nor a "
                "fragment atom; dropped.",
                lp.lp_atom, ligand.name, host,
            )
    return core_lp_names


def decompose_ligand_set(
    ligands: Sequence[LigandRecord],
    core_query_mol: Chem.Mol,
    reference_name: str,
    *,
    core_id: int = 1,
) -> CoreDecomposition:
    """Decompose every ligand in `ligands` into core + per-site fragments, against `core_query_mol`.

    `core_id` is stamped onto the result as-is; it has no effect on the
    decomposition itself. Set it when `ligands` is one group out of several
    produced by `grouping.group_ligands_by_similarity`, so the resulting
    `CoreDecomposition` records which group it is (see io_mcs_text's "CORE
    ID" sections).
    """
    reference = next((lig for lig in ligands if lig.name == reference_name), None)
    if reference is None:
        msg = f"Reference ligand {reference_name} not found in ligand set."
        logger.error(msg)
        raise DecompositionError(msg)

    # `core_query_mol` has heavy atoms only; its matches are heavy-atom indices
    # in query order, and each core heavy atom's hydrogens are added after it.
    ref_heavy, _ = match_core(reference.mol, core_query_mol, ligand_name=reference.name)
    ref_indices = _with_core_hydrogens(reference.mol, ref_heavy, reference.mol, ref_heavy, ligand_name=reference.name)
    ref_anchor_positions = _anchor_query_positions(reference.mol, ref_heavy, set(ref_indices))
    ref_off_core_pattern = tuple(
        any(nbr.GetIdx() not in set(ref_heavy) for nbr in reference.mol.GetAtomWithIdx(idx).GetNeighbors())
        for idx in ref_heavy
    )
    core_atom_order = tuple(
        reference.mol.GetAtomWithIdx(idx).GetProp("_TriposAtomName") for idx in ref_indices
    )

    core_atom_names: dict[str, tuple[str, ...]] = {}
    core_lp_atom_names: dict[str, dict[str, str]] = {}
    tied_core_matches: list[str] = []

    # First pass: each ligand's core, and which anchors each of its non-core pieces joins.
    matched: list[tuple[LigandRecord, tuple[int, ...], tuple[int, ...]]] = []
    joined: list[set[int]] = []
    for ligand in ligands:
        heavy, tied = match_core(
            ligand.mol, core_query_mol, ligand_name=ligand.name, reference_pattern=ref_off_core_pattern
        )
        if tied:
            tied_core_matches.append(ligand.name)
        indices = _with_core_hydrogens(ligand.mol, heavy, reference.mol, ref_heavy, ligand_name=ligand.name)
        anchor_positions = _anchor_query_positions(ligand.mol, heavy, set(indices))
        if anchor_positions != ref_anchor_positions:
            msg = (
                f"Ligand {ligand.name}'s substitution pattern ({sorted(anchor_positions)}) "
                f"doesn't match the reference's ({sorted(ref_anchor_positions)}); the shared "
                "core is inconsistent for this ligand."
            )
            logger.error(msg)
            raise DecompositionError(msg)
        matched.append((ligand, heavy, indices))
        joined.extend(
            _attached_positions(ligand.mol, component, heavy)
            for component in _fragment_components(ligand.mol, set(indices))
        )

    anchor_groups_query = _group_anchor_positions(ref_anchor_positions, joined)
    nsites = len(anchor_groups_query)
    sites = [Site(index=i) for i in range(nsites)]

    # Second pass: each ligand's fragments, one per site.
    for ligand, heavy, indices in matched:
        names = tuple(
            ligand.mol.GetAtomWithIdx(idx).GetProp("_TriposAtomName") for idx in indices
        )
        core_atom_names[ligand.name] = names

        anchor_groups_local = [
            [heavy[pos] for pos in group] for group in anchor_groups_query
        ]
        site_atom_indices = assign_fragment_atoms(ligand.mol, set(indices), heavy, anchor_groups_query)

        ligand_fragments = []
        for site_idx, atom_indices in enumerate(site_atom_indices):
            atom_names = [
                ligand.mol.GetAtomWithIdx(idx).GetProp("_TriposAtomName") for idx in atom_indices
            ]
            anchor_names = [
                ligand.mol.GetAtomWithIdx(idx).GetProp("_TriposAtomName")
                for idx in anchor_groups_local[site_idx]
            ]
            fragment = Fragment(
                ligand_name=ligand.name,
                site_index=site_idx,
                atom_indices=list(atom_indices),
                atom_names=atom_names,
                anchor_atom_names=anchor_names,
            )
            sites[site_idx].fragments[ligand.name] = fragment
            ligand_fragments.append(fragment)

        core_lps = _fold_in_lone_pairs(ligand, core_atom_order, names, ligand_fragments)
        for ref_name, lp_name in core_lps.items():
            core_lp_atom_names.setdefault(ref_name, {})[ligand.name] = lp_name

    logger.info(
        "Decomposed %d ligands into a %d-atom core across %d sites",
        len(ligands), len(core_atom_order), nsites,
    )

    return CoreDecomposition(
        reference_name=reference_name,
        core_query_mol=core_query_mol,
        core_atom_order=core_atom_order,
        core_atom_names=core_atom_names,
        core_lp_atom_names=core_lp_atom_names,
        sites=sites,
        tied_core_matches=tied_core_matches,
        core_id=core_id,
    )


def _validated_site_lists(
    value: list[list[str]] | None, nsites: int, param_name: str
) -> list[list[str]]:
    if value is not None and len(value) == nsites:
        return value
    if value is not None:
        warnings.warn(
            f"{param_name} has {len(value)} entries but there are {nsites} sites; "
            "falling back to an empty default for all sites.",
            stacklevel=3,
        )
    return [[] for _ in range(nsites)]


def apply_manual_overrides(
    decomposition: CoreDecomposition,
    ligands: Sequence[LigandRecord],
    in_frag: list[list[str]] | None = None,
    an_core: list[list[str]] | None = None,
) -> CoreDecomposition:
    """Move atoms between core and fragments per site.

    By default, every site's anchor atom(s) (in the reference ligand) are
    moved OUT of the core and into each ligand's fragment at that site --
    matching legacy's default behavior (`msld_crn.py:250-254`), since the
    anchor position is exactly where chemistry is most likely to vary
    between substituents. `an_core` keeps specific reference-named atoms in
    the core instead; `in_frag` moves additional reference-named atoms out.
    Bonded hydrogens of any core atom being moved are dragged along with it.
    """
    nsites = len(decomposition.sites)
    in_frag = _validated_site_lists(in_frag, nsites, "in_frag")
    an_core = _validated_site_lists(an_core, nsites, "an_core")

    ref_mol = next(lig.mol for lig in ligands if lig.name == decomposition.reference_name)
    name_to_idx = {
        atom.GetProp("_TriposAtomName"): atom.GetIdx()
        for atom in ref_mol.GetAtoms()
        if atom.HasProp("_TriposAtomName")
    }

    droplist: list[list[str]] = []
    for site in decomposition.sites:
        ref_frag = site.fragments.get(decomposition.reference_name)
        droplist.append(list(ref_frag.anchor_atom_names) if ref_frag else [])

    for site_idx in range(nsites):
        for name in in_frag[site_idx]:
            if name not in droplist[site_idx]:
                droplist[site_idx].append(name)
        for name in an_core[site_idx]:
            if name in droplist[site_idx]:
                droplist[site_idx].remove(name)

    for site_idx in range(nsites):
        for name in list(droplist[site_idx]):
            if name.startswith("H"):
                continue
            idx = name_to_idx.get(name)
            if idx is None:
                continue
            for neighbor in ref_mol.GetAtomWithIdx(idx).GetNeighbors():
                if neighbor.GetSymbol() != "H" or not neighbor.HasProp("_TriposAtomName"):
                    continue
                h_name = neighbor.GetProp("_TriposAtomName")
                if h_name in decomposition.core_atom_order and h_name not in droplist[site_idx]:
                    droplist[site_idx].append(h_name)

    drop_positions_by_site = [
        [decomposition.core_atom_order.index(name) for name in names] for names in droplist
    ]
    all_drop_positions = sorted({pos for positions in drop_positions_by_site for pos in positions})
    if not all_drop_positions:
        return decomposition

    new_core_atom_order = tuple(
        name for i, name in enumerate(decomposition.core_atom_order) if i not in all_drop_positions
    )
    new_core_atom_names = {
        ligand_name: tuple(
            name for i, name in enumerate(names) if i not in all_drop_positions
        )
        for ligand_name, names in decomposition.core_atom_names.items()
    }

    ligand_by_name = {lig.name: lig for lig in ligands}
    new_sites = []
    for site in decomposition.sites:
        new_fragments = {}
        dropped_names = droplist[site.index]
        for ligand_name, fragment in site.fragments.items():
            extra_names = [
                decomposition.core_atom_names[ligand_name][pos]
                for pos in drop_positions_by_site[site.index]
            ]
            mol = ligand_by_name[ligand_name].mol
            name_to_local_idx = {
                atom.GetProp("_TriposAtomName"): atom.GetIdx()
                for atom in mol.GetAtoms()
                if atom.HasProp("_TriposAtomName")
            }
            extra_indices = [name_to_local_idx[name] for name in extra_names]
            extra_lp_names = [
                decomposition.core_lp_atom_names[ref_name][ligand_name]
                for ref_name in dropped_names
                if ligand_name in decomposition.core_lp_atom_names.get(ref_name, {})
            ]
            new_fragments[ligand_name] = Fragment(
                ligand_name=ligand_name,
                site_index=fragment.site_index,
                atom_indices=[*extra_indices, *fragment.atom_indices],
                atom_names=[*extra_names, *fragment.atom_names],
                lp_atom_names=[*extra_lp_names, *fragment.lp_atom_names],
                anchor_atom_names=list(fragment.anchor_atom_names),
            )
        new_sites.append(Site(index=site.index, fragments=new_fragments))

    all_dropped_names = {name for names in droplist for name in names}
    new_core_lp_atom_names = {
        ref_name: lp_map
        for ref_name, lp_map in decomposition.core_lp_atom_names.items()
        if ref_name not in all_dropped_names
    }

    logger.info(
        "apply_manual_overrides moved %d core atom(s) to fragments: %s",
        len(all_drop_positions), droplist,
    )

    return CoreDecomposition(
        reference_name=decomposition.reference_name,
        core_query_mol=decomposition.core_query_mol,
        core_atom_order=new_core_atom_order,
        core_atom_names=new_core_atom_names,
        core_lp_atom_names=new_core_lp_atom_names,
        sites=new_sites,
        tied_core_matches=decomposition.tied_core_matches,
        core_id=decomposition.core_id,
    )
