"""Split a ligand set into groups that each have a sane shared core, instead
of forcing one global `find_shared_core` MCS across every ligand at once.

Why this exists: `find_shared_core` needs a substructure common to literally
every ligand passed to it (isotope/CGenFF-type-locked, complete rings only).
If the set mixes simple substituent analogs with even a few structurally
different ("scaffold-hopped") analogs, the global MCS collapses to whatever
tiny fragment survives everywhere -- often a single anchor atom, i.e. one
degenerate "site" -- even though most of the set would happily share a much
larger, multi-site core on its own. Grouping by whole-molecule similarity
*before* running MCS keeps structurally compatible ligands together and
gives incompatible ones their own group, each with its own local core.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Sequence

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from .datamodel import LigandRecord

logger = logging.getLogger(__name__)

DEFAULT_SIMILARITY_THRESHOLD = 0.7
_MORGAN_RADIUS = 2

_fingerprint_generator = rdFingerprintGenerator.GetMorganGenerator(radius=_MORGAN_RADIUS)


@dataclass
class LigandGroup:
    """One core group: `reference` anchors core-finding/site identity for
    `analogs`, exactly like the top-level reference does for the whole set
    in the ungrouped case.
    """

    core_id: int
    reference: LigandRecord
    analogs: list[LigandRecord] = field(default_factory=list)


def _fingerprint(mol: Chem.Mol):
    # Morgan fingerprints are conventionally over the heavy-atom graph;
    # ligand Mols here carry explicit Hs (mol2 convention), so strip them
    # first -- same reasoning as align.find_shared_core's own Chem.RemoveHs.
    return _fingerprint_generator.GetFingerprint(Chem.RemoveHs(mol))


def group_ligands_by_similarity(
    reference: LigandRecord,
    analogs: Sequence[LigandRecord],
    *,
    threshold: float = DEFAULT_SIMILARITY_THRESHOLD,
) -> tuple[list[LigandGroup], list[LigandRecord]]:
    """Partition `analogs` into one or more `LigandGroup`s by Tanimoto
    similarity over Morgan (radius-2) fingerprints.

    `threshold <= 0` disables grouping entirely: returns a single group
    with `reference` and every analog, unchanged -- today's one-global-core
    behavior.

    Group 1 always pairs `reference` with every analog whose similarity to
    it is >= `threshold`. Remaining ("leftover") analogs are then clustered
    among themselves: repeatedly, whichever leftover analog is within
    threshold of the most *other* leftover analogs seeds a new group --
    becoming that group's own local reference for core-finding and anchor
    identity, since the original reference is by definition not a good
    anchor for ligands dissimilar to it -- pulling in every other leftover
    analog within `threshold` of the seed. This repeats until every analog
    belongs to some group.

    An analog left with no partner within threshold of anything remaining
    ends up alone; a single ligand has no "shared core" to find, so it's
    returned separately in the second tuple element (excluded from every
    group) rather than as a degenerate one-ligand group.
    """
    if threshold <= 0:
        return [LigandGroup(core_id=1, reference=reference, analogs=list(analogs))], []

    by_name = {lig.name: lig for lig in [reference, *analogs]}
    fps = {name: _fingerprint(lig.mol) for name, lig in by_name.items()}

    def similarity(name_a: str, name_b: str) -> float:
        return DataStructs.TanimotoSimilarity(fps[name_a], fps[name_b])

    pool = [lig.name for lig in analogs]
    primary_names = [name for name in pool if similarity(name, reference.name) >= threshold]
    primary_set = set(primary_names)
    pool = [name for name in pool if name not in primary_set]

    groups = [LigandGroup(core_id=1, reference=reference, analogs=[by_name[n] for n in primary_names])]
    excluded: list[LigandRecord] = []

    while pool:
        def _connectivity(name: str) -> int:
            return sum(1 for other in pool if other != name and similarity(name, other) >= threshold)

        seed_name = max(pool, key=_connectivity)
        pool.remove(seed_name)
        member_names = [name for name in pool if similarity(name, seed_name) >= threshold]
        for name in member_names:
            pool.remove(name)

        if member_names:
            groups.append(LigandGroup(
                core_id=len(groups) + 1,
                reference=by_name[seed_name],
                analogs=[by_name[n] for n in member_names],
            ))
        else:
            excluded.append(by_name[seed_name])

    if excluded:
        logger.warning(
            "%d analog(s) had no partner within the %.2f similarity threshold and were "
            "excluded from every core group (a single ligand has no shared core to find): %s",
            len(excluded), threshold, [lig.name for lig in excluded],
        )
    logger.info(
        "Grouped %d analog(s) + reference into %d core group(s) at similarity threshold %.2f: %s",
        len(analogs), len(groups), threshold,
        {g.core_id: (g.reference.name, len(g.analogs)) for g in groups},
    )
    return groups, excluded
