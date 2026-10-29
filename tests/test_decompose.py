import pytest

from msldprep.errors import DecompositionError


def test_partition_property(all_ligands, raw_decomposition):
    """Every atom of every ligand is assigned to exactly the core or one fragment, never both/neither."""
    for ligand in all_ligands:
        all_idx = {a.GetIdx() for a in ligand.mol.GetAtoms()}
        core_names = set(raw_decomposition.core_atom_names[ligand.name])
        name_to_idx = {
            a.GetProp("_TriposAtomName"): a.GetIdx()
            for a in ligand.mol.GetAtoms()
            if a.HasProp("_TriposAtomName")
        }
        core_idx = {name_to_idx[n] for n in core_names if n in name_to_idx}

        frag_idx = set()
        for site in raw_decomposition.sites:
            frag_idx.update(site.fragments[ligand.name].atom_indices)

        assert all_idx - core_idx - frag_idx == set(), f"{ligand.name} has unassigned atoms"
        assert core_idx & frag_idx == set(), f"{ligand.name} has atoms in both core and a fragment"


def test_site_and_anchor_consistency_across_ligands(raw_decomposition):
    nsites = len(raw_decomposition.sites)
    assert nsites > 0
    for site in raw_decomposition.sites:
        anchor_sizes = {
            len(frag.anchor_atom_names) for frag in site.fragments.values()
        }
        assert len(anchor_sizes) == 1, f"site {site.index} has inconsistent anchor-group sizes: {anchor_sizes}"


def test_charge_sum_matches_ligand_net_charge(all_ligands, raw_decomposition):
    """Sum of (core + own fragments + own core-hosted LPs) original charges == this ligand's net charge."""
    for ligand in all_ligands:
        total = 0.0
        for name in raw_decomposition.core_atom_names[ligand.name]:
            total += ligand.atom_charges[name]
        for lp_by_ligand in raw_decomposition.core_lp_atom_names.values():
            if ligand.name in lp_by_ligand:
                total += ligand.atom_charges[lp_by_ligand[ligand.name]]
        for site in raw_decomposition.sites:
            fragment = site.fragments[ligand.name]
            for name in [*fragment.atom_names, *fragment.lp_atom_names]:
                total += ligand.atom_charges[name]
        assert round(total, 4) == round(ligand.net_charge, 4), (
            f"{ligand.name}: assembled charge {total} != net_charge {ligand.net_charge}"
        )


def test_apply_manual_overrides_moves_anchors_to_fragments(raw_decomposition, all_ligands, decomposition):
    # decomposition fixture = apply_manual_overrides(raw_decomposition, ...) with defaults
    assert len(decomposition.core_atom_order) < len(raw_decomposition.core_atom_order)
    ref_name = decomposition.reference_name
    raw_ref_anchors = {
        name for site in raw_decomposition.sites
        for name in site.fragments[ref_name].anchor_atom_names
    }
    assert raw_ref_anchors, "expected at least one anchor atom to have been identified"
    assert not (raw_ref_anchors & set(decomposition.core_atom_order)), (
        "default override should have moved every anchor atom out of the core"
    )


def test_decompose_raises_on_inconsistent_substitution_pattern(all_ligands, prepared_ligand_set, monkeypatch):
    """A ligand whose anchor pattern doesn't match the reference's should raise, not silently misalign."""
    import msldprep.decompose as decompose_mod

    real_positions = decompose_mod._anchor_query_positions
    call_count = {"n": 0}

    def fake_positions(*args):
        call_count["n"] += 1
        # corrupt the result for the 2nd ligand processed (not the reference)
        if call_count["n"] == 2:
            return set()
        return real_positions(*args)

    monkeypatch.setattr(decompose_mod, "_anchor_query_positions", fake_positions)
    with pytest.raises(DecompositionError):
        decompose_mod.decompose_ligand_set(all_ligands, prepared_ligand_set.core_query_mol, "11896")


def _heavy(names):
    return {name for name in names if not name.startswith("H")}


def test_the_core_follows_the_legacy_neighbor_type_rule(raw_decomposition):
    # An atom is core only if its CGenFF type and its neighbors' types match in every
    # ligand (design note (e)). The legacy msld_mcs.MsldMCS run on these four ligands
    # (2026-10-29) gives exactly this heavy-atom core for the reference; the benzo-ring
    # carbons whose substituents differ between ligands are left out.
    assert _heavy(raw_decomposition.core_atom_order) == {"N009", "C010", "N011", "C012", "C016", "C017"}


def test_core_hydrogens_follow_their_heavy_atom_and_no_site_is_only_hydrogens(all_ligands, raw_decomposition):
    # design note (i), #15
    for ligand in all_ligands:
        core = set(raw_decomposition.core_atom_names[ligand.name])
        for atom in ligand.mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                heavy = atom.GetNeighbors()[0].GetProp("_TriposAtomName")
                assert (atom.GetProp("_TriposAtomName") in core) == (heavy in core), ligand.name
        for site in raw_decomposition.sites:
            assert _heavy(site.fragments[ligand.name].atom_names), f"site {site.index} of {ligand.name} has only hydrogens"


def test_a_core_position_has_the_same_atom_type_in_every_ligand(all_ligands, raw_decomposition):
    by_name = {ligand.name: ligand for ligand in all_ligands}
    for position in range(len(raw_decomposition.core_atom_order)):
        types = {
            by_name[name].atom_types[names[position]] for name, names in raw_decomposition.core_atom_names.items()
        }
        assert len(types) == 1, (raw_decomposition.core_atom_order[position], types)


def test_sites_are_the_connected_fragments_legacy_finds(raw_decomposition):
    # The core cuts the benzo ring of the quinazoline, so C013-C015 is one piece attached at
    # C012 and C016: one site, not split by anchor. The legacy msld_mcs.MsldMCS run on these
    # four ligands (2026-10-29) gives exactly these two reference fragments.
    fragments = [set(site.fragments["11896"].atom_names) for site in raw_decomposition.sites]
    assert fragments == [
        {"C008", "N007", "H022", "C006", "C018", "H028", "C002", "C003", "H019", "C004", "H020", "C005", "H021", "I001"},
        {"C013", "H024", "C014", "H025", "C015", "H026"},
    ]
