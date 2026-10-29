import os

import pytest

from msldprep import io_mcs_text


def test_export_import_round_trip(tmp_path, raw_decomposition, all_ligands):
    path = str(tmp_path / "mcs_checkpoint.txt")
    io_mcs_text.export_mcs_text([raw_decomposition], path)
    assert os.path.exists(path)

    reimported = io_mcs_text.import_mcs_text(path, all_ligands)

    assert reimported.core_atom_order == raw_decomposition.core_atom_order
    assert reimported.core_atom_names == raw_decomposition.core_atom_names
    assert len(reimported.sites) == len(raw_decomposition.sites)


def test_round_trip_recovers_multi_atom_anchor(raw_decomposition, tmp_path, all_ligands):
    """Legacy's ANCHOR ATOMS section collapses multi-atom anchors to the
    literal string "DUM" -- confirm import_mcs_text still recovers the real
    names structurally rather than trusting that lossy field.
    """
    multi_atom_site = next(
        site for site in raw_decomposition.sites
        if any(len(f.anchor_atom_names) > 1 for f in site.fragments.values())
    )
    path = str(tmp_path / "mcs_checkpoint.txt")
    io_mcs_text.export_mcs_text([raw_decomposition], path)

    with open(path) as fp:
        text = fp.read()
    assert "DUM" in text  # confirms the lossy path is actually exercised

    reimported = io_mcs_text.import_mcs_text(path, all_ligands)
    reimported_site = reimported.sites[multi_atom_site.index]
    for ligand_name, fragment in multi_atom_site.fragments.items():
        assert set(reimported_site.fragments[ligand_name].anchor_atom_names) == set(fragment.anchor_atom_names)
        assert len(fragment.anchor_atom_names) > 1  # sanity: this really is the multi-atom case


def test_import_raises_on_unknown_atom_name(raw_decomposition, all_ligands, tmp_path):
    path = str(tmp_path / "mcs_checkpoint.txt")
    io_mcs_text.export_mcs_text([raw_decomposition], path)
    with open(path) as fp:
        text = fp.read()
    text = text.replace(raw_decomposition.core_atom_order[0], "ZZZZ999", 1)
    with open(path, "w") as fp:
        fp.write(text)

    from msldprep.errors import DecompositionError
    import pytest

    with pytest.raises(DecompositionError):
        io_mcs_text.import_mcs_text(path, all_ligands)


def test_tied_and_far_ligands_are_flagged_on_their_core_line_and_reimport_unchanged(raw_decomposition, all_ligands, tmp_path):
    import dataclasses

    analog = next(name for name in raw_decomposition.core_atom_names if name != raw_decomposition.reference_name)
    far_names = list(raw_decomposition.core_atom_names[analog][:2])
    flagged = dataclasses.replace(raw_decomposition, tied_core_matches=[analog], far_core_atoms={analog: far_names})
    path = str(tmp_path / "mcs_checkpoint.txt")
    io_mcs_text.export_mcs_text([flagged], path)

    with open(path) as fp:
        core_section = fp.read().split("\nCORE\n", 1)[1].split("\nANCHOR ATOMS\n", 1)[0]
    core_lines = {line.split()[0]: line for line in core_section.splitlines() if line.strip()}
    assert "# TIE: the core fits this ligand equally well in more than one way in 3D" in core_lines[analog]
    assert f"# FAR: {' '.join(far_names)} more than 1.1 A from the reference's atoms" in core_lines[analog]
    assert "#" not in core_lines[raw_decomposition.reference_name]

    reimported = io_mcs_text.import_mcs_text(path, all_ligands)
    assert reimported.core_atom_names == raw_decomposition.core_atom_names


@pytest.mark.xfail(
    strict=True,
    reason="export_mcs_text writes no lone-pair names and import_mcs_text reads them only from SITE lines, "
    "so the build step loses every lone pair (#18)",
)
def test_lone_pairs_survive_the_checkpoint_round_trip(raw_decomposition, all_ligands, tmp_path):
    def lone_pairs(decomposition):
        fragments = {(site.index, name): sorted(f.lp_atom_names) for site in decomposition.sites for name, f in site.fragments.items()}
        return fragments, decomposition.core_lp_atom_names

    before = lone_pairs(raw_decomposition)
    assert any(names for names in before[0].values())  # the ChEMBL halogens do carry lone pairs
    path = str(tmp_path / "mcs_checkpoint.txt")
    io_mcs_text.export_mcs_text([raw_decomposition], path)
    assert lone_pairs(io_mcs_text.import_mcs_text(path, all_ligands)) == before
