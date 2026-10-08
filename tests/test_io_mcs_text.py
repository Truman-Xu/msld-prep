import os

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
