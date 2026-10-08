import filecmp
import os

from msldprep.writers import fragment_json, pdb, rtf as rtf_writer


def test_json_round_trip_matches_direct_writers(tmp_path, decomposition, all_ligands, crn_result, translation):
    """The JSON detour must be a pure serialization step -- writing files
    directly and writing files from a JSON round-trip of the same records
    must produce byte-identical output.
    """
    direct_dir = tmp_path / "direct"
    json_dir = tmp_path / "fromjson"
    direct_dir.mkdir()
    json_dir.mkdir()

    pdb.write_all_pdb_files(decomposition, all_ligands, translation, str(direct_dir), write_large_lig=False)
    rtf_writer.write_core_rtf(decomposition, all_ligands, crn_result, translation, str(direct_dir))
    rtf_writer.write_all_pres_rtf(decomposition, all_ligands, crn_result, translation, str(direct_dir))

    core_data, fragment_data = fragment_json.compute_all_file_data(decomposition, all_ligands, crn_result, translation)
    json_path = fragment_json.write_records_json(core_data, fragment_data, str(json_dir))
    fragment_json.write_files_from_json(json_path, str(json_dir))

    direct_files = {f for f in os.listdir(direct_dir)}
    json_files = {f for f in os.listdir(json_dir) if f != "records.json"}
    assert direct_files == json_files

    mismatches = []
    for name in direct_files:
        if not filecmp.cmp(direct_dir / name, json_dir / name, shallow=False):
            mismatches.append(name)
    assert mismatches == []


def test_records_json_is_valid_json_with_expected_keys(tmp_path, decomposition, all_ligands, crn_result, translation):
    import json

    core_data, fragment_data = fragment_json.compute_all_file_data(decomposition, all_ligands, crn_result, translation)
    path = fragment_json.write_records_json(core_data, fragment_data, str(tmp_path))
    with open(path) as fp:
        payload = json.load(fp)

    assert "core" in payload
    assert "fragments" in payload
    assert payload["core"]["source_name"] == decomposition.reference_name
    for site in decomposition.sites:
        for frag_number in range(1, len(site.fragments) + 1):
            key = f"site{site.index + 1}_sub{frag_number}"
            assert key in payload["fragments"]
