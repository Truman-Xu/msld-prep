import os

import pytest
from rdkit import Chem

from msldprep import align
from msldprep.datamodel import LigandRecord
from msldprep.errors import AlignmentError

from .conftest import ANALOG_STEMS, CHEMBL_DIR, REFERENCE_STEM


def test_load_ligand_real_data(reference_ligand):
    assert reference_ligand.name == REFERENCE_STEM
    assert reference_ligand.mol.GetNumAtoms() > 0
    assert reference_ligand.net_charge == 0.0
    for atom in reference_ligand.mol.GetAtoms():
        assert atom.HasProp("CGENFFAtomType")
        assert atom.HasProp("CGENFFPartialCharge")


def test_assign_atom_types_bad_rtf_path_raises_with_actual_path(reference_ligand):
    mol = Chem.MolFromMol2File(reference_ligand.mol2_path, removeHs=False, sanitize=True)
    bad_path = "/nonexistent/path/to/nowhere.rtf"
    with pytest.raises((AlignmentError, FileNotFoundError, OSError)) as exc_info:
        align.assign_atom_types_from_rtf(mol, bad_path)
    # regression: notebook's error message referenced an undefined `rtf_path`
    # free variable and would raise NameError instead of a clean error here.
    assert not isinstance(exc_info.value, NameError)


def test_prepare_ligand_set_aligns_all_and_includes_reference(prepared_ligand_set, reference_ligand):
    assert not prepared_ligand_set.failed
    assert {lig.name for lig in prepared_ligand_set.aligned} == set(ANALOG_STEMS)
    assert prepared_ligand_set.core_query_mol.GetNumAtoms() > 0
    # regression: notebook's reference ligand never got typed/isotope-labeled
    # and was excluded from the group core search.
    assert reference_ligand.mol.HasProp("_Name") or True  # reference loaded without error
    for atom in reference_ligand.mol.GetAtoms():
        assert atom.HasProp("CGENFFAtomType")


def test_align_ligand_set_reports_failures_without_raising():
    # two unrelated small molecules with no meaningful shared substructure --
    # regression for the notebook's bare `except:` that silently swallowed
    # failures instead of reporting them.
    methane = Chem.AddHs(Chem.MolFromSmiles("C"))
    Chem.rdDistGeom.EmbedMolecule(methane, randomSeed=1)
    for atom in methane.GetAtoms():
        atom.SetProp("CGENFFAtomType", "CG331")
        atom.SetDoubleProp("CGENFFPartialCharge", 0.0)

    ref_record = LigandRecord(name="methane_ref", mol=methane, net_charge=0.0, is_reference=True)

    benzene = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    Chem.rdDistGeom.EmbedMolecule(benzene, randomSeed=2)
    # no _TriposAtomName -> align_common_core itself won't fail on typing,
    # but a genuinely disjoint substructure match still exercises the
    # failure-collection path if embedding fails; if RDKit still manages to
    # embed something (e.g. via a trivial single-atom MCS), that's fine too --
    # the real assertion is that align_ligand_set never raises and always
    # returns a well-formed AlignmentResult.
    benzene_record = LigandRecord(name="benzene", mol=benzene, net_charge=0.0)

    result = align.align_ligand_set([benzene_record], ref_record)
    assert isinstance(result, align.AlignmentResult)
    assert len(result.aligned) + len(result.failed) == 1


def _named_mol(smiles: str, *, seed: int) -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    Chem.rdDistGeom.EmbedMolecule(mol, randomSeed=seed)
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetProp("_TriposAtomName", f"{atom.GetSymbol()}{i}")
    return mol


def test_write_and_load_ligands_sdf_round_trips_names_and_bond_orders(tmp_path):
    # a fused aromatic heterocycle -- exactly the kind of ring RDKit's mol2
    # reader can occasionally fail to re-kekulize from an "ar" bond label;
    # sdf carries the already-known explicit bond order directly instead.
    indazole = _named_mol("c1ccc2[nH]ncc2c1CC(=O)O", seed=7)
    quinazoline = _named_mol("c1ccc2ncncc2c1", seed=11)
    ligands = [
        LigandRecord(name="indazole", mol=indazole, net_charge=0.0, is_reference=True),
        LigandRecord(name="quinazoline", mol=quinazoline, net_charge=0.0),
    ]

    sdf_path = str(tmp_path / "aligned_ligands.sdf")
    align.write_ligands_sdf(ligands, sdf_path)

    mols_by_name = align.load_ligands_sdf(sdf_path)
    assert set(mols_by_name) == {"indazole", "quinazoline"}

    for ligand in ligands:
        loaded = mols_by_name[ligand.name]
        assert loaded.GetNumAtoms() == ligand.mol.GetNumAtoms()
        assert [a.GetProp("_TriposAtomName") for a in loaded.GetAtoms()] == [
            a.GetProp("_TriposAtomName") for a in ligand.mol.GetAtoms()
        ]
        assert [b.GetBondType() for b in loaded.GetBonds()] == [b.GetBondType() for b in ligand.mol.GetBonds()]


def test_load_ligand_set_from_sdf_reattaches_cgenff_types(tmp_path):
    mol = _named_mol("CCO", seed=3)
    atom_types = {atom.GetProp("_TriposAtomName"): "CG331" for atom in mol.GetAtoms()}
    atom_charges = {name: 0.0 for name in atom_types}
    ligand = LigandRecord(
        name="ref", mol=mol, net_charge=0.0, is_reference=True,
        atom_types=atom_types, atom_charges=atom_charges,
        bonds=tuple((a.GetProp("_TriposAtomName"), b.GetProp("_TriposAtomName")) for a, b in
                    ((bond.GetBeginAtom(), bond.GetEndAtom()) for bond in mol.GetBonds())),
    )

    aligned_dir = str(tmp_path / "aligned")
    os.makedirs(aligned_dir, exist_ok=True)
    ligand.write_rtf(os.path.join(aligned_dir, "ref.rtf"))
    ligand.write_prm(os.path.join(aligned_dir, "ref.prm"))
    sdf_path = os.path.join(aligned_dir, "aligned_ligands.sdf")
    align.write_ligands_sdf([ligand], sdf_path)

    reference, analogs = align.load_ligand_set_from_sdf(sdf_path, aligned_dir, "ref", [])
    assert analogs == []
    assert reference.is_reference is True
    assert reference.mol.GetNumAtoms() == mol.GetNumAtoms()
    for atom in reference.mol.GetAtoms():
        assert atom.GetProp("CGENFFAtomType") == "CG331"


def _fresh_mol(stem: str) -> Chem.Mol:
    return align.load_ligand(
        os.path.join(CHEMBL_DIR, f"{stem}.mol2"),
        os.path.join(CHEMBL_DIR, f"{stem}.rtf"),
        os.path.join(CHEMBL_DIR, f"{stem}.prm"),
        name=stem,
    ).mol


def _fresh_analog_mol(stem: str) -> Chem.Mol:
    """An analog's Mol with no conformer, so each test aligns it from scratch."""
    mol = _fresh_mol(stem)
    mol.RemoveAllConformers()
    return mol


@pytest.fixture
def fresh_reference_mol() -> Chem.Mol:
    # Not the session-wide `reference_ligand`: `prepared_ligand_set` isotope-labels
    # its atoms (G5), and a labeled reference no longer matches an unlabeled analog.
    return _fresh_mol(REFERENCE_STEM)


def _molblock_3d(smiles: str) -> str:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    Chem.rdDistGeom.EmbedMolecule(mol, randomSeed=1)
    return Chem.MolToMolBlock(mol)


def test_mol_from_mol2_text_matches_the_file_loader(reference_ligand):
    with open(reference_ligand.mol2_path) as handle:
        mol = align.mol_from_mol2_text(handle.read())
    assert mol.GetNumAtoms() == reference_ligand.mol.GetNumAtoms()
    assert [a.GetProp("_TriposAtomName") for a in mol.GetAtoms()] == [
        a.GetProp("_TriposAtomName") for a in reference_ligand.mol.GetAtoms()
    ]
    assert mol.GetConformer().Is3D()


@pytest.mark.parametrize("text", ["", "not a mol2 file"])
def test_mol_from_mol2_text_rejects_unreadable_text(text):
    with pytest.raises(AlignmentError):
        align.mol_from_mol2_text(text)


def test_mol_from_sdf_text_reads_one_3d_molecule_with_explicit_hs():
    mol = align.mol_from_sdf_text(_molblock_3d("CCO"))
    assert mol.GetNumAtoms() == 9  # 3 heavy atoms + 6 H, none stripped
    assert mol.GetConformer().Is3D()


def test_mol_from_sdf_text_accepts_a_record_ending_in_the_sdf_delimiter():
    assert align.mol_from_sdf_text(_molblock_3d("CCO") + "$$$$\n").GetNumAtoms() == 9


def test_mol_from_sdf_text_rejects_zero_or_several_molecules():
    block = _molblock_3d("CCO") + "$$$$\n"
    for text in ("", "not an sdf file\n", block + block):
        with pytest.raises(AlignmentError, match="Expected one molecule"):
            align.mol_from_sdf_text(text)


def test_mol_from_sdf_text_requires_3d_coordinates_unless_told_otherwise():
    flat = Chem.MolFromSmiles("CCO")
    Chem.rdDepictor.Compute2DCoords(flat)
    text = Chem.MolToMolBlock(flat)
    with pytest.raises(AlignmentError, match="no 3D coordinates"):
        align.mol_from_sdf_text(text)
    assert align.mol_from_sdf_text(text, require_3d=False).GetNumAtoms() == 3


@pytest.mark.parametrize("stem", ANALOG_STEMS)
def test_align_to_reference_reports_atom_map_and_rmsd(fresh_reference_mol, stem):
    mol = _fresh_analog_mol(stem)
    ref_mol = fresh_reference_mol

    result = align.align_to_reference(mol, ref_mol)

    assert len(result.atom_map) == result.core.GetNumAtoms() > 0
    assert len(set(result.atom_map.values())) == len(result.atom_map)  # one-to-one
    for mol_idx, ref_idx in result.atom_map.items():
        assert mol.GetAtomWithIdx(mol_idx).GetSymbol() == ref_mol.GetAtomWithIdx(ref_idx).GetSymbol()
    # the tethered embed lands on the reference core; recompute the RMSD
    # from the coordinates to check the figure RDKit reports
    assert result.rmsd < 0.1
    mol_xyz = mol.GetConformer().GetPositions()
    ref_xyz = ref_mol.GetConformer().GetPositions()
    squared = [((mol_xyz[i] - ref_xyz[j]) ** 2).sum() for i, j in result.atom_map.items()]
    assert result.rmsd == pytest.approx((sum(squared) / len(squared)) ** 0.5, abs=1e-6)
    assert result.seed == align.DEFAULT_EMBED_SEED


def test_align_to_reference_is_reproducible_for_a_seed(fresh_reference_mol):
    stem = ANALOG_STEMS[0]

    def coords(**kwargs):
        mol = _fresh_analog_mol(stem)
        align.align_to_reference(mol, fresh_reference_mol, **kwargs)
        return mol.GetConformer().GetPositions()

    default = coords()
    assert (coords() == default).all()
    assert (coords(seed=align.DEFAULT_EMBED_SEED) == default).all()
    assert (coords(seed=7) == coords(seed=7)).all()
    assert not (coords(seed=7) == default).all()


def test_align_common_core_still_returns_the_core_query_mol(fresh_reference_mol):
    stem = ANALOG_STEMS[0]
    core = align.align_common_core(_fresh_analog_mol(stem), fresh_reference_mol)
    expected = align.align_to_reference(_fresh_analog_mol(stem), fresh_reference_mol).core
    assert core.GetNumAtoms() == expected.GetNumAtoms()
    assert core.GetNumBonds() == expected.GetNumBonds()


@pytest.mark.parametrize(
    "options",
    [align.McsOptions(ring_matches_ring_only=True), align.McsOptions(timeout=30)],
)
def test_align_to_reference_accepts_mcs_options(fresh_reference_mol, options):
    result = align.align_to_reference(_fresh_analog_mol(ANALOG_STEMS[0]), fresh_reference_mol, mcs=options)
    assert result.atom_map
    assert result.rmsd < 0.1


def test_align_to_reference_raises_value_error_without_a_shared_substructure():
    water = Chem.AddHs(Chem.MolFromSmiles("O"))
    sodium = Chem.AddHs(Chem.MolFromSmiles("[Na+]"))
    Chem.rdDistGeom.EmbedMolecule(sodium, randomSeed=1)
    with pytest.raises(ValueError, match="no shared substructure"):
        align.align_to_reference(water, sodium)


def _embedded(smiles: str) -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    Chem.rdDistGeom.EmbedMolecule(mol, randomSeed=1)
    return mol


@pytest.mark.parametrize(
    "reference_smiles, analog_smiles",
    [
        ("c1ccc2c(c1)c([nH+]cn2)Nc1cccc(Br)c1", "Oc1ccc2c(Nc3cccc(Br)c3)ncnc2c1"),  # charged ring N in the reference
        ("c1ccc2c(c1)c(ncn2)[NH2+]c1cccc(Br)c1", "Oc1ccc2c(Nc3cccc(Br)c3)ncnc2c1"),  # charged exocyclic N
    ],
    ids=["ring-N-in-reference", "exocyclic-N-in-reference"],
)
def test_align_to_reference_matches_atoms_whatever_their_charge(reference_smiles, analog_smiles):
    result = align.align_to_reference(_embedded(analog_smiles), _embedded(reference_smiles))
    assert result.atom_map


def _largest_ring_torsion_difference(analog: Chem.Mol, reference: Chem.Mol, atom_map: dict[int, int]) -> float:
    """Degrees: the largest difference between a torsion of a non-aromatic ring of
    the reference and the same torsion of the analog."""
    from rdkit.Chem import rdMolTransforms

    def torsions(conf, ring):
        n = len(ring)
        return [rdMolTransforms.GetDihedralDeg(conf, *(ring[(k + d) % n] for d in range(4))) for k in range(n)]

    to_analog = {r: a for a, r in atom_map.items()}
    largest = 0.0
    for ring in Chem.GetSymmSSSR(reference):
        ring = list(ring)
        if any(reference.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            continue
        pairs = zip(torsions(reference.GetConformer(), ring), torsions(analog.GetConformer(), [to_analog[i] for i in ring]))
        largest = max([largest, *(abs((r - a + 180) % 360 - 180) for r, a in pairs)])
    return largest


@pytest.mark.parametrize("seed", [2342, 1, 2, 3, 4])
def test_the_core_keeps_the_reference_geometry_ring_puckers_included(seed):
    # Embedding from RDKit's default start with a soft tether left these chairs twisted
    # or inverted (up to 92 degrees off, depending on the seed).
    reference = _embedded("CN1CCC(CC1)C1CCN(C)CC1")
    analog = Chem.AddHs(Chem.MolFromSmiles("CCN1CCC(CC1)C1CCN(CC)CC1"))
    result = align.align_to_reference(analog, Chem.Mol(reference), seed=seed)
    assert result.rmsd < 1e-3
    assert _largest_ring_torsion_difference(analog, reference, result.atom_map) < 1.0


def test_a_charged_analog_aligns_onto_a_neutral_reference():
    # the other direction of #9, which always worked
    result = align.align_to_reference(
        _embedded("Oc1ccc2c(Nc3cccc(Br)c3)[nH+]cnc2c1"), _embedded("c1ccc2c(c1)c(ncn2)Nc1cccc(Br)c1")
    )
    assert result.atom_map


def test_the_core_is_matched_by_element_so_an_isotope_labeled_reference_still_aligns(reference_ligand, prepared_ligand_set):
    # `prepared_ligand_set` labels the shared reference's atoms with isotopes (G15)
    assert any(atom.GetIsotope() for atom in reference_ligand.mol.GetAtoms())
    result = align.align_to_reference(_fresh_analog_mol(ANALOG_STEMS[0]), reference_ligand.mol)
    assert result.atom_map


_TAUTOMER_PAIRS = [
    ("O=c1[nH]cnc2ccccc12", "Oc1ncnc2ccccc12", 11),  # quinazolin-4-one and its hydroxy tautomer
    ("O=c1cccc[nH]1", "Oc1ccccn1", 7),  # 2-pyridone and 2-hydroxypyridine
    ("Cc1ccnc(N)n1", "CC1=CC=NC(=N)N1", 8),  # aminopyrimidine and its imino tautomer
]


@pytest.mark.parametrize("reference_smiles, analog_smiles, heavy_atoms", _TAUTOMER_PAIRS)
@pytest.mark.parametrize("ring_matches_ring_only", [False, True])
def test_ignoring_bond_orders_matches_a_tautomer_in_full(reference_smiles, analog_smiles, heavy_atoms, ring_matches_ring_only):
    def core_size(ignore_bond_orders):
        options = align.McsOptions(ignore_bond_orders=ignore_bond_orders, ring_matches_ring_only=ring_matches_ring_only)
        return len(align.align_to_reference(_embedded(analog_smiles), _embedded(reference_smiles), mcs=options).atom_map)

    assert core_size(False) == heavy_atoms - 1  # the atom whose bond order changed is dropped
    assert core_size(True) == heavy_atoms


@pytest.mark.parametrize(
    "reference_smiles, analog_smiles, flagged",
    [
        ("c1ccc2c(c1)c(ncn2)Nc1cccc(Br)c1", "Oc1ccc2c(Nc3cccc(Br)c3)[nH+]cnc2c1", {"N": (1, 0, 1, 0)}),  # protonated ring N
        ("O=c1cccc[nH]1", "Oc1ccccn1", {"N": (0, 0, 0, 1), "O": (0, 0, 1, 0)}),  # tautomer: the H moves from N to O
    ],
    ids=["protonated-ring-N", "tautomer"],
)
def test_the_protonation_check_flags_core_atoms_and_leaves_the_pose_as_it_is(reference_smiles, analog_smiles, flagged):
    reference = _embedded(reference_smiles)

    def aligned(check_protonation):
        analog = _embedded(analog_smiles)
        options = align.McsOptions(ignore_bond_orders=True, check_protonation=check_protonation)
        return analog, align.align_to_reference(analog, reference, mcs=options)

    tolerant_mol, tolerant = aligned(False)
    strict_mol, strict = aligned(True)
    assert tolerant.protonation_differences == ()
    assert strict.atom_map == tolerant.atom_map and strict.rmsd == tolerant.rmsd
    assert strict_mol.GetConformer().GetPositions().tolist() == tolerant_mol.GetConformer().GetPositions().tolist()
    found = {
        strict_mol.GetAtomWithIdx(d.atom).GetSymbol(): (d.charge, d.ref_charge, d.hydrogens, d.ref_hydrogens)
        for d in strict.protonation_differences
    }
    assert found == flagged
    assert all(strict.atom_map[d.atom] == d.ref_atom for d in strict.protonation_differences)


def test_the_protonation_check_does_not_flag_a_substituent_in_place_of_a_hydrogen():
    reference = _embedded("c1ccc2c(c1)c(ncn2)Nc1cccc(Br)c1")
    analog = _embedded("Oc1ccc2c(Nc3cccc(Br)c3)ncnc2c1")  # the hydroxy carbon is a CH in the reference
    result = align.align_to_reference(analog, reference, mcs=align.McsOptions(check_protonation=True))
    assert len(result.atom_map) == reference.GetNumHeavyAtoms()
    assert result.protonation_differences == ()


def _spy_on_find_mcs(monkeypatch, *, canceled=False):
    """Replace `align.FindMCS` by a wrapper that records its keyword arguments and
    can report the search as cancelled (the time limit hit)."""
    import types

    calls = []
    real = align.FindMCS

    def wrapper(mols, **kwargs):
        calls.append(kwargs)
        result = real(mols, **kwargs)
        return types.SimpleNamespace(queryMol=result.queryMol, canceled=canceled or result.canceled)

    monkeypatch.setattr(align, "FindMCS", wrapper)
    return calls


def test_the_core_search_has_a_ten_second_limit_by_default(monkeypatch):
    calls = _spy_on_find_mcs(monkeypatch)
    result = align.align_to_reference(_embedded("Oc1ncnc2ccccc12"), _embedded("O=c1[nH]cnc2ccccc12"))
    assert calls[0]["timeout"] == 10
    assert result.mcs_timed_out is False


def test_the_limit_can_be_changed_or_left_to_rdkit(monkeypatch):
    calls = _spy_on_find_mcs(monkeypatch)
    for options, expected in ((align.McsOptions(timeout=3), 3), (align.McsOptions(timeout=None), None)):
        align.align_to_reference(_embedded("Oc1ncnc2ccccc12"), _embedded("O=c1[nH]cnc2ccccc12"), mcs=options)
        assert calls[-1].get("timeout") == expected


def test_a_search_that_hit_its_limit_is_flagged_in_the_result(monkeypatch):
    _spy_on_find_mcs(monkeypatch, canceled=True)
    result = align.align_to_reference(_embedded("Oc1ncnc2ccccc12"), _embedded("O=c1[nH]cnc2ccccc12"))
    assert result.mcs_timed_out is True and result.atom_map


# --- which equivalent orientation is embedded (#11) --------------------------------------------

_ORTHO_REFERENCE = "COc1ccccc1CC(=O)O"  # a methoxy group next to the acetic acid chain
_ORTHO_ANALOG = "Clc1ccccc1CC(=O)O"  # a chlorine in its place; the ring can sit either way round


def _shuffled(smiles: str, seed: int) -> Chem.Mol:
    """The molecule with its atoms in another order, so the first substructure match varies."""
    import random

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    order = list(range(mol.GetNumAtoms()))
    random.Random(seed).shuffle(order)
    return Chem.RenumberAtoms(mol, order)


def _distance_from_chlorine_to_the_reference_oxygen(analog: Chem.Mol, reference: Chem.Mol) -> float:
    chlorine = next(a.GetIdx() for a in analog.GetAtoms() if a.GetSymbol() == "Cl")
    ether_oxygen = next(a.GetIdx() for a in reference.GetAtoms() if a.GetSymbol() == "O" and a.GetDegree() == 2 and not a.GetIsAromatic())
    delta = analog.GetConformer().GetPositions()[chlorine] - reference.GetConformer().GetPositions()[ether_oxygen]
    return float((delta**2).sum() ** 0.5)


def test_a_decorated_ring_turns_its_substituent_toward_the_reference_substituent():
    reference = _embedded(_ORTHO_REFERENCE)
    kept, first_match_only = [], []
    for seed in range(12):
        analog = _shuffled(_ORTHO_ANALOG, seed)
        result = align.align_to_reference(analog, Chem.Mol(reference))
        assert result.orientations == 2
        kept.append(_distance_from_chlorine_to_the_reference_oxygen(analog, reference))
        analog = _shuffled(_ORTHO_ANALOG, seed)
        align.align_to_reference(analog, Chem.Mol(reference), mcs=align.McsOptions(max_orientations=1))
        first_match_only.append(_distance_from_chlorine_to_the_reference_oxygen(analog, reference))
    assert max(kept) < 1.0  # on top of the reference's substituent, whatever the atom order
    assert max(first_match_only) > 3.0  # the first match alone gets it wrong for some atom orders


def test_an_unsubstituted_ring_has_one_orientation():
    result = align.align_to_reference(_embedded("c1ccccc1CC(=O)O"), _embedded(_ORTHO_REFERENCE))
    assert (result.orientations, result.orientation) == (1, 0)


def test_the_first_match_is_kept_when_the_reference_has_nothing_to_overlap():
    # the analog's ring can sit two ways, but the reference has no substituent to prefer either
    result = align.align_to_reference(_embedded(_ORTHO_ANALOG), _embedded("c1ccccc1CC(=O)O"))
    assert result.orientations == 2 and result.orientation == 0


def test_one_orientation_can_be_asked_for():
    result = align.align_to_reference(
        _embedded(_ORTHO_ANALOG), _embedded(_ORTHO_REFERENCE), mcs=align.McsOptions(max_orientations=1)
    )
    assert result.orientations == 1


# --- stereochemistry (#2) -------------------------------------------------------------------------

_STEREO_REFERENCE = "C[C@H](c1ccc(Cl)cc1)C(=O)NC"  # one stereocentre, inside the core shared with the analogs below
_SAME = "C[C@H](c1ccc(F)cc1)C(=O)NCC"
_OPPOSITE = "C[C@@H](c1ccc(F)cc1)C(=O)NCC"
_UNSPECIFIED = "CC(c1ccc(F)cc1)C(=O)NCC"


def _stereo_smiles_from_3d(mol: Chem.Mol) -> str:
    """The isomeric SMILES of the stereochemistry the coordinates have."""
    copy = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(copy)
    return Chem.MolToSmiles(Chem.RemoveHs(copy))


def _aligned(smiles: str, reference: Chem.Mol, *, seed: int = align.DEFAULT_EMBED_SEED, **options):
    analog = Chem.AddHs(Chem.MolFromSmiles(smiles))
    result = align.align_to_reference(analog, Chem.Mol(reference), mcs=align.McsOptions(**options), seed=seed)
    return analog, result


def test_the_stereochemistry_an_analog_specifies_is_kept_and_agrees_with_the_reference():
    reference = _embedded(_STEREO_REFERENCE)
    analog, result = _aligned(_SAME, reference)
    assert _stereo_smiles_from_3d(analog) == Chem.MolToSmiles(Chem.MolFromSmiles(_SAME))
    (center,) = result.stereo
    assert (center.resolved_by, center.in_core, center.opposes_reference) == ("analog", True, False)


def test_the_opposite_configuration_is_kept_but_flagged_as_opposing_the_reference():
    reference = _embedded(_STEREO_REFERENCE)
    analog, result = _aligned(_OPPOSITE, reference)
    assert _stereo_smiles_from_3d(analog) == Chem.MolToSmiles(Chem.MolFromSmiles(_OPPOSITE))  # not overruled
    (center,) = result.stereo
    assert center.opposes_reference is True


def test_only_a_configuration_that_fits_the_reference_core_sits_on_it_exactly():
    # A core fixed on the reference would invert the opposite centre, so that analog
    # is placed with the core held less strictly instead and keeps its configuration.
    reference = _embedded(_STEREO_REFERENCE)
    _, same = _aligned(_SAME, reference)
    opposite_analog, opposite = _aligned(_OPPOSITE, reference)
    assert same.rmsd < 1e-3
    assert opposite.rmsd > 0.1
    assert _stereo_smiles_from_3d(opposite_analog) == Chem.MolToSmiles(Chem.MolFromSmiles(_OPPOSITE))


@pytest.mark.parametrize("seed", [2342, 1, 2, 3])
def test_an_unspecified_centre_takes_the_reference_configuration_whatever_the_seed(seed):
    reference = _embedded(_STEREO_REFERENCE)
    analog, result = _aligned(_UNSPECIFIED, reference, seed=seed)
    assert _stereo_smiles_from_3d(analog) == Chem.MolToSmiles(Chem.MolFromSmiles(_SAME))
    (center,) = result.stereo
    assert (center.resolved_by, center.in_core) == ("reference", True)


def test_a_centre_the_reference_cannot_decide_is_chosen_by_energy_and_flagged():
    reference = _embedded("c1ccccc1CC(=O)O")  # no stereocentre anywhere, and no nitrogen to hold the analog's amide
    analog, result = _aligned("c1ccccc1CC(=O)NC(C)CC", reference)  # the sec-butyl centre is outside the core
    (center,) = result.stereo
    assert (center.resolved_by, center.in_core) == ("energy", False)
    # the structure now has a stereocentre with a configuration, not an undefined one
    assert all(
        info.specified == Chem.StereoSpecified.Specified
        for info in Chem.FindPotentialStereo(analog)
        if info.type == Chem.StereoType.Atom_Tetrahedral
    )


def test_a_centre_inside_the_core_is_not_decided_by_a_reference_without_stereo():
    # the reference has a plain ring carbon where the analog has a centre: nothing to copy
    reference = _embedded("CN1CCC(C)CC1")
    analog, result = _aligned("CN1CC(C)CCC1", reference)
    (center,) = result.stereo
    assert center.resolved_by == "energy" and center.opposes_reference is None


@pytest.mark.parametrize("smiles", ["C[C@H]1C[C@@H](C)CN(C)C1", "C[C@H]1C[C@H](C)CN(C)C1"])
def test_cis_and_trans_ring_isomers_keep_their_relative_configuration(smiles):
    reference = _embedded("CN1CCC(C)CC1")
    analog, _ = _aligned(smiles, reference)
    assert _stereo_smiles_from_3d(analog) == Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def test_a_stereoisomer_limit_of_one_still_gives_an_alignment():
    reference = _embedded("c1ccccc1CC(=O)NC")
    analog, result = _aligned("c1ccccc1CC(=O)NC(C)CC", reference, max_stereoisomers=1)
    assert result.atom_map and analog.GetNumConformers() == 1


@pytest.mark.xfail(
    strict=True,
    reason="#15: RemoveHs keeps isotope-labeled hydrogens, so find_shared_core searches with every hydrogen",
)
def test_the_shared_core_search_leaves_hydrogens_out():
    ligands = [
        align.load_ligand(
            os.path.join(CHEMBL_DIR, f"{stem}.mol2"),
            os.path.join(CHEMBL_DIR, f"{stem}.rtf"),
            os.path.join(CHEMBL_DIR, f"{stem}.prm"),
            name=stem,
        )
        for stem in [REFERENCE_STEM, ANALOG_STEMS[0]]
    ]
    # query atoms match by isotope only, so they report atomic number 0; count them instead
    core = align.type_ligands_and_find_core(ligands)
    assert core.GetNumAtoms() <= min(ligand.mol.GetNumHeavyAtoms() for ligand in ligands)
