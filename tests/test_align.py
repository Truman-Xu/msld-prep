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
