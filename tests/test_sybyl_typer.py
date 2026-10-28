import os
import shutil
import subprocess

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from msldprep import sybyl_typer

# The obabel program (used here only as a validation oracle) is found through
# MSLDPREP_OBABEL_PATH, else on PATH.
OBABEL = os.environ.get("MSLDPREP_OBABEL_PATH") or shutil.which("obabel") or ""
requires_obabel = pytest.mark.skipif(
    not os.path.isfile(OBABEL),
    reason="obabel not found: set MSLDPREP_OBABEL_PATH or put obabel on PATH (used here only as a validation oracle)",
)


def _embed(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    cid = AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE)
    if cid < 0:
        AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE, useRandomCoords=True)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def _internal_and_sybyl_types(mol):
    internal = sybyl_typer.assign_internal_types(mol)
    sybyl = sybyl_typer.assign_sybyl_types(mol)
    return internal, sybyl


class TestSyntheticMolecules:
    """No external tools needed -- checks the rule engine directly against
    hand-derived expected values (also cross-checked against real obabel
    for the same molecules in TestAgainstRealObabel below)."""

    def test_ethane_carbons_are_sp3(self):
        mol = _embed("CC")
        _, sybyl = _internal_and_sybyl_types(mol)
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                assert sybyl[atom.GetIdx()] == "C.3"

    def test_benzene_carbons_are_aromatic(self):
        mol = _embed("c1ccccc1")
        _, sybyl = _internal_and_sybyl_types(mol)
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                assert sybyl[atom.GetIdx()] == "C.ar"

    def test_carboxylate_oxygens_are_o_co2(self):
        mol = _embed("CC(=O)[O-]")
        _, sybyl = _internal_and_sybyl_types(mol)
        oxygen_types = sorted(sybyl[a.GetIdx()] for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
        assert oxygen_types == ["O.co2", "O.co2"]

    def test_amide_carbon_and_nitrogen(self):
        mol = _embed("CC(=O)NC")
        _, sybyl = _internal_and_sybyl_types(mol)
        carbonyl_c = next(
            a for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and any(n.GetAtomicNum() == 8 for n in a.GetNeighbors())
        )
        amide_n = next(a for a in mol.GetAtoms() if a.GetAtomicNum() == 7)
        assert sybyl[carbonyl_c.GetIdx()] == "C.2"
        assert sybyl[amide_n.GetIdx()] == "N.am"

    def test_amide_bond_typed_am(self):
        mol = _embed("CC(=O)NC")
        internal, _ = _internal_and_sybyl_types(mol)
        carbonyl_c = next(
            a for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and any(n.GetAtomicNum() == 8 for n in a.GetNeighbors())
        )
        amide_n = next(a for a in mol.GetAtoms() if a.GetAtomicNum() == 7)
        bond = mol.GetBondBetweenAtoms(carbonyl_c.GetIdx(), amide_n.GetIdx())
        assert sybyl_typer.sybyl_bond_type(bond, internal) == "am"

    def test_aromatic_ring_bonds_typed_ar(self):
        mol = _embed("c1ccccc1")
        internal, _ = _internal_and_sybyl_types(mol)
        ring_bonds = [b for b in mol.GetBonds() if b.GetIsAromatic()]
        assert len(ring_bonds) == 6
        for bond in ring_bonds:
            assert sybyl_typer.sybyl_bond_type(bond, internal) == "ar"

    def test_guanidinium_central_carbon_is_cationic(self):
        mol = _embed("NC(=[NH2+])N")
        _, sybyl = _internal_and_sybyl_types(mol)
        central_c = next(a for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
        assert sybyl[central_c.GetIdx()] == "C.cat"

    def test_halogens_fall_back_to_element_symbol(self):
        mol = _embed("Clc1ccccc1")
        _, sybyl = _internal_and_sybyl_types(mol)
        cl = next(a for a in mol.GetAtoms() if a.GetAtomicNum() == 17)
        assert sybyl[cl.GetIdx()] == "Cl"

    def test_generic_hydrogen_is_h(self):
        mol = _embed("CC")
        _, sybyl = _internal_and_sybyl_types(mol)
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                assert sybyl[atom.GetIdx()] == "H"


class TestAnnotateSybylTypes:
    def test_sets_atom_and_bond_props(self):
        mol = _embed("c1ccccc1")
        sybyl_typer.annotate_sybyl_types(mol)
        for atom in mol.GetAtoms():
            assert atom.HasProp(sybyl_typer.ATOM_TYPE_PROP)
            expected = "C.ar" if atom.GetAtomicNum() == 6 else "H"
            assert atom.GetProp(sybyl_typer.ATOM_TYPE_PROP) == expected
        for bond in mol.GetBonds():
            assert bond.HasProp(sybyl_typer.BOND_TYPE_PROP)
            expected = "ar" if bond.GetIsAromatic() else "1"
            assert bond.GetProp(sybyl_typer.BOND_TYPE_PROP) == expected

    def test_props_match_the_pure_dict_functions(self):
        mol = _embed("CC(=O)NC")
        internal, sybyl = _internal_and_sybyl_types(mol)
        sybyl_typer.annotate_sybyl_types(mol)
        for atom in mol.GetAtoms():
            assert atom.GetProp(sybyl_typer.ATOM_TYPE_PROP) == sybyl[atom.GetIdx()]
        for bond in mol.GetBonds():
            assert bond.GetProp(sybyl_typer.BOND_TYPE_PROP) == sybyl_typer.sybyl_bond_type(bond, internal)


class TestWriteMol2File:
    def test_writes_valid_parseable_mol2_with_unique_names(self, tmp_path):
        mol = _embed("CCO")
        sybyl_typer.annotate_sybyl_types(mol)
        path = str(tmp_path / "ethanol.mol2")
        names = sybyl_typer.write_mol2_file(mol, path, mol_name="ethanol")

        assert len(set(names.values())) == mol.GetNumAtoms()
        text = open(path).read()
        assert "@<TRIPOS>MOLECULE" in text
        assert "@<TRIPOS>ATOM" in text
        assert "@<TRIPOS>BOND" in text

        reloaded = Chem.MolFromMol2File(path, removeHs=False)
        assert reloaded is not None
        assert reloaded.GetNumAtoms() == mol.GetNumAtoms()
        reloaded_names = {a.GetIdx(): a.GetProp("_TriposAtomName") for a in reloaded.GetAtoms()}
        assert set(reloaded_names.values()) == set(names.values())

    def test_raises_without_conformer(self, tmp_path):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        sybyl_typer.annotate_sybyl_types(mol)
        with pytest.raises(sybyl_typer.SybylTypingError):
            sybyl_typer.write_mol2_file(mol, str(tmp_path / "out.mol2"))

    def test_raises_if_not_annotated_first(self, tmp_path):
        mol = _embed("CCO")
        with pytest.raises(sybyl_typer.SybylTypingError):
            sybyl_typer.write_mol2_file(mol, str(tmp_path / "out.mol2"))


# Battery spanning the functional groups relevant to CGenFF-parameterizable
# drug-like small molecules: simple/aromatic/heteroaromatic carbons,
# carboxylic acids/esters/amides, guanidinium, sulfoxides/sulfones,
# phosphates, halogens, nitro/N-oxide/azide/isocyanate, boronic/sulfonic
# acids. Mirrors the battery used to develop this module against real
# obabel output during implementation.
BATTERY = {
    "ethane": "CC",
    "ethene": "C=C",
    "ethyne": "C#C",
    "benzene": "c1ccccc1",
    "toluene": "Cc1ccccc1",
    "phenol": "Oc1ccccc1",
    "aniline": "Nc1ccccc1",
    "pyridine": "c1ccncc1",
    "pyrimidine": "c1ccncn1",
    "furan": "c1ccoc1",
    "thiophene": "c1ccsc1",
    "pyrrole": "c1cc[nH]c1",
    "imidazole": "c1cnc[nH]1",
    "indole": "c1ccc2[nH]ccc2c1",
    "quinazoline": "c1ccc2ncncc2c1",
    "acetic_acid": "CC(=O)O",
    "acetate": "CC(=O)[O-]",
    "acetamide": "CC(=O)N",
    "n_methylacetamide": "CC(=O)NC",
    "methyl_acetate": "CC(=O)OC",
    "acetonitrile": "CC#N",
    "nitrobenzene": "c1ccccc1[N+](=O)[O-]",
    "guanidine": "NC(=N)N",
    "guanidinium": "NC(=[NH2+])N",
    "dimethyl_sulfoxide": "CS(=O)C",
    "dimethyl_sulfone": "CS(=O)(=O)C",
    "methanethiol": "CS",
    "thioanisole": "CSc1ccccc1",
    "trimethylamine": "CN(C)C",
    "ammonium": "[NH4+]",
    "trimethylphosphine_oxide": "CP(=O)(C)C",
    "phosphate": "OP(=O)(O)O",
    "chlorobenzene": "Clc1ccccc1",
    "bromobenzene": "Brc1ccccc1",
    "fluorobenzene": "Fc1ccccc1",
    "iodobenzene": "Ic1ccccc1",
    "dichloromethane": "ClCCl",
    "tert_butanol": "CC(C)(C)O",
    "diethyl_ether": "CCOCC",
    "cyclohexane": "C1CCCCC1",
    "acetone": "CC(=O)C",
    "acetaldehyde": "CC=O",
    "urea": "NC(=O)N",
    "carbamate": "COC(=O)N",
    "quinazoline_f_analog": "Fc1ccc(Nc2ncnc3ccccc23)cc1",
    "quinazoline_cl_analog": "Clc1ccc(Nc2ncnc3ccccc23)cc1",
    "quinazoline_me_analog": "Cc1ccc(Nc2ncnc3ccccc23)cc1",
    "pyridine_n_oxide": "[O-][n+]1ccccc1",
    "boronic_acid": "OB(O)c1ccccc1",
    "sulfonic_acid": "OS(=O)(=O)c1ccccc1",
    "isocyanate": "O=C=Nc1ccccc1",
    "azide": "CN=[N+]=[N-]",
    "carboxylate_amide_mix": "CC(=O)Nc1ccc(C(=O)O)cc1",
}


def _obabel_env():
    env = dict(os.environ)
    prefix = os.path.dirname(os.path.dirname(os.path.abspath(OBABEL)))
    lib_dirs = [os.path.join(prefix, "lib"), os.path.join(prefix, "lib64")]
    existing = env.get("LD_LIBRARY_PATH", "")
    env["LD_LIBRARY_PATH"] = ":".join([*lib_dirs, existing]) if existing else ":".join(lib_dirs)
    return env


def _obabel_mol2_types(mol: Chem.Mol, tmp_path, name: str):
    sdf_path = tmp_path / f"{name}.sdf"
    mol2_path = tmp_path / f"{name}_ob.mol2"
    with Chem.SDWriter(str(sdf_path)) as writer:
        writer.write(mol)
    result = subprocess.run(
        [OBABEL, str(sdf_path), "-O", str(mol2_path)], capture_output=True, text=True, env=_obabel_env(), timeout=60
    )
    assert result.returncode == 0, result.stderr

    lines = mol2_path.read_text().splitlines()
    atom_start = lines.index("@<TRIPOS>ATOM") + 1
    atom_types = []
    for line in lines[atom_start:]:
        if line.startswith("@<TRIPOS>"):
            break
        if line.strip():
            atom_types.append(line.split()[5])

    bond_start = lines.index("@<TRIPOS>BOND") + 1
    bond_types = {}
    for line in lines[bond_start:]:
        if not line.strip() or line.startswith("@"):
            break
        fields = line.split()
        i, j = int(fields[1]), int(fields[2])
        bond_types[(min(i, j), max(i, j))] = fields[3]
    return atom_types, bond_types


@requires_obabel
@pytest.mark.parametrize("name,smiles", sorted(BATTERY.items()))
def test_matches_real_obabel_typing(name, smiles, tmp_path):
    """Regression guard: this custom typer must keep reproducing real
    obabel's atom/bond typing across a broad functional-group battery,
    since it fully replaces the obabel CLI dependency in analog_prep.py.
    """
    mol = _embed(smiles)
    ob_atom_types, ob_bond_types = _obabel_mol2_types(mol, tmp_path, name)

    internal, sybyl = _internal_and_sybyl_types(mol)
    my_atom_types = [sybyl[i] for i in range(mol.GetNumAtoms())]
    assert my_atom_types == ob_atom_types

    my_bond_types = {}
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx() + 1, bond.GetEndAtomIdx() + 1
        my_bond_types[(min(i, j), max(i, j))] = sybyl_typer.sybyl_bond_type(bond, internal)
    assert my_bond_types == ob_bond_types


def _named(smiles: str, names: dict[int, str]) -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    for idx, name in names.items():
        mol.GetAtomWithIdx(idx).SetProp("_TriposAtomName", name)
    return mol


def test_assign_unique_atom_names_replaces_every_name_by_default():
    mol = _named("CCO", {0: "CA", 1: "CB", 2: "OG"})
    names = sybyl_typer.assign_unique_atom_names(mol)
    assert [names[i] for i in range(3)] == ["C001", "C002", "O001"]


def test_kept_atom_names_stay_and_new_ones_never_clash_with_them():
    # "C001" and "h001" are what fresh naming would hand out first
    mol = _named("CCO", {0: "C001", 1: "CB", 2: "h001"})
    names = sybyl_typer.assign_unique_atom_names(mol, keep_existing=True)
    assert [names[i] for i in range(3)] == ["C001", "CB", "h001"]
    hydrogens = [names[a.GetIdx()] for a in mol.GetAtoms() if a.GetSymbol() == "H"]
    assert hydrogens == ["H002", "H003", "H004", "H005", "H006", "H007"]
    assert all(a.GetProp("_TriposAtomName") == names[a.GetIdx()] for a in mol.GetAtoms())


@pytest.mark.parametrize(
    "names, message",
    [
        ({0: "CA", 1: "ca"}, r"used more than once: C1 \('CA'\), C2 \('ca'\)"),
        ({0: "CALPHA"}, r"longer than 4 characters: C1 \('CALPHA'\)"),
    ],
    ids=["duplicate-ignoring-case", "too-long"],
)
def test_kept_atom_names_that_cgenff_would_renumber_are_refused_with_the_atoms_named(names, message):
    with pytest.raises(sybyl_typer.SybylTypingError, match=message):
        sybyl_typer.assign_unique_atom_names(_named("CCO", names), keep_existing=True)
