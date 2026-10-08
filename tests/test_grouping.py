from rdkit import Chem

from msldprep import grouping
from msldprep.datamodel import LigandRecord

# A biphenyl-sulfonamide-piperazine scaffold (mirrors the shape of a real
# MSLD analog set: one shared core, small substituent differences) plus two
# structurally unrelated families, chosen so their pairwise Tanimoto/Morgan
# similarities straddle the default 0.7 threshold in known ways (see the
# assertions below) -- confirmed against real RDKit output, not guessed.
_REF_SMILES = "c1ccc(cc1)c1ccc(cc1)S(=O)(=O)N1CCN(C)CC1"
_A_ME_SMILES = "Cc1ccc(cc1)c1ccc(cc1)S(=O)(=O)N1CCN(C)CC1"  # sim to ref: 0.758
_A_CL_SMILES = "Clc1ccc(cc1)c1ccc(cc1)S(=O)(=O)N1CCN(C)CC1"  # sim to ref: 0.735
_A_ET_SMILES = "c1ccc(cc1)c1ccc(cc1)S(=O)(=O)N1CCN(CC)CC1"  # sim to ref: 0.694 (just below)
_B1_SMILES = "CCCCCCCCCCCC"  # dodecane; sim to ref: 0.028, sim to _B2: 1.000
_B2_SMILES = "CCCCCCCCCCCCC"  # tridecane


def _ligand(name: str, smiles: str) -> LigandRecord:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    return LigandRecord(name=name, mol=mol, net_charge=0.0)


def _names(ligands):
    return {lig.name for lig in ligands}


def test_threshold_zero_disables_grouping():
    reference = _ligand("ref", _REF_SMILES)
    analogs = [_ligand("a_et", _A_ET_SMILES), _ligand("b1", _B1_SMILES), _ligand("b2", _B2_SMILES)]

    groups, excluded = grouping.group_ligands_by_similarity(reference, analogs, threshold=0)

    assert len(groups) == 1
    assert groups[0].core_id == 1
    assert groups[0].reference is reference
    assert _names(groups[0].analogs) == {"a_et", "b1", "b2"}
    assert excluded == []


def test_similar_analogs_stay_with_reference():
    reference = _ligand("ref", _REF_SMILES)
    analogs = [_ligand("a_me", _A_ME_SMILES), _ligand("a_cl", _A_CL_SMILES)]

    groups, excluded = grouping.group_ligands_by_similarity(reference, analogs, threshold=0.7)

    assert len(groups) == 1
    assert groups[0].reference is reference
    assert _names(groups[0].analogs) == {"a_me", "a_cl"}
    assert excluded == []


def test_leftover_analogs_form_their_own_group():
    reference = _ligand("ref", _REF_SMILES)
    analogs = [
        _ligand("a_me", _A_ME_SMILES), _ligand("a_cl", _A_CL_SMILES),
        _ligand("b1", _B1_SMILES), _ligand("b2", _B2_SMILES),
    ]

    groups, excluded = grouping.group_ligands_by_similarity(reference, analogs, threshold=0.7)

    assert excluded == []
    assert len(groups) == 2
    primary, secondary = groups
    assert primary.core_id == 1
    assert primary.reference is reference
    assert _names(primary.analogs) == {"a_me", "a_cl"}
    # b1/b2 are mutually identical fingerprints (both alkanes), unrelated to
    # ref -- they must cluster into a second group, seeded by one of them.
    assert secondary.core_id == 2
    assert {secondary.reference.name, *_names(secondary.analogs)} == {"b1", "b2"}


def test_lone_analog_with_no_partner_is_excluded():
    reference = _ligand("ref", _REF_SMILES)
    analogs = [_ligand("a_me", _A_ME_SMILES), _ligand("a_cl", _A_CL_SMILES), _ligand("a_et", _A_ET_SMILES)]

    groups, excluded = grouping.group_ligands_by_similarity(reference, analogs, threshold=0.7)

    assert len(groups) == 1
    assert _names(groups[0].analogs) == {"a_me", "a_cl"}
    # a_et is below threshold against both the reference and every other
    # leftover analog -- a single ligand has no shared core to find, so it's
    # excluded rather than becoming its own degenerate one-ligand group.
    assert _names(excluded) == {"a_et"}
