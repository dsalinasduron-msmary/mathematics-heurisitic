from rdkit import Chem
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
from rdkit.Chem import Draw, Recap, BRICS
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import Mol

from vina import Vina

mol = Chem.MolFromSmiles(smiles)
fragmented = BRICS.BRICSDecompose(mol,returnMols=True)

mk_prep = MoleculePreparation(
    merge_these_atom_types=("H"),
    charge_model="gasteiger",
)

def remove_wildcard(mol):
    mol = Chem.RWMol(mol)
    for a in [ a for a in mol.GetAtoms() if a.GetAtomicNum() == 0 ] :
        mol.RemoveAtom(a.GetIdx())
    return mol.GetMol()

def prepare_for_docking(molf):
    mol = remove_wildcard(molf)
    Chem.SanitizeMol(mol)
    mol = Chem.rdmolops.AddHs(mol)
    rdDistGeom.EmbedMolecule(mol)
    molsetup_list = mk_prep(mol)
    molsetup = molsetup_list[0]
    pdbqt_string = PDBQTWriterLegacy.write_string(molsetup)
    return mol, pdbqt_string[0]

test_mol = list(fragmented)[0]
print(Chem.MolToSmiles(test_mol))

#final_mol, docking_qts = prepare_for_docking(test_mol)
## dock!
#v = Vina(sf_name='vina', seed=99, verbosity=0)
#v.set_receptor('box_0.pdbqt')
#
#box_center = [16.6141,56.2758,12.7343]
#box_size = [19.1097,19.4189,16.1378]
#v.compute_vina_maps(center=box_center, box_size=box_size)
#v.set_ligand_from_string(docking_qts)
#v.dock(exhaustiveness=32)
#print(v.energies())
