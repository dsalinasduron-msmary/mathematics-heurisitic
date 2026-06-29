import json
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit.Chem import rdDistGeom


mk_prep = MoleculePreparation(
    merge_these_atom_types=("H"),
    charge_model="gasteiger",
)


def prepare_for_docking(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    rdDistGeom.EmbedMolecule(mol)
    molsetup_list = mk_prep(mol)
    molsetup = molsetup_list[0]
    pdbqt_string = PDBQTWriterLegacy.write_string(molsetup)
    return mol, pdbqt_string[0]


def main():
    with open("job.jsonl") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            job = json.loads(line)
            smiles = job["SMILES"]
            _, pdbqt = prepare_for_docking(smiles)
            print(pdbqt)


if __name__ == "__main__":
    main()
