import json
from rdkit import Chem


def write_sdf(smiles, path):
    mol = Chem.MolFromSmiles(smiles)
    mol.SetProp("_SMILES", smiles)
    name = f"{smiles[:8]}.sdf"
    out_path = path if path.endswith(".sdf") else f"{path}/{name}"
    with Chem.SDWriter(out_path) as writer:
        writer.write(mol)


if __name__ == "__main__":
    with open("tree.jsonl") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            record = json.loads(line)
            smiles = record["SMILES"]
            mol = Chem.MolFromSmiles(smiles)
            mol.SetProp("_SMILES", smiles)
            name = f"{record['functional_group']}_site{record['site']}"
            out_path = f"output_{name}.sdf"
            with Chem.SDWriter(out_path) as writer:
                writer.write(mol)
            print(f"Wrote {out_path}")

