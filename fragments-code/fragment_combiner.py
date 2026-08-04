"""Combine two fragments by forming a bond at H-bearing attachment points.

X1-H  +  H-X2  =>  X1-X2
"""

from itertools import product

from rdkit import Chem


def _attachment_atoms(mol):
    """Return indices of heavy atoms that have at least one available hydrogen."""
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetTotalNumHs() > 0]


def combine(mol1, mol2, idx1, idx2):
    """Bond atom idx1 of mol1 to atom idx2 of mol2, consuming one H from each.

    Returns the product Mol, or None if the resulting valence is invalid.
    """
    rw = Chem.RWMol(Chem.CombineMols(mol1, mol2))
    rw.AddBond(idx1, idx2 + mol1.GetNumAtoms(), Chem.BondType.SINGLE)
    try:
        Chem.SanitizeMol(rw)
    except (Chem.AtomValenceException, Chem.KekulizeException):
        return None
    return rw.GetMol()


def all_combinations(mol1, mol2):
    """Yield all unique products from coupling mol1 and mol2 at every H-bearing atom pair."""
    sites1 = _attachment_atoms(mol1)
    sites2 = _attachment_atoms(mol2)
    seen = set()
    for i, j in product(sites1, sites2):
        prod = combine(mol1, mol2, i, j)
        if prod is None:
            continue
        smi = Chem.MolToSmiles(prod)
        if smi not in seen:
            seen.add(smi)
            yield prod


if __name__ == "__main__":
    import sys

    smi1 = sys.argv[1] if len(sys.argv) > 1 else "c1cnccc1"   # pyridine
    smi2 = sys.argv[2] if len(sys.argv) > 2 else "c1cncnc1"   # pyrimidine

    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)
    if mol1 is None or mol2 is None:
        print("Invalid SMILES", file=sys.stderr)
        sys.exit(1)

    products = list(all_combinations(mol1, mol2))
    print(f"{len(products)} unique product(s) from {smi1!r} x {smi2!r}:")
    for p in products:
        print(" ", Chem.MolToSmiles(p))
