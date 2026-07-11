"""Iterate over every atom in a SMILES string using RDKit."""

from __future__ import annotations

import sys
from typing import Dict, List

from rdkit import Chem


def atoms_by_element(smiles: str) -> Dict[str, int]:
    """Count the number of each element type in the given SMILES string.

    Example::

        >>> atoms_by_element("CCO")
        {'C': 2, 'O': 1}

    Args:
        smiles: A valid SMILES string representing a molecule.

    Returns:
        Dictionary mapping element name to its count in the molecule.

    Raises:
        ValueError: If the SMILES string is invalid or RDKit cannot parse it.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES string: {smiles!r}")

    counts: Dict[str, int] = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        counts[symbol] = counts.get(symbol, 0) + 1
    return counts


def iter_atoms(smiles: str) -> List[Dict]:
    """Iterate over every atom and return a list of per-atom dicts.

    For each atom the returned dict includes:
        - index:     RDKit 0-based atom index
        - symbol:    Element symbol (e.g. 'C', 'N', 'O')
        - degree:    Number of bonded neighbours
        - formal_charge: Formal charge on the atom
        - aromatic:  Whether the atom is in an aromatic ring
        - hybridization: String representation of hybridization state

    Args:
        smiles: A valid SMILES string representing a molecule.

    Returns:
        List of dicts, one per atom.

    Raises:
        ValueError: If the SMILES string is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        pretty_smiles = smiles[:60] + ("..." if len(smiles) > 60 else "")
        raise ValueError(f"Could not parse SMILES string: {pretty_smiles!r}")

    atoms = []
    for atom in mol.GetAtoms():
        atoms.append(
            {
                "index": atom.GetIdx(),
                "symbol": atom.GetSymbol(),
                "degree": atom.GetDegree(),
                "formal_charge": atom.GetFormalCharge(),
                "aromatic": atom.GetIsAromatic(),
                "hybridization": str(atom.GetHybridization()),
            }
        )
    return atoms


def _print_atoms(smiles: str) -> None:
    """Print a human-readable table of all atoms in *smiles*.

    Exits with code 1 when the SMILES is invalid.
    """
    counts = atoms_by_element(smiles)
    print(f"SMILES : {smiles}")
    print(f"Elements: {counts}")
    print()

    rows = iter_atoms(smiles)
    if not rows:
        print("No atoms found.")
        return

    header = f"{'Idx':>4}  {'Symbol':>6}  {'Degree':>6}  {'Charge':>7}  {'Aromatic':>8}  {'Hybridization':>15}"
    sep = "-" * len(header)
    print(sep)
    print(header)
    print(sep)
    for row in rows:
        print(
            f"{row['index']:4d}  {row['symbol']:>6s}  {row['degree']:6d}  "
            f"{row['formal_charge']:7d}  {str(row['aromatic']):>8s}  "
            f"{row['hybridization']:>15s}"
        )


def expand_ring(mol: Chem.RWMol, carbon_idx_1: int, carbon_idx_2: int) -> Chem.Mol:
    """Expand an aromatic ring by inserting N–C–N between two adjacent carbons.

    Given two carbon atoms that are bonded to each other inside an aromatic ring,
    this function removes their direct bond and inserts a three-atom chain::

        C₁ — Nₐ — C_new — N_b — C₂

    The original carbons (``carbon_idx_1`` and ``carbon_idx_2``) must belong to an
    aromatic ring and be directly bonded.  Aromaticity is allowed to relax – the
    function does **not** re-aromatise the ring afterwards.

    Args:
        mol:       An RDKit RWMol object containing an aromatic ring.
        carbon_idx_1: Index of the first ring carbon.
        carbon_idx_2: Index of the second ring carbon (must be bonded to ``carbon_idx_1``).

    Returns:
        The original Mol object with the expanded ring.  The original molecule is modified.

    Raises:
        ValueError: If the indices are invalid or if the two carbons share no bond.
    """
    # Work on a copy so the caller's molecule is untouched.
    Chem.SanitizeMol(mol)  # catch any latent issues early

    # ---- validate inputs --------------------------------------------------
    if carbon_idx_1 < 0 or carbon_idx_1 >= mol.GetNumAtoms():
        raise ValueError(
            f"carbon_idx_1={carbon_idx_1} is out of range "
            f"(molecule has {mol.GetNumAtoms()} atoms)"
        )
    if carbon_idx_2 < 0 or carbon_idx_2 >= mol.GetNumAtoms():
        raise ValueError(
            f"carbon_idx_2={carbon_idx_2} is out of range "
            f"(molecule has {mol.GetNumAtoms()} atoms)"
        )

    c1 = mol.GetAtomWithIdx(carbon_idx_1)
    c2 = mol.GetAtomWithIdx(carbon_idx_2)
    if not (c1.GetAtomicNum() == 6 and c2.GetAtomicNum() == 6):
        raise ValueError("Both target atoms must be carbons")

    # ---- insert Nₐ attached to C₁ ----------------------------------------
    na_idx = mol.AddAtom(Chem.Atom(7))  # nitrogen
    mol.AddBond(carbon_idx_1, na_idx, Chem.BondType.SINGLE)
    na = mol.GetAtomWithIdx(na_idx)
    na.SetFormalCharge(0)

    # ---- insert N_b attached to C₂ ---------------------------------------
    nb_idx = mol.AddAtom(Chem.Atom(7))  # nitrogen
    mol.AddBond(carbon_idx_2, nb_idx, Chem.BondType.SINGLE)
    nb = mol.GetAtomWithIdx(nb_idx)
    nb.SetFormalCharge(0)

    # ---- insert C_new attached to both nitrogens --------------------------
    cn_idx = mol.AddAtom(Chem.Atom(6))  # carbon
    mol.AddBond(na_idx, cn_idx, Chem.BondType.DOUBLE)
    mol.AddBond(nb_idx, cn_idx, Chem.BondType.SINGLE)
    c_new = mol.GetAtomWithIdx(cn_idx)
    c_new.SetFormalCharge(0)

    # Let RDKit assign implicit hydrogens and clean up.
    mol.UpdatePropertyCache()
    mol = Chem.AddHs(mol, addCoords=True)
    Chem.SanitizeMol(mol)
    return mol

def _render_mols(mols: list[Chem.Mol], name_prefix: str) -> None:
    """Render a list of RDKit Mols as PNG images using MolDraw2DCairo."""
    try:
        from rdkit.Chem import Draw
    except ImportError:
        print("(skipping rendering — rdkit.Chem.Draw not available)")
        return

    if not mols:
        return

    if len(mols) == 1:
        img = Draw.MolToImage(mols[0], size=(400, 300))
    else:
        img = Draw.MolsToGridImage(
            mols, molsPerRow=min(len(mols), 3),
            subImgSize=(400, 300), legends=[f"expansion #{i + 1}" for i in range(len(mols))],
        )

    out = f"{name_prefix}.png"
    img.save(out)
    print(f"Saved {out}")


def main() -> None:
    if len(sys.argv) < 2:
        # Default example molecule: aspirin
        smiles = "CC(=O)Oc1ccccc1C(=O)O"
    else:
        smiles = sys.argv[1]

    _print_atoms(smiles)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")

    out = expand_ring(Chem.RWMol(mol),0,1)
    _render_mols([out],"out")
    #expanded = add_imidazole_to_aromatic(Chem.RWMol(mol))
    #print(f"\nFound {len(expanded)} aromatic C–C bond(s) to expand.\n")

    #_render_mols(expanded, "expanded")


if __name__ == "__main__":
    main()
