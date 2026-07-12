from __future__ import annotations

import sys
from typing import Dict, List

from rdkit import Chem

def add_imidazole_to_aromatic(mol: Chem.RWMol, carbon_idx_1: int, carbon_idx_2: int) -> Chem.Mol:
    """Expand an aromatic ring by inserting N–C–N between two adjacent carbons.

    Given two carbon atoms that are bonded to each other inside an aromatic ring,
    this function inserts a three-atom chain:

        C₁ — Nₐ — C_new — N_b — C₂

    The original carbons (``carbon_idx_1`` and ``carbon_idx_2``) must belong to an
    aromatic ring and be directly bonded.

    Args:
        mol:       An RDKit RWMol object containing an aromatic ring.
        carbon_idx_1: Index of the first ring carbon.
        carbon_idx_2: Index of the second ring carbon (must be bonded to ``carbon_idx_1``).

    Returns:
        The original Mol object with the expanded ring. The input molecule is modified.

    Raises:
        ValueError: If the indices are invalid or if the two carbons share no bond.
    """
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
