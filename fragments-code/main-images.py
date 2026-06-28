"""Fragment-based design: attach fragments from a library to molecular scaffolds."""

from rdkit import Chem
from rdkit.Chem import Draw, rdFreeSASA
from collections import defaultdict
from pathlib import Path
from typing import Optional

# ---------------------------------------------------------------------------
# Fragment library I/O
# ---------------------------------------------------------------------------

class FragmentLibrary:
    """Read and store a collection of fragments from a SMILES file.

    Expected file formats (one fragment per line):
        csv:  name,smiles          (first column = label, second = SMILES)
        tsv:  <tab>separated labels and smiles
        plain:  one SMILES per line, auto-numbered as frag_0, frag_1, ...
    """

    def __init__(self, path: str):
        self.path = path
        self.fragments: dict[str, Chem.Mol] = {}  # name -> mol

        p = Path(path)
        text = p.read_text().strip()

        if "," in text.splitlines()[0]:
            self._parse_csv(text)
        elif "\t" in text.splitlines()[0]:
            self._parse_tsv(text)
        else:
            self._parse_plain(text)

    # -- parsers -----------------------------------------------------------

    def _parse_csv(self, text: str) -> None:
        for lineno, line in enumerate(text.splitlines(), start=1):
            parts = [c.strip() for c in line.split(",")]
            if len(parts) < 2:
                continue
            name, smiles = parts[0], parts[1]
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                self.fragments[name] = mol

    def _parse_tsv(self, text: str) -> None:
        for lineno, line in enumerate(text.splitlines(), start=1):
            parts = [c.strip() for c in line.split("\t")]
            if len(parts) < 2:
                continue
            name, smiles = parts[0], parts[1]
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                self.fragments[name] = mol

    def _parse_plain(self, text: str) -> None:
        for lineno, line in enumerate(text.splitlines(), start=1):
            smiles = line.strip()
            if not smiles or smiles.startswith("#"):
                continue
            name = f"frag_{lineno}"
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                self.fragments[name] = mol

    # -- accessor ----------------------------------------------------------

    def get(self, name: str) -> Optional[Chem.Mol]:
        """Return the RDKit Mol for *name*, or ``None``."""
        return self.fragments.get(name)

    def names(self) -> list[str]:
        return list(self.fragments.keys())

    @property
    def smiles_dict(self) -> dict[str, str]:
        """{name: canonical SMILES}"""
        return {n: Chem.MolToSmiles(m) for n, m in self.fragments.items()}


def read_fragment_library(path: str) -> FragmentLibrary:
    """Convenience wrapper that returns a :class:`FragmentLibrary`."""
    return FragmentLibrary(path)


# ---------------------------------------------------------------------------
# Scaffold / fragment attachment helpers
# ---------------------------------------------------------------------------

def find_linker_sites(mol: Chem.Mol) -> list[int]:
    """Return atom indices that are good candidates for linker/fragment attachment.

    Heuristic (matches the style of generate_fragments.py): any heavy, non-ring
    atom with at least one hydrogen that can be replaced.
    """
    sites: list[int] = []
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 0:          # wildcard — skip
            continue
        if a.IsInRing():                    # ring atoms are less flexible
            continue
        if a.GetTotalNumHs() < 1:           # no H to replace
            continue
        sites.append(a.GetIdx())
    return sites


def add_fragment_to_scaffold(
    scaffold: Chem.Mol,
    fragment: Chem.Mol,
    frag_name: str,
    scaffold_atom_idx: int,
) -> tuple[Chem.Mol, str]:
    """Attach *fragment* to *scaffold* at *scaffold_atom_idx*.

    The attachment replaces a hydrogen on the scaffold with the first non-hydrogen
    atom of the fragment.  Returns ``(product_mol, product_smiles)``.

    Raises ``ValueError`` if the fragment has no suitable linker atom.
    """
    # --- pick attachment atom on fragment (first heavy, non-ring atom) ----
    frag_linker: int | None = None
    for a in fragment.GetAtoms():
        if a.GetAtomicNum() == 0:
            continue
        if not a.IsInRing():
            frag_linker = a.GetIdx()
            break

    if frag_linker is None:
        raise ValueError(
            f"Fragment '{frag_name}' has no non-ring atom to serve as a linker."
        )

    # --- replace H on scaffold with fragment ------------------------------
    rw = Chem.RWMol(scaffold)

    # Mark the scaffold attachment atom
    scaffold_atom = rw.GetAtomWithIdx(scaffold_atom_idx)

    # Remove any H on the scaffold site
    bonds_to_delete: list[tuple[int, int]] = []
    for bond in scaffold_atom.GetBonds():
        partner = bond.GetOtherAtomIdx(scaffold_atom_idx)
        partner_a = rw.GetAtomWithIdx(partner)
        if partner_a.GetAtomicNum() == 1:          # hydrogen
            bonds_to_delete.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))

    for b1, b2 in bonds_to_delete:
        rw.RemoveBond(b1, b2)

    # Insert fragment as a subgraph — new atoms start at scaffold.GetNumAtoms()
    num_before_insert = rw.GetNumAtoms()
    rw.InsertMol(fragment)

    # Attach first heavy non-ring atom of fragment to the scaffold site
    attach_atom_idx: int | None = None
    for a in fragment.GetAtoms():
        if a.GetAtomicNum() != 0 and not a.IsInRing():
            attach_atom_idx = num_before_insert + a.GetIdx()
            break

    if attach_atom_idx is None:
        raise ValueError(
            f"Fragment '{frag_name}' has no non-ring atom to serve as a linker."
        )

    rw.AddBond(scaffold_atom_idx, attach_atom_idx, Chem.BondType.SINGLE)

    product = rw.GetMol()
    Chem.SanitizeMol(product)
    smiles = Chem.MolToSmiles(product)
    return product, smiles


def scaffold_with_all_fragments(
    scaffold_smiles: str,
    frag_lib: FragmentLibrary,
) -> list[tuple[str, Chem.Mol, str]]:
    """For every fragment in *frag_lib*, try attaching it at every linker site on *scaffold*.

    Returns a list of ``(fragment_name, product_mol, product_smiles)`` tuples.
    Duplicates (same SMILES) are skipped.
    """
    scaffold = Chem.MolFromSmiles(scaffold_smiles)
    if scaffold is None:
        raise ValueError(f"Cannot parse scaffold SMILES: {scaffold_smiles}")

    sites = find_linker_sites(scaffold)
    seen: set[str] = set()
    results: list[tuple[str, Chem.Mol, str]] = []

    for name, frag in frag_lib.fragments.items():
        for site in sites:
            try:
                product, smiles = add_fragment_to_scaffold(
                    scaffold, frag, name, site
                )
            except ValueError:
                continue
            if smiles in seen:
                continue
            seen.add(smiles)
            results.append((name, product, smiles))

    return results


# ---------------------------------------------------------------------------
# CLI entry-point (run as `python main.py <scaffold_smiles> <fragments.tsv>`)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <scaffold_smiles> <fragment_file>")
        sys.exit(1)

    scaffold_smiles, frag_path = sys.argv[1], sys.argv[2]

    lib = read_fragment_library(frag_path)
    results = scaffold_with_all_fragments(scaffold_smiles, lib)

    print(f"Scaffold:  {scaffold_smiles}")
    print(f"Library:   {len(lib.fragments)} fragments ({frag_path})")
    print(f"Products:  {len(results)} unique adducts\n")

    for frag_name, product_mol, smiles in results[:20]:
        sa = find_linker_sites(product_mol)
        print(f"  {frag_name:20s} -> {smiles}")

    if len(results) > 20:
        print(f"  ... and {len(results) - 20} more")

# ---------------------------------------------------------------------------
# Export a grid image of all product molecules
# ---------------------------------------------------------------------------

    # Build list of (mol, label) pairs for rendering
    mols_with_labels: list[tuple[str, Chem.Mol]] = []
    for frag_name, product_mol, smiles in results:
        mols_with_labels.append((frag_name, product_mol))

    if mols_with_labels:
        img = Draw.MolsToGridImage(
            [m for _, m in mols_with_labels],
            molsPerRow=4,
            subImgSize=(300, 250),
            legends=[name for name, _ in mols_with_labels],
            returnPNG=False,          # return a PIL image instead of raw bytes
        )
        out_path = Path("products.png")
        img.save(out_path)
        print(f"\nImage saved to {out_path} ({len(mols_with_labels)} molecules)")
