"""Fragment-based design: attach fragments from a library to molecular scaffolds."""

from rdkit import Chem
from rdkit.Chem import Draw
from collections import defaultdict
from pathlib import Path
from typing import Optional

from fragment_combiner import combine

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
    """Return indices of heavy atoms with at least one H available for attachment."""
    return [a.GetIdx() for a in mol.GetAtoms()
            if a.GetAtomicNum() != 0 and a.GetTotalNumHs() > 0]


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
    frag_attach = next(
        (a.GetIdx() for a in fragment.GetAtoms() if a.GetTotalNumHs() > 0), None
    )
    if frag_attach is None:
        raise ValueError(f"Fragment '{frag_name}' has no H-bearing atom for attachment.")

    product = combine(scaffold, fragment, scaffold_atom_idx, frag_attach)
    if product is None:
        raise ValueError(
            f"Fragment '{frag_name}': invalid valence at scaffold atom {scaffold_atom_idx}."
        )

    return product, Chem.MolToSmiles(product)


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
        print(f"Usage: {sys.argv[0]} <scaffold_file> <fragment_file>")
        sys.exit(1)

    scaffold_path, frag_path = sys.argv[1], sys.argv[2]
    scaffold_smiles = Path(scaffold_path).read_text().strip()

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
