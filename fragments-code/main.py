"""Main entry-point: scaffold with fragments, then write each product to its own folder."""

import importlib.util
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Import from attach-fragments (hyphenated filename can't be a module)
# ---------------------------------------------------------------------------

spec = importlib.util.spec_from_file_location(
    "attach_fragments",
    str(Path(__file__).parent / "attach-fragments.py"),
)
_attach = importlib.util.module_from_spec(spec)
sys.modules["attach_fragments"] = _attach
spec.loader.exec_module(_attach)

FragmentLibrary = _attach.FragmentLibrary
scaffold_with_all_fragments = _attach.scaffold_with_all_fragments


def scaffold_and_save(
    scaffold_path: str,
    frag_lib_path: str,
    out_dir: str = "products",
) -> None:
    """Attach every fragment to every linker site of the scaffold, then write each product."""

    scaffold_smiles = Path(scaffold_path).read_text().strip()
    lib = FragmentLibrary(frag_lib_path)
    results = scaffold_with_all_fragments(scaffold_smiles, lib)

    out_root = Path(out_dir)
    out_root.mkdir(exist_ok=True)

    # --- deduplicate across all (fragment, site) pairs --------------------
    # Different fragments can accidentally produce the same SMILES.
    # We keep only the first occurrence of each unique SMILES.
    seen: dict[str, str] = {}  # smiles -> source fragment name
    for frag_name, product_mol, product_smiles in results:
        if product_smiles in seen:
            continue
        seen[product_smiles] = frag_name

    # --- build unique folder names ----------------------------------------
    smilies_list = list(seen.keys())  # preserve insertion order
    base_names = [seen[s] for s in smilies_list]

    # Assign a clean integer suffix (0-based) to each unique SMILES,
    # regardless of which fragment produced it.
    name_counts: dict[str, int] = {}
    folders: list[Path] = []
    for base_name, smiles in zip(base_names, smilies_list):
        if base_name not in name_counts:
            name_counts[base_name] = 0
        else:
            name_counts[base_name] += 1

        folder = out_root / f"{base_name}_{name_counts[base_name]}"
        folders.append(folder)

    # --- write products ---------------------------------------------------
    for smiles, folder in zip(smilies_list, folders):
        folder.mkdir(exist_ok=True)
        (folder / "product.smi").write_text(smiles + "\n")

    for folder in sorted(folders):
        smiles = Path(folder / "product.smi").read_text().strip()
        print(f"  {folder.name:20s} -> product.smi ({smiles})")

    print(f"\n{len(folders)} unique product(s) written to {out_root}/")


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        print(f"Usage: {sys.argv[0]} <scaffold_file> <fragment_file> [output_dir]")
        sys.exit(1)

    scaffold_path = sys.argv[1]
    frag_lib_path = sys.argv[2]
    out_dir = sys.argv[3] if len(sys.argv) == 4 else "products"

    scaffold_and_save(scaffold_path, frag_lib_path, out_dir)
