"""Render 2D molecular structures for each entry in an input SMILES file.

Usage:
    python draw_rings.py <input.smi> [output_dir]
"""

import argparse
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

BASE = Path(__file__).resolve().parent

W, H = 300, 250


def _validate_input(path: Path) -> None:
    if not path.exists():
        print(f"Error: input file not found: {path}", file=sys.stderr)
        sys.exit(1)
    if not path.is_file():
        print(f"Error: input must be a regular file: {path}", file=sys.stderr)
        sys.exit(1)


def main(input_file: str = "rings_6-aromatic.smi", outdir_name: str = "figures"):
    INPUT = Path(input_file) if Path(input_file).is_absolute() else BASE / input_file
    _validate_input(INPUT)

    OUTDIR = BASE / outdir_name

    lines = [l.strip() for l in INPUT.read_text().splitlines() if l.strip()]
    parts = [r.rsplit(",", 1) for r in lines]
    names, smileses = zip(*parts)

    OUTDIR.mkdir(exist_ok=True)

    for name, smiles in zip(names, smileses):
        mol = _make_mol(smiles)
        svg = _render_svg(mol, legend=name)
        fname = OUTDIR / f"{name}.svg"
        fname.write_text(svg)
        print(f"Wrote {fname}")

    # Also make a single combined HTML file with all structures inline
    styles = ("body{font-family:sans-serif;margin:2em;} "
              "div{display:inline-block;margin:1.5em;vertical-align:top;text-align:center}"
              "p{margin:6px 0;font-size:14px;font-weight:bold}")
    html_lines = [f'<html><head><style>{styles}</style></head><body>']
    for name, smiles in zip(names, smileses):
        mol = _make_mol(smiles)
        svg = _render_svg(mol, legend="")
        html_lines.append(f'<p>{name}</p>')
        html_lines.append(svg)
    html_lines.append("</body></html>")

    fname = OUTDIR / "all_molecules.html"
    fname.write_text("\n".join(html_lines))
    print(f"Wrote {fname} (combined view)")


def _make_mol(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")
    rdDepictor.Compute2DCoords(mol)
    return mol


def _render_svg(mol: Chem.Mol, legend: str = "") -> str:
    drawer = rdMolDraw2D.MolDraw2DSVG(W, H)
    drawer.SetLineWidth(1.5)
    drawer.DrawMolecule(mol, legend=legend)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Render 2D molecular structures from a SMILES file.",
    )
    parser.add_argument(
        "input",
        help="Path to the input SMILES file (default: rings_6-aromatic.smi)",
        default="rings_6-aromatic.smi",
    )
    parser.add_argument(
        "outdir",
        nargs="?",
        help="Output directory name (default: figures)",
        default="figures",
    )
    args = parser.parse_args()
    main(input_file=args.input, outdir_name=args.outdir)
