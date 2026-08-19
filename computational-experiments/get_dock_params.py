#!/usr/bin/env python3
"""Compute the pocket-sphere centre and ligand bounding box.

Usage:
    get_dock_params.py -p <pocket> -s <sphere> [-t <pocket-spheres.tsv>] ligand.pdbqt

Reads the PDBQT file via pdbqt_box.read_pdbqt() and the pocket-sphere entry
via sphere_box (subprocess), then prints a combined report.
"""

import argparse
import subprocess
from pathlib import Path

from pdbqt_box import read_pdbqt


def main():
    parser = argparse.ArgumentParser(
        description="Compute pocket centre and ligand box from PDBQT / pocket-spheres data"
    )
    parser.add_argument("pdbqt", help="Path to the ligand PDBQT file")
    parser.add_argument("-p", "--pocket", type=int, required=True,
                        help="Pocket ID (matches the first column of pocket-spheres.tsv)")
    parser.add_argument("-s", "--sphere", type=int, required=True,
                        help="Sphere ID (matches the second column of pocket-spheres.tsv)")
    parser.add_argument("-t", "--tsv-file", type=Path, default=None,
                        help="Path to pocket-spheres.tsv (default: pocket-spheres.tsv in CWD)")
    args = parser.parse_args()

    # --- ligand box from PDBQT via pdbqt_box ---
    xs, ys, zs = read_pdbqt(args.pdbqt)
    if not xs:
        parser.error("no ATOM records found in the PDBQT file")

    lig_center_x = (min(xs) + max(xs)) / 2
    lig_center_y = (min(ys) + max(ys)) / 2
    lig_center_z = (min(zs) + max(zs)) / 2
    lig_size = max(max(xs) - min(xs),
                   max(ys) - min(ys),
                   max(zs) - min(zs)) + 1.0
    lig_size_x = lig_size_y = lig_size_z = lig_size

    # --- pocket sphere centre via subprocess (sphere_box has no importable API) ---
    tsv_arg = [str(args.tsv_file)] if args.tsv_file else []
    result = subprocess.run(
        ["python3", str(Path(__file__).with_name("sphere_box.py")),
         str(args.pocket), str(args.sphere), *tsv_arg],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        parser.error(result.stderr.strip())

    parts = result.stdout.strip().split()       # cx cy cz sx sy sz (tab-separated)
    sphere_cx = float(parts[0])
    sphere_cy = float(parts[1])
    sphere_cz = float(parts[2])

    # --- report (matching example.txt format) ---
    print(f"center_x = {sphere_cx:.3f}")
    print(f"center_y = {sphere_cy:.3f}")
    print(f"center_z = {sphere_cz:.3f}")
    print(f"size_x = {lig_size_x:.1f}")
    print(f"size_y = {lig_size_y:.1f}")
    print(f"size_z = {lig_size_z:.1f}")


if __name__ == "__main__":
    main()
