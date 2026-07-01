#!/usr/bin/env python3
"""Read a PDBQT file and print an AutoDock box centered on the ligand."""

import argparse


def read_pdbqt(filename: str):
    """Return (x, y, z) lists of all ATOM coordinates."""
    xs, ys, zs = [], [], []
    with open(filename) as f:
        for line in f:
            if line.startswith("ATOM"):
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
    return xs, ys, zs


def main():
    parser = argparse.ArgumentParser(description="Compute ligand bounding box from a PDBQT file")
    parser.add_argument("pdbqt", help="Path to the PDBQT file")
    args = parser.parse_args()

    xs, ys, zs = read_pdbqt(args.pdbqt)
    if not xs:
        print("ERROR: no ATOM records found", flush=True)
        return 1

    center_x = (min(xs) + max(xs)) / 2
    center_y = (min(ys) + max(ys)) / 2
    center_z = (min(zs) + max(zs)) / 2
    size_x = max(xs) - min(xs) + 1.0  # pad 1 A so ligand fits snugly
    size_y = max(ys) - min(ys) + 1.0
    size_z = max(zs) - min(zs) + 1.0

    print(f"Center_X {center_x:8.3f}")
    print(f"Center_Y {center_y:8.3f}")
    print(f"Center_Z {center_z:8.3f}")
    print(f"Size_X   {size_x:8.3f}")
    print(f"Size_Y   {size_y:8.3f}")
    print(f"Size_Z   {size_z:8.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
