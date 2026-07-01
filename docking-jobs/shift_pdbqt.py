#!/usr/bin/env python3
"""Shift all atoms in a PDBQT file so their center moves to a given centroid."""

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
    parser = argparse.ArgumentParser(description="Shift all atoms in a PDBQT file to a given centroid")
    parser.add_argument("pdbqt", help="Path to the PDBQT file")
    parser.add_argument("centroid", help="New center (x,y,z)")
    args = parser.parse_args()

    cx, cy, cz = (float(x) for x in args.centroid.split(","))
    xs, ys, zs = read_pdbqt(args.pdbqt)
    if not xs:
        print("ERROR: no ATOM records found", flush=True)
        return 1

    center_x = (min(xs) + max(xs)) / 2
    center_y = (min(ys) + max(ys)) / 2
    center_z = (min(zs) + max(zs)) / 2
    dx = cx - center_x
    dy = cy - center_y
    dz = cz - center_z

    with open(args.pdbqt) as f:
        for line in f:
            if line.startswith("ATOM"):
                x = float(line[30:38]) + dx
                y = float(line[38:46]) + dy
                z = float(line[46:54]) + dz
                print(f"{line[:30]}{x:>8.3f}{y:>8.3f}{z:>8.3f}{line[54:]}")
            else:
                print(line, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
