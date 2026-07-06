#!/usr/bin/env python3
"""Print a bounding box for a single pocket-sphere entry.

Usage: sphere-box.py <pocket> <sphere> [tsv_file]

Reads pocket-spheres.tsv and outputs six values (tab-separated):
  center_x center_y center_z side_x side_y side_z

The box is the smallest cube enclosing the selected sphere, so
the sides all equal 2*r and the centre is the sphere centre."""

import sys
from pathlib import Path


def main():
    if len(sys.argv) < 3:
        print("Usage: sphere-box.py <pocket> <sphere> [tsv_file]", file=sys.stderr)
        sys.exit(1)

    pocket = int(sys.argv[1])
    sphere = int(sys.argv[2])
    tsv_path = Path(sys.argv[3]) if len(sys.argv) > 3 else Path(__file__).with_name("pocket-spheres.tsv")

    with open(tsv_path) as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue  # skip header
            parts = line.strip().split("\t")
            if int(parts[0]) == pocket and int(parts[1]) == sphere:
                x, y, z, r = (float(v) for v in parts[2:])
                side = 2 * r
                print(f"{x:.4f}\t{y:.4f}\t{z:.4f}\t{side:.4f}\t{side:.4f}\t{side:.4f}")
                return

    print(f"Error: no entry with pocket={pocket}, sphere={sphere} in {tsv_path}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
