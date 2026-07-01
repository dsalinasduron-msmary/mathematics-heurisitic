#!/usr/bin/env python3
"""Extract pocket information from a Pocketeer output JSON."""

import json
import sys


def main():
    if len(sys.argv) < 2:
        print("Usage: pockets.py <pocketer.json>", file=sys.stderr)
        sys.exit(1)

    with open(sys.argv[1]) as f:
        data = json.load(f)

    pockets = data["object_data"]["pockets"]

    print("pocket\tcenter_x\tcenter_y\tcenter_z\tside_x\tside_y\tside_z")
    for i, pocket in enumerate(pockets):
        center = pocket["pocket_center"]
        sides = pocket["pocket_sides"]
        print(f"{i}\t{center[0]:.4f}\t{center[1]:.4f}\t{center[2]:.4f}\t{sides[0]:.4f}\t{sides[1]:.4f}\t{sides[2]:.4f}")


if __name__ == "__main__":
    main()
