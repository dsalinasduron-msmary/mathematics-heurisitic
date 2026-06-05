#!/usr/bin/env python3
"""Extract sphere center and radius for every pocket from a Pocketeer output JSON."""

import json
import sys


def main():
    if len(sys.argv) < 2:
        print("Usage: pocket_spheres.py <pocketer.json>", file=sys.stderr)
        sys.exit(1)

    with open(sys.argv[1]) as f:
        data = json.load(f)

    pockets = data["object_data"]["pockets"]

    for i, pocket in enumerate(pockets):
        centers = pocket["sphere_centers"]
        radii = pocket["sphere_radii"]
        print(f"# Pocket {i} ({len(centers)} spheres)")
        for j, (c, r) in enumerate(zip(centers, radii)):
            print(f"{i}\t{j}\t{c[0]:.4f}\t{c[1]:.4f}\t{c[2]:.4f}\t{r:.4f}")


if __name__ == "__main__":
    main()
