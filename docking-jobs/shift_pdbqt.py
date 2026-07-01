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


def main(pdbqt_string: str, cx: float = 0.0, cy: float = 0.0, cz: float = 0.0):
    """Shift all atoms in a PDBQT string so their centroid moves to (cx, cy, cz).

    Args:
        pdbqt_string: PDBQT contents as a string.
        cx, cy, cz: Target center coordinates (default origin).

    Returns:
        The shifted PDBQT as a string.
    """
    # Parse all ATOM coordinates and the original lines
    xs, ys, zs, lines = [], [], [], pdbqt_string.splitlines(keepends=True)
    for line in lines:
        if line.startswith("ATOM"):
            xs.append(float(line[30:38]))
            ys.append(float(line[38:46]))
            zs.append(float(line[46:54]))

    if not xs:
        raise ValueError("no ATOM records found")

    # Compute centroid (average of all atom coordinates)
    center_x = sum(xs) / len(xs)
    center_y = sum(ys) / len(ys)
    center_z = sum(zs) / len(zs)

    dx = cx - center_x
    dy = cy - center_y
    dz = cz - center_z

    out = []
    for line in lines:
        if line.startswith("ATOM"):
            x = float(line[30:38]) + dx
            y = float(line[38:46]) + dy
            z = float(line[46:54]) + dz
            out.append(f"{line[:30]}{x:>8.3f}{y:>8.3f}{z:>8.3f}{line[54:]}")
        else:
            out.append(line)
    return "".join(out)


if __name__ == "__main__":
    # Usage: python shift_pdbqt.py
    # Called programmatically:  result = main(pdbqt_string, cx, cy, cz)
    pass
