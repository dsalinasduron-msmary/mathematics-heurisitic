#!/usr/bin/env python3
"""Search a folder tree for Vina output PDBQT files and print a TSV of results.

Usage: parse_vina.py <parent_folder>

Output columns (tab-separated):
  smiles    affinity    inter    pocket    sphere
Only MODEL 1 is parsed from each file.
"""
import os
import re
import sys

TARGET = "1iep_ligand_vina_out.pdbqt"


def parse_model1(path):
    """Return (smiles, affinity, inter) from MODEL 1 of a Vina PDBQT, or None."""
    smiles = affinity = inter = None
    in_model1 = False
    with open(path) as f:
        for line in f:
            if line.startswith("MODEL 1"):
                in_model1 = True
                continue
            if in_model1:
                if line.startswith("ENDMDL") or line.startswith("MODEL "):
                    break
                if line.startswith("REMARK SMILES ") and not line.startswith("REMARK SMILES IDX"):
                    smiles = line.split(None, 2)[2].strip()
                elif line.startswith("REMARK VINA RESULT:"):
                    affinity = line.split()[3]
                elif re.match(r"REMARK INTER:\s", line):
                    inter = line.split()[2]
    if smiles and affinity and inter:
        return smiles, affinity, inter
    return None


if len(sys.argv) != 2:
    sys.exit("Usage: parse_vina.py <parent_folder>")

root = sys.argv[1]
if not os.path.isdir(root):
    sys.exit(f"Not a directory: {root}")

for dirpath, _, filenames in os.walk(root):
    if TARGET in filenames:
        result = parse_model1(os.path.join(dirpath, TARGET))
        if result:
            folder = os.path.basename(os.path.normpath(dirpath))
            m = re.search(r"(pocket\d+)_(sphere\d+)", folder)
            pocket = m.group(1) if m else ""
            sphere = m.group(2) if m else ""
            print("\t".join((*result, pocket, sphere)))
