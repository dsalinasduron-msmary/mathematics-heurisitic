#!/usr/bin/env python3
"""Compare INTER docking forces between scaffolds and their fragment-grown children.

Usage: compare_inter.py <docking_folder> [tree.jsonl]

Scans <docking_folder> for ligX_pocketY_sphereZ subdirectories containing
1iep_ligand_vina_out.pdbqt files, then for every scaffold→child pair in
tree.jsonl emits one row per pocket-sphere where both have docking results.

Output columns (tab-separated, no header):
  pocket  sphere  scaffold_smiles  inter_scaffold  child_smiles  inter_child
"""
import os
import re
import sys


TARGET = "1iep_ligand_vina_out.pdbqt"


def parse_inter(path):
    """Return INTER value from MODEL 1 of a Vina PDBQT, or None."""
    in_model1 = False
    with open(path) as f:
        for line in f:
            if line.startswith("MODEL 1"):
                in_model1 = True
                continue
            if in_model1:
                if line.startswith("ENDMDL") or line.startswith("MODEL "):
                    break
                if re.match(r"REMARK INTER:\s", line):
                    return line.split()[2]
    return None


if len(sys.argv) < 2:
    sys.exit("Usage: compare_inter.py <docking_folder> [tree.jsonl]")

docking_root = sys.argv[1]
jsonl_path   = sys.argv[2] if len(sys.argv) >= 3 else os.path.join(docking_root, "tree.jsonl")

# -- load tree.jsonl ----------------------------------------------------------─
smiles_to_idx      = {}   # smiles -> line index (= ligX number)
child_to_scaffold  = {}   # child_smiles -> scaffold_smiles

with open(jsonl_path) as f:
    import json
    for idx, line in enumerate(f):
        record = json.loads(line.strip())
        smiles   = record["SMILES"]
        scaffold = record.get("scaffold", "")
        smiles_to_idx[smiles] = idx
        child_to_scaffold[smiles] = scaffold

# -- scan docking folders ------------------------------------------------------
# inter_map[(lig_idx, pocket, sphere)] = inter_value
inter_map = {}
folder_re = re.compile(r"lig(\d+)_(pocket\d+)_(sphere\d+)$")

for entry in os.scandir(docking_root):
    m = folder_re.match(entry.name)
    if not m:
        continue
    pdbqt = os.path.join(entry.path, TARGET)
    if not os.path.exists(pdbqt):
        continue
    inter = parse_inter(pdbqt)
    if inter is not None:
        key = (int(m.group(1)), m.group(2), m.group(3))
        inter_map[key] = inter

# -- emit rows ----------------------------------------------------------------─
for child_smiles, scaffold_smiles in child_to_scaffold.items():
    if child_smiles == scaffold_smiles:
        continue                      # skip root scaffold pointing to itself
    scaffold_idx = smiles_to_idx.get(scaffold_smiles)
    child_idx    = smiles_to_idx.get(child_smiles)
    if scaffold_idx is None or child_idx is None:
        continue
    # find all pocket-sphere combos where both have results
    scaffold_keys = {(p, s): v for (i, p, s), v in inter_map.items() if i == scaffold_idx}
    child_keys    = {(p, s): v for (i, p, s), v in inter_map.items() if i == child_idx}
    for (pocket, sphere), inter_scaffold in scaffold_keys.items():
        inter_child = child_keys.get((pocket, sphere))
        if inter_child is not None:
            print("\t".join([pocket, sphere, scaffold_smiles,
                             inter_scaffold, child_smiles, inter_child]))
