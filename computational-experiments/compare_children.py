#!/usr/bin/env python3
"""Compare docking results between scaffolds and their fragment-grown children.

Usage: compare_children.py <docking_folder> [tree.jsonl]

Scans <docking_folder> for ligX_pocketY_sphereZ subdirectories containing
1iep_ligand_vina_out.pdbqt files, then for every scaffold->child pair in
tree.jsonl emits one row per pocket-sphere where both have docking results.

Output columns (tab-separated, no header):
  pocket  sphere
  scaffold_smiles  scaffold_affinity  scaffold_inter  scaffold_intra
  child_smiles     child_affinity     child_inter     child_intra
"""
import json
import os
import re
import sys


TARGET = "1iep_ligand_vina_out.pdbqt"


def parse_model1(path):
    """Return (affinity, inter, intra) from MODEL 1 of a Vina PDBQT, or None."""
    affinity = inter = intra = None
    in_model1 = False
    with open(path) as f:
        for line in f:
            if line.startswith("MODEL 1"):
                in_model1 = True
                continue
            if in_model1:
                if line.startswith("ENDMDL") or line.startswith("MODEL "):
                    break
                if line.startswith("REMARK VINA RESULT:"):
                    affinity = line.split()[3]
                elif re.match(r"REMARK INTER:\s", line):
                    inter = line.split()[2]
                elif line.startswith("REMARK INTRA:"):
                    intra = line.split()[2]
    if affinity and inter and intra:
        return affinity, inter, intra
    return None


if len(sys.argv) < 2:
    sys.exit("Usage: compare_children.py <docking_folder> [tree.jsonl]")

docking_root = sys.argv[1]
jsonl_path   = sys.argv[2] if len(sys.argv) >= 3 else os.path.join(docking_root, "tree.jsonl")

# -- load tree.jsonl -----------------------------------------------------------
smiles_to_idx     = {}   # smiles -> line index (= ligX number)
child_to_scaffold = {}   # child_smiles -> scaffold_smiles

with open(jsonl_path) as f:
    for idx, line in enumerate(f):
        record = json.loads(line.strip())
        smiles   = record["SMILES"]
        scaffold = record.get("scaffold", "")
        smiles_to_idx[smiles] = idx
        child_to_scaffold[smiles] = scaffold

# -- scan docking folders ------------------------------------------------------
# results_map[(lig_idx, pocket, sphere)] = (affinity, inter, intra)
results_map = {}
folder_re = re.compile(r"lig(\d+)_(pocket\d+)_(sphere\d+)$")

for entry in os.scandir(docking_root):
    m = folder_re.match(entry.name)
    if not m:
        continue
    pdbqt = os.path.join(entry.path, TARGET)
    if not os.path.exists(pdbqt):
        continue
    result = parse_model1(pdbqt)
    if result is not None:
        key = (int(m.group(1)), m.group(2), m.group(3))
        results_map[key] = result

# -- emit rows -----------------------------------------------------------------
for child_smiles, scaffold_smiles in child_to_scaffold.items():
    if child_smiles == scaffold_smiles:
        continue                      # skip root scaffold pointing to itself
    scaffold_idx = smiles_to_idx.get(scaffold_smiles)
    child_idx    = smiles_to_idx.get(child_smiles)
    if scaffold_idx is None or child_idx is None:
        continue
    scaffold_keys = {(p, s): v for (i, p, s), v in results_map.items() if i == scaffold_idx}
    child_keys    = {(p, s): v for (i, p, s), v in results_map.items() if i == child_idx}
    for (pocket, sphere), (s_aff, s_inter, s_intra) in scaffold_keys.items():
        child_result = child_keys.get((pocket, sphere))
        if child_result is not None:
            c_aff, c_inter, c_intra = child_result
            print("\t".join([pocket, sphere,
                             scaffold_smiles, s_aff, s_inter, s_intra,
                             child_smiles,    c_aff, c_inter, c_intra]))
