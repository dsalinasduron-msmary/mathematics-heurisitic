#!/usr/bin/env python3
"""Cartesian product of every ligand × every pocket sphere, with PDBQT."""

import json
from ligand_pdbqt import prepare_for_docking
from shift_pdbqt import main as shift_main

with open("tree.jsonl") as f:
    ligands = [json.loads(line) for line in f]

spheres = []
with open("pocket-spheres_head.tsv") as f:
    header = f.readline()          # skip header
    for line in f:
        pocket, sphere, x, y, z, r = line.split()
        spheres.append((int(pocket), int(sphere), float(x), float(y), float(z), float(r)))

for lig in ligands:
    mol, pdbqt_string = prepare_for_docking(lig["SMILES"])
    for pot, sph, x, y, z, r in spheres:
        print(f"--- SMILES {lig['SMILES']}\t pocket-{pot} sphere-{sph} ({x:.4f}, {y:.4f}, {z:.4f}) radius={r}")
        shifted = shift_main(pdbqt_string, x, y, z)
        print(shifted)
