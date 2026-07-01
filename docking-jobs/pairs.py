#!/usr/bin/env python3
"""Cartesian product of every ligand × every pocket sphere; writes shifted PDBQT files."""

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

for i, lig in enumerate(ligands):
    _, pdbqt_string = prepare_for_docking(lig["SMILES"])
    for j, (pot, sph, x, y, z, r) in enumerate(spheres):
        outfile = f"pocket-{pot}_ligand-{i}_sphere-{j}.pdbqt"
        print(f"--- SMILES {lig['SMILES']}\t pocket-{pot} sphere-{sph} ({x:.4f}, {y:.4f}, {z:.4f}) radius={r} → {outfile}")
        shifted = shift_main(pdbqt_string, x, y, z)
        with open(outfile, "w") as f:
            f.write(shifted)
