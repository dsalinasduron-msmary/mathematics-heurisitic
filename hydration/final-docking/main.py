#!/usr/bin/env python3
import json
import os
import subprocess

from ligand_sdf import write_sdf

SCRIPT = "./prepare-ligand_hydration.sh"

with open("tree.jsonl") as f:
    entries = [json.loads(line) for line in f]

with open("pocket-spheres.tsv") as f:
    header = next(f).strip()  # "pocket\tsphere\tx\ty\tz\tr"
    for pair_idx, line in enumerate(f):
        parts = line.strip().split("\t")
        pocket, sphere = parts[0], parts[1]
        for ligand_idx in range(len(entries)):
            folder_name = f"lig{ligand_idx}_pocket{pocket}_sphere{sphere}"
            os.makedirs(folder_name, exist_ok=True)
            # Write the full entry JSONL
            filepath = os.path.join(folder_name, "entry.jsonl")
            with open(filepath, "w") as out:
                out.write(json.dumps(entries[ligand_idx]) + "\n")
            # Write the pocket-sphere row as TSV with headers
            sphpath = os.path.join(folder_name, "sphere.tsv")
            with open(sphpath, "w") as out:
                out.write(header + "\n")
                out.write(line.strip() + "\n")
            # Write the ligand SDF file to this pair's folder
            sdf_path = os.path.join(folder_name, "ligand.sdf")
            write_sdf(entries[ligand_idx]["SMILES"], sdf_path)
            # Protonate and add hydration shell via prepare-ligand_hydration.sh
            subprocess.run(
                [SCRIPT, sdf_path],
                check=True,
            )

            # Generate docking box parameters from the hydrated ligand PDBQT
            pdbqt_path = os.path.join(folder_name, "ligand_H.pdbqt")
            result = subprocess.run(
                ["python3", "get_dock_params.py", "-p", pocket, "-s", sphere, pdbqt_path],
                capture_output=True, text=True, check=True,
            )
            params_path = os.path.join(folder_name, "dock_params.txt")
            with open(params_path, "w") as out:
                out.write(result.stdout)
