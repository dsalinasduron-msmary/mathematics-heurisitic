#!/bin/bash
# Run a docking job using files from a (ligand x sphere-pocket) pair folder.
# Usage: run-docking-job.sh <pair_folder>

if [ $# -lt 1 ]; then
    echo "ERROR: usage: $0 <pair_folder>" >&2
    exit 1
fi
DIR="$1"

vina --receptor "$DIR/1iep_receptor.pdbqt" \
     --ligand "$DIR/ligand_H.pdbqt" \
     --config "$DIR/1iep_receptor.box.txt" \
     --exhaustiveness=32 \
     --out "$DIR/1iep_ligand_vina_out.pdbqt"
