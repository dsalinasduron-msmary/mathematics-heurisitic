#!/bin/bash
# Prepare receptor box info by reading center and size from a boxinfo file.
# Usage: prepare-receptor_boxinfo.sh <boxinfo.txt> <output_dir>

if [ $# -lt 2 ] || [ ! -f "$1" ]; then
    echo "ERROR: usage: $0 <boxinfo.txt> <output_dir>" >&2
    exit 1
fi
INFO_FILE="$1"
OUT_DIR="$2"

# Read values from the provided boxinfo file (Python key=value format)
center_x=$(grep '^center_x' "$INFO_FILE" | cut -d= -f2 | tr -d ' ')
center_y=$(grep '^center_y' "$INFO_FILE" | cut -d= -f2 | tr -d ' ')
center_z=$(grep '^center_z' "$INFO_FILE" | cut -d= -f2 | tr -d ' ')
size_x=$(grep '^size_x' "$INFO_FILE" | cut -d= -f2 | tr -d ' ')
size_y=$(grep '^size_y' "$INFO_FILE" | cut -d= -f2 | tr -d ' ')
size_z=$(grep '^size_z' "$INFO_FILE" | cut -d= -f2 | tr -d ' ')

# Validate that all values were read
for var in center_x center_y center_z size_x size_y size_z; do
    if [ -z "${!var}" ]; then
        echo "ERROR: $var not found or empty in $INFO_FILE" >&2
        exit 1
    fi
done

# Use the original PDB input (unchanged) and pass box info from file
mk_prepare_receptor.py -i 1iep_receptorH.pdb -o "$OUT_DIR/1iep_receptor" -p -v \
	--box_size "$size_x" "$size_y" "$size_z" \
	--box_center "$center_x" "$center_y" "$center_z"
