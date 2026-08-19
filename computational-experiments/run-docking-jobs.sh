#!/bin/bash
# Run docking jobs for a list of (ligand x sphere-pocket) pair folders.
# Usage: run-docking-jobs.sh <folder1> [<folder2> ...]

if [ $# -lt 1 ]; then
    echo "ERROR: usage: $0 <folder1> [<folder2> ...]" >&2
    exit 1
fi

SCRIPT_DIR="$(dirname "$0")"

for DIR in "$@"; do
    echo "=== Docking: $DIR ==="
    bash "$SCRIPT_DIR/run-docking-job.sh" "$DIR"
done
