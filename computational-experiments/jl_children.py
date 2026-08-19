#!/usr/bin/env python3
"""Print SMILES of all ligands built on a given scaffold.

Usage: jl_children.py <scaffold_smiles> [file.jsonl]

Reads from stdin if piped, otherwise reads the optional file argument,
falling back to tree.jsonl.

Examples:
  jl_children.py 'c1ccccc1'
  jl_children.py 'Cc1ccccc1' docking-tree/playout_stage9/tree.jsonl
  jl_children.py 'c1ccccc1' | jl_children.py 'Cc1ccccc1'
"""
import json
import os
import stat
import sys

if len(sys.argv) < 2:
    sys.exit("Usage: jl_children.py <scaffold_smiles> [file.jsonl]")

scaffold = sys.argv[1]

stdin_mode = os.fstat(sys.stdin.fileno()).st_mode
stdin_is_pipe = stat.S_ISFIFO(stdin_mode) or stat.S_ISREG(stdin_mode)

if stdin_is_pipe:
    source = sys.stdin
elif len(sys.argv) >= 3:
    source = open(sys.argv[2])
else:
    source = open("tree.jsonl")

with source:
    for line in source:
        line = line.strip()
        if not line:
            continue
        record = json.loads(line)
        if record.get("scaffold") == scaffold:
            print(record["SMILES"])
