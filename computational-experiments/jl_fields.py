#!/usr/bin/env python3
"""Extract fields from JSONL and print as TSV.

Usage: jl_fields.py FIELD [FIELD ...]

Reads from stdin if piped, otherwise reads tree.jsonl.

Examples:
  jl_fields.py SMILES scaffold
  jl_fields.py SMILES | sort -u
  cat other.jsonl | jl_fields.py functional_group site
"""
import json
import os
import stat
import sys

if len(sys.argv) < 2:
    sys.exit("Usage: jl_fields.py FIELD [FIELD ...]")

fields = sys.argv[1:]

stdin_mode = os.fstat(sys.stdin.fileno()).st_mode
stdin_is_pipe = stat.S_ISFIFO(stdin_mode) or stat.S_ISREG(stdin_mode)
source = sys.stdin if stdin_is_pipe else open("tree.jsonl")

with source:
    for line in source:
        line = line.strip()
        if not line:
            continue
        record = json.loads(line)
        print("\t".join(str(record.get(f, "")) for f in fields))
