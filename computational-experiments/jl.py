#!/usr/bin/env python3
"""Filter and select fields from JSONL records.

Reads from stdin if data is piped in, otherwise reads tree.jsonl.

Examples:
  jl.py --eq functional_group=scaffold
  jl.py --ne site=null
  jl.py --select SMILES,functional_group
  jl.py --eq functional_group=scaffold | jl.py --select SMILES
"""
import argparse
import json
import os
import stat
import sys


def parse_value(v):
    """Convert string 'null'/'true'/'false'/numbers to Python types."""
    if v == "null":  return None
    if v == "true":  return True
    if v == "false": return False
    try:             return int(v)
    except ValueError: pass
    try:             return float(v)
    except ValueError: pass
    return v


def matches(record, eq, ne, has, lacks):
    for kv in eq:
        k, _, v = kv.partition("=")
        if record.get(k) != parse_value(v):
            return False
    for kv in ne:
        k, _, v = kv.partition("=")
        if record.get(k) == parse_value(v):
            return False
    for k in has:
        if k not in record:
            return False
    for k in lacks:
        if k in record:
            return False
    return True


parser = argparse.ArgumentParser(
    description="Filter and select fields from JSONL records.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)
parser.add_argument("--eq",     metavar="KEY=VAL", action="append", default=[],
                    help="keep records where KEY equals VAL")
parser.add_argument("--ne",     metavar="KEY=VAL", action="append", default=[],
                    help="keep records where KEY does not equal VAL")
parser.add_argument("--has",    metavar="KEY",     action="append", default=[],
                    help="keep records that contain KEY")
parser.add_argument("--lacks",  metavar="KEY",     action="append", default=[],
                    help="keep records that do not contain KEY")
parser.add_argument("--select", metavar="KEY,...",
                    help="output only these comma-separated fields")
parser.add_argument("file", nargs="?", default="tree.jsonl",
                    help="JSONL file to read (default: tree.jsonl); ignored when stdin is piped")
args = parser.parse_args()

fields = [f.strip() for f in args.select.split(",")] if args.select else None

stdin_mode = os.fstat(sys.stdin.fileno()).st_mode
stdin_is_pipe = stat.S_ISFIFO(stdin_mode) or stat.S_ISREG(stdin_mode)
source = sys.stdin if stdin_is_pipe else open(args.file)

with source:
    for line in source:
        line = line.strip()
        if not line:
            continue
        record = json.loads(line)
        if not matches(record, args.eq, args.ne, args.has, args.lacks):
            continue
        if fields:
            record = {k: record[k] for k in fields if k in record}
        print(json.dumps(record))
