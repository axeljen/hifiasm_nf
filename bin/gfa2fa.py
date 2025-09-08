#!/usr/bin/env python3

import sys

gfa = sys.argv[1]

out = sys.argv[2]

def parse_line(line):
    if line.startswith("S"):
        parts = line.strip().split("\t")
        seq_id = parts[1]
        sequence = parts[2]
        # wrap sequence every 60th character
        wrapped_sequence = "\n".join([sequence[i:i+60] for i in range(0, len(sequence), 60)])
        return f">{seq_id}\n{wrapped_sequence}\n"
    return ""

with open(gfa, "r") as gfa_file, open(out, "w") as fasta_file:
    for line in gfa_file:
        fasta_line = parse_line(line)
        fasta_file.write(fasta_line)
