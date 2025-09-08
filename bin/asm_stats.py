#! /usr/bin/env python

import sys
# from Bio import bgzf

asm = sys.argv[1]

def open_asm(asm):
    return open(asm)
    #if asm.endswith("gz"):
     #   return bgzf.open(asm)
    #else:
     #   return open(asm)

def parse_asm(asm):
    scaffold_lengths = []
    contig_lengths = []
    current_scaffold = ""
    with open_asm(asm) as f:
        for line in f:
            if line.startswith(">"):
                if current_scaffold != "":
                    scaffold_lengths.append(len(current_scaffold))
                    contigs = [seq for seq in current_scaffold.upper().split("N") if len(seq) > 0 and not seq == "N"]
                    for contig in contigs:
                        contig_lengths.append(len(contig))
                current_scaffold = ""
            else:
                current_scaffold += line.strip()
    return scaffold_lengths, contig_lengths


def get_lengths(asm):
    return sorted([len(seq) for seq in asm.values()], reverse=True)

def get_totlen(lengths):
    return sum(lengths)

def get_n50(lengths):
    total = get_totlen(lengths)
    half = total / 2
    running_total = 0
    for length in lengths:
        running_total += length
        if running_total >= half:
            return length
    return 0

def get_l50(lengths):
    total = get_totlen(lengths)
    half = total / 2
    running_total = 0
    for i, length in enumerate(lengths):
        running_total += length
        if running_total >= half:
            return i + 1
    return 0

def get_count(lengths):
    return len(lengths)

def get_longest(lengths):
    return max(lengths) if lengths else 0

# prep a dict for holding the results
results = {'contigs': {'N50': 0, 'L50': 0, 'count': 0, 'longest': 0, 'totlen': 0},
           'scaffolds': {'N50': 0, 'L50': 0, 'count': 0, 'longest': 0, 'totlen': 0}}

scaffolds,contigs = parse_asm(asm)

results['contigs']['N50'] = get_n50(contigs)
results['contigs']['L50'] = get_l50(contigs)
results['contigs']['count'] = get_count(contigs)
results['contigs']['longest'] = get_longest(contigs)
results['contigs']['totlen'] = get_totlen(contigs)

results['scaffolds']['N50'] = get_n50(scaffolds)
results['scaffolds']['L50'] = get_l50(scaffolds)
results['scaffolds']['count'] = get_count(scaffolds)
results['scaffolds']['longest'] = get_longest(scaffolds)
results['scaffolds']['totlen'] = get_totlen(scaffolds)

print("level\ttotal_length\tlongest_sequence\tnum_sequence\tN50\tL50")
for level in ['contigs', 'scaffolds']:
    print(f"{level}\t{results[level]['totlen']}\t{results[level]['longest']}\t{results[level]['count']}\t{results[level]['N50']}\t{results[level]['L50']}")