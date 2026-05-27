#!/usr/bin/env python3

import sys
from collections import Counter

if len(sys.argv) != 4:
    print(f"Usage: {sys.argv[0]} <sequence_file_mapping.tsv> <taxonomy.tsv> <output.txt>")
    sys.exit(1)

file1, file2, outfile = sys.argv[1], sys.argv[2], sys.argv[3]

taxonomy = {}
with open(file2) as f2:
    for line in f2:
        parts = line.strip().split(None, 1)  # split only on first whitespace
        if len(parts) == 2:
            key, tax = parts
            taxonomy[key] = tax

tax_list = []  # to collect taxonomy strings for counting
with open(file1) as f1, open(outfile+'.all', "w") as out:
    out.write("#READ\tTAXONOMY\n")  
    for line in f1:
        parts = line.strip().split()
        if len(parts) >= 4:
            key = parts[0]
            fasta_name = parts[3].removesuffix(".fasta")
            if key in taxonomy:
                tax_str = taxonomy[key]
                out.write(f"{fasta_name}\t{tax_str}\n")
                tax_list.append(tax_str)

counts = Counter(tax_list)
with open(outfile, "w") as s:
    s.write("#LEMMI16s\nGroup_ID\tSize\tTaxonomy\n")
    for i, (tax_str, count) in enumerate(counts.items(), start=1):
        s.write(f"{i}\t{count}\t{tax_str}\n")
