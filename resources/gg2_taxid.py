#!/usr/bin/env python3
"""
Map Greengenes2 taxonomy to NCBI TaxIDs.
Generates a tab-delimited file with sequence_id and mapped taxid.

Usage:
    python3 gg2_taxid_map.py taxonomy.tsv names.dmp tax_id.txt
"""

import sys

if len(sys.argv) != 4:
    print("Usage: python3 gg2_to_taxid_map.py <taxonomy.tsv> <names.dmp> <output.tsv>")
    sys.exit(1)

taxonomy_file = sys.argv[1]
names_dmp = sys.argv[2]
output_file = sys.argv[3]

# 1. Load NCBI taxonomy names
print("[INFO] Loading NCBI taxonomy names...")
taxid_map = {}
with open(names_dmp) as f:
    for line in f:
        parts = [p.strip() for p in line.split("|")]
        taxid, name_txt, _, name_class = parts[0], parts[1], parts[2], parts[3]
        if name_class == "scientific name":
            taxid_map[name_txt] = taxid
print(f"[INFO] Loaded {len(taxid_map)} scientific names from NCBI")

# 2. Process Greengenes2 taxonomy file
print("[INFO] Mapping Greengenes2 taxonomy...")
out = open(output_file, "w")
not_mapped = 0
mapped = 0

with open(taxonomy_file) as tf:
    for line in tf:
        seqid, taxonomy = line.strip().split("\t")
        # remove prefixes and underscores
        ranks = [t.split("__")[1].replace("_", " ") for t in taxonomy.split(";") if "__" in t]

        # try match from species -> genus -> higher ranks
        taxid = None
        for name in reversed(ranks):
            if name in taxid_map:
                taxid = taxid_map[name]
                break
        if taxid is None:
            taxid = "1"  # root if no match
            not_mapped += 1
        else:
            mapped += 1

        out.write(f"{seqid}\t{taxid}\n")

out.close()

print(f"[INFO] Done. Mapped: {mapped}, Unmapped: {not_mapped}")
print(f"[INFO] Output written to {output_file}")
