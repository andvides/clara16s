#!/usr/bin/env python3
import sys, gzip
import os
import pandas as pd
import numpy as np
import csv

pd.options.mode.chained_assignment = None  # default='warn'

taxo_result=sys.argv[1]
TOTAL_READS = int (sys.argv[2])   
out1       =sys.argv[3]

rank_prefixes = [
    "d__", "p__", "c__", "o__", "f__", "g__","s__"
]

with open(taxo_result) as infile, open(out1, "w", newline="") as outfile:
    reader = csv.DictReader(infile, delimiter="\t")
    writer = csv.writer(outfile, delimiter="\t")
    writer.writerow(["#LEMMI16s"])
    writer.writerow(["Group_ID", "Size", "Taxonomy"])

    for i, row in enumerate(reader, start=1):
        otu_id = f"OTU_{i}"

        # abundance → count
        size = int(float(row["abundance"]) * TOTAL_READS)

        lineage = row["lineage"].strip(";").split(";")

        taxonomy = ";".join(
            f"{rank_prefixes[j]}{lineage[j]}" if j < len(lineage) else f"{rank_prefixes[j]}na"
            for j in range(len(rank_prefixes))
        )

        #taxonomy += ";s__na"

        writer.writerow([otu_id, size, taxonomy])


