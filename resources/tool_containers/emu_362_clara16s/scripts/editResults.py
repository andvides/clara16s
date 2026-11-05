#!/usr/bin/env python3

import sys
import csv
import re
from decimal import Decimal, InvalidOperation

# --------------------------------------------------
# Arguments
# --------------------------------------------------

if len(sys.argv) != 4:
    sys.exit(
        "Usage: editResults.py <emu_results.tsv> <total_reads> <output.tsv>"
    )

taxo_result = sys.argv[1]
total_reads = int(sys.argv[2])
output_file = sys.argv[3]

# --------------------------------------------------
# Taxonomic ranks
# --------------------------------------------------

rank_prefixes = [
    "d__",
    "p__",
    "c__",
    "o__",
    "f__",
    "g__",
    "s__"
]

# --------------------------------------------------
# Process Emu results
# --------------------------------------------------

otu_counter = 0
total_estimated_counts = Decimal("0")
total_abundance = Decimal("0")

with open(
    taxo_result, "r", newline="", encoding="utf-8-sig"
) as infile:

    reader = csv.DictReader(infile, delimiter="\t")

    if not reader.fieldnames:
        sys.exit("ERROR: Empty or invalid Emu results file")

    # Normalize column names:
    # "estimated counts" -> "estimated_counts"
    column_mapping = {
        re.sub(r"\s+", "_", name.strip().lower()): name
        for name in reader.fieldnames
    }

    required_columns = {
        "tax_id",
        "abundance",
        "lineage",
        "estimated_counts"
    }

    missing_columns = required_columns - set(column_mapping)

    if missing_columns:
        sys.exit(
            "ERROR: Missing columns in Emu results: "
            + ", ".join(sorted(missing_columns))
        )

    tax_id_col = column_mapping["tax_id"]
    abundance_col = column_mapping["abundance"]
    lineage_col = column_mapping["lineage"]
    counts_col = column_mapping["estimated_counts"]

    with open(
        output_file, "w", newline="", encoding="utf-8"
    ) as outfile:

        writer = csv.writer(outfile, delimiter="\t")

        # LEMMI16s output header
        writer.writerow(["#LEMMI16s"])
        writer.writerow(["Group_ID", "Size", "Taxonomy"])

        for row in reader:

            tax_id = (row[tax_id_col] or "").strip()

            # Exclude unclassified reads
            if tax_id.lower() == "unclassified":
                continue

            abundance_value = (row[abundance_col] or "").strip()
            counts_value = (row[counts_col] or "").strip()
            lineage_value = (row[lineage_col] or "").strip()

            if not abundance_value or not counts_value:
                continue

            try:
                abundance = Decimal(abundance_value)
                estimated_counts = Decimal(counts_value)

            except InvalidOperation:
                sys.exit(
                    f"ERROR: Invalid numeric value for tax_id {tax_id}"
                )

            if abundance <= 0 or estimated_counts <= 0:
                continue

            if not lineage_value:
                continue

            # Parse taxonomic lineage
            lineage = [
                value.strip()
                for value in lineage_value.strip(";").split(";")
            ]

            # Build standardized LEMMI16s taxonomy
            taxonomy = []

            for i, prefix in enumerate(rank_prefixes):

                value = lineage[i] if i < len(lineage) else ""

                if not value or value.lower() in ("nan", "none"):
                    value = "na"

                taxonomy.append(f"{prefix}{value}")

            taxonomy_string = ";".join(taxonomy)

            # Use Emu estimated counts directly
            size = format(estimated_counts, "f")

            otu_counter += 1
            total_estimated_counts += estimated_counts
            total_abundance += abundance

            writer.writerow([
                f"OTU_{otu_counter}",
                size,
                taxonomy_string
            ])

# --------------------------------------------------
# Summary
# --------------------------------------------------

print("----------------------------------------")
print("LEMMI16s results generated successfully")
print("----------------------------------------")
print(f"Input reads             : {total_reads}")
print(f"Number of taxa          : {otu_counter}")
print(f"Sum of abundance        : {total_abundance}")
print(f"Sum of estimated counts : {total_estimated_counts}")
print(f"Output file             : {output_file}")
print("----------------------------------------")
