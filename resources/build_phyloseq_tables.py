#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import sys


def parse_taxonomy(tax_string):
    """
    Convert:
    d__Bacteria;p__Firmicutes;...
    to dictionary with taxonomic ranks.
    """

    ranks = {
        "d": "Kingdom",
        "p": "Phylum",
        "c": "Class",
        "o": "Order",
        "f": "Family",
        "g": "Genus",
        "s": "Species",
    }

    result = {
        "Kingdom": "NA",
        "Phylum": "NA",
        "Class": "NA",
        "Order": "NA",
        "Family": "NA",
        "Genus": "NA",
        "Species": "NA",
    }

    for item in tax_string.split(";"):

        item = item.strip()

        if "__" not in item:
            continue

        prefix, value = item.split("__", 1)

        if prefix in ranks:
            result[ranks[prefix]] = value

    return result


def read_lemmi_file(filename):

    df = pd.read_csv(
        filename,
        sep="\t",
        comment="#"
    )

    expected = {"Group_ID", "Size", "Taxonomy"}

    if not expected.issubset(df.columns):
        raise ValueError(
            f"{filename}: expected columns "
            f"{expected}, found {list(df.columns)}"
        )

    return df


def main():

    if len(sys.argv) < 5:
        sys.exit(
            f"Usage:\n"
            f"{sys.argv[0]} "
            f"otu_table.tsv "
            f"taxonomy.tsv "
            f"feature_mapping.tsv "
            f"sample1.tab sample2.tab ..."
        )

    otu_output = sys.argv[1]
    taxonomy_output = sys.argv[2]
    mapping_output = sys.argv[3]

    input_files = sys.argv[4:]

    abundance = {}

    for infile in input_files:

        sample = (
            Path(infile)
            .stem
            .replace(".predictions", "")
        )
        sample=sample.split("-")[-1]
        print(f"Processing {infile}")

        df = read_lemmi_file(infile)

        for _, row in df.iterrows():

            tax = str(row["Taxonomy"]).strip()
            size = int(row["Size"])

            if tax not in abundance:
                abundance[tax] = {}

            abundance[tax][sample] = (
                abundance[tax].get(sample, 0)
                + size
            )

    # ---------------------------------------------------------
    # Feature IDs
    # ---------------------------------------------------------

    taxonomies = sorted(abundance.keys())

    feature_map = {
        tax: f"TAX_{i:06d}"
        for i, tax in enumerate(taxonomies, start=1)
    }

    # ---------------------------------------------------------
    # OTU TABLE
    # ---------------------------------------------------------

    otu = pd.DataFrame.from_dict(
        abundance,
        orient="index"
    ).fillna(0)

    otu = otu.astype(int)

    otu.index = [
        feature_map[tax]
        for tax in otu.index
    ]

    otu.index.name = "FeatureID"

    otu = otu.reset_index()

    otu.to_csv(
        otu_output,
        sep="\t",
        index=False
    )

    print(f"Wrote {otu_output}")

    # ---------------------------------------------------------
    # TAXONOMY TABLE
    # ---------------------------------------------------------

    taxonomy_rows = []

    for tax in taxonomies:

        row = parse_taxonomy(tax)

        row["FeatureID"] = feature_map[tax]

        taxonomy_rows.append(row)

    taxonomy_df = pd.DataFrame(taxonomy_rows)

    taxonomy_df = taxonomy_df[
        [
            "FeatureID",
            "Kingdom",
            "Phylum",
            "Class",
            "Order",
            "Family",
            "Genus",
            "Species",
        ]
    ]

    taxonomy_df.to_csv(
        taxonomy_output,
        sep="\t",
        index=False
    )

    print(f"Wrote {taxonomy_output}")

    # ---------------------------------------------------------
    # FEATURE MAPPING
    # ---------------------------------------------------------

    mapping_df = pd.DataFrame({
        "FeatureID": [
            feature_map[t]
            for t in taxonomies
        ],
        "FullTaxonomy": taxonomies
    })

    mapping_df.to_csv(
        mapping_output,
        sep="\t",
        index=False
    )

    print(f"Wrote {mapping_output}")


if __name__ == "__main__":
    main()
