#!/usr/bin/env python3
import csv
import sys
from ete3 import NCBITaxa
from tqdm import tqdm

# --------------------
# Input files
# --------------------
TAX_FILE=sys.argv[2]
FASTA_IN=sys.argv[1]
#TAX_FILE = "reference.tsv"
#FASTA_IN = "reference.fasta"
FASTA_OUT = "reference_kraken_format.fasta"

# --------------------
# Load NCBI taxonomy
# --------------------
ncbi = NCBITaxa()

# --------------------
# Load taxonomy table
# --------------------
taxonomy = {}
with open(TAX_FILE) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        taxonomy[row["Feature ID"]] = row["Taxon"]

# --------------------
# Helper functions
# --------------------
def extract_tax_names(lineage):
    """
    Converts:
    'd__Bacteria; p__Bacteroidota; ...; g__Aureispira; s__Aureispira marina'
    into:
    ['Aureispira marina', 'Aureispira', 'Saprospiraceae', ...]
    ordered deepest-first.
    """
    parts = [x.strip() for x in lineage.split(";")]
    names = []

    # Reverse (species → genus → family → ...)
    for p in reversed(parts):
        if "__" not in p:
            continue
        prefix, name = p.split("__", 1)
        name = name.strip()
        if name not in ["", "uncultured", "metagenome"]:
            names.append(name)

    return names


def get_taxid(seqid):
    lineage = taxonomy.get(seqid)
    if not lineage:
        return 1   # root

    names = extract_tax_names(lineage)

    # Try each name until one exists in NCBI
    for name in names:
        try:
            hit = ncbi.get_name_translator([name])
            if name in hit:
                return hit[name][0]
        except:
            pass

    # fallback: assign "root"
    return 1


# --------------------
# Rewrite FASTA
# --------------------
print("Rewriting FASTA with NCBI taxids...")

with open(FASTA_IN) as fin, open(FASTA_OUT, "w") as fout:
    for line in fin:
        if line.startswith(">"):
            seqid = line[1:].strip()
            taxid = get_taxid(seqid)
            fout.write(f">{seqid}|kraken:taxid|{taxid}\n")
        else:
            fout.write(line)

print("\nDone!")
print(f"Output FASTA: {FASTA_OUT}")
