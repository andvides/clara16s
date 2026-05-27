#!/usr/bin/env python3

import argparse
import pandas as pd


def find_start_line(file):
    """
    Detecta la línea donde empiezan los datos reales (OTU_/ASV_)
    """
    with open(file) as f:
        for i, line in enumerate(f):
            if line.startswith("OTU_") or line.startswith("ASV_"):
                return i
    return 0


def parse_tax(tax):
    """
    Convierte taxonomía tipo:
    d__Bacteria;p__Firmicutes;...
    a lista limpia para Krona
    """
    levels = []
    for t in str(tax).split(";"):
        if "__" in t:
            val = t.split("__", 1)[1]
        else:
            val = t

        if val.strip() == "" or val.lower() == "na":
            val = "Unclassified"

        levels.append(val)

    return levels


def main():
    parser = argparse.ArgumentParser(
        description="Convert taxonomy table (OTU/ASV) to Krona format"
    )
    parser.add_argument("input", help="Input TSV file")
    parser.add_argument("-o", "--output", default="krona_input.txt",
                        help="Output file (default: krona_input.txt)")

    args = parser.parse_args()

    # Detectar inicio de datos
    start = find_start_line(args.input)

    # Leer archivo
    df = pd.read_csv(
        args.input,
        sep="\t",
        skiprows=start,
        names=["Group_ID", "Size", "Taxonomy"]
    )

    # Procesar filas
    rows = []
    for _, row in df.iterrows():
        try:
            abundance = int(row["Size"])
        except:
            continue  # saltar filas problemáticas

        tax_levels = parse_tax(row["Taxonomy"])
        rows.append([abundance] + tax_levels)

    # Guardar salida
    out = pd.DataFrame(rows)
    out.to_csv(args.output, sep="\t", index=False, header=False)

    print(f"[OK] Krona input generado: {args.output}")
    print(f"[INFO] Filas procesadas: {len(rows)}")


if __name__ == "__main__":
    main()