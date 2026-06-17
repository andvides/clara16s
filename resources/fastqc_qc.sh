#!/usr/bin/env bash

set -euo pipefail

########################################
# Usage
########################################
usage() {
cat << EOF
Usage:
  $(basename "$0") \
      -i <reads_dir> \
      -t <illumina|pacbio|nanopore> \
      -o <output_dir> \
      -p <threads>

Example:
  $(basename "$0") \
      -i ./fastq \
      -t illumina \
      -o ./qc \
      -p 16

EOF
exit 1
}

########################################
# Parse arguments
########################################
READ_DIR=""
TECH=""
OUTDIR=""
THREADS=""

while getopts "i:t:o:p:h" opt; do
    case "$opt" in
        i) READ_DIR="$OPTARG" ;;
        t) TECH="${OPTARG,,}" ;;
        o) OUTDIR="$OPTARG" ;;
        p) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

[[ -z "$READ_DIR" ]] && usage
[[ -z "$TECH" ]] && usage
[[ -z "$OUTDIR" ]] && usage
[[ -z "$THREADS" ]] && usage

########################################
# Validate technology
########################################
case "$TECH" in
    illumina|pacbio|nanopore)
        ;;
    *)
        echo "ERROR ${TECH}: Technology must be one of:"
        echo "  illumina"
        echo "  pacbio"
        echo "  nanopore"
        exit 1
        ;;
esac

########################################
# Check dependencies
########################################
command -v fastqc >/dev/null 2>&1 || {
    echo "ERROR: fastqc not found"
    exit 1
}

mkdir -p "$OUTDIR"
FASTQC_DIR="${OUTDIR}/fastqc"

mkdir -p "$FASTQC_DIR"

########################################
# Find reads
########################################

if [[ -f "$READ_DIR" ]]; then

    READS=("$READ_DIR")

elif [[ -d "$READ_DIR" ]]; then

    mapfile -t READS < <(
    find "$READ_DIR" -type f \
    \( -name "*.fastq" \
    -o -name "*.fq" \
    -o -name "*.fastq.gz" \
    -o -name "*.fq.gz" \)
    )

else
    echo "Input does not exist: $READ_DIR"
    exit 1
fi

########################################
# Run FastQC
########################################
echo "======================================="
echo "Running FastQC"
echo "Technology: $TECH"
echo "Files found: ${#READS[@]}"
echo "======================================="

fastqc \
    --threads "$THREADS" \
    --outdir "$FASTQC_DIR" \
    "${READS[@]}"

########################################
# Metrics file
########################################
METRICS="${OUTDIR}/qc_metrics.tsv"

echo -e \
"Sample\tTechnology\tTotal_Sequences\tSequence_Length\tPercent_GC" \
> "$METRICS"

########################################
# Summary file
########################################
SUMMARY="${OUTDIR}/qc_summary.tsv"

echo -e \
"Sample\tModule\tStatus" \
> "$SUMMARY"

########################################
# Parse FastQC output
########################################
for ZIP in "$FASTQC_DIR"/*_fastqc.zip
do
    SAMPLE=$(basename "$ZIP" _fastqc.zip)

    TMPDIR=$(mktemp -d)

    unzip -qq "$ZIP" -d "$TMPDIR"

    DATA=$(find "$TMPDIR" -name fastqc_data.txt)
    SUM=$(find "$TMPDIR" -name summary.txt)

    TOTAL=$(grep "^Total Sequences" "$DATA" | cut -f2)
    LEN=$(grep "^Sequence length" "$DATA" | cut -f2)
    GC=$(grep "^%GC" "$DATA" | cut -f2)

    echo -e \
    "${SAMPLE}\t${TECH}\t${TOTAL}\t${LEN}\t${GC}" \
    >> "$METRICS"

    while IFS=$'\t' read -r STATUS MODULE FILE
    do
        echo -e \
        "${SAMPLE}\t${MODULE}\t${STATUS}" \
        >> "$SUMMARY"
    done < "$SUM"

    rm -rf "$TMPDIR"
done


########################################
# Technology-specific notes
########################################
case "$TECH" in

illumina)
cat << EOF

Recommended review:
  - Per base sequence quality
  - Adapter content
  - Duplication levels

EOF
;;

pacbio)
cat << EOF

Recommended review:
  - Read length distribution
  - GC content
  - HiFi read yield

FastQC quality plots are less informative
for PacBio HiFi data.

EOF
;;

nanopore)
cat << EOF

Recommended review:
  - Read length distribution
  - GC content

Consider:
  NanoPlot
  NanoStat
  PycoQC

for deeper Nanopore QC.

EOF
;;

esac

echo
echo "======================================="
echo "FastQC completed"
echo "======================================="
echo
echo "FastQC reports:"
echo "  ${FASTQC_DIR}"
echo
echo "Metrics:"
echo "  ${METRICS}"
echo
echo "Summary:"
echo "  ${SUMMARY}"
echo