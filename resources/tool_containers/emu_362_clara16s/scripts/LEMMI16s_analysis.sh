
#!/usr/bin/env bash

set -e
set -o xtrace

# --------------------------------------------------
# Arguments
# --------------------------------------------------

technology=$1
sample_folder=$2
query=$3
model_path=$4
output_file=$5
aux_parameters=${6:-none=none}

# --------------------------------------------------
# Configuration
# --------------------------------------------------

echo "Starting Emu abundance estimation"
echo "Technology    : $technology"
echo "Sample        : $sample_folder"
echo "Query         : $query"
echo "Model path    : $model_path"
echo "Output file   : $output_file"
echo "Aux parameters: $aux_parameters"

# --------------------------------------------------
# Activate environment
# --------------------------------------------------

source activate emu

if [ -z "${cpus+x}" ]; then
    cpus=$(nproc)
fi

# --------------------------------------------------
# Resolve paths before changing directory
# --------------------------------------------------

query=$(realpath "$query")
model_path=$(realpath "$model_path")

mkdir -p "$(dirname "$output_file")"
output_file="$(realpath -m "$output_file")"

# --------------------------------------------------
# Identify Emu database
# --------------------------------------------------

if [[ -f "$model_path/species_taxid.fasta" &&
      -f "$model_path/taxonomy.tsv" ]]; then

    emu_database="$model_path"

elif [[ -f "$model_path/database/species_taxid.fasta" &&
        -f "$model_path/database/taxonomy.tsv" ]]; then

    emu_database="$model_path/database"

else
    echo "ERROR: Emu database not found in: $model_path"
    exit 1
fi

echo "Emu database: $emu_database"

# --------------------------------------------------
# Map sequencing technology to Emu type
# --------------------------------------------------

case "${technology,,}" in

    illumina)
        emu_type="sr"
        ;;

    nanopore|ont)
        emu_type="map-ont"
        ;;

    pacbio)
        emu_type="map-pb"
        ;;

    pacbio-hifi|hifi)
        emu_type="map-hifi"
        ;;

    *)
        echo "ERROR: Unsupported sequencing technology: $technology"
        exit 1
        ;;
esac

echo "Emu sequencing type: $emu_type"

# --------------------------------------------------
# Prepare working directory
# --------------------------------------------------

work_dir="tmp_${sample_folder}"

mkdir -p "$work_dir/results"
cd "$work_dir"

# --------------------------------------------------
# Parse additional parameters
# --------------------------------------------------

aux_array=()

if [[ -n "$aux_parameters" &&
      "$aux_parameters" != "none" &&
      "$aux_parameters" != "none=none" ]]; then

    read -r -a raw_array <<< "$aux_parameters"

    for arg in "${raw_array[@]}"; do
        [[ "$arg" == "none" ]] && continue
        aux_array+=("$arg")
    done
fi

echo "Additional Emu parameters:"
if [[ ${#aux_array[@]} -eq 0 ]]; then
    echo "  none"
else
    printf '  %s\n' "${aux_array[@]}"
fi

# --------------------------------------------------
# Identify query reads
# --------------------------------------------------

query_files=()

if [[ "$emu_type" == "sr" ]]; then

    if [[ -f "$query/queryReads_1.fastq" &&
          -f "$query/queryReads_2.fastq" ]]; then

        query_files=(
            "$query/queryReads_1.fastq"
            "$query/queryReads_2.fastq"
        )

    elif [[ -f "$query/queryReads.fq" ]]; then

        query_files=("$query/queryReads.fq")

    else
        echo "ERROR: Illumina query FASTQ files not found"
        exit 1
    fi

else

    if [[ -f "$query/queryReads.fq" ]]; then
        query_files=("$query/queryReads.fq")

    elif [[ -f "$query/queryReads.fastq" ]]; then
        query_files=("$query/queryReads.fastq")

    else
        echo "ERROR: Query FASTQ file not found in $query"
        exit 1
    fi

fi

echo "Input FASTQ files:"
printf '  %s\n' "${query_files[@]}"

# --------------------------------------------------
# Count input reads
# --------------------------------------------------

total=$(awk '
    NR % 4 == 1 && /^@/ { count++ }
    END { print count+0 }
' "${query_files[0]}")

echo "Input reads: $total"

# --------------------------------------------------
# Run Emu abundance estimation
# --------------------------------------------------

echo "Running Emu abundance estimation"

emu abundance "${query_files[@]}" \
    --db "$emu_database" \
    --type "$emu_type" \
    --threads "$cpus" \
    --output-dir results \
    --output-basename queryReads \
    --keep-counts \
    "${aux_array[@]}"

# --------------------------------------------------
# Format results
# --------------------------------------------------

emu_results="results/queryReads_rel-abundance.tsv"

if [[ ! -f "$emu_results" ]]; then
    echo "ERROR: Emu output file not found: $emu_results"
    exit 1
fi

echo "Formatting Emu results"

editResults.py \
    "$emu_results" \
    "$total" \
    results.tsv

# --------------------------------------------------
# Copy final output
# --------------------------------------------------

cp results.tsv "$output_file"

echo "Analysis completed successfully"
echo "Final output: $output_file"
