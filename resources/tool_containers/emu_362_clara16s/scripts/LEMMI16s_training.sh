
#!/usr/bin/env bash

set -o xtrace
set -e

# --------------------------------------------------
# Arguments
# --------------------------------------------------

sample_folder=$1
dataset=$2
database=$3
load_model=$4
model=$5
aux_parameters=${6:-none=none}

# --------------------------------------------------
# Environment
# --------------------------------------------------

echo "Starting Emu database training"
echo "Sample       : $sample_folder"
echo "Dataset      : $dataset"
echo "Database     : $database"
echo "Load model   : $load_model"
echo "Model        : $model"
echo "Aux params   : $aux_parameters"

source activate emu

if [ -z "${cpus+x}" ]; then
    cpus=$(nproc)
fi

# --------------------------------------------------
# Prepare working directory
# --------------------------------------------------

work_dir="tmp_${sample_folder}"

mkdir -p "$work_dir/tmp"
mkdir -p "$work_dir/database"
mkdir -p "$work_dir/database/taxonomy"

# Resolve paths before changing directory
dataset=$(realpath "$dataset")
database=$(realpath "$database")

cd "$work_dir"

# --------------------------------------------------
# Parse additional parameters
# --------------------------------------------------

aux_array=()

if [[ -n "$aux_parameters" &&
      "$aux_parameters" != "none" &&
      "$aux_parameters" != "none=none" ]]; then

    read -r -a aux_array <<< "$aux_parameters"
fi

echo "Additional Emu parameters:"
printf '  %s\n' "${aux_array[@]:-none}"

# --------------------------------------------------
# Identify database
# --------------------------------------------------

database_name=$(basename "$database" | cut -f1 -d ':')

echo "Database type: $database_name"

# --------------------------------------------------
# Load existing model or build a new one
# --------------------------------------------------

if [[ "$load_model" == "1" ]]; then

    echo "Using pre-built Emu database"

    if [[ "$model" == "none" || ! -d "$model" ]]; then
        echo "ERROR: Invalid pre-built model path: $model"
        exit 1
    fi

    model=$(realpath "$model")

    echo "Model path: $model"

elif [[ "$load_model" == "0" ]]; then

    echo "Building Emu database"

    if [[ "$database_name" == "SILVA" ]]; then

        echo "Preparing SILVA taxonomy"

        echo -e "id\tlineage" > ex-taxonomy2.tsv

        awk '{print $2"\t"$1}' \
            "$database"/tax_slv_ssu_*.txt \
            >> ex-taxonomy2.tsv

        emu build-database database \
            --sequences "$dataset/reference.fasta" \
            --seq2tax "$database/tax_id.txt" \
            --taxonomy-list ex-taxonomy2.tsv \
            "${aux_array[@]}"

    elif [[ "$database_name" == "GREENGENES2" ]]; then

        echo "Preparing Greengenes2 database"

        emu build-database database \
            --sequences "$dataset/reference.fasta" \
            --seq2tax "$database/tax_id.txt" \
            --taxonomy-list "$dataset/reference.tsv" \
            "${aux_array[@]}"

    else

        echo "Preparing NCBI taxonomy"

        tar -zxvf "$database/taxdump.tar.gz" -C taxonomy/

        emu build-database database \
            --sequences "$dataset/reference.fasta" \
            --seq2tax "$database/tax_id.txt" \
            --ncbi-taxonomy taxonomy \
            "${aux_array[@]}"

    fi

    model="$(pwd)/database"

    echo "Emu database successfully built"
    echo "Model path: $model"

else

    echo "ERROR: Invalid load_model value: $load_model"
    echo "Expected 0 (build model) or 1 (load model)"
    exit 1

fi

echo "Training completed successfully"
