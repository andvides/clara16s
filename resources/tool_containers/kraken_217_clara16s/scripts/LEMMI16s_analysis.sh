#!/usr/bin/env bash
set -Eeuo pipefail

# Arguments:
#   $1 technology: "illumina" for paired-end; otherwise single-end/PacBio
#   $2 sample_folder
#   $3 query directory
#   $4 model_path
#   $5 output_file
#   $6 aux_parameters
#
# Examples:
#   kraken2.threads=20,bracken.l="G"
#   kraken2.threads=16,kraken2.confidence=0.1,bracken.r=150
#   kraken2.quick=true,bracken.t=10

technology=$1
sample_folder=$2
query=$3
model_path=$4
output_file=$5
aux_parameters="${6:-none=none}"

query="$(realpath "$query")"
model_path="$(realpath "$model_path")"

if [[ "$output_file" != /* ]]; then
    output_file="$(realpath -m "$output_file")"
fi

cpus="${cpus:-$(nproc)}"

mkdir -p "tmp_${sample_folder}"
cd "tmp_${sample_folder}"

declare -A KRAKEN2_OPTIONS=(
    [threads]="$cpus"
)

declare -A BRACKEN_OPTIONS=(
    [l]="G"
    [r]="150"
)

trim() {
    local value="$1"
    value="${value#"${value%%[![:space:]]*}"}"
    value="${value%"${value##*[![:space:]]}"}"
    printf '%s' "$value"
}

normalize_option() {
    local option="$1"
    option="${option#--}"
    option="${option#-}"
    printf '%s' "$option"
}

set_aux_option() {
    local tool="$1"
    local option="$2"
    local value="$3"

    option="$(normalize_option "$option")"

    [[ -n "$option" ]] || {
        echo "ERROR: empty option in aux_parameters." >&2
        exit 2
    }

    case "$tool:$option" in
        kraken2:db|\
        kraken2:output|\
        kraken2:report|\
        kraken2:paired|\
        kraken2:classified-out|\
        kraken2:unclassified-out|\
        bracken:d|\
        bracken:i|\
        bracken:o|\
        bracken:w)
            echo "ERROR: '$tool.$option' is managed by LEMMI16s_analysis.sh." >&2
            exit 2
            ;;
    esac

    case "$tool" in
        kraken2)
            KRAKEN2_OPTIONS["$option"]="$value"
            ;;
        bracken)
            BRACKEN_OPTIONS["$option"]="$value"
            ;;
        *)
            echo "ERROR: unsupported tool '$tool'." >&2
            echo "Use kraken2.<option> or bracken.<option>." >&2
            exit 2
            ;;
    esac
}

parse_aux_parameters() {
    local entry key value tool option
    local first_char last_char
    local -a entries=()

    [[ -n "$aux_parameters" && "$aux_parameters" != "none=none" ]] || return 0

    IFS=',' read -r -a entries <<< "$aux_parameters"

    for entry in "${entries[@]}"; do
        entry="$(trim "$entry")"
        [[ -n "$entry" ]] || continue

        if [[ "$entry" != *=* ]]; then
            echo "ERROR: invalid parameter '$entry'; expected key=value." >&2
            exit 2
        fi

        key="$(trim "${entry%%=*}")"
        value="$(trim "${entry#*=}")"

        if [[ "$key" == kraken2.* ]]; then
            tool="kraken2"
            option="${key#kraken2.}"
        elif [[ "$key" == bracken.* ]]; then
            tool="bracken"
            option="${key#bracken.}"
        elif [[ "$key" == kraken2-build.* || "$key" == bracken-build.* ]]; then
            # Training-only parameters are ignored during analysis.
            continue
        else
            echo "ERROR: unknown parameter namespace in '$key'." >&2
            echo "Analysis accepts kraken2.* and bracken.*." >&2
            exit 2
        fi

        # Remove optional quotes around values, e.g. bracken.l="G".
        if [[ ${#value} -ge 2 ]]; then
            first_char="${value:0:1}"
            last_char="${value: -1}"

            if [[ ( "$first_char" == '"' && "$last_char" == '"' ) ||
                  ( "$first_char" == "'" && "$last_char" == "'" ) ]]; then
                value="${value:1:${#value}-2}"
            fi
        fi

        set_aux_option "$tool" "$option" "$value"
    done
}

append_options() {
    local tool="$1"
    local map_name="$2"
    local -n options_ref="$map_name"
    local -a keys=()
    local key value flag

    for key in "${!options_ref[@]}"; do
        keys+=("$key")
    done

    if ((${#keys[@]})); then
        while IFS= read -r key; do
            value="${options_ref[$key]}"

            [[ "$value" == "false" ]] && continue

            if [[ "$tool" == "kraken2" ]]; then
                flag="--$key"
            elif ((${#key} == 1)); then
                flag="-$key"
            else
                flag="--$key"
            fi

            if [[ "$value" == "true" ]]; then
                CLI_ARGS+=("$flag")
            elif [[ -n "$value" ]]; then
                CLI_ARGS+=("$flag" "$value")
            else
                echo "ERROR: empty value for $tool option '$key'." >&2
                exit 2
            fi
        done < <(printf '%s\n' "${keys[@]}" | sort)
    fi
}

parse_aux_parameters

echo "Technology: $technology"
echo "Sample: $sample_folder"
echo "Model: $model_path"
echo "Auxiliary parameters: $aux_parameters"

if [[ "$technology" == "illumina" ]]; then
    echo "Running Kraken2 for paired-end Illumina reads"

    read1="$query/queryReads1.fq"
    read2="$query/queryReads2.fq"

    [[ -f "$read1" && -f "$read2" ]] || {
        echo "ERROR: paired FASTQ files not found: $read1 and/or $read2" >&2
        exit 1
    }

    CLI_ARGS=()
    append_options "kraken2" KRAKEN2_OPTIONS

    kraken2 "${CLI_ARGS[@]}" \
        --db "$model_path/database" \
        --output query.out \
        --report query.report \
        --paired "$read1" "$read2"

else
    echo "Running Kraken2 for single-end/PacBio reads"

    reads="$query/queryReads.fq"

    [[ -f "$reads" ]] || {
        echo "ERROR: FASTQ file not found: $reads" >&2
        exit 1
    }

    CLI_ARGS=()
    append_options "kraken2" KRAKEN2_OPTIONS

    kraken2 "${CLI_ARGS[@]}" \
        --db "$model_path/database" \
        --output query.out \
        --report query.report \
        "$reads"
fi

CLI_ARGS=()
append_options "bracken" BRACKEN_OPTIONS

bracken "${CLI_ARGS[@]}" \
    -d "$model_path/database" \
    -i query.report \
    -o query.bracken \
    -w query_bracken_genuses.report

editResult.py \
    query_bracken_genuses.report \
    results.tsv \
    summary_taxonomy.tsv

cp results.tsv "$output_file"
