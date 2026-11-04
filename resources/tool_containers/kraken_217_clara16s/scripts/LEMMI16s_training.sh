#!/usr/bin/env bash
set -Eeuo pipefail

# Arguments:
#   $1 sample_folder
#   $2 dataset directory
#   $3 database/reference directory
#   $4 load_model: 1 = use pre-built database, 0 = build database
#   $5 model: pre-built Kraken2 database path
#   $6 aux_parameters
#
# Examples:
#   kraken2-build.threads=12
#   bracken-build.k=35
#   bracken-build.l=1000
#   kraken2-build.threads=12,bracken-build.k=35,bracken-build.l=1000

sample_folder=$1
dataset=$2
database=$3
load_model="${4:-0}"
model="${5:-}"
aux_parameters="${6:-none=none}"

if [[ "$load_model" != "0" && "$load_model" != "1" ]]; then
    echo "ERROR: load_model must be 0 or 1." >&2
    exit 2
fi

if [[ "$load_model" == "1" && -z "$model" ]]; then
    echo "ERROR: model path is required when load_model=1." >&2
    exit 2
fi

cpus="${cpus:-$(nproc)}"
SILVA_TAXONOMY_DIR="${SILVA_TAXONOMY_DIR:-/silvaDB/taxonomy}"

if [[ "$load_model" == "0" ]]; then
    dataset="$(realpath "$dataset")"
    database="$(realpath "$database")"
else
    model="$(realpath "$model")"
    [[ -d "$model" ]] || {
        echo "ERROR: model directory not found: $model" >&2
        exit 2
    }
fi

work_root="tmp_${sample_folder}"
mkdir -p "$work_root"
cd "$work_root"

database_name=""
if [[ "$load_model" == "0" ]]; then
    database_name="$(basename "$database" | cut -f1 -d ':')"
fi

declare -A K2_COMMON=(
    [threads]="$cpus"
    [db]="database"
)

declare -A K2_ADD=()
declare -A K2_BUILD=()

declare -A BRACKEN=(
    [d]="database"
    [t]="$cpus"
    [k]="35"
    [l]="150"
)

bracken_db_overridden=false

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

set_option() {
    local namespace="$1"
    local option="$2"
    local value="$3"

    option="$(normalize_option "$option")"

    [[ -n "$option" ]] || {
        echo "ERROR: empty option in aux_parameters." >&2
        exit 2
    }

    case "$namespace" in
        kraken2-build)
            K2_COMMON["$option"]="$value"
            ;;
        kraken2-build.add-to-library)
            K2_ADD["$option"]="$value"
            ;;
        kraken2-build.build)
            K2_BUILD["$option"]="$value"
            ;;
        bracken-build)
            BRACKEN["$option"]="$value"

            if [[ "$option" == "d" || "$option" == "database" ]]; then
                bracken_db_overridden=true
            fi
            ;;
        *)
            echo "ERROR: unsupported namespace '$namespace'." >&2
            exit 2
            ;;
    esac
}

parse_aux_parameters() {
    local entry key value namespace option
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

        if [[ "$key" == kraken2-build.add-to-library.* ]]; then
            namespace="kraken2-build.add-to-library"
            option="${key#kraken2-build.add-to-library.}"
        elif [[ "$key" == kraken2-build.build.* ]]; then
            namespace="kraken2-build.build"
            option="${key#kraken2-build.build.}"
        elif [[ "$key" == kraken2-build.* ]]; then
            namespace="kraken2-build"
            option="${key#kraken2-build.}"
        elif [[ "$key" == bracken-build.* ]]; then
            namespace="bracken-build"
            option="${key#bracken-build.}"
        elif [[ "$key" == kraken2.* || "$key" == bracken.* ]]; then
            # Analysis-only parameters are ignored during training.
            continue
        else
            echo "ERROR: unknown parameter namespace in '$key'." >&2
            echo "Training accepts kraken2-build.* and bracken-build.*." >&2
            exit 2
        fi

        # Remove optional quotes around values.
        if [[ ${#value} -ge 2 ]]; then
            first_char="${value:0:1}"
            last_char="${value: -1}"

            if [[ ( "$first_char" == '"' && "$last_char" == '"' ) ||
                  ( "$first_char" == "'" && "$last_char" == "'" ) ]]; then
                value="${value:1:${#value}-2}"
            fi
        fi

        # Prevent overriding workflow-managed actions.
        case "$namespace:$option" in
            kraken2-build:add-to-library|\
            kraken2-build:build|\
            kraken2-build:clean|\
            kraken2-build:download-taxonomy|\
            kraken2-build:download-library|\
            kraken2-build:standard|\
            kraken2-build:special|\
            kraken2-build.add-to-library:add-to-library|\
            kraken2-build.build:build)
                echo "ERROR: '$key' is managed by the workflow." >&2
                exit 2
                ;;
        esac

        set_option "$namespace" "$option" "$value"
    done
}

append_options() {
    local tool="$1"
    local map_name="$2"
    local key value flag
    local -n options_ref="$map_name"
    local -a keys=()
    local k

    for k in "${!options_ref[@]}"; do
        keys+=("$k")
    done

    if ((${#keys[@]})); then
        while IFS= read -r k; do
            key="$k"
            value="${options_ref[$key]}"

            [[ "$value" == "false" ]] && continue

            if [[ "$tool" == "bracken-build" ]]; then
                if [[ "$key" == -* ]]; then
                    flag="$key"
                elif ((${#key} == 1)); then
                    flag="-$key"
                else
                    flag="--$key"
                fi
            else
                flag="--$key"
            fi

            if [[ "$value" == "true" ]]; then
                CLI_ARGS+=("$flag")
            elif [[ -n "$value" ]]; then
                CLI_ARGS+=("$flag" "$value")
            else
                echo "ERROR: empty value for '$tool' option '$key'." >&2
                exit 2
            fi
        done < <(printf '%s\n' "${keys[@]}" | sort)
    fi
}

parse_aux_parameters

db_dir="${K2_COMMON[db]}"

if [[ "$load_model" == "1" ]]; then

    # Validate the pre-built Kraken2 database.
    for index_file in hash.k2d opts.k2d taxo.k2d; do
        [[ -s "$model/$index_file" ]] || {
            echo "ERROR: incomplete Kraken2 database '$model'; missing $index_file." >&2
            exit 2
        }
    done

    # Create symbolic link to the expected database location.
    mkdir -p "$(dirname "$db_dir")"

    if [[ -L "$db_dir" ]]; then
        current_target="$(readlink "$db_dir")"
        if [[ "$current_target" != "$model" ]]; then
            rm -- "$db_dir"
        fi
    elif [[ -e "$db_dir" ]]; then
        if [[ -d "$db_dir" && -z "$(ls -A "$db_dir")" ]]; then
            rmdir -- "$db_dir"
        else
            echo "ERROR: '$db_dir' already exists and is not a symbolic link." >&2
            echo "Move or remove it before using load_model=1." >&2
            exit 2
        fi
    fi

    if [[ ! -L "$db_dir" ]]; then
        ln -s -- "$model" "$db_dir"
    fi

    echo "Using pre-built Kraken2 database: $model"
    echo "Model link created: $(pwd)/$db_dir -> $model"
    echo "Training skipped because load_model=1."
    exit 0
fi

db_dir="${K2_COMMON[db]}"
mkdir -p "$db_dir" "$db_dir/taxonomy" "$db_dir/library" tmp

if [[ "$bracken_db_overridden" == false ]]; then
    BRACKEN[d]="$db_dir"
fi

echo "Working directory: $(pwd)"
echo "CPUs: $cpus"
echo "Kraken2 database: $db_dir"
echo "Auxiliary parameters: $aux_parameters"

if [[ "$database_name" == "GREENGENES2" ]]; then
    gg2kdb.r --path_in "$dataset" --n_cores "$cpus"

    outdir="tmp"

    cp "$outdir"/*nodes.dmp "$db_dir/taxonomy/nodes.dmp"
    cp "$outdir"/*names.dmp "$db_dir/taxonomy/names.dmp"
    cp "$outdir"/*seqs.fna reference_kraken_format.fasta

    sed -i -e 's/$/\t|/g' "$db_dir/taxonomy/names.dmp"
    sed -i -e 's/$/\t|/g' "$db_dir/taxonomy/nodes.dmp"

else
    [[ -f "$database/taxdump.tar.gz" ]] || {
        echo "ERROR: taxonomy dump not found: $database/taxdump.tar.gz" >&2
        exit 1
    }

    tar -xzf "$database/taxdump.tar.gz" -C "$db_dir/taxonomy/"

    rm -f reference_kraken_format.fasta reference_taxID.txt

    editReference.py \
        "$database/tax_id.txt" \
        "$dataset/reference.fasta" \
        "$dataset/reference.tsv" > reference_taxID.txt

    [[ -s reference_kraken_format.fasta ]] || {
        echo "ERROR: reference_kraken_format.fasta is missing or empty." >&2
        exit 2
    }

    [[ -s reference_taxID.txt ]] || {
        echo "ERROR: reference_taxID.txt is missing or empty." >&2
        exit 2
    }
fi

CLI_ARGS=()
append_options "kraken2-build" K2_COMMON
append_options "kraken2-build" K2_ADD

kraken2-build "${CLI_ARGS[@]}" \
    --add-to-library reference_kraken_format.fasta

CLI_ARGS=()
append_options "kraken2-build" K2_COMMON
append_options "kraken2-build" K2_BUILD

kraken2-build "${CLI_ARGS[@]}" --build

CLI_ARGS=()
append_options "bracken-build" BRACKEN

bracken-build "${CLI_ARGS[@]}"
