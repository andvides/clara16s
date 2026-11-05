#!/usr/bin/env bash
set -o xtrace
set -e
technology=$1
sample_folder=$2
query=$3
model_path=$4
output_file=$5

mkdir -p tmp_$sample_folder
cd tmp_$sample_folder
source activate custom-emu 

if [ -z ${cpus+x} ]; then cpus=$(nproc);fi

if [ "$technology" = "illumina" ]; then
    echo "running illumina"
else
    echo "running PacBio or Oxford"
    emu abundance $query/queryReads.fq  --db $model_path/database --threads $cpus
fi

total=$(( $(wc -l < ${query}/queryReads.fq) / 4 ))
editResults.py results/queryReads_rel-abundance.tsv $total results.tsv
cp results.tsv $output_file 




