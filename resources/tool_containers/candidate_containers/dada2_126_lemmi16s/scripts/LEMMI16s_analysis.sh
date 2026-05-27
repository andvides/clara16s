#!/usr/bin/env bash
set -o xtrace
set -e

technology=$1
sample_folder=$2
query=$3
model_path=$4
output_file=$5
aux_parameters=$6


mkdir -p tmp_$sample_folder
cd tmp_$sample_folder
directory=`pwd`

if [ -z ${cpus+x} ]; then cpus=$(nproc);fi

#Tutorial from https://benjjneb.github.io/dada2/tutorial.html
if [ "$technology" = "illumina" ]; then
    echo "running illumina"
    cp $query/queryReads1.fq .
    cp $query/queryReads2.fq .
    runClassifier.R queryReads_$i.fasta $model_path/reference.fasta $cpus $directory "$aux_parameters"
else
    echo "running PacBio"
    cp $query/queryReads.fq .
    runClassifier-longReads.R queryReads.fasta $model_path/reference.fasta $cpus $directory "$aux_parameters"
fi

editResult.py results.txt results.tsv summary_taxonomy.tsv
cp results.tsv  $output_file

