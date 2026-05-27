#!/usr/bin/env bash
set -o xtrace
set -e
sample_folder=$1
query=$2
model_path=$3
output_file=$4

mkdir -p tmp_$sample_folder
cd tmp_$sample_folder

if [ -z ${cpus+x} ]; then cpus=$(nproc);fi


kraken2 --threads $cpus --db $model_path/database --output query.out --report query.report --paired $query/queryReads1.fq $query/queryReads2.fq 
/bracken2/bracken -d $model_path/database -i  query.report -o query.bracken -l 'G' -r 150

editResult.py query_bracken_genuses.report results.tsv summary_taxonomy.tsv
cp results.tsv $output_file 




