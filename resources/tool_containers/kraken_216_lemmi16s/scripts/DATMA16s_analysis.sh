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

if [ -z ${cpus+x} ]; then cpus=$(nproc);fi

if [ "$technology" = "illumina" ]; then
    echo "running illumina"
    cat $query/1* > queryReads1
    cat $query/2* > queryReads2
    kraken2 --threads $cpus --db $model_path/database --output query.out --report query.report --paired queryReads1 queryReads2 
else
    echo "running PacBio"
    cat $query/* > queryReads
    kraken2 --threads $cpus --db $model_path/database --output query.out --report query.report queryReads
fi
/bracken2/bracken -d $model_path/database -i  query.report -o query.bracken -l 'G' -r 150

editResult.py query_bracken_genuses.report results.tsv summary_taxonomy.tsv
cp results.tsv $output_file 




