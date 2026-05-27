#!/usr/bin/env bash
set -o xtrace
set -e
cpus=$1
classifier=$2
query=$3
output=$4

mkdir -p results_tax
qiime feature-classifier classify-sklearn \
    --i-classifier $classifier \
    --i-reads $query \
    --p-n-jobs $cpus \
    --o-classification $output.qza 

qiime metadata tabulate \
    --m-input-file $output.qza \
    --o-visualization $output.qzv

qiime tools export \
  --input-path $output.qzv \
  --output-path ./results_tax/

cat results_tax/metadata.tsv > $output.tsv


#unzip -q $output.qzv -d results_aba 
#mv results_aba/*/data/metadata.tsv $output.tsv
#rm -rf results_aba
