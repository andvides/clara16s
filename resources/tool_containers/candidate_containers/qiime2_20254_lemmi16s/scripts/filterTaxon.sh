#!/usr/bin/env bash
set -o xtrace
set -e
silvaRef=$1
silvaTaxo=$2
filter=$3
len=$4
outName=$5
qiime rescript filter-seqs-length-by-taxon \
    --i-sequences $silvaRef \
    --i-taxonomy $silvaTaxo \
    --p-labels $filter \
    --p-max-lens $len \
    --o-filtered-seqs other$outName.qza \
    --o-discarded-seqs $outName.qza 
