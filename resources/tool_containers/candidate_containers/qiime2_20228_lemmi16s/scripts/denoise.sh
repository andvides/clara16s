#!/usr/bin/env bash
set -o xtrace
set -e

inputFile=$1
outputFile=$2
  --type 'SampleData[PairedEndSequencesWithQuality]'  \
  --input-path $manifestFile \
  --output-path $outputFile \
  --input-format PairedEndFastqManifestPhred33 \

qiime cutadapt trim-paired \
    --i-demultiplexed-sequences $inputFile  \
    --p-cores 1 \
    --p-error-rate 0.1 \
    --p-overlap 3 \
    --verbose \
    --o-trimmed-sequences paired-end-demux-trimmed.qza

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs paired-end-demux-trimmed.qza \
    --p-n-threads 0 \
    --p-trunc-q 2 \
    --p-trunc-len-f 219 \
    --p-trunc-len-r 194 \
    --p-max-ee-f 2 \
    --p-max-ee-r 4 \
    --p-n-reads-learn 1000000 \
    --p-chimera-method none \
    --o-table table-dada2.qza \
    --o-representative-sequences $ $outputFile \
    --o-denoising-stats stats-dada2.qza

