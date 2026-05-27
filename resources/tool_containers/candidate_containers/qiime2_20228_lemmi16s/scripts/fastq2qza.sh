#!/usr/bin/env bash
set -o xtrace
set -e

queryPath=$1
outputFile=$2
manifestFile="manifest.tsv"

#manifest file
echo "sample-id,absolute-filepath,direction" > $manifestFile
echo "demo1,"$queryPath"/queryReads1.fq,forward">> $manifestFile
echo "demo1,"$queryPath"/queryReads2.fq,reverse">> $manifestFile
qiime tools import  \
  --type 'SampleData[PairedEndSequencesWithQuality]'  \
  --input-path $manifestFile \
  --output-path $outputFile \
  --input-format PairedEndFastqManifestPhred33 \

#qiime tools import \
#  --type 'FeatureData[Sequence]' \
#  --input-path $1 \
#  --output-path $2
