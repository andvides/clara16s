#!/usr/bin/env bash
set -o xtrace
set -e
sample_folder=$1
dataset=$2
database=$3
use_full_ref=$4
load_model=$5

mkdir -p tmp_$sample_folder
cd tmp_$sample_folder

if [ -z ${cpus+x} ]; then cpus=$(nproc);fi

if [ $load_model -gt 0 ]; then
  wget -c https://data.qiime2.org/2023.5/common/silva-138-99-nb-classifier.qza
  mv silva-138-99-nb-classifier.qza referenceClassificator.qza
else
  if [ $use_full_ref -gt 0 ]; then
    fasta2qza.sh  $database/dna-sequences.fasta reference.qza
    taxo2qza.sh $database/taxonomy.tsv  referenceTaxo.qza
    trainClass.sh reference.qza referenceTaxo.qza referenceClassificator
  else
    fasta2qza.sh  $dataset/reference.fasta reference.qza
    taxo2qza.sh $dataset/reference.tsv  referenceTaxo.qza
    trainClass.sh reference.qza referenceTaxo.qza referenceClassificator
  fi
fi

#sleep 10
 
