#!/usr/bin/env bash
set -o xtrace
set -e
sample_folder=$1
dataset=$2
database=$3
load_model=$4
model=$5
aux_parameters=$6

mkdir -p tmp_$sample_folder
cd tmp_$sample_folder

if [ -z ${cpus+x} ]; then cpus=$(nproc);fi

if [ $load_model -gt 0 ]; then
  wget -c -q $model -O model.tar.gz
  tar xvzf model.tar.gz
  fasta2qza.sh  custom_reference.fasta reference.qza
  taxo2qza.sh custom_reference.tsv  referenceTaxo.qza
  trainClass.sh reference.qza referenceTaxo.qza referenceClassificator
else
  fasta2qza.sh  $dataset/reference.fasta reference.qza
  taxo2qza.sh $dataset/reference.tsv  referenceTaxo.qza
  trainClass.sh reference.qza referenceTaxo.qza referenceClassificator
fi
