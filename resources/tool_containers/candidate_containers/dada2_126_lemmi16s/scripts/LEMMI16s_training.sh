#!/usr/bin/env bash
set -o xtrace
set -e
sample_folder=$1
dataset=$2
database=$3
load_model=$4
model=$5
aux_parameters=$6

database_name=`basename $database | cut -f1 -d ':'` 

mkdir -p tmp_$sample_folder
cd tmp_$sample_folder

if [ -z ${cpus+x} ]; then cpus=$(nproc);fi


if [ $load_model -gt 0 ]; then
  wget -c -q $model -O model.tar.gz
  tar xvzf model.tar.gz 
  editReference.py custom_reference.fasta custom_reference.tsv reference.fasta
else
  editReference.py $dataset/reference.fasta $dataset/reference.tsv reference.fasta
fi

