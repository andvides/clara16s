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
echo "intersect_bp" >  targets.out

if [ -z ${cpus+x} ]; then cpus=$(nproc);fi

#editReference.py ../../../datasets/$sample_folder/reference.fasta ../../../datasets/$sample_folder/reference.tsv  reference.fasta > reference.tsv
#sleep 30

if [ $load_model -gt 0 ]; then
  #Example wget https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset18_062020.pds.tgz
  #wget -c -q $model -O model.tgz
  #fasta=`tar -tf model.tgz | grep fasta`
  #taxo=`tar -tf model.tgz | grep tax`
  wget -c -q $model -O model.tar.gz
  tar xvzf model.tar.gz
  editReference.py custom_reference.fasta custom_reference.tsv reference.fasta > reference.tsv
else
  editReference.py $dataset/reference.fasta $dataset/reference.tsv  reference.fasta > reference.tsv    
fi
