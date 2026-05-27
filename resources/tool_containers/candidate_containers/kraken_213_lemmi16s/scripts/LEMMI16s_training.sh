#!/usr/bin/env bash
set -o xtrace
set -e
sample_folder=$1
dataset=$2
database=$3

mkdir -p tmp_$sample_folder/database
mkdir -p tmp_$sample_folder/database/taxonomy

cd tmp_$sample_folder

database_name=`basename $database | cut -f1 -d ':'`

pwd
echo $USER

if [ -z ${cpus+x} ]; then cpus=$(nproc);fi


if [ $database_name=="SILVA" ]; then
   cp -r /silvaDB/taxonomy database/
else
   tar zxvf $database/taxdump.tar.gz -C  database/taxonomy/
fi

editReference.py $database/tax_id.txt $dataset/reference.fasta $dataset/reference.tsv > reference_taxID.txt

kraken2-build --threads $cpus --add-to-library reference_kraken_format.fasta --db database
kraken2-build --threads $cpus --build --db database
/bracken2/bracken-build -d database -t $cpus -k 35 -l 150

