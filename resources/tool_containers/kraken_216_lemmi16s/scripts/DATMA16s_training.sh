#!/usr/bin/env bash
set -o xtrace
set -e
sample_folder=$1
dataset=$2
database=$3

mkdir -p tmp_$sample_folder/tmp
mkdir -p tmp_$sample_folder/database
mkdir -p tmp_$sample_folder/database/taxonomy

cd tmp_$sample_folder

database_name=`basename $database | cut -f1 -d ':'`
database_name="${database_name^^}"

pwd
echo $USER

if [ -z ${cpus+x} ]; then cpus=$(nproc);fi


if [ "$database_name" = "GREENGENES2" ]; then

   gg2kdb.r --path_in $database --n_cores $cpus
   outdir=tmp
   k2db=database
   # rm -rf $k2db
   mkdir -p $k2db $k2db/taxonomy $k2db/library

   cp $outdir/*nodes.dmp $k2db/taxonomy/nodes.dmp
   cp $outdir/*names.dmp $k2db/taxonomy/names.dmp
   cp $outdir/*seqs.fna  reference_kraken_format.fasta

   # fix delimiters! 
   sed -i -e 's/$/\t|/g' $k2db/taxonomy/names.dmp
   sed -i -e 's/$/\t|/g' $k2db/taxonomy/nodes.dmp
else
   tar zxvf $database/taxdump.tar.gz -C  database/taxonomy/
   editReference.py $database/tax_id.txt $database/dna-sequences.fasta $database/taxonomy.tsv > reference_taxID.txt
fi
kraken2-build --threads $cpus --add-to-library reference_kraken_format.fasta --db database
kraken2-build --threads $cpus --build --db database
/bracken2/bracken-build -d database -t $cpus -k 35 -l 150

