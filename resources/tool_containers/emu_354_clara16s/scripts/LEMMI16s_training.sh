#!/usr/bin/env bash
set -o xtrace
set -e
sample_folder=$1
dataset=$2
database=$3
model_path=$4
output_file=$5

mkdir -p tmp_$sample_folder/tmp
mkdir -p tmp_$sample_folder/database
mkdir -p tmp_$sample_folder/database/taxonomy

cd tmp_$sample_folder

database_name=`basename $database | cut -f1 -d ':'`

echo 'emu version 2025'
source activate custom-emu

if [ -z ${cpus+x} ]; then cpus=$(nproc);fi

if [ "$database_name" = "SILVA" ]; then
  echo -e "id\tlineage" > ex-taxonomy2.tsv
  awk '{print $2"\t"$1}' $database/tax_slv_ssu_*.txt >> ex-taxonomy2.tsv
    emu build-database database --sequences $dataset/reference.fasta \
    --seq2tax $database/tax_id.txt \
    --taxonomy-list ex-taxonomy2.tsv
elif [ "$database_name" = "GREENGENES2" ]; then
  # gg2kdb-2025.r --path_in $dataset --n_cores $cpus
  # outdir=taxonomy
  # k2db=database
  
  # # fix delimiters! 
  # sed -i -e 's/$/\t|/g' $outdir/names.dmp
  # sed -i -e 's/$/\t|/g' $outdir/nodes.dmp
  
  # emu build-database  $k2db --sequences $dataset/reference.fasta \
  #   --seq2tax $outdir/tax_id.txt \
  #   --ncbi-taxonomy $outdir

  emu build-database  database --sequences $dataset/reference.fasta \
  --seq2tax $database/tax_id.txt \
  --taxonomy-list  $dataset/reference.tsv

else
   tar zxvf $database/taxdump.tar.gz -C  taxonomy/
   emu build-database  database --sequences $dataset/reference.fasta \
    --seq2tax $database/tax_id.txt \
    --ncbi-taxonomy taxonomy  
fi
