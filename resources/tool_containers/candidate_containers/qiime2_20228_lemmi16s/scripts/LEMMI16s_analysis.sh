#!/usr/bin/env bash
set -o xtrace
set -e
sample_folder=$1
query=$2
model_path=$3
output_file=$4
aux_parameters=$5

#Default tool parameters
declare -A parameters
parameters['type']='SampleData[PairedEndSequencesWithQuality]'
parameters['format']='PairedEndFastqManifestPhred33'
parameters['trim_length']=120
parameters["tax_level"]=7
parameters["q_score"]=20

#updating parameters fom the config
IFS=, read -a fields <<< $aux_parameters
for i in "${fields[@]}"; do 
  key=`echo "$i" | cut -d'=' -f1`
  value=`echo "$i" | cut -d'=' -f2`
  parameters[$key]=$value
done

mkdir -p tmp_$sample_folder
cd tmp_$sample_folder

if [ -z ${cpus+x} ]; then cpus=$(nproc);fi

#Tutorial from: https://docs.qiime2.org/2022.2/tutorials/moving-pictures/
echo "running"
echo "sample-id,absolute-filepath,direction" > manifest.tsv
#Deblur cannot operate on sample IDs that contain underscores.
echo -e "${sample_folder//_/-},$query/queryReads1.fq,forward" >> manifest.tsv
echo -e "${sample_folder//_/-},$query/queryReads2.fq,reverse" >> manifest.tsv

qiime tools import \
  --type ${parameters["type"]} \
  --input-path manifest.tsv \
  --output-path demux.qza  \
  --input-format ${parameters["format"]}

qiime quality-filter q-score \
  --i-demux demux.qza \
  --p-min-quality ${parameters["q_score"]} \
  --o-filtered-sequences demux-filtered.qza \
  --o-filter-stats demux-filter-stats.qza

qiime deblur denoise-16S \
  --i-demultiplexed-seqs demux-filtered.qza   \
  --p-trim-length ${parameters["trim_length"]}   \
  --o-representative-sequences rep-seqs-deblur.qza \
  --o-table table-deblur.qza \
  --p-sample-stats \
  --o-stats deblur-stats.qza

mv rep-seqs-deblur.qza rep-seqs.qza
mv table-deblur.qza table.qza


qiime feature-classifier classify-sklearn  \
  --i-classifier $model_path/referenceClassificator.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level ${parameters["tax_level"]}  \
  --o-collapsed-table phyla-table.qza

qiime tools export --input-path phyla-table.qza --output-path .
biom convert -i feature-table.biom -o results.txt --to-tsv

editResult.py results.txt results.tsv summary_taxonomy.tsv
cp results.tsv $output_file

