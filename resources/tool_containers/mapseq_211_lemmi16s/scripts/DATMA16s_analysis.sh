#!/usr/bin/env bash
set -o xtrace
set -e
technology=$1
sample_folder=$2
query=$3
model_path=$4
output_file=$5
aux_parameters=$6

#Default tool parameters
declare -A parameters
parameters['AVGQUAL']=20
parameters['HEADCROP']=12
parameters['LEADING']=3
parameters['TRAILING']=3
parameters['MINLEN']=75

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

if [ "$technology" = "illumina" ]; then
  echo "running illumina"
  #Tutorial from https://research.csc.fi/metagenomics_quality
  #Clean and merge 

  #java -jar /usr/src/Trimmomatic/0.38/Trimmomatic-0.38/trimmomatic-0.38.jar PE $query/queryReads1.fq $query/queryReads2.fq output_forward_paired.fastq \
    #output_forward_unpaired.fastq output_reverse_paired.fastq output_reverse_unpaired.fastq \
    #-phred33 -threads $cpus AVGQUAL:${parameters['AVGQUAL']} HEADCROP:${parameters['HEADCROP']} LEADING:${parameters['LEADING']} TRAILING:${parameters['TRAILING']} MINLEN:${parameters['MINLEN']}
    #AVGQUAL:20 HEADCROP:12 LEADING:3 TRAILING:3 MINLEN:75 -phred33 -threads $cpus

  #seqprep -f output_forward_paired.fastq -r output_reverse_paired.fastq -1 sample1_seqprep_R1.fastq.gz -2 sample1_seqprep_R2.fastq.gz -s sample1_seqprep_merged.fastq.gz
else
  echo "running PacBio"
  cat $query/* > sample1_seqprep_merged.fastq.gz
  #gzip -c $query/queryReads.fq > sample1_seqprep_merged.fastq.gz
fi

#Classification
REF=$model_path/reference_mapseq.fna
if [ -f "$REF" ]; then
  echo "Using local reference"
  mapseq -fastq sample1_seqprep_merged.fastq.gz -nthreads $cpus $model_path/reference_mapseq.fna $model_path/reference_mapseq.tax > results.mseq
else
  echo "Using default MAPseq reference"
  mapseq -fastq sample1_seqprep_merged.fastq.gz -nthreads $cpus > results.mseq
fi

#Write the results in the LEMMI16s format
echo -e "Taxonomy\tTaxonomyLevel\tLabel\tCounts" > results.txt
mapseq -otucounts results.mseq | grep "^0" | awk '{if($2==6) print}' >> results.txt 
editResult.py results.txt  results.tsv summary_taxonomy.tsv

cp results.tsv $output_file 




