#!/usr/bin/env bash
set -o xtrace
set -e
technology=$1
sample_folder=$2
query=$3
model_path=$4
output_file=$5
aux_parameters=$6

source activate LotuS3

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
  mkdir -p reads
  cp $query/queryReads1.fq reads/.
  cp $query/queryReads2.fq reads/.  
  
  touch reads/miSeqMap.sm.txt
  echo -e'#SampleID\tfastqFile\tSequencingRun' >> reads/miSeqMap.sm.txt
  echo -e'Sample1\tqueryReads1.fq,queryReads2.fq\tRun1' >> reads/miSeqMap.sm.txt

  lotus3 -i reads/ -m reads/miSeqMap.sm.txt  \
         -o reads-output/ \
         -refDB $model_path/database/reference.fasta \
         -tax4refDB $model_path/database/reference.tsv
else
  echo "running PacBio"

  mkdir -p reads
  cp $query/queryReads.fq reads/
  touch reads/pacbioMap.sm.txt
  echo -e '#SampleID\tfastqFile\tSequencingRun' > reads/pacbioMap.sm.txt
  echo -e 'Sample1\tqueryReads.fq\tRun1' >> reads/pacbioMap.sm.txt
  cat /data/sdm_PacBio_LSSU.txt > reads/sdm_PacBio_LSSU.txt

  lotus3 -i reads/ -m reads/pacbioMap.sm.txt \
         -o reads-output \
         -refDB  $model_path/database/reference.fasta \
         -tax4refDB  $model_path/database/reference.tsv \
         -s reads/sdm_PacBio_LSSU.txt -p PacBio -t $cpus
fi
biom convert -i reads-output/OTU.biom -o results.txt --to-tsv --header-key taxonomy
source deactivate 

#Write the results in the LEMMI16s format 
editResult.py results.txt  results.tsv summary_taxonomy.tsv
cp results.tsv $output_file 



