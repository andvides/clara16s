#!/usr/bin/env bash
set -o xtrace
set -e

queryfasta=$1
queryTaxo=$2
outdir=$3
numReads=$4
dbRef=$5
cpus=$6
params=$7
source activate mhass

mhass --amplicon-fasta $queryfasta --amplicon-genome-labels $queryTaxo --output-dir $outdir --num-reads $4  --threads $cpus $params
mv $outdir/combined_reads.fastq $outdir/queryReads.fq 
mhass_longRead_taxoRef.py $outdir/sequence_file_mapping.tsv $dbRef  $outdir/queryTaxo.tsv
source deactivate 
