#!/usr/bin/env bash
set -o xtrace
set -e

taxonIDS_unknown=$1
baseDir=$2
taxonomy=$3
fasta=$4
output_fasta=$5
output_tsv=$6


cat $baseDir/*/*.reg.list | cut -f1 -d' ' | cut -d'>' -f2 >> $baseDir/usedReads.txt 
for organisms in $1
do
 grep $organisms  $taxonomy | cut -f1 | cut -f1 -d'-' >> $baseDir/usedReads.txt
done
selectFasta -list $baseDir/usedReads.txt -fasta $fasta > $output_fasta
#hyperex -p $output --quiet $region $output_fasta > $output_fasta.hyperex || true
writeRefTaxo.py  $output_fasta  $taxonomy $output_tsv
