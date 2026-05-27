#!/usr/bin/env bash
set -o xtrace
set -e

fasta=$1
train_rate=$2
random_seed=$3

index=`basename   $fasta .reg`
#query_size=`grep -F "$index" $taxonDist | awk -F"\t" '{print $4}' | bc` #number of amplicons for the query dataset
query_size=` grep -c '>' $fasta | awk -v train_rate=$train_rate '{printf "%.0f\n", $1*train_rate/100}' | bc`

if [[ $query_size -gt 0 ]]
then
  selectFasta -seed $random_seed -random $query_size -fasta $fasta > $fasta.query
  grep ">" $fasta.query > $fasta.list
  selectFasta -list  $fasta.list -fasta $fasta > $fasta.ref
else
  cp $fasta $fasta.ref
  touch $fasta.query
fi


