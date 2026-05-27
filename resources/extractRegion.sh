#!/usr/bin/env bash
set -o xtrace
set -e
fasta=$1
output=$2
region=$3

min_seq_length=10 #To skip hyperex panicked at 'slice index starts at ' error

if [ -f $output.hyperex ]; then
  rm $output.*
fi
hyperex -p $output --quiet $region $fasta > $output.hyperex || true
touch $output.fa
seqkit seq -g -m $min_seq_length -o $output $output.fa  
#mv $output.fa  $output

