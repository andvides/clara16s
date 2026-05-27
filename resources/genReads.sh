#!/usr/bin/env bash
set -o xtrace
set -e

taxonDist=$1
read_length_distribution=$2
chp=$3
rs=$4
fasta=$5
output=$6

DIR=`echo $output |xargs dirname`
#index=`basename   $output .query `
index=`basename   $output .sim1.fq `

coverage=`grep  "$index" $taxonDist | awk -F"\t" '{print $(NF)}' | bc`

if [[ $coverage -gt 0 ]]
then
  #primer_forward=`head -1 $fasta | cut -f3 -d" " |cut -f2 -d"="`
  #primer_reverse=`head -1 $fasta | cut -f4 -d" " |cut -f2 -d"="`
  #cutadapt -a $primer_forward...$primer_reverse -o $output.trimm $fasta
  #art_illumina -nf 0 -na -q -p $read_length_distribution -rs $rs -c $coverage -qL 30 -i $output.trimm -o $DIR/$index.sim
  art_illumina -amp -nf 0 -na -q -p $read_length_distribution -rs $rs -c $coverage -i $fasta -o $DIR/$index.sim

  #grinder -rd $read_length_distribution -od $DIR -bn $index.query -tr $tr -cp $chp -rs $rs -rf $fasta
  #mv $output-reads.fa $output
  #mv $output-ranks.txt $output.txt
else
  touch $DIR/${index}.sim1.fq
  touch $DIR/${index}.sim2.fq
  #touch $output.txt
fi
