#!/usr/bin/env bash
set -o xtrace
set -e

#To choose a specific organism, the name in the config file (taxo option) must coincide with the name in the repository/Database_name/taxonomy.tsv file                  
#Te scripts lists all the reads that contains the taxonIDS in the taxonomy.tsv.
#SelectFasta recovers the fasta file from the original database for each taxonIDS

taxonIDS=$1
taxonomy=$2
output_list=$3
dbs=$4
output_fasta=$5

selectTaxo.py $taxonIDS $taxonomy $output_list
selectFasta -list  $output_list -fasta $dbs -fasta_sel > $output_fasta
