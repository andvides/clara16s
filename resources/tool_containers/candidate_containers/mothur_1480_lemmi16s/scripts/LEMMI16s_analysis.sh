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
parameters['maxambig']=0
parameters['minlength']=50
parameters['maxhomop']=8
parameters['diffs']=2
parameters['taxlevel']=4
parameters['cutoff']=0.03
parameters['label']=0.03

#updating parameters fom the config
IFS=, read -a fields <<< $aux_parameters
for i in "${fields[@]}"; do 
  key=`echo "$i" | cut -d'=' -f1`
  value=`echo "$i" | cut -d'=' -f2`
  parameters[$key]=$value
done

mkdir -p tmp_$sample_folder
cd tmp_$sample_folder
directory=`pwd`
if [ -z ${cpus+x} ]; then cpus=$(nproc);fi

echo "running"
cp $query/queryReads1.fq queryReads1.fastq
cp $query/queryReads2.fq queryReads2.fastq

reference=$model_path'/reference.fasta'
taxonomy=$model_path'/reference.tsv'
export PROC=$cpus

#Clean seqs
mothur "#make.file(inputdir=$directory, type=fastq, prefix=stability);
        make.contigs(file=stability.files,processors=$PROC);
        screen.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table, processors=$PROC, 
                    maxambig=${parameters['maxambig']}, minlength=${parameters['minlength']}, maxhomop=${parameters['maxhomop']});
        unique.seqs(fasta=stability.trim.contigs.good.fasta, count=stability.contigs.good.count_table);
        pre.cluster(fasta=stability.trim.contigs.good.unique.fasta, count=stability.trim.contigs.good.count_table, diffs=${parameters['diffs']},processors=$PROC);
        chimera.vsearch(fasta=stability.trim.contigs.good.unique.precluster.fasta, count=stability.trim.contigs.good.unique.precluster.count_table, dereplicate=t, processors=$PROC);
        rename.file(fasta=current, count=current, taxonomy=current, prefix=clean)"

#Classify seq
mothur "#classify.seqs(fasta=clean.fasta, count=clean.count_table, reference=$reference, taxonomy=$taxonomy, processors=$PROC);
        rename.file(fasta=current, count=current, taxonomy=current, prefix=final)"
#OTU
mothur "#cluster.split(fasta=final.fasta, count=final.count_table, taxonomy=final.taxonomy, taxlevel=${parameters['taxlevel']}, cutoff=${parameters['cutoff']}, processors=$PROC);
        make.shared(list=final.opti_mcc.list, count=final.count_table, label=${parameters['level']});
        classify.otu(list=final.opti_mcc.list, count=final.count_table, taxonomy=final.taxonomy, label=${parameters['level']})"

editResult.py final.opti_mcc.0.03.cons.taxonomy results.tsv summary_taxonomy.tsv
cp results.tsv $output_file



