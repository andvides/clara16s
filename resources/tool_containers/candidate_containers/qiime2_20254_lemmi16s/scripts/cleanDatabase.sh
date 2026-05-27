#!/usr/bin/env bash
set -o xtrace
set -e

#Script to remove sequences from the original dataset
#Because Qiime2 is used, this file must be part of the Dockerfile for Qiime2
#See https://github.com/bokulich-lab/RESCRIPt for details

repository=$1
Database=`basename   $repository `
release=$2
version=`basename   $release | cut -f2 -d '_' ` 

cd $repository
echo "Preparing the reference database"

if [ -z ${cpus+x} ]; then cpus=$(nproc);fi

if [ "$Database" = "SILVA" ]; then
  #Fasta format to qiime structure
  qiime tools import \
    --type 'FeatureData[SILVATaxonomy]' \
    --input-path tax_slv_ssu_${version}.txt \
    --output-path taxranks-silva-${version}-ssu-nr99.qza \
    
  qiime tools import \
    --type 'FeatureData[SILVATaxidMap]' \
    --input-path taxmap_slv_ssu_ref_nr_${version}.txt \
    --output-path taxmap-silva-${version}-ssu-nr99.qza \
    
  qiime tools import \
    --type 'Phylogeny[Rooted]' \
    --input-path tax_slv_ssu_${version}.tre \
    --output-path taxtree-silva-${version}-nr99.qza \
    
  qiime tools import \
    --type 'FeatureData[RNASequence]' \
    --input-path SILVA_${version}_SSURef_NR99_tax_silva_trunc.fasta \
    --output-path silva-${version}-ssu-nr99-rna-seqs.qza
    
  qiime rescript reverse-transcribe \
    --i-rna-sequences silva-${version}-ssu-nr99-rna-seqs.qza \
    --o-dna-sequences raw_sequences.qza
    
  qiime rescript parse-silva-taxonomy \
    --i-taxonomy-tree taxtree-silva-${version}-nr99.qza \
    --i-taxonomy-map taxmap-silva-${version}-ssu-nr99.qza \
    --i-taxonomy-ranks taxranks-silva-${version}-ssu-nr99.qza \
    --p-include-species-labels \
    --o-taxonomy raw_sequences_taxonomy.qza
    
elif [ "$Database" = "NCBI" ]; then
  BioProjects=`echo $version | sed "s/,/[BioProject] OR /g"`
  BioProjects+=[BioProject]
  qiime rescript get-ncbi-data \
    --p-query "$BioProjects"\
    --o-sequences raw_sequences.qza \
    --o-taxonomy raw_sequences_taxonomy.qza
else
  qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path raw_sequences.fasta\
    --output-path raw_sequences.qza

  qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path raw_sequences.taxonomy \
    --output-path raw_sequences_taxonomy.qza
fi

#Culling low-quality sequences with cull-seqs
qiime rescript cull-seqs \
  --i-sequences raw_sequences.qza \
  --o-clean-sequences raw_sequences-cleaned.qza \
  --p-n-jobs $cpus
    
#Filtering sequences by length and taxonomy
qiime rescript filter-seqs-length-by-taxon \
  --i-sequences raw_sequences-cleaned.qza \
  --i-taxonomy  raw_sequences_taxonomy.qza \
  --p-labels Archaea Bacteria Eukaryota \
  --p-min-lens 900 1200 1400 \
  --o-filtered-seqs raw_sequences-filt.qza \
  --o-discarded-seqs raw_sequences-discard.qza
    
#Dereplication of sequences and taxonomy
qiime rescript dereplicate \
  --i-sequences raw_sequences-filt.qza  \
  --i-taxa raw_sequences_taxonomy.qza \
  --p-mode 'uniq' \
  --o-dereplicated-sequences raw_sequences-derep-uniq.qza \
  --o-dereplicated-taxa raw_sequences-tax-derep-uniq.qza
    
#Export from qiime structure to fasta and tsv format 
#The dna-sequences.fasta and taxonomy.tsv files are automatically generated
qiime tools export --input-path raw_sequences-derep-uniq.qza --output-path .
qiime tools export --input-path raw_sequences-tax-derep-uniq.qza --output-path .   






