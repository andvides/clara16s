#!/usr/bin/env bash
set -o xtrace
set -e

repository=$1
Database=`basename   $repository `
Database="${Database^^}"

release=$2
version=`basename   $release | cut -f2 -d '_' ` 

output=$3

mkdir -p $repository
echo "Version 2025:Download files for : "$Database

if [[ "$Database" = "SILVA" ]]; then
  #Download NCBI tax ID
  wget -q -c -P $repository https://ftp.arb-silva.de//$release/Exports/taxonomy/tax_slv_ssu_${version}.acc_taxid.gz
  gunzip -f $repository/tax_slv_ssu_${version}.acc_taxid.gz 
  cp $repository/tax_slv_ssu_${version}.acc_taxid $repository/tax_id.txt 

  #Download the Taxonomy Rank file. This maps the taxonomic rank and taxonomy to the taxid.
  wget -q -c -P $repository https://www.arb-silva.de/fileadmin/silva_databases/$release/Exports/taxonomy/tax_slv_ssu_${version}.txt.gz
  gunzip -f $repository/tax_slv_ssu_${version}.txt.gz

  #Download the Taxonomy Map file. This maps the sequence Accessions to the Organism Name and Taxonomy IDs.
  wget -q -c -P $repository https://www.arb-silva.de/fileadmin/silva_databases/$release/Exports/taxonomy/taxmap_slv_ssu_ref_nr_${version}.txt.gz
  gunzip -f $repository/taxmap_slv_ssu_ref_nr_${version}.txt.gz

  #Download the Taxonomy Tree file. This file contains the hierarchical relationship of the taxonomy IDs in tree form.
  wget -q -c -P $repository https://www.arb-silva.de/fileadmin/silva_databases/$release/Exports/taxonomy/tax_slv_ssu_${version}.tre.gz
  gunzip -f $repository/tax_slv_ssu_${version}.tre.gz

  #Download the SILVA NR99 sequences (non-redundant and unaligned)
  wget -q -c -P $repository https://www.arb-silva.de/fileadmin/silva_databases/$release/Exports/SILVA_${version}_SSURef_NR99_tax_silva_trunc.fasta.gz
  gunzip -f $repository/SILVA_${version}_SSURef_NR99_tax_silva_trunc.fasta.gz
  
  #tax dump
  wget -P $repository https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
  #Rename the output file to standard name
  #mv $repository/SILVA_${version}_SSURef_NR99_tax_silva_trunc.fasta $repository/$output.fasta
  
elif [[ "$Database" = "NCBI" ]]; then
  #Download NCBI tax ID the other files are downloaded using qiime 2 plugin
  wget -q -c -P $repository https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
  gunzip -f $repository/nucl_gb.accession2taxid.gz
  #Temporal file to skip the rule
  touch $repository/$output.fasta
  
elif [[ "$Database" = "GREENGENES2" ]]; then

  #version=`basename   $release | sed -r 's/[.]+/_/g' ` 
  BASE_URL='https://ftp.microbio.me/greengenes_release/'${version}

  # Download taxonomy map
  #wget -q -c -P $repository ${BASE_URL}/${version}.taxonomy.asv.tsv.qza
  wget -q -c -P $repository ${BASE_URL}/${version}.backbone.tax.qza
  mv $repository/${version}.backbone.tax.qza  $repository/raw_sequences_taxonomy.qza
  touch $repository/raw_sequences_taxonomy.tsv

  # Download representative sequences (full-length 16S + fragments)
  #wget -q -c -P $repository ${BASE_URL}/${version}.seqs.fna.qza
  wget -q -c -P $repository ${BASE_URL}/${version}.backbone.full-length.fna.qza
  mv $repository/${version}.backbone.full-length.fna.qza  $repository/raw_sequences.qza
  touch $repository/raw_sequences.fasta

  # Download phylogeny files
  wget -q -c -P $repository ${BASE_URL}/${version}.taxonomy.id.nwk
  mv $repository/${version}.taxonomy.id.nwk  $repository/raw_taxonomy.id.nwk

  #tax dump
  wget -P $repository https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
  tar xvzf $repository/taxdump.tar.gz -C $repository

elif [[ "$Database" = "GTDB" ]]; then
  version=`basename   $release | cut -f2 -d '_' ` 
  release=`basename   $release | cut -f1 -d '_' ` 
  rversion=`basename  $version | cut -f1 -d '.' `

  #tax dump
  BASE_URL="https://data.gtdb.ecogenomic.org/releases/$release/$version/auxillary_files/"
  html_content=$(curl -s $BASE_URL)
  dates=$(echo "$html_content" | grep -oP 'ncbi_taxdump_\d{8}.tar.gz' | sed -r 's/ncbi_taxdump_([0-9]{8}).tar.gz/\1/' | sort | uniq)
  latest_date=$(echo "$dates" | tail -n 1)
  if [[ -n "$latest_date" ]]; then
    echo "Downloaded ncbi_taxdump_${latest_date}.tar.gz"
    wget -c -P $repository "${BASE_URL}ncbi_taxdump_${latest_date}.tar.gz"
    mv $repository/ncbi_taxdump_${latest_date}.tar.gz  $repository/taxdump.tar.gz
  else
    echo "No date found in the repository."
  fi

  #Download the sequences example "release214_214.0" or "release220_220.0" 
  ext1=tar
  ext2=tar 
  ext3="-C $repository"
  tar_gunzip="tar xvzf"
  
  if (( $(echo "$version >= 220" | bc -l) )); then
   ext1=fna
   ext2=tsv
   ext3=""
   tar_gunzip="gunzip -f "
  fi 
  
  wget -c -P $repository https://data.gtdb.ecogenomic.org/releases/$release/$version/genomic_files_all/ssu_all_r$rversion.$ext1.gz
  $tar_gunzip $repository/ssu_all_r$rversion.$ext1.gz $ext3

  #Download metadata
  wget -c -P $repository https://data.gtdb.ecogenomic.org/releases/$release/$version/bac120_metadata_r$rversion.$ext2.gz
  $tar_gunzip $repository/bac120_metadata_r$rversion.$ext2.gz $ext3

  wget -c -P $repository https://data.gtdb.ecogenomic.org/releases/$release/$version/ar53_metadata_r$rversion.$ext2.gz
  $tar_gunzip $repository/ar53_metadata_r$rversion.$ext2.gz $ext3

  #Generate taxonomy
  editReference.py  $Database $repository/ssu_all_r$rversion.fna $repository/bac120_metadata_r$rversion.tsv,$repository/ar53_metadata_r$rversion.tsv \
                    $repository/$output.fasta $repository/$output.taxonomy

  #TaxID 
  grep '>' $repository/$output.fasta | awk -F "\t" '{print $1"\t"$2}' | cut -f2 -d'>' > $repository/tax_id.txt 
      
elif [[ "$Database" = "rrnDB" ]]; then
  #Download the sequences 
  wget -q -c -P $repository https://rrndb.umms.med.umich.edu/static/download/rrnDB-${release}_16S_rRNA.fasta.zip
  zcat $repository/rrnDB-${release}_16S_rRNA.fasta.zip  > $repository/$output.tmp1
  
  #Download metadata
  wget -q -c -P $repository https://rrndb.umms.med.umich.edu/static/download/rrnDB-${release}.tsv.zip
  zcat $repository/rrnDB-${release}.tsv.zip  > $repository/$output.tmp2
  
  #Generate taxonomy
  editReference.py  $Database $repository/$output.tmp1 $repository/$output.tmp2  $repository/$output.fasta $repository/$output.taxonomy

  #TaxID 
  grep '>' $repository/$output.fasta | awk -F "\t" '{print $1"\t"$2}' | cut -f2 -d'>' > $repository/tax_id.txt 

  #tax dump
  wget -q -c -P $repository ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz  
else
    echo "Failed, please select a valid database error 2025"
    exit 1
fi
