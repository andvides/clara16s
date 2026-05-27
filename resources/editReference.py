#!/usr/bin/env python3
import sys, gzip
import os
import pandas as pd
import tarfile
import numpy as np 

database =sys.argv[1]
fastaFile=sys.argv[2]
taxonomy =sys.argv[3]
outFasta =open(sys.argv[4], 'w+')
outTaxa  =open(sys.argv[5], 'w+')


if database == "rrnDB":
  fastaFile_data=open(fastaFile, 'r')
  taxonomy_data =open(taxonomy, 'r')

  #Reading taxo file
  taxa = pd.read_table(taxonomy_data)
  taxa = taxa[['Data source record id','NCBI tax id','RDP taxonomic lineage','Data source organism name']]
  taxa.set_index('Data source record id',inplace=True)
  #taxa.to_csv('borrar.txt')  
  #taxa = taxa.to_dict()
  #print(taxa)

  #Reading fasta file
  #>Pseudomonas aeruginosa|GCF_009676885.1|NZ_CP046069.1|Chromosome: ANONYMOUS|711785..713295 +
  i=1
  for line in fastaFile_data:
    if line.startswith('>'):
      read=line.split("|")[1]
      taxid=str(taxa.loc[read]['NCBI tax id'])
      lineage=taxa.loc[read]['RDP taxonomic lineage']
      lineage= lineage.split("|||")[1] if "|||" in lineage else lineage
      specie=taxa.loc[read]['Data source organism name']
      taxo='d__'+lineage+'; '+specie
      taxo=taxo.replace("|domain; ", ";p__")
      taxo=taxo.replace("|phylum; ", ";c__")
      taxo=taxo.replace("|class; ",  ";o__")
      taxo=taxo.replace("|order; ",  ";f__")
      taxo=taxo.replace("|family; ", ";g__")
      taxo=taxo.replace("|genus; ",  ";s__")
      taxo=taxo.replace(" ",  "_")
      taxo=taxo.replace("|",  "_") 
      line='>R'+str(i)+'_'+read+" "+taxid+' '+taxo+'\n'
      printTaxa='R'+str(i)+'_'+read+'\t'+taxo+'\n'
      outTaxa.write(printTaxa)
      i+=1
    #print(line.split('\n')[0])
    outFasta.write(line)

elif database == "GTDB":
  #Reading metadata bacteria
  taxID={}
  metadata=open(taxonomy.split(',')[0], 'r')
  for line in metadata:
    if line.startswith('accession'):
      continue
    line=line.split("\t")
    read=  line[0]
    taxid= line[77]
    taxID[read]=taxid
    
    
  #Reading metadata arc
  metadata=open(taxonomy.split(',')[1], 'r')
  for line in metadata:
    if line.startswith('accession'):
      continue
    line=line.split("\t")
    read=  line[0]
    taxid= line[77]
    taxID[read]=taxid

  #Reading fasta file
  fastaFile_data=open(fastaFile, 'r')
  i=1
  for line in fastaFile_data:
    if line.startswith('>'):
      read=line.split("~")[0].split('>')[1]
      taxo=line.split()[1]
      line='>R'+str(i)+'_'+read+" "+str(taxID[read])+" "+taxo.replace(" ",  "_")+'\n'
      printTaxa='R'+str(i)+'_'+read+'\t'+taxo+'\n'
      outTaxa.write(printTaxa)
      i+=1
    outFasta.write(line)
    

    
