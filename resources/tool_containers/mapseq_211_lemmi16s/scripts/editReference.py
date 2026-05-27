#!/usr/bin/env python3
import sys, gzip
import os
import re

#Input and output files
fastaFile=open(sys.argv[1], 'r')
taxonomy=open(sys.argv[2], 'r')
outp =open(sys.argv[3], 'w+')

#Reading taxo file
print("#cutoff: 0.00:0.08 0.70:0.35 0.70:0.35 0.70:0.35 0.80:0.25 0.92:0.08 0.95:0.05")
print("#name: NCBI")
print("#levels: Kingdom Phylum Class Order Family Genus Species")
taxo={}
for line in taxonomy:
  line=line.split('\n')[0]
  """read=line.split()[0]
  taxonomy=line.split(read)[1]"""
  read, taxonomy = line.split(maxsplit=1)
  taxonomy=taxonomy.replace("; ", ";")
  taxonomy=taxonomy.replace("\t", "")
  taxonomy=re.sub("d__|p__|c__|o__|f__|g__|s__","",taxonomy)
  print(read+"\t"+taxonomy)

#Reading fasta file
for line in fastaFile:
  if line.startswith('>'):
    line=line.split()
    read=line[0].split(">")[1]
    line='>'+read+'\n'
  outp.write(line) 
