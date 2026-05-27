#!/usr/bin/env python3
import sys, gzip
import os
import re

#Input and output files
fastaFile=open(sys.argv[1], 'r')
taxonomy=open(sys.argv[2], 'r')
outp =open(sys.argv[3], 'w+')

#Reading taxo file
taxo={}
for line in taxonomy:
  line=line.split('\n')[0]
  read=line.split()[0]
  taxonomy=line.split(read)[1]
  taxonomy=taxonomy.replace(" ", "")
  taxonomy=taxonomy.replace("\t", "")
  #taxonomy=re.sub("d__|p__|c__|o__|f__|g__|s__","",taxonomy)
  taxo[read]=taxonomy

#Reading fasta file
for line in fastaFile:
  if line.startswith('>'):
    line=line.split()
    read=line[0].split(">")[1]
    line='>'+taxo[read]+'\n'
  outp.write(line) 
