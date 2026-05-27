#!/usr/bin/env python3
import sys, gzip
import os
import numpy as np
from collections import Counter

fasta=open(sys.argv[1], 'r')
taxid=open(sys.argv[2], 'r')

#KF697569.1.1227 d__Bacteria; p__Chloroflexi; c__Anaerolineae; o__Caldilineales; f__Caldilineaceae; g__uncultured; s__uncultured_bacterium
reads={}
for line in taxid:
  if line.startswith("Feature"):
    continue
  else:
    line=line.split('\n')[0]
    reference=line.split('\t')[0]
    taxonomy =line.split('\t')[1]
    reads[reference]=taxonomy

#@CP007206.2008927.2010444-1800/1
#GTGCCAGAAACCGAGATAATACTGTGGATCCAAGGGTTATCCGGAATCATTGGGTTTAAAGGGTCCGTAGGCGGCCAGATAAGTCAGTGGTGAAAGCCCA
#+
#F6"EEFF"2"FFF"F"FEE=FF"="EFDFFFFFF"EF*FF+EFFF8F4FBED>FFDF0FFFCFFFFFEFF*)FFFDFFE;FFEFFFFFFFDFFB8F@FF;
outp=open(sys.argv[3]+".all", 'w+')
toprint="#READ\tTAXONOMY\n"
outp.write(toprint)
i=0
taxa=[]
for line in fasta:
  if i%4==0:
    name="-".join(line.split("-")[:-1]).split('@')[1]
    #name=line.split()[0].split('-')[0].split('@')[1]
    toprint=line.split()[0].split('/')[0].split('@')[1]+'\t'+reads[name]+'\n'
    taxa.append(reads[name])
    outp.write(toprint)
  i+=1
outp.close()

outp=open(sys.argv[3], 'w+')
toprint="#LEMMI16s\nGroup_ID\tSize\tTaxonomy\n"
outp.write(toprint)
count = Counter(taxa)
i=1
for items in count.keys():
  toprint=str(i)+"\t"+str(count[items])+"\t"+items+"\n"
  i+=1
  outp.write(toprint)
outp.close()
