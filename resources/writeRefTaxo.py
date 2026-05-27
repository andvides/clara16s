#!/usr/bin/env python3
import sys, gzip
import os

fastaFile=open(sys.argv[1], 'r')
taxonomy=open(sys.argv[2], 'r')
outp =open(sys.argv[3], 'w+')
header = sys.argv[4]+'\n' if len(sys.argv) > 4 else "Feature ID\tTaxon\n"
#Reading fasta file
seq=set()
for line in fastaFile:
  if line.startswith('>'):
    line=line.split()
    read=line[0].split(">")[1]
    seq.add(read)
#print(seq)


#Reading taxo file
taxo={}
for line in taxonomy:
  read=line.split()[0]
  taxo[read]=line
#print(seq)

#Output
outp.write(str(header))
for read in seq:
  outp.write(taxo[read])
outp.close()
