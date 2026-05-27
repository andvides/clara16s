#!/usr/bin/env python3
import sys, gzip
import os

fastaFile=open(sys.argv[1], 'r')
taxonomy=open(sys.argv[2], 'r')
outp =open(sys.argv[3], 'w+')

#Reading taxo file
taxo={}
for line in taxonomy:
  if "Feature" not in line:
    line=line.replace(" ", "")
    line=line.split('\n')[0]
    read=line.split()[0]
    taxo[read]=line.split(read)[1]
    print(read+'\t'+taxo[read].strip()+';') 
#print(seq)

#Reading fasta file
seq=set()
for line in fastaFile:
  if line.startswith('>'):
    line=line.split()
    read=line[0].split(">")[1]
    line='>'+read+taxo[read]+'\n'
    #seq.add(read)
  #print(line.split('\n')[0])
  outp.write(line)
#print(seq)



#Output
#for read in seq:
  #line=read+'\t'+taxo[read]+'\n'
#  outp.write(taxo[read])
#outp.close()
