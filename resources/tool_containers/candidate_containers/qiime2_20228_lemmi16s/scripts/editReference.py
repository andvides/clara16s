#!/usr/bin/env python3
import sys, gzip
import os

fastaFile=open(sys.argv[1], 'r')
taxonomy=open(sys.argv[2], 'r')
#outp =open(sys.argv[3], 'w+')

#Reading taxo file
taxo={}
for line in taxonomy:
  line=line.split('\n')[0]
  #line=line.replace(';', '\t')
  read=line.split()[0]
  taxo[read]=line.split(read)[1]
#print(seq)

#Reading fasta file
seq=set()
i=1
for line in fastaFile:
  if line.startswith('>'):
    line=line.split('-')[0]
    read=line.split(">")[1]
    line='>R'+str(i)+'_'+read+"\t"+taxo[read]
    i+=1
    #seq.add(read)
  print(line.split('\n')[0])
#print(seq)



#Output
#for read in seq:
  #line=read+'\t'+taxo[read]+'\n'
#  outp.write(taxo[read])
#outp.close()
