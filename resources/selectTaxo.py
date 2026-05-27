#!/usr/bin/env python3
import sys, gzip
import re 
import os

organism=sys.argv[1].split()
taxoFile=open(sys.argv[2], 'r')
outputList=open(sys.argv[3], 'w+')
outputTaxo=open(sys.argv[3]+".tsv", 'w+')

for line in taxoFile:
  if line.startswith('Feature'):
    continue
  else:
    for items in organism:
      if re.search(items, line, re.I): 
        data=line.split()[0]+"\n"
        outputTaxo.write(line)
        outputList.write(data)

outputList.close()
outputTaxo.close()
