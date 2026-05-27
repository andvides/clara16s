#!/usr/bin/env python3
import sys, gzip
open('reference_kraken_format.fasta', 'w').close()
outp=open('reference_kraken_format.fasta', 'a')

taxid_file=sys.argv[1]
fasta_file=sys.argv[2]
taxo_file =sys.argv[3]

#Read taxID
taxID={}
for line in open(taxid_file):
    line=line.split('\n')[0]
    read =line.split()[0]
    taxid=line.split()[1]
    taxID[read]=taxid

for line in open(fasta_file):
  if line.startswith('>'):
    read =line.split()[0]
    read=read.split('>')[1]
    if read in taxID:
      taxid=taxID[read]
    else:
      taxid='0'
    line = '>{}|kraken:taxid|{}\n'.format(read, taxid)
  outp.write(line)
outp.close()

for line in open(taxo_file):
  line=line.split('\n')[0]
  read =line.split()[0]
  taxonomy=str(line.split('\t')[-1])
  if read in taxID:
    taxid=taxID[read]
  else:
    taxid='0'
    taxonomy='d__unknow;p__unknow;c__unknow;o__unknow;f__unknow;g__unknow;s__unknow'
  print(read+'\t'+taxid+'\t'+taxonomy)
