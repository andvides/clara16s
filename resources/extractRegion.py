#!/usr/bin/env python3
import argparse
from Bio import SeqIO
import os
import sys

fasta_file=sys.argv[1]
output_file=sys.argv[2]
region=sys.argv[3]

with open(output_file, 'w+') as file:
    pass

try:
  with open(fasta_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
      file = open(output_file+"_hyperex_tmp.fna", "w+")
      file.write('>'+record.description+'\n')
      file.write(str(record.seq))
      file.close()
      command = "hyperex -p "+output_file+"_hyperex_out --quiet "+region+" "+output_file+"_hyperex_tmp.fna > "+output_file+"_hyperex_log || true"
      os.system(command)
      command = "cat "+output_file+"_hyperex_out.fa >> "+output_file+" && rm "+output_file+"_hyperex_*"
      os.system(command)
except FileNotFoundError:
  print(f"File path does not exist")
except Exception as e:
  print(f"An error occurred: {str(e)}")
