#!/usr/bin/env python3
import sys, gzip
import os
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'

taxo_result=sys.argv[1]
out1       =open(sys.argv[2], 'w+')
out2       =open(sys.argv[3], 'w+')


#Reading taxo file
"""
# Constructed from biom file
#OTU ID Sample1 taxonomy
OTU1    258.0   k__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria; o__Pseudomonadales_A_650611; f__Pseudomonadaceae; g__Pseudomonas_E_647464; s__?

"""
taxonomy = pd.read_csv(str(taxo_result), header=1, sep="\t")
taxonomy["taxonomy"]=taxonomy["taxonomy"].str.replace('; ',';')
taxonomy["taxonomy"]=taxonomy["taxonomy"].str.replace('k__','d__')
taxonomy["taxonomy"]=taxonomy["taxonomy"].str.replace('?','')
otus = taxonomy["taxonomy"].str.split(';',expand=True)
otus.insert(0, "Size", taxonomy["Sample1"].astype(int))
#otus.insert(0, "OTU", taxonomy["OTU"])
otus.columns = ["Size","domain","phylum","class","order","family","genus","species"]
#print(otus)

#Write results
out1.write("#LEMMI16s\n")
out1.write("Group_ID\tSize\tTaxonomy\n")
for index, row in otus.iterrows():
  toprint='group_'+str(index)+"\t"+str(row["Size"])+'\t'+ \
          row["domain"]+';'+row["phylum"]+';'+row["class"]+';'+ \
          row["order"]+';'+row["family"]+';'+row["genus"]+';'+row["species"]+'\n'
  out1.write(toprint)

#Write summary taxonomy
out2.write("#LEMMI16s\norganism\treads\tabundance\n")
for taxo in ["domain","phylum","class","order","family","genus","species"]:
  df=otus[[taxo,'Size']]
  df[taxo]=df[taxo]
  total=df['Size'].sum()
  reads_abundance = df.groupby([taxo])['Size'].sum()
  abundance=reads_abundance/total
  df=pd.DataFrame({"reads":reads_abundance,"abundance":abundance})
  df.to_csv(out2, header=False,sep='\t')

