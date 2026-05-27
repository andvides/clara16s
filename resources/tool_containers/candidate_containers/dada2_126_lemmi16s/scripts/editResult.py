#!/usr/bin/env python3
import sys, gzip
import os
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'

taxo_result=sys.argv[1]
out1       =open(sys.argv[2], 'w+')
out2       =open(sys.argv[3], 'w+')

#Reading Otu file
otus = pd.read_csv(str(taxo_result), header=0, sep=",")
otus.columns = ['Size',"domain","phylum","class","order","family","genus","species"]
otus.replace(np.nan,"na__na",inplace=True)
#print(otus)

#Write results
out1.write("#LEMMI16s\n")
out1.write("Group_ID\tSize\tTaxonomy\n")
for index, row in otus.iterrows():
  toprint='group_'+str(index)+"\t"+str(row["Size"])+'\t'+row["domain"]+';'+row["phylum"]+';'+row["class"]+';'+row["order"]+';'+row["family"]+';'+row["genus"]+';'+row["species"]+'\n'
  out1.write(toprint)

#Write summary taxonomy
out2.write("#LEMMI16s\norganism\treads\tabundance\n")
for taxo in ["domain","phylum","class","order","family","genus","species"]:
  df=otus[[taxo,'Size']]
  df[taxo].replace("na__na",taxo[0]+"__",inplace=True)
  total=df['Size'].sum()
  reads_abundance = df.groupby([taxo])['Size'].sum()
  abundance=reads_abundance/total
  df=pd.DataFrame({"reads":reads_abundance,"abundance":abundance})
  df.to_csv(out2, header=False,sep='\t')
