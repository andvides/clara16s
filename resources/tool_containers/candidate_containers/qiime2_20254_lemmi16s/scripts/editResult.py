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
otus = pd.read_csv(str(taxo_result), header=1, sep="\t")
otus.columns = ['Taxonomy',"Size"]


#Write results
out1.write("#LEMMI16s\n")
out1.write("Group_ID\tSize\tTaxonomy\n")
for index, row in otus.iterrows():
  size=int(row["Size"])
  if size>0:
    toprint='group_'+str(index)+"\t"+str(size)+'\t'+row["Taxonomy"].replace(';__',';na__na')+'\n'
    out1.write(toprint)

#Write summary taxonomy
out2.write("#LEMMI16s\norganism\treads\tabundance\n")
df2 = otus["Taxonomy"].str.split(';',expand=True)
df2["Size"]=otus["Size"].astype('int')
df2.columns = ["domain","phylum","class","order","family","genus","species","Size"]

for taxo in ["domain","phylum","class","order","family","genus","species"]:
  df=df2[[taxo,'Size']]
  df[taxo].replace("__",taxo[0]+"__",inplace=True)
  total=df['Size'].sum()
  reads_abundance = df.groupby([taxo])['Size'].sum()
  abundance=reads_abundance/total
  df=pd.DataFrame({"reads":reads_abundance,"abundance":abundance})
  df=df[~(df == 0).all(axis=1)]
  df.to_csv(out2, header=False,sep='\t')
