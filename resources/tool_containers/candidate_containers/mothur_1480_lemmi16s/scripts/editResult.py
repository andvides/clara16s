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
taxonomy = pd.read_csv(str(taxo_result), header=0, sep="\t")
otus = taxonomy["Taxonomy"].str.split(';',expand=True)
otus = otus.drop([7], axis=1)
otus.insert(0, "Size", taxonomy["Size"])
otus.insert(0, "OTU", taxonomy["OTU"])
otus.columns = ['Group_ID',"Size","domain","phylum","class","order","family","genus","species"]
#print(otus)

#Write results
out1.write("#LEMMI16s\n")
out1.write("Group_ID\tSize\tTaxonomy\n")
for index, row in otus.iterrows():
  toprint='group_'+str(index)+"\t"+str(row["Size"])+'\t'+ \
          row["domain"].split('(')[0]+';'+row["phylum"].split('(')[0]+';'+row["class"].split('(')[0]+';'+ \
          row["order"].split('(')[0]+';'+row["family"].split('(')[0]+';'+row["genus"].split('(')[0]+';'+row["species"].split('(')[0]+'\n'
  out1.write(toprint)

#Write summary taxonomy
out2.write("#LEMMI16s\norganism\treads\tabundance\n")
for taxo in ["domain","phylum","class","order","family","genus","species"]:
  df=otus[[taxo,'Size']]
  df[taxo]=df[taxo].str.replace("__",taxo[0]+"__").str[1:]
  total=df['Size'].sum()
  reads_abundance = df.groupby([taxo])['Size'].sum()
  abundance=reads_abundance/total
  df=pd.DataFrame({"reads":reads_abundance,"abundance":abundance})
  df.to_csv(out2, header=False,sep='\t')
