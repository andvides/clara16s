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
Taxonomy	TaxonomyLevel	Label	Counts
0	0	Bacteria	9954
0	1	Bacteria;Proteobacteria	9954
0	2	Bacteria;Proteobacteria;Gammaproteobacteria	9954
"""
taxonomy = pd.read_csv(str(taxo_result), header=0, sep="\t")
taxonomy["Label"]=taxonomy["Label"].str.replace(' ','_')
otus = taxonomy["Label"].str.split(';',expand=True)
otus.insert(0, "Size", taxonomy["Counts"])
#otus.insert(0, "OTU", taxonomy["OTU"])
otus.columns = ["Size","domain","phylum","class","order","family","genus","species"]
#print(otus)

#Write results
out1.write("#LEMMI16s\n")
out1.write("Group_ID\tSize\tTaxonomy\n")
for index, row in otus.iterrows():
  toprint='group_'+str(index)+"\t"+str(row["Size"])+'\td__'+ \
          row["domain"]+';p__'+row["phylum"]+';c__'+row["class"]+';o__'+ \
          row["order"]+';f__'+row["family"]+';g__'+row["genus"]+';s__'+row["species"]+'\n'
  out1.write(toprint)

#Write summary taxonomy
out2.write("#LEMMI16s\norganism\treads\tabundance\n")
for taxo in ["domain","phylum","class","order","family","genus","species"]:
  df=otus[[taxo,'Size']]
  df[taxo]=taxo[0]+"__"+df[taxo]
  total=df['Size'].sum()
  reads_abundance = df.groupby([taxo])['Size'].sum()
  abundance=reads_abundance/total
  df=pd.DataFrame({"reads":reads_abundance,"abundance":abundance})
  df.to_csv(out2, header=False,sep='\t')

