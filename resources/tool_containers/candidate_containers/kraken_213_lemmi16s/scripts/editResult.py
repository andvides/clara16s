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
otus = pd.read_csv(str(taxo_result), header=0, sep="\t")
otus.columns = ['abundance',"Size","Size2","Taxo","TaxId","name"]

#Write results
out1.write("#LEMMI16s\n")
out1.write("Group_ID\tSize\tTaxonomy\n")
taxo={'D':'d__na','P':'p__na','C':'c__na','O':'o__na','F':'f__na','G':'g__na','S':'s__na',}
i=1
for index, row in otus.iterrows():
  row["name"]=row["Taxo"].lower()+'__'+row["name"].replace(" ","")
  taxo[row["Taxo"]]=row["name"]
  if row["Taxo"]=='G':
    toprint='OTU_'+str(i)+"\t"+str(row["Size"])+'\t'+taxo['D']+';'+taxo['P']+';'+taxo['C']+';'+taxo['O']+';'+taxo['F']+";"+taxo['G']+";"+taxo['S']+"\n"
    out1.write(toprint)
    i+=1
  delete=False
  for key in taxo:
    if key==row["Taxo"]:
      delete=True
    elif delete:
      taxo[key]=key.lower()+'__na'
    

#Write summary taxonomy
out2.write("#LEMMI16s\norganism\treads\tabundance\n")
for taxo in ["D","P","C","O","F","G"]:
  df=otus.loc[otus['Taxo'] == taxo]
  df["name"]=taxo.lower()+'__'+df["name"].str.replace(" ","")
  df["abundance"]=df["abundance"]/100
  df=df[["name","Size","abundance"]]
  df.to_csv(out2, header=False,index=False,sep='\t')

#Reply genus to species
df=otus.loc[otus['Taxo'] == 'G']
df["name"]='s__'+df["name"].str.replace(" ","")
df["abundance"]=df["abundance"]/100
df=df[["name","Size","abundance"]]
df.to_csv(out2, header=False,index=False,sep='\t')

