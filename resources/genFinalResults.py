#!/usr/bin/env python3
import sys#!/usr/bin/env python3
import sys
import os.path
import pandas as pd
import numpy as np


inputFiles  = sys.argv[1].split()
outputFile  = sys.argv[2]
project=os.path.basename(outputFile)
project=project.split('.tsv')[0]

tool_list=[]
dataset_list=[]
taxonomy_list=[]
value_list=[]
for file in inputFiles:
  metric=os.path.basename(file).split('.')[0]
  tool = os.path.basename(file).split('.')[1]
  dataset = os.path.basename(file).split('.')[2]
  taxonomy = os.path.basename(file).split('.')[3]
  for line in open(file, 'r'):
      value=line.split('\n')[0]
      tool_list.append(tool)
      dataset_list.append(project.split('.')[0]+'-e-'+dataset)
      taxonomy_list.append(taxonomy)
      value_list.append(float(value))

d = {'toolname':tool_list, 'sample':dataset_list,'rank':taxonomy_list,'value':value_list}
df = pd.DataFrame(data=d)

sets=df.groupby(['toolname', 'rank']).groups
avg=df.groupby(['toolname', 'rank'])['value'].mean()
i=0
for group in sets:
   tool_list.append(group[0])
   dataset_list.append(project.split('.')[0])
   taxonomy_list.append(group[1])
   value_list.append(avg[i])
   i+=1

d = {'toolname':tool_list, 'sample':dataset_list,'rank':taxonomy_list,'value':value_list}
df2 = pd.DataFrame(data=d)
df2 = df2.reset_index(drop=True)
#df2 = df2.drop_duplicates()
df2.to_csv(outputFile, sep="\t")
outputFile=outputFile.split('tsv')[0]+'json'
df2.to_json(outputFile)
