#!/usr/bin/env python3
import warnings
import sys
import pandas as pd
import numpy as np
import os.path
import glob
import subprocess

totalReads=int(sys.argv[1])
sd_lognormal=int(sys.argv[2])
amplicon_ratio=int(sys.argv[3])
train_rate=int(sys.argv[4])
random_seed=int(sys.argv[5])
directory=sys.argv[6]
taxonDist=open(sys.argv[7],'w+')
gnc_ref=sys.argv[8]


def mean_by_category(gnc_ref):
  gcn=pd.read_csv(str(gnc_ref), header=0, sep="\t")
  gcn = gcn[['rank','mean']]
  gcn["mean"]=gcn["mean"].round().astype('int')
  grouped = gcn.groupby('rank')
  mean_by_category = grouped['mean'].mean().round().astype('int')
  df_dict = mean_by_category.to_dict()
  df_dict['d']=df_dict['superkingdom']
  del df_dict['superkingdom']
  df_dict['p']=df_dict['phylum']
  del df_dict['phylum']
  df_dict['c']=df_dict['class']
  del df_dict['class']
  df_dict['o']=df_dict['order']
  del df_dict['order']
  df_dict['f']=df_dict['family']
  del df_dict['family']
  df_dict['g']=df_dict['genus']
  del df_dict['genus']
  df_dict['s']=df_dict['species']
  del df_dict['species']
  #print(df_dict)
  return df_dict

def GCN_reference(gnc_ref):
  gcn=pd.read_csv(str(gnc_ref), header=0, sep="\t")
  gcn = gcn[['taxid','name','mean']]
  gcn["mean"]=gcn["mean"].round().astype('int')
  gcn.columns = ['taxid','name','GCN_mean']
  #print(gcn)
  return gcn


def loadRecollectedSeqs(directory):
  files = glob.glob(directory+"/*.reg")
  taxonIDS={}
  for file in files:
    basename = os.path.basename(file)
    basename = basename.split('.reg')[0]
    if os.path.getsize(file) > 0: #to avoid empty files from previous stages
      numAmplicons=subprocess.check_output("grep -c '>' "+file, shell=True)
      numAmplicons=int(numAmplicons.decode("utf-8"))
      taxonIDS[basename]=numAmplicons
    else:
      toprint=basename+'\t0\t0\t0\t0\t0\t0\t0\t0\n'
      taxonDist.write(toprint)
  df=pd.DataFrame.from_dict(taxonIDS,orient='index')
  df.sort_index(inplace=True)
  return df

def genDistribution(num_samples,totalReads,amplicon_ratio,seed_base):
  #Draw samples from a log-normal distribution
  #num_samples=len(taxonIDS)
  np.random.seed(seed_base)
  p_abu = np.random.lognormal(1, sd_lognormal, num_samples)
  p_reads = [
            int(p)
            for p in p_abu
                     * totalReads
                     * amplicon_ratio
                     / sum(p_abu)
       ]
  #reads_per_sample = {taxonIDS[i]: p_reads[i] for i in range(len(taxonIDS))}
  return p_reads

toprint="name\ttaxid\tmarkers_in_DB\tmarkers_used\tGCN_mean\tamplicons\treads\tamp_coverage\tfull_coverage\n"
taxonDist.write(toprint)

gcn=GCN_reference(gnc_ref)

df = pd.DataFrame()  
df=loadRecollectedSeqs(directory)
df.reset_index(inplace=True)
df.columns=["id","markers in DB"]

mean_by_rank=mean_by_category(gnc_ref)

df['name'] = df['id'].str.split('__').str.get(1)
df['name'] = df['name'].str.replace('_',' ')
df = df.merge(gcn, how='left', on='name')
#df.fillna(1, inplace=True)
nan_indices = np.where(pd.isna(df))
nan_locations = list(zip(nan_indices[0], nan_indices[1]))
for row, col in nan_locations:
  df.iloc[row, 4]=mean_by_rank[df.iloc[row, 0][0]]

df["markers used"]= df["markers in DB"]*(train_rate/100)
df["markers used"]=df["markers used"].astype('int')

df["amplicons"]=df["markers used"]*df["GCN_mean"]
df['amplicons'] = df['amplicons'].round()
df["reads"]=genDistribution(len(df),totalReads,amplicon_ratio,random_seed)
df["coverage"]=df["reads"]/df["amplicons"]
df.replace([np.inf, -np.inf], np.nan, inplace=True)
df.fillna(0, inplace=True)

#df['coverage'] = df['coverage'].round()
#df["coverage"]=df["coverage"].astype('int')

df["Full coverage"]= df['coverage']*df["GCN_mean"]
df['Full coverage']= df['Full coverage'].round()
df["Full coverage"]= df["Full coverage"].astype('int')

#df["PCR"]=df["reads"]/(df["markers used"]*df["GCN_mean"])
#df["coverage"]=df["reads"]/df["markers used"]
#cols=["name","taxid","markers in DB","markers used","GCN_mean","amplicons","PCR","coverage","reads"]
cols=["id","taxid","markers in DB","markers used","GCN_mean","amplicons","reads","coverage","Full coverage"]
df = df[cols]
df = df.set_index('id')
df.to_csv(taxonDist, sep="\t", header=False, mode='a', index=True)


