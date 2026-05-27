#Snakemake file to create the instances
#Use the yam/instances/*.yaml config file to set the parameters

import toolbox
import os
import re
import numpy as np

root = toolbox.set_root()
tmp = toolbox.set_tmp(root, config)
os.chdir(tmp)

container_runner, docker_container, clean = toolbox.define_docker_command(
    "LEMMI16s_master", "lemmi16s_master",root, tmp, config, int(workflow.cores / config["distribute_cores_to_n_tasks"])
    if workflow.cores > 1
    else 1,
)

#DATABASE parameters (the same names should be used in repository.smk)
#Database_name=config["Database"].split(':')[0]
Database_name=config["targets_taxonomy"].upper()
sequences = Database_name+'/dna-sequences.fasta' 
taxonomy  = Database_name+'/taxonomy.tsv'
GCN_ref   = Database_name+'/GCN_ref.tsv'

#READS CONFIG options
taxonIDS = [s.split(":")[0] for s in config["taxon"]] 
taxonIDS_unknown = [s.split(":")[0] if s.split(":")[-1]=='unknown' else '' for s in config["taxon"] ] 
#print(taxonIDS_unknown)
    
#regIDS = config["region"]
totalReads= config["totalReads"]
train_rate=config["QR_ratio"]  if config["QR_ratio"] else 20
sd_lognormal=2
amplicon_ratio=1
regions = config.get("region", [])
regIDS = ['--region ' + r for r in regions] if regions else [' ']

chimera= config.get("chimera", 0)
random_seed= config.get("random_seed", 1)
read_length= config.get("read_length", 200) 
read_mean= config.get("read_mean", 250) 
read_std= config.get("read_std", 2)
read_seqSys= config.get("read_seqSys", 'MSv3') 
max_base_quality=config.get("max_base_quality", 35) 
min_base_quality=config.get("min_base_quality", 10) 
#Pass the reads parameters as a string because some tools takes them in the same parameter
read_length_distribution='-ss '+read_seqSys+' -l '+str(read_length)+' -m '+str(read_mean)+' -s '+str(read_std)+' -qL '+str(min_base_quality)+' -qU '+str(max_base_quality)

#Number of calibration and evaluation samples to generate
calibration_samples = ["c{}".format(s) for s in config["calibration"]]
evaluation_samples  = ["e{}".format(s) for s in config["evaluation"]]
all_sample=[
             "{}-{}".format(config["instance"], s)
             for s in evaluation_samples + calibration_samples
            ]
            
#Using the same seed to reproduce the results
np.random.seed(random_seed)
random_samples =np.random.choice(range(len(all_sample)), len(all_sample), replace=False)
random_instance = {all_sample[i].split("-")[-1]: random_samples[i] for i in range(len(all_sample))}

#Output files and directories
sample=config["instance"]
OUT_DIR = "instances/"+sample


rule all:
  input: refTaxo=expand("{root}{OUT_DIR}/reference.tsv",OUT_DIR=OUT_DIR,root=root),
         ref=expand("{root}{OUT_DIR}/reference.fasta",OUT_DIR=OUT_DIR,root=root),
         sample_dis=expand("{root}{OUT_DIR}/{all_sample}/filterTaxon.dist",OUT_DIR=OUT_DIR,root=root,all_sample=all_sample),
         sample_tax=expand("{root}{OUT_DIR}/{all_sample}/queryTaxo.tsv",OUT_DIR=OUT_DIR,root=root,all_sample=all_sample),
         query_forward=expand("{root}{OUT_DIR}/{all_sample}/queryReads.fq",OUT_DIR=OUT_DIR,root=root,all_sample=all_sample),
         query_reverse=expand("{root}{OUT_DIR}/{all_sample}/queryTaxo.tsv",OUT_DIR=OUT_DIR,root=root,all_sample=all_sample),

rule selecTax:  
  input:  taxo=expand("{root}repository/{taxonomy}",root=root,taxonomy=taxonomy), dbs=expand("{root}repository/{sequences}",root=root,sequences=sequences)
  output: list="{root}{OUT_DIR}/{sample}/{taxonIDS}.list", seqs="{root}{OUT_DIR}/{sample}/{taxonIDS}.fna"
  container:"{}sif/{}.sif".format(root, config["lemmi16s_master"].split("/")[-1].split(":")[0])
  params: container_suffix=(lambda wildcards: re.sub(r'[\W_]', '',wildcards.taxonIDS.strip()+'selecTax') if container_runner != "" else "")
  shell: """
         {container_runner}{params.container_suffix} {docker_container} filterTaxon.sh '{wildcards.taxonIDS}' {input.taxo} {output.list} {input.dbs} {output.seqs} 
         """

rule selecReg:
  input:  "{root}{OUT_DIR}/{sample}/{taxonIDS}.fna"
  output: "{root}{OUT_DIR}/{sample}/{taxonIDS}.reg"
  container:"{}sif/{}.sif".format(root, config["lemmi16s_master"].split("/")[-1].split(":")[0])
  params: reg=regIDS, container_suffix=(lambda wildcards: re.sub(r'[\W_]', '',wildcards.taxonIDS.strip()+'selecReg') if container_runner != "" else "")
  shell: """
         {container_runner}{params.container_suffix} {docker_container} extractRegion.py {root}{OUT_DIR}/{sample}/{wildcards.taxonIDS}.fna {output} '{params.reg}'
         """

rule splitDataset:
  input:  taxon_reg=expand("{root}{OUT_DIR}/{sample}/{taxonIDS}.reg",OUT_DIR=OUT_DIR,root=root,taxonIDS=taxonIDS,sample=sample),
  output: taxon_query="{root}{OUT_DIR}/{sample}/{taxonIDS}.reg.query", taxon_ref="{root}{OUT_DIR}/{sample}/{taxonIDS}.reg.ref"
  container:"{}sif/{}.sif".format(root, config["lemmi16s_master"].split("/")[-1].split(":")[0])  
  params: train_rate=train_rate, rs=random_seed,
          container_suffix=(lambda wildcards: re.sub(r'[\W_]', '',wildcards.taxonIDS.strip()+'splitDataset') if container_runner != "" else "")
  threads: 1
  shell: """
         {container_runner}{params.container_suffix} {docker_container} splitDataset.sh {root}{OUT_DIR}/{sample}/{wildcards.taxonIDS}.reg {params.train_rate} {params.rs}
         """
         
rule do_sampling:
  input:  taxon_query=expand("{root}{OUT_DIR}/{sample}/{taxonIDS}.reg.query",OUT_DIR=OUT_DIR,root=root,taxonIDS=taxonIDS,sample=sample)
  output: sample_dis="{root}{OUT_DIR}/{all_sample}/filterTaxon.dist"
  container:"{}sif/{}.sif".format(root, config["lemmi16s_master"].split("/")[-1].split(":")[0])
  params: tr=totalReads, sd_lognormal=sd_lognormal, amplicon_ratio=amplicon_ratio, train_rate=train_rate, 
          rs=lambda wildcards: [random_instance[wildcards.all_sample.split("-")[-1]]],
          container_suffix=(lambda wildcards: re.sub(r'[\W_]', '',wildcards.all_sample.strip()+'do_sampling') if container_runner != "" else "")
  shell: """
         {container_runner}{params.container_suffix} {docker_container} \
         computeDistribution.py {params.tr} {params.sd_lognormal} {params.amplicon_ratio} {params.train_rate} {params.rs} \
                                {root}{OUT_DIR}/{sample} {output.sample_dis} {root}repository/{GCN_ref}
         """

rule genReads:
  input:  sample_dis=expand("{root}{OUT_DIR}/{all_sample}/filterTaxon.dist",OUT_DIR=OUT_DIR,root=root,all_sample=all_sample),
          taxon_query=expand("{root}{OUT_DIR}/{sample}/{taxonIDS}.reg.query",OUT_DIR=OUT_DIR,root=root,taxonIDS=taxonIDS,sample=sample)
  output: taxonForward="{root}{OUT_DIR}/{all_sample}/{taxonIDS}.sim1.fq",
          taxonReverse="{root}{OUT_DIR}/{all_sample}/{taxonIDS}.sim2.fq"
  container:"{}sif/{}.sif".format(root, config["lemmi16s_master"].split("/")[-1].split(":")[0])          
  params: read_length_distribution=read_length_distribution, cp=chimera, rs=random_seed,  
          container_suffix=(lambda wildcards: re.sub(r'[\W_]', '',wildcards.taxonIDS.strip()+wildcards.all_sample.strip()) if container_runner != "" else "")
  shell: "{container_runner}{params.container_suffix} {docker_container} \
          genReads.sh {root}{OUT_DIR}/{wildcards.all_sample}/filterTaxon.dist '{params.read_length_distribution}' {params.cp} {params.rs} {root}{OUT_DIR}/{sample}/{wildcards.taxonIDS}.reg.query {output.taxonForward} {output.taxonReverse}"


rule join:
  input:  taxonForward=expand("{root}{OUT_DIR}/{all_sample}/{taxonIDS}.sim1.fq",OUT_DIR=OUT_DIR,root=root,all_sample=all_sample,taxonIDS=taxonIDS),
          taxonReverse=expand("{root}{OUT_DIR}/{all_sample}/{taxonIDS}.sim2.fq",OUT_DIR=OUT_DIR,root=root,all_sample=all_sample,taxonIDS=taxonIDS),
  output: query_forward="{root}{OUT_DIR}/{all_sample}/queryReads1.fq", query_reverse="{root}{OUT_DIR}/{all_sample}/queryReads2.fq"
  container:"{}sif/{}.sif".format(root, config["lemmi16s_master"].split("/")[-1].split(":")[0])  
  params: container_suffix=(lambda wildcards: re.sub(r'[\W_]', '',wildcards.all_sample.strip()+'join') if container_runner != "" else "")
  threads: 1  
  shell: "cat {root}{OUT_DIR}/{wildcards.all_sample}/*.sim1.fq > {output.query_forward} | \
          cat {root}{OUT_DIR}/{wildcards.all_sample}/*.sim2.fq > {output.query_reverse} "             

rule writeTaxo:
  input:  query_forward=expand("{root}{OUT_DIR}/{all_sample}/queryReads1.fq",OUT_DIR=OUT_DIR,root=root,all_sample=all_sample,taxonIDS=taxonIDS),
          query_reverse=expand("{root}{OUT_DIR}/{all_sample}/queryReads2.fq",OUT_DIR=OUT_DIR,root=root,all_sample=all_sample,taxonIDS=taxonIDS),
  output: sample_tax="{root}{OUT_DIR}/{all_sample}/queryTaxo.tsv"
  container:"{}sif/{}.sif".format(root, config["lemmi16s_master"].split("/")[-1].split(":")[0])  
  params: container_suffix=(lambda wildcards: re.sub(r'[\W_]', '',wildcards.all_sample.strip()+'write_calibrarion') if container_runner != "" else "")
  shell: "{container_runner}{params.container_suffix} {docker_container} writeQueryTaxo.py {root}{OUT_DIR}/{wildcards.all_sample}/queryReads1.fq {root}repository/{taxonomy} {output.sample_tax}"             


rule genRef:
  input:  taxon_ref=expand("{root}{OUT_DIR}/{sample}/{taxonIDS}.reg.ref",OUT_DIR=OUT_DIR,root=root,taxonIDS=taxonIDS,sample=sample),
          sample_tax=expand("{root}{OUT_DIR}/{all_sample}/queryTaxo.tsv",OUT_DIR=OUT_DIR,root=root,all_sample=all_sample),
  output: refTaxo="{root}{OUT_DIR}/reference.tsv", ref="{root}{OUT_DIR}/reference.fasta"
  container:"{}sif/{}.sif".format(root, config["lemmi16s_master"].split("/")[-1].split(":")[0])  
  params: container_suffix=("genRef" if container_runner != "" else "")
  threads: 1  
  shell: """
         {container_runner}{params.container_suffix} {docker_container} writeRefTaxo.sh '{taxonIDS_unknown}' {root}{OUT_DIR} {root}repository/{taxonomy} {root}repository/{sequences} {output.ref} {output.refTaxo}
         {clean}
         """
