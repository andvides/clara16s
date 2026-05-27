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
train_rate= config.get("QR_ratio", 20)  
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

#mhass Command Options
listParams = [' --' + seqParams for seqParams in config["sequencerParams"]]  if config["sequencerParams"] else ''
sequencerParams = " ".join(s.strip() for s in listParams)
print(f'sequencerParams:{sequencerParams}')
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
         base_query_fasta=expand("{root}{OUT_DIR}/{sample}/base_query.fasta",OUT_DIR=OUT_DIR,root=root,sample=sample),
         base_query_taxo=expand("{root}{OUT_DIR}/{sample}/base_queryTaxo.tsv",OUT_DIR=OUT_DIR,root=root,sample=sample),
         query_reads=expand("{root}{OUT_DIR}/{all_sample}/queryReads.fq",OUT_DIR=OUT_DIR,root=root,all_sample=all_sample),
         query_taxo=expand("{root}{OUT_DIR}/{all_sample}/queryTaxo.tsv",OUT_DIR=OUT_DIR,root=root,all_sample=all_sample),

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
        if [ -z "$(echo "{params.reg}" | xargs)" ]; then
           echo "Skipping selecReg"
           ln -sf {root}{OUT_DIR}/{sample}/{wildcards.taxonIDS}.fna {output} 
        else
           {container_runner}{params.container_suffix} {docker_container} extractRegion.py {root}{OUT_DIR}/{sample}/{wildcards.taxonIDS}.fna {output} '{params.reg}'
        fi 
         """

rule splitDataset:
  input:  taxon_reg=expand("{root}{OUT_DIR}/{sample}/{taxonIDS}.reg",OUT_DIR=OUT_DIR,root=root,taxonIDS=taxonIDS,sample=sample),
  output: taxon_query="{root}{OUT_DIR}/{sample}/{taxonIDS}.reg.query", taxon_ref="{root}{OUT_DIR}/{sample}/{taxonIDS}.reg.ref"
  container:"{}sif/{}.sif".format(root, config["lemmi16s_master"].split("/")[-1].split(":")[0])  
  params: train_rate=train_rate, rs=random_seed,
          container_suffix=(lambda wildcards: re.sub(r'[\W_]', '',wildcards.taxonIDS.strip()+'splitDataset') if container_runner != "" else "")
  threads: 1
  shell: """
         {container_runner}{params.container_suffix} {docker_container} splitDataset.sh {root}{OUT_DIR}/{sample}/{wildcards.taxonIDS}.reg {params.train_rate} '{params.rs}'
         """

rule join:
  input:  taxon_query=expand("{root}{OUT_DIR}/{sample}/{taxonIDS}.reg.query",OUT_DIR=OUT_DIR,root=root,taxonIDS=taxonIDS,sample=sample)
  output: base_query_fasta="{root}{OUT_DIR}/{sample}/base_query.fasta"
  threads: 1  
  shell: " cat {root}{OUT_DIR}/{sample}/*.reg.query > {output.base_query_fasta} && \
           sed -E 's/^(>[^ ]+).*/\\1/' {output.base_query_fasta} > {output.base_query_fasta}.tmp && \
           mv {output.base_query_fasta}.tmp {output.base_query_fasta}"

rule writeTaxo:
  input:  base_query_fasta=expand("{root}{OUT_DIR}/{sample}/base_query.fasta",OUT_DIR=OUT_DIR,root=root,sample=sample)
  output: base_query_taxo="{root}{OUT_DIR}/{sample}/base_queryTaxo.tsv",
  container: "{}sif/{}.sif".format(root, config["lemmi16s_master"].split("/")[-1].split(":")[0]) 
  params: container_suffix=("writeTaxo" if container_runner != "" else ""),
          header="asvid\tgenomeid"
  threads: 1        
  shell: """
      {container_runner}{params.container_suffix} {docker_container} writeRefTaxo.py {input.base_query_fasta} {root}repository/{taxonomy}  {output.base_query_taxo} '{params.header}'
      """             

rule genReads:
  input:  base_query_fasta=expand("{root}{OUT_DIR}/{sample}/base_query.fasta",OUT_DIR=OUT_DIR,root=root,sample=sample),
          base_query_taxo=expand("{root}{OUT_DIR}/{sample}/base_queryTaxo.tsv",OUT_DIR=OUT_DIR,root=root,sample=sample),
  output: query_read="{root}{OUT_DIR}/{all_sample}/queryReads.fq",
          query_taxo="{root}{OUT_DIR}/{all_sample}/queryTaxo.tsv",
  container:"{}sif/{}.sif".format(root, config["lemmi16s_master"].split("/")[-1].split(":")[0])          
  params: read_length_distribution=read_length_distribution, cp=chimera, rs=random_seed, tr=totalReads, sp=sequencerParams,
          container_suffix=(lambda wildcards: re.sub(r'[\W_]', '',wildcards.all_sample.strip()) if container_runner != "" else "")
  shell: "{container_runner}{params.container_suffix} {docker_container} \
        genLongReads.sh '{root}{OUT_DIR}/{sample}/base_query.fasta' '{root}{OUT_DIR}/{sample}/base_queryTaxo.tsv' '{root}{OUT_DIR}/{wildcards.all_sample}'  {params.tr} {root}repository/{taxonomy} {workflow.cores} '{params.sp}'"  
#genReads.sh {root}{OUT_DIR}/{wildcards.all_sample}/filterTaxon.dist '{params.read_length_distribution}' {params.cp} {params.rs} {root}{OUT_DIR}/{sample}/{wildcards.taxonIDS}.reg.query {output.taxonForward} {output.taxonReverse}

rule genRef:
  input:  taxon_ref=expand("{root}{OUT_DIR}/{sample}/{taxonIDS}.reg.ref",OUT_DIR=OUT_DIR,root=root,taxonIDS=taxonIDS,sample=sample),
  output: refTaxo="{root}{OUT_DIR}/reference.tsv", ref="{root}{OUT_DIR}/reference.fasta"
  container:"{}sif/{}.sif".format(root, config["lemmi16s_master"].split("/")[-1].split(":")[0])  
  params: container_suffix=("genRef" if container_runner != "" else "")
  threads: 1  
  shell: """
         {container_runner}{params.container_suffix} {docker_container} writeRefTaxo.sh '{taxonIDS_unknown}' {root}{OUT_DIR} {root}repository/{taxonomy} {root}repository/{sequences} {output.ref} {output.refTaxo}
         {clean}
         """
