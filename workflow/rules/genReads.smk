import toolbox
import os
import re

root = toolbox.set_root()
tmp = toolbox.set_tmp(root, config)

container_runner, docker_container, clean = toolbox.define_docker_command(
    "LEMMI16s_master", "lemmi16s_master",root, tmp, config, int(workflow.cores)
)

#DATABASE
Database_name=config["Database"].split(':')[0]
sequences = Database_name+'/dna-sequences.fasta'
taxonomy = Database_name+'/taxonomy.tsv'
  
#DIRECTORIES
OUT_DIR="datasets/"+config["dataset"]
OUT_DIR_CALIBRATION = "datasets/"+config["dataset"]+"/calibration"
OUT_DIR_EVALUATION  = "datasets/"+config["dataset"]+"/evaluation"

isExist = os.path.exists(root+OUT_DIR_CALIBRATION)
if not isExist:
  os.makedirs(root+OUT_DIR_CALIBRATION)

isExist = os.path.exists(root+OUT_DIR_EVALUATION)
if not isExist:
  os.makedirs(root+OUT_DIR_EVALUATION)


#READS CONFIG
taxonIDS = config["taxon"]
regIDS = config["region"]
chimera= config["chimera"]
totalReads= int(config["totalReads"]/len(taxonIDS))
random_seed= config["random_seed"]
read_length= config["read_length"]
read_mean= config["read_mean"]
read_std= config["read_std"]
read_seqSys= config["read_seqSys"]
read_length_distribution='-ss '+read_seqSys+' -l '+str(read_length)+' -m '+str(read_mean)+' -s '+str(read_std) #because some tools takes them in the same parameter

#sd_lognormal=2
#amplicon_ratio=1


rule all:
  input: tax_calibration=expand("{root}{OUT_DIR_CALIBRATION}/calibration.queryTaxo.tsv",OUT_DIR_CALIBRATION=OUT_DIR_CALIBRATION,root=root),
         tax_evaluation =expand("{root}{OUT_DIR_EVALUATION}/evaluation.queryTaxo.tsv",OUT_DIR_EVALUATION=OUT_DIR_EVALUATION,root=root)

#CALIBRATION reads
rule calibration_genReads:
  input:  taxon_dist=expand("{root}{OUT_DIR}/filterTaxon.dist", root=root, OUT_DIR=OUT_DIR),
          taxon_query=expand("{root}{OUT_DIR}/{taxonIDS}.reg.query",root=root, OUT_DIR=OUT_DIR,taxonIDS=taxonIDS)
  output: seqF="{root}{OUT_DIR_CALIBRATION}/{taxonIDS}.calibration.query1.fq", 
          seqR="{root}{OUT_DIR_CALIBRATION}/{taxonIDS}.calibration.query2.fq"
  params: read_length_distribution=read_length_distribution, cp=chimera, rs=random_seed,  
          container_suffix=(lambda wildcards: re.sub(r'[\W_]', '',wildcards.taxonIDS.strip()) if container_runner != "" else "")
  shell: "{container_runner}{params.container_suffix}_calibration --rm {docker_container} \
          genReads.sh {input.taxon_dist} '{params.read_length_distribution}' {params.cp} {params.rs} {root}{OUT_DIR}/{wildcards.taxonIDS}.reg.query  {root}{OUT_DIR_CALIBRATION}/{wildcards.taxonIDS}.calibration.query"

rule calibration_join:
  input:  tid_fq1=expand("{root}{OUT_DIR_CALIBRATION}/{taxonIDS}.calibration.query1.fq",OUT_DIR_CALIBRATION=OUT_DIR_CALIBRATION,root=root,taxonIDS=taxonIDS),
          tid_fq2=expand("{root}{OUT_DIR_CALIBRATION}/{taxonIDS}.calibration.query2.fq",OUT_DIR_CALIBRATION=OUT_DIR_CALIBRATION,root=root,taxonIDS=taxonIDS)
  output: tax=expand("{root}{OUT_DIR_CALIBRATION}/calibration.queryTaxo.tsv",OUT_DIR_CALIBRATION=OUT_DIR_CALIBRATION,root=root)
  params: container_suffix=("join_calibration" if container_runner != "" else "")
  shell: "cat {root}{OUT_DIR_CALIBRATION}/*.calibration.query1.fq > {root}{OUT_DIR_CALIBRATION}/queryReads1.fq | \
          cat {root}{OUT_DIR_CALIBRATION}/*.calibration.query2.fq > {root}{OUT_DIR_CALIBRATION}/queryReads2.fq | \
          {container_runner}{params.container_suffix}2 --rm {docker_container} writeQueryTaxo.py {root}{OUT_DIR_CALIBRATION}/queryReads1.fq {root}repository/{taxonomy} {output.tax}"

#EVALUATION  reads
rule evaluation_genReads:
  input:  taxon_dist=expand("{root}{OUT_DIR}/filterTaxon.dist", root=root, OUT_DIR=OUT_DIR),
          taxon_query=expand("{root}{OUT_DIR}/{taxonIDS}.reg.query",root=root, OUT_DIR=OUT_DIR,taxonIDS=taxonIDS)
  output: seqF="{root}{OUT_DIR_EVALUATION}/{taxonIDS}.evaluation.query1.fq", 
          seqR="{root}{OUT_DIR_EVALUATION}/{taxonIDS}.evaluation.query2.fq"
  params: read_length_distribution=read_length_distribution, cp=chimera, rs=random_seed+1,  
          container_suffix=(lambda wildcards: re.sub(r'[\W_]', '',wildcards.taxonIDS.strip()) if container_runner != "" else "")
  shell: "{container_runner}{params.container_suffix}_evaluation --rm {docker_container} \
          genReads.sh {input.taxon_dist} '{params.read_length_distribution}' {params.cp} {params.rs} {root}{OUT_DIR}/{wildcards.taxonIDS}.reg.query  {root}{OUT_DIR_EVALUATION}/{wildcards.taxonIDS}.evaluation.query"

rule evaluation_join:
  input:  tid_fq1=expand("{root}{OUT_DIR_EVALUATION}/{taxonIDS}.evaluation.query1.fq",OUT_DIR_EVALUATION=OUT_DIR_EVALUATION,root=root,taxonIDS=taxonIDS),
          tid_fq2=expand("{root}{OUT_DIR_EVALUATION}/{taxonIDS}.evaluation.query2.fq",OUT_DIR_EVALUATION=OUT_DIR_EVALUATION,root=root,taxonIDS=taxonIDS)
  output: tax=expand("{root}{OUT_DIR_EVALUATION}/evaluation.queryTaxo.tsv",OUT_DIR_EVALUATION=OUT_DIR_EVALUATION,root=root)
  params: container_suffix=("join_evaluation" if container_runner != "" else "")
  shell: "cat {root}{OUT_DIR_EVALUATION}/*.evaluation.query1.fq > {root}{OUT_DIR_EVALUATION}/queryReads1.fq | \
          cat {root}{OUT_DIR_EVALUATION}/*.evaluation.query2.fq > {root}{OUT_DIR_EVALUATION}/queryReads2.fq | \
          {container_runner}{params.container_suffix}2 --rm {docker_container} writeQueryTaxo.py {root}{OUT_DIR_EVALUATION}/queryReads1.fq {root}repository/{taxonomy} {output.tax}"


