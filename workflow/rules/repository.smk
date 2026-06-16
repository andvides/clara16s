#Snakemake file to download and prepare the ribosomal RNA gene database
#Sequences are removed from the downloaded dataset using Qiime2 and RESCRIPt plugin (qiime2 docker image is used).

import toolbox
import os

root = toolbox.set_root()
tmp = toolbox.set_tmp(root, config)
os.chdir(tmp)

#CLARA16s_master container used for basic process
container_runner, docker_container, clean = toolbox.define_docker_command(
    "clara16s_master", "clara16s_master",root, tmp, config, int(workflow.cores)
)


#Config option to indicate the reference database (see /config/config.yaml)
Database_name=config["Database"].split(':')[0].upper()
release=config["Database"].split(':')[-1]
version=release.split('_')[-1]

#Config option to indicate the version of rrnDB16S gene copy number (GCN) database
#https://rrndb.umms.med.umich.edu/
GCN_ref=config["GCN_ref"]

#Output files for this script
crudeSequences='raw_sequences'
if Database_name=='SILVA':
  crudeSequences='SILVA_'+version+'_SSURef_NR99_tax_silva_trunc'

#Final files
sequences = 'dna-sequences.fasta'
taxonomy =  'taxonomy.tsv'


rule all:
  input: expand("{root}repository/{Database_name}/{sequences}",root=root,Database_name=Database_name,sequences=sequences), 
         expand("{root}repository/{Database_name}/{taxonomy}",root=root,Database_name=Database_name,taxonomy=taxonomy),
         expand("{root}repository/{Database_name}/GCN_ref.tsv",root=root,Database_name=Database_name)

rule download_seq:
  output: expand("{root}repository/{Database_name}/{crudeSequences}.fasta",root=root,Database_name=Database_name,crudeSequences=crudeSequences)
  container:"{}sif/{}.sif".format(root, config["clara16s_master"].split("/")[-1].split(":")[0])
  params: container_suffix=("downloadDatabase" if container_runner != "" else "")
  shell:  """
          {container_runner}{params.container_suffix} {docker_container} downloadDatabase.sh {root}repository/{Database_name} {release} {crudeSequences}
          {clean}
          """

rule curation:
  input:  expand("{root}repository/{Database_name}/{crudeSequences}.fasta",root=root,Database_name=Database_name,crudeSequences=crudeSequences)
  container:"{}sif/{}.sif".format(root, config["clara16s_master"].split("/")[-1].split(":")[0])
  output: seq=expand("{root}repository/{Database_name}/{sequences}",root=root,Database_name=Database_name,sequences=sequences), 
          tax=expand("{root}repository/{Database_name}/{taxonomy}",root=root,Database_name=Database_name,taxonomy=taxonomy)
  params: container_suffix=("curation" if container_runner != "" else "")
  shell:  """
          {container_runner}{params.container_suffix}  {docker_container} cleanDatabase.sh {root}repository/{Database_name} {release} {output.seq} {output.tax}  
          {clean}
          """

rule download_GCN_ref:
  output: expand("{root}repository/{Database_name}/GCN_ref.tsv",root=root,Database_name=Database_name)
  shell: "wget --no-check-certificate -q -c {GCN_ref} --output-document={output}.zip && zcat {output}.zip >  {output}"

