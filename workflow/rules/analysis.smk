import glob
import shutil
import numpy as np
import pandas as pd
import toolbox
import os
root = toolbox.set_root()
tmp = toolbox.set_tmp(root, config)

os.chdir(tmp)

Database_name=config["Database"].split(':')[0].upper()
#Database_name=config["targets_taxonomy"].upper()
dataset=config["Dataset"] 
extract = 'zcat' if config["format"] == 'gz' else 'cat'
toolname = config["container"].split("/")[-1].split(":")[0]
instances = config["instance"]

load_model=config["load_model"]
technology=config["technology"]
model=config["model"] if load_model != 0 else "none"

container_runner, docker_container, clean = toolbox.define_docker_command(
    toolname, "container", root, tmp, config, int(workflow.cores / config["distribute_cores_to_n_tasks"])
    if workflow.cores > 1
    else 1,
)
toolbox.validate_config(config)

#READS CONFIG
all_samples = config["sample"]
print(f"{dataset}+{' '.join(all_samples)}")
for sample in config["sample"]:
    dir_name = f"{root}/instances/{instances}/{sample}"
    os.makedirs(dir_name, exist_ok=True)

#TOOL config
aux_parameters=",".join(config["aux_parameters"]) if ("aux_parameters" in config and config["aux_parameters"] != None) else "none=none"
aux_parameters=aux_parameters.replace(" ", "")
print(aux_parameters)

source = Path(f"{root}/repository/{Database_name}/dna-sequences.fasta")
link = Path(f"{root}/instances/{instances}/reference.fasta")

if link.exists() or link.is_symlink():
    link.unlink()
link.symlink_to(source.resolve())


source = Path(f"{root}/repository/{Database_name}/taxonomy.tsv")
link = Path(f"{root}/instances/{instances}/reference.tsv")

if link.exists() or link.is_symlink():
    link.unlink()
link.symlink_to(source.resolve())

localrules:
    target,


rule target:
    input:
        predictions2=expand(
            "{root}analysis_outputs/{toolname}.{instance}-{e}.predictions.tsv",
            root=root,
            toolname=toolname,
            instance=instances,
            e=all_samples,
        ),
        runtime_analysis2=expand(
            "{root}analysis_outputs/{toolname}.{instance}-{e}.runtime_analysis.txt",
            root=root,
            toolname=toolname,
            instance=instances,
            e=all_samples,
        ),
        memory_analysis2=expand(
            "{root}analysis_outputs/{toolname}.{instance}-{e}.memory_analysis.txt",
            root=root,
            toolname=toolname,
            instance=instances,
            e=all_samples,
        ),    
        runtime_training=expand(
            "{root}analysis_outputs/{toolname}.{instance}.runtime_training.txt",
            root=root,
            toolname=toolname,
            instance=instances,
        ),
        memory_training=expand(
            "{root}analysis_outputs/{toolname}.{instance}.memory_training.txt",
            root=root,
            toolname=toolname,
            instance=instances,
        ),


rule cleanReads:
    input:
        data = lambda wc: f"{dataset}/{wc.all_samples}"
    output:
        cleanReads="{root}instances/{instance}/{all_samples}/queryReads.fq"
    params: docker_suffix="clean" if container_runner != "" else "", 
    message:
        "Running cleanReads for sample {wildcards.all_samples}..."
    container:
        "{}sif/{}.sif".format(root, toolname)
    threads: int(workflow.cores)   
    shell:
        """
       {extract} {input.data}/* > {output.cleanReads}_clean
        """

rule training_model:
    input:
        instance=expand("{root}instances/{instance}", root=root, instance=instances),
    output:
        runtime_training=expand(
            "{root}analysis_outputs/{toolname}.{{instance}}.runtime_training.txt",
            root=root,
            toolname=toolname,
        ),
        memory_training=expand(
            "{root}analysis_outputs/{toolname}.{{instance}}.memory_training.txt",
            root=root,
            toolname=toolname,
        ),
    params: docker_suffix="training" if container_runner != "" else "", 
    message:
        "run training..."
    container:
        "{}sif/{}.sif".format(root, toolname)
    threads: int(workflow.cores)   
    shell:
        """
        # run the tool        
        {container_runner}{params.docker_suffix} {docker_container} {root}../workflow/scripts/task_wrapper.sh \
        LEMMI16s_training.sh \
        {wildcards.instance} {input.instance} {root}repository/{Database_name} {load_model} {model} '{aux_parameters}' training_{wildcards.instance} 

        mv {tmp}memory.training_{wildcards.instance}.txt {output.memory_training}
        mv {tmp}runtime.training_{wildcards.instance}.txt {output.runtime_training}

        # remove stopped docker containers
        {clean}

        """
        
rule run_analysis:
    input:
        cleanReads=expand(
            "{root}instances/{instance}/{all_samples}/queryReads.fq",
            root=root,
            instance=instances,
            all_samples=all_samples,
        ),
        runtime_training=expand(
            "{root}analysis_outputs/{toolname}.{instance}.memory_training.txt",
            root=root,
            toolname=toolname,
            instance=instances,
        ),
        memory_training=expand(
            "{root}analysis_outputs/{toolname}.{instance}.memory_training.txt",
            root=root,
            toolname=toolname,
            instance=instances,
        )
    output:
        predictions=expand(
            "{root}analysis_outputs/{toolname}.{{instance}}-{{all_samples}}.predictions.tsv",
            root=root,
            toolname=toolname,
        ),
        runtime_analysis=expand(
            "{root}analysis_outputs/{toolname}.{{instance}}-{{all_samples}}.runtime_analysis.txt",
            root=root,
            toolname=toolname,
        ),
        memory_analysis=expand(
            "{root}analysis_outputs/{toolname}.{{instance}}-{{all_samples}}.memory_analysis.txt",
            root=root,
            toolname=toolname,
        ),
    params: docker_suffix=(lambda wildcards: wildcards.all_samples if container_runner != "" else ""),
            tech=technology if technology == "PacBio" else "illumina", 
    message:
        "run analysis on calibrarion and evaluation sets..."
    container:
        "{}sif/{}.sif".format(root, toolname)
    threads: int(workflow.cores)
    shell:
        """
        {container_runner}{params.docker_suffix} {docker_container} {root}../workflow/scripts/task_wrapper.sh \
        LEMMI16s_analysis.sh {params.tech}\
        {wildcards.instance}-{wildcards.all_samples} {input.instance}/{wildcards.all_samples} \
        {tmp}tmp_{wildcards.instance} {output.predictions} '{aux_parameters}' analysis_{wildcards.instance}-{wildcards.all_samples}

        mv {tmp}memory.analysis_{wildcards.instance}-{wildcards.all_samples}.txt {output.memory_analysis}
        mv {tmp}runtime.analysis_{wildcards.instance}-{wildcards.all_samples}.txt {output.runtime_analysis}

        # remove stopped docker containers
        {clean}
        """
        
        

