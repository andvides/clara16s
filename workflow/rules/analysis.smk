import glob
import shutil
import numpy as np
import pandas as pd
import toolbox
import reads_layout
import os
import sys
root = toolbox.set_root()
tmp = toolbox.set_tmp(root, config)

os.chdir(tmp)

Database_name=config["Database"].split(':')[0].upper()
#Database_name=config["targets_taxonomy"].upper()
dataset=config["Dataset"] 
toolname = config["container"].split("/")[-1].split(":")[0]
instances = config["instance"]

load_model=config["load_model"]
model=config["model"] if load_model != 0 else "none"

# SEQUENCING TECHNOLOGY: illumina | nanopore | pacbio (case insensitive), same values as cleaning.smk
technology = str(config["technology"]).strip().lower()
if technology not in ("illumina", "nanopore", "pacbio"):
    sys.exit(f"ERROR: technology must be illumina, nanopore or pacbio (got '{config['technology']}')")

# LIBRARY LAYOUT (only used with illumina, optional): auto | paired | interleaved | single
# cleaning.smk delivers paired and interleaved reads as two files (forward/reverse) and single-end
# reads as one file, so for the analysis "interleaved" is the same as "paired".
library_layout = str(config.get("library_layout") or "auto").strip().lower()
if library_layout not in reads_layout.LAYOUTS:
    sys.exit(f"ERROR: library_layout must be one of {', '.join(reads_layout.LAYOUTS)} (got '{library_layout}')")

# UNCOMPRESS (optional, default false): how the trimmed reads are handed to the tool scripts
#   false -> symbolic links (queryReads*.fq) to the gzip files produced by cleaning.smk
#   true  -> uncompressed copies (plain-text queryReads*.fq), for tools that cannot read gzip files
#            (needs about 4-5 times the disk space of the compressed reads for every sample)
uncompress_value = config.get("uncompress")
if uncompress_value is None or str(uncompress_value).strip() == "":
    uncompress = False
elif isinstance(uncompress_value, bool):
    uncompress = uncompress_value
elif str(uncompress_value).strip().lower() in ("true", "yes", "y", "1"):
    uncompress = True
elif str(uncompress_value).strip().lower() in ("false", "no", "n", "0"):
    uncompress = False
else:
    sys.exit(f"ERROR: uncompress must be true or false (got '{uncompress_value}')")

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

# READS produced by cleaning.smk (trimmed reads of every sample)
def trimmed_reads(sample):
    folder = f"{root}instances/{instances}/{sample}/trimmed"
    return {"paired": [f"{folder}/queryReads1.trimmed.fastq.gz", f"{folder}/queryReads2.trimmed.fastq.gz"],
            "single": [f"{folder}/queryReads.trimmed.fastq.gz"]}


def read_structure(sample):
    """'paired' (forward + reverse files) or 'single' (one file) for the trimmed reads of a sample."""
    if technology != "illumina" or library_layout == "single":
        return "single"
    if library_layout in ("paired", "interleaved"):
        return "paired"
    # auto: use the files that cleaning.smk produced
    reads = trimmed_reads(sample)
    has_paired = all(os.path.exists(path) for path in reads["paired"])
    has_single = os.path.exists(reads["single"][0])
    if has_paired and not has_single:
        return "paired"
    if has_single and not has_paired:
        return "single"
    if has_paired and has_single:
        raise ValueError("both paired (queryReads1/2) and single (queryReads) trimmed reads exist; "
                         "set library_layout to paired or single")
    raise ValueError(f"no trimmed reads found in {os.path.dirname(reads['single'][0])}; "
                     "run the cleaning step first (or set library_layout)")


READ_STRUCTURE = {}
for sample in all_samples:
    try:
        READ_STRUCTURE[sample] = read_structure(sample)
    except ValueError as error:
        sys.exit(f"ERROR in sample {sample}: {error}")
    print(f"Sample {sample}: {technology} {READ_STRUCTURE[sample]} reads")


def sample_reads(wildcards):
    return trimmed_reads(wildcards.all_samples)[READ_STRUCTURE[wildcards.all_samples]]


def tool_technology(wildcards):
    # The tool scripts (LEMMI16s_analysis.sh) only compare this value with "illumina":
    #   illumina        -> paired reads: queryReads1.fq + queryReads2.fq
    #   anything else   -> one read file: queryReads.fq (long-read branch of the script)
    if technology == "illumina":
        return "illumina" if READ_STRUCTURE[wildcards.all_samples] == "paired" else "illumina_single"
    return {"pacbio": "PacBio", "nanopore": "Nanopore"}[technology]

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
        instance=expand("{root}instances/{instance}", root=root, instance=instances),
        reads=sample_reads,
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
            tech=tool_technology,
            structure=lambda wildcards: READ_STRUCTURE[wildcards.all_samples],
            uncompress="1" if uncompress else "0",
            reads_forward=lambda wildcards, input: input.reads[0],
            reads_reverse=lambda wildcards, input: input.reads[1] if len(input.reads) > 1 else "",
    message:
        "run analysis on calibrarion and evaluation sets..."
    container:
        "{}sif/{}.sif".format(root, toolname)
    threads: int(workflow.cores)
    shell:
        """
        # link the trimmed reads with the names expected by the tool scripts:
        #   paired reads (illumina paired / interleaved) -> queryReads1.fq + queryReads2.fq
        #   one file (illumina single, PacBio, Nanopore)  -> queryReads.fq
        # uncompress=false: symbolic links to the gzip files; uncompress=true: uncompressed copies
        # (the old files/links are removed first so a copy is never written through an old link)
        rm -f {input.instance}/{wildcards.all_samples}/queryReads.fq {input.instance}/{wildcards.all_samples}/queryReads1.fq {input.instance}/{wildcards.all_samples}/queryReads2.fq
        if [ "{params.structure}" = "paired" ]; then
            if [ "{params.uncompress}" = "1" ]; then
                gzip -dc {params.reads_forward} > {input.instance}/{wildcards.all_samples}/queryReads1.fq
                gzip -dc {params.reads_reverse} > {input.instance}/{wildcards.all_samples}/queryReads2.fq
            else
                ln -s {params.reads_forward} {input.instance}/{wildcards.all_samples}/queryReads1.fq
                ln -s {params.reads_reverse} {input.instance}/{wildcards.all_samples}/queryReads2.fq
            fi
        else
            if [ "{params.uncompress}" = "1" ]; then
                gzip -dc {params.reads_forward} > {input.instance}/{wildcards.all_samples}/queryReads.fq
            else
                ln -s {params.reads_forward} {input.instance}/{wildcards.all_samples}/queryReads.fq
            fi
        fi

        {container_runner}{params.docker_suffix} {docker_container} {root}../workflow/scripts/task_wrapper.sh \
        LEMMI16s_analysis.sh {params.tech} \
        {wildcards.instance}-{wildcards.all_samples} {input.instance}/{wildcards.all_samples} \
        {tmp}tmp_{wildcards.instance} {output.predictions} '{aux_parameters}' analysis_{wildcards.instance}-{wildcards.all_samples}

        mv {tmp}memory.analysis_{wildcards.instance}-{wildcards.all_samples}.txt {output.memory_analysis}
        mv {tmp}runtime.analysis_{wildcards.instance}-{wildcards.all_samples}.txt {output.runtime_analysis}

        # remove stopped docker containers
        {clean}

        """
        
        

