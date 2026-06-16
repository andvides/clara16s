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
dataset=config["Dataset"] 
extract = 'zcat' if config["format"] == 'gz' else 'cat'
instances = config["instance"]
technology=config["technology"]

#CLARA16s_master container used for basic process
container_runner, docker_container, clean = toolbox.define_docker_command(
    "clara16s_master", "clara16s_master",root, tmp, config, int(workflow.cores)
)
toolbox.validate_config(config)

#READS CONFIG
all_samples = config["sample"]
print(f"{dataset}+{' '.join(all_samples)}")
for sample in config["sample"]:
    dir_name = f"{root}/instances/{instances}/{sample}"
    os.makedirs(dir_name, exist_ok=True)


localrules:
    target,

rule target:
    input:
        qc_metrics_raw=expand("{root}instances/{instance}/{all_samples}/raw/qc_metrics.tsv",
            root=root,
            instance=instances,
            all_samples=all_samples,),
        qc_summary_raw=expand("{root}instances/{instance}/{all_samples}/raw/qc_summary.tsv",
            root=root,
            instance=instances,
            all_samples=all_samples,),
        merged_reads=expand("{root}instances/{instance}/{all_samples}/merged/queryReads.concat.fastq.gz",
            root=root,
            instance=instances,
            all_samples=all_samples,),
        trimmed_reads=expand("{root}instances/{instance}/{all_samples}/trimmed/queryReads.trimmed.fastq.gz",
            root=root,
            instance=instances,
            all_samples=all_samples,),
        json=expand("{root}instances/{instance}/{all_samples}/trimmed/cutadapt.json",
            root=root,
            instance=instances,
            all_samples=all_samples,),
        log=expand("{root}instances/{instance}/{all_samples}/trimmed/cutadapt.log",
            root=root,
            instance=instances,
            all_samples=all_samples,),
        qc_metrics=expand("{root}instances/{instance}/{all_samples}/trimmed/qc_metrics.tsv",
            root=root,
            instance=instances,
            all_samples=all_samples,),
        qc_summary=expand("{root}instances/{instance}/{all_samples}/trimmed/qc_summary.tsv",
            root=root,
            instance=instances,
            all_samples=all_samples,),
        report=expand("{root}instances/{instance}/multiqc/multiqc_report.html",
            root=root,
            instance=instances,
            all_samples=all_samples,),
        
   
rule amplicon_qc_raw:
    input:
        reads=lambda wc: f"{dataset}/{wc.all_samples}"
    output:
        qc_metrics="{root}instances/{instance}/{all_samples}/raw/qc_metrics.tsv",
        qc_summary="{root}instances/{instance}/{all_samples}/raw/qc_summary.tsv"
    message:
        "Raw QC: {wildcards.all_samples}"
    params:docker_suffix=lambda wc: f"amplicon_qc_raw{wc.all_samples}" if container_runner != "" else "",
    shell:
        """
        {container_runner}{params.container_suffix} {docker_container} 
        ./fastqc_qc.sh \
            -i {input.reads} \
            -t {technology} \
            -o {root}instances/{instance}/{wildcards.all_samples}/raw
        {clean}
        """

rule mergeReads:
    input:
        reads=lambda wc: f"{dataset}/{wc.all_samples}",
        qc=rules.amplicon_qc_raw.output.qc_metrics
    output:
        reads="{root}instances/{instance}/{all_samples}/merged/queryReads.concat.fastq.gz"
    message:
        "Merging reads: {wildcards.all_samples}"
    params: docker_suffix=lambda wc: f"mergeReads{wc.all_samples}" if container_runner != "" else "",
    container:
        "{}sif/{}.sif".format(root, toolname)
    threads: int(workflow.cores)
    shell:
        """
        mkdir -p $(dirname {output.reads})
        cat {input.reads}/*.fastq.gz > {output.reads}
        """

rule trim_primers:
    input:
        reads=rules.mergeReads.output.reads
    output:
        reads="{root}instances/{instance}/{all_samples}/trimmed/queryReads.trimmed.fastq.gz",
        json="{root}instances/{instance}/{all_samples}/trimmed/cutadapt.json",
        log="{root}instances/{instance}/{all_samples}/trimmed/cutadapt.log"
    params:docker_suffix=lambda wc: f"trim_primers{wc.all_samples}" if container_runner != "" else "",
            tech=technology if technology == "PacBio" else "illumina", 
            fwd=forward_primer,
            rev=reverse_primer,
            min_len=100,
    message:
        "Cutadapt: {wildcards.all_samples}"
    container:
        "{}sif/{}.sif".format(root, toolname)
    threads: int(workflow.cores)
    shell:
        """
        mkdir -p $(dirname {output.reads})

        cutadapt \
            -g {params.fwd} \
            -a {params.rev} \
            --discard-untrimmed \
            --minimum-length {params.min_len} \
            --json {output.json} \
            -o {output.reads} \
            {input.reads} \
            > {output.log} 2>&1
        """

rule amplicon_qc_trimmed:
    input:
        reads=rules.trim_primers.output.reads
    output:
        qc_metrics="{root}instances/{instance}/{all_samples}/trimmed/qc_metrics.tsv",
        qc_summary="{root}instances/{instance}/{all_samples}/trimmed/qc_summary.tsv"
    params:docker_suffix=lambda wc: f"amplicon_qc_trimmed{wc.all_samples}" if container_runner != "" else "",
            tech=technology if technology == "PacBio" else "illumina", 
    message:
        "Trimmed QC: {wildcards.all_samples}"
    container:
        "{}sif/{}.sif".format(root, toolname)
    threads: int(workflow.cores)
    shell:
        """
            ./fastqc_qc.sh \
                -i {input.reads} \
                -t {technology} \
                -o {params.outdir} \
                -p {threads}
        """
        

rule multiqc:
    input:
        expand(
            "{root}instances/{instance}/{all_samples}/trimmed/cutadapt.json",
            all_samples=all_samples
        )
    output:
        report="{root}instances/{instance}/multiqc/multiqc_report.html"
    params: docker_suffix=lambda wc: f"multiqc{wc.all_samples}" if container_runner != "" else "",
            tech=technology if technology == "PacBio" else "illumina", 
    message:
        "run multiqc sets..."
    container:
        "{}sif/{}.sif".format(root, toolname)
    threads: int(workflow.cores)
    shell:
        """
        mkdir -p {root}instances/{instance}/multiqc

        multiqc \
            {root}instances/{instance} \
            --force \
            -o {root}instances/{instance}/multiqc
        """

        
        

