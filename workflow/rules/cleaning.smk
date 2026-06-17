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
instances = config["instance"]
technology=config["technology"]

#CLARA16s_master container used for basic process
container_runner, docker_container, clean = toolbox.define_docker_command(
    "clara16s_master", "clara16s_master",root, tmp, config, int(workflow.cores)
)
toolbox.validate_config(config)

#READS CONFIG
all_samples = config["sample"]
print(f"{dataset}{' '.join(all_samples)}")
for sample in config["sample"]:
    dir_name = f"{root}/instances/{instances}/{sample}"
    os.makedirs(dir_name, exist_ok=True)

forward_primer = config["forward_primer"] if config["forward_primer"] != '' else 'TTTCTGTTGGTGCTGATATTGC'
reverse_primer = config["reverse_primer"] if config["reverse_primer"] != '' else 'ACTTGCCTGTCGCTCTATCTTC'
min_len = int(config["min_len"]) if config["min_len"] != '' else 100
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
    params:container_suffix=lambda wc: f"amplicon_qc_raw{wc.all_samples}" if container_runner != "" else "",
    threads: int(workflow.cores)
    shell:
        """
        {container_runner}{params.container_suffix} -v {dataset}:{dataset} {docker_container} fastqc_qc.sh \
            -i {input.reads} \
            -t {technology} \
            -o {root}instances/{wildcards.instance}/{wildcards.all_samples}/raw \
            -p {threads}
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
    params: container_suffix=lambda wc: f"mergeReads{wc.all_samples}" if container_runner != "" else "",
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
    params:container_suffix=lambda wc: f"trim_primers{wc.all_samples}" if container_runner != "" else "",
            tech=technology if technology == "PacBio" else "illumina", 
            fwd=forward_primer,
            rev=reverse_primer,
            min_len=min_len,
    message:
        "Cutadapt: {wildcards.all_samples}"
    threads: int(workflow.cores)
    shell:
        """
        mkdir -p $(dirname {output.reads})

        {container_runner}{params.container_suffix} {docker_container} cutadapt \
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
    params:container_suffix=lambda wc: f"amplicon_qc_trimmed{wc.all_samples}" if container_runner != "" else "",
            tech=technology if technology == "PacBio" else "illumina", 
    message:
        "Trimmed QC: {wildcards.all_samples}"
    threads: int(workflow.cores)
    shell:
        """
            {container_runner}{params.container_suffix} {docker_container} fastqc_qc.sh \
                -i {input.reads} \
                -t {technology} \
                -o {root}instances/{wildcards.instance}/{wildcards.all_samples}/trimmed \
                -p {threads}
        """
        

rule multiqc:
    input:
        expand(
            "{root}instances/{instance}/{all_samples}/trimmed/cutadapt.json",
            root=root,
            instance=instances,
            all_samples=all_samples
        )
    output:
        report="{root}instances/{instance}/multiqc/multiqc_report.html"
    params: container_suffix=lambda wc: f"multiqc{wc.instance}" if container_runner != "" else "",
            tech=technology if technology == "PacBio" else "illumina", 
    message:
        "run multiqc sets..."
    threads: int(workflow.cores)
    shell:
        """
        mkdir -p {root}instances/{wildcards.instance}/multiqc

        {container_runner}{params.container_suffix} {docker_container} multiqc \
            {root}instances/{wildcards.instance} \
            --force \
            -o {root}instances/{wildcards.instance}/multiqc
        """

        
        

