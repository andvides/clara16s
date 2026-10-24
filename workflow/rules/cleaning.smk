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
dataset=config["Dataset"]
instances = config["instance"]

# SEQUENCING TECHNOLOGY: illumina | nanopore | pacbio (case insensitive)
technology = str(config["technology"]).strip().lower()
if technology not in ("illumina", "nanopore", "pacbio"):
    sys.exit(f"ERROR: technology must be illumina, nanopore or pacbio (got '{config['technology']}')")

# LIBRARY LAYOUT (only used with illumina, optional): auto | paired | interleaved | single
#   paired      -> two files per sample (forward/reverse, e.g. sample_R1.fastq.gz + sample_R2.fastq.gz)
#   interleaved -> one file with forward and reverse reads alternating
#   single      -> single-end reads (no mates); processed like nanopore/pacbio reads
#   auto        -> detected from the file names and read headers (default)
library_layout = str(config.get("library_layout") or "auto").strip().lower()
if library_layout not in reads_layout.LAYOUTS:
    sys.exit(f"ERROR: library_layout must be one of {', '.join(reads_layout.LAYOUTS)} (got '{library_layout}')")

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

# Paired Illumina samples (paired or interleaved files) keep two files (forward/reverse) from the
# merge step to the trimmed reads. Single-end Illumina, Nanopore and PacBio reads stay in one file.
ILLUMINA_INPUT = {}
if technology == "illumina":
    for sample in all_samples:
        try:
            ILLUMINA_INPUT[sample] = reads_layout.resolve_illumina_sample(f"{dataset}/{sample}", library_layout)
        except (OSError, ValueError) as error:
            sys.exit(f"ERROR in sample {sample}: {error}")
        print(f"Sample {sample}: illumina {ILLUMINA_INPUT[sample]['layout']} reads")

    # The output files differ between paired (two files) and single-end (one file) reads, and the
    # rules below are defined once for the whole run, so all the samples must share the structure.
    read_structures = {"single" if reads["layout"] == "single" else "paired" for reads in ILLUMINA_INPUT.values()}
    if len(read_structures) > 1:
        detail = ", ".join(f"{sample}: {ILLUMINA_INPUT[sample]['layout']}" for sample in all_samples)
        sys.exit("ERROR: single-end and paired-end Illumina samples cannot be processed in the same run; "
                 f"run them separately ({detail})")
    PAIRED = read_structures == {"paired"}
else:
    PAIRED = False
    if library_layout != "auto":
        print(f"library_layout '{library_layout}' is only used with illumina reads; ignored for {technology}")

forward_primer = config["forward_primer"] if config["forward_primer"] != '' else 'TTTCTGTTGGTGCTGATATTGC'
reverse_primer = config["reverse_primer"] if config["reverse_primer"] != '' else 'ACTTGCCTGTCGCTCTATCTTC'
min_len = int(config["min_len"]) if config["min_len"] != '' else 100

SAMPLE_DIR = "{root}instances/{instance}/{all_samples}"
if PAIRED:
    MERGED_READS = [SAMPLE_DIR + "/merged/queryReads1.concat.fastq.gz",
                    SAMPLE_DIR + "/merged/queryReads2.concat.fastq.gz"]
    TRIMMED_READS = [SAMPLE_DIR + "/trimmed/queryReads1.trimmed.fastq.gz",
                     SAMPLE_DIR + "/trimmed/queryReads2.trimmed.fastq.gz"]
else:
    MERGED_READS = [SAMPLE_DIR + "/merged/queryReads.concat.fastq.gz"]
    TRIMMED_READS = [SAMPLE_DIR + "/trimmed/queryReads.trimmed.fastq.gz"]


def is_interleaved_sample(wc):
    return PAIRED and ILLUMINA_INPUT[wc.all_samples]["layout"] == "interleaved"


def raw_qc_input(wc):
    # Interleaved files mix both mates, so their raw QC runs on the separated (merged) files
    if is_interleaved_sample(wc):
        return [f"{wc.root}instances/{wc.instance}/{wc.all_samples}/merged/queryReads1.concat.fastq.gz",
                f"{wc.root}instances/{wc.instance}/{wc.all_samples}/merged/queryReads2.concat.fastq.gz"]
    return f"{dataset}/{wc.all_samples}"


def raw_qc_target(wc):
    if is_interleaved_sample(wc):
        return f"{wc.root}instances/{wc.instance}/{wc.all_samples}/merged"
    return f"{dataset}/{wc.all_samples}"


def merge_after_qc(wc):
    # merging waits for the raw QC, except for interleaved samples (their raw QC needs the merged files)
    if is_interleaved_sample(wc):
        return []
    return f"{wc.root}instances/{wc.instance}/{wc.all_samples}/raw/qc_metrics.tsv"


def trimmed_qc_target(wc):
    # paired reads: QC the folder with both mates; otherwise the single trimmed file
    if PAIRED:
        return f"{wc.root}instances/{wc.instance}/{wc.all_samples}/trimmed"
    return f"{wc.root}instances/{wc.instance}/{wc.all_samples}/trimmed/queryReads.trimmed.fastq.gz"


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
        merged_reads=expand(MERGED_READS,
            root=root,
            instance=instances,
            all_samples=all_samples,),
        trimmed_reads=expand(TRIMMED_READS,
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
        reads=raw_qc_input
    output:
        qc_metrics="{root}instances/{instance}/{all_samples}/raw/qc_metrics.tsv",
        qc_summary="{root}instances/{instance}/{all_samples}/raw/qc_summary.tsv"
    message:
        "Raw QC: {wildcards.all_samples}"
    params:container_suffix=lambda wc: f"amplicon_qc_raw{wc.all_samples}" if container_runner != "" else "",
            qc_input=raw_qc_target,
    threads: int(workflow.cores)
    shell:
        """
        {container_runner}{params.container_suffix} -v {dataset}:{dataset} {docker_container} fastqc_qc.sh \
            -i {params.qc_input} \
            -t {technology} \
            -o {root}instances/{wildcards.instance}/{wildcards.all_samples}/raw \
            -p {threads}
        {clean}
        """

if PAIRED:

    # Illumina: one forward and one reverse file per sample.
    # Paired files (several lanes allowed) are concatenated per mate; interleaved files are split.
    rule mergeReads:
        input:
            reads=lambda wc: (ILLUMINA_INPUT[wc.all_samples]["files"] if is_interleaved_sample(wc)
                              else ILLUMINA_INPUT[wc.all_samples]["r1"] + ILLUMINA_INPUT[wc.all_samples]["r2"]),
            qc=merge_after_qc
        output:
            reads1="{root}instances/{instance}/{all_samples}/merged/queryReads1.concat.fastq.gz",
            reads2="{root}instances/{instance}/{all_samples}/merged/queryReads2.concat.fastq.gz"
        message:
            "Merging reads: {wildcards.all_samples}"
        threads: int(workflow.cores)
        run:
            os.makedirs(os.path.dirname(output.reads1), exist_ok=True)
            sample_reads = ILLUMINA_INPUT[wildcards.all_samples]
            if sample_reads["layout"] == "interleaved":
                n_pairs = reads_layout.deinterleave(sample_reads["files"], output.reads1, output.reads2)
                print(f"{wildcards.all_samples}: {n_pairs} read pairs split from interleaved files")
            else:
                reads_layout.concat_fastq(sample_reads["r1"], output.reads1)
                reads_layout.concat_fastq(sample_reads["r2"], output.reads2)

    # Forward primer is expected at the 5' end of the forward reads and the reverse primer
    # (as ordered, 5'->3') at the 5' end of the reverse reads. Pairs are kept only when
    # both reads contain their primer and are at least min_len long.
    rule trim_primers:
        input:
            reads1=rules.mergeReads.output.reads1,
            reads2=rules.mergeReads.output.reads2
        output:
            reads1="{root}instances/{instance}/{all_samples}/trimmed/queryReads1.trimmed.fastq.gz",
            reads2="{root}instances/{instance}/{all_samples}/trimmed/queryReads2.trimmed.fastq.gz",
            json="{root}instances/{instance}/{all_samples}/trimmed/cutadapt.json",
            log="{root}instances/{instance}/{all_samples}/trimmed/cutadapt.log"
        params:container_suffix=lambda wc: f"trim_primers{wc.all_samples}" if container_runner != "" else "",
                fwd=forward_primer,
                rev=reverse_primer,
                min_len=min_len,
        message:
            "Cutadapt (paired): {wildcards.all_samples}"
        threads: int(workflow.cores)
        shell:
            """
            mkdir -p $(dirname {output.reads1})

            {container_runner}{params.container_suffix} {docker_container} cutadapt \
                -g {params.fwd} \
                -G {params.rev} \
                --pair-filter=any \
                --discard-untrimmed \
                --minimum-length {params.min_len} \
                --json {output.json} \
                -o {output.reads1} \
                -p {output.reads2} \
                {input.reads1} {input.reads2} \
                > {output.log} 2>&1
            """

else:

    # Nanopore / PacBio / single-end Illumina: all the files of the sample are concatenated into a single file
    rule mergeReads:
        input:
            reads=lambda wc: f"{dataset}/{wc.all_samples}",
            qc=rules.amplicon_qc_raw.output.qc_metrics
        output:
            reads="{root}instances/{instance}/{all_samples}/merged/queryReads.concat.fastq.gz"
        message:
            "Merging reads: {wildcards.all_samples}"
        threads: int(workflow.cores)
        run:
            os.makedirs(os.path.dirname(output.reads), exist_ok=True)
            reads_layout.concat_fastq(reads_layout.list_fastq(input.reads), output.reads)

    rule trim_primers:
        input:
            reads=rules.mergeReads.output.reads
        output:
            reads="{root}instances/{instance}/{all_samples}/trimmed/queryReads.trimmed.fastq.gz",
            json="{root}instances/{instance}/{all_samples}/trimmed/cutadapt.json",
            log="{root}instances/{instance}/{all_samples}/trimmed/cutadapt.log"
        params:container_suffix=lambda wc: f"trim_primers{wc.all_samples}" if container_runner != "" else "",
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
        reads=TRIMMED_READS
    output:
        qc_metrics="{root}instances/{instance}/{all_samples}/trimmed/qc_metrics.tsv",
        qc_summary="{root}instances/{instance}/{all_samples}/trimmed/qc_summary.tsv"
    params:container_suffix=lambda wc: f"amplicon_qc_trimmed{wc.all_samples}" if container_runner != "" else "",
            qc_input=trimmed_qc_target,
    message:
        "Trimmed QC: {wildcards.all_samples}"
    threads: int(workflow.cores)
    shell:
        """
            {container_runner}{params.container_suffix} {docker_container} fastqc_qc.sh \
                -i {params.qc_input} \
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




