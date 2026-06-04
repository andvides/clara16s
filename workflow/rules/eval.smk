import glob
import shutil

import numpy as np
import pandas as pd
import toolbox
import yaml

root = toolbox.set_root()
tmp = toolbox.set_tmp(root, config)
os.chdir(tmp)

container_runner, docker_container, clean = toolbox.define_docker_command(
    "CLARA16s_master", "clara16s_master",root, tmp, config, int(workflow.cores)
)


instance_name_nosuffix = config["instance"].split("___")[0]
try:
    suffix = "___" + config["instance"].split("___")[1]
except IndexError:
    suffix = ""

toolnames = set(
    [
        predictions_path.split("/")[-1].split(".")[0]
        for predictions_path in glob.glob(
            "{}analysis_outputs/*.{}*.predictions.tsv".format(
                root, instance_name_nosuffix
            )
        )
        
        if predictions_path.split("/")[-1].split(".")[0] not in config["no_evaluation"]
    ]
)


Database_name=config["Database"].split(':')[0].upper()
dataset=config["Dataset"] 
instances = config["instance"]
gcn_normalization = "yes" if config.get("gcn_normalization") == 1 else "no"
#print(gcn_normalization)

all_samples = config["sample"]
print(f"{dataset}+{' '.join(all_samples)}")
for sample in config["sample"]:
    dir_name = f"{root}/instances/{instances}/{sample}"
    os.makedirs(dir_name, exist_ok=True)

localrules:
    target,

rule target:
    input:
        expand(
            "{root}evaluations/{tool}.{instance}-{sample}.tab",
            root=root,
            tool=toolnames,
            instance=instances,
            sample=all_samples,
        ),
        expand(
            "{root}evaluations/{tool}.{instance}.krona.html",
            root=root,
            tool=toolnames,
            instance=instances,
        ),
        expand(
            "{root}evaluations/{tool}.{instance}.otu_table.tsv",
            root=root,
            tool=toolnames,
            instance=instances,
        ),
        expand(
            "{root}evaluations/{tool}.{instance}.taxonomy.tsv",
            root=root,
            tool=toolnames,
            instance=instances,
        ),
        expand(
            "{root}evaluations/{tool}.{instance}.feature_mapping.tsv",
            root=root,
            tool=toolnames,
            instance=instances,
        ),
        expand(
            "{root}evaluations/{tool}.{instance}.report.html",
            root=root,
            tool=toolnames,
            instance=instances,
        ),

rule tables:
    input:
        predictions=(
            "{root}/analysis_outputs/"
            "{tool}.{instance_name}-{sample}.predictions.tsv"
        ),
        #gcn_ref=(
        #    "{root}/repository/{Database_name}/GCN_ref.tsv"
        #),

    output:
        table=(
            "{root}/evaluations/"
            "{tool}.{instance_name}" + suffix + "-{sample}.tab"
        ),

    params:
        container_suffix=lambda wc: (
            re.sub(
                r"[\W_]",
                "",
                f"{wc.tool}{wc.instance_name}{wc.sample}"
            )
            if container_runner != ""
            else ""
        ),

    wildcard_constraints:
        sample="(?!real_sample_).*",

    message:
        "Generating Krona tables",

    container:
        f"{root}/sif/{config['clara16s_master'].split('/')[-1].split(':')[0]}.sif",

    shell:
        """
        {container_runner}{params.container_suffix} \
        {docker_container} \
        krona_convert.py \
        {input.predictions} \
        -o {output.table}
        """


rule krona:
    input:
        lambda wc: expand(
            "{root}evaluations/{tool}.{instance_name}" + suffix + "-{sample}.tab",
            root=root,
            tool=wc.tool,
            instance_name=instance_name_nosuffix,
            sample=all_samples,
        ),

    output:
        "{root}/evaluations/{tool}.{instance}.krona.html",

    params:
        container_suffix=lambda wc: (
            f"krona_{wc.tool}"
            if container_runner != ""
            else ""
        ),

    message:
        "Generating Krona HTML report",

    container:
        f"{root}/sif/{config['clara16s_master'].split('/')[-1].split(':')[0]}.sif",

    shell:
        """
        {container_runner}{params.container_suffix} \
        {docker_container} \
        ktImportText \
        {input} \
        -o {output}
        """

rule otu_table:
    input:
        lambda wc: expand(
            "{root}analysis_outputs/{tool}.{instance_name}" + suffix + "-{sample}.predictions.tsv",
            root=root,
            tool=wc.tool,
            instance_name=instance_name_nosuffix,
            sample=all_samples,
        ),

    output:
        otu_tables="{root}/evaluations/{tool}.{instance}.otu_table.tsv",
        taxonomy="{root}/evaluations/{tool}.{instance}.taxonomy.tsv",
        feature="{root}/evaluations/{tool}.{instance}.feature_mapping.tsv",
    params:
        container_suffix=lambda wc: (
            f"otu_table_{wc.tool}"
            if container_runner != ""
            else ""
        ),

    message:
        "Generating otu_table report",
    container:
        f"{root}/sif/{config['clara16s_master'].split('/')[-1].split(':')[0]}.sif",
    shell:
        """
        {container_runner}{params.container_suffix} \
        {docker_container} \
        build_phyloseq_tables.py  \
        {output.otu_tables} {output.taxonomy} {output.feature}  \
        {input} 
        """

rule report:
    input:
        otu_table="{root}/evaluations/{tool}.{instance}.otu_table.tsv",
        taxonomy="{root}/evaluations/{tool}.{instance}.taxonomy.tsv",
        feature_mapping="{root}/evaluations/{tool}.{instance}.feature_mapping.tsv",
        metadata="{root}/evaluations/metadata.tsv",
        qmd="{root}/evaluations/microbiome_report.qmd",

    output:
        html="{root}/evaluations/{tool}.{instance}.report.html"
    
    params:
        container_suffix=lambda wc: (
            f"html_report_{wc.tool}_{wc.instance}"
            if container_runner != ""
            else ""
        ),
    
    message:
        "Generating html report",
    container:
        f"{root}/sif/{config['clara16s_master'].split('/')[-1].split(':')[0]}.sif",
    threads:
        workflow.cores
    shell:
        """
        {container_runner}{params.container_suffix} \
        {docker_container} report.sh  {input.qmd} \
            {input.otu_table} \
           {input.taxonomy} \
           {input.metadata} \
           {wildcards.tool}.{wildcards.instance}.report.html
        """