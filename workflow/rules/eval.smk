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
print(f"dataset:{dataset}")

all_samples = config["sample"]
print(f"{dataset}+{' '.join(all_samples)}")
for sample in config["sample"]:
    dir_name = f"{root}/instances/{instances}/{sample}"
    os.makedirs(dir_name, exist_ok=True)

shutil.copy(f"{root}/final_results/microbiome_report.qmd", f"{root}/evaluations/microbiome_report.qmd")

# ---------------------------------------------------------------------------
# Web (clara16s_web): HTML and instances.json
# ---------------------------------------------------------------------------
import json
import re as _re
from datetime import date as _date

portal_data_dir = config.get("portal_data_dir", f"{root}final_results")
portal_instance_dir = f"{portal_data_dir}/{instances}"
portal_manifest = f"{portal_data_dir}/instances.json"


def _first_number(path):
    try:
        with open(path) as fh:
            text = fh.read()
    except OSError:
        return None
    m = _re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text)
    return float(m.group()) if m else None


def _sum_numbers(paths):
    vals = [v for v in (_first_number(p) for p in paths) if v is not None]
    return sum(vals) if vals else None


def collect_benchmark(tool):
    runtime_files = [f"{root}analysis_outputs/{tool}.{instances}.runtime_training.txt"]
    memory_files = [f"{root}analysis_outputs/{tool}.{instances}.memory_training.txt"]
    for s in all_samples:
        runtime_files.append(
            f"{root}analysis_outputs/{tool}.{instances}-{s}.runtime_analysis.txt"
        )
        memory_files.append(
            f"{root}analysis_outputs/{tool}.{instances}-{s}.memory_analysis.txt"
        )

    runtime = _sum_numbers(runtime_files)
    memory = _sum_numbers(memory_files)
    bench = {}
    if runtime is not None:
        bench["runtime_seconds"] = round(runtime, 2)
    if memory is not None:
        # los archivos de memoria se asumen en MB; si vienen en KB, ajustar aqui
        bench["peak_memory_mb"] = round(memory, 2)
    return bench or None



localrules:
    target,
    portal,


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
        f"{portal_instance_dir}/.portal_done",


rule portal:
    """Copia Krona/Quarto/MultiQC al portal y actualiza instances.json."""
    input:
        krona=expand(
            "{root}evaluations/{tool}.{instance}.krona.html",
            root=root,
            tool=toolnames,
            instance=instances,
        ),
        krona_combined=expand(
            "{root}evaluations/{tool}.{instance}.krona_combined.html",
            root=root,
            tool=toolnames,
            instance=instances,
        ),
        report=expand(
            "{root}evaluations/{tool}.{instance}.report.html",
            root=root,
            tool=toolnames,
            instance=instances,
        ),
    output:
        marker=f"{portal_instance_dir}/.portal_done",
        manifest=portal_manifest,
    message:
        "Publicando resultados en el portal web",
    run:
        os.makedirs(portal_instance_dir, exist_ok=True)

        tools_entry = []
        for tool in sorted(toolnames):
            krona_src = f"{root}/evaluations/{tool}.{instances}.krona.html"
            report_src = f"{root}/evaluations/{tool}.{instances}.report.html"
            krona_dst = f"{portal_instance_dir}/{tool}_krona.html"
            report_dst = f"{portal_instance_dir}/{tool}_report.html"
            shutil.copyfile(krona_src, krona_dst)
            shutil.copyfile(
                f"{root}/evaluations/{tool}.{instances}.krona_combined.html",
                f"{portal_instance_dir}/{tool}_krona_combined.html",
            )
            shutil.copyfile(report_src, report_dst)

            entry = {
                "id": tool,
                "name": tool,
                "krona": f"data/{instances}/{tool}_krona.html",
                "krona_combined": f"data/{instances}/{tool}_krona_combined.html",
                "quarto": f"data/{instances}/{tool}_report.html",
            }
            bench = collect_benchmark(tool)
            if bench:
                entry["benchmark"] = bench
            tools_entry.append(entry)

        instance_entry = {
            "id": instances,
            "name": config.get("instance_name", instances),
            "description": config.get(
                "instance_description",
                f"Amplicon analysis of dataset {os.path.basename(str(dataset).rstrip('/'))}.",
            ),
            "date": config.get("date", _date.today().isoformat()),
            "tags": config.get("tags", [Database_name, os.path.basename(str(dataset).rstrip("/"))]),
            "tools": tools_entry,
            "metadata": {
                "Dataset": str(dataset),
                "Database": config["Database"],
                "Samples": len(all_samples),
                "Tools": ", ".join(sorted(toolnames)),
                "GCN normalization": gcn_normalization,
            },
        }

        multiqc_src = f"{root}/instances/{instances}/multiqc/multiqc_report.html"
        if os.path.exists(multiqc_src):
            shutil.copyfile(multiqc_src, f"{portal_instance_dir}/multiqc.html")
            instance_entry["multiqc"] = f"data/{instances}/multiqc.html"
        else:
            print(f"Advertencia: El archivo no existe en la ruta {multiqc_src}")

        # --- merge acumulativo del manifiesto ---
        manifest = {
            "title": config.get("portal_title", "Amplicon Sequencing Analysis"),
            "subtitle": config.get(
                "portal_subtitle",
                "Interactive results from amplicon sequencing analysis experiments",
            ),
            "instances": [],
        }
        if os.path.exists(portal_manifest):
            try:
                with open(portal_manifest) as fh:
                    previous = json.load(fh)
                manifest["title"] = previous.get("title", manifest["title"])
                manifest["subtitle"] = previous.get("subtitle", manifest["subtitle"])
                manifest["instances"] = [
                    i for i in previous.get("instances", []) if i.get("id") != instances
                ]
            except (OSError, ValueError):
                print("instances.json ilegible: se reconstruye desde cero")

        manifest["instances"].append(instance_entry)
        manifest["instances"].sort(key=lambda i: (i.get("date", ""), i.get("id", "")))

        with open(portal_manifest, "w") as fh:
            json.dump(manifest, fh, indent=2, ensure_ascii=False)

        with open(output.marker, "w") as fh:
            fh.write(instances + "\n")

        print(f"Portal actualizado: {portal_manifest}")


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


rule krona_combined:
    input:
        lambda wc: expand(
            "{root}evaluations/{tool}.{instance_name}" + suffix + "-{sample}.tab",
            root=root,
            tool=wc.tool,
            instance_name=instance_name_nosuffix,
            sample=all_samples,
        ),

    output:
        "{root}/evaluations/{tool}.{instance}.krona_combined.html",

    params:
        container_suffix=lambda wc: (
            f"krona_combined_{wc.tool}"
            if container_runner != ""
            else ""
        ),

    message:
        "Generating combined Krona HTML report",

    container:
        f"{root}/sif/{config['clara16s_master'].split('/')[-1].split(':')[0]}.sif",

    shell:
        """
        {container_runner}{params.container_suffix} \
        {docker_container} \
        ktImportText -c \
        {input:q} \
        -o {output:q}
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
        metadata=expand("{dataset}/metadata.tsv",dataset=dataset),
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
        {container_runner}{params.container_suffix} -v {dataset}:{dataset} \
        {docker_container} report.sh  {input.qmd} \
            {input.otu_table} \
            {input.taxonomy} \
            {input.metadata} \
            {wildcards.tool}.{wildcards.instance}.report.html
        """
