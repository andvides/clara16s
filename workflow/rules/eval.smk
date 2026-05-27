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


config["instance_name_nosuffix"] = config["instance"].split("___")[0]
try:
    suffix = "___" + config["instance"].split("___")[1]
except IndexError:
    suffix = ""

toolnames = set(
    [
        predictions_path.split("/")[-1].split(".")[0]
        for predictions_path in glob.glob(
            "{}analysis_outputs/*.{}*.predictions.tsv".format(
                root, config["instance_name_nosuffix"]
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
            "{root}evaluations/{toolname}.{instance}-{all_sample}.tab",
            root=root,
            toolname=toolnames,
            instance=instances,
            all_sample=all_samples,
        ),
        expand(
            "{root}evaluations/{toolname}.{instance}.krona.html",
            root=root,
            toolname=toolnames,
            instance=instances,
        ),
 

rule tables:
    input:
        tabs=expand(
            "{root}analysis_outputs/{{toolname}}.{{instance_name_nosuffix}}-{{all_sample}}.predictions.tsv",
            root=root,
        ),
        gnc_ref=expand(
            "{root}repository/{Database_name}/GCN_ref.tsv",
            root=root,
            Database_name=Database_name,
        ),
    output:
        tables_out=expand(
            "{root}evaluations/{{toolname}}.{{instance_name_nosuffix}}{suffix}-{{all_sample}}.tab",
            root=root,
            suffix=suffix,
        ),
    wildcard_constraints:
        all_sample="(?!real_sample_).*",
    message:
        "Tables"
    container:"{}sif/{}.sif".format(root, config["clara16s_master"].split("/")[-1].split(":")[0])
    params: container_suffix=(lambda wildcards: re.sub(r'[\W_]', '',wildcards.toolname.strip()+wildcards.instance_name_nosuffix.strip()+wildcards.all_sample.strip()) if container_runner != "" else "")

    shell:
        """
        {container_runner}{params.container_suffix} {docker_container} krona_convert.py \
        {input.predictions} -o {output.tables_out}
        """


rule krona:
    input:
        tables=expand(
            "{root}evaluations/{{toolname}}.{{instance_name_nosuffix}}-{{all_sample}}.tab",
            root=root,
            toolname=toolnames,
            all_sample=all_samples,
        ),
    output:
        krona_out=expand(
            "{root}evaluations/{toolname}.{instance}.krona.html",
            root=root,
            toolname=toolnames,
            instance=instances,
        ),
      
    message:
        "krona"
    container:"{}sif/{}.sif".format(root, config["clara16s_master"].split("/")[-1].split(":")[0])
    params: container_suffix=(lambda wildcards: re.sub(r'[\W_]', '',wildcards.taxonIDS.strip()+'selecTax') if container_runner != "" else "")
    #threads: 1
    shell:
        """
        {container_runner}{params.container_suffix} {docker_container} ktImportText \
        {input} -o {output}
        """