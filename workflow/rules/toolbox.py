import glob
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd
import yaml
from Bio import Entrez, SeqIO
from ete3 import NCBITaxa
pd.options.mode.chained_assignment = None  # default='warn'

pd.set_option("display.max_rows", 20)
pd.set_option("display.max_columns", 999)
pd.set_option("display.max_colwidth", 9999)

Entrez.email = "A.N.Other@example.com"

GTDB_LINEAGE_MAPPING = {'d__': 'domain', 'p__': 'phylum', 'c__': 'class',
                        'o__': 'order', 'f__': 'family', 'g__': 'genus', 's__': 'species'}

GTDB_LINEAGE_MAPPING_REV = {GTDB_LINEAGE_MAPPING[v]: v for v in GTDB_LINEAGE_MAPPING}

ORDERED_LINEAGE = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'domain']


# function for all smk

def validate_config(config):
    error_msg = ''
    try:
        if '-' in config['instance_name']:
            error_msg = 'The character - is not allowed in instance names'
        if '-' in config['container']:
            error_msg = 'The character - is not allowed in tool names'

    except KeyError:
        pass
    if error_msg != '':
        exit(error_msg)


def set_root():
    """
    The root is where the data content is produced
    :return: the root path
    """
    root = os.environ["CLARA16s_ROOT"]
    if not root.endswith("/"):
        root += "/"
    root += "data/"
    return root


def set_tmp(root, config):
    """
    The tmp should be set as the current working directory in the
    snakemake file. If a value is mentioned manually in the config
    with the tmp parameter, use it (useful in analysis to recover previous content)
    else use the datetime (timestamp) to provide a clean folder.

    :param root: the path to the root folder
    :param config: the config dictionnary
    :return: the path to the tmp folder
    """
    if "tmp" in config:
        tmp = "{}tmp/{}/".format(root, config["tmp"])
    else:
        tmp = "{}tmp/{}/".format(root, config["datetime"])
    try:
        os.mkdir(tmp)
    except FileExistsError:
        pass
    return tmp


def get_genome_size(g):
    """
    TODO
    :param g:
    :return:
    """
    with gzip.open(g, "rt") as handle:
        return sum(
            [len(record.seq) for record in SeqIO.parse(handle, "fasta")]
        )


def define_docker_command(name, config_key, root, tmp, config, cores):
    """
    TODO
    :param name:
    :param config_key:
    :param root:
    :param tmp:
    :param config:
    :param cores:
    :return:
    """

    if "--use-singularity" not in sys.argv:
        container_runner = "docker run --memory='{}' --memory-swap='{}' --cpus {} -e cpus={} " \
                           "-v {}../:{}../ -w {} --name {}_{}_".format(
            config['docker_memory'], config['docker_memory'], cores, cores, root, root, tmp, name, config["datetime"]
        )
        docker_container = config[config_key]
        clean = 'docker rm $(docker ps -a | grep -v "Up " | grep "{}_" | cut -f1 -d" ") || true'.format(
            name
        )
    else:
        container_runner = ""
        clean = ""
        docker_container = ""

    return container_runner, docker_container, clean


# function for genomes.download


def define_accessions(root):
    """
    TODO
    :param root:
    :return:
    """
    try:
        v_acc = [
            a.split("\t")[0]
            for a in open("{}repository/vir_list.tsv".format(root))
            if ("ncbi_genbank_assembly_accession" not in a and not a.startswith('#'))
        ]
        e_acc = [
            a.split("\t")[0]
            for a in open("{}repository/euk_list.tsv".format(root))
            if ("ncbi_genbank_assembly_accession" not in a and not a.startswith('#'))
        ]
        p_acc = [
            a.split("\t")[3]
            for a in open("{}repository/prok_list.tsv".format(root))
            if ("ncbi_genbank_assembly_accession" not in a and not a.startswith('#'))
        ]
        pn_acc = [
            a.split("\t")[0]
            for a in open("{}repository/prokncbi_list.tsv".format(root))
            if ("ncbi_genbank_assembly_accession" not in a and not a.startswith('#'))
        ]
        all_acc = v_acc + e_acc + p_acc + pn_acc
        for element in v_acc + e_acc + p_acc + pn_acc:
            if os.path.exists("{}repository/genomes/{}.fna.gz".format(root, element)):
                all_acc.remove(element)
        return all_acc
    except FileNotFoundError:
        exit("missing genome list file")


# function for repository.setup

def _extract_size(xml):
    """
    TODO
    :param xml:
    :return:
    """
    try:
        return int(
            re.findall("[0-9]+", re.findall("total_length.*?</Stat>", xml["Meta"])[0])[
                0
            ]
        )
    except IndexError:
        return 0


def _get_primary_id(acc):
    """
    TODO
    :param acc:
    :return:
    """
    log=open('primary_id', 'a')
    log.write('{}\n'.format(acc))
    log.close()
    try:
        result = Entrez.read(Entrez.esearch(db="assembly", term=acc))["IdList"]
        log=open('primary_id', 'a')
        log.write('{}\n'.format(result))
        log.close()
    except RuntimeError:
        return 0
    try:
        return result[0]
    except IndexError:
        return 0


def create_list(clade, input, output, config):
    """
    TODO
    :param clade:
    :param input:
    :param output:
    :param config:
    :return:
    """
    data = pd.read_csv(
        input.genbank_assembly[0],
        header=0,
        sep="\t",
    )
    # filter
    data = data[
        data.apply(
            lambda row: len(
                set(row["__lineage"].split(";"))
                & set([str(v) for v in config["{}_clade".format(clade)]])
            )
                        > 0,
            axis=1,
        )
    ]
    # keep those with at least two genomes for one taxid
    data = data[data.taxid.groupby(data.taxid).transform(func=len) > 1]
    # reset column index
    data.reset_index(inplace=True, drop=True)
    # rename assembly accession
    data.rename(
        columns={"# assembly_accession": "ncbi_genbank_assembly_accession"},
        inplace=True,
    )
    data = data.groupby('taxid').tail(int(config['max_genomes_per_organism']))
    # fetch the size
    # first get all primary ids in separate requests. biopython should allow 3 requests/sec max.
    # to avoid more (will likely block the requests), But only one rule should call biopython at any moment
    # I prevent create_euk_list and create_vir_list and create_prokncbi_list to run together
    # by creating an artificial dependence of one to another in the rules
    data["__primary_id"] = data.apply(
        lambda row: _get_primary_id(row["ncbi_genbank_assembly_accession"]), axis=1
    )
    data.drop(data[data["__primary_id"] == 0].index, axis=0, inplace=True)
    # get all info in one requests... but splitting by 10000 as crashes seem to happen otherwise
    all_summaries = []
    for i in range(0, len(list(data["__primary_id"])), 10000):
        a = i
        if i + 10000 < len(list(data["__primary_id"])):
            z = i + 10000
        else:
            z = len(list(data["__primary_id"]))
        all_summaries += Entrez.read(
            Entrez.esummary(
                db="assembly", id=",".join(list(data["__primary_id"])[a:z])
            )
        )["DocumentSummarySet"]["DocumentSummary"]
    data.reset_index(drop=True, inplace=True)
    data["__size"] = data.apply(
        lambda row: _extract_size(all_summaries[row.name]), axis=1
    )
    data.drop(["__primary_id"], inplace=True, axis=1)
    data.drop(data[data["__size"] == 0].index, axis=0, inplace=True)
    # recheck to keep only those with two genome per taxid still in
    data = data[data.taxid.groupby(data.taxid).transform(func=len) > 1]
    # keep a given nb of genome per organism according to the defined limit
    if clade == 'euk':
        outp = open(output.euk_list[0], 'w')
        outp.write('# CLARA16s content selected from NCBI data\n')
        outp.close()
        data.to_csv(
            output.euk_list[0], sep="\t", header=True, mode='a', index=False,
        )
    elif clade == 'vir':
        outp = open(output.vir_list[0], 'w')
        outp.write('# CLARA16s content selected from NCBI data\n')
        outp.close()
        data.to_csv(
            output.vir_list[0], sep="\t", header=True, mode='a', index=False,
        )
    else:
        outp = open(output.prokncbi_list[0], 'w')
        outp.write('# CLARA16s content selected from NCBI data\n')
        outp.close()
        data.to_csv(
            output.prokncbi_list[0], sep="\t", header=True, mode='a', index=False,
        )



# function for instance.definition

def _update_selected_ids(selection, taxon_field, selected_ids):
    """
    TODO
    :param selection:
    :param taxon_field:
    :param selected_ids:
    :return:
    """
    if len(selection) > 0:
        for acc in selection["ncbi_genbank_assembly_accession"]:
            taxon = selection[selection["ncbi_genbank_assembly_accession"] == acc][
                taxon_field
            ].values
            if taxon[0] not in selected_ids:
                selected_ids[taxon[0]] = []
            selection_acc = list(
                selection[selection[taxon_field] == taxon[0]][
                    "ncbi_genbank_assembly_accession"
                ].values
            )
            selected_ids[taxon[0]] += selection_acc


def do_selection(output_file, prok_list, euk_list, vir_list, selected_ids, host_size, config, negative=False):
    """
    TODO
    :param output_file:
    :param prok_list:
    :param euk_list:
    :param vir_list:
    :param selected_ids:
    :param host_size:
    :param config:
    :param negative:
    :return:
    """

    # percent of reads * reads gives your reads
    # use a distribution to decide abundance of each of the x
    # normalize this to reads
    # figure out how value relate to genome size, and voila

    # get an empty version of each list as default
    selection_prok = prok_list[
        prok_list["ncbi_genbank_assembly_accession"] == "nothing"
        ]
    selection_euk = euk_list[
        euk_list["ncbi_genbank_assembly_accession"] == "nothing"
        ]
    selection_vir = vir_list[
        vir_list["ncbi_genbank_assembly_accession"] == "nothing"
        ]

    if not negative or 'prok' not in config['targets']:
        if config['prok_taxa'] and len(config['prok_taxa']) > 0:
            prok_list_custom = None
            for entry in config['prok_taxa']:
                prok_list_tmp = prok_list[prok_list['__lineage'].str.contains(entry)]
                if type(prok_list_custom) == pd.core.frame.DataFrame:
                    prok_list_custom = prok_list_custom.append(prok_list_tmp)
                else:
                    prok_list_custom = prok_list_tmp
            prok_list = prok_list_custom
        if config['use_prokncbi'] == 1:
            if len(prok_list['__lineage']) < config["prok_nb"]:
                n_prok = len(prok_list['__lineage'])
            else:
                n_prok = config["prok_nb"]

            selection_prok = prok_list.drop_duplicates().sample(n=n_prok)
            _update_selected_ids(selection_prok, 'taxid', selected_ids)

        else:
            if len(prok_list['gtdb_taxonomy']) < config["prok_nb"]:
                n_prok = len(prok_list['gtdb_taxonomy'])
            else:
                n_prok = config["prok_nb"]
            selection_prok = prok_list.drop_duplicates().sample(n=n_prok)
            _update_selected_ids(selection_prok, 'gtdb_taxonomy', selected_ids)

    if not negative or 'euk' not in config['targets']:
        if config['euk_taxa'] and len(config['euk_taxa']) > 0:
            euk_list_custom = None
            for entry in config['euk_taxa']:
                euk_list_tmp = euk_list[euk_list['__lineage'].str.contains(entry)]
                if type(euk_list_custom) == pd.core.frame.DataFrame:
                    euk_list_custom = euk_list_custom.append(euk_list_tmp)
                else:
                    euk_list_custom = euk_list_tmp
            euk_list = euk_list_custom
        if len(euk_list['__lineage']) < config["euk_nb"]:
            n_euk = len(euk_list['__lineage'])
        else:
            n_euk = config["euk_nb"]
        selection_euk = euk_list.drop_duplicates().sample(n=n_euk)
        _update_selected_ids(selection_euk, 'taxid', selected_ids)

    if not negative or 'vir' not in config['targets']:
        if config['vir_taxa'] and len(config['vir_taxa']) > 0:
            vir_list_custom = None
            for entry in config['vir_taxa']:
                vir_list_tmp = vir_list[vir_list['__lineage'].str.contains(entry)]
                if type(vir_list_custom) == pd.core.frame.DataFrame:
                    vir_list_custom = vir_list_custom.append(vir_list_tmp)
                else:
                    vir_list_custom = vir_list_tmp
            vir_list = vir_list_custom
        if len(vir_list['__lineage']) < config["vir_nb"]:
            n_vir = len(vir_list['__lineage'])
        else:
            n_vir = config["vir_nb"]
        selection_vir = vir_list.drop_duplicates().sample(n=n_vir)
        _update_selected_ids(selection_vir, 'taxid', selected_ids)

    # define the relative abundance from a lognormal distrib for each group
    # apply this to the proportion of reads dedicated to this group
    # remove those with 0 reads
    # for each genome, produce the coverage that matches this nb of reads:w
    if not negative or 'vir' not in config['targets']:
        v_abu = np.random.lognormal(1, config["sd_lognormal_vir"], n_vir)
        v_reads = [
            int(v)
            for v in v_abu
                     * config["total_reads_nb"]
                     * config["vir_ratio"]
                     / sum(v_abu)
        ]
    else:
        v_reads = []

    selection_vir["reads"] = v_reads
    # drop from selected_ids and selection_vir if reads = 0
    selection_vir[selection_vir["reads"] == 0].apply(
        lambda row: selected_ids[row["taxid"]].remove(
            row["ncbi_genbank_assembly_accession"]
        ),
        axis=1,
    )
    selection_vir.drop(
        selection_vir[selection_vir["reads"] == 0].index, axis=0, inplace=True
    )
    selection_vir["coverage"] = (
            config["read_length"] * selection_vir["reads"] / selection_vir["__size"]
    )

    if not negative or 'prok' not in config['targets']:
        p_abu = np.random.lognormal(1, config["sd_lognormal_prok"], n_prok)
        p_reads = [
            int(p)
            for p in p_abu
                     * config["total_reads_nb"]
                     * config["prok_ratio"]
                     / sum(p_abu)
        ]
    else:
        p_reads = []

    selection_prok["reads"] = p_reads
    if config['use_prokncbi'] == 1:
        selection_prok[selection_prok["reads"] == 0].apply(
            lambda row: selected_ids[row['taxid']].remove(
                row["ncbi_genbank_assembly_accession"]
            ),
            axis=1,
        )
    else:
        selection_prok[selection_prok["reads"] == 0].apply(
            lambda row: selected_ids[row["gtdb_taxonomy"]].remove(
                row["ncbi_genbank_assembly_accession"]
            ),
            axis=1,
        )
    selection_prok.drop(
        selection_prok[selection_prok["reads"] == 0].index, axis=0, inplace=True
    )
    if config['use_prokncbi']:
        selection_prok["coverage"] = (
                config["read_length"]
                * selection_prok["reads"]
                / selection_prok["__size"]
        )
    else:
        selection_prok["coverage"] = (
                config["read_length"]
                * selection_prok["reads"]
                / selection_prok["genome_size"]
        )
    if not negative or 'euk' not in config['targets']:
        e_abu = np.random.lognormal(1, config["sd_lognormal_euk"], n_euk)
        e_reads = [
            int(e)
            for e in e_abu
                     * config["total_reads_nb"]
                     * config["euk_ratio"]
                     / sum(e_abu)
        ]
    else:
        e_reads = []
    selection_euk["reads"] = e_reads
    selection_euk[selection_euk["reads"] == 0].apply(
        lambda row: selected_ids[row["taxid"]].remove(
            row["ncbi_genbank_assembly_accession"]
        ),
        axis=1,
    )
    selection_euk.drop(
        selection_euk[selection_euk["reads"] == 0].index, axis=0, inplace=True
    )
    selection_euk["coverage"] = (
            config["read_length"] * selection_euk["reads"] / selection_euk["__size"]
    )
    value_to_use="coverage"
    host_value = config["read_length"]* config["total_reads_nb"]* config["host_ratio"] / host_size
    if "long_reads" in config and config["long_reads"] == 1:
        value_to_use="reads"
        host_value = config["total_reads_nb"] * config["host_ratio"]

    yaml_content = {
        "instance_name": output_file.split("/")[-1].split(".")[0].split('-')[0],
        "sample_name": "-".join(output_file.split("/")[-1].split(".")[0:-1]),
        "euk": [
            {
                k: float(
                    selection_euk[
                        selection_euk["ncbi_genbank_assembly_accession"] == k
                        ][value_to_use].values[0]
                )
            }
            for k in selection_euk["ncbi_genbank_assembly_accession"].values
        ],
        "prok": [
            {
                k: float(
                    selection_prok[
                        selection_prok["ncbi_genbank_assembly_accession"] == k
                        ][value_to_use].values[0]
                )
            }
            for k in selection_prok["ncbi_genbank_assembly_accession"].values
        ],
        "vir": [
            {
                k: float(
                    selection_vir[
                        selection_vir["ncbi_genbank_assembly_accession"] == k
                        ][value_to_use].values[0]
                )
            }
            for k in selection_vir["ncbi_genbank_assembly_accession"].values
        ],
        "host": host_value,
        "tech": config['tech'],
        "read_length": config['read_length'],
        "total_reads_nb": config['total_reads_nb'],
        "fragment_mean_size": config['fragment_mean_size'],
        "fragment_sd": config['fragment_sd'],
        "max_base_quality": config['max_base_quality'],
        "core_simul": config['core_simul'],
        "long_reads": config['long_reads'] if 'long_reads' in config else 0
    }
    if 'long_reads' in config and config['long_reads'] == 1:
        for e in yaml_content["prok"]:
            for k in e:
                e[k] = int(e[k])
        for e in yaml_content["euk"]:
            for k in e:
                e[k] = int(e[k])
        for e in yaml_content["vir"]:
            for k in e:
                e[k] = int(e[k])
        yaml_content["host"] = int(yaml_content["host"])
    with open(output_file, "w") as yaml_f:
        yaml.dump(yaml_content, yaml_f)

    return selected_ids


def create_full_ref(full_list, selected_ids, output_file, columns, unknown_lineages):
    """
    TODO
    :param full_list:
    :param selected_ids:
    :param output_file:
    :param columns:
    :return:
    """
    list_ref = full_list[
        ~full_list["ncbi_genbank_assembly_accession"].isin(
            [
                item
                for sublist in list(selected_ids.values())
                for item in sublist
            ]
        )
    ]
    outp = open(output_file, 'w')
    outp.write('# CLARA16s content selected from NCBI+GTDB data\n')
    outp.close()
    list_ref.to_csv(
        output_file,
        sep="\t",
        header=True,
        index=False,
        columns=columns,
        mode='a'
    )

    return list_ref


def create_best_ref(full_list, output_file, sort_columns, write_columns):
    """
    TODO
    :param full_list:
    :param output_file:
    :param sort_columns:
    :param write_columns:
    :return:
    """
    # the best_ref, one for each taxid or gtdb entry,
    # ref genome then complete genome then chromosome then contig then
    # rename assembly level so alphabetic order will do. 1=Complete Genome, 2=Chromosome, 4=Contig, 3=Scaffold

    full_list["assembly_level"][
        full_list["assembly_level"] == "Complete Genome"
        ] = "1"
    full_list["assembly_level"][full_list["assembly_level"] == "Contig"] = "4"
    full_list["assembly_level"][
        full_list["assembly_level"] == "Scaffold"
        ] = "3"
    full_list["assembly_level"][
        full_list["assembly_level"] == "Chromosome"
        ] = "2"
    full_list.sort_values(
        by=sort_columns,
        ascending=[True, False, True, False],
        inplace=True,
    )
    full_list.drop_duplicates(inplace=True, subset=sort_columns[0], keep="first")
    outp = open(output_file, 'w')
    outp.write('# CLARA16s content selected from NCBI+GTDB data\n')
    outp.close()
    full_list.to_csv(
        output_file,
        sep="\t",
        header=True,
        index=False,
        columns=write_columns,
        mode='a'
    )


# function for analysis


def get_samples(root, config):
    """
    TODO
    :param root:
    :param config:
    :return:
    """
    samples = []
    for f in glob.glob(
            "{}instances/{}/{}-*".format(root, config["instance_name"], config["instance_name"])
    ):
        if os.path.isdir(f) and "real_sample" not in f:
            samples.append(f.split("/")[-1])
    return samples


def define_tmp_folder_rules_inputs(actual_file, mock_file, real_target_files):
    """
    The goal here is to allow the user to delete the tmp folders without having snakemake rerunning
    even though tmp are inputs of rules. So given the target files that have to be really absent for the pipeline to rerun
    test if one is missing. if yes, return the real input file for the rule (will trigger it)
    else return the mock file, which is supposed to be present and untouched, it will never trigger the rule no matter
    whether the tmp files are there or not
    Note: this goes against the modification of how to rerun introduced in smk 7.8.0
    See
    https://github.com/snakemake/snakemake/issues/1694
    fortunately, the option --rerun-triggers mtime restore the previous behavior so we can keep working like this
    :param actual_file:
    :param mock_file:
    :param real_target_files:
    :return:
    """
    for f_list in real_target_files:
        for f in f_list:
            if not os.path.exists(f): # test whether the important targets are there
                return actual_file # if one is missing, return the real input of the rule
    return mock_file # if none of the real outputs of the pipeline are missing, return this mock untouched file


# function for evaluation

def _is_keep_gtdb(row, keep_idx):
    """
    TODO
    :param row:
    :return:
    """
    try:
        return ORDERED_LINEAGE.index(GTDB_LINEAGE_MAPPING[row['organism'][0:3]]) <= keep_idx
    except KeyError:
        # if this is the case, it is because it is the lowest level (accession) so we keep it
        # as it is below the targeted level
        return True


def _is_keep_ncbi(row, keep_idx):
    """
    TODO
    :param row:
    :return:
    """
    try:
        ncbi = NCBITaxa()
        taxid = ncbi.get_name_translator([row['organism']])
        rank = ncbi.get_rank(taxid[row['organism']])
        if rank[taxid[row['organism']][0]] not in ORDERED_LINEAGE:
            return True
        return ORDERED_LINEAGE.index(rank[taxid[row['organism']][0]]) <= keep_idx
    except KeyError:
        return True
    except IndexError:
        return True


def _define_evaluation_taxa_gtdb(row, ref_genomes, rank, dontsum):
    """
    TODO
    :param row:
    :return:
    """
    if row['organism'].startswith(GTDB_LINEAGE_MAPPING_REV[rank]):
        return row['organism']
    elif dontsum:
        return np.nan
    else:

        try:
            return re.findall('{}[^;]*'.format(GTDB_LINEAGE_MAPPING_REV[rank]),
                              ref_genomes[ref_genomes['__gtdb_accession'] == row['organism']].gtdb_taxonomy.values[0])[
                0]
        except IndexError:
            try:
                return re.findall('{}[^;]*'.format(GTDB_LINEAGE_MAPPING_REV[rank]), ref_genomes[
                    ref_genomes['gtdb_taxonomy'].str.contains(row['organism'], na=False)].gtdb_taxonomy.values[0])[0]
            except IndexError:
                # this is the case a ncbi accession was passed but gtdb is used to assess.
                # The function load_truth below uses that
                return re.findall('{}[^;]*'.format(GTDB_LINEAGE_MAPPING_REV[rank]), ref_genomes[
                    ref_genomes['ncbi_genbank_assembly_accession'].str.contains(
                        row['organism'])].gtdb_taxonomy.values[
                    0])[0]


def _define_evaluation_taxa_ncbi(row, ref_genomes, rank, dontsum):
    """
    TODO
    :param row:
    :return:
    """
    ncbi = NCBITaxa()
    if dontsum:
        try:
            taxid = ncbi.get_name_translator([row['organism']])[row['organism']][0]
            if ncbi.get_rank([taxid])[taxid] != rank:
                return np.nan
        except KeyError:
            return np.nan
    try:
        taxid = ncbi.get_name_translator([row['organism']])[row['organism']][0]
        lineage = ncbi.get_lineage(taxid)
        for taxid in lineage:
            if ncbi.get_rank([taxid])[taxid] == rank:
                return ncbi.get_taxid_translator([taxid])[taxid]
    except KeyError:
        try:
            lineage = ref_genomes[ref_genomes['ncbi_genbank_assembly_accession'] == row['organism']].__lineage.values[0]
        except IndexError:
            return ''
        for name in lineage.split(';'):
            taxid = ncbi.get_name_translator([name])[name][0]
            if ncbi.get_rank([taxid])[taxid] == rank:
                return name


def sum_to_taxlevel(predictions_input, ref_genomes_files, taxonomy, rank, config, outp_file=None):
    """
    TODO note explain that the tool should not report a read twice
    :param predictions_input:
    :param ref_genomes_files:
    :param taxonomy:
    :param rank:
    :param outp_file:
    :return:
    """
    try:
        tool=predictions_input.split('/')[-1].split('.')[0]
    except AttributeError:
        tool=np.nan

    dontsum=False
    if tool in config['no_sum_tools']:
        dontsum=True

    if type(predictions_input) == pd.core.frame.DataFrame:
        predictions = predictions_input
    else:
        predictions = pd.read_csv(predictions_input, header=1, sep="\t")
    predictions['evaluation_taxa'] = np.nan
    if outp_file:
        outp = open(outp_file, 'w')
        outp.write('#CLARA16s_{}\n'.format(taxonomy))
        outp.close()
        if len(predictions) == 0:
            predictions.to_csv(outp_file, index=False, sep="\t",
                               columns=['evaluation_taxa',
                                        'reads',
                                        'abundance'], mode='a')
            return

    ref_genomes_df = []

    for filename in ref_genomes_files:
        df = pd.read_csv(filename, header=1, sep="\t")
        ref_genomes_df.append(df)

    ref_genomes = pd.concat(ref_genomes_df, axis=0, ignore_index=True)
    if taxonomy == 'gtdb':
        keep_idx = ORDERED_LINEAGE.index(rank)
        if len(predictions) > 0:
            predictions['keep'] = \
                predictions.apply(lambda row: _is_keep_gtdb(row, keep_idx), axis=1)
        else:
            predictions['keep'] = None
    else:
        if len(predictions) > 0:
            keep_idx = ORDERED_LINEAGE.index(rank)
            predictions['keep'] = \
                predictions.apply(lambda row: _is_keep_ncbi(row, keep_idx), axis=1)
        else:
            predictions['keep'] = None

    predictions_filtered = predictions[predictions['keep']]
    if len(predictions_filtered) == 0:
        if outp_file:
            predictions_filtered.to_csv(outp_file, index=False,sep="\t",columns=['evaluation_taxa','reads','abundance'],mode='a')

        return predictions_filtered

    if taxonomy == 'gtdb':
        predictions_filtered.loc[:, 'evaluation_taxa'] = predictions_filtered.apply(
            lambda row: _define_evaluation_taxa_gtdb(row, ref_genomes, rank, dontsum),
            axis=1)
    else:
        predictions_filtered.loc[:, 'evaluation_taxa'] = predictions_filtered.apply(
            lambda row: _define_evaluation_taxa_ncbi(row, ref_genomes, rank, dontsum),
            axis=1)
    if outp_file:
        predictions_filtered.groupby(by=["evaluation_taxa"], as_index=False).sum().to_csv(outp_file, index=False,
                                                                                          sep="\t",
                                                                                          columns=['evaluation_taxa',
                                                                                                   'reads',
                                                                                                   'abundance'],
                                                                                          mode='a')
    else:
        return predictions_filtered.groupby(by=["evaluation_taxa"], as_index=False).sum()


def load_truth(yaml_file, ref_genomes_files, config, rank, unknown_lineages):
    """
    TODO
    :param yaml_file:
    :param ref_genomes_files:
    :param taxonomy:
    :param rank:
    :param targets:
    :param unknown_lineages:
    :return:
    """
    taxonomy = config["targets_taxonomy"]
    targets = config["targets"]
    with open(yaml_file, "r") as file:
        sample_composition = yaml.safe_load(file)
    sample_composition_df = pd.DataFrame(columns=['organism', 'abundance'])
    to_ignore = [line.strip() for line in open(unknown_lineages)]
    for target in targets:
        for acc in sample_composition[target]:
            if list(acc.keys())[0] in to_ignore:
                continue
            sample_composition_df = sample_composition_df.append(
                pd.DataFrame([[list(acc.keys())[0], list(acc.values())[0]]], columns=['organism', 'abundance']))
    return sum_to_taxlevel(sample_composition_df, ref_genomes_files, taxonomy, rank, config)


def load_classified_reads_truth(yaml_file, targets):
    """
    TODO
    :param yaml_file:
    :param targets:
    :return:
    """
    with open(yaml_file, "r") as file:
        instance_composition = yaml.safe_load(file)
    # sample_composition_df = pd.DataFrame(columns=['organism', 'abundance'])
    #
    # for target in targets:
    #     for acc in sample_composition[target]:
    #         sample_composition_df = sample_composition_df.append(
    #             pd.DataFrame([[list(acc.keys())[0], list(acc.values())[0]]], columns=['organism', 'abundance']))
    # return sum_to_taxlevel(sample_composition_df, ref_genomes_files, taxonomy, rank)
    total_ratio = 0
    if 'prok' in instance_composition['targets']:
        total_ratio += instance_composition['prok_ratio']
    if 'euk' in instance_composition['targets']:
        total_ratio += instance_composition['euk_ratio']
    if 'vir' in instance_composition['targets']:
        total_ratio += instance_composition['vir_ratio']
    return total_ratio * instance_composition['total_reads_nb']


def evaluation_results(predictions, truth, threshold):
    """
    TODO
    :param predictions:
    :param truth:
    :param threshold:
    :return:
    """
    predictions = predictions[predictions["reads"] >= threshold]
    if len(predictions) == 0:
        return predictions
    if len(truth) == 0:
        return predictions
    predictions.loc[:, "evaluation"] = predictions.apply(
        lambda row: row["evaluation_taxa"]
                    in truth["evaluation_taxa"].values,
        axis=1,
    )
    return predictions


def get_precision(predictions, truth, threshold):
    """
    TODO
    :param predictions:
    :param truth:
    :param threshold:
    :return:
    """
    try:
        predictions = evaluation_results(predictions, truth, threshold)
        #print(predictions)
        precision = predictions[predictions["evaluation"]].evaluation.count() / predictions["evaluation"].count()
        if precision >= 0:
            return precision
        return 0
    except KeyError:
        return 0


def get_recall(predictions, truth, threshold):
    """
    TODO
    :param predictions:
    :param truth:
    :param threshold:
    :return:
    """
    try:
        predictions = evaluation_results(predictions, truth, threshold)
        recall = predictions[predictions["evaluation"]].evaluation.count() / truth["evaluation_taxa"].count()
        if recall >= 0:
            return recall
        return recall
    except KeyError:
        return 0


def get_tp(predictions, truth, threshold):
    """
    TODO
    :param predictions:
    :param truth:
    :param threshold:
    :return:
    """
    try:
        predictions = evaluation_results(predictions, truth, threshold)
        return predictions[predictions["evaluation"]].evaluation.count()
    except KeyError:
        return 0


def get_fp(predictions, truth, threshold):
    """
    TODO
    :param predictions:
    :param truth:
    :param threshold:
    :return:
    """
    try:
        predictions = evaluation_results(predictions, truth, threshold)
        return predictions[~predictions["evaluation"]].evaluation.count()
    except KeyError:
        return 0


def get_f1(predictions, truth, threshold):
    """
    TODO
    :param predictions:
    :param truth:
    :param threshold:
    :return:
    """
    try:
        predictions = evaluation_results(predictions, truth, threshold)
        precision = predictions[predictions["evaluation"]].evaluation.count() / predictions["evaluation"].count()
        recall = predictions[predictions["evaluation"]].evaluation.count() / truth["evaluation_taxa"].count()
        if (2 * (precision * recall) / (precision + recall)) >= 0:
            return 2 * (precision * recall) / (precision + recall)
        return 0
    except KeyError:
        return 0


def by_row_squared_distance(row, predictions):
    """
    TODO
    :param row:
    :param predictions:
    :return:
    """
    return (row['norm_abundance'] - predictions[predictions['evaluation_taxa'] == row["evaluation_taxa"]][
        'abundance']) ** 2


def get_l2(predictions, truth, threshold):
    """
    TODO
    :param predictions:
    :param truth:
    :param threshold:
    :return:
    """
    try:
        predictions = evaluation_results(predictions, truth, threshold)
        truth.loc[:, 'norm_abundance'] = truth.apply(lambda row: row['abundance'] / truth['abundance'].sum(), axis=1)
        squared_distances = truth.apply(lambda row: by_row_squared_distance(row, predictions), axis=1)
        return np.sqrt(squared_distances.sum().values[0])
    except ValueError:
        return -1
    except KeyError:
        return -1
    except IndexError:
        return -1


def to_json(inp, outp, config,  orient="columns", sort_with=False, reverse=False):
    """
    TODO
    :param inp:
    :param outp:
    :return:
    """
    data = pd.read_csv(inp, sep="\t", header=1)
    if sort_with == 'self':
        data = data.assign(sort_with=1)
        data_sort = data.loc[data['rank'] == config['rank_to_sort_by']]
        data_sort = data_sort.loc[data['sample'] == config['instance_name']]
        for tool in data_sort['toolname']:
            data.loc[data['toolname']==tool, 'sort_with'] = data_sort.loc[data['toolname']==tool, 'value'].values[0]
    elif sort_with:
        data = data.assign(sort_with=1)
        for line in open(sort_with):
            if line.startswith('#'):
                continue
            if line.split('\t')[1] == config['instance_name'] and line.split('\t')[2] == config['rank_to_sort_by']:
                value = float(line.strip().split('\t')[3])
                tool = line.split('\t')[0]
                data.loc[data['toolname']==tool, 'sort_with'] = value
    try:
        data.sort_values(by=['sort_with'], inplace=True, ignore_index=True, ascending=True*(not reverse))
        data.drop(["sort_with"], inplace=True, axis=1)
    except KeyError:
        pass
    data.to_json(
        outp, index=True, indent=4, orient=orient
    )
#CLARA16s functions
def getMax(row):
 maxGCN=1
 for search in ["domain_GCN","phylum_GCN","class_GCN","order_GCN","family_GCN","genus_GCN","species_GCN"]:
  if row[search] > maxGCN:
      maxGCN=row[search]
 return maxGCN

def load_classified_reads_predictions(gcn, predictionFile, taxo, output, gcn_normalization):
  predictions_reads = pd.read_csv(predictionFile, header=1, sep="\t")
  df2 = predictions_reads['Taxonomy'].str.replace('; ',';')
  df2 = df2.str.split(';',expand=True)
  df2.insert(0, "size", predictions_reads["Size"])
  df2.insert(0, "id", predictions_reads["Group_ID"])
  df2.columns = ['id',"size","domain","phylum","class","order","family","genus","species"]
  df2.replace('Unassigned','na__Unassigned', inplace=True)
  df2.replace(np.nan,'na__Unassigned', inplace=True)

  for search_tax_level in ["domain","phylum","class","order","family","genus","species"]:
     df = pd.DataFrame()
     df['id'] = df2[["id"]] 
     df['name'] = df2[[search_tax_level]]
     df['name'] = df['name'].str.split('__').str.get(1)
     df['name'] = df['name'].str.replace('_',' ')
     df = df.merge(gcn, how='left', on='name')
     df2[search_tax_level+'_GCN']=df["GCN_mean"]
  df2.fillna(0, inplace=True)
  #df2=df2.sort_values(by=['id'])

  sumatory={}
  gcn_value={}
  for values in df2[taxo].unique():
    sumatory[values]=0

  for index, row in df2.iterrows():
    values=row[taxo]
    reported_gcn=float(getMax(row)) if gcn_normalization=='yes' else 1
    sumatory[values]=sumatory[values]+(row["size"]/reported_gcn)
    gcn_value[values]=reported_gcn

  df_tool_abud = pd.DataFrame([sumatory])
  df_tool_abud = df_tool_abud.T
  df_tool_abud.reset_index(inplace=True)
  df_tool_abud.columns =['evaluation_taxa','reads']
  df_tool_abud['abundance']=df_tool_abud['reads']/df_tool_abud['reads'].sum()
  df_tool_abud['GCN_mean']= df_tool_abud['evaluation_taxa'].map(gcn_value)
  #print(df_tool_abud)

  #Savave the results
  file=open(output,'w+')
  toprint='#CLARA16s\n'
  file.write(toprint)
  df_tool_abud.to_csv(file, sep="\t", header=True, mode='a', index=True)
  file.close()

def load_truth_reads_predictions(gcn, predictionFile, taxo, gcn_normalization):
  predictions_reads = pd.read_csv(predictionFile, header=1, sep="\t")
  df2 = predictions_reads['Taxonomy'].str.replace('; ',';')
  df2 = df2.str.split(';',expand=True)
  df2.insert(0, "size", predictions_reads["Size"])
  df2.insert(0, "id", predictions_reads["Group_ID"])
  df2.columns = ['id',"size","domain","phylum","class","order","family","genus","species"]
  df2.replace('Unassigned','na__Unassigned', inplace=True)
  df2.replace(np.nan,'na__Unassigned', inplace=True)

  for search_tax_level in ["domain","phylum","class","order","family","genus","species"]:
     df = pd.DataFrame()
     df['id'] = df2[["id"]] 
     df['name'] = df2[[search_tax_level]]
     df['name'] = df['name'].str.split('__').str.get(1)
     df['name'] = df['name'].str.replace('_',' ')
     df = df.merge(gcn, how='left', on='name')
     df2[search_tax_level+'_GCN']=df["GCN_mean"]
  df2.fillna(0, inplace=True)
  #df2=df2.sort_values(by=['id'])

  sumatory={}
  gcn_value={}
  for values in df2[taxo].unique():
    sumatory[values]=0

  for index, row in df2.iterrows():
    values=row[taxo]
    reported_gcn=float(getMax(row)) if gcn_normalization=='yes' else 1
    sumatory[values]=sumatory[values]+(row["size"]/reported_gcn)
    gcn_value[values]=reported_gcn

  df_tool_abud = pd.DataFrame([sumatory])
  df_tool_abud = df_tool_abud.T
  df_tool_abud.reset_index(inplace=True)
  df_tool_abud.columns =['evaluation_taxa','reads']
  df_tool_abud['abundance']=df_tool_abud['reads']/df_tool_abud['reads'].sum()
  df_tool_abud['GCN_mean']= df_tool_abud['evaluation_taxa'].map(gcn_value)
  #print(df_tool_abud)
  return df_tool_abud

def GCN_reference(gnc_ref):
  gcn=pd.read_csv(str(gnc_ref), header=0, sep="\t")
  gcn = gcn[['taxid','name','mean']]
  gcn["mean"]=gcn["mean"].round().astype('int')
  gcn.columns = ['taxid','name','GCN_mean']
  #print(gcn)
  return gcn
