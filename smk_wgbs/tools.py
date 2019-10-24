import itertools
import os
import re
import glob
from pathlib import Path
from pkg_resources import resource_filename
import pandas as pd
import itertools
from pandas.api.types import CategoricalDtype
from typing import Optional, Dict


def sel_expand(template, **kwargs):
    fields = kwargs.keys()
    values = [kwargs[f] for f in fields]
    values = [[val] if isinstance(val, (int, str)) else val for val in values]
    value_combinations = itertools.product(*values)

    def get_expanded_template(template, fields, comb):
        for field, value in zip(fields, comb):
            template = template.replace("{" + field + "}", value)
        return template

    res = [get_expanded_template(template, fields, comb) for comb in value_combinations]
    if len(res) == 1:
        return res[0]
    return res


def create_metadata_table_from_file_pattern(
    wildcard_pattern: str, field_constraints: Optional[Dict] = None
):
    """Create metadata table for all files matching a snakemake-like pattern

    Output: metadata table with these columns:
      - one column per wildcard field
      - 'path' contains the full match for the snakemake-like pattern
      - fields in the output table are in the order of appearance from the filepath pattern

    Details:
      - fields may occur multiple times in the pattern
      - files are found by replacing each field with a '*' and using glob.glob
      - metadata are extracted using a regex constructed as follows:
        - the first occurence of each field is replaced with the default regex ('.+'),
          or the regex supplied via field_constraints
        - all following occurences of the field are replace with a backreference
          to the first match for the field


    """
    if field_constraints is None:
        field_constraints = {}

    field_names_set = set()
    all_field_name_occurences = re.findall(r"{(.*?)}", wildcard_pattern)
    field_names_in_order_of_appearance = [
        x
        for x in all_field_name_occurences
        if not (x in field_names_set or field_names_set.add(x))
    ]

    assert set(field_constraints.keys()) <= field_names_set

    glob_pattern = re.sub(r"{(.+?)}", r"*", wildcard_pattern)
    glob_results = glob.glob(glob_pattern)
    if not glob_results:
        raise ValueError(f"Could not find any file matching:\n{glob_pattern}")

    regex_pattern = wildcard_pattern
    for field_name in field_names_set:
        if field_name in field_constraints:
            # replace first field with regex
            regex_pattern = regex_pattern.replace(
                "{" + field_name + "}",
                f"(?P<{field_name}>{field_constraints[field_name]})",
                1,
            )
            # replace following fields with named backreference, if there are any
            regex_pattern = regex_pattern.replace(
                "{" + field_name + "}", f"(?P={field_name})"
            )
        else:
            regex_pattern = regex_pattern.replace(
                "{" + field_name + "}", f"(?P<{field_name}>.+)", 1
            )
            # replace following fields with named backreference, if there are any
            regex_pattern = regex_pattern.replace(
                "{" + field_name + "}", f"(?P={field_name})"
            )

    metadata_df = pd.Series(glob_results).str.extract(regex_pattern)
    metadata_df["path"] = glob_results
    metadata_df = metadata_df[["path"] + field_names_in_order_of_appearance]

    return metadata_df


def get_snakefile_fp() -> str:
    return resource_filename("smk_wgbs", "wgbs_alignment.smk")

def print_snakefile_fp():
    print(get_snakefile_fp())


def get_demo_config_fp():
    return resource_filename("smk_wgbs", "demo_config.yaml")

def print_demo_config_fp():
    print(get_demo_config_fp())


def find_workflow_version():
    return "v0.1.0"


def fp_to_parquet(fp, suffix="bed"):
    return re.sub(f".{suffix}$", ".parquet", fp)


index_pattern_cg = "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mcall_qc/sandbox/genomes/GRCm38mm10_PhiX_Lambda/GRCm38mm10_PhiX_Lambda_CG-nome_{chromosome}.bed.gz"
index_pattern_ch = "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mcall_qc/sandbox/genomes/GRCm38mm10_PhiX_Lambda/GRCm38mm10_PhiX_Lambda_CH-nome_{chromosome}.bed.gz"
index_pattern = "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mcall_qc/sandbox/genomes/GRCm38mm10_PhiX_Lambda/GRCm38mm10_PhiX_Lambda_CG-CHG-CHH_{chromosome}.bed.gz"
nome_index_pattern = "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mcall_qc/sandbox/genomes/GRCm38mm10_PhiX_Lambda/GRCm38mm10_PhiX_Lambda_CG-CHG-CHH_{chromosome}_NOMe-context.bed.gz"


def prepare_index_files():

    # fmt: off
    chromosomes = [
        "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
        "11", "12", "13", "14", "15", "16", "17", "18", "19",
        "X", "Y", "MT", "phix", "L",
    ]
    # fmt: on
    chrom_dtype = CategoricalDtype(categories=chromosomes, ordered=True)
    strand_dtype = CategoricalDtype(categories=["+", "-"], ordered=True)
    motif_dtype = pd.CategoricalDtype(["CG", "CHG", "CHH"], ordered=True)
    seq_contexts_l = sorted(
        "".join(t)
        for t in itertools.product(
            list("ACGT"), list("ACGT"), ["C"], list("ACGT"), list("ACGT")
        )
    )
    seq_context_dtype = pd.CategoricalDtype(seq_contexts_l, ordered=True)

    seq_context_to_nomeseq_d = {seq: seq[1:3] + ('G' if seq[3] is 'G' else 'H')
                                for seq in seq_contexts_l}


    nome_seq_context_l = sorted(
        "".join(t) for t in itertools.product(list("ACGT"), "C", ["G", "H"])
    )
    assert set(seq_context_to_nomeseq_d.values()) == set(nome_seq_context_l)
    # nome_seq_context_dtype = pd.CategoricalDtype(nome_seq_context_l, ordered=True)

    for chromosome in chromosomes:
        print(chromosome)

        index_df = pd.read_csv(
            index_pattern.format(chromosome=chromosome),
            sep="\t",
            # fmt: off
            names=["Chromosome", "Start", "End", "motif", "score", "strand", "seq_context"],
            usecols=["Chromosome", "Start", "End", "motif", "strand", "seq_context"],
            header=0,
            # fmt: on
            dtype={
                "Chromosome": chrom_dtype,
                "strand": strand_dtype,
                "seq_context": seq_context_dtype,
                "motif": motif_dtype,
            },
        )

        index_df['nome_seq_context'] = index_df['seq_context'].str[1:4]

        index_df = index_df[
            [
                "Chromosome",
                "Start",
                "End",
                "motif",
                "strand",
                "seq_context",
                "nome_seq_context",
            ]
        ]

        out_bed = nome_index_pattern.format(chromosome=chromosome)
        # csv writer too slow
        # index_df.rename(columns={"Chromosome": "#Chromosome"}).to_csv(
        #     out_bed, sep="\t", index=False
        # )
        index_df.to_parquet(fp_to_parquet(out_bed))


def nome_filtering(input, output, params):

    print('start nome filtering')

    # class input:
    #     CG = '/icgc/dkfzlsdf/analysis/hs_ontogeny/results/scnmt/results-per-entity/lsk150_pr-hiera_ex-scnmt-1_cls-1/blood67/smk-wgbs/default-params_v0.1.0/meth-calls/partial-meth-calls/SE_lsk150_pr-hiera_ex-scnmt-1_cls-1_blood67_per-cytosine_CG.parquet'
    #     CHG = '/icgc/dkfzlsdf/analysis/hs_ontogeny/results/scnmt/results-per-entity/lsk150_pr-hiera_ex-scnmt-1_cls-1/blood67/smk-wgbs/default-params_v0.1.0/meth-calls/partial-meth-calls/SE_lsk150_pr-hiera_ex-scnmt-1_cls-1_blood67_per-cytosine_CHG.parquet'
    #     CHH = '/icgc/dkfzlsdf/analysis/hs_ontogeny/results/scnmt/results-per-entity/lsk150_pr-hiera_ex-scnmt-1_cls-1/blood67/smk-wgbs/default-params_v0.1.0/meth-calls/partial-meth-calls/SE_lsk150_pr-hiera_ex-scnmt-1_cls-1_blood67_per-cytosine_CHH.parquet'

    columns = ["Chromosome", "Start", "End", "beta_value", "n_meth", "n_total"]
    cg_df = pd.read_parquet(input.CG, columns=columns)
    chg_df = pd.read_parquet(input.CHG, columns=columns)
    chh_df = pd.read_parquet(input.CHH, columns=columns)

    cg_chrom_dfs = []
    ch_chrom_dfs = []
    for chromosome in params.chromosomes:
        print('Chromosome', chromosome)
        print('CG merge')
        chrom_index = pd.read_parquet(
            nome_index_pattern.format(chromosome=chromosome),
                columns=[
                "Chromosome",
                "Start",
                "End",
                "motif",
                "seq_context",
                "nome_seq_context",
            ],
        )
        cg_chrom_df = cg_df.loc[cg_df["Chromosome"] == chromosome]
        cg_chrom_df_annotated = pd.merge(cg_chrom_df,
                                         chrom_index.loc[chrom_index['motif'] == 'CG'],
                                         how="inner")
        if cg_chrom_df.shape[0] != cg_chrom_df_annotated.shape[0]:
            print(f'WARNING: cg_calls {cg_chrom_df.shape[0]}, inner merge: {cg_chrom_df_annotated.shape[0]}')
        if not cg_chrom_df_annotated.empty:
            cg_chrom_dfs.append(cg_chrom_df_annotated)

        print('CH merge')
        ch_chrom_df = pd.concat(
            [
                chg_df.loc[chg_df["Chromosome"] == chromosome],
                chh_df.loc[chh_df["Chromosome"] == chromosome],
            ],
            axis=0,
        ).sort_values(["Chromosome", "Start", "End"])
        ch_chrom_df_annotated = pd.merge(
                ch_chrom_df,
                chrom_index.loc[chrom_index['motif'] != 'CG'],
                how="inner")
        if ch_chrom_df.shape[0] != ch_chrom_df_annotated.shape[0]:
            print(f'WARNING: ch_calls {ch_chrom_df.shape[0]}, inner merge: {ch_chrom_df_annotated.shape[0]}')
        if not ch_chrom_df_annotated.empty:
            ch_chrom_dfs.append(ch_chrom_df_annotated)

    print('Saving results')
    cg_full_df = pd.concat(cg_chrom_dfs, axis=0)
    cg_full_df.to_parquet(output.cg_full_parquet)

    cg_meth_df = cg_full_df.query('nome_seq_context in ["ACG", "TCG"]').reset_index()
    cg_meth_df.to_parquet(output.cg_meth_parquet)

    ch_full_df = pd.concat(ch_chrom_dfs, axis=0)
    ch_full_df.to_parquet(output.ch_full_parquet)

    ch_acc_df = ch_full_df.query('nome_seq_context in ["GCA", "GCT", "GCC"]')
    ch_acc_df.to_parquet(output.ch_acc_parquet)

    ch_meth_df = ch_full_df.query('nome_seq_context in ["ACA", "ACT", "ACC", "TCA", "TCT", "TCC"]')

    chg_meth_df = ch_meth_df.query('motif == "CHG"')
    chg_meth_df.to_parquet(output.chg_meth_parquet)

    chh_meth_df = ch_meth_df.query('motif == "CHH"')
    chh_meth_df.to_parquet(output.chh_meth_parquet)
