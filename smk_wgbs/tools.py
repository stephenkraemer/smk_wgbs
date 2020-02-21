import glob
from collections import defaultdict
import pickle
import itertools
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Dict
from snakemake.io import expand

import pandas as pd
import numpy as np
from pandas.api.types import CategoricalDtype
from pkg_resources import resource_filename


def to_pickle(obj, fp, protocol=4):
    with open(fp, "wb") as fout:
        pickle.dump(obj, fout, protocol=protocol)


def from_pickle(fp):
    with open(fp, "rb") as fin:
        return pickle.load(fin)


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


def fp_to_tsv(fp, suffix="bed"):
    return re.sub(f".{suffix}$", ".tsv", fp)


def fp_to_pickle(fp, suffix="bed"):
    return re.sub(f".{suffix}$", ".p", fp)


def fp_to_bedgraph(fp, suffix="bed"):
    return re.sub(f".{suffix}$", ".bedGraph", fp)


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

    seq_context_to_nomeseq_d = {
        seq: seq[1:3] + ("G" if seq[3] is "G" else "H") for seq in seq_contexts_l
    }

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

        index_df["nome_seq_context"] = index_df["seq_context"].str[1:4]

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
        index_df.to_pickle(fp_to_pickle(out_bed))


def nome_filtering(input, output, params):

    print("start nome filtering")

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
        print("Chromosome", chromosome)
        print("CG merge")

        # TEST

        chrom_index = pd.read_csv(
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

        chrom_index = pd.read_pickle(
            nome_index_pattern.format(chromosome=chromosome) + ".p"
        )

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
            engine="pyarrow",
        )

        # /TEST

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
        assert cg_chrom_df["Chromosome"].eq(chromosome).all()
        assert chrom_index["Chromosome"].eq(chromosome).all()
        cg_chrom_df_annotated = pd.merge(
            cg_chrom_df,
            chrom_index.loc[chrom_index["motif"] == "CG"],
            on=["Start", "End"],
            how="left",
        )
        if cg_chrom_df.shape[0] != cg_chrom_df_annotated.shape[0]:
            print(
                f"WARNING: cg_calls {cg_chrom_df.shape[0]}, inner merge: {cg_chrom_df_annotated.shape[0]}"
            )
        if not cg_chrom_df_annotated.empty:
            cg_chrom_dfs.append(cg_chrom_df_annotated)

        print("CH merge")
        ch_chrom_df = pd.concat(
            [
                chg_df.loc[chg_df["Chromosome"] == chromosome],
                chh_df.loc[chh_df["Chromosome"] == chromosome],
            ],
            axis=0,
        ).sort_values(["Chromosome", "Start", "End"])
        assert ch_chrom_df["Chromosome"].eq(chromosome).all()

        # TEST
        chrom_index = chrom_index.loc[chrom_index["motif"] != "CG"][
            ["Start", "End"]
        ].copy()
        chrom_index["value"] = np.arange(len(chrom_index))
        ch_chrom_df = ch_chrom_df.drop("Chromosome", axis=1)

        ch_chrom_df_annotated = pd.merge(
            ch_chrom_df,
            chrom_index.loc[chrom_index["motif"] != "CG"],
            on=["Start", "End"],
            how="left",
        )
        # \ TEST

        ch_chrom_df_annotated = pd.merge(
            ch_chrom_df,
            chrom_index.loc[chrom_index["motif"] != "CG"],
            on=["Start", "End"],
            how="left",
        )

        if ch_chrom_df.shape[0] != ch_chrom_df_annotated.shape[0]:
            print(
                f"WARNING: ch_calls {ch_chrom_df.shape[0]}, inner merge: {ch_chrom_df_annotated.shape[0]}"
            )
        if not ch_chrom_df_annotated.empty:
            ch_chrom_dfs.append(ch_chrom_df_annotated)

    print("Saving results")
    cg_full_df = pd.concat(cg_chrom_dfs, axis=0)
    cg_full_df.to_parquet(output.cg_full_parquet)

    cg_meth_df = cg_full_df.query('nome_seq_context in ["ACG", "TCG"]').reset_index()
    cg_meth_df.to_parquet(output.cg_meth_parquet)

    ch_full_df = pd.concat(ch_chrom_dfs, axis=0)
    ch_full_df.to_parquet(output.ch_full_parquet)

    ch_acc_df = ch_full_df.query('nome_seq_context in ["GCA", "GCT", "GCC"]')
    ch_acc_df.to_parquet(output.ch_acc_parquet)

    ch_meth_df = ch_full_df.query(
        'nome_seq_context in ["ACA", "ACT", "ACC", "TCA", "TCT", "TCC"]'
    )

    chg_meth_df = ch_meth_df.query('motif == "CHG"')
    chg_meth_df.to_parquet(output.chg_meth_parquet)

    chh_meth_df = ch_meth_df.query('motif == "CHH"')
    chh_meth_df.to_parquet(output.chh_meth_parquet)


def run_methyldackel(input, output, params, threads):

    print("Running MethylDackel")
    subprocess.run(
        f"""
            MethylDackel extract \
            {params.calling_options} \
            {params.motif_params} \
            {params.sampling_options} \
            -o {params.prefix} \
            -@ {threads} \
            {input.ref_genome_unconverted} \
            {input.bam}
            """.split()
    )
    print("Done with MethylDackel")

    chrom_dtype = CategoricalDtype(categories=params.chromosomes, ordered=True)

    assert len(output) % 5 == 0
    for i in range(int(len(output) / 4)):
        bed, bedgraph, parquet, pickle, tmp_bedgraph = output[5 * i : (5 * i) + 5]
        # convert bedgraph to bed and parquet
        shutil.copy(tmp_bedgraph, bedgraph)
        methyldackel_bedgraph_to_bed_and_parquet(
            bed, bedgraph, chrom_dtype, parquet, pickle
        )


def methyldackel_bedgraph_to_bed_and_parquet(
    bed, bedgraph, chrom_dtype, parquet, pickle
):
    print("Create parquet and BED file")
    df = pd.read_csv(
        bedgraph,
        sep="\t",
        skiprows=1,
        names=["#Chromosome", "Start", "End", "beta_value", "n_meth", "n_unmeth"],
        dtype={"#Chromosome": chrom_dtype},
    )
    # Discard unwanted chromosomes/or more commonly: contigs (not in config['chromosomes'])
    df = df.dropna(subset=["#Chromosome"])
    df["n_total"] = df.eval("n_meth + n_unmeth")
    df["beta_value"] = df.eval("n_meth / n_total")
    # Add motif column
    # TODO: improve
    if "CG" in bed:
        motif = "CG"
    elif "CHG" in bed:
        motif = "CHG"
    elif "CHH" in bed:
        motif = "CHH"
    else:
        motif = "NA"
    df["motif"] = motif
    # Final column order
    df = df[["#Chromosome", "Start", "End", "motif", "beta_value", "n_total", "n_meth"]]
    # Save to BED and parquet
    df.to_csv(bed, sep="\t", header=True, index=False)
    df.rename(columns={"#Chromosome": "Chromosome"}).to_parquet(parquet)
    df.rename(columns={"#Chromosome": "Chromosome"}).to_pickle(pickle)
    print("done")


def run_bismark_alignment(input, output, log, params, is_paired, upto=-1):
    # temp_dir
    # - need different temp_dir per atomic alignment
    # - rule removes the atomic temp_dir at the end of the rule using the temp()
    # directive. This is more robust than removing the directory inside the rule
    # (in case of crashed jobs)
    # - if bismark changes its behavior to completely remove the temp_dir at the end
    # of its run, this rule will fail, because the temp_dir is an expected (temp) output

    output_dir = Path(output.bam).parent
    if is_paired:
        input_option = f"-1 {input.fq_r1} -2 {input.fq_r2}"
        unmapped_option = "--unmapped"
        stem = re.sub(".fastq(.gz)?$", "", Path(input.fq_r1).name)
    else:
        input_option = f"--single_end {input.fq}"
        unmapped_option = ""
        stem = re.sub(".fastq(.gz)?$", "", Path(input.fq).name)
    if upto > 0:
        upto_option = f"--upto {upto}"
        print(f"WARNING: only performing the first {upto} alignments")
    else:
        upto_option = ""
    subprocess.run(
        f"""
            bismark \
            {input_option} \
            {unmapped_option} \
            {upto_option} \
            --non_directional \
            --local \
            --output_dir {output_dir} \
            --temp_dir {output.temp_dir}\
            {params.genome_dir} \
            > {log} 2>&1
            """.split()
    )

    # obtain paths where bismark will place files
    # bismark hardcodes output files basenames
    # (--basename cannot be used together with multithreading atm)

    if is_paired:
        report_name = stem + "_bismark_bt2_PE_report.txt"
        bismark_bam_fp = output_dir.joinpath(stem + "_bismark_bt2_pe.bam")
    else:
        report_name = stem + "_bismark_bt2_SE_report.txt"
        bismark_bam_fp = output_dir.joinpath(stem + "_bismark_bt2.bam")
    bismark_report_fp = output_dir.joinpath(report_name)
    # Move results from hard coded paths to user-specified paths
    bismark_bam_fp.rename(output.bam)
    bismark_report_fp.rename(output.bismark_report)
    if is_paired:

        bismark_unmapped_r1 = output_dir.joinpath(
            Path(input.fq_r1).name + "_unmapped_reads_1.fq.gz"
        )
        bismark_unmapped_r1.rename(output.unmapped_fq1)

        bismark_unmapped_r2 = output_dir.joinpath(
            Path(input.fq_r2).name + "_unmapped_reads_2.fq.gz"
        )
        bismark_unmapped_r2.rename(output.unmapped_fq2)


def merge_or_symlink(input, output):
    # hard link if only one bam, else merge
    if len(input) == 1:
        os.link(input[0], output.bam)
        os.link(input[0] + ".bai", output.bam + ".bai")
    else:
        subprocess.run("samtools merge -f -r".split() + output + input)
        subprocess.run(["samtools", "index", output.bam])


def find_fastqs(wildcards, metadata_table):
    df = metadata_table
    is_in_selected_pair = (
        (df["protocol"] == wildcards.protocol)
        & (df["entity"] == wildcards.entity)
        & (df["sample"] == wildcards.sample)
        & (df["lib"] == wildcards.lib)
        & (df["uid"] == wildcards.uid)
    )
    fastq_r1_ser = df.loc[is_in_selected_pair & (df["read_number"] == "1"), "path"]
    assert fastq_r1_ser.shape[0] == 1
    fastq_r2_ser = df.loc[is_in_selected_pair & (df["read_number"] == "2"), "path"]
    assert fastq_r2_ser.shape[0] == 1
    return {"fq_r1": fastq_r1_ser.iloc[0], "fq_r2": fastq_r2_ser.iloc[0]}


def find_individual_fastqs(wildcards, metadata_table):
    df = metadata_table
    is_in_selected_pair = (
        (df["protocol"] == wildcards.protocol)
        & (df["entity"] == wildcards.entity)
        & (df["sample"] == wildcards.sample)
        & (df["lib"] == wildcards.lib)
        & (df["uid"] == wildcards.uid)
        & (df["read_number"] == wildcards.read_number)
    )
    fastq_ser = df.loc[is_in_selected_pair, "path"]
    assert fastq_ser.shape[0] == 1
    return fastq_ser.iloc[0]


def merge_mcalls(input, output, wildcards):

    if len(input) == 1:
        input_parquet = Path(input[0])
        os.link(input_parquet, output.parquet)
        os.link(input_parquet.with_suffix(".bed"), output.bed)
        os.link(input_parquet.with_suffix(".bedGraph"), output.bedgraph)
        os.link(input_parquet.with_suffix(".p"), output.pickle)
        return

    assert len(input) == 2
    df1 = pd.read_pickle(Path(input[0]).with_suffix(".p"))
    assert df1["motif"].nunique() == 1
    df2 = pd.read_pickle(Path(input[1]).with_suffix(".p"))
    assert df2["motif"].nunique() == 1
    # outer leads to lexicographic sort, which ignores categorical order as of 0.25.1
    res = pd.merge(df1, df2, on=["Chromosome", "Start", "End"], how="outer")
    res = res.sort_values(["Chromosome", "Start", "End"])
    res["n_meth"] = res["n_meth_x"].add(res["n_meth_y"], fill_value=0)
    res["n_total"] = res["n_total_x"].add(res["n_total_y"], fill_value=0)
    res["beta_value"] = res["n_meth"] / res["n_total"]
    res["motif"] = df1["motif"].iloc[0]
    res = res[
        ["Chromosome", "Start", "End", "motif", "beta_value", "n_total", "n_meth"]
    ]
    assert res.isna().sum().sum() == 0
    res.to_parquet(output.parquet)
    res.to_pickle(output.pickle)
    res.rename(columns={"Chromosome": "#Chromosome"}).to_csv(
        output.bed, sep="\t", header=True, index=False
    )
    res["n_unmeth"] = res["n_total"] - res["n_meth"]
    with open(output.bedgraph, "wt") as fout:
        fout.write(
            f'track type="bedGraph" description="{wildcards.entity}_{wildcards.sample}"'
        )
    res[["Chromosome", "Start", "End", "beta_value", "n_meth", "n_unmeth"]].to_csv(
        output.bedgraph, mode="a", sep="\t", header=False, index=False
    )


def metadata_table_from_fastq_pattern(config):
    global_samples = config.get("samples", [])
    metadata_tables = []
    assert len(config["entities"]) > 0
    for elem in config["entities"]:
        if isinstance(elem, str):
            entity = elem
            samples = global_samples
        else:
            entity, *samples = elem

        # create metadata table per entity, using only specified entities
        if samples:
            for sample in samples:
                curr_metadata_table = create_metadata_table_from_file_pattern(
                    sel_expand(config["fastq_pattern"], entity=entity, sample=sample)
                ).assign(entity=entity, sample=sample)
                metadata_tables.append(curr_metadata_table)
        else:
            curr_metadata_table = create_metadata_table_from_file_pattern(
                sel_expand(config["fastq_pattern"], entity=entity)
            ).assign(entity=entity)
            metadata_tables.append(curr_metadata_table)

    metadata_table = pd.concat(metadata_tables, axis=0)

    # fill in optional pattern fields
    if "lib" not in metadata_table:
        metadata_table["lib"] = "1"
    if "protocol" not in metadata_table:
        metadata_table["protocol"] = config["protocol"]

    metadata_table = metadata_table.astype(str)

    # logging of metadata table not yet implemented
    # timestamp = time.strftime('%d-%m-%y_%H:%M:%S')
    # metadata_table_tsv = Path(config['metadata_table_tsv']).with_suffix(f'.{timestamp}.tsv')
    # Path(metadata_table_tsv).parent.mkdir(parents=True, exist_ok=True)
    # metadata_table.to_csv(metadata_table_tsv, header=True, index=False, sep='\t')

    return metadata_table


def test_meth_calling_qc_stats(config, sampling_name, wildcards, keys, names):

    import yaml

    config_yaml = get_demo_config_fp()
    with open(config_yaml) as fin:
        config = yaml.load(fin, Loader=yaml.FullLoader)
    config["protocol"] = "sctm"

    # Fill prefix with alignment config_name and workflow_version
    results_prefix = sel_expand(
        config["results_prefix"],
        config_name=config["alignment_config_name"],
        workflow_version=find_workflow_version(),
    )

    # concatenate paths
    config["log_dir"] = os.path.join(
        config["results_dir"], results_prefix, config["log_dir"]
    )
    for pattern_name, pattern in config["result_patterns"].items():
        pattern = pattern.replace(
            "{mcalling_config_name}", config["mcalling_config_name"]
        )
        config["result_patterns"][pattern_name] = os.path.join(
            config["results_dir"], results_prefix, pattern
        )

    class wildcards:
        entity = "lsk150_pr-hiera_ex-sctm-1_cls-1"
        sample = "blood1"

    sampling_name = "_unused"

    metadata_table_with_wildcards = create_mcalls_metadata_table(config, sampling_name)
    metadata_table = expand_mcalls_metadata_table(
        metadata_table_with_wildcards, wildcards.entity, wildcards.sample
    )
    # print(metadata_table.fp.iloc[0])
    compute_meth_calling_qc_stats(
        metadata_table=metadata_table,
        chrom_regions=config["chrom_regions"],
        prefix="/home/kraemers/temp/test_f",
        keys=[wildcards.entity, wildcards.sample],
        names=["entity", "sample"],
    )


def create_mcalls_metadata_table(config, sampling_name):

    # TODO-protocol
    patterns = {
        "SE": {
            "merged": sel_expand(
                config["result_patterns"]["pe_or_se_mcalls_merged_per_motif"],
                pairing="SE",
                protocol=config["protocol"],
            ),
            "single": sel_expand(
                config["result_patterns"]["pe_or_se_mcalls_cytosines_per_motif"],
                pairing="SE",
                protocol=config["protocol"],
            ),
        },
        "PE": {
            "merged": sel_expand(
                config["result_patterns"]["pe_or_se_mcalls_merged_per_motif"],
                pairing="PE",
                protocol=config["protocol"],
            ),
            "single": sel_expand(
                config["result_patterns"]["pe_or_se_mcalls_cytosines_per_motif"],
                pairing="PE",
                protocol=config["protocol"],
            ),
        },
        "SE_PE": {
            "merged": sel_expand(
                config["result_patterns"]["final_mcalls_merged_per_motif"],
                protocol=config["protocol"],
                alignment_combi="SE_PE",
            ),
            "single": sel_expand(
                config["result_patterns"]["final_mcalls_cytosines_per_motif"],
                protocol=config["protocol"],
                alignment_combi="SE_PE",
            ),
        },
    }
    for alignment_mode, d1 in patterns.items():
        for merge_mode, d2 in d1.items():
            d1[merge_mode] = fp_to_pickle(d2)

    alignment_combis = (
        ["PE", "SE", "SE_PE"]
        if config["alignment_combi"] == "SE_PE"
        else [config["alignment_combi"]]
    )
    input_rows = []
    for alignment_mode in alignment_combis:
        for motif, user_request in config["protocol_config"][config["protocol"]][
            "motifs"
        ].items():
            if not user_request:
                continue
            if user_request == "full":
                region_str = ""
            else:
                region_str = sampling_name
            for merging_mode in ["merged", "single"]:
                if merging_mode == "merged" and motif == "CHH":
                    continue
                input_rows.append(
                    dict(
                        pairing=alignment_mode,
                        merging_mode=merging_mode,
                        motif=motif,
                        fp=sel_expand(
                            patterns[alignment_mode][merging_mode],
                            motif=motif,
                            region=region_str,
                        ),
                    )
                )
    metadata_table = pd.DataFrame(input_rows)

    return metadata_table


def expand_mcalls_metadata_table(metadata_table, entity, sample):
    metadata_table = metadata_table.copy()
    metadata_table["fp"] = metadata_table["fp"].apply(
        lambda fp: expand(fp, entity=entity, sample=sample)[0]
    )
    return metadata_table


def save_df(df, out_p, keys, names, add_as="data"):
    df = df.copy()
    if add_as == "data":
        data_cols = df.columns.to_list()
        for name, key in zip(names, keys):
            df[name] = key
        df = df[names + data_cols]
    elif add_as == "index":
        df = pd.concat([df], keys=[tuple(keys)], names=names)
    else:
        raise ValueError
    df.to_pickle(out_p)
    df.to_csv(fp_to_tsv(out_p, "p"), sep="\t", header=True, index=False)
    df.columns = df.columns.astype(str)
    df.to_parquet(fp_to_parquet(out_p, "p"))


def compute_meth_calling_qc_stats(metadata_table, chrom_regions, prefix, keys, names):
    coverage_df_rows = []
    beta_value_dist_index_keys = []
    metadata_field_names = metadata_table.columns[0:-1].to_list()
    beta_value_dist_index_names = metadata_field_names + "is_corrected region".split()
    beta_value_dist_dfs = []
    global_meth_df_rows = []
    coverage_dist_dfs = []
    coverage_dist_index_keys = []
    coverage_dist_index_names = metadata_field_names + ["region"]
    for _unused, row_ser in metadata_table.iterrows():
        metadata_ser = row_ser.drop("fp")
        df = pd.read_pickle(row_ser["fp"])

        coverage_df_rows.append(
            dict(**metadata_ser, region="global", n_covered=df.shape[0])
        )

        df["n_total_capped"] = df["n_total"].mask(lambda ser: ser.gt(30), 30)
        coverage_dist_dfs.append(
            df["n_total_capped"]
            .value_counts()
            .reindex(np.arange(1, 31))
            .to_frame("frequency")
            .rename_axis("n_total")
        )
        coverage_dist_index_keys.append((*metadata_ser.to_list(), "global"))

        for region_name, chrom_or_chroms in chrom_regions.items():
            if isinstance(chrom_or_chroms, str):
                chrom_df = df.query("Chromosome == @chrom_or_chroms").copy()
            else:
                chrom_df = df.query("Chromosome in @chrom_or_chroms").copy()
            chrom_df["beta_value_corrected"] = chrom_df["beta_value"].mask(
                lambda ser: ser.between(0, 1, inclusive=False)
            )
            for is_corrected, beta_value_col in zip(
                [False, True], ["beta_value", "beta_value_corrected"]
            ):
                if chrom_df.empty:
                    if isinstance(chrom_or_chroms, str):
                        assert chrom_or_chroms in df["Chromosome"].cat.categories
                    else:
                        for c in chrom_or_chroms:
                            assert c in df["Chromosomes"].cat.categories
                    continue
                # print(region_name, chrom_or_chroms, beta_value_col)
                beta_value_dist_dfs.append(
                    pd.cut(
                        chrom_df[beta_value_col],
                        bins=np.linspace(0, 1.01, 102),
                        right=False,
                        include_lowest=False,
                    )
                    .value_counts()
                    .sort_index()
                    .to_frame("frequency")
                    .T
                )
                beta_value_dist_index_keys.append(
                    (*metadata_ser.to_list(), is_corrected, region_name)
                )

                if is_corrected:
                    n_meth = chrom_df[beta_value_col].eq(1).sum()
                    n_total = chrom_df[beta_value_col].notnull().sum()
                else:
                    n_meth = chrom_df["n_meth"].sum()
                    n_total = chrom_df["n_total"].sum()
                beta_value = n_meth / n_total
                n_intermediate_beta = (
                    chrom_df[beta_value_col].between(0, 1, inclusive=False).sum()
                )

                global_meth_df_rows.append(
                    dict(
                        **metadata_ser,
                        region=region_name,
                        is_corrected=is_corrected,
                        n_meth=n_meth,
                        n_total=n_total,
                        n_intermediate_beta=n_intermediate_beta,
                        beta_value=beta_value,
                    )
                )

    coverage_df = pd.DataFrame(coverage_df_rows)
    merge_motifs_groupby_cols = metadata_table.columns.drop(
        ["motif", "fp"]
    ).to_list() + ["region"]
    new = (
        coverage_df.query("not merging_mode == 'merged'")
        .groupby(merge_motifs_groupby_cols)
        .sum()
        .assign(motif="C")
        .reset_index()
    )
    coverage_df = pd.concat([coverage_df, new], sort=False).reset_index(drop=True)
    save_df(coverage_df, prefix + "_coverage-long.p", keys=keys, names=names)

    coverage_df_wide = (
        coverage_df.pivot_table(
            index=coverage_df.columns.drop(["pairing", "n_covered"]).to_list(),
            columns=["pairing"],
            values="n_covered",
        )
        .assign(
            overlap=lambda df: df.eval("SE + PE - SE_PE"),
            overlap_perc=lambda df: df.eval("overlap / SE_PE"),
        )
        .reset_index()
    )
    save_df(coverage_df_wide, prefix + "_coverage-wide.p", keys=keys, names=names)

    global_meth_df = pd.DataFrame(global_meth_df_rows)
    save_df(global_meth_df, prefix + "_global-meth.p", keys=keys, names=names)

    beta_value_dist_df = pd.concat(
        beta_value_dist_dfs,
        keys=beta_value_dist_index_keys,
        names=beta_value_dist_index_names,
    )
    beta_value_dist_df.index = beta_value_dist_df.index.droplevel(-1)
    save_df(
        beta_value_dist_df,
        prefix + "_beta-value-dist.p",
        keys=keys,
        names=names,
        add_as="index",
    )

    motif_coverage_dist_df = pd.concat(
        coverage_dist_dfs,
        keys=coverage_dist_index_keys,
        names=coverage_dist_index_names,
    )
    merge_motifs_groupby_cols = metadata_table.columns.drop(
        ["motif", "fp"]
    ).to_list() + ["region", "n_total"]
    agg_motifs_coverage_df = pd.concat(
        [motif_coverage_dist_df.groupby(merge_motifs_groupby_cols).sum()],
        keys=["C"],
        names=["motif"],
    )
    coverage_dist_df = pd.concat(
        [motif_coverage_dist_df.reset_index(), agg_motifs_coverage_df.reset_index()],
        sort=False,
    ).reset_index(drop=True)
    save_df(coverage_dist_df, prefix + "_coverage-dist.p", keys=keys, names=names)


def concat_meth_calling_qc_data(metadata_table, prefix, output):
    metadata_table = metadata_table[["entity", "sample"]].drop_duplicates()

    def concat_and_save(suffix, out_p):
        print(suffix)
        print(out_p)
        fps_ser = metadata_table.apply(
            lambda ser: expand(prefix + suffix, **ser)[0], axis=1
        )
        concat_df = pd.concat([pd.read_pickle(fp) for fp in fps_ser])
        concat_df.to_pickle(out_p)
        concat_df.to_csv(fp_to_tsv(out_p, "p"), sep="\t", header=True, index=False)
        concat_df.columns = concat_df.columns.astype(str)
        concat_df.to_parquet(fp_to_parquet(out_p, "p"))

    concat_and_save("_coverage-long.p", output.meth_calling_exp_qc_coverage_long)
    concat_and_save("_coverage-wide.p", output.meth_calling_exp_qc_coverage_wide)
    concat_and_save("_global-meth.p", output.meth_calling_exp_qc_global_meth)
    concat_and_save("_beta-value-dist.p", output.meth_calling_exp_qc_beta_value_dist)
    concat_and_save("_coverage-dist.p", output.meth_calling_exp_qc_coverage_dist)


def fn():

    import seaborn as sns
    import pandas as pd

    coverage_long_p = "/icgc/dkfzlsdf/analysis/hs_ontogeny/results/sctm/cohort-results/experiment-qc/sctm-1-meth-calling-qc_coverage-long.p"
    coverage_wide_p = "/icgc/dkfzlsdf/analysis/hs_ontogeny/results/sctm/cohort-results/experiment-qc/sctm-1-meth-calling-qc_coverage-wide.p"
    global_meth_p = "/icgc/dkfzlsdf/analysis/hs_ontogeny/results/sctm/cohort-results/experiment-qc/sctm-1-meth-calling-qc_global-meth.p"
    beta_dist_p = "/icgc/dkfzlsdf/analysis/hs_ontogeny/results/sctm/cohort-results/experiment-qc/sctm-1-meth-calling-qc_beta-value-dist.p"
    coverage_dist_p = "/icgc/dkfzlsdf/analysis/hs_ontogeny/results/sctm/cohort-results/experiment-qc/sctm-1-meth-calling-qc_coverage-dist.p"

    # DF // 'alignment_config_name', 'mcalling_config_name', 'entity', 'sample',
    #        'pairing', 'merging_mode', 'motif', 'region', 'n_covered'
    coverage_long_df = pd.read_pickle(coverage_long_p)

    # DF // 'alignment_config_name', 'mcalling_config_name', 'entity', 'sample',
    #        'merging_mode', 'motif', 'region', 'PE', 'SE', 'SE_PE', 'overlap',
    #        'overlap_perc'
    coverage_wide_df = pd.read_pickle(coverage_wide_p)

    global_meth_df = pd.read_pickle(global_meth_p)
    beta_dist_df = pd.read_pickle(beta_dist_p)
    coverage_dist_p = pd.read_pickle(coverage_dist_p)


def parse_samtools_stats(fp, keys):

    stats_in = defaultdict(list)
    with open(fp) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            stats_in[fields[0]].append(fields[1:])
    dfs = {
        "IS": {
            "columns": [
                "insert_size",
                "pairs_total",
                "inward_oriented_pairs",
                "outward_oriented_pairs",
                "other_pairs",
            ],
            "key": "insert_size",
        }
    }
    stats_out = {}
    summary_stats_l = stats_in.pop("SN")

    df = (
        pd.DataFrame(
            [
                (name[:-1], int(n) if "." not in n else float(n))
                for (name, n, *comment) in summary_stats_l
            ],
            columns=["name", "count"],
        )
        .assign(**keys)[list(keys.keys()) + ["name", "count"]]
        .pivot_table(index=list(keys.keys()), columns="name", values="count")
        .reset_index(drop=False)
    )
    df.columns.name = None

    stats_out["summary_stats"] = df

    for k, v in stats_in.items():
        stat_df = pd.DataFrame(v)
        if k in dfs:
            stat_df.columns = dfs[k]["columns"]
            stat_name = dfs[k]["key"]
        else:
            stat_name = k
        orig_cols = stat_df.columns.to_list()
        stat_df = stat_df.assign(**keys)[list(keys.keys()) + orig_cols]
        stats_out[stat_name] = stat_df

    return stats_out


def collect_samtools_stats(input_fps, output_p, keys, metadata_table):
    stats = defaultdict(list)
    for fp, (_, row_ser) in zip(
        input_fps, metadata_table[["entity", "sample"]].iterrows()
    ):
        curr_stats = parse_samtools_stats(fp, dict(**keys, **row_ser.to_dict()))
        for k, v in curr_stats.items():
            stats[k].append(v)
    for k, v in stats.items():
        stats[k] = pd.concat(v)
    to_pickle(stats, output_p)


def test_collect_samtools_stats():

    metadata_table = pd.DataFrame({"entity": "entA", "sample": ["a", "b"]})

    collect_samtools_stats(
        input_fps=[
            "/icgc/dkfzlsdf/analysis/hs_ontogeny/results/sctm/results-per-entity/lsk150_pr-hiera_ex-sctm-1_cls-1/blood1/smk-wgbs/default-params_v0.1.0/alignments/lsk150_pr-hiera_ex-sctm-1_cls-1_blood1_SE_PE_samtools-stats.txt",
            "/icgc/dkfzlsdf/analysis/hs_ontogeny/results/sctm/results-per-entity/lsk150_pr-hiera_ex-sctm-1_cls-1/blood2/smk-wgbs/default-params_v0.1.0/alignments/lsk150_pr-hiera_ex-sctm-1_cls-1_blood2_SE_PE_samtools-stats.txt",
        ],
        output_p="/home/kraemers/temp/test.p",
        keys=dict(e="exp1", a="al1", m="meth1"),
        metadata_table=metadata_table,
    )
