import glob
import itertools
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Dict

import pandas as pd
from pandas.api.types import CategoricalDtype
from pkg_resources import resource_filename


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
        cg_chrom_df_annotated = pd.merge(
            cg_chrom_df, chrom_index.loc[chrom_index["motif"] == "CG"], how="inner"
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
        ch_chrom_df_annotated = pd.merge(
            ch_chrom_df, chrom_index.loc[chrom_index["motif"] != "CG"], how="inner"
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
            --ignoreFlags 3840 \
            --requireFlags 0 \
            {params.motif_params} \
            {params.sampling_options} \
            -o {params.prefix} \
            -q 15 \
            -p 15 \
            -@ {threads} \
            {input.ref_genome_unconverted} \
            {input.bam}
            """.split()
    )
    print("Done with MethylDackel")

    chrom_dtype = CategoricalDtype(categories=params.chromosomes, ordered=True)

    assert len(output) % 4 == 0
    for i in range(int(len(output) / 4)):
        bed, bedgraph, parquet, tmp_bedgraph = output[4 * i : (4 * i) + 4]
        # convert bedgraph to bed and parquet
        shutil.copy(tmp_bedgraph, bedgraph)
        methyldackel_bedgraph_to_bed_and_parquet(bed, bedgraph, chrom_dtype, parquet)


def methyldackel_bedgraph_to_bed_and_parquet(bed, bedgraph, chrom_dtype, parquet):
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
    if "CpG" in bed:
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
        upto_option = f'--upto {upto}'
        print(f'WARNING: only performing the first {upto} alignments')
    else:
        upto_option = ''
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

        bismark_unmapped_r1 = output_dir.joinpath(Path(input.fq_r1).name + "_unmapped_reads_1.fq.gz")
        bismark_unmapped_r1.rename(output.unmapped_fq1)

        bismark_unmapped_r2 = output_dir.joinpath(Path(input.fq_r2).name + "_unmapped_reads_2.fq.gz")
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


def merge_mcalls(input, output, wildcards):
    if len(input) == 1:
        os.link(input.parquet, output.parquet)
        os.link(input.bed, output.bed)
        os.link(input.bedgraph, output.bedgraph)
        return

    else:
        assert len(input) == 2
        df1 = pd.read_parquet(input[0])
        assert df1["motif"].nunique() == 1
        df2 = pd.read_parquet(input[1])
        assert df2["motif"].nunique() == 1
        res = pd.merge(df1, df2, on=["Chromosome", "Start", "End"])
        res["n_meth"] = res["n_meth_x"] + res["n_meth_y"]
        res["n_total"] = res["n_total_x"] + res["n_total_y"]
        res["beta_value"] = res["n_meth"] / res["n_total"]
        res["motif"] = res["motif_x"]
        res = res[
            ["Chromosome", "Start", "End", "motif", "beta_value", "n_total", "n_meth"]
        ]
        res.to_parquet(output.parquet)
        res.rename(columns={"Chromosome": "#Chromosome"}).to_csv(
            output.bed, sep="\t", header=True, index=False
        )
        res['n_unmeth'] = res['n_total'] - res['n_meth']
        with open(output.bedgraph, 'wt') as fout:
            fout.write(f'track type="bedGraph" description="{wildcards.entity}_{wildcards.sample}"')
        res[["Chromosome", "Start", "End", "beta_value", "n_meth", "n_unmeth"]].to_csv(
            output.bedgraph, mode='a', sep="\t", header=False, index=False
        )
