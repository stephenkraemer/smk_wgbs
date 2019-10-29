"""
export PATH=/home/kraemers/projects/Bismark:$PATH

cd /icgc/dkfzlsdf/analysis/hs_ontogeny

snakemake \
--snakefile $(smk_wgbs_snakefile) \
--latency-wait 60 \
--configfile $(smk_wgbs_demo_config) \
--cluster "bsub -R rusage[mem={params.avg_mem}] -M {params.max_mem} -n {threads} -J {params.name} -W {params.walltime} -o /home/kraemers/temp/logs/" \
--jobs 1000 \
--keep-going \
--forcerun trim_reads_pe \
--config \
  experiment_name=scnmt-1 \
  protocol=scnmt \
  'entities=[
    "lsk150_pr-hiera_ex-scnmt-1_cls-25"
  ]' \
  upto=10000 \
--dryrun \



# scm
--config \
  experiment_name=pilot-scM \
  protocol=sctm \
  'entities=[
     "hsc_pr-scbs_ex-test-1_cls-1_2",
  ]' \
  upto=10000 \

# scmt
--config \
  experiment_name=sctm-1 \
  protocol=sctm \
  'entities=[
     "lsk150_pr-hiera_ex-sctm-1_cls-0",
     "lsk150_pr-hiera_ex-sctm-1_cls-1",
     "lsk150_pr-hiera_ex-sctm-1_cls-25"
  ]' \

# scnmt
--config \
  experiment_name=scnmt-1 \
  protocol=scnmt \
  'entities=[
    "lsk150_pr-hiera_ex-scnmt-1_cls-1",
    "lsk150_pr-hiera_ex-scnmt-1_cls-0",
    "lsk150_pr-hiera_ex-scnmt-1_cls-1_methneg",
    "lsk150_pr-hiera_ex-scnmt-1_cls-25"
  ]' \


--forcerun multiqc \
--forcerun nome_filtering \
--rerun-incomplete \
-p \

"""

import time
import os
import re
from pathlib import Path
import pandas as pd
from smk_wgbs import create_metadata_table_from_file_pattern, sel_expand, find_workflow_version
from smk_wgbs.tools import fp_to_parquet, fp_to_bedgraph, run_methyldackel, methyldackel_bedgraph_to_bed_and_parquet
from functools import partial
import smk_wgbs.tools


# For development
# import yaml
# config_yaml = '/home/kraemers/projects/smk_wgbs/doc/demo_config.yaml'
# with open(config_yaml) as fin:
#     config = yaml.load(fin, Loader=yaml.FullLoader)

# Fill prefix with config_name and workflow_version
results_prefix = sel_expand(
        config['results_prefix'],
        config_name=config['name'],
        workflow_version=find_workflow_version()
)

# concatenate paths
config['log_dir'] = os.path.join(config['results_dir'], results_prefix, config['log_dir'])
for pattern_name, pattern in config['result_patterns'].items():
    config['result_patterns'][pattern_name] = (
        os.path.join(config['results_dir'], results_prefix, pattern))

# Prepare metadata table and adjust uid fields
# ======================================================================

# if a fastq_pattern is provided, it must have the fields: entity, sample, read_number and it must additionally have at least one additional field which serves to distinguish read pairs.
# The fastq_pattern is then used to create a metadata table with one column per field (using smk_wgbs.create_metadata_table_from_file_pattern).
# The metadata table is placed at metadata_table_tsv (with a timestamp before the suffix)
# without fastq_pattern: metadata table is read from metadata_table_tsv
if 'fastq_pattern' in config:
    if 'entities' in config:
        # create metadata table per entity, using only specified entities
        metadata_table = pd.concat([
            create_metadata_table_from_file_pattern(
                    sel_expand(config['fastq_pattern'], entity=entity)).assign(entity=entity)
            for entity in config['entities']], axis=0)
    else:
        raise NotImplementedError

    # fill in optional pattern fields
    if 'lib' not in metadata_table:
        metadata_table['lib'] = '1'
    if 'protocol' not in metadata_table:
        metadata_table['protocol'] = config['protocol']

    metadata_table = metadata_table.astype(str)

    timestamp = time.strftime('%d-%m-%y_%H:%M:%S')
    metadata_table_tsv = Path(config['metadata_table_tsv']).with_suffix(f'.{timestamp}.tsv')
    Path(metadata_table_tsv).parent.mkdir(parents=True, exist_ok=True)
    metadata_table.to_csv(metadata_table_tsv, header=True, index=False, sep='\t')

    # REMOVE
    # metadata_table = metadata_table.iloc[0:10].copy()
    # \REMOVE

else:
    metadata_table = pd.read_csv(config['metadata_table_tsv'], sep='\t', header=0, index=False)

uid_columns = metadata_table.columns[
    ~metadata_table.columns.isin(['entity', 'sample', 'read_number', 'path', 'lib', 'protocol'])]

# add read pair UID field to metadata table
metadata_table['uid'] = metadata_table[uid_columns].apply(
        lambda ser: '_'.join(f'{uid_col}-{ser[uid_col]}' for uid_col in uid_columns), axis=1)


# Targets
# ======================================================================
wildcard_constraints:
    lib = '[^_]+',
    motif = '(CG|CHG|CHH)',
    region = '.*',


is_nome = (config['protocol_config'][config['protocol']].get('is_nome_seq', False))
if is_nome:
    for motif, user_request in config['protocol_config'][config['protocol']]['motifs'].items():
        if not user_request:
            raise ValueError(f'NOMe-Seq requires calls for all motifs, but motif {motif} is turned off')
        elif user_request == 'sampled':
            print(f'WARNING: NOMe seq, but {motif} is set to "sampled"')

sampling_name = '1:40e-80e'

targets = []
for _, row_ser in metadata_table.iterrows():
    for motif, user_request in config['protocol_config'][config['protocol']]['motifs'].items():

        if not user_request:
            continue
        elif user_request == 'sampled':
            region = f'_sampled-{sampling_name}'
        elif user_request == 'full':
            region = ''
        else:
            raise ValueError()

        if config['protocol_config'][config['protocol']].get('is_nome_seq', False):
            patterns = [sel_expand(config['result_patterns']['cg_full_parquet'],
                                   alignment_combi=config['alignment_combi'])]
        else:
            patterns = [sel_expand(config['result_patterns']['final_mcalls_merged_per_motif'],
                                   alignment_combi=config['alignment_combi'])]
            if motif == 'CHH':
                pass
            elif motif in ['CG', 'CHG']:
                pairings = config['alignment_combi'].split('_')
                for pairing in pairings:
                    patterns.append(
                            sel_expand(
                                    config['result_patterns']['pe_or_se_mcalls_merged_per_motif'],
                                    pairing=pairing
                            )
                    )
            else:
                ValueError()
        for pattern in patterns:
            targets.append(sel_expand(
                    pattern,
                    **row_ser,
                    motif=motif,
                    region=region,
            ))
print(sorted(targets), sep='\n')
multiqc_report_fp = sel_expand(config['multiqc_by_experiment'],
                               protocol=config['protocol'],
                               experiment=config['experiment_name'])


fastqc_reports = (
    metadata_table[['protocol', 'entity', 'sample', 'lib', 'uid', 'read_number']]
        .apply(lambda ser: expand(config['result_patterns']['fastq_fastqc'], **ser)[0], axis=1)
)


rule all:
    input:
         targets,
         multiqc_report_fp,
         # fastqc_reports,

# SE Alignment
# ======================================================================

# on --parallel
# If system resources are plentiful this is a viable option to speed up the alignment process
# (we observed a near linear speed increase for up to --parallel 8 tested). However, please note
# that a typical Bismark run will use several cores already (Bismark itself, 2 or 4 threads of
# Bowtie2/HISAT2, Samtools, gzip etc...) and ~10-16GB of memory depending on the choice of aligner
# and genome. WARNING: Bismark Parallel (BP?) is resource hungry! Each value of --parallel specified
# will effectively lead to a linear increase in compute and memory requirements, so --parallel 4 for
#     e.g. the GRCm38 mouse genome will probably use ~20 cores and eat ~40GB or RAM, but at the same time
# reduce the alignment time to ~25-30%. You have been warned.

"""
-o/--output_dir <dir>    Write all output files into this directory. By default the output files will be written into
                         the same folder as the input file(s). If the specified folder does not exist, Bismark will attempt
                         to create it first. The path to the output folder can be either relative or absolute.

--temp_dir <dir>         Write temporary files to this directory instead of into the same directory as the input files. If
                         the specified folder does not exist, Bismark will attempt to create it first. The path to the
                         temporary folder can be either relative or absolute.
"""
if config['alignment_combi'] == 'SE_PE':
    se_fq_input = config['result_patterns']['unmapped_fastq']
else:
    se_fq_input = config['result_patterns']['trimmed_fastq']

rule bismark_se_atomic_alignment:
    input:
         fq = se_fq_input,
         # not directly used, but still import files (parent directory path is bismark arg)
         refconvert_CT = config['genome_dir'] + "/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
         refconvert_GA = config['genome_dir'] + "/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
    output:
          bam = config['result_patterns']['se_atomic_bam_unsorted'],
          bismark_report = config['result_patterns']['se_atomic_bam_unsorted'].replace('.bam', '_SE_report.txt'),
          temp_dir = temp(directory(config['result_patterns']['se_atomic_bam_unsorted'] + '___bismark_tempfiles'))
    # non specified output files: will not be auto-removed
    params:
        genome_dir = config['genome_dir'],
        avg_mem = 12000,
        max_mem = 16000,
        name = 'bismark_alignment_SE_{entity}_{sample}',
        walltime = '08:00',
        # walltime = '00:30',
    threads: 6
    log:
        config['log_dir'] + '/bismark_atomic-alignment_SE_{protocol}_{entity}_{sample}_{uid}_{lib}_R{read_number}.log'
    # basename and multicore together are not possible
    # --score_min G,10,2 \
    # --local
    # --parallel 4 \
    run:
        smk_wgbs.tools.run_bismark_alignment(
                input, output, log, params, is_paired=False, upto=config['upto'],
        )


rule bismark_pe_atomic_alignment:
    input:
         fq_r1 = sel_expand(config['result_patterns']['trimmed_fastq'], read_number='1'),
         fq_r2 = sel_expand(config['result_patterns']['trimmed_fastq'], read_number='2'),
         # not directly used, but still import files (parent directory path is bismark arg)
         refconvert_CT = config['genome_dir'] + "/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
         refconvert_GA = config['genome_dir'] + "/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
    output:
          bam = config['result_patterns']['pe_atomic_bam_unsorted'],
          unmapped_fq1 = sel_expand(config['result_patterns']['unmapped_fastq'],
                                    read_number='1'),
          unmapped_fq2 = sel_expand(config['result_patterns']['unmapped_fastq'],
                                    read_number='2'),
          bismark_report = config['result_patterns']['pe_atomic_bam_unsorted'].replace('.bam', '_PE_report.txt'),
          temp_dir = temp(directory(config['result_patterns']['pe_atomic_bam_unsorted'] + '___bismark_tempfiles'))
    # non specified output files: will not be auto-removed
    params:
          genome_dir = config['genome_dir'],
          avg_mem = 24000,
          max_mem = 32000,
          name = 'bismark_alignment_PE_{entity}_{sample}',
          walltime = '08:00',
          # walltime = '00:30',
    threads: 6
    log:
       config['log_dir'] + '/bismark_atomic-alignment_PE_{protocol}_{entity}_{sample}_{uid}_{lib}.log'
    # basename and multicore together are not possible
    # --score_min G,10,2 \
    # --local
    # --parallel 4 \
    run:
        smk_wgbs.tools.run_bismark_alignment(
                input, output, log, params, is_paired=True, upto=config['upto'],
        )

def find_final_alignment_input(wildcards):
    if wildcards.alignment_combi == 'SE':
        return expand(config['result_patterns']['pe_or_se_bam'], **wildcards, pairing='SE')
    elif wildcards.alignment_combi == 'PE':
        return expand(config['result_patterns']['pe_or_se_bam'], **wildcards, pairing='PE')
    else:
        return expand(config['result_patterns']['pe_or_se_bam'], **wildcards, pairing=['PE', 'SE'])

rule final_alignment:
    input:
         find_final_alignment_input,
    output:
          bam = config['result_patterns']['final_bam'],
    params:
          avg_mem = 6000,
          max_mem = 8000,
          name = 'merge-final-alignment_{entity}_{sample}',
          walltime = '02:00',
    run:
        smk_wgbs.tools.merge_or_symlink(input, output)


rule bismark_se_sort_atomic_bam:
    input:
         bam = config['result_patterns']['se_atomic_bam_unsorted'],
    output:
          bam = config['result_patterns']['se_atomic_bam_sorted'],
    params:
          # default -m: 768 Mb
          avg_mem = 8 * 1000,
          max_mem = 8 * 1600,
          walltime = '00:30',
          name = 'sort_atomic_bam_{entity}_{sample}_{lib}_{uid}'
    threads: 8
    shell:
        "samtools sort -o {output} -T {output} -@ 8 {input}"


rule bismark_pe_sort_atomic_bam:
    input:
         bam = config['result_patterns']['pe_atomic_bam_unsorted'],
    output:
          bam = config['result_patterns']['pe_atomic_bam_sorted'],
    params:
          # default -m: 768 Mb
          avg_mem = 8 * 1000,
          max_mem = 8 * 1600,
          walltime = '00:30',
          name = 'sort_atomic_bam_PE_{entity}_{sample}_{lib}_{uid}'
    threads: 8
    shell:
         "samtools sort -o {output} -T {output} -@ 8 {input}"


def find_atomic_bams_for_library(wildcards):
    # import pudb; pudb.set_trace()
    df = metadata_table
    matches_library_and_read = (
            (df['protocol'] == wildcards.protocol)
            & (df['entity'] == wildcards.entity)
            & (df['sample'] == wildcards.sample)
            & (df['lib'] == wildcards.lib)
            # & (df['read_number'] == wildcards.read_number)
    )
    wildcard_cols = ['protocol', 'entity', 'sample', 'lib', 'uid']
    if wildcards.pairing == 'SE':
        atomic_bam_pattern = config['result_patterns']['se_atomic_bam_sorted']
        wildcard_cols += ['read_number']
    else:
        atomic_bam_pattern = config['result_patterns']['pe_atomic_bam_sorted']
    # duplicates only allowed for PE
    wildcard_df = df.loc[matches_library_and_read, wildcard_cols]
    if wildcards.pairing == 'SE':
        assert wildcard_df.duplicated().sum() == 0
    else:
        wildcard_df = wildcard_df.drop_duplicates()
    ser = wildcard_df.apply(lambda ser: expand(atomic_bam_pattern, **ser)[0], axis=1)
    l = ser.to_list()
    print(l)
    return l


rule bismark_merge_library_mdup:
    input:
        find_atomic_bams_for_library
    output:
         bam = config['result_patterns']['pe_or_se_library_bam'],
         bai = config['result_patterns']['pe_or_se_library_bam'] + '.bai',
         metrics = re.sub(
                 '.bam$',
                 '_mdup-metrics.txt',
                 config['result_patterns']['pe_or_se_library_bam']),
    params:
          # default -m: 768 Mb
          avg_mem = 8 * 1000,
          max_mem = 8 * 1600,
          max_mem_gb = '13G',
          walltime = '00:30',
          name = 'mdup_merge_{pairing}_{entity}_{sample}_{lib}',
          input_spec = lambda wildcards, input: ' '.join(f'INPUT={i}' for i in input)
    # java  -jar picard.jar
    # OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \ #changed from default of 100
    shell:
        """
        picard -Xmx{params.max_mem_gb} MarkDuplicates \
        {params.input_spec} \
        OUTPUT={output.bam} \
        METRICS_FILE={output.metrics} \
        CREATE_INDEX=false \
        TMP_DIR=/tmp
        
        samtools index {output.bam}
        """


def find_library_bams(wildcards):
    df = metadata_table
    matches_library = (
            (df['protocol'] == wildcards.protocol)
            & (df['entity'] == wildcards.entity)
            & (df['sample'] == wildcards.sample)
    )
    ser = (df
           .loc[matches_library, ['protocol', 'entity', 'sample', 'lib']]
           .drop_duplicates()
            # expand always returns list, get string out of list of length 1
           .apply(lambda ser: expand(
                    config['result_patterns']['pe_or_se_library_bam'],
                    **ser,
                    pairing=wildcards.pairing)[0],
                  axis=1)
           )
    l = ser.to_list()
    return l


rule bismark_merge_libraries:
    input:
        find_library_bams
    output:
        bam = config['result_patterns']['pe_or_se_bam'],
        bai = config['result_patterns']['pe_or_se_bam'] + '.bai',
    # TODO: function which returns **resources and returns minimal resources if only hard linking is required
    params:
          avg_mem = 6000,
          max_mem = 8000,
          name = 'merge-lib_{entity}_{sample}_{pairing}',
          walltime = '00:30',
    run:
        smk_wgbs.tools.merge_or_symlink(input, output)



# ==============================================================================

def get_motif_params_str(wildcards):
    chromosome1_name = config['chromosomes'][0]
    if wildcards.region:
        sample_command = f'-r {chromosome1_name}:40000000-100000000'
    else:
        sample_command = ''
    if wildcards.motif == 'CG':
        return ''
    elif wildcards.motif == 'CHH':
        return f'--CHH --noCpG {sample_command}'
    elif wildcards.motif == 'CHH':
        return f'--CHG --noCpG {sample_command}'
    else:
        raise ValueError()

def get_methyldackel_walltime(wildcards):
    if wildcards.region:
        sampled_str = 'sampled'
    else:
        sampled_str = 'full'
    walltime_d = {
        'CHH': {
            'sampled': '02:00',
            'full': '10:00',
        },
        'CHG': {
            'sampled': '02:00',
            'full': '10:00',
        },
        'CG': {
            'sampled': '01:00',
            'full': '02:00',
        },
    }
    return walltime_d[wildcards.motif][sampled_str]


methyldackel_input = dict(
    bam = config['result_patterns']['pe_or_se_bam'],
    ref_genome_unconverted = ancient(config['genome_fa']),
)

# rule for full runs
def get_methyldackel_output_params(
        cg_walltime,
        chh_walltime,
        sampling_name,
        calling_mode
):
    methyldackel_output = []
    motif_params_l = ['--noCpG']
    methyldackel_walltime = cg_walltime
    for motif, motif_calling_mode in config['protocol_config'][config['protocol']]['motifs'].items():
        if motif_calling_mode == calling_mode:
            if motif == 'CG':
                motif_params_l.remove('--noCpG')
            else:  # CHH or CHG
                assert motif in ['CHH', 'CHG']
                motif_params_l.append(f'--{motif}')
                methyldackel_walltime = chh_walltime  # walltime for full run including CHH and/or CHG
            bed = Path(sel_expand(
                    config['result_patterns']['pe_or_se_mcalls_cytosines_per_motif'],
                    motif=motif,
                    region=f'_sampled-{sampling_name}' if sampling_name else '',
            ))
            methyldackel_motif = motif if motif != 'CG' else 'CpG'
            methyldackel_output += [
                bed,
                bed.with_suffix('.bedGraph'),
                bed.with_suffix('.parquet'),
                temp(bed.parent.joinpath('{pairing}_{entity}_{sample}' + f'_{methyldackel_motif}.bedGraph')),
            ]
    motif_params_str = ' '.join(motif_params_l)
    if not methyldackel_output:
        methyldackel_output = [sel_expand(config['result_patterns']['pe_or_se_mcalls_cytosines_per_motif'],
                               motif='UNUSED',
                               region='UNUSED')]
    return methyldackel_output, methyldackel_walltime, motif_params_str

methyldackel_output, methyldackel_walltime, motif_params_str = get_methyldackel_output_params(
        cg_walltime='02:00',
        chh_walltime='10:00',
        sampling_name='',
        calling_mode='full'
)
rule methyldackel_per_cytosine_per_motif:
    input:
          **methyldackel_input,
    output:
          methyldackel_output,
    params:
          avg_mem = 6000,
          max_mem = 10000,
          name = 'methyldackel_per_cytosine_{pairing}_{entity}_{sample}',
          walltime = methyldackel_walltime,
          prefix = lambda wildcards, output: Path(output[0]).parent.joinpath(
                  f'{wildcards.pairing}_{wildcards.entity}_{wildcards.sample}'),
          # TODO: is config directly available in rule?
          chromosomes = config['chromosomes'],
          motif_params = motif_params_str,
          sampling_options = '',
    # TODO: defaults for maxVariantFrac and minOppositeDepth
    threads: 4
    log:
       config['log_dir'] + '/methyldackel_per_cytosine_per_motif:{pairing}_{protocol}_{entity}_{sample}.log'
    run:
        run_methyldackel(input, output, params, threads)

methyldackel_output, methyldackel_walltime, motif_params_str = get_methyldackel_output_params(
        cg_walltime='00:30',
        chh_walltime='02:00',
        sampling_name=sampling_name,
        calling_mode='sampled')

rule methyldackel_per_cytosine_per_motif_sampled:
    input:
         **methyldackel_input,
    output:
          methyldackel_output,
    params:
          avg_mem = 6000,
          max_mem = 10000,
          name = 'methyldackel_se_per_cytosine_{entity}_{sample}' + sampling_name,
          walltime = methyldackel_walltime,
          prefix = lambda wildcards, output: Path(output[0]).parent.joinpath(
                  f'{wildcards.pairing}_{wildcards.entity}_{wildcards.sample}'),
          # TODO: is config directly available in rule?
          chromosomes = config['chromosomes'],
          motif_params = motif_params_str,
          sampling_options = '-r 1:40000000-100000000',
    # TODO: defaults for maxVariantFrac and minOppositeDepth
    threads: 4
    log:
       config['log_dir'] + '/methyldackel_per_cytosine_per_motif_sampled:{pairing}_{protocol}_{entity}_{sample}.log'
    run:
        run_methyldackel(input, output, params, threads)


rule methyldackel_merge_context:
    input:
        ref_genome_unconverted = ancient(config['genome_fa']),
        bedgraph = fp_to_bedgraph(config['result_patterns']['pe_or_se_mcalls_cytosines_per_motif']),
    output:
        bed = config['result_patterns']['pe_or_se_mcalls_merged_per_motif'],
        bedgraph = fp_to_bedgraph(config['result_patterns']['pe_or_se_mcalls_merged_per_motif']),
        parquet = fp_to_parquet(config['result_patterns']['pe_or_se_mcalls_merged_per_motif']),
    params:
        chromosomes = config['chromosomes'],
        avg_mem = 8000,
        max_mem = 12000,
        walltime = '00:20',
        name = 'methyldackel_SE_mergeContext_{entity}_{sample}_{motif}{region}',
    log:
        config['log_dir'] + '/methyldackel_mergeContext_{pairing}_{entity}_{sample}_{motif}{region}'
    run:
        from pandas.api.types import CategoricalDtype

        print('run MethylDackel MergeContext')
        shell("""
        MethylDackel mergeContext \
        -o {output.bedgraph} \
        {input.ref_genome_unconverted} \
        {input.bedgraph}
        """)

        print('Convert BedGraph to other file formats')
        methyldackel_bedgraph_to_bed_and_parquet(
                bed=output.bed,
                bedgraph=output.bedgraph,
                chrom_dtype = CategoricalDtype(categories=params.chromosomes, ordered=True),
                parquet=output.parquet
        )

def find_merge_cytosine_level_calls_input(wildcards):
    if wildcards.alignment_combi in ['SE', 'PE']:
        return expand(fp_to_parquet(config['result_patterns']['pe_or_se_mcalls_cytosines_per_motif']),
                      **wildcards,
                      pairing=wildcards.alignment_combi)
    else:
        return expand(fp_to_parquet(config['result_patterns']['pe_or_se_mcalls_cytosines_per_motif']),
                      **wildcards,
                      pairing=['PE', 'SE'])

final_mcalls_cyt_bed = config['result_patterns']['final_mcalls_cytosines_per_motif']
rule final_mcalls:
    input:
         find_merge_cytosine_level_calls_input
    output:
         bed = final_mcalls_cyt_bed,
         parquet = fp_to_parquet(final_mcalls_cyt_bed),
         bedgraph = fp_to_bedgraph(final_mcalls_cyt_bed),
    params:
          avg_mem = 8000,
          max_mem = 12000,
          walltime = '00:30',
          name = 'merge-final-calls_{entity}_{sample}'
    run:
        smk_wgbs.tools.merge_mcalls(input, output, wildcards)

final_mcalls_merged_bed = config['result_patterns']['final_mcalls_merged_per_motif']
rule final_mcalls_merged_motifs:
    input:
         bedgraph = fp_to_bedgraph(final_mcalls_cyt_bed),
         ref_genome_unconverted = ancient(config['genome_fa']),
    output:
         bed = final_mcalls_merged_bed,
         bedgraph = fp_to_bedgraph(final_mcalls_merged_bed),
         parquet = fp_to_parquet(final_mcalls_merged_bed)
    params:
          avg_mem = 2000,
          max_mem = 4000,
          walltime = '01:00',
          name = 'merge-context_final-mcalls_{entity}_{sample}_{motif}',
          chromosomes = config['chromosomes'],
    run:
        from pandas.api.types import CategoricalDtype

        print('run MethylDackel MergeContext')
        shell("""
        MethylDackel mergeContext \
        -o {output.bedgraph} \
        {input.ref_genome_unconverted} \
        {input.bedgraph}
        """)

        print('Convert BedGraph to other file formats')
        methyldackel_bedgraph_to_bed_and_parquet(
                bed=output.bed,
                bedgraph=output.bedgraph,
                chrom_dtype = CategoricalDtype(categories=params.chromosomes, ordered=True),
                parquet=output.parquet
        )

# ==========================================================================================
# Generate methyl-converted version of the reference genome, if necessary:
# indexing is run twice in parallel, so this will use twice the specified amount of threads
# code (with help at bottom: https://github.com/FelixKrueger/Bismark/blob/master/bismark_genome_preparation)
# mouse genome run took a little more than 2h
# ~ 10GB max memory for mouse genome with --parallel 32
rule bismark_genome_preparation:
    input:
         ancient(config['genome_dir'])
    output:
          config['genome_dir'] + "/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
          config['genome_dir'] + "/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"
    params:
        avg_mem = 20000,
        max_mem = 24000,
        name = 'ref_conversion',
        walltime = '08:00',
    threads: 64
    # message: "Converting Genome into Bisulfite analogue"
    shell:
         """
         bismark_genome_preparation \
         --bowtie2 \
         --parallel 32 \
         --genomic_composition \
         --verbose \
         {input} \
         > {log} 2>&1
         """


# ==========================================================================================
# Trim the reads for adapter-ends and quality

# cutadapt speed up linear with cores
# parallelization requires pigz python package to be installed (by default if bismark is installed via conda)

rule trim_reads_pe:
    input:
        unpack(partial(smk_wgbs.tools.find_fastqs, metadata_table=metadata_table)),
    output:
          fq_r1 = sel_expand(config['result_patterns']['trimmed_fastq'], read_number='1'),
          fq_r2 = sel_expand(config['result_patterns']['trimmed_fastq'], read_number='2'),
    log:
       config['log_dir'] + '/cutadapt_{protocol}_{entity}_{sample}_{uid}_{lib}.log'
    # message:
    #        "Trimming raw paired-end read data\n{input.fq_r1}\n{input.fq_r2}"
    params:
          avg_mem = 4000,
          max_mem = 6000,
          name = 'trim_reads_{entity}_{sample}',
          walltime = '00:30',
          phred_filter_command = '--nextseq-trim 10' if config['illumina_machine'] == 'nextseq' else '-q 10',
    threads: 16
    # currently, do not trim random hexamers priming sites \
    # note that this anyway does not take care of priming sites at the end of short reads \
    # -u 6
    # -U 6
    # discard read pair if any read is bad; given that we do SE alignment this may not be ideal
    # --pair-filter=any \
    shell:
         """
         cutadapt \
         --cores {threads} \
         -a AGATCGGAAGAGCG \
         -A AGATCGGAAGAGCG \
         -o {output.fq_r1} \
         -p {output.fq_r2} \
         {params.phred_filter_command} \
         -u 6 \
         -U 6 \
         --minimum-length 30 \
         --max-n 0.3 \
         --pair-filter=any \
         {input.fq_r1} \
         {input.fq_r2} \
         > {log} 2>&1
         """


input_d = {}
for motif, calling_mode in config['protocol_config'][config['protocol']]['motifs'].items():
    if calling_mode == 'full':
        region = ''
    else:
        region = sampling_name
    input_d[motif] = fp_to_parquet(
            sel_expand(
                    config['result_patterns']['final_mcalls_cytosines_per_motif'],
                    motif=motif,
                    region=region,
                    alignment_combi=config['alignment_combi']
            )
    )


RESPAT = config['result_patterns']
rule nome_filtering:
    input:
         **input_d
    output:
          cg_full_parquet = RESPAT['cg_full_parquet'],
          cg_meth_parquet = RESPAT['cg_meth_parquet'],
          ch_full_parquet = RESPAT['ch_full_parquet'],
          ch_acc_parquet = RESPAT['ch_acc_parquet'],
          chg_meth_parquet = RESPAT['chg_meth_parquet'],
          chh_meth_parquet = RESPAT['chh_meth_parquet'],
    params:
          avg_mem = 60000,
          max_mem = 80000,
          walltime = '04:00',
          name = 'nome_filtering_{entity}_{sample}',
          chromosomes = config['chromosomes'],
    log:
       config['log_dir'] + '/nome-filtering_{protocol}_{alignment_combi}_{entity}_{sample}{region}.log'
    run:
        smk_wgbs.tools.nome_filtering(input, output, params)


rule samtools_stats:
    input:
         bam = config['result_patterns']['final_bam'],
         ref_genome_unconverted = ancient(config['genome_fa']),
    output:
          config['result_patterns']['final_bam_samtools_stats'],
    params:
          avg_mem = 500,
          max_mem = 1000,
          walltime = '00:20',
          name = 'samtools_stats_{entity}_{sample}',
    shell:
         """
         samtools stats \
         --insert-size 3000 \
         --ref-seq {input.ref_genome_unconverted} \
         --remove-dups \
         {input.bam} > {output[0]}
         """


rule samtools_flagstats:
    input:
         bam = config['result_patterns']['final_bam'],
    output:
          config['result_patterns']['final_bam_samtools_flagstats'],
    params:
          avg_mem = 500,
          max_mem = 1000,
          walltime = '00:20',
          name = 'samtools_flagstats_{entity}_{sample}',
    shell:
         """
         samtools flagstat {input.bam} > {output[0]}
         """


samtools_stat_files = (
    metadata_table[['entity', 'sample', 'protocol']]
        .drop_duplicates()
        .apply(
            lambda ser: sel_expand(
                    config['result_patterns']['final_bam_samtools_stats'],
                    **ser,
                    alignment_combi=config['alignment_combi'],
            ),
            axis=1
    )
)

samtools_flagstat_files = (
    metadata_table[['entity', 'sample', 'protocol']]
        .drop_duplicates()
        .apply(
            lambda ser: sel_expand(
                    config['result_patterns']['final_bam_samtools_flagstats'],
                    **ser,
                    alignment_combi=config['alignment_combi'],
            ),
            axis=1
    )
)


bismark_report_pattern = (config['result_patterns']['se_atomic_bam_unsorted']
                          .replace('.bam', '_SE_report.txt'))
bismark_se_reports = (metadata_table[['protocol', 'entity', 'sample', 'lib', 'uid', 'read_number']]
    .apply(lambda ser: expand(bismark_report_pattern, **ser)[0], axis=1)
)

rule multiqc:
    input:
         samtools_stat_files,
         # bismark_se_reports,
         # fastqc_reports,
         samtools_flagstat_files,
    output:
          report_html=multiqc_report_fp,
          file_list=multiqc_report_fp + '_file-list.txt',
    params:
          avg_mem = 5000,
          max_mem = 10000,
          walltime = '00:15',
          name = f'multiqc_{config["name"]}',
    run:
        with open(output.file_list, 'wt') as fout:
            fout.write('\n'.join(input) + '\n')
        shell(f"""
              multiqc \
              --filename {output.report_html} \
              --force \
              --file-list {output.file_list}
              """
        )


rule fastqc:
    input:
         unpack(partial(smk_wgbs.tools.find_fastqs, metadata_table=metadata_table)),
    output:
          qc_r1 = sel_expand(config['result_patterns']['fastq_fastqc'], read_number='1'),
          fq_r1 = temp(sel_expand(config['result_patterns']['fastq_fastqc'], read_number='1').replace('_fastqc.zip', '.fastq.gz')),
          qc_r2 = sel_expand(config['result_patterns']['fastq_fastqc'], read_number='2'),
          fq_r2 = temp(sel_expand(config['result_patterns']['fastq_fastqc'], read_number='2').replace('_fastqc.zip', '.fastq.gz')),
    params:
          output_dir = lambda wildcards, output: Path(output[0]).parent,
          avg_mem = 2000,
          max_mem = 4000,
          walltime = '00:30',
          name = 'fastqc',
    threads: 2
    # -c {input.adapter_sequences_tsv} \
    run:
        import os
        from pathlib import Path
        os.symlink(input.fq_r1, output.fq_r1)
        os.symlink(input.fq_r2, output.fq_r2)
        # -o {params.output_dir} \
        shell(
                """
                fastqc \
                --noextract \
                --threads {threads} \
                {output.fq_r1} {output.fq_r2}
                """
         )
         # fastqc_r1_fp = Path(params.output_dir).joinpath(Path(input.fq_r1).name.replace('.fastq.gz', '_fastqc.zip'))
         # fastqc_r2_fp = Path(params.output_dir).joinpath(Path(input.fq_r2).name.replace('.fastq.gz', '_fastqc.zip'))
         # fastqc_r1_fp.rename(output.qc_r1)
         # fastqc_r2_fp.rename(output.qc_r2)


