"""
snakemake \
--snakefile $(smk_wgbs_snakefile) \
--latency-wait 60 \
--configfile $(smk_wgbs_demo_config) \
--cluster "bsub -R rusage[mem={params.avg_mem}] -M {params.mem_mb} -n {threads} -J {params.name} -W {params.walltime} -o /home/kraemers/temp/logs/" \
--jobs 1000 \
--keep-going \
--config \
  experiment_name=sctm-1 \
  protocol=sctm \
  'entities=[
     "lsk150_pr-hiera_ex-sctm-1_cls-1",
  ]' \
  'samples=["blood2"]' \
  upto=100000 \
--forcerun trim_reads_pe \

--forcerun methyldackel_per_cytosine_per_motif \
--dryrun \






--forcerun trim_reads_pe \

# scm
--config \
  experiment_name=pilot-scM \
  protocol=sctm \
  'entities=[
     "hsc_pr-scbs_ex-test-1_cls-1_2",
  ]' \
  upto=10000 \

# sctm
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

import os
import re
from pathlib import Path
import pandas as pd
import numpy as np
from smk_wgbs import sel_expand, find_workflow_version
from smk_wgbs.tools import fp_to_parquet, fp_to_bedgraph, run_methyldackel, methyldackel_bedgraph_to_bed_and_parquet, fp_to_pickle
from functools import partial
import smk_wgbs.tools


# For development
# import yaml
# config_yaml = '/home/kraemers/projects/smk_wgbs/smk_wgbs/demo_config.yaml'
# with open(config_yaml) as fin:
#     config = yaml.load(fin, Loader=yaml.FullLoader)

# Fill prefix with alignment config_name and workflow_version
results_prefix = sel_expand(
        config['results_prefix'],
        config_name=config['alignment_config_name'],
        workflow_version=find_workflow_version()
)

# concatenate paths
# results per pid
config['log_dir'] = os.path.join(config['results_dir'], results_prefix, config['log_dir'])
for pattern_name, pattern in config['result_patterns'].items():
    pattern = pattern.replace('{mcalling_config_name}', config['mcalling_config_name'])
    config['result_patterns'][pattern_name] = (
        os.path.join(config['results_dir'], results_prefix, pattern))

# results per experiment
for pattern_name, pattern in config['experiment_result_patterns'].items():
    config['experiment_result_patterns'][pattern_name] = (
        os.path.join(
                config['results_dir'],
                expand(
                        pattern,
                        protocol=config['protocol'],
                        experiment=config['experiment_name'],
                        alignment_config_name=config['alignment_config_name'],
                        mcalling_config_name=config['mcalling_config_name'],
                )[0]
        )
    )

# databases
for pattern_name, pattern in config['database_patterns'].items():
    config['database_patterns'][pattern_name] = (
        os.path.join(config['database_dir'], pattern))



# Prepare metadata table and adjust uid fields
# ======================================================================

# if a fastq_pattern is provided, it must have the fields: entity, sample, read_number and it must additionally have at least one additional field which serves to distinguish read pairs.
# The fastq_pattern is then used to create a metadata table with one column per field (using smk_wgbs.create_metadata_table_from_file_pattern).
# The metadata table is placed at metadata_table_tsv (with a timestamp before the suffix)
# without fastq_pattern: metadata table is read from metadata_table_tsv
if config.get('fastq_pattern') and config.get('metadata_table_tsv'):
    raise ValueError()

if 'fastq_pattern' in config:
    metadata_table = smk_wgbs.tools.metadata_table_from_fastq_pattern(config)
else:  # metadata_table_tsv given
    metadata_table = pd.read_csv(
            config['metadata_table_tsv'], sep='\t', header=0, index=False
    )

uid_columns = metadata_table.columns[
    ~metadata_table.columns.isin(
            ['entity', 'sample', 'read_number', 'path', 'lib', 'protocol']
    )
]

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
            if motif == 'CHH':
                patterns = [sel_expand(config['result_patterns']['final_mcalls_cytosines_per_motif'],
                                   alignment_combi=config['alignment_combi'])]
            elif motif in ['CG', 'CHG']:
                patterns = [sel_expand(config['result_patterns']['final_mcalls_merged_per_motif'],
                                   alignment_combi=config['alignment_combi'])]
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
# print(sorted(targets), sep='\n')
multiqc_report_fp = sel_expand(config['multiqc_by_experiment'],
                               protocol=config['protocol'],
                               experiment=config['experiment_name'])

fastqc_reports = (
    metadata_table[['protocol', 'entity', 'sample', 'lib', 'uid', 'read_number']]
        .apply(lambda ser: expand(config['result_patterns']['fastq_fastqc'], **ser)[0], axis=1)
)


fastq_screen_dirs = (
    metadata_table[['protocol', 'entity', 'sample', 'lib', 'uid', 'read_number']]
        .apply(lambda ser: expand(config['result_patterns']['fastq_screen'], **ser)[0], axis=1)
)

# untrimmed
def get_fastq_screen_txt_fp(row_ser):
    wildcards = row_ser[['protocol', 'entity', 'sample', 'lib', 'uid', 'read_number']]
    out_dir = expand(config['result_patterns']['fastq_screen'], **wildcards)[0]
    fastq_stem = os.path.basename(row_ser['path']).replace('.fastq.gz', '')
    return out_dir + f'/{fastq_stem}_screen.txt'
fastq_screen_txts_untrimmed = metadata_table.apply(get_fastq_screen_txt_fp, axis=1)

# initial solution to allow turning some specialized functionality on and off
# implementation will change in the future
# currently, this is only targeted at fastq_screen from the experiment results
experiment_result_targets_d = config['experiment_result_patterns'].copy()
if not config['modules']['fastq_screen']:
    experiment_result_targets_d.pop('multiqc_read_level')
experiment_result_targets = list(experiment_result_targets_d.values())

rule all:
    input:
         targets,
         multiqc_report_fp,
         experiment_result_targets,
         # config['database_patterns']['fastq_screen_index_dir'] + '/fastq_screen.conf',
         # fastq_screen_dirs,
         # config['experiment_result_patterns']['multiqc_read_level'],
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


def bismark_se_alignment_get_walltime(wildcards, input, attempt):
    lane_size_gb = os.path.getsize(input.fq) / 1e9 * 2
    # 30 min minimum processing time, plus 2h per 4 GB of data
    walltime_min = min(360 + 60 * lane_size_gb, 590)
    # convert to string HH:MM
    # return f'{np.floor(walltime_min / 60):02.0f}:{walltime_min % 60:02.0f}'
    return int(walltime_min * attempt)


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
    resources:
        avg_mem = lambda wildcards, attempt: 18 * 1024,
        mem_mb = lambda wildcards, attempt: 18 * 1024,
        walltime_min = bismark_se_alignment_get_walltime,
        attempt = lambda wildcards, attempt: attempt,
    params:
        genome_dir = config['genome_dir'],
        name = 'bismark_alignment_SE_{entity}_{sample}',
        # walltime_min = '00:30',
    threads: 4
    log:
        config['log_dir'] + '/bismark_atomic-alignment_SE_{protocol}_{entity}_{sample}_{uid}_{lib}_R{read_number}.log'
    # basename and multicore together are not possible
    # --score_min G,10,2 \
    # --local
    # --parallel 4 \
    run:
        smk_wgbs.tools.abort_on_nth_attempt(resources.attempt, 4)
        smk_wgbs.tools.run_bismark_alignment(
                input, output, log, params, is_paired=False, upto=config['upto'],
        )


def bismark_pe_alignment_get_walltime(wildcards, input, attempt):
    lane_size_gb = os.path.getsize(input.fq_r1) / 1e9 * 2
    # 30 min minimum processing time, plus 2h per 4 GB of data
    walltime_min = min(420 + 60 * lane_size_gb, 590)
    # convert to string HH:MM
    # return f'{np.floor(walltime_min / 60):02.0f}:{walltime_min % 60:02.0f}'
    return int(walltime_min * attempt)


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
    resources:
         avg_mem = lambda wildcards, attempt: 18 * 1024,
         mem_mb = lambda wildcards, attempt: 18 * 1024,
         walltime_min = bismark_pe_alignment_get_walltime,
         # walltime_min = '00:30',
         attempt = lambda wildcards, attempt: attempt,
    params:
          genome_dir = config['genome_dir'],
          name = 'bismark_alignment_PE_{entity}_{sample}',
    threads: 4
    log:
       config['log_dir'] + '/bismark_atomic-alignment_PE_{protocol}_{entity}_{sample}_{uid}_{lib}.log'
    # basename and multicore together are not possible
    # --score_min G,10,2 \
    # --local
    # --parallel 4 \
    run:
        smk_wgbs.tools.abort_on_nth_attempt(resources.attempt, 4)
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
    resources:
         walltime_min = lambda wildcards, input, attempt: (20 if len(input) > 1 else 3) * attempt,
         avg_mem = lambda wildcards, input, attempt: (600 if len(input) > 1 else 100) * attempt,
         mem_mb = lambda wildcards, input, attempt: (800 if len(input) > 1 else 200) * attempt,
         attempt = lambda wildcards, attempt: attempt,
    params:
         name = 'merge-final-alignment_{entity}_{sample}',
    run:
        smk_wgbs.tools.abort_on_nth_attempt(resources.attempt, 4)
        smk_wgbs.tools.merge_or_symlink(input, output)


rule bismark_se_sort_atomic_bam:
    input:
         bam = config['result_patterns']['se_atomic_bam_unsorted'],
    output:
          bam = config['result_patterns']['se_atomic_bam_sorted'],
    resources:
        # default -m: 768 Mb
        avg_mem = lambda wildcards, attempt: 1000 * attempt,
        mem_mb = lambda wildcards, attempt: 1500 * attempt,
        walltime_min = lambda  wildcards, attempt: 10 * attempt,
        attempt = lambda wildcards, attempt: attempt,
    params:
          name = 'sort_atomic_bam_{entity}_{sample}_{lib}_{uid}',
    threads: 8
    shell:
        """        
        if (( {resources.attempt} >= 4 )) ; then
            exit 1
        else
            samtools sort -o {output} -T {output} -@ 8 {input}
        fi
        """


rule bismark_pe_sort_atomic_bam:
    input:
         bam = config['result_patterns']['pe_atomic_bam_unsorted'],
    output:
          bam = config['result_patterns']['pe_atomic_bam_sorted'],
    resources:
         avg_mem = lambda wildcards, attempt: 1600 * attempt,
         mem_mb = lambda wildcards, attempt: 2000 * attempt,
         walltime_min = lambda wildcards, attempt: 20 * attempt,
         attempt = lambda wildcards, attempt: attempt,
    params:
          name = 'sort_atomic_bam_PE_{entity}_{sample}_{lib}_{uid}',
    threads: 8
    shell:
        """        
        if (( {resources.attempt} >= 4 )) ; then
            exit 1
        else
            samtools sort -o {output} -T {output} -@ 8 {input}
        fi
        """


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
    resources:
         # default -m: 768 Mb
         avg_mem = lambda wildcards, attempt: 6000,
         mem_mb = lambda wildcards, attempt: 8000,
         walltime_min = lambda wildcards, attempt: 20 * attempt,
         attempt = lambda wildcards, attempt: attempt,
    params:
          max_mem_gb = '8G',
          name = 'mdup_merge_{pairing}_{entity}_{sample}_{lib}',
          input_spec = lambda wildcards, input: ' '.join(f'INPUT={i}' for i in input)
    # java  -jar picard.jar
    # OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \ #changed from default of 100
    shell:
        """
        if (( {resources.attempt} >= 4 )) ; then
            exit 1
        else
            picard -Xmx{params.max_mem_gb} MarkDuplicates \
            {params.input_spec} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metrics} \
            CREATE_INDEX=false \
            TMP_DIR=/tmp
            
            samtools index {output.bam}
        fi
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
    resources:
         walltime_min = lambda wildcards, input, attempt: (30 if len(input) > 1 else 3) * attempt,
         avg_mem = lambda wildcards, input, attempt: (6000 if len(input) > 1 else 100),
         mem_mb = lambda wildcards, input, attempt: (8000 if len(input) > 1 else 200),
         attempt = lambda wildcards, attempt: attempt,
    params:
          name = 'merge-lib_{entity}_{sample}_{pairing}',
          # unpacking for resources / params: https://bitbucket.org/snakemake/snakemake/issues/471/unpack-for-params
    run:
        smk_wgbs.tools.abort_on_nth_attempt(resources.attempt, 4)
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
                bed.with_suffix('.p'),
                temp(bed.parent.joinpath('{pairing}_{entity}_{sample}' + f'_{methyldackel_motif}.bedGraph')),
            ]
    motif_params_str = ' '.join(motif_params_l)
    if not methyldackel_output:
        methyldackel_output = [sel_expand(config['result_patterns']['pe_or_se_mcalls_cytosines_per_motif'],
                               motif='UNUSED',
                               region='UNUSED')]
    return methyldackel_output, methyldackel_walltime, motif_params_str

methyldackel_output, methyldackel_walltime, motif_params_str = get_methyldackel_output_params(
        cg_walltime=60,
        chh_walltime=120,
        sampling_name='',
        calling_mode='full'
)
rule methyldackel_per_cytosine_per_motif:
    input:
          **methyldackel_input,
    output:
          methyldackel_output,
    resources:
         avg_mem = lambda wildcards, attempt: 4000,
         mem_mb = lambda wildcards, attempt: 8000,
         walltime_min = lambda wildcards, attempt: methyldackel_walltime * attempt,
         attempt = lambda wildcards, attempt: attempt,
    params:
          name = 'methyldackel_per_cytosine_{pairing}_{entity}_{sample}',
          prefix = lambda wildcards, output: Path(output[0]).parent.joinpath(
                  f'{wildcards.pairing}_{wildcards.entity}_{wildcards.sample}'),
          # TODO: is config directly available in rule?
          chromosomes = config['chromosomes'],
          motif_params = motif_params_str,
          sampling_options = '',
          calling_options = config['methyldackel_options'],
    # TODO: defaults for maxVariantFrac and minOppositeDepth
    threads: 4
    log:
       config['log_dir'] + '/methyldackel_per_cytosine_per_motif:{pairing}_{protocol}_{entity}_{sample}.log'
    run:
        smk_wgbs.tools.abort_on_nth_attempt(resources.attempt, 4)
        run_methyldackel(input, output, params, threads)

methyldackel_output, methyldackel_walltime, motif_params_str = get_methyldackel_output_params(
        cg_walltime=60,
        chh_walltime=120,
        sampling_name=sampling_name,
        calling_mode='sampled')

rule methyldackel_per_cytosine_per_motif_sampled:
    input:
         **methyldackel_input,
    output:
          methyldackel_output,
    resources:
             avg_mem = lambda wildcards, attempt: 3000,
             mem_mb = lambda wildcards, attempt: 4000,
             walltime_min = lambda wildcards, attempt: methyldackel_walltime * attempt,
             attempt = lambda wildcards, attempt: attempt,
    params:
          name = 'methyldackel_se_per_cytosine_{entity}_{sample}' + sampling_name,
          prefix = lambda wildcards, output: Path(output[0]).parent.joinpath(
                  f'{wildcards.pairing}_{wildcards.entity}_{wildcards.sample}'),
          # TODO: is config directly available in rule?
          chromosomes = config['chromosomes'],
          motif_params = motif_params_str,
          sampling_options = '-r 1:40000000-100000000',
          calling_options = config['methyldackel_options'],
    # TODO: defaults for maxVariantFrac and minOppositeDepth
    threads: 4
    log:
       config['log_dir'] + '/methyldackel_per_cytosine_per_motif_sampled:{pairing}_{protocol}_{entity}_{sample}.log'
    run:
        smk_wgbs.tools.abort_on_nth_attempt(resources.attempt, 4)
        run_methyldackel(input, output, params, threads)


rule methyldackel_merge_context:
    input:
        ref_genome_unconverted = ancient(config['genome_fa']),
        bedgraph = fp_to_bedgraph(config['result_patterns']['pe_or_se_mcalls_cytosines_per_motif']),
    output:
        bed = config['result_patterns']['pe_or_se_mcalls_merged_per_motif'],
        bedgraph = fp_to_bedgraph(config['result_patterns']['pe_or_se_mcalls_merged_per_motif']),
        parquet = fp_to_parquet(config['result_patterns']['pe_or_se_mcalls_merged_per_motif']),
        pickle = fp_to_pickle(config['result_patterns']['pe_or_se_mcalls_merged_per_motif'])
    resources:
         avg_mem = lambda wildcards, attempt: 800 * attempt,
         mem_mb = lambda wildcards, attempt: 1200 * attempt,
         walltime_min = lambda wildcards, attempt: 20 * attempt,
         attempt = lambda wildcards, attempt: attempt,
    params:
        chromosomes = config['chromosomes'],
        name = 'methyldackel_SE_mergeContext_{entity}_{sample}_{motif}{region}',
    log:
        config['log_dir'] + '/methyldackel_mergeContext_{pairing}_{entity}_{sample}_{motif}{region}'
    run:
        smk_wgbs.tools.abort_on_nth_attempt(resources.attempt, 4)
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
                parquet=output.parquet,
                pickle=output.pickle,
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
         pickle = fp_to_pickle(final_mcalls_cyt_bed),
    resources:
         avg_mem = lambda wildcards, attempt: 3000 * attempt,
         mem_mb = lambda wildcards, attempt: 5000 * attempt,
         walltime_min = lambda wildcards, attempt: 20 * attempt,
         attempt = lambda wildcards, attempt: attempt,
    params:
          name = 'merge-final-calls_{entity}_{sample}'
    run:
        smk_wgbs.tools.abort_on_nth_attempt(resources.attempt, 4)
        smk_wgbs.tools.merge_mcalls(input, output, wildcards)

final_mcalls_merged_bed = config['result_patterns']['final_mcalls_merged_per_motif']
rule final_mcalls_merged_motifs:
    input:
         bedgraph = fp_to_bedgraph(final_mcalls_cyt_bed),
         ref_genome_unconverted = ancient(config['genome_fa']),
    output:
         bed = final_mcalls_merged_bed,
         bedgraph = fp_to_bedgraph(final_mcalls_merged_bed),
         parquet = fp_to_parquet(final_mcalls_merged_bed),
         pickle = fp_to_pickle(final_mcalls_merged_bed),
    resources:
         avg_mem = lambda wildcards, attempt: 600 * attempt,
         mem_mb = lambda wildcards, attempt: 1000 * attempt,
         walltime_min = lambda wildcards, attempt: 20 * attempt,
         attempt = lambda wildcards, attempt: attempt,
    params:
          name = 'merge-context_final-mcalls_{entity}_{sample}_{motif}',
          chromosomes = config['chromosomes'],
    run:
        smk_wgbs.tools.abort_on_nth_attempt(resources.attempt, 4)
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
                parquet=output.parquet,
                pickle=output.pickle,
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
    resources:
         avg_mem = lambda wildcards, attempt: 20000,
         mem_mb = lambda wildcards, attempt: 24000,
         walltime_min = lambda wildcards, attempt: 24 * 60 * attempt,
         attempt = lambda wildcards, attempt: attempt,
    threads: 32
    # message: "Converting Genome into Bisulfite analogue"
    shell:
         """
         if (( {resources.attempt} >= 4 )) ; then
             exit 1
         else
             bismark_genome_preparation \
             --bowtie2 \
             --parallel 16 \
             --genomic_composition \
             --verbose \
             {input}
         fi
         """

# ==========================================================================================
# Trim the reads for adapter-ends and quality

# cutadapt speed up linear with cores
# parallelization requires pigz python package to be installed (by default if bismark is installed via conda)

cutadapt_phred_filter = config["cutadapt_phred_filter"]
if config['illumina_machine'] == 'nextseq':
    phred_filter_command = f'--nextseq-trim {cutadapt_phred_filter}'
else:
    phred_filter_command = f'-q {cutadapt_phred_filter}'


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
    resources:
         avg_mem = lambda wildcards, attempt: 800 * attempt,
         mem_mb = lambda wildcards, attempt: 1500 * attempt,
         walltime_min = lambda wildcards, attempt: 20 * attempt,
         attempt = lambda wildcards, attempt: attempt,
    params:
          name = 'trim_reads_{entity}_{sample}',
          phred_filter_command = phred_filter_command,
    threads: 8
    # currently, do not trim random hexamers priming sites \
    # note that this anyway does not take care of priming sites at the end of short reads \
    # -u 6
    # -U 6
    # discard read pair if any read is bad; given that we do SE alignment this may not be ideal
    # --pair-filter=any \
    shell:
         """
         if (( {resources.attempt} >= 4 )) ; then
             exit 1
         else
             cutadapt \
             {config[cutadapt_options]} \
             {params.phred_filter_command} \
             --cores {threads} \
             -a {config[read1_adapter_seq]} \
             -A {config[read2_adapter_seq]} \
             -o {output.fq_r1} \
             -p {output.fq_r2} \
             --pair-filter=any \
             {input.fq_r1} \
             {input.fq_r2} \
             > {log} 2>&1
         fi
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
    resources:
         avg_mem = lambda wildcards, attempt: 60000 * attempt,
         mem_mb = lambda wildcards, attempt: 80000 * attempt,
         walltime_min = lambda wildcards, attempt: 4 * 60 * attempt,
         attempt = lambda wildcards, attempt: attempt,
    params:
          name = 'nome_filtering_{entity}_{sample}',
          chromosomes = config['chromosomes'],
    log:
       config['log_dir'] + '/nome-filtering_{protocol}_{alignment_combi}_{entity}_{sample}{region}.log'
    run:
        smk_wgbs.tools.abort_on_nth_attempt(resources.attempt, 4)
        smk_wgbs.tools.nome_filtering(input, output, params)


rule samtools_stats:
    input:
         bam = config['result_patterns']['final_bam'],
         ref_genome_unconverted = ancient(config['genome_fa']),
    output:
          config['result_patterns']['final_bam_samtools_stats'],
    resources:
         avg_mem = lambda wildcards, attempt: 100 * attempt,
         mem_mb = lambda wildcards, attempt: 150 * attempt,
         walltime_min = lambda wildcards, attempt: 20 * attempt,
         attempt = lambda wildcards, attempt: attempt,
    params:
          name = 'samtools_stats_{entity}_{sample}',
    shell:
         """
         if (( {resources.attempt} >= 4 )) ; then
             exit 1
         else
             samtools stats \
             --insert-size 3000 \
             --ref-seq {input.ref_genome_unconverted} \
             {input.bam} \
             > {output[0]}
         fi
         """


rule samtools_flagstats:
    input:
         bam = config['result_patterns']['final_bam'],
    output:
          config['result_patterns']['final_bam_samtools_flagstats'],
    resources:
         avg_mem = lambda wildcards, attempt: 500 * attempt,
         mem_mb = lambda wildcards, attempt: 1000 * attempt,
         walltime_min = lambda wildcards, attempt: 20 * attempt,
         attempt = lambda wildcards, attempt: attempt,
    params:
          name = 'samtools_flagstats_{entity}_{sample}',
    shell:
         """
         if (( {resources.attempt} >= 4 )) ; then
             exit 1
         else
             samtools flagstat {input.bam} > {output[0]}
         fi
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
    resources:
         avg_mem = lambda wildcards, attempt: 5000,
         mem_mb = lambda wildcards, attempt: 10000,
         walltime_min = lambda wildcards, attempt: 20 * attempt,
         attempt = lambda wildcards, attempt: attempt,
    params:
          name = f'multiqc_{config["alignment_config_name"]}',
    run:
        smk_wgbs.tools.abort_on_nth_attempt(resources.attempt, 4)
        with open(output.file_list, 'wt') as fout:
            fout.write('\n'.join(input) + '\n')
        shell(f"""
              multiqc \
              --filename {output.report_html} \
              --force \
              --file-list {output.file_list}
              """
        )


fastq_screen_txts_trimmed = metadata_table[
    ['protocol', 'entity', 'sample', 'lib', 'uid', 'read_number']
].apply(
        lambda ser: (
            expand(config['result_patterns']['fastq_screen'], **ser)[0]
            + '/'
            + os.path.basename(expand(config['result_patterns']['trimmed_fastq'], **ser)[0])
              .replace('.fastq.gz', '_screen.txt'),
        ),
        axis=1
).to_list()
rule multiqc_read_level:
    input:
         # fastqc_reports,
         fastq_screen_txts_trimmed,
    output:
          report_html=config['experiment_result_patterns']['multiqc_read_level'],
          file_list= config['experiment_result_patterns']['multiqc_read_level'] + '_file-list.txt',
    resources:
         avg_mem = lambda wildcards, attempt: 2000 * attempt,
         mem_mb = lambda wildcards, attempt: 3000 * attempt,
         walltime_min = lambda wildcards, attempt: 20 * attempt,
         attempt = lambda wildcards, attempt: attempt,
    params:
          name = f'multiqc_{config["alignment_config_name"]}_read-level',
    run:
        smk_wgbs.tools.abort_on_nth_attempt(resources.attempt, 4)
        with open(output.file_list, 'wt') as fout:
            fout.write('\n'.join(input) + '\n')
        shell(f"""
              multiqc \
              --filename {output.report_html} \
              --force \
              --file-list {output.file_list}
              """
              )


meth_calling_qc_metadata_table = smk_wgbs.tools.create_mcalls_pattern_table(
        config, sampling_name
)


# TODO-protocol
mcalling_qc_report_done_files = (
    metadata_table[['entity', 'sample']]
        .drop_duplicates()
        .apply(lambda ser: expand(
            config['result_patterns']['meth_calling_qc_prefix'] + '.qc-done',
            **ser,
            protocol=config['protocol'],
            )[0],
               axis=1
               )
        .to_list()
)


rule meth_calling_qc:
    input:
        meth_calling_qc_metadata_table['fp'].to_list()
    output:
        done = touch(config['result_patterns']['meth_calling_qc_prefix'] + '.qc-done'),
    resources:
         avg_mem = lambda wildcards, attempt: 4000 * attempt,
         mem_mb = lambda wildcards, attempt: 4000 * attempt,
         walltime_min = lambda wildcards, attempt: 20 * attempt,
         attempt = lambda wildcards, attempt: attempt,
    params:
        metadata_table_expanded = lambda wildcards: smk_wgbs.tools.expand_mcalls_metadata_table(
                meth_calling_qc_metadata_table,
                entity=wildcards.entity,
                sample=wildcards.sample
        ),
        prefix=config['result_patterns']['meth_calling_qc_prefix'],
        name = 'mcalling_qc_{entity}_{sample}',
    run:
        smk_wgbs.tools.abort_on_nth_attempt(resources.attempt, 4)
        smk_wgbs.tools.compute_meth_calling_qc_stats(
                metadata_table=params.metadata_table_expanded,
                chrom_regions=config['chrom_regions'],
                prefix=params.prefix,
                keys=[
                    config['experiment_name'],
                    config['alignment_config_name'],
                    config['mcalling_config_name'],
                    wildcards.entity,
                    wildcards.sample
                ],
                names=[
                    'experiment_name',
                    'alignment_config_name',
                    'mcalling_config_name',
                    'entity',
                    'sample'
                ]
        )


rule concat_meth_calling_qc_data:
    input:
        mcalling_qc_report_done_files
    output:
        **{k: v for k, v in config['experiment_result_patterns'].items()
         if k.startswith('meth_calling')}
    resources:
             avg_mem = lambda wildcards, attempt: 24000 * attempt,
             mem_mb = lambda  wildcards, attempt: 32000 * attempt,
             walltime_min = lambda  wildcards, attempt: 59 * attempt,
             attempt = lambda wildcards, attempt: attempt,
    params:
          # TODO-protocol
          # Protect against expansion attempt for entity and sample
          prefix=lambda wildcards: sel_expand(
                  config['result_patterns']['meth_calling_qc_prefix'],
                  protocol=config['protocol']),
          name = 'concat_meth_calling_qc_data',
    run:
        smk_wgbs.tools.abort_on_nth_attempt(resources.attempt, 4)
        smk_wgbs.tools.concat_meth_calling_qc_data(
                metadata_table, params.prefix, output
        )


# TODO-protocol
# for each individual sample, ie unique (entity, sample), we take the final bam
# for samtools stats computation
unique_entity_sample_df = metadata_table[['entity', 'sample']].drop_duplicates()
samtools_stats_files = (
        unique_entity_sample_df
        .apply(lambda ser: expand(
            config['result_patterns']['final_bam_samtools_stats'],
            **ser,
            protocol=config['protocol'],
            alignment_combi=config['alignment_combi'],
            )[0],
               axis=1
               )
        .to_list()
)

rule collect_samtools_stats:
    input:
        samtools_stats_files,
    output:
        config['experiment_result_patterns']['samtools_stats']
    resources:
             avg_mem = lambda  wildcards, attempt: 4000 * attempt,
             mem_mb = lambda  wildcards, attempt: 5000 * attempt,
             walltime_min = lambda  wildcards, attempt: 20 * attempt,
             attempt = lambda wildcards, attempt: attempt,
    params:
          name = 'collect_samtools_stats',
    run:
        smk_wgbs.tools.abort_on_nth_attempt(resources.attempt, 4)
        smk_wgbs.tools.collect_samtools_stats(
                input_fps=input,
                output_p=output[0],
                keys=dict(experiment_name = config['experiment_name'],
                          alignment_config_name=config['alignment_config_name']),
                unique_entity_sample_df=unique_entity_sample_df,
        )


# TODO: can I rely on snakemake creating the folder for the done file
#   before executing the rule, so that fastq_screen can place files in this folder?
rule download_fastq_screen_indices:
    output:
        fastq_screen_index_conf = config['database_patterns']['fastq_screen_index_dir'] + '/fastq_screen.conf',
    resources:
             avg_mem = lambda wildcards, attempt: 4000 * attempt,
             mem_mb = lambda wildcards, attempt: 4000 * attempt,
             walltime_min = lambda wildcards, attempt: 2 * 60 * attempt,
             attempt = lambda wildcards, attempt: attempt,
    params:
        outdir=os.path.dirname(config['database_patterns']['fastq_screen_index_dir']),
        name = 'fastq-screen_get-genomes'
    shell:
        """
        if (( {resources.attempt} >= 4 )) ; then
            exit 1
        else
            if [ -d "{params.outdir}/FastQ_Screen_Genomes_Bisulfite" ]; then
                rm -r {params.outdir}/FastQ_Screen_Genomes_Bisulfite
            fi
            fastq_screen \
            --bisulfite \
            --get_genomes \
            --outdir {params.outdir}
        fi
        """

# TODO: conf file had to be manually edited and we use the copy in the smk repo for now
#    the fastq_screen.conf file on the fastq server is not correct currently (obviously a development version)

rule fastq_screen:
    input:
         # raw
         # fq = partial(smk_wgbs.tools.find_individual_fastqs, metadata_table=metadata_table),
         # trimmed
         fq = config['result_patterns']['trimmed_fastq'],
         fastq_screen_index_conf = config['database_patterns']['fastq_screen_index_dir'] + '/fastq_screen.conf',
    output:
        txt = config['result_patterns']['fastq_screen'] + '/{fastq_stem}_screen.txt'
    resources:
             avg_mem = lambda wildcards, attempt: 16000,
             mem_mb = lambda wildcards, attempt: 17000,
             walltime_min = lambda wildcards, attempt: 90 * attempt,
             attempt = lambda wildcards, attempt: attempt,
    threads: 6
    params:
        name = 'fastq-screen',
        outdir = lambda wildcards, output: os.path.dirname(output.txt)
    # --bowtie2 '--trim5 6' \
    # threads is ignored, put number of threads for one bismark run
    #   --threads {threads} \
    shell:
        """
        if (( {resources.attempt} >= 4 )) ; then
            exit 1
        else
            mkdir -p {params.outdir}/fastq-screen-tmp-files
            cd {params.outdir}/fastq-screen-tmp-files
            fastq_screen \
            --aligner bowtie2 \
            --bisulfite \
            --conf /home/kraemers/projects/smk_wgbs/smk_wgbs/fastq_screen.conf \
            --force \
            --nohits \
            --outdir {params.outdir} \
            --subset 100000 \
            {input.fq}
        fi
        """

# rule fastqc:
#     input:
#          unpack(partial(smk_wgbs.tools.find_fastqs, metadata_table=metadata_table)),
#     output:
#           qc_r1 = sel_expand(config['result_patterns']['fastq_fastqc'], read_number='1'),
#           fq_r1 = temp(sel_expand(config['result_patterns']['fastq_fastqc'], read_number='1').replace('_fastqc.zip', '.fastq.gz')),
#           qc_r2 = sel_expand(config['result_patterns']['fastq_fastqc'], read_number='2'),
#           fq_r2 = temp(sel_expand(config['result_patterns']['fastq_fastqc'], read_number='2').replace('_fastqc.zip', '.fastq.gz')),
#     params:
#           output_dir = lambda wildcards, output: Path(output[0]).parent,
#           avg_mem = lambda wildcards, attempt: 2000 * attempt,
#           mem_mb = lambda wildcards, attempt: 4000 * attempt,
#           walltime_min = lambda wildcards, attempt: 30 * attempt,
#           name = 'fastqc',
#     threads: 2
#     # -c {input.adapter_sequences_tsv} \
#     run:
#         smk_wgbs.tools.abort_on_nth_attempt(resources.attempt, 4)
#         import os
#         from pathlib import Path
#         os.symlink(input.fq_r1, output.fq_r1)
#         os.symlink(input.fq_r2, output.fq_r2)
#         # -o {params.output_dir} \
#         shell(
#                 """
#                 fastqc \
#                 --noextract \
#                 --threads {threads} \
#                 {output.fq_r1} {output.fq_r2}
#                 """
#          )
#          # fastqc_r1_fp = Path(params.output_dir).joinpath(Path(input.fq_r1).name.replace('.fastq.gz', '_fastqc.zip'))
#          # fastqc_r2_fp = Path(params.output_dir).joinpath(Path(input.fq_r2).name.replace('.fastq.gz', '_fastqc.zip'))
#          # fastqc_r1_fp.rename(output.qc_r1)
#          # fastqc_r2_fp.rename(output.qc_r2)


