"""
export PATH=/home/kraemers/projects/Bismark:$PATH

snakemake \
--snakefile /home/kraemers/projects/smk_wgbs/smk_wgbs/wgbs_alignment.smk \
--latency-wait 120 \
--configfile /home/kraemers/projects/smk_wgbs/doc/demo_config.yaml \
--cluster "bsub -R rusage[mem={params.avg_mem}] -M {params.max_mem} -n {threads} -J {params.name} -W {params.walltime} -o /home/kraemers/temp/logs/" \
--jobs 1000 \
--keep-going \
--rerun-incomplete \
--forcerun trim_reads_pe \
--dryrun \


--forcerun bismark_se_local_alignment_per_lane \
--forcerun nome_filtering \
"""

import time
import re
import shutil
import os
from pathlib import Path
import pandas as pd
import subprocess
import pandas as pd
from pandas.api.types import CategoricalDtype
from smk_wgbs import create_metadata_table_from_file_pattern, sel_expand, find_workflow_version
from smk_wgbs.tools import fp_to_parquet
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
    samples = ['blood41', 'blood48', 'blood79', 'blood12', 'blood22', 'blood65', 'blood59', 'blood63', 'blood35', 'blood30', 'blood80', 'blood49', 'blood78', 'blood28', 'blood76', 'blood32', 'blood72', 'blood73', 'blood91', 'blood40', 'blood90', 'blood9', 'blood74', 'blood71', 'blood57', 'blood75', 'blood87', 'blood18', 'blood56', 'blood68', 'blood10', 'blood70', 'blood84', 'blood62', 'blood26', 'blood42', 'blood85', 'blood89', 'blood29', 'blood50', 'blood45', 'blood60', 'blood83', 'bloodmethneg1']
    metadata_table = metadata_table.query('sample in @samples')
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
            pattern = config['result_patterns']['cg_full_parquet']
        else:
            pattern = config['result_patterns']['se_mcalls_cytosines_per_motif']
        targets.append(sel_expand(
                pattern,
                **row_ser,
                motif=motif,
                region=region,
        ))


rule all:
    input:
         targets,
         # trimmed_fastq_files = trimmed_fastqs,

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
rule bismark_se_local_alignment_per_lane:
    input:
         fq = config['result_patterns']['trimmed_fastq'],
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
        name = 'bismark_alignment_{entity}_{sample}',
        walltime = '08:00',
        # walltime = '00:30',
    threads: 6
    log:
        config['log_dir'] + '/bismark_atomic-alignment_{protocol}_{entity}_{sample}_{uid}_{lib}_R{read_number}.log'
    # basename and multicore together are not possible
    # --score_min G,10,2 \
    # --local
    # --parallel 4 \
    run:
        # temp_dir
        # - need different temp_dir per atomic alignment
        # - rule removes the atomic temp_dir at the end of the rule using the temp()
        # directive. This is more robust than removing the directory inside the rule
        # (in case of crashed jobs)
        # - if bismark changes its behavior to completely remove the temp_dir at the end
        # of its run, this rule will fail, because the temp_dir is an expected (temp) output
        from pathlib import Path
        import re

        output_dir = Path(output.bam).parent
        shell(
           """
           bismark \
           --single_end {input.fq} \
           --non_directional \
           --local \
           --output_dir {output_dir} \
           --temp_dir {output.temp_dir}\
           {params.genome_dir} \
           > {log} 2>&1
           """
        )

        # obtain paths where bismark will place files
        # bismark hardcodes output files basenames
        # (--basename cannot be used together with multithreading atm)
        stem = re.sub('.fastq(.gz)?$', '', Path(input.fq).name)
        bismark_bam_fp = output_dir.joinpath(stem + '_bismark_bt2.bam')
        bismark_se_report_fp = output_dir.joinpath(stem + '_bismark_bt2_SE_report.txt')
        # Move results from hard coded paths to user-specified paths
        bismark_bam_fp.rename(output.bam)
        bismark_se_report_fp.rename(output.bismark_report)


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

def find_atomic_bams_for_library(wildcards):
    df = metadata_table
    matches_library_and_read = (
            (df['protocol'] == wildcards.protocol)
            & (df['entity'] == wildcards.entity)
            & (df['sample'] == wildcards.sample)
            & (df['lib'] == wildcards.lib)
            # & (df['read_number'] == wildcards.read_number)
    )
    ser = (df
           .loc[matches_library_and_read, ['protocol', 'entity', 'sample', 'lib', 'uid', 'read_number']]
           .apply(lambda ser: expand(config['result_patterns']['se_atomic_bam_sorted'], **ser)[0],
                  axis=1)
           )
    l = ser.to_list()
    return l


rule bismark_se_merge_library_mdup:
    input:
        find_atomic_bams_for_library
    output:
         bam = config['result_patterns']['se_library_bam'],
         bai = config['result_patterns']['se_library_bam'] + '.bai',
         metrics = re.sub(
                 '.bam$',
                 '_mdup-metrics.txt',
                 config['result_patterns']['se_library_bam']),
    params:
          # default -m: 768 Mb
          avg_mem = 8 * 1000,
          max_mem = 8 * 1600,
          max_mem_gb = '13G',
          walltime = '00:30',
          name = 'mdup_merge_{entity}_{sample}_{lib}',
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
           .apply(lambda ser: expand(config['result_patterns']['se_library_bam'], **ser)[0],
                  axis=1)
           )
    l = ser.to_list()
    return l
# TODO: index
rule bismark_se_merge_libraries:
    input:
        find_library_bams
    output:
        bam = config['result_patterns']['se_bam'],
        bai = config['result_patterns']['se_bam'] + '.bai',
    # TODO: function which returns **resources and returns minimal resources if only hard linking is required
    params:
          avg_mem = 6000,
          max_mem = 8000,
          name = 'merge-lib_{entity}_{sample}',
          walltime = '00:30',
    run:
        # hard link if only one library, else merge
        import os
        import subprocess
        if len(input) == 1:
            os.link(input[0], output.bam)
            os.link(input[0] + '.bai', output.bam + '.bai')
        else:
            subprocess.run('samtools merge -f -r'.split() + output + input)
            subprocess.run(['samtools', 'index', output.bam])



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


def run_methyldackel(input, output, params, threads):

    print('Running MethylDackel')
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
    print('Done with MethylDackel')

    chrom_dtype = CategoricalDtype(categories=params.chromosomes, ordered=True)

    assert len(output) % 4 == 0
    for i in range(int(len(output) / 4)):
        bed, bedgraph, parquet, tmp_bedgraph  = output[4*i:(4*i)+4]
        # convert bedgraph to bed and parquet
        shutil.copy(tmp_bedgraph, bedgraph)
        print('Create parquet and BED file')
        df = pd.read_csv(
                bedgraph,
                # '/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/scwgbs_alignment/results-per-entity/hsc_pr-scbs_ex-test-1_cls-10_1/blood/meth-calls/partial-meth-calls/SE_hsc_pr-scbs_ex-test-1_cls-10_1_blood_per-cytosine_CpG.bedGraph',
                sep='\t',
                skiprows=1,
                names=['#Chromosome', 'Start', 'End', 'beta_value', 'n_meth', 'n_unmeth'],
                dtype={'#Chromosome': chrom_dtype}
        )
        # Discard unwanted chromosomes/or more commonly: contigs (not in config['chromosomes'])
        df = df.dropna(subset=['#Chromosome'])
        df['n_total'] = df.eval('n_meth + n_unmeth')
        df['beta_value'] = df.eval('n_meth / n_total')
        # Add motif column
        # TODO: improve
        if 'CpG' in bed:
            motif = 'CG'
        elif 'CHG' in bed:
            motif  = 'CHG'
        elif 'CHH' in bed:
            motif = 'CHH'
        else:
            motif = 'NA'
        df['motif'] = motif

        # Final column order
        df = df[['#Chromosome', 'Start', 'End', 'motif', 'beta_value', 'n_total', 'n_meth']]

        # Save to BED and parquet
        df.to_csv(bed, sep='\t', header=True, index=False)
        df.rename(columns={'#Chromosome': 'Chromosome'}).to_parquet(parquet)
        print('done')
        # else: continue

methyldackel_input = dict(
    bam = config['result_patterns']['se_bam'],
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
                    config['result_patterns']['se_mcalls_cytosines_per_motif'],
                    motif=motif,
                    region=f'_sampled-{sampling_name}' if sampling_name else '',
            ))
            methyldackel_motif = motif if motif != 'CG' else 'CpG'
            methyldackel_output += [
                bed,
                bed.with_suffix('.bedGraph'),
                bed.with_suffix('.parquet'),
                temp(bed.parent.joinpath('SE_{entity}_{sample}' + f'_{methyldackel_motif}.bedGraph')),
            ]
    motif_params_str = ' '.join(motif_params_l)
    if not methyldackel_output:
        methyldackel_output = [sel_expand(config['result_patterns']['se_mcalls_cytosines_per_motif'],
                               motif='UNUSED',
                               region='UNUSED')]
    return methyldackel_output, methyldackel_walltime, motif_params_str

methyldackel_output, methyldackel_walltime, motif_params_str = get_methyldackel_output_params(
        cg_walltime='02:00',
        chh_walltime='10:00',
        sampling_name='',
        calling_mode='full'
)
rule methyldackel_se_CG_per_cytosine:
    input:
          **methyldackel_input,
    output:
          methyldackel_output,
    params:
          avg_mem = 6000,
          max_mem = 10000,
          name = 'methyldackel_se_per_cytosine_{entity}_{sample}',
          walltime = methyldackel_walltime,
          prefix = lambda wildcards, output: Path(output[0]).parent.joinpath(
                  f'SE_{wildcards.entity}_{wildcards.sample}'),
          # TODO: is config directly available in rule?
          chromosomes = config['chromosomes'],
          motif_params = motif_params_str,
          sampling_options = '',
    # TODO: defaults for maxVariantFrac and minOppositeDepth
    threads: 4
    log:
       config['log_dir'] + '/methyldackel_se_CG_per_cytosine:{protocol}_{entity}_{sample}.log'
    run:
        run_methyldackel(input, output, params, threads)

methyldackel_output, methyldackel_walltime, motif_params_str = get_methyldackel_output_params(
        cg_walltime='00:30',
        chh_walltime='02:00',
        sampling_name=sampling_name,
        calling_mode='sampled')

rule methyldackel_se_CG_per_cytosine_sampled:
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
                  f'SE_{wildcards.entity}_{wildcards.sample}'),
          # TODO: is config directly available in rule?
          chromosomes = config['chromosomes'],
          motif_params = motif_params_str,
          sampling_options = '-r 1:40000000-100000000',
    # TODO: defaults for maxVariantFrac and minOppositeDepth
    threads: 4
    log:
       config['log_dir'] + '/methyldackel_se_CG_per_cytosine:{protocol}_{entity}_{sample}.log'
    run:
        run_methyldackel(input, output, params, threads)




# --nOT   6,0,6,0 \
# --nOB   6,0,6,0 \
# --nCTOT 6,0,6,0 \
# --nCTOB 6,0,6,0 \

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
def find_fastqs(wildcards):
    df = metadata_table
    is_in_selected_pair = (
            (df['protocol'] == wildcards.protocol)
            & (df['entity'] == wildcards.entity)
            & (df['sample'] == wildcards.sample)
            & (df['lib'] == wildcards.lib)
            & (df['uid'] == wildcards.uid)
    )
    fastq_r1_ser = df.loc[is_in_selected_pair & (df['read_number'] == '1'), 'path']
    assert fastq_r1_ser.shape[0] == 1
    fastq_r2_ser = df.loc[is_in_selected_pair & (df['read_number'] == '2'), 'path']
    assert fastq_r2_ser.shape[0] == 1
    return {'fq_r1': fastq_r1_ser.iloc[0], 'fq_r2': fastq_r2_ser.iloc[0]}

rule trim_reads_pe:
    input:
        unpack(find_fastqs)
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
         -q 10 \
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
                    config['result_patterns']['se_mcalls_cytosines_per_motif'],
                    motif=motif,
                    region=region
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
       config['log_dir'] + '/nome-filtering_{protocol}_{entity}_{sample}{region}.log'
    run:
        smk_wgbs.tools.nome_filtering(input, output, params)
