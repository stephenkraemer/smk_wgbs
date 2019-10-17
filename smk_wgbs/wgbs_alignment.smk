"""
reference genome creation	20

trimming	120

chromosome length creation	20
se alignment	20
se duplicate removal	20
sorting	20
flagstats	20
fastqc before trimming	20
fastqc after trimming	10
methylation calling with bismark	45
reindex methylation calls	90
save as pickle / parquet / hdf5	10
7.25

snakemake \
--snakefile /home/kraemers/projects/smk_wgbs/smk_wgbs/wgbs_alignment.smk \
--latency-wait 60 \
--cluster "bsub -R rusage[mem={params.avg_mem}] -M {params.max_mem} -n {threads} -J {params.name} -W {params.walltime} -o /home/kraemers/temp/logs/" \
--configfile /home/kraemers/projects/smk_wgbs/doc/demo_config.yaml \
--forcerun bismark_se_local_alignment_per_lane \
--until bismark_se_local_alignment_per_lane \
--jobs 1000 \

--dryrun \
--forcerun trim_reads_pe \


--forcerun methyldackel_se_CG_per_cytosine \

--forcerun bismark_se_merge_library_mdup \
--forcerun bismark_se_merge_libraries \



--forcerun bismark_se_local_alignment_per_lane \




"""

import time
import re
import os
from pathlib import Path
import pandas as pd
from smk_wgbs import create_metadata_table_from_file_pattern, sel_expand


# For development
# import yaml
# config_yaml = '/home/kraemers/projects/smk_wgbs/doc/demo_config.yaml'
# with open(config_yaml) as fin:
#     config = yaml.load(fin, Loader=yaml.FullLoader)


# convert relative to absolute paths
for pattern_name, pattern in config['result_patterns'].items():
    config['result_patterns'][pattern_name] = (
        pattern if os.path.isabs(pattern)
        else os.path.join(config['rpe_dir'], pattern))

# for pattern_name, pattern in config['result_patterns'].items():
#     config['result_patterns'][pattern_name] = pattern.replace('{uid}', '_'.join(colname + '-{' + colname + '}' for colname in uid_columns))

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
    if 'lib' not in metadata_table:
        metadata_table['lib'] = '1'
    metadata_table = metadata_table.astype(str)
    timestamp = time.strftime('%d-%m-%y_%H:%M:%S')
    metadata_table.to_csv(
            Path(config['metadata_table_tsv']).with_suffix(f'.{timestamp}.tsv'),
            header=True, index=False, sep='\t')
    # REMOVE
    metadata_table = metadata_table.iloc[0:2].copy()
    # \REMOVE
else:
    metadata_table = pd.read_csv(config['metadata_table_tsv'], sep='\t', header=0, index=False)
uid_columns = metadata_table.columns[
    ~metadata_table.columns.isin(['entity', 'sample', 'read_number', 'path', 'lib'])]

# add read pair UID field to metadata table
metadata_table['uid'] = metadata_table[uid_columns].apply(
        lambda ser: '_'.join(f'{uid_col}-{ser[uid_col]}' for uid_col in uid_columns), axis=1)


# Targets
# ======================================================================
wildcard_constraints:
    lib = '[^_]+'

lane_alignments = []
for _, row_ser in metadata_table.iterrows():
    lane_alignments.append(sel_expand(config['result_patterns']['mcalls_se_cg_per_cyt'].replace('.bed', '.bedGraph'), **row_ser))


rule all:
    input:
         lane_alignments,
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
          bam = config['result_patterns']['atomic_bam_unsorted'],
          report = config['result_patterns']['atomic_bam_unsorted'].replace('.bam', '_SE_report.txt'),
    # non specified output files: will not be auto-removed
    params:
        genome_dir = config['genome_dir'],
        avg_mem = 50000,
        max_mem = 60000,
        name = 'bismark_alignment_{entity}_{sample}',
        walltime = '04:00',
        output_dir = lambda wildcards, output: str(Path(output.bam).parent),
        temp_dir = lambda wildcards, output: str(Path(output.bam).parent.joinpath('bismark_tempfiles')),
    threads: 24
    log:
        config['log_dir'] + '/bismark_atomic-alignment_{entity}_{sample}_{uid}_{lib}_R{read_number}.log'
    # basename and multicore together are not possible
    # --score_min G,10,2 \
    # --local
    shell:
        """
        bismark \
        --single_end {input.fq} \
        --parallel 4 \
        --non_directional \
        --local \
        --upto 100000 \
        --output_dir {params.output_dir} \
        --temp_dir {params.temp_dir}\
        {params.genome_dir} \
        > {log} 2>&1
        """

rule bismark_se_sort_atomic_bam:
    input:
         bam = config['result_patterns']['atomic_bam_unsorted'],
    output:
          bam = config['result_patterns']['atomic_bam_sorted'],
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
            (df['entity'] == wildcards.entity)
            & (df['sample'] == wildcards.sample)
            & (df['lib'] == wildcards.lib)
            # & (df['read_number'] == wildcards.read_number)
    )
    ser = (df
           .loc[matches_library_and_read, ['entity', 'sample', 'lib', 'uid', 'read_number']]
           .apply(lambda ser: expand(config['result_patterns']['atomic_bam_sorted'], **ser)[0],
                  axis=1)
           )
    l = ser.to_list()
    return l


rule bismark_se_merge_library_mdup:
    input:
        find_atomic_bams_for_library
    output:
         bam = config['result_patterns']['library_bam'],
         bai = config['result_patterns']['library_bam'] + '.bai',
         metrics = re.sub(
                 '.bam$',
                 '_mdup-metrics.txt',
                 config['result_patterns']['library_bam']),
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
            (df['entity'] == wildcards.entity)
            & (df['sample'] == wildcards.sample)
    )
    ser = (df
           .loc[matches_library, ['entity', 'sample', 'lib']]
           .drop_duplicates()
            # expand always returns list, get string out of list of length 1
           .apply(lambda ser: expand(config['result_patterns']['library_bam'], **ser)[0],
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

rule methyldackel_se_CG_per_cytosine:
    input:
         bam = config['result_patterns']['se_bam'],
         ref_genome_unconverted = ancient(config['genome_fa']),
    output:
        bedgraph = str(Path(config['result_patterns']['mcalls_se_cg_per_cyt']).with_suffix('.bedGraph')),
        bed = config['result_patterns']['mcalls_se_cg_per_cyt'],
        parquet = str(Path(config['result_patterns']['mcalls_se_cg_per_cyt']).with_suffix('.parquet')),
    params:
          avg_mem = 6000,
          max_mem = 10000,
          name = 'methyldackel_se_CG_per_cytosine_{entity}_{sample}',
          walltime = '02:00',
          prefix = lambda wildcards, output: output.bedgraph.replace('_CpG.bedGraph', ''),
          # TODO: is config directly available in rule?
          chromosomes = config['chromosomes'],
    # TODO: defaults for maxVariantFrac and minOppositeDepth
    threads: 4
    log:
       config['log_dir'] + '/methyldackel_se_CG_per_cytosine:_{entity}_{sample}.log'
    run:
        import subprocess
        import pandas as pd
        from pandas.api.types import CategoricalDtype

        print('Running MethylDackel')
        subprocess.run(
            f"""
            MethylDackel extract \
            --ignoreFlags 3840 \
            --requireFlags 0 \
            -o {params.prefix} \
            -q 15 \
            -p 15 \
            -@ {threads} \
            {input.ref_genome_unconverted} \
            {input.bam}
            """.split()
        )
        print('Done with MethylDackel')

        # NOTE: parquet does not retain the categorical dtype atm
        print('Create parquet and BED file')
        chrom_dtype = CategoricalDtype(categories=params.chromosomes, ordered=True)
        df = pd.read_csv(
                # output.bedGraph,
                '/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/scwgbs_alignment/results-per-entity/lsk150_pr-hiera_ex-scnmt-1_cls-1/blood31/meth-calls/partial-meth-calls/SE_lsk150_pr-hiera_ex-scnmt-1_cls-1_blood31_per-cytosine_CpG.bedGraph',
                sep='\t',
                skiprows=1,
                names=['#Chromosome', 'Start', 'End', 'beta_value', 'n_meth', 'n_unmeth'],
                dtype={'#Chromosome': chrom_dtype}
        )
        # Discard unwanted chromosomes/or more commonly: contigs (not in config['chromosomes'])
        df = df.dropna(subset=['#Chromosome'])
        df['n_total'] = df.eval('n_meth + n_unmeth')
        df['beta_value'] = df.eval('n_meth / n_total')
        df = df[['#Chromosome', 'Start', 'End', 'beta_value', 'n_total', 'n_meth']]

        df.to_csv(output.bed, sep='\t', header=True, index=False)
        df.rename(columns={'#Chromosome': 'Chromosome'}).to_parquet(output.parquet)
        print('done')




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
    log:
       config['log_dir'] + '/bismark_genome_preparation.log'
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
    is_in_selected_pair = ((df['entity'] == wildcards.entity)
                           & (df['sample'] == wildcards.sample)
                           & (df['lib'] == wildcards.lib)
                           & (df['uid'] == wildcards.uid))
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
       config['log_dir'] + '/cutadapt_{entity}_{sample}_{uid}_{lib}.log'
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
