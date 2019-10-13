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
--snakefile /home/kraemers/projects/smk_wgbs/smk_wgbs.smk \
--latency-wait 60 \
--cluster "bsub -R rusage[mem={params.avg_mem}] -M {params.max_mem} -n {threads} -J {params.name} -W {params.walltime} -o /home/kraemers/temp/logs/" \
--jobs 1000 \
--dryrun \

"""

import itertools
import os
import re
import glob
from pathlib import Path

import pandas as pd
from typing import Optional, Dict



# ==========================================================================================
# Generate methyl-converted version of the reference genome, if necessary:

# The path to the folder containing the genome to be bisulfite converted.
# The Bismark Genome Preparation expects one or more fastA files in the folder
# (with the file extension: .fa or .fasta (also ending in .gz)).
# Specifying this path is mandatory.


metadata_table_tsv = ''

lane_alignments = []
for _, row_ser in metadata_table.iterrows():
    lane_alignments.append(sel_expand(RESULT_PATHS['bam_per_lane_read'], **row_ser))


rule all:
    input:
         lane_alignments,
         # trimmed_fastq_files = trimmed_fastqs,


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
         fq = RESULT_PATHS['trimmed_fastq'],
         # not directly used, but still import files (parent directory path is bismark arg)
         refconvert_CT = config['genome_dir'] + "/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
         refconvert_GA = config['genome_dir'] + "/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
    output:
          bam = RESULT_PATHS['bam_per_lane_read'],
          report = RESULT_PATHS['bam_per_lane_read'].replace('.bam', '_SE_report.txt'),
    # non specified output files: will not be auto-removed
    params:
        genome_dir = config['genome_dir'],
        avg_mem = 50000,
        max_mem = 60000,
        name = 'bismark_alignment_{pid}_{cell_id}',
        walltime = '04:00',
        output_dir = lambda wildcards, output: str(Path(output.bam).parent),
        temp_dir = lambda wildcards, output: str(Path(output.bam).parent.joinpath('bismark_tempfiles')),
    threads: 24
    log:
        config['log_dir'] + '/bismark_lane-alignment_{pid}_{cell_id}_{run_id}_{as_id}_{lane_id}_{read_number}.log'
    # basename and multicore together are not possible
    shell:
        """
        bismark \
        --single_end {input.fq} \
        --parallel 4 \
        --local \
        --non_directional \
        --output_dir {params.output_dir} \
        --temp_dir {params.temp_dir}\
        {params.genome_dir} \
        > {log} 2>&1
        """


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
    message: "Converting Genome into Bisulfite analogue"
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
         fq_r1 = sel_expand(fastq_pattern, read_number='1'),
         fq_r2 = sel_expand(fastq_pattern, read_number='2'),
    output:
          fq_r1 = sel_expand(RESULT_PATHS['trimmed_fastq'], read_number='1'),
          fq_r2 = sel_expand(RESULT_PATHS['trimmed_fastq'], read_number='2'),
    log:
       config['log_dir'] + '/cutadapt_{pid}_{cell_id}_{run_id}_{as_id}_{lane_id}.log'
    message:
           "Trimming raw paired-end read data\n{input.fq_r1}\n{input.fq_r2}"
    params:
          avg_mem = 4000,
          max_mem = 6000,
          name = 'trim_reads_{pid}_{cell_id}',
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
         -a AGATCGGAAGAGCGGTTCAGCA \
         -A AGATCGGAAGAGCGTCGTGTAG \
         -o {output.fq_r1} \
         -p {output.fq_r2} \
         -q 10 \
         --minimum-length 30 \
         --max-n 0.3 \
         --pair-filter=any \
         {input.fq_r1} \
         {input.fq_r2} \
         > {log} 2>&1
         """
