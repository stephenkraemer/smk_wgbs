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
--snakefile /home/kraemers/projects/scwgbs_alignment/bismark_scwgbs_alignment.smk \
--latency-wait 60 \
--cluster "bsub -R rusage[mem={params.avg_mem}] -M {params.max_mem} -n {threads} -J {params.name} -W {params.walltime} -o /home/kraemers/temp/logs/" \
--jobs 1000 \
--dryrun \

"""



# ==========================================================================================
# Generate methyl-converted version of the reference genome, if necessary:

# The path to the folder containing the genome to be bisulfite converted.
# The Bismark Genome Preparation expects one or more fastA files in the folder
# (with the file extension: .fa or .fasta (also ending in .gz)).
# Specifying this path is mandatory.
GENOMEPATH = "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mbias/sandbox/genomes/GRCm38mm10_PhiX_Lambda_bismark"
LOGDIR = '/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/scwgbs_alignment/logs'

rule all:
    input:
         refconvert_CT = GENOMEPATH + "/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
         refconvert_GA = GENOMEPATH + "/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",

# indexing is run twice in parallel, so this will use twice the specified amount of threads
# code (with help at bottom: https://github.com/FelixKrueger/Bismark/blob/master/bismark_genome_preparation)
rule bismark_genome_preparation:
    input:
         ancient(GENOMEPATH)
    output:
          GENOMEPATH + "/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
          GENOMEPATH + "/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"
    params:
        avg_mem = 64000,
        max_mem = 125000,
        name = 'ref_conversion',
        walltime = '32:00',
    threads: 64
    log:
       LOGDIR + '/bismark_genome_preparation.log'
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





