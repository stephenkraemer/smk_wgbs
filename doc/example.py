# * Example from python

import pandas as pd
import snakemake
import yaml
from smk_wgbs import create_metadata_table_from_file_pattern, get_snakefile_path, sel_expand

# ** Generate config file
# 1. From a new dict
# 2. Or based on smk_wgbs/doc/demo_config.yaml
# - copy to appropriate location and edit as needed
# - or load and modify in python as needed
# For the example, we use the provided demo configfile


# ** Generate metadata table and save to file


metadata_table = pd.concat([
    create_metadata_table_from_file_pattern(sel_expand(fastq_pattern, pid=pid)).assign(pid=pid)
    for pid in pids
], axis=0)

metadata_table = metadata_table.iloc[0:2].copy()

# as specified in demo_config.yaml
metadata_table_tsv = "/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/scwgbs_alignment/alignment-metadata/bravo-pilot.tsv"
metadata_table.to_csv(metadata_table_tsv,
                      header=True, index=False, sep='\t')

# ** Run from python

snakemake.snakemake(
        snakefile=get_snakefile_path(),
        configfile='/home/kraemers/projects/smk_wgbs/doc/demo_config.yaml',
)
