# * Example from python

import pandas as pd
import snakemake
import yaml
from smk_wgbs import get_files_df, get_snakefile_path, sel_expand

# ** Generate config file
# Can be based on smk_wgbs/doc/demo_config.yaml
# - copy to appropriate location and edit as needed
# - or load and modify in python as needed

config_yaml = 'smk_wgbs/doc/demo_config.yaml'
with open(config_yaml) as fin:
    config = yaml.load(fin)

# ** Automatically generate metadata table

view_by_pid_dir = "/icgc/dkfzlsdf/project/mouse_hematopoiesis_PBAT/sequencing/whole_genome_bisulfite_tagmentation_sequencing/view-by-pid"

fastq_pattern = view_by_pid_dir + "/{pid}/{cell_id}/paired/run{run_id}/sequence/AS-{as_id}-LR-{lane_id}_R{read_number}.fastq.gz"

pids = ['lsk150_pr-hiera_ex-scnmt-1_cls-1']

metadata_table = pd.concat([
    get_files_df(sel_expand(fastq_pattern, pid=pid)).assign(pid=pid)
    for pid in pids
], axis=0)

metadata_table = metadata_table.iloc[0:2].copy()

metadata_table.to_csv(config['metadata_table_tsv'],
                      header=True, index=False, sep='\t')

snakemake.snakemake(
        snakefile=get_snakefile_path(),
        configfile=config_yaml,
)
