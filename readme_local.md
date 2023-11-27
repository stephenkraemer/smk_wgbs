# On the environment files

- smk_wgbs is known to fail with more recent versions of htslib, methyldackel and samtools, this hasn't been adressed yet, so usable versions of these libraries are fixed also in the dev_env.yaml
   - see also issue: "Improve current env file: problems with methyldackel 0.5.x, htslib ~1.10, samtools ~1.10; had to fix the methyldackel=0.4, samtools=1.9 and htslib=1.9 versions"
   - more deprectation issues could have accumulated by now, the codebase hasn't been updated for a while now, eg. newer snakemake versions would also be a candidate for problems
- the lock file smk_wgbs/smk-wgbs_0.2-env_locked.yaml and the env_dev.yaml file smk_wgbs/smk_wgbs_dev_env.yaml have been created at the same day, so the lock file is current and matched to the dev_env.yaml file (esp. in terms of the htslib, samtools and methyldackel versions)

# Installation

The easiest way to use the workflow is as a python package. To set this up:

1. Create a conda environment
    - install miniconda as described here:
      https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
    - using mamba is recommended, to install it: `conda install mamba`
    - then create the workflow environment: `mamba env create -f smk_wgbs/smk-wgbs_0.2-env_locked.yaml`
      - for development, consider installing smk_wgbs/smk_wgbs_dev_env.yaml (see section on environments above)
2. There is no pypi/conda package yet (coming soon). To install the package
    - clone the repository
    - run `pip install /path/to/smk_wgbs`
