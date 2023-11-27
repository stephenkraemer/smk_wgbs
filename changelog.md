0.2.0

- fix: log directory for bismark genome preparation is not defined, resulting in bash error
  when trying to execute the rule. Fixed by removing redirection to log file from bash code running bismark genome preparation.
- use 32 cores for bismark genome preparation instead of 64, to make queuing easier - the high core number was 
  just used during development. The walltime has consequently been increased to 24h. Memory has been left as is for now,
  can probably be reduced in the future.
- deleted apparently unused params section from bismark genome preparation rule: 
  params: name = 'ref_conversion'
- fix: add check=True to run_methyldackel
- fix: in tools.prepare_index_files, replace 'is' with '=='
- fix: increase mem for bismark to 18 GiB (from 16 GB), just to be safe
- add smk_wgbs_dev_env.yaml
- add env lockfile smk-wgbs_0.2-env_locked.yaml