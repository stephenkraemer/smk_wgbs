"""Tools for smk_wgbs users eg to get tables of produced files etc."""
import itertools
import os
from copy import deepcopy
from typing import Dict, Any, Optional, List
import pandas as pd
import smk_wgbs.utils as ut
import smk_wgbs.tools


# TODO: fix non-consecutive index
def create_mcalls_metadata_table(
    config: Dict,
    motifs: Optional[List[str]] = None,
    protocols: Optional[List[str]] = None,
    alignment_config_names: Optional[List[str]] = None,
    mcalling_config_names: Optional[List[str]] = None,
    entities: Optional[List[str]] = None,
    samples: Optional[List[str]] = None,
    workflow_versions: Optional[List[str]] = None,
    alignment_combis: Optional[List[str]] = None,
    merging_modes: Optional[List[str]] = None,
    sampling_names: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Create table of methylation calling results

    Note that the methylation calling results are placed in a results-per-sample folder,
    which does not contain the experiment name in its path as of now. Therefore, it is
    not possible to specify experiments for this function, and the metadata table has
    no column with the corresponding experiment. Maybe this should be changed in the future.


    Parameters
    ----------
    config: smk_wgbs base config dict, the config variables which are passed as args to this function are directly controlled by this function and will not be taken from the config dict. (Config dict is only used to look up path patterns actually).
    protocols: smk_wgbs config var, eg 'sctm'
    motifs: CG, CHG, CHH
    alignment_config_names: smk_wgbs config var
    mcalling_config_names: smk_wgbs config var
    entities: smk_wgbs config var, if defined: entities and samples must be List[str], ie nested entity lists are not allowed.
    samples: smk_wgbs config var, must be List[str] if defined
    workflow_versions: In the future, the workflow_version will likely be removed from the filepaths, and will then not need to be used here any more.
    alignment_combis: smk_wgbs config var, SE, PE, SE_PE or final. If final, points to the final workflow result, which may be either of the pairings. if SE_PE was used, SE and PE calls will also be available and can be retrieved via this arg.
    merging_modes: 'merged' or 'single'. Note: 'merged' is only available for CG and CHG
    sampling_names: sampling name, currently defined in wgbs_alignment.smk, so must be passed
        in addition to config dict

    Returns
    -------
    metadata table with fields:
                    "protocol",
                    "pairing",
                    "merging_mode",
                    "alignment_config_name",
                    "mcalling_config_name",
                    "entities",
                    "samples",
                    "workflow_version",
                    "path"
    """

    print(
        """
        WARNING
        This currently is working, but it produces a non-consecutive, non-sorted
        index. Consider fixing this before re-using this function.
        """
    )

    config = deepcopy(config)

    # currently, entity/sample definition by passing nested lists through entities is not allowed - only plain List[str]
    if "entities" in config:
        assert len(config["entities"] >= 1)
        assert all(isinstance(elem, str) for elem in config["entities"])
    if "samples" in config:
        assert len(config["samples"] >= 1)

    # provide mapping field name in path patterns -> List[field values]
    fields_ser = pd.Series(
        dict(
            zip(
                [
                    "protocol",
                    "motif",
                    "alignment_combi",
                    "merging_mode",
                    # alignment_config_name has field {config_name} atm, should be changed for clarity
                    "config_name",
                    "mcalling_config_name",
                    "entity",
                    "samples",
                    "workflow_version",
                    # sampling name is internally called region
                    "region",
                ],
                [
                    protocols,
                    motifs,
                    alignment_combis,
                    merging_modes,
                    alignment_config_names,
                    mcalling_config_names,
                    entities,
                    samples,
                    workflow_versions,
                    sampling_names,
                ],
            )
        )
    )

    # full region has empty field. If 'full' is specified, the user means ""
    fields_ser["region"] = [r if r != "full" else "" for r in fields_ser["region"]]
    fields_ser = fields_ser.dropna()

    # patterns are given as relative paths, prefix the results dir
    for pattern_name, pattern in config["result_patterns"].items():
        config["result_patterns"][pattern_name] = os.path.join(
            config["results_dir"], config["results_prefix"], pattern
        )

    # look up table for the appropriate result patterns (we retrieve the pickle path)
    patterns = {
        "SE": {
            "merged": config["result_patterns"]["pe_or_se_mcalls_merged_per_motif"],
            "single": config["result_patterns"]["pe_or_se_mcalls_cytosines_per_motif"],
        },
        "PE": {
            "merged": config["result_patterns"]["pe_or_se_mcalls_merged_per_motif"],
            "single": config["result_patterns"]["pe_or_se_mcalls_cytosines_per_motif"],
        },
        "SE_PE": {
            "merged": config["result_patterns"]["final_mcalls_merged_per_motif"],
            "single": config["result_patterns"]["final_mcalls_cytosines_per_motif"],
        },
        "final": {
            "merged": config["result_patterns"]["final_mcalls_merged_per_motif"],
            "single": config["result_patterns"]["final_mcalls_cytosines_per_motif"],
        },
    }
    # convert the bed file patterns taken from the config file to pickle
    # only pickle patterns are listed in the output df
    for alignment_mode, d1 in patterns.items():
        for merge_mode, d2 in d1.items():
            d1[merge_mode] = ut.fp_to_pickle(d2)

    # loop over all possible combinations of the given field values
    # for each combination, find all matching files via globbing if necessary
    # then concat and return
    df_l = []
    for combi in itertools.product(*fields_ser):
        curr_fields_ser = pd.Series(dict(zip(fields_ser.index, combi)))
        pattern = patterns[curr_fields_ser["alignment_combi"]][
            curr_fields_ser["merging_mode"]
        ]
        # alignment combi still needs to be expanded, don't drop it
        expanded_pattern = ut.sel_expand(
            pattern, **curr_fields_ser.drop(["merging_mode"])
        )
        if "{" in expanded_pattern:
            files_df = smk_wgbs.tools.create_metadata_table_from_file_pattern(
                expanded_pattern
            )
        else:
            files_df = pd.DataFrame({"path": expanded_pattern})
        for field_name, field_value in curr_fields_ser.iteritems():
            files_df[field_name] = field_value
        df_l.append(files_df)
    metadata_table = pd.concat(df_l)

    return metadata_table


# TODO-refactor: this is copy-pasted, then edited, from wgbs_alignment.smk.
#   this will fail if the path expansion logic in the workflow changes
def expand_config_relative_path_patterns(config: Dict[str, Any]) -> Dict:
    """Expand relative path patterns (*only results per pid currently*)

    Parameters
    ----------
    config: complete config dict used in one workflow run

    Returns
    -------
    Dict, deepcopy of config dict with expanded paths

    """

    config = deepcopy(config)

    # Fill prefix with alignment config_name and workflow_version
    results_prefix = ut.sel_expand(
        config["results_prefix"],
        config_name=config["alignment_config_name"],
        workflow_version=ut.find_workflow_version(),
    )

    # concatenate paths
    # currently only for results per pid
    for pattern_name, pattern in config["result_patterns"].items():
        pattern = pattern.replace(
            "{mcalling_config_name}", config["mcalling_config_name"]
        )
        config["result_patterns"][pattern_name] = os.path.join(
            config["results_dir"], results_prefix, pattern
        )

    return config
