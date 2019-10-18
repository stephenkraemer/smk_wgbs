import itertools
import os
import re
import glob
from pathlib import Path
from pkg_resources import resource_filename
import pandas as pd
from typing import Optional, Dict


def sel_expand(template, **kwargs):
    fields = kwargs.keys()
    values = [kwargs[f] for f in fields]
    values = [[val] if isinstance(val, (int, str)) else val
              for val in values]
    value_combinations = itertools.product(*values)
    def get_expanded_template(template, fields, comb):
        for field, value in zip(fields, comb):
            template = template.replace('{' + field + '}', value)
        return template
    res = [get_expanded_template(template, fields, comb) for comb in value_combinations]
    if len(res) == 1:
        return res[0]
    return res


def create_metadata_table_from_file_pattern(wildcard_pattern: str, field_constraints: Optional[Dict] = None):
    """Create metadata table for all files matching a snakemake-like pattern

    Output: metadata table with these columns:
      - one column per wildcard field
      - 'path' contains the full match for the snakemake-like pattern
      - fields in the output table are in the order of appearance from the filepath pattern

    Details:
      - fields may occur multiple times in the pattern
      - files are found by replacing each field with a '*' and using glob.glob
      - metadata are extracted using a regex constructed as follows:
        - the first occurence of each field is replaced with the default regex ('.+'),
          or the regex supplied via field_constraints
        - all following occurences of the field are replace with a backreference
          to the first match for the field


    """
    if field_constraints is None:
        field_constraints = {}

    field_names_set = set()
    all_field_name_occurences = re.findall(r'{(.*?)}', wildcard_pattern)
    field_names_in_order_of_appearance = [x for x in all_field_name_occurences
                                          if not (x in field_names_set or field_names_set.add(x))]

    assert set(field_constraints.keys()) <= field_names_set

    glob_pattern = re.sub(r'{(.+?)}', r'*', wildcard_pattern)
    glob_results = glob.glob(glob_pattern)
    if not glob_results:
        raise ValueError(f'Could not find any file matching:\n{glob_pattern}')

    regex_pattern = wildcard_pattern
    for field_name in field_names_set:
        if field_name in field_constraints:
            # replace first field with regex
            regex_pattern = regex_pattern.replace('{' + field_name + '}', f'(?P<{field_name}>{field_constraints[field_name]})', 1)
            # replace following fields with named backreference, if there are any
            regex_pattern = regex_pattern.replace('{' + field_name + '}', f'(?P={field_name})')
        else:
            regex_pattern = regex_pattern.replace('{' + field_name + '}', f'(?P<{field_name}>.+)', 1)
            # replace following fields with named backreference, if there are any
            regex_pattern = regex_pattern.replace('{' + field_name + '}', f'(?P={field_name})')

    metadata_df = pd.Series(glob_results).str.extract(regex_pattern)
    metadata_df['path'] = glob_results
    metadata_df = metadata_df[['path'] + field_names_in_order_of_appearance]

    return metadata_df


def get_snakefile_path() -> str:
    return resource_filename('smk_wgbs', 'smk_wgbs.smk')

# def get_demo_config_dict() -> dict:
#     return resource_filename

def find_workflow_version():
    return 'v0.1.0'

