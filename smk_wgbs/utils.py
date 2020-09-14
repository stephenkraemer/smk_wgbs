import re
import pickle
import itertools

def sel_expand(template, **kwargs):
    fields = kwargs.keys()
    values = [kwargs[f] for f in fields]
    values = [[val] if isinstance(val, (int, str)) else val for val in values]
    value_combinations = itertools.product(*values)

    def get_expanded_template(template, fields, comb):
        for field, value in zip(fields, comb):
            template = template.replace("{" + field + "}", value)
        return template

    res = [get_expanded_template(template, fields, comb) for comb in value_combinations]
    if len(res) == 1:
        return res[0]
    return res


def to_pickle(obj, fp, protocol=4):
    with open(fp, "wb") as fout:
        pickle.dump(obj, fout, protocol=protocol)


def from_pickle(fp):
    with open(fp, "rb") as fin:
        return pickle.load(fin)

def fp_to_pickle(fp, suffix="bed"):
    return re.sub(f".{suffix}$", ".p", fp)

def find_workflow_version():
    return "v0.1.0"
