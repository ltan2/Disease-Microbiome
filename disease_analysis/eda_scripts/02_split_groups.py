"""Split the labeled samples into the three groups for one disease (see common.split_groups).

Only samples with data for KIND are kept: a human-free sketch (sketch), a YACHT result (taxa) or
a KO profile (func).
Writes results/<disease>/<kind>/groups.tsv (acc, group), in the row order used by later steps.
Usage: 02_split_groups.py DISEASE sketch|taxa|func
"""
import os
import sys

import pandas as pd

import common

disease, kind = sys.argv[1], sys.argv[2]
groups = common.split_groups(disease)

if kind == "sketch":
    has_data = lambda a: os.path.exists(common.nohuman_path(a))
elif kind == "taxa":
    has_data = set(open(common.YACHT_SAMPLES).read().split()).__contains__
else:
    has_data = lambda a: os.path.exists(common.ko_path(a))

rows = [(a, g) for g in common.GROUPS for a in sorted(groups[g]) if has_data(a)]
df = pd.DataFrame(rows, columns=["acc", "group"])
df.to_csv(f"{common.out_dir(disease, kind)}/groups.tsv", sep="\t", index=False)
print(df.group.value_counts().reindex(common.GROUPS).to_string())
