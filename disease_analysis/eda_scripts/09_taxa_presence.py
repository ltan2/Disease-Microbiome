"""Taxa step 2: per-group presence/absence table (sample x taxon x coverage).

Writes results/<disease>/taxa/presence_<group>.tsv.gz: acc, taxon, one 0/1 column per coverage
(the format read by common.read_presence).
Samples with no taxa have no rows; group sizes come from groups.tsv.
Usage: 09_taxa_presence.py DISEASE
"""
import sys

import pandas as pd

import common

out = common.out_dir(sys.argv[1], "taxa")
groups = pd.read_csv(f"{out}/groups.tsv", sep="\t")
pres = pd.read_csv(common.YACHT_PRESENCE, sep="\t").merge(groups, on="acc")
for grp in common.GROUPS:
    pres[pres.group == grp].drop(columns="group").to_csv(f"{out}/presence_{grp}.tsv.gz", sep="\t", index=False)
print(pres.groupby("group").acc.nunique().reindex(common.GROUPS).to_string())
