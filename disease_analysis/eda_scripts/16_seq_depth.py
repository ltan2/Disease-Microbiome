"""Sequencing depth (SRA mbases, before human removal) of every sample in one data kind's groups.

Samples without an mbases value are left out.
Writes results/<disease>/<kind>/depth.tsv (acc, group, mbases).
Usage: 16_seq_depth.py DISEASE sketch|taxa|func
"""
import sys

import duckdb
import pandas as pd

import common

disease, kind = sys.argv[1], sys.argv[2]
out = common.out_dir(disease, kind)
groups = pd.read_csv(f"{out}/groups.tsv", sep="\t")

con = duckdb.connect(common.METADATA_DB, read_only=True)
con.register("groups", groups)
df = con.sql("""SELECT acc, "group", mbases FROM groups
                JOIN (SELECT DISTINCT ON (acc) acc, mbases FROM sra_metadata) USING (acc)
                WHERE mbases IS NOT NULL""").df()
df.to_csv(f"{out}/depth.tsv", sep="\t", index=False)
print(f"{len(df)} of {len(groups)} samples have mbases")
