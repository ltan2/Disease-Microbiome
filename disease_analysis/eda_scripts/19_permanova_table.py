"""One summary table of the PERMANOVA/PERMDISP results of a disease (KO Aitchison + taxa Jaccard per coverage).

R^2 (share of variation explained by group) is computed from the pseudo-F:
R^2 = F (k-1) / (F (k-1) + n - k), with k groups and n samples.
Writes results/<disease>/permanova_table.tsv and prints it as a markdown table.
Usage: 19_permanova_table.py DISEASE
"""
import os
import sys

import pandas as pd

import common

out = os.path.dirname(common.out_dir(sys.argv[1], "func"))
files = {"KO (Aitchison)": f"{out}/func/aitchison_permanova.tsv"}
files |= {f"Taxa (Jaccard, cov {c})": f"{out}/taxa/permanova_{c}.tsv" for c in common.COVERAGES}

rows = []
for data, f in files.items():
    if not os.path.exists(f):
        continue
    r = pd.read_csv(f, sep="\t")
    k = r.groups.str.count(" vs ") + 1
    r["R2"] = r.F * (k - 1) / (r.F * (k - 1) + r["sample size"] - k)
    a, b = r[r.test == "PERMANOVA"].set_index("groups"), r[r.test == "PERMDISP"].set_index("groups")
    for g in a.index:
        rows.append({"data": data, "comparison": "all 3 groups" if g.count(" vs ") == 2 else g,
                     "n": a.loc[g, "sample size"], "PERMANOVA F": round(a.loc[g, "F"], 1),
                     "R2": round(a.loc[g, "R2"], 3), "PERMANOVA p": a.loc[g, "p-value"],
                     "PERMDISP F": round(b.loc[g, "F"], 2), "PERMDISP p": round(b.loc[g, "p-value"], 3)})
df = pd.DataFrame(rows)
df.to_csv(f"{out}/permanova_table.tsv", sep="\t", index=False)

print("| " + " | ".join(df.columns) + " |\n|" + "---|" * df.shape[1])
for r in df.itertuples(index=False):
    print("| " + " | ".join(f"{v:,}" if isinstance(v, int) else str(v) for v in r) + " |")
