"""One group's samples, a per-sample value against sequencing depth, colored by BioProject.

The 8 projects with the most runs get their own color (legend: body site, runs, within-project Spearman
rho with depth); the rest are gray. The title gives the group's overall rho, the rho between projects
(project medians) and the pooled rho within projects (depth and value ranked inside each colored project).
Also writes OUT_projects.tsv (bioproject, site, runs, median depth, median value, rho).

Usage: plot_depth_projects.py IN.tsv COLUMN GROUP OUT.png TITLE [--ylabel TEXT] [--xlabel TEXT]
  IN.tsv  from 36_per_mbp.py (acc, group, mbases, ...); COLUMN is any of its columns
"""
import argparse

import duckdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import common

# categorical palette, fixed order (dataviz reference palette, light mode)
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, MUTED, GRID, OTHER = "#2b2b28", "#6b6a64", "#e5e4df", "#b5b4ad"

ap = argparse.ArgumentParser()
for a in ["tsv", "column", "group", "out", "title"]:
    ap.add_argument(a)
ap.add_argument("--ylabel")
ap.add_argument("--xlabel", default="sequencing depth (Mbases, before human removal)")
args = ap.parse_args()
col = args.column

df = pd.read_csv(args.tsv, sep="\t", dtype={"acc": str})
df = df[(df.group == args.group) & (df[col] > 0) & (df.mbases > 0)]
con = duckdb.connect(common.METADATA_DB, read_only=True)
con.register("df", df)
df = con.sql("""SELECT df.*, bioproject, organism FROM df
                JOIN (SELECT DISTINCT ON (acc) * FROM sra_metadata) USING (acc)""").df()
df["site"] = df.organism.str.replace("human ", "").str.replace(" metagenome", "").replace("metagenome", "unknown")


def rho(d):
    return spearmanr(d.mbases, d[col]).statistic if len(d) >= 5 and d.mbases.nunique() > 2 else np.nan


proj = (df.groupby("bioproject").agg(site=("site", "first"), runs=("acc", "size"), depth_median=("mbases", "median"),
                                     value_median=(col, "median"))
        .assign(rho=df.groupby("bioproject").apply(rho)).reset_index())
colored = proj.nlargest(len(COLORS), "runs").sort_values("depth_median")
proj.sort_values("depth_median").to_csv(args.out.removesuffix(".png") + "_projects.tsv", sep="\t", index=False,
                                        float_format="%.4g")

big = df[df.bioproject.isin(colored.bioproject)]
ranks = big.groupby("bioproject")[["mbases", col]].rank(pct=True)
r_within = spearmanr(ranks.mbases, ranks[col]).statistic
r_between = spearmanr(proj.depth_median, proj.value_median).statistic

fig, ax = plt.subplots(figsize=(8, 5))
rest = df[~df.bioproject.isin(colored.bioproject)]
if len(rest):
    ax.scatter(rest.mbases, rest[col], s=12, color=OTHER, linewidths=0,
               label=f"other ({rest.bioproject.nunique()} projects, {len(rest)} runs)")
for c, r in zip(COLORS, colored.itertuples()):
    d = df[df.bioproject == r.bioproject]
    lab = f"{r.bioproject} ({r.site}, {r.runs} runs" + (f", rho={r.rho:.2f})" if not np.isnan(r.rho) else ")")
    ax.scatter(d.mbases, d[col], s=12, color=c, alpha=0.8, linewidths=0, label=lab)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(args.xlabel, color=INK)
ax.set_ylabel(args.ylabel or col, color=INK)
ax.grid(color=GRID, lw=0.6)
ax.set_axisbelow(True)
ax.tick_params(colors=MUTED)
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
ax.legend(fontsize=7, frameon=False, labelcolor=INK, loc="upper left", bbox_to_anchor=(1.01, 1))
ax.set_title(f"{args.title}\n{args.group}: overall rho={rho(df):.2f}; between projects {r_between:.2f}; "
             f"within projects {r_within:.2f}", color=INK, fontsize=10)
fig.tight_layout()
fig.savefig(args.out, dpi=200, bbox_inches="tight")
