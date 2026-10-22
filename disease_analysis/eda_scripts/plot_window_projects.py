"""Disease and control samples within one sequencing-depth window, split by BioProject, one value per person.

Left panel: COLUMN per Mbase per person, one column of points per BioProject, labeled with its library
selection and body site (SRA libraryselection, organism), so a study with a different protocol stands out.
Right panel: COLUMN against depth for the same people; --highlight projects are drawn as open markers,
showing whether they differ at equal depth. Runs of one biosample are collapsed to their medians.
Also writes OUT_projects.tsv (bioproject, group, people, median depth, COLUMN, COLUMN per Mbase).

Usage: plot_window_projects.py IN.tsv COLUMN LO HI OUT.png TITLE [--highlight PRJ ...] [--ylabel TEXT] [--xlabel TEXT]
  IN.tsv  from 36_per_mbp.py: acc, group, mbases, COLUMN, COLUMN_per_mbp
  LO, HI  depth window in Mbases (inclusive)
"""
import argparse

import duckdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import common

# same group colors as plot_group_hist.py and plot_depth_curve.py
GROUP_COLOR = dict(zip(common.GROUPS, ["#2a78d6", "#eb6834", "#1baf7a"]))
INK, MUTED, GRID = "#2b2b28", "#6b6a64", "#e5e4df"

ap = argparse.ArgumentParser()
for a in ["tsv", "column", "lo", "hi", "out", "title"]:
    ap.add_argument(a)
ap.add_argument("--highlight", nargs="*", default=[])
ap.add_argument("--ylabel")
ap.add_argument("--xlabel", default="sequencing depth (Mbases, before human removal)")
args = ap.parse_args()
col, per = args.column, f"{args.column}_per_mbp"
ylabel = args.ylabel or col

df = pd.read_csv(args.tsv, sep="\t", dtype={"acc": str})
df = df[df.group.isin(["disease", "control"]) & df.mbases.between(float(args.lo), float(args.hi)) & (df[col] > 0)]
con = duckdb.connect(common.METADATA_DB, read_only=True)
con.register("df", df)
df = con.sql("""SELECT df.*, bioproject, biosample, libraryselection, organism FROM df
                JOIN (SELECT DISTINCT ON (acc) * FROM sra_metadata) USING (acc)""").df()
people = (df.groupby(["group", "bioproject", "biosample"])
          .agg(mbases=("mbases", "median"), n=(col, "median"), r=(per, "median"),
               sel=("libraryselection", "first"), org=("organism", "first")).reset_index())

proj = (people.groupby(["bioproject", "group"])
        .agg(people=("n", "size"), mbases=("mbases", "median"), n=("n", "median"), r=("r", "median")).reset_index()
        .rename(columns={"n": col, "r": per}))
proj.to_csv(args.out.removesuffix(".png") + "_projects.tsv", sep="\t", index=False, float_format="%.4g")

# projects ordered by their median per Mbase; tick label = project, protocol, site
order = people.groupby("bioproject").r.median().sort_values().index.tolist()
info = people.groupby("bioproject").agg(sel=("sel", "first"), org=("org", "first"))
ticks = [f"{p}\n{info.sel[p]}\n{info.org[p].replace('human ', '').replace(' metagenome', '') or 'unknown site'}" for p in order]

rng = np.random.default_rng(0)
fig, axes = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw={"width_ratios": [1.3, 1]})
ax = axes[0]
for i, p in enumerate(order):
    for j, grp in enumerate(["disease", "control"]):
        d = people[(people.bioproject == p) & (people.group == grp)]
        if d.empty:
            continue
        x = i + (j - 0.5) * 0.3 + rng.uniform(-0.08, 0.08, len(d))
        ax.scatter(x, d.r, s=14, color=GROUP_COLOR[grp], alpha=0.7, linewidths=0)
        ax.plot([i + (j - 0.5) * 0.3 - 0.12, i + (j - 0.5) * 0.3 + 0.12], [d.r.median()] * 2, color=INK, lw=2)
ax.axhline(people.r.median(), color=MUTED, lw=1.2, ls="--")
ax.set_xticks(range(len(order)), ticks, fontsize=7)
for t, p in zip(ax.get_xticklabels(), order):
    t.set_color(INK if p in args.highlight else MUTED)
    t.set_fontweight("bold" if p in args.highlight else "normal")
ax.set_ylabel(f"{ylabel} / Mbase (per person)", color=INK)
ax.set_title("by BioProject (bars: median; dashed: all people)", color=INK, fontsize=10)

ax = axes[1]
for grp in ["disease", "control"]:
    for hl in [False, True]:
        d = people[(people.group == grp) & (people.bioproject.isin(args.highlight) == hl)]
        if d.empty:
            continue
        ax.scatter(d.mbases, d.n, s=22 if hl else 14, linewidths=1.2 if hl else 0, alpha=0.9 if hl else 0.6,
                   facecolors="none" if hl else GROUP_COLOR[grp], edgecolors=GROUP_COLOR[grp],
                   label=f"{grp}, {', '.join(args.highlight)} (n={len(d)})" if hl else f"{grp}, other projects (n={len(d)})")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(args.xlabel, color=INK)
ax.set_ylabel(f"{ylabel} (per person)", color=INK)
ax.set_title("same people against depth (median of each person's runs)", color=INK, fontsize=10)
ax.legend(fontsize=7, frameon=False, labelcolor=INK)

for ax in axes:
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTED)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
fig.suptitle(f"{args.title} ({float(args.lo):g}-{float(args.hi):g} Mbases)", color=INK)
fig.tight_layout()
fig.savefig(args.out, dpi=200, bbox_inches="tight")
