"""Group differences in a per-sample count at matched sequencing depth, instead of dividing by depth.

Counts grow slower than depth (log-log slope < 1), so count / Mbase favours shallow samples. Two checks:
  bins        samples cut into depth quantile bins (pooled over groups); per bin and group the median
              count, and disease vs control within the bin (median ratio, Mann-Whitney p)
  regression  OLS of log(count) on group (reference: control), three ways: group only; + log(depth);
              + log(depth) + BioProject fixed effects (so only within-project contrasts count;
              NA for a group that shares no BioProject with control).
              exp(coef) is the fold change at equal depth; SEs are clustered by BioProject.
Samples with count <= 0 or mbases <= 0 are left out. With --facet, everything is run per facet value.
Writes OUT_bins.tsv, OUT_regression.tsv and OUT.png (bin medians vs depth, log-log).
Usage: 39_depth_adjusted.py IN.tsv COLUMN OUT TITLE [--facet COL] [--bins N] [--ylabel TEXT] [--xlabel TEXT]
"""
import argparse
import warnings

import duckdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu

import common

ap = argparse.ArgumentParser()
for a in ["tsv", "column", "out", "title"]:
    ap.add_argument(a)
ap.add_argument("--facet")
ap.add_argument("--bins", type=int, default=10)
ap.add_argument("--ylabel")
ap.add_argument("--xlabel", default="sequencing depth (Mbases)")
args = ap.parse_args()
col = args.column
COLORS = {"disease": "#eb6834", "control": "#2a78d6", "other_control": "#1baf7a"}
MIN_N = 5  # fewer samples of a group in a bin: no median or test for it

df = pd.read_csv(args.tsv, sep="\t", dtype={"acc": str, **({args.facet: str} if args.facet else {})})
df = df[(df[col] > 0) & (df.mbases > 0)]
con = duckdb.connect(common.METADATA_DB, read_only=True)
con.register("df", df)
df = con.sql("""SELECT df.*, bioproject FROM df
                JOIN (SELECT DISTINCT ON (acc) acc, bioproject FROM sra_metadata) USING (acc)""").df()
df["bioproject"] = df.bioproject.fillna("unknown")
df["y"], df["x"] = np.log(df[col]), np.log(df.mbases)
df["group"] = pd.Categorical(df.group, categories=["control", "disease", "other_control"])


def bins(d):
    d = d.assign(bin=pd.qcut(d.mbases, args.bins, labels=False, duplicates="drop"))
    rows = []
    for b, e in d.groupby("bin"):
        r = {"bin": b, "depth_lo": e.mbases.min(), "depth_hi": e.mbases.max(), "depth_median": e.mbases.median()}
        x = {}
        for g in common.GROUPS:
            v = e.loc[e.group == g, col]
            r[f"n_{g}"] = len(v)
            r[f"median_{g}"] = v.median() if len(v) >= MIN_N else np.nan
            x[g] = v
        ok = min(len(x["disease"]), len(x["control"])) >= MIN_N
        r["ratio_disease_control"] = r["median_disease"] / r["median_control"] if ok else np.nan
        r["p_mannwhitney"] = mannwhitneyu(x["disease"], x["control"]).pvalue if ok else np.nan
        rows.append(r)
    return pd.DataFrame(rows)


def regression(d):
    both = {g: d.groupby("bioproject").group.agg(lambda x: {g, "control"} <= set(x)).sum()
            for g in ["disease", "other_control"]}
    rows = []
    for model, f in [("group only", "y ~ group"), ("+ depth", "y ~ group + x"),
                     ("+ depth + BioProject", "y ~ group + x + C(bioproject)")]:
        with warnings.catch_warnings():  # rank-deficient when a group shares no project with control
            warnings.simplefilter("ignore")
            fit = smf.ols(f, d).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(d.bioproject)[0]})
        ci = fit.conf_int()
        for g in ["disease", "other_control"]:
            t = f"group[T.{g}]"
            na = "BioProject" in model and both[g] == 0  # no within-project contrast: not identified
            rows.append({"model": model, "contrast": f"{g} vs control", "n": int(fit.nobs),
                         "n_projects": d.bioproject.nunique(), "n_projects_with_both": both[g],
                         **{k: np.nan if na else v for k, v in {
                             "fold": np.exp(fit.params[t]), "fold_lo": np.exp(ci.loc[t, 0]),
                             "fold_hi": np.exp(ci.loc[t, 1]), "p-value": fit.pvalues[t]}.items()},
                         "depth_slope": fit.params.get("x", np.nan)})
    return pd.DataFrame(rows)


panels = list(df.groupby(args.facet, sort=False)) if args.facet else [(None, df)]
binned, reg = [], []
for key, d in panels:
    b, r = bins(d), regression(d)
    if args.facet:
        b.insert(0, args.facet, key)
        r.insert(0, args.facet, key)
    binned.append(b)
    reg.append(r)
binned, reg = pd.concat(binned), pd.concat(reg)
binned.to_csv(f"{args.out}_bins.tsv", sep="\t", index=False, float_format="%.6g")
reg.to_csv(f"{args.out}_regression.tsv", sep="\t", index=False, float_format="%.6g")

ncol = min(len(panels), 4)
nrow = -(-len(panels) // ncol)
fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 3.8 * nrow), squeeze=False)
for ax, (key, _) in zip(axes.flat, panels):
    b = binned[binned[args.facet] == key] if args.facet else binned
    r = reg[reg[args.facet] == key] if args.facet else reg
    for g in common.GROUPS:
        ax.plot(b.depth_median, b[f"median_{g}"], "-o", color=COLORS[g], lw=2, ms=5, label=g)
    ax.set(xscale="log", yscale="log", xlabel=args.xlabel, ylabel=args.ylabel or col)
    dp = r[(r.model == "+ depth") & (r.contrast == "disease vs control")].iloc[0]
    ax.text(0.02, 0.98, f"disease/control at equal depth: {dp.fold:.2f}x (p={dp['p-value']:.2g})",
            transform=ax.transAxes, va="top", fontsize=8, color="0.25")
    ax.grid(alpha=0.3, lw=0.5)
    if args.facet:
        ax.set_title(f"{args.facet} {key}", fontsize=9)
for ax in axes.flat[len(panels):]:
    ax.axis("off")
axes.flat[0].legend(fontsize=8, loc="lower right", frameon=False)
fig.suptitle(f"{args.title}\nmedian per depth bin ({args.bins} quantile bins; groups with < {MIN_N} samples left out)")
fig.tight_layout()
fig.savefig(f"{args.out}.png", dpi=200, bbox_inches="tight")
print(reg[reg.contrast == "disease vs control"].round(4).to_string(index=False))
