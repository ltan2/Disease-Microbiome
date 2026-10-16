"""Per-study disease vs control differential abundance of the literature-marker species, as a heatmap
(after Wirbel et al. 2019, Nat Med, Fig. 1b).

Studies are BioProjects (SRA metadata). Samples: disease and control (common.split_groups) with
YACHT taxa; samples with no taxa at all are left out (relative abundance undefined), as are
samples without a BioProject. A study is kept when it has >= MIN_PER_GROUP samples in each group.

Per study and GTDB r232 species (sum of its genomes' abundances, as in 25_taxa_assoc.py):
  effect  log2 fold change of geometric means, mean(log2(a + pc)) disease - control, where a = 0
          for samples where the species was not detected and pc = half the smallest nonzero
          abundance in results/<disease>/taxa/abundance.tsv.gz (same pc as 25_taxa_assoc.py)
  p       two-sided Wilcoxon rank-sum (Mann-Whitney U) on a; zeros are tied and get mid-ranks
          (the test is rank-based, so the log and pc do not change it)
  q       Benjamini-Hochberg over all tested species of that study
A species is tested in a study when present in >= MIN_PREV of the study's disease or control
samples and in >= MIN_DETECT of its samples overall (so 1 sample of a small study is not enough);
otherwise its cell is untestable (NaN, gray), never a zero effect.

Colors are capped at the 98th percentile of |effect| (symmetric).
Heatmap columns: the species of the literature_markers.tsv clades for DISEASE (GTDB suffixed
splits included, e.g. g__Fusobacterium_A), tested in at least one kept study, grouped by marker.
Writes results/<disease>/taxa/study_da_<COLUMN>.tsv (every tested species x study),
study_heatmap_<COLUMN>.png and study_report_<COLUMN>.txt (samples/studies kept and excluded).
Usage: 26_study_heatmap.py DISEASE COLUMN [MIN_PER_GROUP=5] [MIN_PREV=0.1] [MIN_DETECT=3]
"""
import sys

import duckdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import common

ALPHA = 0.05
INK, MUTED, MISSING = "#2b2b28", "#6b6a64", "#c8c7c1"

disease, column = sys.argv[1], sys.argv[2]
min_n = int(sys.argv[3]) if len(sys.argv) > 3 else 5
min_prev = float(sys.argv[4]) if len(sys.argv) > 4 else 0.1
min_detect = int(sys.argv[5]) if len(sys.argv) > 5 else 3
out = common.out_dir(disease, "taxa")
report = []
log = lambda s: (report.append(s), print(s))

# ---- samples and studies ----
ab = pd.read_csv(f"{out}/abundance.tsv.gz", sep="\t", usecols=["acc", "group", "taxon", column], dtype={"acc": str})
pc = ab[column][ab[column] > 0].min() / 2
groups = pd.read_csv(f"{out}/groups.tsv", sep="\t", dtype={"acc": str})
groups = groups[groups.group.isin(["disease", "control"])]
bioproject = duckdb.connect(common.METADATA_DB, read_only=True).sql(
    "SELECT DISTINCT ON (acc) acc, bioproject FROM sra_metadata").df().set_index("acc").bioproject
groups["study"] = bioproject.reindex(groups.acc).to_numpy()

log(f"{disease}: disease vs control by BioProject ({column}, pseudocount {pc:.3g})")
log(f"  labeled samples with a YACHT result: {groups.group.value_counts().to_dict()}")
no_taxa = ~groups.acc.isin(set(ab.acc))
log(f"  excluded, no taxa detected: {groups[no_taxa].group.value_counts().to_dict()}")
no_study = ~no_taxa & groups.study.isna()
log(f"  excluded, no BioProject: {groups[no_study].group.value_counts().to_dict()}")
groups = groups[~no_taxa & ~no_study]
counts = pd.crosstab(groups.study, groups.group).reindex(columns=["disease", "control"], fill_value=0)
keep = counts[(counts.disease >= min_n) & (counts.control >= min_n)].sort_values("disease", ascending=False)
log(f"  studies kept (>= {min_n} disease and >= {min_n} control): {len(keep)}")
for s, r in keep.iterrows():
    log(f"    {s}: {r.disease} disease, {r.control} control")
dropped = counts.drop(keep.index)
log(f"  studies excluded: {len(dropped)} ({dropped.disease.sum()} disease, {dropped.control.sum()} control samples)")
for s, r in dropped.sort_values("disease", ascending=False).iterrows():
    why = "no control" if r.control == 0 else "no disease" if r.disease == 0 else f"< {min_n} in a group"
    log(f"    {s}: {r.disease} disease, {r.control} control ({why})")
if keep.empty:
    raise SystemExit("no study passes; lower MIN_PER_GROUP")

# ---- per-study tests on every species ----
tax = common.read_gtdb().set_index("taxon")
ab = ab[ab.acc.isin(groups.acc)]
ab["species"] = tax.species.reindex(ab.taxon).to_numpy()
rows = []
for study in keep.index:
    g = groups[groups.study == study].set_index("acc").group
    m = (ab[ab.acc.isin(g.index)].pivot_table(index="acc", columns="species", values=column, aggfunc="sum")
         .reindex(g.index).fillna(0))  # samples x species, 0 = not detected
    dis, ctl = m[g == "disease"], m[g == "control"]
    prev = pd.DataFrame({"prev_disease": (dis > 0).mean(), "prev_control": (ctl > 0).mean()})
    tested = (prev.max(axis=1) >= min_prev) & ((m > 0).sum() >= min_detect)
    dis, ctl, prev = dis.loc[:, tested], ctl.loc[:, tested], prev[tested]
    eff = np.log2(dis + pc).mean() - np.log2(ctl + pc).mean()
    p = stats.mannwhitneyu(dis, ctl, alternative="two-sided", axis=0).pvalue
    r = prev.assign(effect=eff, p=p, q=stats.false_discovery_control(p), n_disease=len(dis), n_control=len(ctl))
    rows.append(r.reset_index().assign(study=study))
    log(f"  {study}: {tested.sum()} species tested (prevalence >= {min_prev} in a group, "
        f"detected in >= {min_detect} samples), "
        f"{(r.q < ALPHA).sum()} with q < {ALPHA}")
res = pd.concat(rows, ignore_index=True)
res = res[["study", "species", "n_disease", "n_control", "effect", "p", "q", "prev_disease", "prev_control"]]
res.to_csv(f"{out}/study_da_{column}.tsv", sep="\t", index=False)

# ---- marker species ----
lit = pd.read_csv("literature_markers.tsv", sep="\t")
lit = lit[lit.disease == disease]
cols, marker_of = [], {}
arrow = dict(zip(lit.clade, lit.expected.map({"up": "↑", "down": "↓"}).fillna("↕")))
for clade in lit.clade:
    rank = common.RANK_PREFIX[clade[:3]]
    in_clade = (tax[rank] == clade) | tax[rank].str.startswith(clade + "_")
    sp = sorted(set(tax.species[in_clade]) & set(res.species))
    marker_of.update({s: clade for s in sp})
    # within a marker, order by effect in the largest study (untestable there last)
    eff = res[res.study == keep.index[0]].set_index("species").effect.reindex(sp)
    cols += list(eff.sort_values(ascending=False, na_position="last").index)
log(f"  marker species tested in >= 1 kept study: {len(cols)} "
    f"({', '.join(f'{c}: {sum(v == c for v in marker_of.values())}' for c in lit.clade)})")
if not cols:
    raise SystemExit("no marker species tested in any kept study")
eff = res.pivot(index="study", columns="species", values="effect").reindex(index=keep.index, columns=cols)
q = res.pivot(index="study", columns="species", values="q").reindex(index=keep.index, columns=cols)

# symmetric limits at the 98th percentile of |effect| so a single extreme cell of a small study
# does not wash out the rest; beyond it the colors saturate (arrows on the colorbar)
vmax = np.nanpercentile(np.abs(eff.to_numpy()), 98)
cmap = plt.get_cmap("RdBu_r").copy()
cmap.set_bad(MISSING)
fig, ax = plt.subplots(figsize=(0.32 * len(cols) + 3.5, 0.45 * len(keep) + 3.2))
im = ax.imshow(np.ma.masked_invalid(eff.to_numpy()), cmap=cmap, norm=mcolors.Normalize(-vmax, vmax), aspect="auto")
for (i, j), v in np.ndenumerate(q.to_numpy()):
    if v < ALPHA:
        ax.text(j, i, "*", ha="center", va="center", color="white" if abs(eff.iat[i, j]) > vmax * 0.6 else INK,
                fontsize=11, fontweight="bold")
ax.set_xticks(range(len(cols)))
ax.set_xticklabels([f"{c.removeprefix('s__')} {arrow[marker_of[c]]}" for c in cols], rotation=90, fontsize=7,
                   color=INK, style="italic")
ax.set_yticks(range(len(keep)))
ax.set_yticklabels([f"{s}  (CRC n={r.disease}, ctrl n={r.control})" if disease == "Colorectal cancer"
                    else f"{s}  (disease n={r.disease}, ctrl n={r.control})" for s, r in keep.iterrows()],
                   fontsize=8, color=INK)
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
# white gaps between cells, marker brackets above the columns
ax.set_xticks(np.arange(len(cols) + 1) - 0.5, minor=True)
ax.set_yticks(np.arange(len(keep) + 1) - 0.5, minor=True)
ax.grid(which="minor", color="white", lw=1.5)
ax.tick_params(which="minor", length=0)
start = 0
for clade in lit.clade:
    n = sum(marker_of.get(c) == clade for c in cols)
    if n and not clade.startswith("s__"):  # species columns are already named
        top = ax.get_xaxis_transform()  # x in columns, y in axes fraction
        ax.plot([start - 0.4, start + n - 0.6], [1.02, 1.02], color=MUTED, lw=1.2, transform=top, clip_on=False)
        ax.text(start + (n - 1) / 2, 1.04, clade.removeprefix("g__").removeprefix("s__").removeprefix("f__"),
                ha="center", va="bottom", fontsize=8, color=INK, style="italic", transform=top)
    if n:
        if start:
            ax.axvline(start - 0.5, color="white", lw=4)
        start += n
cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02, extend="both")
cb.set_label("log2 fold change, disease vs control\n(red: enriched in disease, blue: depleted)", fontsize=8,
             color=INK)
cb.ax.tick_params(labelsize=7, colors=MUTED, length=0)
cb.outline.set_visible(False)
ax.set_title(f"{disease}: per-study differential abundance of literature-marker species ({column})\n"
             f"* q < {ALPHA} (Wilcoxon rank-sum, BH within study); gray: not tested in that study; "
             f"arrows: expected direction (↕ variable)",
             fontsize=9, color=INK, pad=30)
fig.savefig(f"{out}/study_heatmap_{column}.png", dpi=200, bbox_inches="tight")
open(f"{out}/study_report_{column}.txt", "w").write("\n".join(report) + "\n")
