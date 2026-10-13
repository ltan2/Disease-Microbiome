"""Metadata EDA of one disease's disease / control / other_control samples (see common.split_groups).

Who the samples are and where they come from: group sizes and data availability (sketch, taxa, KO),
BioProjects (study design: case+control vs case-only), countries / continents, sequencing
(platform, layout, library selection, organism = body site, release year, depth), sample material,
and host covariates (age, sex, BMI) with how often each is filled. Also flags how confounded
disease vs control is with BioProject and country (Cramer's V).

Writes results/<disease>/metadata/:
  samples.tsv        one row per labeled sample: group, data availability, SRA + host metadata
  counts.tsv         per group: samples, data availability, BioProjects, countries, biosamples, subjects
  bioprojects.tsv    per BioProject: n per group, design, countries, centers, instruments, median Mbases
  <field>.tsv        group x value counts for country, continent, platform, instrument, librarylayout,
                     libraryselection, organism, material, release_year, sex
  host.tsv           per group: coverage and median / IQR of age and BMI
  summary.txt        headline numbers and confounding checks
  groups.png geography.png bioprojects.png sequencing.png host.png
Usage: 18_metadata_eda.py DISEASE
"""
import os
import sys

import duckdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

import common

# categorical slots 1-8 (blue, orange, aqua, yellow, magenta, green, violet, red); "Other"/missing in grey
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
GREY = "#b4b2a9"
GROUP_COLOR = dict(zip(common.GROUPS, PALETTE))
# BioSample attribute keys per host field, in order of preference (first filled one wins)
ATTR = {
    "subject": ["host_subject_id_sam", "subject_id_sam", "host_subject_id_sam_s_dpl", "id_sam"],
    "age": ["host_age_sam", "age_sam", "age_at_collection_sam", "host_age_years_sam", "age_years_sam"],
    "sex": ["sex_calc", "host_sex_sam", "sex_sam", "gender_sam", "host_gender_sam"],
    "bmi": ["host_body_mass_index_sam_s_dpl230", "bmi_sam", "host_bmi_sam", "host_body_mass_index_sam"],
    "material": ["env_material_sam", "isolation_source_sam", "sample_type_sam", "body_site_sam",
                 "host_body_site_sam", "tissue_sam", "env_medium_sam"],
}
MISSING = {"", "na", "n/a", "nan", "none", "missing", "not applicable", "not collected", "not provided",
           "unknown", "not available", "uncalculated", "-", "missing: not provided"}

disease = sys.argv[1]
out = common.out_dir(disease, "metadata")
groups = common.split_groups(disease)
lab = pd.DataFrame([(a, g) for g in common.GROUPS for a in sorted(groups[g])], columns=["acc", "group"])

# ---- data availability per sample ----
yacht = set(open(common.YACHT_SAMPLES).read().split())
lab["sketch"] = lab.acc.map(lambda a: os.path.exists(common.nohuman_path(a)))
lab["taxa"] = lab.acc.isin(yacht)
lab["func"] = lab.acc.map(lambda a: os.path.exists(common.ko_path(a)))

# ---- SRA run metadata + selected BioSample attributes ----
con = duckdb.connect(common.METADATA_DB, read_only=True)
con.register("lab", lab[["acc"]])
run = con.sql("""SELECT acc, bioproject, sra_study, biosample, center_name, organism, platform, instrument,
                        librarylayout, libraryselection, mbases, avgspotlen, year(releasedate) AS release_year,
                        geo_loc_name_country_calc AS country, geo_loc_name_country_continent_calc AS continent
                 FROM (SELECT DISTINCT ON (acc) * FROM sra_metadata) JOIN lab USING (acc)""").df()
keys = sorted({k for ks in ATTR.values() for k in ks})
attr = con.sql(f"""SELECT acc, k, any_value(v) AS v FROM attributes JOIN lab USING (acc)
                   WHERE k IN ({', '.join(f"'{k}'" for k in keys)}) GROUP BY acc, k""").df()
attr = attr[~attr.v.str.strip().str.lower().isin(MISSING)].pivot(index="acc", columns="k", values="v")
host = pd.DataFrame(index=attr.index)
for field, ks in ATTR.items():
    ks = [k for k in ks if k in attr]
    host[field] = attr[ks].bfill(axis=1).iloc[:, 0] if ks else np.nan

df = lab.merge(run, on="acc", how="left").merge(host, left_on="acc", right_index=True, how="left")
for c in ["country", "continent"]:
    df[c] = df[c].where(~df[c].astype(str).str.lower().isin(MISSING))
df["age"] = pd.to_numeric(df.age.astype(str).str.extract(r"^\s*(\d+(?:\.\d+)?)", expand=False), errors="coerce")
df.loc[~df.age.between(0, 120), "age"] = np.nan
df["bmi"] = pd.to_numeric(df.bmi, errors="coerce")
df.loc[~df.bmi.between(10, 80), "bmi"] = np.nan
sex = df.sex.astype(str).str.strip().str.lower()
df["sex"] = np.select([sex.isin(["male", "m", "man"]), sex.isin(["female", "f", "woman"])], ["male", "female"], None)
df["material_raw"] = df.material.str.strip().str.lower()
# free-text material -> a few broad classes (first match wins); unmatched values are kept as they are
MATERIAL = [("stool", r"stool|fec|faec|feces|distal gut|caec|cecum"), ("biopsy / tissue", r"biops|mucos|tissue|intestin|colon|ileum|rect"),
            ("cell culture", r"culture"), ("oral", r"oral|saliva|tongue|plaque|mouth"), ("respiratory", r"respir|nasal|sputum|lung")]
df["material"] = np.select([df.material_raw.str.contains(rx, na=False) for _, rx in MATERIAL], [m for m, _ in MATERIAL],
                           None)
df["material"] = df.material.where(df.material.notna(), df.material_raw)
df["subject"] = np.where(df.subject.notna(), df.bioproject + ":" + df.subject.astype(str), None)
df.to_csv(f"{out}/samples.tsv", sep="\t", index=False)

# ---- per-group counts ----
g = df.groupby("group", sort=False)
counts = pd.DataFrame({
    "samples": g.size(), "with_sketch": g.sketch.sum(), "with_taxa": g.taxa.sum(), "with_func": g.func.sum(),
    "bioprojects": g.bioproject.nunique(), "countries": g.country.nunique(), "biosamples": g.biosample.nunique(),
    "with_subject_id": g.subject.count(), "subjects": g.subject.nunique(),
}).reindex(common.GROUPS)
counts.to_csv(f"{out}/counts.tsv", sep="\t")

# ---- BioProjects ----
bp = df.pivot_table(index="bioproject", columns="group", values="acc", aggfunc="size", fill_value=0)
bp = bp.reindex(columns=common.GROUPS, fill_value=0)
bp["design"] = np.select([(bp.disease > 0) & (bp.control > 0), bp.disease > 0, bp.control > 0],
                         ["case+control", "case only", "control only"], "other_control only")
agg = lambda c: df.groupby("bioproject")[c].agg(lambda v: "; ".join(sorted(v.dropna().astype(str).unique())))
bp["countries"], bp["centers"], bp["instruments"] = agg("country"), agg("center_name"), agg("instrument")
bp["median_mbases"] = df.groupby("bioproject").mbases.median()
bp = bp.sort_values(["disease", "control", "other_control"], ascending=False)
bp.to_csv(f"{out}/bioprojects.tsv", sep="\t")


def crosstab(col, top=None):
    """group x value counts (columns = groups); values past `top` (by total) fold into 'Other'."""
    v = df[col].astype("string").fillna("missing")
    if top is not None:
        keep = v[v != "missing"].value_counts().index[:top]
        v = v.where(v.isin(keep) | (v == "missing"), "Other")
    ct = pd.crosstab(v, df.group).reindex(columns=common.GROUPS, fill_value=0)
    return ct.loc[ct.sum(axis=1).sort_values(ascending=False).index]


FIELDS = ["country", "continent", "platform", "instrument", "librarylayout", "libraryselection", "organism",
          "material", "release_year", "sex"]
for f in FIELDS:
    crosstab(f).to_csv(f"{out}/{f}.tsv", sep="\t")


# ---- host covariates ----
def q(v, p):
    return v.quantile(p) if v.notna().any() else np.nan


rows = []
for grp in common.GROUPS:
    d = df[df.group == grp]
    r = {"group": grp, "n": len(d)}
    for f in ["age", "bmi"]:
        r |= {f"{f}_filled": d[f].notna().mean(), f"{f}_median": q(d[f], .5), f"{f}_q1": q(d[f], .25), f"{f}_q3": q(d[f], .75)}
    r |= {"sex_filled": d.sex.notna().mean(), "frac_female": (d.sex == "female").sum() / max(d.sex.notna().sum(), 1)}
    rows.append(r)
hosttab = pd.DataFrame(rows)
hosttab.to_csv(f"{out}/host.tsv", sep="\t", index=False, float_format="%.3g")


# ---- confounding: disease vs control only ----
def cramers_v(col):
    d = df[df.group.isin(["disease", "control"]) & df[col].notna()]
    ct = pd.crosstab(d[col], d.group)
    if ct.shape[0] < 2 or ct.shape[1] < 2:
        return np.nan
    return np.sqrt(chi2_contingency(ct)[0] / ct.values.sum() / (min(ct.shape) - 1))


dis = df[df.group == "disease"]
paired = bp.index[bp.design == "case+control"]
lines = [
    f"Disease: {disease}",
    "",
    counts.to_string(),
    "",
    f"BioProjects with disease or control samples: {(bp.design != 'other_control only').sum()}",
    bp.design.value_counts().to_string(),
    f"disease samples in a case+control BioProject: {dis.bioproject.isin(paired).sum()} / {len(dis)}"
    f" ({dis.bioproject.isin(paired).mean():.1%})",
    f"largest BioProject share of disease samples: {bp.disease.max() / max(len(dis), 1):.1%} ({bp.disease.idxmax()})",
    "",
    "Cramer's V between disease/control label and (0 = independent, 1 = label fully determined):",
    *(f"  {c}: {cramers_v(c):.3f}" for c in ["bioproject", "country", "platform", "librarylayout", "organism"]),
    "",
    "Top countries (disease / control / other_control):",
    crosstab("country", 10).to_string(),
    "",
    "Host covariates:",
    hosttab.to_string(index=False, float_format=lambda x: f"{x:.3g}"),
]
open(f"{out}/summary.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))


# ---- figures ----
def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#e5e4df", lw=0.8)
    ax.set_axisbelow(True)


def stacked_counts(ax, ct, title, groups=common.GROUPS):
    """Horizontal bars of counts per value, stacked by group, largest value on top."""
    ct = ct[groups].iloc[::-1]
    left = np.zeros(len(ct))
    for grp in groups:
        ax.barh(ct.index.astype(str), ct[grp], left=left, color=GROUP_COLOR[grp], label=grp,
                edgecolor="white", linewidth=1)
        left += ct[grp].to_numpy()
    ax.set_title(title, loc="left")
    ax.set_xlabel("samples")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    style(ax)


def share_bars(ax, col, title, top=7):
    """One 100%-stacked bar per group showing the share of each value of `col` (top values + Other + missing)."""
    ct = crosstab(col, top)
    order = [v for v in ct.index if v not in ("Other", "missing")]
    order = (sorted(order) if col == "release_year" else order) + [v for v in ("Other", "missing") if v in ct.index]
    share = (ct.loc[order] / ct.sum()).T.reindex(common.GROUPS[::-1])
    left = np.zeros(len(share))
    colors = PALETTE + [GREY] * 2
    for i, v in enumerate(order):
        c = GREY if v == "Other" else "#dddcd5" if v == "missing" else colors[i]
        ax.barh(share.index, share[v], left=left, color=c, label=v, edgecolor="white", linewidth=1)
        left += share[v].to_numpy()
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1))
    ax.set_title(title, loc="left")
    ax.legend(frameon=False, fontsize=7, bbox_to_anchor=(1, 0.5), loc="center left")
    style(ax)


# groups: labeled samples and how many have each data kind
fig, ax = plt.subplots(figsize=(9, 4))
kinds = ["samples", "with_sketch", "with_taxa", "with_func"]
y = np.arange(len(common.GROUPS))
for i, k in enumerate(kinds):
    ax.barh(y + (i - 1.5) * 0.2, counts[k], height=0.18, color=PALETTE[i], label=k, edgecolor="white")
    for yy, v in zip(y + (i - 1.5) * 0.2, counts[k]):
        ax.text(v, yy, f" {v:,}", va="center", fontsize=7, color="#3d3d3a")
ax.set_yticks(y, [f"{grp}\n{counts.bioprojects[grp]} BioProjects" for grp in common.GROUPS])
ax.invert_yaxis()
ax.set_xlabel("samples")
ax.set_title(f"{disease}: samples per group and data availability", loc="left")
ax.legend(frameon=False, fontsize=8, loc="lower right")
style(ax)
fig.tight_layout()
fig.savefig(f"{out}/groups.png", dpi=200, bbox_inches="tight")

# geography: countries of disease + control samples, continent share for every group
fig, (a1, a2) = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={"width_ratios": [1.3, 1]})
dc = crosstab("country")[["disease", "control"]]
dc = dc[dc.sum(axis=1) > 0].loc[lambda d: d.sum(axis=1).sort_values(ascending=False).index]
dc = pd.concat([dc[:20], dc[20:].sum().to_frame("Other").T]) if len(dc) > 21 else dc
stacked_counts(a1, dc, f"countries with disease/control samples ({len(dc)} shown)", ["disease", "control"])
share_bars(a2, "continent", "continent share per group")
fig.suptitle(f"{disease}: geographic distribution")
fig.tight_layout()
fig.savefig(f"{out}/geography.png", dpi=200, bbox_inches="tight")

# BioProjects: top 30 by disease + control, labeled with country
top_bp = bp[bp.design != "other_control only"].assign(n=lambda b: b.disease + b.control).nlargest(30, "n")
fig, ax = plt.subplots(figsize=(11, max(4, 0.32 * len(top_bp) + 1.5)))
labels = top_bp.index + "  (" + top_bp.countries.str.slice(0, 25).replace("", "?") + ")"
stacked_counts(ax, top_bp.set_index(labels)[["disease", "control"]],
               f"top {len(top_bp)} of {(bp.design != 'other_control only').sum()} BioProjects with disease/control samples",
               ["disease", "control"])
fig.suptitle(f"{disease}: samples per BioProject")
fig.tight_layout()
fig.savefig(f"{out}/bioprojects.png", dpi=200, bbox_inches="tight")

# sequencing / sample type
fig, axes = plt.subplots(3, 2, figsize=(16, 11))
for ax, (col, title) in zip(axes.flat, [("platform", "sequencing platform"), ("librarylayout", "library layout"),
                                         ("libraryselection", "library selection"), ("organism", "SRA organism (body site)"),
                                         ("material", "sample material (BioSample)"), ("release_year", "SRA release year")]):
    share_bars(ax, col, title)
fig.suptitle(f"{disease}: sequencing and sample type, share of each group")
fig.tight_layout()
fig.savefig(f"{out}/sequencing.png", dpi=200, bbox_inches="tight")

# host: age and BMI distributions, sex share
fig, axes = plt.subplots(1, 3, figsize=(17, 4.5))
for ax, f, bins, xl in [(axes[0], "age", np.arange(0, 101, 5), "age (years)"), (axes[1], "bmi", np.arange(10, 61, 2), "BMI")]:
    for grp in common.GROUPS:
        v = df[f][df.group == grp].dropna()
        if len(v):
            ax.hist(v, bins=bins, density=True, histtype="step", lw=2, color=GROUP_COLOR[grp],
                    label=f"{grp} (n={len(v):,} of {(df.group == grp).sum():,}, median {v.median():.0f})")
    ax.set_xlabel(xl)
    ax.set_ylabel("density")
    ax.set_title(f"{f} (only samples that report it)", loc="left")
    ax.legend(frameon=False, fontsize=7)
    style(ax)
share_bars(axes[2], "sex", "sex")
fig.suptitle(f"{disease}: host covariates")
fig.tight_layout()
fig.savefig(f"{out}/host.png", dpi=200, bbox_inches="tight")
