"""LaTeX tables of the per-sample count tests for several diseases.

  all samples     Kruskal-Wallis + Dunn (from <kind>/kruskal_dunn.tsv), with n and median count per group
  within studies  per-BioProject Mann-Whitney + van Elteren (from within_project.tsv)
Writes results/count_tests_all.tex and results/count_tests_within.tex (longtable fragments sized for a
single-column portrait page; need the booktabs and longtable packages and must not be put inside a table float) and results/count_tests.tex, a standalone document with both.
With --depth nonhuman|raw, only the all-samples table is written, on counts per Mbase of sequencing depth
(from <kind>/kruskal_dunn_per_mbp<suffix>.tsv, made by make count-tests-per-mbp), to
results/count_tests_all_per_mbp<suffix>.tex and the standalone results/count_tests_per_mbp<suffix>.tex.
Usage: 22_latex_tables.py [--depth nonhuman|raw] DISEASE [DISEASE ...]
"""
import argparse
import os

import pandas as pd

import common

# kind -> (data label, count file, count column)
KINDS = {"sketch": ("Sketch (unique hashes)", "hash_counts.tsv", "n_unique_hashes"),
         "func": ("KO (distinct KOs)", "distinct_per_sample.tsv", "n_distinct"),
         "taxa": ("Taxa (distinct taxa)", "distinct_per_sample.tsv", "n_distinct")}
# per-Mbase versions: (data label, count file, count column), file and test suffix appended
KINDS_PER_MBP = {"sketch": ("Sketch (unique hashes / Mbase)", "hash_counts_per_mbp", "n_unique_hashes_per_mbp"),
                 "func": ("KO (distinct KOs / Mbase)", "distinct_per_mbp", "n_distinct_per_mbp"),
                 "taxa": ("Taxa (distinct taxa / Mbase)", "distinct_per_mbp", "n_distinct_per_mbp")}
LABEL = {"disease": "D", "control": "C", "other_control": "O"}


def tex(s):
    return str(s).replace("_", r"\_")


def num(x, digits=2):
    """Fixed point with a real minus sign; -0.00 shows as 0.00."""
    if pd.isna(x):
        return "--"
    s = f"{x:,.{digits}f}"
    return f"$-${s[1:]}" if s.startswith("-") and float(s.replace(",", "")) != 0 else s.lstrip("-")


def med(x):
    """Whole counts as integers; per-Mbase rates to 3 significant figures."""
    if pd.isna(x):
        return "--"
    if float(x).is_integer() or abs(x) >= 1000:
        return f"{x:,.0f}"
    return f"{x:.3g}"


def pval(p):
    if pd.isna(p):
        return "--"
    if p >= 0.01:
        return f"{p:.3f}"
    if p >= 0.001:
        return f"{p:.4f}"
    if p == 0:  # below double precision
        return r"$<10^{-300}$"
    m, e = f"{p:.1e}".split("e")
    return rf"${m}\times10^{{{int(e)}}}$"


def longtable(cols, header, rows, caption, label):
    """Footnote-size longtable with tight columns, so it fits a portrait single-column page."""
    ncol = len(cols.replace("@{}", ""))
    lines = [r"\begingroup\footnotesize\setlength{\tabcolsep}{2.5pt}", rf"\begin{{longtable}}{{{cols}}}",
             rf"\caption{{{caption}}}\label{{{label}}}\\", r"\toprule", header + r" \\", r"\midrule", r"\endfirsthead",
             rf"\multicolumn{{{ncol}}}{{@{{}}l}}{{\emph{{(continued)}}}}\\", r"\toprule",
             header + r" \\", r"\midrule", r"\endhead", r"\bottomrule", r"\endfoot"]
    return "\n".join(lines + rows + [r"\end{longtable}", r"\endgroup", ""])


def section(text, ncol, rule=True):
    """Full-width heading row: disease (bold, between rules) or data / coverage block (italic)."""
    if rule:
        return rf"\midrule\multicolumn{{{ncol}}}{{@{{}}l}}{{\textbf{{{tex(text)}}}}} \\ \midrule"
    return rf"\multicolumn{{{ncol}}}{{@{{}}l}}{{\emph{{{tex(text)}}}}} \\"


def block_name(data, level):
    return f"{data}, coverage {level}" if level else data


def all_samples(diseases, suffix=None):
    """suffix None: raw counts; otherwise per-Mbase counts, suffix "" (raw depth) or "_nonhuman"."""
    rows = []
    for disease in diseases:
        rows.append(section(disease, 7))
        if suffix is None:
            kinds = {k: (data, file, col, "kruskal_dunn.tsv") for k, (data, file, col) in KINDS.items()}
        else:
            kinds = {k: (data, f"{file}{suffix}.tsv", col, f"kruskal_dunn_per_mbp{suffix}.tsv")
                     for k, (data, file, col) in KINDS_PER_MBP.items()}
        for kind, (data, file, col, test_file) in kinds.items():
            out = common.out_dir(disease, kind)
            if not os.path.exists(f"{out}/{test_file}"):
                print(f"skipping {disease} {kind}: no {out}/{test_file}")
                continue
            res = pd.read_csv(f"{out}/{test_file}", sep="\t", dtype={"coverage": str})
            counts = pd.read_csv(f"{out}/{file}", sep="\t", dtype={"coverage": str})
            facet = "coverage" if "coverage" in res else None
            for level, r in res.groupby(facet, sort=False) if facet else [("", res)]:
                c = counts[counts.coverage == level] if facet else counts
                m, n = c.groupby("group")[col].median(), c.groupby("group").size()
                rows.append(section(block_name(data, level), 7, rule=False))
                for _, t in r.iterrows():
                    groups = t.groups.split(" vs ")
                    kw = t.test == "Kruskal-Wallis"
                    rows.append(" & ".join([
                        "KW" if kw else "Dunn", ", ".join(LABEL[g] for g in groups) if kw else
                        f"{LABEL[groups[0]]}--{LABEL[groups[1]]}",
                        "/".join(f"{n[g]}" for g in groups), "/".join(med(m[g]) for g in groups),
                        num(t.statistic), pval(t["p-value"]), pval(t["p-adj (Holm)"])]) + r" \\")
                rows.append(r"\addlinespace")
    header = r"Test & Groups & $n$ & Median & $H$ or $z$ & $p$ & $p_{\mathrm{Holm}}$"
    caption = (r"Per-sample counts compared across all samples: unique hashes of the human-free sketch, "
               r"distinct KOs, and distinct YACHT taxa at each minimum coverage. D = disease, C = control, "
               r"O = other control (controls of other diseases). KW = Kruskal--Wallis ($H$) across D, C and O; "
               r"Dunn = pairwise Dunn test ($z$, tie-corrected; $p_{\mathrm{Holm}}$ adjusted over the three "
               r"pairs), $z>0$ meaning the first group ranks higher. $n$ and medians are given in the order "
               r"of the groups.")
    label = "tab:count-tests-all"
    if suffix is not None:
        depth = ("human-free sequencing depth (SRA Mbases $\\times$ the sample's human-free share of sketch hashes)"
                 if suffix else "sequencing depth (SRA Mbases, before human removal)")
        caption = (caption.replace("Per-sample counts compared", "Per-sample counts per Mbase compared")
                   + f" Each count is divided by the sample's {depth}; samples without depth are left out.")
        label = f"tab:count-tests-all-per-mbp{suffix.replace('_', '-')}"
    return longtable("@{}llrrrrr@{}", header, rows, caption, label)


def within(diseases):
    rows = []
    for disease in diseases:
        rows.append(section(disease, 9))
        path = f"{os.path.dirname(common.out_dir(disease, 'func'))}/within_project.tsv"
        res = pd.read_csv(path, sep="\t", dtype={"coverage": str})
        res["coverage"] = res.coverage.fillna("")
        for (data, level), r in res.groupby(["data", "coverage"], sort=False):
            rows.append(section(block_name(data, level), 9, rule=False))
            for _, t in r.iterrows():
                combined = t.test == "van Elteren"
                rows.append((r"\cmidrule(l){1-9}" if combined else "") + " & ".join([
                    r"\emph{combined}" if combined else tex(t.bioproject),
                    f"{t.disease:,}", f"{t.control:,}", f"{t['median disease']:,.0f}",
                    f"{t['median control']:,.0f}", num(t.statistic, 2 if combined else 1),
                    num(t["cliffs delta"]), pval(t["p-value"]), pval(t["p-adj (Holm)"])]) + r" \\")
            rows.append(r"\addlinespace")
    header = (r"BioProject & $n_\mathrm{D}$ & $n_\mathrm{C}$ & Med.\ D & Med.\ C & "
              r"$U$ or $z$ & $\delta$ & $p$ & $p_{\mathrm{Holm}}$")
    caption = (r"Disease (D) vs control (C) within BioProjects with $\geq$10 samples of each. Per BioProject: "
               r"two-sided Mann--Whitney $U$, $p_{\mathrm{Holm}}$ adjusted over the BioProjects of the block. "
               r"\emph{combined}: van Elteren stratified rank-sum test ($z$) over those BioProjects, $z>0$ "
               r"meaning disease higher; its medians pool all BioProjects. $\delta$ = Cliff's delta, "
               r"$\delta>0$ meaning disease higher.")
    return longtable("@{}lrrrrrrrr@{}", header, rows, caption, "tab:count-tests-within")


ap = argparse.ArgumentParser()
ap.add_argument("--depth", choices=["nonhuman", "raw"])
ap.add_argument("diseases", nargs="+")
args = ap.parse_args()
diseases = args.diseases
if args.depth:
    suffix = "_nonhuman" if args.depth == "nonhuman" else ""
    name = f"count_tests_all_per_mbp{suffix}"
    open(f"{common.RESULTS}/{name}.tex", "w").write(all_samples(diseases, suffix))
    open(f"{common.RESULTS}/count_tests_per_mbp{suffix}.tex", "w").write("\n".join([
        r"\documentclass[10pt]{article}", r"\usepackage{booktabs,longtable}", r"\begin{document}",
        rf"\input{{{name}}}", r"\end{document}", ""]))
    print(f"wrote {common.RESULTS}/{{{name},count_tests_per_mbp{suffix}}}.tex")
    raise SystemExit
open(f"{common.RESULTS}/count_tests_all.tex", "w").write(all_samples(diseases))
open(f"{common.RESULTS}/count_tests_within.tex", "w").write(within(diseases))
open(f"{common.RESULTS}/count_tests.tex", "w").write("\n".join([
    r"\documentclass[10pt]{article}",
    r"\usepackage{booktabs,longtable}",
    r"\begin{document}",
    r"\input{count_tests_all}", r"\clearpage", r"\input{count_tests_within}", r"\end{document}", ""]))
print(f"wrote {common.RESULTS}/count_tests{{,_all,_within}}.tex")
