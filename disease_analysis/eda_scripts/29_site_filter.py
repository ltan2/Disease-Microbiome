"""Restrict one disease's taxa samples to one body site, so the taxa and classifier scripts can be rerun on it.

Body site per sample from the SRA organism (common.body_site): gut (human gut / gut / feces / human feces
metagenome), oral (human oral / human saliva / oral metagenome), unknown (human metagenome, metagenome or
no organism; these do not say where the sample came from) or other (any other site or host). Samples of
SITE are kept from results/<disease>/taxa/{groups.tsv,abundance.tsv.gz}, every group (disease, control,
other_control) alike.

Writes results/<disease>/site_<SITE>/taxa/:
  groups.tsv, abundance.tsv.gz  the unfiltered files restricted to SITE (same columns), read by
                                25_taxa_assoc.py, 26_study_heatmap.py, 27_classify.py and 28_rank_taxa.py
                                when run with SITE=<SITE> in the environment
  site_composition.tsv          samples per group and body site before filtering, plus the organisms
                                behind each site (the before/after table)
Usage: 29_site_filter.py DISEASE SITE   (SITE: gut or oral)
"""
import sys

import duckdb
import pandas as pd

import common

disease, site = sys.argv[1], sys.argv[2]
if site not in common.SITE_ORGANISMS:
    raise SystemExit(f"SITE must be one of {list(common.SITE_ORGANISMS)}")
base, out = common.out_dir(disease, "taxa", site=""), common.out_dir(disease, "taxa", site=site)

groups = pd.read_csv(f"{base}/groups.tsv", sep="\t", dtype={"acc": str})
organism = duckdb.connect(common.METADATA_DB, read_only=True).sql(
    "SELECT DISTINCT ON (acc) acc, organism FROM sra_metadata").df().set_index("acc").organism
groups["organism"] = organism.reindex(groups.acc).to_numpy()
groups["site"] = groups.organism.map(common.body_site)

ab = pd.read_csv(f"{base}/abundance.tsv.gz", sep="\t", dtype={"acc": str})
with_taxa = groups.acc.isin(set(ab.acc))
comp = (groups[with_taxa].groupby(["site", "organism"], dropna=False).group.value_counts().unstack(fill_value=0)
        .reindex(columns=common.GROUPS, fill_value=0))
comp.to_csv(f"{out}/site_composition.tsv", sep="\t")
print(f"{disease}: samples with taxa by body site and organism, before filtering")
print("  " + comp.to_string().replace("\n", "\n  "))
by_site = comp.groupby(level="site").sum()
print("  by site:\n    " + by_site.to_string().replace("\n", "\n    "))

keep = groups[groups.site == site]
keep[["acc", "group"]].to_csv(f"{out}/groups.tsv", sep="\t", index=False)
ab[ab.acc.isin(set(keep.acc))].to_csv(f"{out}/abundance.tsv.gz", sep="\t", index=False)
after = keep[keep.acc.isin(set(ab.acc))].group.value_counts().reindex(common.GROUPS, fill_value=0)
print(f"  kept {site}: {after.to_dict()} -> {out}")
