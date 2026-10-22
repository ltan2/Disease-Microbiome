# EDA of one disease vs its controls vs other-disease controls.
#
#   make sample-sketch-similar-viz Colorectal cancer    # sketch Jaccard similarity
#   make taxa-pres-abs Colorectal cancer                # YACHT taxa presence/absence
#   make funprofile-pres-abs Colorectal cancer          # KO (function) presence/absence
#   make funprofile-abundance Colorectal cancer         # KO CLR / Aitchison distance
#   make seq-depth Colorectal cancer                    # sequencing depth per group (taxa and KO samples)
#   make per-mbp Colorectal cancer                      # unique hashes, distinct KOs and distinct taxa per Mbase of sequencing
#   make depth-curve Colorectal cancer                  # those counts and counts per Mbase vs depth: per-Mbase does not remove depth
#   make window-projects Colorectal cancer              # disease vs control unique hashes in one depth window, by BioProject
#   make hash-depth Colorectal cancer                   # unique hashes vs depth and per Mbase vs depth, two separate plots
#   make hash-depth-projects Colorectal cancer          # disease samples' unique hashes per Mbase vs depth, colored by BioProject
#   make depth-curve Colorectal cancer DEPTH=raw        # per-mbp, depth-curve, hash-depth on SRA depth (default: human-free)
#   make metadata-eda Colorectal cancer                 # group sizes, BioProjects, geography, sequencing, host covariates
#   make stratified-permanova Colorectal cancer         # disease vs control PERMANOVA before/after accounting for BioProject
#   make count-tests Colorectal cancer                  # Kruskal-Wallis + Dunn on per-sample counts, then disease vs control within BioProjects
#   make taxa-abundance Colorectal cancer               # YACHT taxa relative abundance (containment- and detection-score-based) + grouped histogram
#   make taxa-forest Colorectal cancer                  # per-taxon abundance association forest plots + literature markers (literature_markers.tsv)
#   make study-heatmap Colorectal cancer                # per-BioProject disease vs control log2FC heatmap of literature-marker species
#   make classify Colorectal cancer                     # disease vs control LASSO + random forest on taxa abundance (stratified CV, studies pooled)
#   make rank-taxa Colorectal cancer                    # one table ranking clades by LASSO, random forest and univariate agreement
#   make rank-taxa Colorectal cancer SITE=gut           # same, on gut samples only (results/<disease>/site_gut/; taxa, ml targets)
#   make site-compare Colorectal cancer                 # rank-taxa for all samples and each of SITES, then before/after comparison
#   make dann Colorectal cancer                         # disease vs control neural net, study adversary off (baseline) and on (DANN)
#   make dann Colorectal cancer ADVERSARY=off           # baseline only (ADVERSARY=on: DANN only); DANN_CV=loso: leave-one-study-out
#   make dann-plot Colorectal cancer                    # baseline vs DANN: ROC, AUROC/balanced accuracy, study leakage, training curves
#   make sample-sketch-similar-viz DISEASE="Crohn's disease" PROCS=128
#   make list-diseases                                  # valid disease labels with case/control counts
#
# Outputs go to results/<disease>/{sketch,taxa,func}/ ; labels, human removal and YACHT parsing are
# disease-independent and run once.

PY       ?= /home/grads/lbt5343/miniconda3/envs/logan/bin/python
# PROCS: parallel processes; PERMS: PERMANOVA permutations; MIN_PREV: keep KOs present in at least
# this fraction of samples before CLR; N_SUB: samples per group in the Aitchison heatmap;
# MIN_PER_GROUP: disease and control samples a BioProject needs for the stratified PERMANOVA
PROCS    ?= 64
PERMS    ?= 999
MIN_PREV ?= 0.1
N_SUB    ?= 200
MIN_PER_GROUP ?= 10
# DETECT_WEIGHTS: detection-score weight per coverage threshold, highest to lowest (one per COVERAGES)
DETECT_WEIGHTS ?= 1,1,1,1,1,1,1
# TAXA_MIN_PREV: test clades present in at least this fraction of some group's samples; TOP_N: species per forest
# plot; MEMBERS: most prevalent species shown under each literature marker
TAXA_MIN_PREV ?= 0.1
TOP_N ?= 30
MEMBERS ?= 3
# STUDY_MIN_N: disease and control samples a BioProject needs for the per-study heatmap; HEATMAP_COL: abundance proxy
STUDY_MIN_N ?= 5
# STUDY_MIN_DETECT: samples of a study a species must be detected in to be tested there
STUDY_MIN_DETECT ?= 3
HEATMAP_COL ?= containment_abundance
# CLASSIFY_RANK: GTDB rank the classifier features are summed to (or genome); REPEATS: repeats of the
# stratified 5-fold CV (disease vs control only, BioProject ignored)
CLASSIFY_RANK ?= species
# SITE: run the taxa/ml targets on one body site (gut, oral; 29_site_filter.py); SITES: compared by site-compare
SITE ?=
SITES ?= gut oral
REPEATS ?= 5
# ADVERSARY: study (BioProject) adversary off = baseline, on = DANN, or both; DANN_CV: stratified or loso
# (leave-one-study-out); DANN_LAMBDA: maximum gradient reversal weight
ADVERSARY ?= off on
DANN_CV ?= stratified
DANN_LAMBDA ?= 1

DISEASE_TARGETS := sample-sketch-similar-viz taxa-pres-abs funprofile-pres-abs funprofile-abundance seq-depth per-mbp depth-curve hash-depth hash-depth-projects window-projects metadata-eda stratified-permanova count-tests taxa-abundance taxa-forest study-heatmap classify rank-taxa paper-tables site-compare dann dann-plot

# The disease may follow the target, plain or quoted; its words (and the quoted name, spaces escaped)
# become no-op goals.
space   := $(subst ,, )
ifneq ($(filter $(firstword $(MAKECMDGOALS)),$(DISEASE_TARGETS)),)
  DISEASE ?= $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(DISEASE) $(subst $(space),\ ,$(DISEASE)):;@:)
  ifeq ($(strip $(DISEASE)),)
    $(error give a disease, e.g. make $(firstword $(MAKECMDGOALS)) Colorectal cancer)
  endif
endif
export DISEASE SITE

BASE    := results/$(subst ',,$(subst $(space),_,$(strip $(DISEASE))))
OUT     := $(BASE)$(if $(SITE),/site_$(SITE))
SKETCH  := $(OUT)/sketch
TAXA    := $(OUT)/taxa
FUNC    := $(OUT)/func
META    := $(OUT)/metadata
ML      := $(OUT)/ml
LABELS  := results/disease_labels.tsv
COVERAGES := 1.0 0.5 0.25 0.125 0.0625 0.03125 0.015625
TAXA_JACCARD   := $(foreach c,$(COVERAGES),$(TAXA)/jaccard_$(c).npy $(TAXA)/groups_$(c).tsv)
TAXA_PERMANOVA := $(foreach c,$(COVERAGES),$(TAXA)/permanova_$(c).tsv)
NOHUMAN := results/nohuman_unitigs.done
YACHT   := results/yacht_presence.tsv.gz
CONTAINMENT := results/yacht_containment.tsv.gz
GTDB    := results/gtdb_r232/bac120_taxonomy_r232.tsv.gz results/gtdb_r232/ar53_taxonomy_r232.tsv.gz
ABUND_COLS := containment_abundance detection_abundance
# DEPTH=nonhuman (default): mbases x the sample's human-free share of sketch hashes (37_nonhuman_depth.py);
# DEPTH=raw: SRA mbases (before human removal). nonhuman outputs carry a _nonhuman suffix
DEPTH ?= nonhuman
D       := $(if $(filter nonhuman,$(DEPTH)),_nonhuman)
DLABEL  := $(if $(D),human-free sequencing depth (Mbases, estimated),sequencing depth (Mbases, before human removal))

.PHONY: $(DISEASE_TARGETS) labels remove-human yacht-presence yacht-containment list-diseases

sample-sketch-similar-viz: $(SKETCH)/jaccard_heatmap.png $(SKETCH)/jaccard_summary.tsv $(SKETCH)/permanova.tsv $(SKETCH)/hash_counts.png
taxa-pres-abs: $(TAXA)/prevalence.tsv $(TAXA)/distinct_per_sample.png $(TAXA_PERMANOVA)
funprofile-pres-abs: $(FUNC)/prevalence.tsv $(FUNC)/distinct_per_sample.png $(FUNC)/kruskal_dunn.tsv
funprofile-abundance: $(FUNC)/aitchison_heatmap.png $(FUNC)/aitchison_summary.tsv $(FUNC)/aitchison_permanova.tsv
seq-depth: $(TAXA)/depth.png $(FUNC)/depth.png
per-mbp: $(SKETCH)/hash_counts_per_mbp$(D).png $(TAXA)/distinct_per_mbp$(D).png $(FUNC)/distinct_per_mbp$(D).png
hash-depth-projects: $(SKETCH)/hash_per_mbp_vs_depth_projects$(D).png
hash-depth: $(SKETCH)/hash_vs_depth$(D).png $(SKETCH)/hash_per_mbp_vs_depth$(D).png
window-projects: $(SKETCH)/hash_counts_window_projects$(D).png
depth-curve: $(SKETCH)/hash_counts_depth_curve$(D).png $(TAXA)/distinct_depth_curve$(D).png $(FUNC)/distinct_depth_curve$(D).png
metadata-eda: $(META)/summary.txt
stratified-permanova: $(OUT)/stratified_permanova.tsv
count-tests: $(SKETCH)/kruskal_dunn.tsv $(TAXA)/kruskal_dunn.tsv $(FUNC)/kruskal_dunn.tsv $(OUT)/within_project.tsv
taxa-abundance: $(TAXA)/abundance.tsv.gz $(TAXA)/abundance_hist.png
study-heatmap: $(TAXA)/study_heatmap_$(HEATMAP_COL).png
taxa-forest: $(foreach c,$(ABUND_COLS),$(TAXA)/forest_top_$(c).png $(TAXA)/forest_markers_$(c).png)
classify: $(foreach c,$(ABUND_COLS),$(ML)/auroc_$(c)_$(CLASSIFY_RANK).tsv)
rank-taxa: $(foreach c,$(ABUND_COLS),$(ML)/ranked_$(c)_$(CLASSIFY_RANK).tsv)
paper-tables: $(foreach c,$(ABUND_COLS),$(ML)/paper/tables_$(c)_$(CLASSIFY_RANK).md)
site-compare: $(foreach c,$(ABUND_COLS),$(BASE)/site_compare/compare_$(c)_$(CLASSIFY_RANK).png)
dann-plot: $(foreach c,$(ABUND_COLS),$(ML)/dann/compare_$(c)_$(CLASSIFY_RANK)_$(DANN_CV).png)
dann: $(foreach c,$(ABUND_COLS),$(foreach a,$(ADVERSARY),$(ML)/dann/auroc_$(c)_$(CLASSIFY_RANK)_$(DANN_CV)_$(if $(filter on,$(a)),dann,baseline).tsv))

# ---- shared, run once ----
labels: $(LABELS)
$(LABELS): 00_labels.sh
	bash $< $@

list-diseases: $(LABELS)
	@$(PY) -c "import pandas as pd; print(pd.read_csv('$<', sep='\t').pivot_table(index='disease', columns='status', aggfunc='size', fill_value=0).to_string())"

# drop human hashes from every sketch (re-run to pick up new sketches: make remove-human)
$(NOHUMAN):
	$(PY) 01_remove_human.py $(PROCS) && touch $@
remove-human:
	$(PY) 01_remove_human.py $(PROCS) && touch $(NOHUMAN)

# parse every YACHT result once (re-run to pick up new results: make yacht-presence)
$(YACHT):
	$(PY) 08_yacht_presence.py $(PROCS)
yacht-presence:
	$(PY) 08_yacht_presence.py $(PROCS)

# containment of every detected taxon, parsed from the YACHT results once (re-run: make yacht-containment)
$(CONTAINMENT):
	$(PY) 23_yacht_containment.py $(PROCS)
yacht-containment:
	$(PY) 23_yacht_containment.py $(PROCS)

# ---- sketch similarity, per disease ----
$(SKETCH)/groups.tsv: 02_split_groups.py $(LABELS) $(NOHUMAN)
	$(PY) $< "$$DISEASE" sketch
$(SKETCH)/jaccard.npy: 03_jaccard.py $(SKETCH)/groups.tsv
	$(PY) $< "$$DISEASE" $(PROCS)
$(SKETCH)/jaccard_heatmap.png: 04_heatmap.py $(SKETCH)/jaccard.npy $(SKETCH)/groups.tsv
	$(PY) $^ $@ "$$DISEASE: Jaccard similarity of human-free unitig sketches" "Jaccard similarity"
$(SKETCH)/jaccard_summary.tsv: 05_summary.py $(SKETCH)/jaccard.npy $(SKETCH)/groups.tsv
	$(PY) $^ $@
$(SKETCH)/permanova.tsv $(SKETCH)/pcoa.png &: 06_pcoa_permanova.py $(SKETCH)/jaccard.npy $(SKETCH)/groups.tsv
	$(PY) $^ $(SKETCH)/permanova.tsv $(SKETCH)/pcoa.png "$$DISEASE: Jaccard distance of human-free unitig sketches" $(PERMS) --similarity
$(SKETCH)/hash_counts.tsv: 07_hash_counts.py $(SKETCH)/groups.tsv
	$(PY) $< "$$DISEASE" $(PROCS)
$(SKETCH)/hash_counts.png: plot_group_density.py $(SKETCH)/hash_counts.tsv
	$(PY) $^ n_unique_hashes $@ "$$DISEASE: unique hashes per sample (human removed)" --log --xlabel "# unique hashes"

# ---- presence/absence, per disease: taxa (YACHT) and func (KO profiles) ----
GROUP_NAMES  := disease control other_control
TAXA_PRES    := $(foreach g,$(GROUP_NAMES),$(TAXA)/presence_$(g).tsv.gz)
FUNC_PRES    := $(foreach g,$(GROUP_NAMES),$(FUNC)/presence_$(g).tsv.gz $(FUNC)/matrix_$(g).tsv.gz)

ifeq ($(SITE),)
$(TAXA)/groups.tsv: 02_split_groups.py $(LABELS) $(YACHT)
	$(PY) $< "$$DISEASE" taxa
else
# one body site: the unfiltered groups and abundance (make taxa-abundance without SITE first), restricted
$(TAXA)/groups.tsv $(TAXA)/abundance.tsv.gz &: 29_site_filter.py $(BASE)/taxa/groups.tsv $(BASE)/taxa/abundance.tsv.gz
	$(PY) $< "$$DISEASE" $(SITE)
endif
$(TAXA_PRES) &: 09_taxa_presence.py $(TAXA)/groups.tsv $(YACHT)
	$(PY) $< "$$DISEASE"

$(FUNC)/groups.tsv: 02_split_groups.py $(LABELS)
	$(PY) $< "$$DISEASE" func
$(FUNC_PRES) &: 12_ko_presence.py $(FUNC)/groups.tsv
	$(PY) $< "$$DISEASE" $(PROCS)

# taxa Jaccard (one matrix per coverage), then PERMANOVA/PERMDISP + PCoA per coverage
$(TAXA_JACCARD) &: 17_taxa_jaccard.py $(TAXA_PRES)
	$(PY) $< "$$DISEASE"
$(TAXA)/permanova_%.tsv $(TAXA)/pcoa_%.png: 06_pcoa_permanova.py $(TAXA)/jaccard_%.npy $(TAXA)/groups_%.tsv
	$(PY) $^ $(TAXA)/permanova_$*.tsv $(TAXA)/pcoa_$*.png "$$DISEASE: Jaccard distance of YACHT taxa (coverage $*)" $(PERMS) --similarity

# shared steps; $* is taxa or func
$(OUT)/%/prevalence.tsv: 10_prevalence.py $(OUT)/%/presence_disease.tsv.gz
	$(PY) $< "$$DISEASE" $*
$(OUT)/%/distinct_per_sample.tsv: 11_distinct_per_sample.py $(OUT)/%/presence_disease.tsv.gz
	$(PY) $< "$$DISEASE" $*
$(TAXA)/distinct_per_sample.png: plot_group_density.py $(TAXA)/distinct_per_sample.tsv
	$(PY) $^ n_distinct $@ "$$DISEASE: distinct YACHT taxa per sample" --facet coverage --xlabel "# distinct taxa" --log
$(FUNC)/distinct_per_sample.png: plot_group_density.py $(FUNC)/distinct_per_sample.tsv
	$(PY) $^ n_distinct $@ "$$DISEASE: distinct KOs per sample" --xlabel "# distinct KOs" --log
# Kruskal-Wallis + pairwise Dunn (Holm) on distinct KOs per sample
$(FUNC)/kruskal_dunn.tsv: 18_kruskal_dunn.py $(FUNC)/distinct_per_sample.tsv
	$(PY) $^ n_distinct $@
# same on distinct taxa per sample, per coverage, and on unique hashes per sample
$(TAXA)/kruskal_dunn.tsv: 18_kruskal_dunn.py $(TAXA)/distinct_per_sample.tsv
	$(PY) $^ n_distinct $@ --facet coverage
$(SKETCH)/kruskal_dunn.tsv: 18_kruskal_dunn.py $(SKETCH)/hash_counts.tsv
	$(PY) $^ n_unique_hashes $@

# ---- KO abundance: CLR + Aitchison distance, per disease (reuses func/groups.tsv) ----
$(FUNC)/abundance.tsv.gz: 13_ko_abundance.py $(FUNC)/groups.tsv
	$(PY) $< "$$DISEASE" $(PROCS)
$(FUNC)/clr.tsv.gz $(FUNC)/clr_groups.tsv &: 14_clr.py $(FUNC)/abundance.tsv.gz
	$(PY) $< "$$DISEASE" $(MIN_PREV)
$(FUNC)/aitchison.npy: 15_aitchison.py $(FUNC)/clr.tsv.gz
	$(PY) $< "$$DISEASE"
$(FUNC)/aitchison_heatmap.png: 04_heatmap.py $(FUNC)/aitchison.npy $(FUNC)/clr_groups.tsv
	$(PY) $^ $@ "$$DISEASE: Aitchison distance of KO profiles (up to $(N_SUB) per group)" "Aitchison distance" --n-sub $(N_SUB)
$(FUNC)/aitchison_summary.tsv: 05_summary.py $(FUNC)/aitchison.npy $(FUNC)/clr_groups.tsv
	$(PY) $^ $@
$(FUNC)/aitchison_permanova.tsv $(FUNC)/aitchison_pcoa.png &: 06_pcoa_permanova.py $(FUNC)/aitchison.npy $(FUNC)/clr_groups.tsv
	$(PY) $^ $(FUNC)/aitchison_permanova.tsv $(FUNC)/aitchison_pcoa.png "$$DISEASE: Aitchison distance of KO profiles" $(PERMS)

# ---- sequencing depth per group, for the samples of each data kind ($* = taxa or func) ----
.PRECIOUS: $(OUT)/%/depth.tsv
$(OUT)/%/depth.tsv: 16_seq_depth.py $(OUT)/%/groups.tsv
	$(PY) $< "$$DISEASE" $*
$(OUT)/%/depth.png: plot_group_density.py $(OUT)/%/depth.tsv
	$(PY) $^ mbases $@ "$$DISEASE: sequencing depth ($* samples)" --log --xlabel "Mbases"

# ---- per-sample counts per Mbase of sequencing depth (samples with mbases > 0) ----
.PRECIOUS: $(OUT)/%/depth_nonhuman.tsv
$(OUT)/%/depth_nonhuman.tsv: 37_nonhuman_depth.py $(OUT)/%/depth.tsv
	$(PY) $< "$$DISEASE" $* $(PROCS)
$(SKETCH)/hash_counts_per_mbp$(D).tsv: 36_per_mbp.py $(SKETCH)/hash_counts.tsv $(SKETCH)/depth$(D).tsv
	$(PY) $^ n_unique_hashes $@
$(OUT)/%/distinct_per_mbp$(D).tsv: 36_per_mbp.py $(OUT)/%/distinct_per_sample.tsv $(OUT)/%/depth$(D).tsv
	$(PY) $^ n_distinct $@
$(SKETCH)/hash_counts_per_mbp$(D).png: plot_group_density.py $(SKETCH)/hash_counts_per_mbp$(D).tsv
	$(PY) $^ n_unique_hashes_per_mbp $@ "$$DISEASE: unique hashes per Mbase (human removed)" --log --pseudo 0 --xlabel "unique hashes / Mbase"
$(TAXA)/distinct_per_mbp$(D).png: plot_group_density.py $(TAXA)/distinct_per_mbp$(D).tsv
	$(PY) $^ n_distinct_per_mbp $@ "$$DISEASE: distinct YACHT taxa per Mbase (samples with >= 1 taxon)" --facet coverage --log --pseudo 0 --xlabel "distinct taxa / Mbase"
$(FUNC)/distinct_per_mbp$(D).png: plot_group_density.py $(FUNC)/distinct_per_mbp$(D).tsv
	$(PY) $^ n_distinct_per_mbp $@ "$$DISEASE: distinct KOs per Mbase" --log --pseudo 0 --xlabel "distinct KOs / Mbase"
# disease and control people within one depth window (WINDOW_LO-WINDOW_HI Mbases), by BioProject;
# WINDOW_HIGHLIGHT projects are marked (default: the PCR-selected PRJEB42155)
WINDOW_LO ?= 30
WINDOW_HI ?= 3000
WINDOW_HIGHLIGHT ?= PRJEB42155
$(SKETCH)/hash_counts_window_projects$(D).png $(SKETCH)/hash_counts_window_projects$(D)_projects.tsv &: plot_window_projects.py $(SKETCH)/hash_counts_per_mbp$(D).tsv
	$(PY) $< $(SKETCH)/hash_counts_per_mbp$(D).tsv n_unique_hashes $(WINDOW_LO) $(WINDOW_HI) $(SKETCH)/hash_counts_window_projects$(D).png "$$DISEASE: unique hashes by BioProject" --highlight $(WINDOW_HIGHLIGHT) --ylabel "unique hashes" --xlabel "$(DLABEL)"
# the same counts against depth: raw counts saturate, counts per Mbase fall (taxa at DEPTH_COVERAGE)
DEPTH_COVERAGE ?= 0.0625
$(SKETCH)/hash_counts_depth_curve$(D).png $(SKETCH)/hash_counts_depth_curve$(D)_spearman.tsv &: plot_depth_curve.py $(SKETCH)/hash_counts_per_mbp$(D).tsv
	$(PY) $^ n_unique_hashes $(SKETCH)/hash_counts_depth_curve$(D).png "$$DISEASE: unique hashes vs sequencing depth (human removed)" --ylabel "unique hashes" --xlabel "$(DLABEL)"
$(SKETCH)/hash_vs_depth$(D).png $(SKETCH)/hash_vs_depth$(D)_spearman.tsv &: plot_depth_curve.py $(SKETCH)/hash_counts_per_mbp$(D).tsv
	$(PY) $^ n_unique_hashes $(SKETCH)/hash_vs_depth$(D).png "$$DISEASE: unique hashes vs sequencing depth (human removed)" --ylabel "unique hashes" --xlabel "$(DLABEL)" --panels raw
$(SKETCH)/hash_per_mbp_vs_depth$(D).png $(SKETCH)/hash_per_mbp_vs_depth$(D)_spearman.tsv &: plot_depth_curve.py $(SKETCH)/hash_counts_per_mbp$(D).tsv
	$(PY) $^ n_unique_hashes $(SKETCH)/hash_per_mbp_vs_depth$(D).png "$$DISEASE: unique hashes per Mbase vs sequencing depth (human removed)" --ylabel "unique hashes" --xlabel "$(DLABEL)" --panels per_mbase
$(SKETCH)/hash_per_mbp_vs_depth_projects$(D).png $(SKETCH)/hash_per_mbp_vs_depth_projects$(D)_projects.tsv &: plot_depth_projects.py $(SKETCH)/hash_counts_per_mbp$(D).tsv
	$(PY) $^ n_unique_hashes_per_mbp disease $(SKETCH)/hash_per_mbp_vs_depth_projects$(D).png "$$DISEASE: unique hashes per Mbase vs sequencing depth, by BioProject" --ylabel "unique hashes / Mbase" --xlabel "$(DLABEL)"
$(TAXA)/distinct_depth_curve$(D).png $(TAXA)/distinct_depth_curve$(D)_spearman.tsv &: plot_depth_curve.py $(TAXA)/distinct_per_mbp$(D).tsv
	$(PY) $^ n_distinct $(TAXA)/distinct_depth_curve$(D).png "$$DISEASE: distinct YACHT taxa vs sequencing depth (coverage $(DEPTH_COVERAGE))" --query "coverage == $(DEPTH_COVERAGE)" --ylabel "distinct taxa" --xlabel "$(DLABEL)"
$(FUNC)/distinct_depth_curve$(D).png $(FUNC)/distinct_depth_curve$(D)_spearman.tsv &: plot_depth_curve.py $(FUNC)/distinct_per_mbp$(D).tsv
	$(PY) $^ n_distinct $(FUNC)/distinct_depth_curve$(D).png "$$DISEASE: distinct KOs vs sequencing depth" --ylabel "distinct KOs" --xlabel "$(DLABEL)"

# ---- metadata EDA, per disease: all labeled samples, whatever data they have ----
$(META)/summary.txt: 18_metadata_eda.py $(LABELS)
	$(PY) $< "$$DISEASE"

# ---- disease vs control PERMANOVA before/after accounting for BioProject (samples that finished every step) ----
$(OUT)/stratified_permanova.tsv $(OUT)/stratified_projects.tsv &: 20_stratified_permanova.py $(SKETCH)/jaccard.npy $(FUNC)/aitchison.npy $(FUNC)/clr_groups.tsv $(TAXA_JACCARD) $(TAXA)/groups.tsv
	$(PY) $< "$$DISEASE" $(PERMS) $(MIN_PER_GROUP)

# ---- disease vs control per-sample counts within BioProjects (Mann-Whitney per project + van Elteren) ----
$(OUT)/within_project.tsv: 21_within_project.py $(SKETCH)/hash_counts.tsv $(FUNC)/distinct_per_sample.tsv $(TAXA)/distinct_per_sample.tsv
	$(PY) $< "$$DISEASE" $(MIN_PER_GROUP)

# ---- taxa relative abundance per sample: containment- and detection-score-based ----
ifeq ($(SITE),)
$(TAXA)/abundance.tsv.gz: 24_taxa_abundance.py $(TAXA_PRES) $(CONTAINMENT)
	$(PY) $< "$$DISEASE" $(DETECT_WEIGHTS)
endif
$(TAXA)/abundance_hist.png: plot_group_hist.py $(TAXA)/abundance.tsv.gz
	$(PY) $< $(TAXA)/abundance.tsv.gz $@ "$$DISEASE: YACHT taxa relative abundance proxies" --columns containment_abundance detection_abundance --labels "containment-based abundance" "detection-score-based abundance" --per-sample

# ---- per-taxon abundance association (species/genus/family, GTDB r232) + forest plots ----
$(GTDB):
	mkdir -p $(@D) && curl -sfL -o $@ https://data.gtdb.ecogenomic.org/releases/release232/232.0/$(@F)
.PRECIOUS: $(TAXA)/assoc_%.tsv.gz $(TAXA)/markers_%.tsv
$(TAXA)/assoc_%.tsv.gz $(TAXA)/markers_%.tsv: 25_taxa_assoc.py $(TAXA)/abundance.tsv.gz $(GTDB) literature_markers.tsv
	$(PY) $< "$$DISEASE" $* $(TAXA_MIN_PREV) $(MEMBERS)
$(TAXA)/forest_top_%.png: plot_forest.py $(TAXA)/assoc_%.tsv.gz
	$(PY) $^ $@ "$$DISEASE: top $(TOP_N) species by disease vs control q ($*)" --query "rank == 'species'" --top $(TOP_N)
$(TAXA)/forest_markers_%.png: plot_forest.py $(TAXA)/markers_%.tsv
	$(PY) $^ $@ "$$DISEASE: literature marker taxa and their most prevalent species ($*)"

# ---- per-study (BioProject) disease vs control differential abundance of literature-marker species ----
$(TAXA)/study_heatmap_%.png $(TAXA)/study_da_%.tsv $(TAXA)/study_report_%.txt: 26_study_heatmap.py $(TAXA)/abundance.tsv.gz $(GTDB) literature_markers.tsv
	$(PY) $< "$$DISEASE" $* $(STUDY_MIN_N) $(TAXA_MIN_PREV) $(STUDY_MIN_DETECT)

# ---- disease vs control classification: LASSO and random forest, stratified CV (studies pooled) ----
$(ML)/auroc_%_$(CLASSIFY_RANK).tsv: 27_classify.py $(TAXA)/abundance.tsv.gz $(GTDB)
	$(PY) $< "$$DISEASE" $* $(CLASSIFY_RANK) $(TAXA_MIN_PREV) $(PROCS) $(REPEATS)

# ---- clades ranked by agreement of LASSO, random forest and the univariate disease vs control test ----
$(ML)/ranked_%_$(CLASSIFY_RANK).tsv: 28_rank_taxa.py $(ML)/auroc_%_$(CLASSIFY_RANK).tsv $(TAXA)/assoc_%.tsv.gz literature_markers.tsv
	$(PY) $< "$$DISEASE" $* $(CLASSIFY_RANK)

# ---- publication tables of the classifiers, recomputed from the per-fold outputs ----
$(ML)/paper/tables_%_$(CLASSIFY_RANK).md: 31_paper_tables.py $(ML)/ranked_%_$(CLASSIFY_RANK).tsv
	$(PY) $< "$$DISEASE" $* $(CLASSIFY_RANK)

# ---- before/after body-site filter: all samples vs each of SITES (rank-taxa run per site by a sub-make) ----
$(BASE)/site_compare/compare_%_$(CLASSIFY_RANK).png: 30_compare_sites.py $(BASE)/ml/ranked_%_$(CLASSIFY_RANK).tsv
	for s in $(SITES); do $(MAKE) --no-print-directory rank-taxa SITE=$$s || exit 1; done
	$(PY) $< "$$DISEASE" $* $(CLASSIFY_RANK) $(SITES)

# ---- disease vs control neural net with a study (BioProject) adversary switched off (baseline) or on (DANN) ----
$(ML)/dann/auroc_%_$(CLASSIFY_RANK)_$(DANN_CV)_baseline.tsv: 34_dann.py $(TAXA)/abundance.tsv.gz $(GTDB)
	$(PY) $< "$$DISEASE" $* $(CLASSIFY_RANK) $(TAXA_MIN_PREV) off $(DANN_CV) $(REPEATS) $(DANN_LAMBDA)
$(ML)/dann/auroc_%_$(CLASSIFY_RANK)_$(DANN_CV)_dann.tsv: 34_dann.py $(TAXA)/abundance.tsv.gz $(GTDB)
	$(PY) $< "$$DISEASE" $* $(CLASSIFY_RANK) $(TAXA_MIN_PREV) on $(DANN_CV) $(REPEATS) $(DANN_LAMBDA)
$(ML)/dann/compare_%_$(CLASSIFY_RANK)_$(DANN_CV).png: 35_plot_dann.py $(foreach a,baseline dann,$(ML)/dann/auroc_%_$(CLASSIFY_RANK)_$(DANN_CV)_$(a).tsv)
	$(PY) $< "$$DISEASE" $* $(CLASSIFY_RANK) $(DANN_CV)
