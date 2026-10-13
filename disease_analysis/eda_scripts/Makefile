# EDA of one disease vs its controls vs other-disease controls.
#
#   make sample-sketch-similar-viz Colorectal cancer    # sketch Jaccard similarity
#   make taxa-pres-abs Colorectal cancer                # YACHT taxa presence/absence
#   make funprofile-pres-abs Colorectal cancer          # KO (function) presence/absence
#   make funprofile-abundance Colorectal cancer         # KO CLR / Aitchison distance
#   make seq-depth Colorectal cancer                    # sequencing depth per group (taxa and KO samples)
#   make metadata-eda Colorectal cancer                 # group sizes, BioProjects, geography, sequencing, host covariates
#   make stratified-permanova Colorectal cancer         # disease vs control PERMANOVA before/after accounting for BioProject
#   make count-tests Colorectal cancer                  # Kruskal-Wallis + Dunn on per-sample counts, then disease vs control within BioProjects
#   make taxa-abundance Colorectal cancer               # YACHT taxa relative abundance: containment- and detection-score-based
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

DISEASE_TARGETS := sample-sketch-similar-viz taxa-pres-abs funprofile-pres-abs funprofile-abundance seq-depth metadata-eda stratified-permanova count-tests taxa-abundance

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
export DISEASE

OUT     := results/$(subst ',,$(subst $(space),_,$(strip $(DISEASE))))
SKETCH  := $(OUT)/sketch
TAXA    := $(OUT)/taxa
FUNC    := $(OUT)/func
META    := $(OUT)/metadata
LABELS  := results/disease_labels.tsv
COVERAGES := 1.0 0.5 0.25 0.125 0.0625 0.03125 0.015625
TAXA_JACCARD   := $(foreach c,$(COVERAGES),$(TAXA)/jaccard_$(c).npy $(TAXA)/groups_$(c).tsv)
TAXA_PERMANOVA := $(foreach c,$(COVERAGES),$(TAXA)/permanova_$(c).tsv)
NOHUMAN := results/nohuman_unitigs.done
YACHT   := results/yacht_presence.tsv.gz
CONTAINMENT := results/yacht_containment.tsv.gz

.PHONY: $(DISEASE_TARGETS) labels remove-human yacht-presence yacht-containment list-diseases

sample-sketch-similar-viz: $(SKETCH)/jaccard_heatmap.png $(SKETCH)/jaccard_summary.tsv $(SKETCH)/permanova.tsv $(SKETCH)/hash_counts.png
taxa-pres-abs: $(TAXA)/prevalence.tsv $(TAXA)/distinct_per_sample.png $(TAXA_PERMANOVA)
funprofile-pres-abs: $(FUNC)/prevalence.tsv $(FUNC)/distinct_per_sample.png $(FUNC)/kruskal_dunn.tsv
funprofile-abundance: $(FUNC)/aitchison_heatmap.png $(FUNC)/aitchison_summary.tsv $(FUNC)/aitchison_permanova.tsv
seq-depth: $(TAXA)/depth.png $(FUNC)/depth.png
metadata-eda: $(META)/summary.txt
stratified-permanova: $(OUT)/stratified_permanova.tsv
count-tests: $(SKETCH)/kruskal_dunn.tsv $(TAXA)/kruskal_dunn.tsv $(FUNC)/kruskal_dunn.tsv $(OUT)/within_project.tsv
taxa-abundance: $(TAXA)/abundance.tsv.gz

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

$(TAXA)/groups.tsv: 02_split_groups.py $(LABELS) $(YACHT)
	$(PY) $< "$$DISEASE" taxa
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
$(TAXA)/abundance.tsv.gz: 24_taxa_abundance.py $(TAXA_PRES) $(CONTAINMENT)
	$(PY) $< "$$DISEASE" $(DETECT_WEIGHTS)
