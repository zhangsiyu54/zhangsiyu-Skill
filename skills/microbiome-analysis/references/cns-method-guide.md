# CNS Method Guide

Use this guide before running 16S, ASV/OTU, microbiome diversity, differential taxa, or functional profile analysis.

## Minimal Inputs

- Require feature table, taxonomy table, sample metadata, group column, feature type, sequencing region/platform, and whether values are raw counts, relative abundance, or transformed.
- For raw reads, require primer/region information and target pipeline such as QIIME2/DADA2.
- For shotgun functional profiles, require pathway/module table and database version.

## QC and Preprocessing

- Track reads/features per sample, rare taxa, unassigned taxonomy, negative/blank controls when available, and sample metadata completeness.
- Do not rarefy by default for differential abundance; rarefaction may be used for alpha/beta diversity comparability if justified.
- Filter contaminants, low-depth samples, and extremely rare features using explicit rules.
- Use compositional data awareness: relative abundance is constrained and should not be treated as absolute biomass.

## Diversity and Community Testing

- Alpha diversity: report metric, test, group sample counts, and multiple-testing correction if multiple metrics are tested.
- Beta diversity: use appropriate distances such as Bray-Curtis, Jaccard, weighted/unweighted UniFrac when phylogeny exists.
- PERMANOVA must include permutations, strata/blocking when relevant, and dispersion check when possible.

## Differential Abundance

- Prefer ANCOM-BC, ALDEx2, MaAsLin2, or DESeq2-style models depending on data and covariates.
- Avoid simple t-tests on relative abundance for formal claims unless labeled exploratory.
- Report taxonomic rank, prevalence filters, structural zero handling when applicable, effect direction, FDR, and covariates.

## Required Figures

- Read-depth/QC plot, alpha diversity boxplot, PCoA/NMDS ordination, top taxa stacked bar, taxa heatmap, differential taxa volcano/bar plot, functional pathway plot when applicable.

## Methods Text Must Include

- Sequencing type/region, denoising/clustering method, taxonomy classifier/database, filtering thresholds, normalization/transformation, diversity metrics, PERMANOVA formula, differential method, FDR, package versions, and limitations.

## Algorithm and Design Decision Tree

Before running an analysis, match the scenario below to the recommended method/tool and parameters, and state the rationale in the methods text. See `references/pitfalls.md` for the failure modes each choice is meant to avoid.

| Scenario | Recommended method/tool + key parameters | Rationale |
| --- | --- | --- |
| 16S ASV/OTU differential abundance with compositionality | ANCOM-BC (or ANCOM-BC2)/ALDEx2 with explicit prevalence filter (e.g., taxa in >=10% of samples), structural-zero handling, CLR/log-ratio bias correction, BH-FDR | Corrects the compositional bias and unequal sampling fraction that distort naive relative-abundance tests; structural-zero test avoids forcing a fold-change on taxa that are truly absent in one group |
| Shotgun metagenomics functional profiling | HUMAnN3 against UniRef90/MetaCyc (or ChocoPhlAn for taxonomy), normalize pathway/gene-family abundance to copies-per-million (CPM), test associations with MaAsLin2 | CPM removes sequencing-depth differences before comparing pathway abundance; stratified HUMAnN output lets taxon-level contributions to a pathway shift be checked before over-interpreting a functional change |
| Low-biomass sample with contamination risk | Sequence extraction/PCR negative controls; run decontam (prevalence or frequency method) before any diversity/differential step | Reagent and kit-derived DNA ("kitome") can dominate true signal in low-biomass samples (Salter et al. 2014); without blanks, contaminants are indistinguishable from biology |
| Longitudinal/repeated-measures design | PERMANOVA with `strata=subject` (or a mixed model); MaAsLin2/lme4 with a random intercept (and slope if warranted) for subject; report within- and between-subject variance | Repeated samples per subject are correlated; treating them as independent inflates degrees of freedom and produces anti-conservative p-values (pseudoreplication) |
| Beta diversity choice by biological question | Unweighted UniFrac for phylogenetic membership/rare-taxon questions; weighted UniFrac for phylogenetically-weighted abundance shifts; Bray-Curtis for compositional/non-phylogenetic shifts; Jaccard for presence/absence without phylogeny | Each metric answers a different hypothesis (Lozupone & Knight 2005); the choice should follow the biological question, not be picked post hoc after seeing which one is significant |
| Cross-cohort/cross-study comparison | Re-classify all cohorts with one shared reference database + version (SILVA, Greengenes2, or GTDB); apply batch correction (e.g., MMUPHin, ComBat-seq) or random-effects meta-analysis across studies | Database/version drift (SILVA vs Greengenes vs GTDB) and technical batch differences are common confounders that make raw cross-study taxon comparisons invalid unless harmonized |
| PERMANOVA with possible dispersion confound | Run `betadisper`/PERMDISP2 alongside `adonis2`; report centroid (PERMANOVA) and dispersion results together; note permutations and any blocking/strata | A significant PERMANOVA can reflect unequal within-group dispersion rather than a true centroid shift; both must be reported to support a compositional-difference claim |
| Raw amplicon reads with run/region variability | DADA2 with per-sequencing-run error learning, primer/region-specific truncation lengths, chimera removal (`removeBimeraDenovo`); do not merge ASVs across different variable regions | The DADA2 error model is run/instrument-specific; different hypervariable regions are not directly comparable at the ASV sequence level |
| Absolute abundance / microbial load question | Pair relative-abundance sequencing with 16S qPCR or spike-in standards to scale relative abundances into absolute/estimated load | Sequencing alone yields compositional (relative) data; a claimed increase in one taxon cannot be distinguished from a decrease elsewhere without an orthogonal absolute measurement |
| Small sample size or sparse compositional differential abundance | ALDEx2 with Monte Carlo Dirichlet instances (default 128+) and CLR transform; interpret effect size (median CLR difference) alongside BH-FDR | Propagates count-estimation uncertainty into the test statistic, giving more conservative, reproducible calls in small-n studies where benchmarking shows many methods over-call significance (Nearing et al. 2022) |
