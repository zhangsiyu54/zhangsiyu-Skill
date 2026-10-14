# CNS Method Guide

Use this guide before CRISPR knockout, activation, inhibition, or drug screen hit analysis.

## Minimal Inputs

- Require sgRNA count matrix or gene-level summary, sample metadata, treatment/control labels, library type, gene/sgRNA annotation, replicate structure, and selection direction.
- For MAGeCK-style outputs, require gene summary, sgRNA summary if available, beta/logFC, P value, FDR, and positive/negative selection columns.
- Confirm whether the screen is viability/dropout, positive selection, drug resistance/sensitization, CRISPRi, or CRISPRa.

## QC and Preprocessing

- Check read depth, mapped sgRNA fraction, missing sgRNAs, replicate correlation, library representation, Gini/skewness, and essential/nonessential control behavior when available.
- Normalize counts before testing and filter extremely low-count sgRNAs using explicit rules.
- Keep replicate-level consistency visible; do not trust hits driven by a single replicate or single sgRNA.

## Hit Calling

- Use MAGeCK RRA/MLE, BAGEL, drugZ, or equivalent tools when raw counts and design are present.
- Report positive and negative selection separately.
- Require multiple sgRNAs per gene where possible; flag genes supported by too few guides.
- Use FDR correction and provide effect size plus statistical evidence.
- For drug screens, model treatment-specific depletion/enrichment relative to control and avoid calling generic essential genes as drug-specific hits without interaction evidence.

## Downstream Interpretation

- Prioritize hits with consistent guide direction, known pathway coherence, control behavior, and enrichment support.
- Enrichment should use screened genes as background, not the whole genome.
- State whether results are hypothesis-generating and require validation.

## Required Figures

- QC depth/correlation, volcano plot, ranked gene plot, top hit barplot, sgRNA consistency plot when available, pathway enrichment dot/bar plot.

## Methods Text Must Include

- Library/screen type, contrast, normalization, hit-calling method, replicate design, sgRNA filtering, FDR threshold, direction definition, background universe for enrichment, and package versions.

## Algorithm and Design Decision Tree

Use this table to pick the method/parameter combination before running scripts. Read `references/pitfalls.md` alongside it — several rows correspond directly to a numbered pitfall.

| Scenario | Tool/method | Key parameter(s) | Rationale |
| --- | --- | --- | --- |
| Positive selection (dropout survivors, drug resistance) screen | `mageck test` (RRA), positive-selection mode | `--gene-lfc-method median` or default RRA positive rank; report `pos|score`/`pos|fdr` columns | RRA aggregates ranks per direction; positive-selection sensitivity is generally lower than negative, so require concordant sgRNAs and treat marginal FDR hits cautiously (pitfall 2). BAGEL/BAGEL2 is not applicable here — it is a depletion-only classifier. |
| Negative selection (essentiality/dropout) screen with a nonessential/essential reference set available | BAGEL2 (preferred when reference gene sets exist) or MAGeCK RRA negative mode | BAGEL2 Bayes Factor threshold set via cross-validation (tool default), or `neg|fdr` from `mageck test` | BAGEL2's reference-anchored Bayesian model gives a calibrated Bayes Factor and lower false positives for depletion than rank-based RRA when curated essential/nonessential gene lists are available (pitfall 3). |
| Negative selection screen with no reference essential/nonessential gene list | MAGeCK RRA (negative mode) | Default `mageck test`, explicit low-count sgRNA filter | Without a calibration set, a purely rank-based method is safer than a Bayesian classifier that needs a reference prior; label output as an uncalibrated rank list, not a validated essential-gene set. |
| Multi-condition or covariate-adjusted design (e.g., multiple timepoints, dose series, paired T0) | MAGeCK MLE | Design matrix with explicit time/passage/dose covariates; inspect per-condition beta scores | MLE's NB GLM can model shared and condition-specific effects simultaneously, letting you subtract generic passage/fitness drift instead of conflating it with treatment effect (pitfall 5). |
| Genome-wide screen in an aneuploid/CNV-rich cancer cell line | CERES or CRISPRcleanR (or Chronos if longitudinal data exist), applied before or alongside MAGeCK/BAGEL2 scoring | Supply cell-line copy-number segments (CERES/Chronos) or rely on CRISPRcleanR's unsupervised sgRNA-level segmentation if CN data are unavailable | Cas9 cutting toxicity scales with local copy number, not gene essentiality; uncorrected knockout screens in amplified regions produce clustered false-positive "essential" calls (pitfall 10). |
| Drug/chemogenetic interaction screen (treatment vs. vehicle, same T0) | drugZ | Paired treatment-vs-control input per replicate; inspect `normZ` and FDR | drugZ was purpose-built for treatment-vs-control interaction scoring and explicitly separates drug-specific hits from generic fitness genes better than repurposing a plain essentiality tool. |
| CRISPRi or CRISPRa screen | MAGeCK MLE or RRA, but with CRISPRi/a-designed library and QC | Check guide-to-TSS distance for top hits; do not assume knockout-style exon-tiling logic | Guide efficacy in CRISPRi/a depends on TSS-proximal positioning rather than frameshift logic; treating guide behavior like a knockout screen inflates apparent sgRNA inconsistency (pitfall 9). |
| Low library representation / coverage QC before any hit calling | `scripts/check_environment.py` + count-distribution QC (Gini index, zero-count sgRNA fraction) at T0 | Compare T0/plasmid-pool coverage against the intended cells-per-sgRNA target for the library size used | Bottleneck-driven sgRNA loss at transduction is indistinguishable from true depletion later; catching it at T0 prevents misattributing a coverage artifact to biology (pitfall 1). |
| Combining hits across several independent contrasts (multiple drugs/cell lines/timepoints) into one summary | Per-contrast FDR from the tool of choice, plus explicit note on family-wise inflation | Do not silently pool per-contrast FDR-significant genes into a single "hit in any condition" claim without adjustment | Per-comparison FDR does not control the false-discovery rate of the pooled claim across N contrasts (pitfall 7). |
