# CNS Method Guide

Use this guide before proteomics, metabolomics, or multi-omics table analysis.

## Minimal Inputs

- Require feature abundance matrix, sample metadata, platform, normalization status, feature annotation, group/contrast, missing-value meaning, and batch/run columns if present.
- For proteomics, distinguish LFQ, TMT/iTRAQ, DIA, spectral counts, peptide-level, and protein-level data.
- For metabolomics, distinguish targeted/untargeted, positive/negative mode, feature ID confidence, and metabolite annotation level.

## QC and Preprocessing

- Summarize missingness by sample and feature, total intensity, detected features, PCA, correlation, and batch/run effects.
- Filter features with excessive missingness using explicit group-aware rules.
- Normalize according to platform: median/quantile/variance stabilization for LFQ-style matrices; TMT channel normalization for multiplexed designs; metabolomics TIC/PQN/log transform/scaling when justified.
- Impute only after filtering and record method. Left-censored missingness should use low-intensity methods; random missingness may use KNN/model-based approaches.

## Differential Analysis

- Use limma/MSstats/DEP-style moderated models for proteomics when possible.
- For metabolomics, use linear models or nonparametric tests with FDR correction; do not rely only on VIP scores.
- Include covariates/batches when known. Keep all tested features in result tables.
- Validate annotation and direction before pathway claims.

## Multi-omics Integration

- Align sample IDs across omics layers before integration.
- Use correlation, pathway-level summaries, MOFA/mixOmics-style methods only when sample size supports them.
- Interpret cross-omics correlations as associations unless externally validated.

## Required Figures

- Missingness heatmap/bar, PCA, sample correlation, volcano plot, top feature heatmap, pathway enrichment plot, multi-omics correlation/network plot when applicable.

## Methods Text Must Include

- Platform, acquisition/quantification type if known, normalization, filtering, imputation, model/test, FDR method, annotation database, pathway method, thresholds, package versions, and limitations.

## Algorithm and Design Decision Tree

Before running an analysis, use this table to pick the method/parameter and justify it in the methods text. See `references/pitfalls.md` for the failure modes each choice is meant to avoid.

| Scenario | Recommended approach | Rationale |
|---|---|---|
| TMT/iTRAQ vs LFQ differential protein workflow | TMT: MSstatsTMT/limma on log2 reporter ratios referenced to a pooled bridge channel, with per-plex median or IRS (internal reference scaling) normalization across plexes; LFQ: MaxLFQ intensities (MaxQuant) or MSstats with peptide-level rollup to protein, moderated (eBayes) t-test | TMT reporter ions suffer co-isolation ratio compression and need a shared reference to bridge plexes; LFQ intensities are directly comparable across runs only after MaxLFQ-style delayed normalization, and moderated variance shrinkage compensates for typically small per-group n |
| Metabolomics untargeted vs targeted normalization | Untargeted (urine/serum, no internal standard per analyte): probabilistic quotient normalization (PQN) after integral normalization, plus log transform and Pareto/auto-scaling for multivariate models; targeted (isotope-labeled internal standards available): analyte/internal-standard ratio, no PQN needed | PQN is robust to a few dominant, biologically irrelevant dilution drivers, which is the main failure mode of simple total-intensity/TIC scaling in dilution-variable biofluids like urine; targeted assays already have an analyte-specific internal correction so a global dilution correction is redundant or can even distort the ratio |
| Missing value pattern MNAR vs MAR | MNAR (missingness increases as abundance falls, or "all-or-nothing" between conditions): left-censored imputation (MinProb, QRILC, deterministic minimum/half-minimum); MAR/MCAR (missingness unrelated to abundance, scattered across the intensity range): KNN, missForest, or maximum-likelihood/model-based imputation | Left-censored methods place imputed values in the low-abundance tail consistent with a true below-detection-limit signal; KNN/model-based methods borrow information from correlated features/samples and are appropriate only when the missingness mechanism is not abundance-driven |
| Multi-omics integration approach choice | Small n (roughly <15/group) or exploratory stage: pairwise correlation and pathway/module-level summary only; larger, well-matched multi-omics cohort: MOFA2 (unsupervised latent factors) or DIABLO/mixOmics (supervised, when a clear outcome variable exists) | Latent-factor/supervised multivariate integration models are prone to overfitting and unstable factor loadings at small n; correlation and pathway-level summaries are more conservative and easier to sanity-check when sample size is limited |
| Pathway enrichment for metabolites | Compound-ID-based enrichment (e.g., MetaboAnalyst MSEA/pathway analysis) using a verified name-to-InChIKey-to-KEGG/HMDB crosswalk, with the measured/detectable feature list (not the whole database) as background | Metabolite pathway databases key on exact compound IDs; using the full theoretical compound universe as background inflates significance, and unresolved ID mapping silently drops or duplicates hits |
| Protein rollup with shared peptides (paralogs/isoforms) | Razor + unique peptide quantification (assign shared peptides to the protein group with strongest evidence) rather than naive peptide-intensity summation across all matching proteins | Naive summation double-counts shared peptides and can fabricate apparent quantification of a low-abundance paralog that was never actually distinguished from a more abundant relative |
| Small sample size with many tested features | Moderated/empirical-Bayes variance shrinkage (limma-style eBayes as in DEP/MSstats) with BH-FDR as the primary threshold, not raw t-test with uncorrected p-values or Bonferroni | Per-feature variance estimates are unstable at small n; borrowing variance information across features improves power while FDR control keeps the false-positive rate calibrated for thousands of simultaneous tests |
| Imputation timing relative to peptide-to-protein aggregation | Impute at the same data level (peptide or protein) at which statistical testing will be performed; do not impute at peptide level then re-aggregate with a different assumption, or vice versa, without re-checking the missingness mechanism at that level | Imputation performance is data-level dependent (Lazar et al. 2016); mixing levels propagates bias from whichever level was imputed into the level being tested |
| TMT bridging across multiple plexes/batches | Common pooled reference/bridge channel in every plex plus internal reference scaling (IRS); avoid direct cross-plex ratio comparison without a bridge | Without a shared reference, plex-to-plex technical variation (channel efficiency, batch) is indistinguishable from the biological effect being tested |
| Reporting identification confidence for named metabolites | Carry an MSI-level (1-4) column through the pipeline; confirm MSI level 1 (matched to authentic standard) for any metabolite central to the paper's conclusion | Pathway/mechanistic claims are only as strong as the weakest identification in the chain; level 2-4 annotations should be hedged accordingly |
