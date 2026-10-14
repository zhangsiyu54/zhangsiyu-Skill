# CNS Method Guide

Use this guide before running co-expression network, module, or pathway activity analysis.

## Minimal Inputs

- Require normalized/log-transformed expression or omics matrix, sample metadata, phenotype columns, organism, gene ID type, and pathway gene sets when pathway scoring is requested.
- WGCNA requires enough independent biological samples. Treat very small cohorts as exploratory and avoid overfitting module-trait claims.
- Confirm whether data are bulk, pseudobulk, proteomics, metabolomics, or already summarized module scores.

## Matrix QC

- Remove low-variance and low-detection features before WGCNA or correlation networks.
- Check sample clustering, outliers, batch effects, missing values, and phenotype distributions.
- Use variance-stabilized/log expression for correlation-based methods.
- Document imputation only when unavoidable; avoid imputing large missing blocks for network inference.

## WGCNA and Network Rules

- Choose soft-thresholding power using scale-free topology diagnostics and mean connectivity, not a fixed default.
- Construct signed or signed-hybrid networks unless the biological question requires unsigned relationships.
- Report module eigengenes, module sizes, module-trait correlations, adjusted P values, and hub gene criteria.
- Validate important modules with enrichment, external data, or literature consistency when possible.

## Pathway Activity Rules

- Use GSVA/ssGSEA for sample-level pathway scores from normalized expression; do not run these on raw counts.
- Use curated gene sets, record collection/version, and filter very small/large gene sets.
- For group comparison of pathway scores, use appropriate tests or linear models with covariates.
- Interpret pathway scores as relative activity signatures, not direct biochemical flux.

## Required Figures

- Sample clustering/PCA, soft-threshold plot, module dendrogram, module-trait heatmap, hub gene/network summary, pathway score heatmap, pathway boxplot.

## Methods Text Must Include

- Matrix preprocessing, feature filtering, correlation type, network type, soft power, module merge threshold, trait model/test, pathway scoring method, gene set source/version, thresholds, and package versions.

## Algorithm and Design Decision Tree

Read `references/pitfalls.md` before finalizing conclusions from any row below.

| Scenario | Recommended tool/method + key parameter | Rationale |
| --- | --- | --- |
| Small sample size unsuitable for full WGCNA (roughly <15 samples) | Do not run formal WGCNA module detection; instead report a correlation heatmap of top variable genes/pathway scores against traits, or a simple hierarchical clustering, and label results exploratory | Correlation-based module detection is unstable below ~15 independent samples; forcing WGCNA invites non-reproducible modules and overfit hub genes (see pitfalls.md #2) |
| Choosing signed vs. unsigned network | Default to signed or signed-hybrid network (`networkType = "signed"` or `"signed hybrid"` in `adjacency()`/`blockwiseModules`) unless the question explicitly needs magnitude-only co-regulation | Signed networks preserve the direction of co-expression so modules and eigengenes are interpretable with respect to up/down regulation; unsigned networks can merge anti-correlated genes into one module (pitfalls.md #12) |
| Single-sample pathway activity scoring for cohort comparison | GSVA (Poisson kernel for counts, Gaussian kernel for log-normalized data) or ssGSEA, run jointly on the full comparison cohort in one call | Rank/kernel-based scores are sample-set-dependent; scoring the full cohort together keeps scores internally comparable within that run |
| Single-sample pathway score for true per-patient/clinical deployment | Prefer a sample-independent scoring method (e.g., singscore) or GSVA/ssGSEA re-run against a fixed, versioned reference cohort each time | GSVA/ssGSEA scores are not portable across independent runs; a genuinely single-sample method avoids recomputation drift (pitfalls.md #14) |
| Multi-trait module correlation with correction | Compute the full module x trait correlation matrix (Pearson or biweight midcorrelation `bicor` for outlier robustness), then apply Benjamini-Hochberg FDR across all cells of the matrix, not per-trait | Testing M modules x T traits is a large multiple-comparison problem; per-trait-only correction understates the true false-discovery rate (pitfalls.md #3) |
| Consensus WGCNA across multiple datasets | Use `blockwiseConsensusModules` (or manual consensus TOM construction per Langfelder & Horvath 2007) with matched gene sets across datasets and comparable normalization per dataset | Consensus modules require genes/probes common to all input datasets and comparable per-dataset scaling; naive pooling of raw matrices from different platforms/normalizations biases the consensus network |
| Hub gene prioritization | Rank by intramodular connectivity (kME, i.e., correlation of gene expression to module eigengene) combined with gene significance (correlation to trait); report both thresholds used | kME alone reflects network centrality, not phenotype relevance; combining kME with GS gives "trait-relevant hub" candidates while keeping causal language out (pitfalls.md #4) |
| Very large gene universe (>20,000-30,000 genes) or memory-limited environment | `blockwiseModules` with `maxBlockSize` set as high as available memory allows; report number of blocks and block-assignment method | Block-wise clustering trades completeness for tractability; genes split across blocks cannot be merged into one module, so block count/size must be disclosed as a limitation (pitfalls.md #9) |
| RNA-seq raw counts as WGCNA input | Variance-stabilizing transform (DESeq2 `vst`/`rlog`) or log2(CPM/TPM+1) before correlation/adjacency calculation; never feed raw counts directly | Pearson/bicor correlation assumes roughly homoskedastic, continuous-scale input; raw counts have strong mean-variance dependence that distorts adjacency estimation |
| Proteomics/metabolomics or mixed-omics network construction | Normalize/scale each omics layer within its own platform before combining; consider per-platform WGCNA plus eigengene-level integration instead of one pooled matrix | Platforms differ in dynamic range, missingness, and noise structure; pooling without harmonization lets one platform's variance dominate clustering (pitfalls.md #11) |
| Differential module-trait association across case/control with covariates (age, sex, batch) | Linear model of module eigengene ~ trait + covariates (e.g., `lm`/`limma` on eigengenes) rather than a single unadjusted correlation | Unadjusted eigengene-trait correlation can be confounded by known covariates; modeling covariates explicitly separates trait signal from nuisance variance |
