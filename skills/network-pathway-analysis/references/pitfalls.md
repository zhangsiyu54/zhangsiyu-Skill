# Network and Pathway Analysis Pitfalls

Read this before finalizing conclusions from WGCNA, module-trait, hub gene, or GSVA/ssGSEA results. Each entry lists the trigger, the mechanism that causes the failure, how it shows up, and the fix.

## 1. Soft-thresholding power chosen without checking scale-free fit

- **Trigger**: picking a fixed/default soft power (e.g., always using 6) instead of running `pickSoftThreshold`.
- **Mechanism**: the soft power is meant to push the adjacency matrix toward approximate scale-free topology; skipping the diagnostic means the resulting network may still be close to random, so module detection reflects noise structure rather than biological co-expression.
- **Symptom**: scale-free topology fit index (signed R^2) never approaches roughly 0.8-0.9 at the chosen power, or mean connectivity is implausibly high/low; modules look arbitrary and enrichment is weak.
- **Fix**: run the fit-index-vs-power curve, choose the lowest power where the fit plateaus near the conventional 0.80-0.90 threshold (relax if the curve genuinely never reaches it, but report the achieved R^2 and treat modules as more exploratory), and always report the chosen power and the R^2 achieved.

## 2. Sample size too small for stable modules

- **Trigger**: running WGCNA on fewer than roughly 15-20 independent biological samples.
- **Mechanism**: pairwise correlation estimates are highly variable at low n; module assignment and eigengenes computed from noisy correlation matrices are not stable across resampling.
- **Symptom**: module structure changes substantially with small changes in filtering or with bootstrap/jackknife resampling; hub gene lists are not reproducible when a few samples are dropped.
- **Fix**: treat n<15 as exploratory only, state the sample count prominently, avoid strong causal or clinical claims, and where possible validate modules by resampling/bootstrap or an independent cohort (module preservation statistics, Langfelder et al. 2011).

## 3. Module-trait correlation without multiple-testing correction

- **Trigger**: reporting a module-trait correlation matrix (M modules x T traits) and highlighting "significant" cells using raw, uncorrected p-values.
- **Mechanism**: a full module-trait matrix is M*T simultaneous tests; uncorrected p<0.05 across dozens of cells produces false positives at the nominal rate per cell, not per matrix.
- **Symptom**: many marginally significant cells scattered across the heatmap, several of which do not replicate or make biological sense.
- **Fix**: apply FDR (Benjamini-Hochberg) or family-wise correction across the full module x trait matrix (not per-trait only), report both raw and adjusted p-values, and annotate the heatmap with adjusted-significance stars only.

## 4. Hub gene / kME conflated with causal driver

- **Trigger**: labeling the gene with highest intramodular connectivity (kME) or module membership as "the driver" or "the key regulator" of the phenotype.
- **Mechanism**: kME and gene significance (GS) describe network centrality and correlation strength within the observed co-expression structure; they are descriptive/statistical properties, not evidence of upstream regulatory or causal position.
- **Symptom**: mechanistic language ("gene X drives module Y", "regulates the pathway") appears in a report with no perturbation, eQTL, ChIP, or knockdown/knockout evidence.
- **Fix**: describe hub genes as "topologically central / highly connected candidates for follow-up," require external validation (perturbation, orthogonal cohort, prior literature, regulatory database) before any causal or "driver" language, and state the hub gene criteria used (kME threshold, GS threshold, or both).

## 5. GSVA/ssGSEA scores treated as directly comparable across independent runs

- **Trigger**: computing GSVA/ssGSEA scores in one batch of samples, then comparing those scores numerically to scores computed in a separate run (different cohort, different batch, added/removed samples later).
- **Mechanism**: both GSVA and ssGSEA are rank-based, sample-set-dependent methods — the score for a given sample depends on the empirical distribution of expression across all samples included in that specific run, so scores are relative to the run's sample composition, not an absolute or portable scale.
- **Symptom**: apparent "shifts" in pathway activity between two studies that actually reflect adding/removing samples or re-running the tool on a different cohort composition, not real biology.
- **Fix**: always compute scores for all samples to be compared in a single joint run; if new samples must be added later, recompute the full matrix rather than appending scores; never merge score tables from separate GSVA/ssGSEA invocations without re-running jointly, and state this limitation explicitly in methods.

## 6. GSVA/ssGSEA run on raw counts

- **Trigger**: feeding raw RNA-seq counts directly into `gsva()`/`ssgsea` without normalization or log transformation.
- **Mechanism**: these methods assume the input reflects each gene's relative expression distribution across samples in a roughly continuous, variance-stabilized scale; raw counts have mean-variance dependence and library-size artifacts that distort the rank/kernel estimation step.
- **Symptom**: pathway scores dominated by library size or a handful of very highly expressed genes; scores correlate with sequencing depth rather than biology.
- **Fix**: use normalized, log-transformed expression (e.g., log2(TPM+1), log2(CPM+1), variance-stabilizing transform, or the `kcdf="Poisson"` option in GSVA only when explicitly using integer counts as documented); state the normalization method in methods text.

## 7. Small gene sets producing unstable enrichment scores

- **Trigger**: scoring pathways/gene sets with fewer than roughly 10 genes (after matching to the expression matrix).
- **Mechanism**: GSVA/ssGSEA kernel density and rank-based statistics need enough genes per set to estimate a stable enrichment score; very small sets are dominated by one or two genes and are sensitive to missing/filtered genes.
- **Symptom**: pathway score for a small set swings drastically when one gene is dropped (different ID mapping, low detection) or when gene set version changes.
- **Fix**: filter gene sets to a minimum size (commonly 5-10 as a hard floor, 15+ preferred) and a maximum size before scoring, report the min/max size filter used, and flag any reported pathway near the size floor as lower-confidence.

## 8. Module merge threshold silently changes the biological narrative

- **Trigger**: choosing `mergeCutHeight` (or equivalent module-merging dissimilarity threshold) without justification, or trying several values until the "expected" modules appear.
- **Mechanism**: merging modules by eigengene correlation distance controls how many nearly-redundant modules get collapsed into one; a small threshold keeps many fine-grained modules, a large threshold aggressively merges modules that may actually carry distinct trait associations.
- **Symptom**: a trait-associated module of interest appears, disappears, or changes composition depending on the merge height; module count varies widely across re-analyses of the same data with different thresholds.
- **Fix**: use a conventional default (commonly 0.20-0.25 corresponding to eigengene correlation ~0.75-0.80) unless there is a stated biological or statistical reason to change it, report the exact value used, and check that the trait association of interest is not an artifact of a particular merge choice by inspecting the pre-merge module structure.

## 9. Block-wise WGCNA losing cross-block module relationships

- **Trigger**: using `blockwiseModules` on very large gene sets (commonly >20,000-30,000 genes or memory-limited environments), which splits genes into blocks before clustering.
- **Mechanism**: TOM and module detection are computed within each block independently (using a preliminary correlation-based pre-clustering to assign genes to blocks); genes that are truly co-expressed but assigned to different blocks cannot be merged into the same module, and cross-block correlation structure is not fully captured.
- **Symptom**: biologically related genes known to co-express end up in different modules with no shared module label, or module count/size looks artificially constrained by block boundaries.
- **Fix**: prefer single-block analysis when memory allows; if block-wise is required, increase `maxBlockSize` as far as resources permit, report the number of blocks and block assignment method used, and treat any conclusion that depends on genes being in different blocks as a limitation to flag explicitly.

## 10. Module eigengene-trait correlation treated as causal

- **Trigger**: phrasing results as "module X causes/drives/regulates trait Y" based solely on eigengene-trait correlation and its p-value.
- **Mechanism**: WGCNA module-trait correlation is an association between a co-expression summary (the eigengene, the first principal component of module gene expression) and a phenotype measured in the same cross-sectional samples; it carries no temporal or interventional information.
- **Symptom**: causal verbs in the report ("module drives," "leads to," "controls") attached only to a correlation coefficient and p-value.
- **Fix**: use association language ("associated with," "correlated with," "co-varies with"), state the correlation method (Pearson/Spearman/biweight midcorrelation) and sample size next to every reported correlation, and reserve causal claims for designs that support them (perturbation, longitudinal, Mendelian randomization, or prior mechanistic literature).

## 11. Combining heterogeneous omics types in one WGCNA network without harmonizing distributions

- **Trigger**: building a single co-expression/correlation network directly from a matrix that concatenates transcriptomics (log counts), proteomics (log intensities), and/or metabolomics (abundances) without checking comparability of dynamic range, missingness, and distribution shape.
- **Mechanism**: WGCNA's correlation and adjacency calculations implicitly assume comparable behavior across features; features with very different dynamic range, sparsity, or measurement noise (e.g., proteomics with many missing values vs. RNA-seq with near-complete detection) can dominate or distort clustering, and cross-omics correlations may reflect measurement-platform artifacts rather than shared biology.
- **Symptom**: modules split cleanly along omics-type boundaries rather than by biological pathway, or a handful of high-variance features from one platform dominate several modules.
- **Fix**: scale/normalize each omics layer appropriately for its own distribution before combining (or run within-platform WGCNA and integrate afterward via eigengene correlation, DIABLO/MOFA-style multi-omics integration, or explicit multi-omics network methods), document missing-data handling per platform, and check whether module composition is confounded by omics type before interpreting cross-omics modules biologically.

## 12. Signed vs. unsigned network choice not stated or not matched to the question

- **Trigger**: defaulting to an unsigned network (or not stating which type was used) when the biological question depends on the direction of co-expression, or vice versa.
- **Mechanism**: unsigned networks treat strong positive and strong negative correlations as equally "connected," which can merge genes that are anti-correlated into the same module and obscure directional relationships (e.g., activator vs. repressor); signed (or signed-hybrid) networks only treat positively correlated genes as connected, which better preserves direction but may split biologically coupled but anti-correlated genes into separate modules.
- **Symptom**: a module contains genes with opposite-direction expression patterns relative to the trait, muddying eigengene interpretation, or expected antagonistic pairs never appear in the same module when they should be analyzed together.
- **Fix**: use signed or signed-hybrid networks as the default for most trait-correlation work (recommended by WGCNA's authors for most applications) unless the question specifically calls for unsigned "magnitude of co-regulation," and always state which network type was used and why.

## 13. Sample outliers and batch structure not inspected before network construction

- **Trigger**: proceeding directly to soft-threshold selection and module detection without first clustering samples (e.g., hierarchical clustering or PCA) to check for outlier arrays/libraries or batch/technical grouping.
- **Mechanism**: WGCNA correlations are sensitive to extreme outlier samples and to strong non-biological (batch/technical) variance components; a single aberrant sample can dominate pairwise correlations and distort adjacency estimation, and unaddressed batch effects can generate modules that track technical batch rather than phenotype.
- **Symptom**: one or two samples sit far from all others on a sample dendrogram or PCA plot; a module's eigengene correlates strongly with a batch/processing-date variable rather than (or in addition to) the trait of interest.
- **Fix**: always run and report a sample clustering/PCA QC step before network construction, remove or justify retaining extreme outliers, and check module eigengenes against known technical covariates (batch, processing date, RIN/DV200, sequencing depth) in addition to biological traits — adjusting or reporting confounding explicitly when found.

## 14. Choice of pathway scoring method not matched to gene-set count or sample count

- **Trigger**: using GSVA/ssGSEA by default regardless of how many samples are available or how the pathway scores will be used downstream (e.g., true single-sample clinical scoring vs. batch cohort comparison).
- **Mechanism**: GSVA and ssGSEA estimate each sample's score using the rank/density of expression relative to the other samples in the same run, so results depend on sample-set composition and are less stable with very few samples; genuinely sample-independent methods (e.g., singscore, z-score-based methods) are constructed to not require a background cohort and may be preferable for true single-sample/clinical deployment.
- **Symptom**: a "single-sample score" is reported for a case with no comparator cohort, or scores are re-derived per patient in a clinical pipeline using GSVA/ssGSEA and drift when the reference cohort composition changes.
- **Fix**: use GSVA/ssGSEA for cohort-level relative pathway activity comparison within one consistent sample set; for genuinely single-patient or per-sample deployment scenarios, prefer scoring methods designed to be sample-independent (e.g., singscore, or GSVA/ssGSEA re-run against a fixed, versioned reference cohort) and document which regime applies.
