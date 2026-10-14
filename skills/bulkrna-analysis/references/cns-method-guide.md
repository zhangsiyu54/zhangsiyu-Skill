# CNS Method Guide

Use this guide before running or reporting bulk RNA-seq analysis. Treat it as the methods and QC backbone for manuscript-style delivery.

## Minimal Inputs

- Raw integer gene counts, sample metadata, group labels, organism, gene ID type, and contrast definition are required for formal differential expression.
- Normalized expression matrices can be used for visualization, PCA, clustering, correlation, GSVA/WGCNA, or exploratory differential summaries, but do not report DESeq2/edgeR statistics from already-normalized values.
- Require batch, paired subject, sex, age, tissue, treatment time, or other covariates when the user expects adjusted modeling.

## QC and Filtering

- Check sample ID identity and ordering before modeling.
- Summarize library sizes, detected genes, zero fraction, group balance, batch columns, and replicate count per group.
- Filter low-count genes before modeling. Prefer a rule tied to library depth and group size, such as retaining genes with counts above a minimal CPM/count threshold in enough samples.
- Use PCA, sample correlation heatmap, and hierarchical clustering to detect outliers, swaps, severe batch effects, and mislabeled replicates.
- Do not remove samples silently. Flag suspicious samples and document the reason.

## Modeling

- Use DESeq2 for standard two-group or multifactor raw-count differential expression; use design formulas that include known covariates.
- Use limma-voom or edgeR for complex designs, paired designs, small cohorts with covariates, or when the requested contrast is easier to express with linear models.
- For paired or repeated-measure designs, include subject/blocking terms rather than treating samples as independent.
- For time-course or multi-level factors, define explicit contrasts or use likelihood-ratio style tests only when the biological question asks for any effect across levels.
- Use Benjamini-Hochberg FDR for multiplicity. Report both nominal P value and adjusted P value.
- Use shrunken log2 fold change for ranking/visualization when available; do not confuse shrinkage estimates with the p-value test model.

## Thresholds and Ranking

- Default reporting threshold: adjusted P value < 0.05 plus an effect-size threshold such as absolute log2FC >= 1, unless the user specifies otherwise.
- For enrichment, pass ranked statistics for GSEA and separated up/down significant gene lists for ORA.
- Preserve all tested genes in the full result table, not only significant genes.

## Required Figures

- QC summary: library size, detected genes, PCA, sample correlation or clustering.
- Differential summary: volcano plot, MA plot when possible, top DEG heatmap, DEG count bar.
- Downstream handoff: ranked list, up/down gene lists, and enrichment-ready ID mapping table.

## Interpretation Boundaries

- Do not infer causality from differential expression alone.
- Do not over-interpret pathways supported by only one or two genes.
- State when the analysis is underpowered, lacks biological replicates, uses normalized-only data, or contains unresolved batch effects.

## Methods Text Must Include

- Input file type, genome/gene annotation version when known, filtering rule, normalization method, model package, design formula, contrast, FDR method, thresholds, package versions, and exact date/databases used for enrichment.

## Algorithm and Design Decision Tree

Use this table to pick the concrete method/function/parameter for a given scenario. See `references/pitfalls.md` for the failure modes behind several of these choices.

| Scenario | Recommended Method / Function / Parameter | Rationale |
|---|---|---|
| Two-group raw counts, n>=3/group | DESeq2 `DESeq(dds)` + `results(dds, contrast=c("condition","case","control"))`, then `lfcShrink(dds, res=res, type="apeglm")` for reported LFC | Explicit contrast avoids picking the wrong default coefficient; apeglm shrinkage stabilizes LFC for low-count genes without altering the Wald p-value |
| 3+ level factor, testing "any change across levels" | DESeq2 LRT: `DESeq(dds, test="LRT", reduced=~1)` for the omnibus p-value; report per-level Wald contrasts separately for effect direction/magnitude | LRT gives one FDR-controlled test per gene for "does this factor matter at all"; its logged LFC column is only one coefficient, so effect size must come from separate pairwise contrasts |
| Interaction design (genotype x treatment) | `design = ~ genotype + treatment + genotype:treatment`; test the interaction via LRT dropping the interaction term (`reduced = ~ genotype + treatment`); simple effects via combined coefficients or `results(dds, contrast=list(...))` | Interaction term tests whether the treatment effect differs by genotype; main-effect coefficients alone cannot answer that without the interaction term |
| Salmon/kallisto transcript-level quantification | `tximport(files, type="salmon"/"kallisto", countsFromAbundance="no")` -> `DESeqDataSetFromTximport()` (or edgeR/limma with the matched length-offset workflow) | Preserves the transcript-length offset per sample, correcting for differential isoform usage that a plain summed count matrix would miss |
| Paired tumor/normal or pre/post same subject | `design = ~ patient + condition` (DESeq2/limma); for limma with technical pseudo-replication also consider `duplicateCorrelation()` | Blocks on subject so within-patient correlation is modeled, not treated as independent noise |
| Multi-factor with batch + condition, batch known | `design = ~ batch + condition`; report `condition` contrast; use `limma::removeBatchEffect()` output only for PCA/heatmap, never for testing | Keeps batch variance in the model for valid inference while still allowing a visually batch-cleaned plot (see pitfall on ComBat-before-testing) |
| Small cohort, unmeasured/unknown confounders suspected | Estimate surrogate variables with SVA (`sva::svaseq()`) or RUVSeq and add them as covariates in the design formula, not as a pre-processing "cleaning" step on the counts | Adds unmeasured structure as design covariates so it's properly accounted for in the same GLM, avoiding the double-dipping problem of adjusting-then-testing |
| GSEA-ready ranking | Rank the full (unfiltered) gene table by shrunken log2FC or by `sign(log2FC) * -log10(pvalue)`/Wald statistic; never rank only the significant subset | GSEA requires a complete, monotonic ranking of all tested genes, not a truncated significant-only list |
| Cross-sample PCA / sample correlation / heatmap input | Variance-stabilizing transform (`vst()`) or `rlog()` counts (DESeq2), or voom's logCPM for limma pipelines — not raw counts, not naive `log2(count+1)` | VST/rlog stabilize variance across the mean-count range so PCA/clustering distances aren't dominated by highly expressed genes |
| Transcriptome-wide shift expected (viral shutoff, global stress response, prokaryotic stringent response) | Do not trust default median-of-ratios/TMM size factors; use spike-in-based normalization (ERCC) or an independently validated scaling factor, and state the assumption violation explicitly | Standard normalization assumes most genes are non-DE; that assumption breaks when the majority of the transcriptome shifts in one direction |
| Already-normalized/microarray-style continuous data (TPM, FPKM, intensities) | `limma::lmFit()` + `eBayes()` on log-scale values; do not run DESeq2/edgeR | DESeq2/edgeR require raw integer counts for their negative-binomial likelihood; continuous/pre-normalized data violates that model |
| n=1-2 per group ("pilot"/no true replication) | Report exploratory fold-change/ranking only; explicitly withhold formal significance claims; if forced to run DESeq2 with no replicates, flag the Poisson-only assumption | Dispersion cannot be reliably estimated from so few replicates; formal p-values would overstate confidence |
| Very large cohort (n>>10/group), public data re-analysis | Keep the effect-size threshold (e.g. abs(log2FC) >= 1) as a hard filter alongside padj, and consider reporting effect-size-ranked top hits rather than the full padj<0.05 list | At large n, padj<0.05 alone captures many biologically trivial effect sizes |
