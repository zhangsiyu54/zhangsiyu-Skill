# bulkRNA Advanced Modules

Load this file for non-trivial RNA-seq projects.

## Input Type and Model Choice

- Raw gene-level counts: prefer DESeq2, edgeR, or limma-voom.
- Transcript-level quantification: summarize to gene level before DE, e.g. tximport-style workflow.
- Already normalized matrix: do not run DESeq2; use limma or exploratory plotting and state limitations.
- Small n or no biological replicates: do not claim reliable differential expression.

## Multi-factor Designs

- Record design formula, e.g. `~ batch + condition`, `~ patient + condition`, `~ time + treatment`.
- Confirm reference level, contrast, paired sample column and batch/covariate columns.
- Deliver model matrix/contrast notes and a design sanity table.

## QC and Exploratory Figures

- Required: library size, detected genes, PCA by group/batch, sample correlation heatmap.
- Optional: density/boxplot of normalized expression, Cook's distance/outlier notes, replicate consistency.
- Flag sample swaps, batch-dominated PCA, outlier samples and tiny group sizes.

## Differential Expression

- Use shrinkage or clear reporting for unstable logFC when possible.
- Deliver full DEG table, significant DEG table, up/down gene lists, volcano, MA plot, top gene heatmap.
- Keep thresholds explicit: padj, log2FC, baseMean/filtering.

## Downstream Handoff

- For GO/KEGG/GSEA, invoke or hand off to `enrichment-analysis`.
- For GSVA/WGCNA/module-trait correlation, invoke or hand off to `network-pathway-analysis`.
- For public cohort validation, invoke or hand off to `public-db-survival`.
