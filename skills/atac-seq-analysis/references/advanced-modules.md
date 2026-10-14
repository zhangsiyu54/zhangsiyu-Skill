# ATAC Advanced Modules

## Bulk ATAC

- QC: library size, FRiP, TSS enrichment if available, replicate correlation, PCA.
- Differential accessibility: DESeq2/edgeR/limma-voom on peak counts.
- Peak annotation: nearest gene, promoter/enhancer category, genomic distribution.
- Motif analysis: enriched motifs among opening/closing peaks; state motif database.

## scATAC

- QC: fragments per cell, TSS enrichment, nucleosome signal, blacklist ratio.
- Dimensionality reduction: TF-IDF, LSI, UMAP, clustering.
- Gene activity/gene score: map accessibility to genes with method note.
- Peak-to-gene links: require sufficient cells and clear correlation method.
- Multiome integration: combine with `singlecell-analysis` when RNA modality exists.

## Interpretation

- Accessibility is not expression. Phrase results as "open/closed chromatin" or "motif enrichment", not direct TF activity unless using additional evidence.
