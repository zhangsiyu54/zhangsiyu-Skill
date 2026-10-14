# ATAC Input Contract

Accepted inputs:
- peak count matrix: first column `peak`, remaining columns sample IDs.
- sample metadata with `sample` and group/covariate columns.
- differential peak table with `peak`, `log2FC`, `padj`.
- motif enrichment table with `motif`, `tf`, `pvalue` or `padj`.
- Signac/Seurat object, ArchR project, fragments, peak BED and genome annotation.

Required parameters:
- genome build, e.g. hg38/mm10.
- group column and contrast.
- peak annotation database or TxDb.
- whether analysis is bulk ATAC or scATAC.

Stop and ask if genome build, peak IDs, or group labels are missing.
