# Network Pathway Input Contract

Accepted inputs:
- expression matrix: first column `gene`, remaining columns sample IDs.
- metadata table: one row per sample with `sample` and phenotype columns.
- pathway score table: first column `pathway`, remaining columns sample IDs.
- module-trait table: columns `module`, `trait`, `correlation`, `pvalue` or `padj`.
- gene set GMT for GSVA/ssGSEA when formal scoring is requested.

Required parameters:
- organism and gene ID type.
- phenotype columns for module-trait correlation.
- whether expression values are raw counts, normalized counts, log expression, protein intensity, or metabolite abundance.

Stop and ask if:
- sample IDs do not match metadata.
- sample size is too small for WGCNA.
- gene ID type or gene set database is missing.
