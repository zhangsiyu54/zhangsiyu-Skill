# Proteomics Metabolomics Input Contract

Accepted inputs:
- abundance matrix: first column `feature`, remaining columns sample IDs.
- metadata table with `sample` and group/covariate columns.
- feature annotation table: `feature`, `gene`, `protein`, `metabolite`, `pathway`, `id`.
- differential table with `feature`, `log2FC`, `padj`.

Required parameters:
- omics type: proteomics, phosphoproteomics, metabolomics, lipidomics, or mixed.
- group column and case/control labels.
- normalization status and missing value handling.
- feature ID type for pathway interpretation.

Stop and ask if sample IDs, contrast or missing-value rule is missing.
