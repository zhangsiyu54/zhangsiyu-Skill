# Spatial Transcriptomics Input Contract

Accepted inputs:
- 10x Visium output folder.
- Seurat spatial `.rds` or AnnData `.h5ad`.
- exported spot table with `spot`, `x`, `y`, optional `domain`, `sample`, `nCount`, `nFeature`.
- expression table with `spot` and target gene columns.
- deconvolution table with spot-by-celltype proportions.

Required parameters:
- target genes or features for spatial maps.
- coordinate columns.
- sample/section column when multiple sections exist.
- domain/cell type columns when comparing regions.

Stop and ask if:
- spot IDs do not match between coordinates and expression.
- coordinate orientation is ambiguous and image overlay is required.
- deconvolution is requested but no reference or deconvolution table is provided.
