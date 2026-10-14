# Single-cell Input Contract

Accepted inputs:
- Seurat `.rds`
- AnnData `.h5ad`
- 10x folder with `matrix.mtx`, `features.tsv`, `barcodes.tsv`
- exported tables: `cell_metadata.csv`, `umap.csv`, `marker_expression.csv`, `markers.csv`

Required metadata:
- cell ID column.
- cluster or cell type column.
- condition/group column when comparing groups.
- target genes for FeaturePlot, VlnPlot, DotPlot, or module score.

Stop and ask if:
- object cannot be loaded locally.
- cell IDs do not match between metadata and expression/embedding tables.
- marker genes or group labels are missing for requested plots.
- annotation is requested but no marker evidence is provided.
