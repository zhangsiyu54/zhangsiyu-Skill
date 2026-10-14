# Paper Practice

Follow Seurat/Scanpy/CellChat-style single-cell papers: quality control first, biological annotation second, advanced interpretation only after metadata and evidence are present.

## Reference Sources

- Seurat GitHub: https://github.com/satijalab/seurat
- Seurat integration Cell paper: https://doi.org/10.1016/j.cell.2019.05.031
- Scanpy paper: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1382-0
- CellChat GitHub: https://github.com/sqjin/CellChat
- CellChat paper: https://www.nature.com/articles/s41467-021-21246-9

## Paper-style Flow

1. Inspect object type and metadata; confirm cell IDs, sample IDs, condition, cluster, and annotation fields.
2. Run or review QC: detected genes, UMIs, mitochondrial ratio, doublet/ambient limitations.
3. Normalize, reduce dimensions, cluster, and annotate with marker evidence; do not invent cell type names.
4. Add composition and pseudobulk differential modules when group metadata exist.
5. Add trajectory, module score, or cell-cell communication only when required fields are present.
6. Write methods with filtering thresholds, integration/batch method, clustering resolution, annotation evidence, and databases used.

## Visualization Pattern

- UMAP: small semi-transparent points, stable cluster/celltype colors, legends outside.
- Marker DotPlot: color = average expression, size = percent expressed; reserve wide right margin for colorbar/size legend.
- Cell communication: use circle/network for overview and bubble/dotplot for ligand-receptor pairs; keep labels sparse.
