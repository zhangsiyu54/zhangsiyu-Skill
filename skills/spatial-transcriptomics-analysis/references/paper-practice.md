# Paper Practice

Follow Squidpy/cell2location/Seurat spatial papers: preserve tissue geometry and link spatial patterns to cell type or neighborhood evidence.

## Reference Sources

- Squidpy GitHub: https://github.com/scverse/squidpy
- Squidpy Nature Methods: https://www.nature.com/articles/s41592-021-01358-2
- cell2location GitHub: https://github.com/BayraktarLab/cell2location
- cell2location Nature Biotechnology: https://www.nature.com/articles/s41587-021-01139-4
- Seurat spatial vignette: https://satijalab.org/seurat/articles/spatial_vignette

## Paper-style Flow

1. Confirm spot/barcode IDs, coordinates, image availability, sample/section labels, and target genes/domains.
2. Run or review spatial QC and spot filtering before feature maps.
3. Generate feature maps, domain maps, spatially variable genes, and region comparisons.
4. Add deconvolution or cell type abundance maps only with a reference and method.
5. Add neighborhood enrichment or spatial ligand-receptor context when requested and supported.

## Visualization Pattern

- Preserve coordinate orientation and equal aspect ratio.
- Use light-to-saturated feature scales and categorical domain colors.
- Keep legends outside tissue panels and avoid text over tissue unless labels are sparse.
