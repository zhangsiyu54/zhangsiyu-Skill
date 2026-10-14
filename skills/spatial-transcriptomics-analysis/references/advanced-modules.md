# Spatial Transcriptomics Advanced Modules

## QC

- Report spots, median counts/features, mitochondrial percentage if available, tissue spots and low-quality regions.
- Use spatial QC maps to identify edge artifacts or tissue folds.

## Spatial Feature Maps

- Preserve spatial orientation.
- Use continuous color scales for gene expression and categorical palettes for domains.
- Avoid comparing raw expression across sections without normalization notes.

## Domains and Neighborhoods

- Spatial domains can come from clustering, pathology regions or user labels.
- Neighborhood enrichment needs domain/cell type labels and spatial neighbor graph.
- Report method, neighbor radius/k and sample section.

## Deconvolution

- Use SPOTlight, RCTD, cell2location or provided deconvolution table depending on environment.
- Deliver cell type proportion maps and spot-level proportion table.
- State reference dataset and species/tissue match.

## Integration

- For reference single-cell objects, combine with `singlecell-analysis`.
- For pathway maps, combine with `network-pathway-analysis`.
