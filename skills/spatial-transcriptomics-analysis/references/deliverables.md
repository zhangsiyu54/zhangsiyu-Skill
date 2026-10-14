# Spatial Transcriptomics Deliverables

Minimum output structure:
- `output/delivery.pptx` when the customer requests a presentation-ready handoff; include project overview, input/QC, key results, figures, methods, limitations, and file index.
- `output/report.md`
- `output/methods.md`
- `output/summary.tsv`
- `output/results/spot_qc.tsv`
- `output/results/spatial_feature_values.tsv`
- `output/results/domain_summary.tsv` when domain labels exist
- `output/results/deconvolution_proportions.tsv` when deconvolution exists
- `output/results/spatially_variable_genes.tsv` when the Seurat route computes spatially variable genes
- `output/figures/spatial_feature_map.png` and `.pdf` (plus `spatial_feature_map_<gene>.png`/`.pdf` for each additional requested gene beyond the first)
- `output/figures/spatial_domain_map.png` and `.pdf`
- `output/figures/spot_qc_map.png` and `.pdf`
- `output/figures/deconvolution_composition_map.png` and `.pdf` when deconvolution exists
- `output/figures/umap_clusters.png` and `.pdf` (optional additional diagnostic from the Seurat route)
- `output/run/qc.json`
- `output/run/environment.json`
- `output/run/methods.json`
- `output/run/artifacts.json`

Report sections:
- input object/table summary
- QC observations
- target gene/domain maps
- deconvolution or neighborhood result when available
- limitations
