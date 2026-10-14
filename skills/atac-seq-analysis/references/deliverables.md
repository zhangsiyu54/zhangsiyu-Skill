# ATAC Deliverables

Minimum output structure:
- `output/delivery.pptx` when the customer requests a presentation-ready handoff; include project overview, input/QC, key results, figures, methods, limitations, and file index.
- `output/report.md`
- `output/methods.md`
- `output/summary.tsv`
- `output/results/differential_peaks.tsv`
- `output/results/motif_enrichment.tsv`
- `output/results/peak_annotation.tsv` when annotation is run
- `output/results/scatac_cluster_markers.tsv` when scATAC is run
- `output/figures/differential_peak_volcano.png` and `.pdf`
- `output/figures/motif_enrichment_barplot.png` and `.pdf`
- `output/figures/peak_heatmap.png` and `.pdf`
- `output/figures/pca.png` and `.pdf` when the formal `run_atac_templates.R` route is used
- `output/figures/sample_correlation_heatmap.png` and `.pdf` when the formal `run_atac_templates.R` route is used
- `output/figures/qc_overview.png` and `.pdf` (library depth and peaks-detected per sample) when `scripts/quickplot_atac.py` is used
- `output/results/normalized_accessibility.tsv` when the formal `run_atac_templates.R` route is used
- `output/run/qc.json`
- `output/run/environment.json`
- `output/run/methods.json`
- `output/run/artifacts.json`
