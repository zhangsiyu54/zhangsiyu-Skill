# Proteomics Metabolomics Deliverables

Minimum output structure:
- `output/delivery.pptx` when the customer requests a presentation-ready handoff; include project overview, input/QC, key results, figures, methods, limitations, and file index.
- `output/report.md`
- `output/methods.md`
- `output/summary.tsv`
- `output/results/differential_features.tsv`
- `output/results/missingness_summary.tsv`
- `output/qc/sample_qc.tsv` and `output/qc/feature_filtering.tsv` when the R route (`run_omics_templates.R`) is used
- `output/results/normalized_matrix.tsv` when generated
- `output/results/pathway_summary.tsv` when pathway interpretation is run
- `output/figures/pca.png` and `.pdf`
- `output/figures/missingness.png` and `.pdf` (Python route) or `output/figures/missingness.png`/`.pdf` (R route, per-sample missingness bar)
- `output/figures/sample_correlation.png` and `.pdf` when the R route is used
- `output/figures/volcano.png` and `.pdf`
- `output/figures/top_feature_heatmap.png` and `.pdf`
- `output/run/qc.json`
- `output/run/environment.json`
- `output/run/methods.json`
- `output/run/artifacts.json`

Report sections:
- input matrix summary
- normalization and missingness
- differential results
- pathway/multi-omics interpretation
- limitations
