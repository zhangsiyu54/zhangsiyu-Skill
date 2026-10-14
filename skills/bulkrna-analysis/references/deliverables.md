# bulkRNA Deliverables

Minimum output structure:
- `output/delivery.pptx` when the customer requests a presentation-ready handoff; include project overview, input/QC, key results, figures, methods, limitations, and file index.
- `output/report.md`
- `output/methods.md`
- `output/summary.tsv`
- `output/results/differential_expression.tsv`
- `output/results/limma_voom_differential_expression.tsv` when limma-voom is used
- `output/results/differential_expression.tsv`, `output/results/normalized_expression_used.tsv`, `output/results/up_genes.txt`, `output/results/down_genes.txt`, `output/results/gsea_ranked_list.rnk` when `run_limma_trend.R` (already-normalized/log-scale expression route) is used
- `output/results/normalized_counts.tsv` when generated
- `output/qc/sample_qc.tsv` (sample, group, library_size, detected_genes)
- `output/qc/design_matrix.tsv` and `output/qc/library_size.tsv` when multi-factor modeling is used
- `output/figures/volcano.png` and `output/figures/volcano.pdf`
- `output/figures/heatmap_top_degs.png` and `output/figures/heatmap_top_degs.pdf`
- `output/figures/sample_correlation_heatmap.png` and `.pdf` (Pearson sample correlation, produced by both `run_deseq2.R` and `scripts/quickplot_bulk.py`)
- `output/figures/pca.png` and `.pdf` (PC1 vs PC2 colored by group, produced by both `run_deseq2.R` and `scripts/quickplot_bulk.py`)
- `output/figures/qc_summary.png` and `.pdf` (library size and detected genes per sample, produced by both `run_deseq2.R` and `scripts/quickplot_bulk.py`)
- `output/figures/limma_voom_volcano.png` and `.pdf` when limma-voom is used
- `output/figures/limma_voom_heatmap.png` and `.pdf` when limma-voom is used
- `output/run/qc.json`
- `output/run/environment.json`
- `output/run/methods.json`
- `output/run/artifacts.json`

Report sections:
- input summary
- contrast and thresholds
- QC observations
- main DEG counts
- methods and package versions
- limitations
- file list
