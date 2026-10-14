# Network Pathway Deliverables

Minimum output structure:
- `output/delivery.pptx` when the customer requests a presentation-ready handoff; include project overview, input/QC, key results, figures, methods, limitations, and file index.
- `output/report.md`
- `output/methods.md`
- `output/summary.tsv`
- `output/results/module_trait_correlation.tsv` (module/gene, trait, correlation, pvalue, padj across the full module x trait matrix)
- `output/results/gene_modules.tsv` when WGCNA is run (module assignment table: gene, module)
- `output/results/hub_genes.tsv` when WGCNA is run (kME, gene significance, hub flag and thresholds used)
- `output/results/pathway_scores.tsv` when GSVA/ssGSEA is run
- `output/results/pathway_group_summary.tsv` when groups are available (group means + Wilcoxon/Mann-Whitney or Kruskal-Wallis test, BH-adjusted)
- `output/figures/sample_clustering.png` and `.pdf` when the formal WGCNA route is run
- `output/figures/soft_threshold_selection.png` and `.pdf` when the formal WGCNA route is run
- `output/figures/module_dendrogram.png` and `.pdf` when the formal WGCNA route is run (skipped in exploratory small-sample mode)
- `output/figures/hub_gene_network_summary.png` and `.pdf` when the formal WGCNA route is run
- `output/figures/module_trait_heatmap.png` and `.pdf`
- `output/figures/pathway_activity_heatmap.png` and `.pdf`
- `output/figures/pathway_activity_barplot.png` and `.pdf`
- `output/figures/pathway_activity_boxplot.png` and `.pdf` when groups are available
- `output/run/qc.json`
- `output/run/environment.json`
- `output/run/methods.json`
- `output/run/artifacts.json`

Report sections:
- input matrix and metadata summary
- WGCNA/GSVA parameters
- module-trait results
- pathway activity results
- hub genes and limitations
