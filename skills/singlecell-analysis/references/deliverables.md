# Single-cell Deliverables

Minimum output structure:
- `output/delivery.pptx` when the customer requests a presentation-ready handoff; include project overview, input/QC, key results, figures, methods, limitations, and file index.
- `output/report.md`
- `output/methods.md`
- `output/summary.tsv`
- `output/results/marker_summary.tsv`
- `output/results/celltype_proportion.tsv`
- `output/results/annotation_evidence.tsv` when annotation is requested
- `output/results/module_score_summary.tsv` and `output/results/module_score_percell.tsv` (R route) when module scoring is requested
- `output/results/pseudobulk_summary.tsv` and per-celltype DEG tables when condition comparison is requested
- `output/results/trajectory/pseudotime.tsv` when trajectory is requested
- `output/results/trajectory/pseudotime_genes.tsv` when trajectory is run (pseudotime-associated/branch genes)
- `output/results/cellchat/ligand_receptor_table.tsv` when communication analysis is requested
- `output/results/cellchat/pathway_summary.tsv` when communication analysis is run
- `output/results/cellchat/group_comparison_summary.tsv` when communication analysis is run with a group/condition column
- `output/results/ligand_receptor_table.tsv` when quick communication plotting is run
- `output/figures/qc_violin.png` and `.pdf`
- `output/figures/umap_clusters.png` and `.pdf`
- `output/figures/umap_celltype.png` and `.pdf` when annotation exists
- `output/figures/umap_sample.png` and `.pdf` when a sample column exists
- `output/figures/umap_condition.png` and `.pdf` when a condition/group column exists
- `output/figures/featureplot.png` and `.pdf` when target genes exist (Python quicklook also writes `featureplot_<gene>.png/.pdf` per gene)
- `output/figures/heatmap_markers.png` and `.pdf` (top marker heatmap)
- `output/figures/dotplot_markers.png` and `.pdf`
- `output/figures/celltype_proportion.png` and `.pdf`
- `output/figures/module_score_umap.png` and `.pdf`, `output/figures/module_score_violin.png` and `.pdf` when module scoring is requested
- `output/figures/pseudobulk/*_volcano.png` and `.pdf` when pseudobulk DE is run
- `output/figures/pseudobulk/*_heatmap.png` and `.pdf` when pseudobulk DE is run (top DE gene heatmap)
- `output/figures/trajectory/pseudotime.png` and `.pdf` when trajectory is run
- `output/figures/trajectory/pseudotime_by_celltype.png` and `.pdf` when trajectory is run (pseudotime distribution by group)
- `output/figures/cellchat/interaction_circle.png` and `.pdf` when communication analysis is run
- `output/figures/cellchat/interaction_bubble.png` and `.pdf` when communication analysis is run
- `output/figures/cellchat/interaction_comparison.png`, `pathway_rank_comparison.png`, and `interaction_diff_circle.png` (with `.pdf`) when communication analysis is run with a group/condition column
- `output/figures/communication_network.png` and `.pdf` when quick communication plotting is run
- `output/figures/ligand_receptor_dotplot.png` and `.pdf` when quick communication plotting is run
- `output/run/qc.json`
- `output/run/environment.json`
- `output/run/methods.json`
- `output/run/artifacts.json`

Report sections:
- input object/table summary
- filtering rules
- normalization and dimensional reduction
- cluster/annotation method
- marker and composition results
- pseudobulk, trajectory, module score, and communication results when requested
- limitations and missing metadata
