# Enrichment Deliverables

Minimum output structure:
- `output/delivery.pptx` when the customer requests a presentation-ready handoff; include project overview, input/QC, key results, figures, methods, limitations, and file index.
- `output/report.md` (from `scripts/quickplot_enrichment.py` or `scripts/run_enrichment.R`)
- `output/report_gsea.md` (from `scripts/run_gsea_fgsea.R` or `scripts/quickplot_gsea.py` — kept separate from `report.md` so an ORA run and a GSEA run into the same `output/` directory don't overwrite each other's report)
- `output/methods.md`
- `output/summary.tsv`
- `output/results/go_enrichment.tsv` (full tested-term table, not only significant terms; includes a `significant` column) or `output/results/kegg_enrichment.tsv`, from `scripts/run_enrichment.R`
- `output/results/go_enrichment_simplified.tsv` (or `kegg_enrichment_simplified.tsv`) when GO redundancy reduction via `clusterProfiler::simplify()` was applied and GOSemSim was available
- `output/results/gsea_results.tsv` when ranked input exists
- `output/results/fgsea_results.tsv` (full tested-pathway table with NES/padj) when ranked GMT GSEA is run
- `output/results/gsea_leading_edge.tsv` (long-format pathway/gene leading-edge table) when ranked GMT GSEA (`scripts/run_gsea_fgsea.R`) is run
- `output/results/unmapped_genes.tsv` when ID mapping is performed (e.g. SYMBOL→ENTREZID for KEGG ORA)
- `output/figures/enrichment_bubble.png` and `output/figures/enrichment_bubble.pdf`
- `output/figures/enrichment_barplot.png` and `output/figures/enrichment_barplot.pdf`
- `output/figures/enrichment_cnetplot.png` and `.pdf` (gene-concept network) when a `geneID` column is available and >=3 terms are significant
- `output/figures/go_enrichment_dotplot.png`/`.pdf` (or `kegg_enrichment_dotplot`) and `output/figures/go_enrichment_barplot.png`/`.pdf` (or `kegg_enrichment_barplot`) when `scripts/run_enrichment.R` is run
- `output/figures/gsea_dotplot.png` and `.pdf` (NES summary dotplot) when ranked GSEA (`scripts/run_gsea_fgsea.R`) is run
- `output/figures/gsea_running_score_top.png` and `.pdf` (running-score plot for the top pathway) when ranked GSEA (`scripts/run_gsea_fgsea.R`) is run
- `output/figures/gsea_running_score.png` and `.pdf` when GSEA quicklook (`scripts/quickplot_gsea.py`) is run
- `output/figures/gsea_summary_dotplot.png` and `.pdf` when GSEA quicklook is run
- `output/run/qc.json`
- `output/run/environment.json`
- `output/run/methods.json`
- `output/run/artifacts.json`

Report sections:
- input gene list summary
- database scope
- ID mapping status
- top enriched terms
- plotting rules
- limitations
