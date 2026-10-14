# Microbiome Deliverables

Minimum output structure:
- `output/delivery.pptx` when the customer requests a presentation-ready handoff; include project overview, input/QC, key results, figures, methods, limitations, and file index.
- `output/report.md`
- `output/methods.md`
- `output/summary.tsv`
- `output/results/alpha_diversity.tsv` (Observed/Shannon/Simpson per sample)
- `output/results/alpha_diversity_stats.tsv` (per-metric group test; formal R route only)
- `output/results/taxa_relative_abundance.tsv`
- `output/results/differential_taxa.tsv` when differential abundance is run (ANCOM-BC2, formal R route)
- `output/results/permanova.tsv` when beta diversity testing is run (PERMANOVA + betadisper/dispersion row; both routes run this by default)
- `output/figures/read_depth.png` and `.pdf`
- `output/figures/alpha_diversity.png` and `.pdf`
- `output/figures/taxa_composition.png` and `.pdf`
- `output/figures/ordination.png` and `.pdf` (PCoA on Bray-Curtis distance)
- `output/figures/differential_taxa.png` and `.pdf` when differential abundance is run (formal R route)
- `output/run/qc.json`
- `output/run/environment.json`
- `output/run/methods.json`
- `output/run/artifacts.json`

Report sections:
- input table summary
- taxonomy rank and filtering
- diversity results
- differential abundance method
- compositional limitations
