# Immune Repertoire Deliverables

Minimum output structure:
- `output/delivery.pptx` when the customer requests a presentation-ready handoff; include project overview, input/QC, key results, figures, methods, limitations, and file index.
- `output/report.md`
- `output/methods.md`
- `output/summary.tsv`
- `output/results/clonotype_summary.tsv`
- `output/results/diversity_summary.tsv` (per-sample Shannon raw + rarefied at common depth, Gini-Simpson, inverse Simpson, Chao1, clonality, total count/depth, observed clonotype richness -- see the depth caveat below)
- `output/results/diversity_caveat.txt` (non-silent depth/cell-number-sensitivity caveat, see `references/pitfalls.md` #1/#4)
- `output/results/overlap_matrix_morisita_horn.tsv` and `output/results/overlap_matrix_jaccard.tsv` (pairwise sample x sample clonotype overlap; only produced when 2+ samples are present)
- `output/figures/clonotype_expansion.png` and `.pdf` (clone-size-class stacked bar; never labels individual clones)
- `output/figures/top_clonotype_frequency.png` and `.pdf` (top-N clonotypes only, "Other" bucket for the rest; also serves as a lightweight cross-sample clone-tracking view)
- `output/figures/diversity_barplot.png` and `.pdf` (raw vs. rarefied Shannon plus a sample-depth panel, so depth and diversity are always shown together)
- `output/figures/v_gene_usage.png` and `.pdf`
- `output/figures/j_gene_usage.png` and `.pdf` (only when a `j_gene` column is present)
- `output/figures/cdr3_length_distribution.png` and `.pdf` (only when a `cdr3` column is present)
- `output/figures/clonotype_overlap_heatmap.png` and `.pdf` (Morisita-Horn; only when 2+ samples are present)
- `output/run/environment.json`
- `output/run/methods.json`
- `output/run/qc.json`
