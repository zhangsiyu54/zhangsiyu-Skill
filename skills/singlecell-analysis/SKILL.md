---
name: singlecell-analysis
description: Portable deep single-cell RNA-seq workflow for Claude Code (and other Agent-Skills-compatible coding agents such as Codex CLI). Use when the user provides Seurat RDS, h5ad, 10x matrix, metadata, marker tables, embeddings, clusters, cell types, conditions, or target genes and asks for QC/filtering, doublet or ambient RNA review, normalization, integration or batch correction, PCA/UMAP/clustering, marker genes, annotation, FeaturePlot/DotPlot/violin plots, cell type proportions, pseudobulk differential expression, trajectory or pseudotime, module scores, cell-cell communication, or a complete single-cell deliverable package.
---

# Single-cell Analysis

## Core Workflow

Before choosing scripts, read `references/paper-practice.md`, `references/cns-method-guide.md`, and `references/top-references.md`. Follow the closest CNS-style workflow, including input eligibility, QC order, statistical model, parameter reporting, classic visualizations, interpretation boundaries, and literature-backed methods wording.
1. Inspect input format first: Seurat RDS, h5ad, 10x folder, matrix+metadata, or exported embedding/marker tables.
2. Validate cell IDs, metadata columns, cluster/celltype fields, group labels, and requested genes.
3. For full workflows, prefer `scripts/run_seurat.R` when Seurat is available.
4. Read `references/advanced-modules.md` when the user requests integration, annotation, pseudobulk DE, trajectory, cell-cell communication, or module scoring.
5. For exported tables or lightweight demos, use `scripts/quickplot_singlecell.py` to generate UMAP, marker dotplot, and composition plots.
6. Use the optional R templates for advanced analyses only when dependencies and metadata are present.
7. Keep raw object handling read-only; write all generated files under `output/`.
8. Do not guess cell type names. Use marker evidence or label ambiguous clusters as `needs_annotation`.
9. Write `output/questions.md` if the object cannot be loaded or required metadata is missing.
10. Before finalizing conclusions, read `references/pitfalls.md` and check the analysis against each failure mode (pseudoreplication, over-regression, over-integration, ambient RNA, single-marker annotation, arbitrary resolution/PCs, doublet order-of-operations, and more).

## Visual Standard

Read `references/visual-style.md` before final plots. Use calm categorical palettes, small point sizes for UMAP, fixed legend order, and non-overlapping labels.

## Scripts

- `scripts/make_demo_data.py`: create metadata, UMAP coordinates, and marker expression demo tables.
- `scripts/quickplot_singlecell.py`: UMAP, composition, and dotplot PNG/PDF generator for exported tables.
- `scripts/quickplot_communication.py`: ligand-receptor network and dotplot PNG/PDF generator from edge tables.
- `scripts/run_seurat.R`: formal Seurat workflow template for R users.
- `scripts/run_pseudobulk_de.R`: pseudobulk differential expression template for cell type by condition contrasts; runs per named cell type found in `celltype_col` (no hardcoded cell type) with configurable `padj_threshold` (default 0.05) and `log2fc_threshold` (default 1) CLI arguments.
- `scripts/run_cellchat.R`: cell-cell communication template using CellChat when available.
- `scripts/run_trajectory_monocle3.R`: trajectory/pseudotime template using monocle3 when available.
- `scripts/check_environment.py`: inspect R/Python availability and write `output/run/environment.json`.
- `scripts/write_methods_doc.py`: write `output/methods.md` and `output/run/methods.json`.
- `scripts/write_delivery_pptx.py`: assemble `output/delivery.pptx` from the already-written report.md/methods.md/figures/summary.tsv (run last, after the scripts above).

## Environment and Methods Rule

Before analysis, run `scripts/check_environment.py --out output` (see `references/r-dependencies.md` for the exact R package list and install command). Use R formal workflows when R and required packages are available; otherwise use Python quickplot scripts. If neither route is available, write setup instructions before analysis. After every analysis, run `scripts/write_methods_doc.py` or write an equivalent `output/methods.md` documenting inputs, parameters, route, thresholds, packages/databases, and limitations. After the core deliverables (report.md, methods.md) exist, ALWAYS also run `scripts/write_delivery_pptx.py --out output --skill-name singlecell-analysis` as a standard final step -- do not wait for the user to ask for a presentation deck; every completed analysis ships with one automatically, in bilingual (Chinese/English) form.

## Required Inputs

See `references/input-contract.md` for object/table formats and required metadata.

## Deliverable Checklist

See `references/deliverables.md` before handoff.
