---
name: spatial-transcriptomics-analysis
description: Portable spatial transcriptomics workflow for Claude Code (and other Agent-Skills-compatible coding agents such as Codex CLI). Use when the user provides 10x Visium output, spatial spot coordinates, expression matrix, Seurat spatial object, AnnData/Squidpy object, tissue image, cell type deconvolution table, spatial domains, marker genes, or target pathways and asks for spatial QC, spot-level feature maps, spatial domains, spatially variable genes, neighborhood enrichment, ligand-receptor spatial context, deconvolution heatmaps, region comparison, or spatial transcriptomics deliverables.
---

# Spatial Transcriptomics Analysis

## Core Workflow

Before choosing scripts, read `references/paper-practice.md`, `references/cns-method-guide.md`, and `references/top-references.md`. Follow the closest CNS-style workflow, including input eligibility, QC order, statistical model, parameter reporting, classic visualizations, interpretation boundaries, and literature-backed methods wording.
1. Identify input type: 10x Visium folder, Seurat spatial RDS, h5ad/Squidpy object, spot table, or exported coordinates and expression.
2. Confirm spot IDs, coordinates, sample/section labels, target genes, domain labels and whether a histology image is available.
3. Read `references/advanced-modules.md` for Visium preprocessing, spatial domains, deconvolution, spatially variable genes or neighborhood analysis.
4. Use `scripts/quickplot_spatial.py` for coordinate-based feature and domain figures from exported tables.
5. Use `scripts/run_spatial_seurat.R` when Seurat spatial dependencies and Visium objects are available.
6. Export every figure as PNG and PDF.
7. Do not infer cell types from spots without reference/deconvolution evidence.
8. Before finalizing conclusions, check `references/pitfalls.md` for common spatial analysis failure modes (spot-mixture cell-type calls, expression-only domain clustering vs. over-smoothing, depth-confounded SVGs, image/coordinate misalignment, over-interpreted ligand-receptor co-localization, batch effects, deconvolution reference mismatch, post-hoc region testing).

## Visual Standard

Read `references/visual-style.md` before plotting. Spatial maps must preserve coordinate orientation unless the user requests a flip/rotation.

## Scripts

- `scripts/make_demo_data.py`: create demo spot coordinate/expression/domain tables.
- `scripts/quickplot_spatial.py`: generate spatial feature (one or more genes via comma-separated `--gene`), domain and QC maps as PNG/PDF, spot QC and domain-summary result tables, and an optional dominant-cell-type map/table from a `--deconvolution` proportions CSV.
- `scripts/run_spatial_seurat.R`: formal Seurat spatial pipeline (accepts a Visium directory via `Load10X_Spatial` or a saved Seurat `.rds`) — SCTransform normalization, PCA/clustering/UMAP, `SpatialFeaturePlot`/`SpatialDimPlot` maps, Moran's I spatially variable gene ranking, and result tables — when Seurat/ggplot2 dependencies are available.
- `scripts/check_environment.py`: inspect R/Python availability and write `output/run/environment.json`.
- `scripts/write_methods_doc.py`: write `output/methods.md` and `output/run/methods.json`.
- `scripts/write_delivery_pptx.py`: assemble `output/delivery.pptx` from the already-written report.md/methods.md/figures/summary.tsv (run last, after the scripts above).

## Environment and Methods Rule

Before analysis, run `scripts/check_environment.py --out output` (see `references/r-dependencies.md` for the exact R package list and install command). Use R formal workflows when R and required packages are available; otherwise use Python quickplot scripts. If neither route is available, write setup instructions before analysis. After every analysis, run `scripts/write_methods_doc.py` or write an equivalent `output/methods.md` documenting inputs, parameters, route, thresholds, packages/databases, and limitations. After the core deliverables (report.md, methods.md) exist, ALWAYS also run `scripts/write_delivery_pptx.py --out output --skill-name spatial-transcriptomics-analysis` as a standard final step -- do not wait for the user to ask for a presentation deck; every completed analysis ships with one automatically, in bilingual (Chinese/English) form.

## Required Inputs

See `references/input-contract.md`.

## Deliverable Checklist

See `references/deliverables.md`.
