---
name: network-pathway-analysis
description: Portable co-expression network and pathway activity workflow for Claude Code (and other Agent-Skills-compatible coding agents such as Codex CLI). Use when the user provides bulk RNA-seq, normalized expression, proteomics, metabolomics, single-cell pseudobulk, phenotype metadata, pathway gene sets, module assignments, WGCNA tables, or GSVA/ssGSEA scores and asks for WGCNA, module-trait correlation, hub genes, GSVA, ssGSEA, pathway activity heatmaps, module enrichment, pathway boxplots, correlation networks, or a reproducible network/pathway deliverable package.
---

# Network Pathway Analysis

## Core Workflow

Before choosing scripts, read `references/paper-practice.md`, `references/cns-method-guide.md`, and `references/top-references.md`. Follow the closest CNS-style workflow, including input eligibility, QC order, statistical model, parameter reporting, classic visualizations, interpretation boundaries, and literature-backed methods wording.
1. Identify input type: expression matrix, normalized matrix, module table, pathway score matrix, or existing WGCNA/GSVA results.
2. Confirm sample metadata, phenotype columns, organism, gene ID type, and whether expression has been normalized/log-transformed.
3. Read `references/advanced-modules.md` when WGCNA, GSVA, ssGSEA, module enrichment, or hub gene interpretation is requested.
4. Use `scripts/quickplot_network_pathway.py` for lightweight module-trait and pathway score figures from existing tables.
5. Use `scripts/run_wgcna_gsva.R` when R/WGCNA/GSVA dependencies and enough samples are available.
6. Export every figure as PNG and PDF and keep all outputs under `output/`.
7. Do not infer causality from correlation networks. State sample size and correlation method.
8. Before finalizing conclusions, check `references/pitfalls.md` for known WGCNA/hub-gene/GSVA failure modes and address any that apply.

## Visual Standard

Read `references/visual-style.md` before plotting. Use centered titles, readable axis labels, blue-white-red correlation heatmaps, and restrained pathway colors.

## Scripts

- `scripts/make_demo_data.py`: create demo module-trait, pathway score, and metadata tables.
- `scripts/quickplot_network_pathway.py`: generate module-trait heatmap (BH-FDR significance stars), pathway activity heatmap/barplot, and (with `--metadata`/`--group-col`) a grouped pathway boxplot plus `pathway_group_summary.tsv`, as PNG/PDF from existing tables. `--padj-threshold` (default 0.05) sets the BH-FDR cutoff used for heatmap stars, the pathway group-comparison flagging, and report.md's significant-cell counts/text.
- `scripts/run_wgcna_gsva.R`: formal WGCNA (soft-power selection with scale-free-fit check, block-wise module detection, module eigengenes, module-trait correlation with full-matrix BH-FDR, kME/GS-based hub genes) + GSVA/ssGSEA pathway scoring pipeline when R dependencies and enough samples are available; falls back to an exploratory correlation-only mode below the sample-size threshold (see `references/pitfalls.md` #2). Optional trailing positional args after `min_samples` set `padj_threshold` (default 0.05, same role as the Python flag) and `min_module_size` (default 30, passed to `blockwiseModules()`).
- `scripts/check_environment.py`: inspect R/Python availability and write `output/run/environment.json`.
- `scripts/write_methods_doc.py`: write `output/methods.md` and `output/run/methods.json`.
- `scripts/write_delivery_pptx.py`: assemble `output/delivery.pptx` from the already-written report.md/methods.md/figures/summary.tsv (run last, after the scripts above).

## Environment and Methods Rule

Before analysis, run `scripts/check_environment.py --out output` (see `references/r-dependencies.md` for the exact R package list and install command). Use R formal workflows when R and required packages are available; otherwise use Python quickplot scripts. If neither route is available, write setup instructions before analysis. After every analysis, run `scripts/write_methods_doc.py` or write an equivalent `output/methods.md` documenting inputs, parameters, route, thresholds, packages/databases, and limitations. After the core deliverables (report.md, methods.md) exist, ALWAYS also run `scripts/write_delivery_pptx.py --out output --skill-name network-pathway-analysis` as a standard final step -- do not wait for the user to ask for a presentation deck; every completed analysis ships with one automatically, in bilingual (Chinese/English) form.

## Required Inputs

See `references/input-contract.md`.

## Deliverable Checklist

See `references/deliverables.md`.
