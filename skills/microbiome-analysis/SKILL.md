---
name: microbiome-analysis
description: Portable microbiome analysis workflow for Claude Code (and other Agent-Skills-compatible coding agents such as Codex CLI). Use when the user provides 16S, ASV/OTU abundance table, taxonomy table, sample metadata, alpha diversity, beta diversity, PCoA coordinates, QIIME2/DADA2/phyloseq objects, differential taxa tables, or metagenome functional profiles and asks for microbiome QC, rarefaction review, alpha diversity, beta diversity, PCoA, PERMANOVA, taxa composition stacked bars, heatmaps, differential abundance, ANCOM-BC/DESeq2-style taxa analysis, functional profile plots, or a reproducible microbiome deliverable package.
---

# Microbiome Analysis

## Core Workflow

Before choosing scripts, read `references/paper-practice.md`, `references/cns-method-guide.md`, and `references/top-references.md`. Follow the closest CNS-style workflow, including input eligibility, QC order, statistical model, parameter reporting, classic visualizations, interpretation boundaries, and literature-backed methods wording.
1. Identify input type: raw reads, QIIME2 artifact, ASV/OTU table, taxonomy table, metadata, or existing diversity results.
2. Confirm sample IDs, group column, taxonomy rank, feature type, sequencing region, and whether data are counts or relative abundance.
3. Read `references/advanced-modules.md` for DADA2/QIIME2, compositional differential abundance, PERMANOVA, or functional metagenome interpretation.
4. Use `scripts/quickplot_microbiome.py` for read depth, alpha diversity, taxa composition, Bray-Curtis PCoA/PERMANOVA figures and tables from plain CSV inputs.
5. Use `scripts/run_phyloseq_ancombc.R` for the same figures plus ANCOM-BC2 differential abundance when phyloseq/ANCOMBC/mia/vegan dependencies are available.
6. Export every figure as PNG and PDF under `output/`.
7. Do not interpret differential taxa from relative abundance t tests as definitive; state compositional limitations.
8. Before finalizing conclusions, read `references/pitfalls.md` and check for rarefaction, compositionality, dispersion, batch/kitome, database-version, and method-choice failure modes that may apply.

## Visual Standard

Read `references/visual-style.md` before plotting. Use readable taxonomy labels, stable group colors and restrained stacked bars.

## Scripts

- `scripts/make_demo_data.py`: create demo abundance, taxonomy and metadata tables.
- `scripts/quickplot_microbiome.py`: dependency-light (matplotlib/numpy only) pipeline from abundance/metadata tables: read depth, alpha diversity (Observed/Shannon/Simpson boxplots), taxa composition (top taxa + Other), Bray-Curtis PCoA ordination, permutation-based PERMANOVA + betadisper-style dispersion check, and result tables/PNG/PDF figures under `output/`. Does not run differential abundance (relative-abundance data is compositional; use the R route for that).
- `scripts/run_phyloseq_ancombc.R`: formal phyloseq/ANCOM-BC2 workflow when Bioconductor dependencies (phyloseq, ANCOMBC, mia, vegan, ggplot2) are available: builds the phyloseq object, alpha diversity (Observed/Shannon/Simpson) with group tests, Bray-Curtis PCoA with `vegan::adonis2` PERMANOVA and `betadisper` dispersion check, taxa composition (top taxa + Other), and ANCOM-BC2 differential abundance (prevalence filter, structural-zero handling, BH-FDR) with result table and volcano-style plot.
- `scripts/check_environment.py`: inspect R/Python availability and write `output/run/environment.json`.
- `scripts/write_methods_doc.py`: write `output/methods.md` and `output/run/methods.json`.
- `scripts/write_delivery_pptx.py`: assemble `output/delivery.pptx` from the already-written report.md/methods.md/figures/summary.tsv (run last, after the scripts above).

## Environment and Methods Rule

Before analysis, run `scripts/check_environment.py --out output` (see `references/r-dependencies.md` for the exact R package list and install command). Use R formal workflows when R and required packages are available; otherwise use Python quickplot scripts. If neither route is available, write setup instructions before analysis. After every analysis, run `scripts/write_methods_doc.py` or write an equivalent `output/methods.md` documenting inputs, parameters, route, thresholds, packages/databases, and limitations. After the core deliverables (report.md, methods.md) exist, ALWAYS also run `scripts/write_delivery_pptx.py --out output --skill-name microbiome-analysis` as a standard final step -- do not wait for the user to ask for a presentation deck; every completed analysis ships with one automatically, in bilingual (Chinese/English) form.

## Required Inputs

See `references/input-contract.md`.

## Deliverable Checklist

See `references/deliverables.md`.
