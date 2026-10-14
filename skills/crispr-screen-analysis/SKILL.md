---
name: crispr-screen-analysis
description: Portable CRISPR screen analysis workflow for Claude Code (and other Agent-Skills-compatible coding agents such as Codex CLI). Use when the user provides MAGeCK-style gene summaries, sgRNA count tables, treatment/control groups, beta scores, logFC, FDR, essential gene lists, pathway annotations, or asks for CRISPR screen hit ranking, volcano plot, gene rank plot, positive/negative selection, pathway enrichment handoff, QC summary, or screen deliverables.
---

# CRISPR Screen Analysis

## Core Workflow

Before choosing scripts, read `references/paper-practice.md`, `references/cns-method-guide.md`, and `references/top-references.md`. Follow the closest CNS-style workflow, including input eligibility, QC order, statistical model, parameter reporting, classic visualizations, interpretation boundaries, and literature-backed methods wording.
1. Identify input type: sgRNA counts, MAGeCK gene summary, BAGEL/drugZ output, or existing hit table.
2. Confirm contrast, selection direction, gene identifiers, FDR threshold and screen design.
3. Run `scripts/check_environment.py` first. Use Python quickplot for existing gene summaries; use MAGeCK/BAGEL/drugZ when installed and raw counts are available.
4. Use `scripts/quickplot_crispr.py` for the volcano, gene-rank, and top-hits-bar figures; add `--sgrna-counts` when an sgRNA-level count table is available to also get the QC and sgRNA-consistency figures.
5. Read `references/visual-style.md` before plotting. Export only PNG/PDF unless explicitly requested otherwise.
6. After every run, write `output/methods.md` with `scripts/write_methods_doc.py`.
7. Before finalizing conclusions, check `references/pitfalls.md` (coverage/MOI, selection-direction method choice, control baseline, low-count filtering, batch/timepoint confounds, sgRNA-level consistency, multi-contrast FDR, orthogonal validation, CRISPRi/a assumptions, copy-number artifacts, normalization, off-target guides, enrichment background).
8. Do not call a hit validated without independent validation evidence.
9. After the core deliverables (report.md, methods.md) exist, ALWAYS also run `scripts/write_delivery_pptx.py --out output --skill-name crispr-screen-analysis` as a standard final step -- do not wait for the user to ask for a presentation deck; every completed analysis ships with one automatically, in bilingual (Chinese/English) form.

## Scripts

- `scripts/make_demo_data.py`: create a demo CRISPR gene summary plus an sgRNA-level count table (`--no-sgrna` to skip the latter).
- `scripts/quickplot_crispr.py`: from `--gene-summary` (required), generate the volcano plot, gene rank plot, and top-hits barplot (PNG/PDF) with staggered, collision-avoiding hit labels, plus `results/gene_summary.tsv`, `results/top_hits.tsv`, and `summary.tsv`/`run/qc.json`. Pass `--sgrna-counts <file>` (sgRNA-level raw counts) to additionally generate the QC distribution figure (read depth, sgRNAs-per-gene, replicate correlation, Gini index) and the sgRNA-consistency plot for top hits, plus `results/sgrna_summary.tsv`. Tune `--fdr-threshold`, `--lfc-threshold`, `--top-n-labels`, `--top-n-consistency` as needed.
- `scripts/check_environment.py`: inspect R/Python availability.
- `scripts/write_methods_doc.py`: write methods details.
- `scripts/write_delivery_pptx.py`: assemble `output/delivery.pptx` from the already-written report.md/methods.md/figures/summary.tsv (run last, after the scripts above).

## Required Inputs

See `references/input-contract.md`.

## Deliverable Checklist

See `references/deliverables.md`.
