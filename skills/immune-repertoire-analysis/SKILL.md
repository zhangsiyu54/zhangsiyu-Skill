---
name: immune-repertoire-analysis
description: Portable immune repertoire workflow for Claude Code (and other Agent-Skills-compatible coding agents such as Codex CLI). Use when the user provides TCR/BCR clonotype tables, VDJ contig annotations, scTCR/scBCR metadata, clonotype frequencies, sample groups, cell types, or asks for clonotype expansion, clonal diversity, V/J usage, CDR3 length, overlap, alluvial-style clone tracking, cell type-clonotype coupling, or immune repertoire deliverables.
---

# Immune Repertoire Analysis

## Core Workflow

Before choosing scripts, read `references/paper-practice.md`, `references/cns-method-guide.md`, and `references/top-references.md`. Follow the closest CNS-style workflow, including input eligibility, QC order, statistical model, parameter reporting, classic visualizations, interpretation boundaries, and literature-backed methods wording.
1. Identify input type: clonotype table, VDJ contig table, scTCR/scBCR metadata, or repertoire summary.
2. Confirm sample IDs, chain type, clonotype ID, count/frequency columns, group/celltype columns.
3. Run `scripts/check_environment.py` first. Use Python quickplots when specialized R packages are unavailable.
4. Use `scripts/quickplot_repertoire.py` for clonal expansion (clone-size class + top-N tracking), diversity (with depth caveat), V/J gene usage, CDR3 length, and clonotype overlap figures/tables.
5. Read `references/visual-style.md` before plotting. Export only PNG/PDF unless explicitly requested otherwise.
6. After every run, write `output/methods.md` with `scripts/write_methods_doc.py`.
7. Do not equate clonotype expansion with antigen specificity unless antigen evidence exists.
8. Before finalizing conclusions, check `references/pitfalls.md` for depth/normalization, clonotype-definition, multiplet, convergent-recombination, and BCR-lineage traps.
9. After the core deliverables (report.md, methods.md) exist, ALWAYS also run `scripts/write_delivery_pptx.py --out output --skill-name immune-repertoire-analysis` as a standard final step -- do not wait for the user to ask for a presentation deck; every completed analysis ships with one automatically, in bilingual (Chinese/English) form.

## Scripts

- `scripts/make_demo_data.py`: create demo clonotype table with sample/group/clonotype/count plus v_gene, j_gene, cdr3, celltype (CDR3/V/J identity is fixed per clonotype ID across samples so overlap/tracking analyses are meaningful).
- `scripts/quickplot_repertoire.py`: generate clonal expansion by clone-size class, top-N clonotype frequency (cross-sample clone tracking), diversity metrics (Shannon raw + rarefied, Gini-Simpson, inverse Simpson, Chao1, clonality) with a printed/written depth-sensitivity caveat, V and J gene usage, CDR3 length distribution, and pairwise clonotype overlap (Morisita-Horn + Jaccard) as PNG/PDF figures and TSV tables.
- `scripts/check_environment.py`: inspect R/Python availability.
- `scripts/write_methods_doc.py`: write methods details.
- `scripts/write_delivery_pptx.py`: assemble `output/delivery.pptx` from the already-written report.md/methods.md/figures/summary.tsv (run last, after the scripts above).

## Required Inputs

See `references/input-contract.md`.

## Deliverable Checklist

See `references/deliverables.md`.
