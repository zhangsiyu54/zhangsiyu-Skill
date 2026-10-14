# CNS Method Guide

Use this guide before running or reporting ORA, GSEA, pathway, or gene set analysis.

## Minimal Inputs

- Require organism, gene ID type, source analysis, test universe/background, and directionality.
- For ORA, require significant up/down gene lists and a valid background universe from the tested genes.
- For GSEA/fgsea, require a complete ranked list, preferably all genes with signed statistic, log2FC, or Wald/t statistic.

## ID Mapping and Databases

- Do not guess gene ID type or organism. Confirm SYMBOL, ENTREZ, ENSEMBL, UniProt, or custom IDs.
- Record one-to-many and unmapped ID rates. Keep an ID mapping table.
- Choose databases by biological question: GO BP/MF/CC for function, KEGG/Reactome for pathways, MSigDB Hallmark for compact pathway summaries, custom GMT for project-specific gene sets.
- Record database name, version/date, organism, and gene set size filters.

## Statistical Rules

- ORA uses a hypergeometric/Fisher-style test and must use a defensible background universe.
- GSEA uses ranked genes and avoids arbitrary DEG cutoffs; use permutation/fgsea parameters and report adjusted P values.
- Separate up/down ORA unless the question explicitly asks for combined gene lists.
- Apply FDR correction within each database/collection and avoid mixing databases without clear labels.

## Interpretation

- Prioritize terms with coherent leading-edge/core genes, reasonable gene ratios, and non-redundant biology.
- Collapse or cluster redundant GO terms when delivering to non-technical users.
- Do not claim pathway activation from ORA alone. Use direction-aware GSEA/GSVA/PROGENy-style evidence if activation is requested.

## Algorithm and Design Decision Tree

| Situation | Signal that identifies it | Recommended approach | Typical tool/function | Key caveat |
|---|---|---|---|---|
| Have significant DEG list only (no full ranking retained) | Only up/down gene symbols or a fixed cutoff table provided | ORA (hypergeometric/Fisher) | `clusterProfiler::enrichGO`/`enrichKEGG` with `universe` set to tested genes | Discards magnitude/rank information; sensitive to the arbitrary cutoff and background choice (see `pitfalls.md` #1). |
| Have a full ranked gene list (all genes, signed statistic or log2FC) | DE results table with all tested genes and a continuous statistic | GSEA/fgsea | `fgsea::fgsea`/`fgseaMultilevel` | Rank by a signed test statistic or shrunken log2FC, not raw/unshrunken log2FC (see `pitfalls.md` #2). |
| Many redundant significant GO terms (parent/child, near-duplicate names) | Long GO BP result table with overlapping gene membership across "significant" terms | Semantic similarity reduction | `clusterProfiler::simplify()` (GOSemSim-based) or external REVIGO pass | State the similarity cutoff used; do not silently drop terms without disclosing the reduction rule. |
| Cross-species gene set comparison (e.g., human MSigDB vs mouse data) | Gene set collection organism differs from the expression data organism | Curated ortholog mapping before matching | `orthogene`, `babelgene`, biomaRt/HGNC-MGI homology tables, or `msigdbr` organism-specific export | Never rely on uppercase symbol string-matching alone; report ortholog mapping success rate. |
| Pathway-level visualization across many comparisons/conditions | Multiple contrasts, clusters, or timepoints each with their own enrichment table | Enrichment map / term-similarity network, or a multi-group dot plot | `enrichplot::emapplot`, `enrichplot::cnetplot`, or `clusterProfiler::compareCluster` + dotplot | Restrict to a curated top-N term subset per comparison; unfiltered networks become unreadable. |
| Need leading-edge/core gene identification for a specific pathway | User asks which genes actually drive a GSEA hit | Extract fgsea leading-edge genes | `fgsea` result's `leadingEdge` column, or GSEA desktop's core enrichment output | Report leading-edge size alongside NES/padj; treat overlaps of 2-3 genes as statistically nominal, not a biological conclusion (see `pitfalls.md` #6). |
| Sparse replicates / noisy per-gene statistics, but full expression matrix available | Few biological replicates, wide confidence intervals on individual gene log2FC | Prefer GSEA on a signed/shrunken statistic over ORA on an arbitrary threshold | `fgsea` on `lfcShrink`-derived or Wald-statistic ranking | ORA cutoffs are unstable under low power; GSEA aggregates weak per-gene signal across a whole set. |
| Need per-sample pathway activity (not just a group-level enrichment call) | User wants a pathway score per sample for downstream correlation/clustering/survival analysis | Single-sample gene set scoring | `GSVA::gsva`, ssGSEA, or `PROGENy` for signaling-pathway activity | Distinct from ORA/GSEA — a per-sample score, not a group-comparison p value; do not conflate with GSEA NES. |
| Comparing enrichment across many clusters/conditions in one figure | scRNA-seq cluster markers, or several contrasts analyzed together | Cluster comparison workflow | `clusterProfiler::compareCluster` | Apply consistent background/database and FDR scope across all groups being compared. |
| Project-specific gene signature not present in public databases | User supplies a custom gene list/GMT (e.g., an in-house signature) | Custom GMT as GSEA/ORA input | `fgsea`/`clusterProfiler::enricher` with a user-supplied `TERM2GENE` table | State the signature's provenance and size; avoid circularity if the signature was itself derived from a public database being re-tested (see `pitfalls.md` #10). |

## Required Figures

- Dot/bubble plot with gene ratio, adjusted P value, and count.
- Bar plot for top non-redundant terms.
- GSEA running-score plot for selected pathways.
- Enrichment map or term-similarity network only when enough significant terms exist.

## Methods Text Must Include

- Input gene list/ranked metric, background universe, ID mapping rule, database collection/version, gene set size filters, enrichment test, FDR method, thresholds, and top-term selection rule.
