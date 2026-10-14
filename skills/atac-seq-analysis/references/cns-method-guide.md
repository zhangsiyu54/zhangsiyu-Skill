# CNS Method Guide

Use this guide before running bulk ATAC-seq, peak-table, motif, or scATAC-seq analysis.

## Minimal Inputs

- For bulk ATAC: require peak count matrix, sample metadata, genome build, peak coordinates, groups, contrast, blacklist/filtering status, and peak caller/annotation if known.
- For scATAC: require fragment files or Signac/ArchR object, cell metadata, clusters, genome build, and gene annotation.
- For existing outputs: require differential peak table, motif enrichment table, or gene activity matrix plus clear column meanings.

## QC and Filtering

- Bulk ATAC QC should include library depth, FRiP, duplicate rate, insert size periodicity, TSS enrichment when available, peak count, and sample correlation/PCA.
- scATAC QC should include fragments per cell, TSS enrichment, nucleosome signal, blacklist ratio, peak count per cell, doublet risk, and sample balance.
- Remove blacklist/problematic peaks and low-count peaks with explicit rules.
- Keep genome build consistent across peaks, annotation, motifs, and gene models.

## Differential Accessibility and Annotation

- Use DESeq2/edgeR/limma-style models on raw peak counts with covariates where possible.
- Annotate peaks to promoters, introns, intergenic regions, nearest genes, and distance to TSS; do not equate nearest gene with proven target gene.
- Motif enrichment requires a matched background peak set; report motif database/version and background strategy.
- Footprinting-style interpretation requires sufficient coverage and should be framed cautiously.

## scATAC Rules

- Use TF-IDF/LSI, depth-aware clustering, UMAP, peak calling or consensus peaks, gene activity scores, and motif deviation analysis when objects support them.
- Peak-to-gene links are correlation-based unless validated; do not claim enhancer-target proof.

## Required Figures

- QC metrics, PCA/UMAP, differential peak volcano/MA, peak heatmap, genomic annotation pie/bar, motif enrichment dot/bar plot, gene activity maps, peak-to-gene examples when requested.

## Methods Text Must Include

- Genome build, peak caller/source, blacklist and filtering, normalization/model, contrast, annotation database, motif database/background, scATAC dimensionality parameters, thresholds, and package versions.

## Algorithm and Design Decision Tree

Use this table to pick the tool/parameter route for a given task and cite the rationale in methods text. Treat rows as defaults to start from, not rigid rules — deviate and record the reason when the input contract or biology demands it.

| Scenario | Recommended tool/parameter | Rationale |
| --- | --- | --- |
| Bulk ATAC, two-group differential peaks | DESeq2 or edgeR on fixed-width (e.g., 500 bp) summit-centered consensus peak counts, with `apeglm`/`ashr` log2FC shrinkage; include peak width or GC as a covariate only if a width/GC-vs-log2FC trend is visible in diagnostics | Fixed-width peaks remove the length confound baked into raw variable-width MACS peaks; shrinkage stabilizes fold-change ranking for low-count peaks the same way it does for lowly-expressed genes in RNA-seq, without pretending peaks behave like genes |
| scATAC clustering from a fragment file | ArchR `IterativeLSI` (genome-wide tile matrix, 500 bp bins) or Signac `RunTFIDF` + `FindTopFeatures` + `RunSVD` on a peak/window matrix; drop or down-weight LSI/SVD components correlated with sequencing depth before UMAP/clustering | TF-IDF + LSI corrects for per-cell depth and per-feature ubiquity, which raw PCA on a sparse binary matrix cannot; iterative variants further remove depth-driven structure by re-selecting variable features from a first-pass clustering |
| Motif enrichment for a differential peak set | HOMER `findMotifsGenome.pl` with a GC/length-matched background peak set, or chromVAR for continuous per-sample/per-cell motif deviation when a background peak set and count matrix already exist | Matched background controls for base-composition-driven false enrichment (GC-rich motifs inflate against a naive genomic background); chromVAR is preferable when the question is "how does motif accessibility vary across cells/samples" rather than a one-shot enrichment test |
| Gene activity score for RNA integration | ArchR `addGeneScoreMatrix` (distance-weighted exponential decay model) for ArchR-native workflows; Signac `GeneActivity()` (simple gene body + promoter window counting) when working in a Seurat/Signac object; either used only as a label-transfer/integration proxy | Both are heuristic accessibility-to-gene proxies with different bias profiles (decay-weighted vs. flat window); pick the one native to the object type in hand rather than converting objects, and always caveat the score as a proxy, not measured expression |
| Peak-to-gene linkage | ArchR `addPeak2GeneLinks` or Signac `LinkPeaks`, restricting search window to a stated distance cutoff (commonly ~500 kb) and requiring a minimum correlation/FDR threshold across cells or metacells | Both are correlation-based across single-cell (or aggregated pseudo-bulk) profiles and require enough cells for a stable correlation estimate; state the distance window and threshold explicitly since they materially change the number of reported links |
| Chromatin footprinting for a specific TF | TOBIAS (`ATACorrect` + `ScoreBigwig` + `BINDetect`) as the general-purpose route with explicit Tn5 sequence-bias correction; HINT-ATAC as an alternative with its own position-dependency bias model; only on (pseudo)bulk aggregates with adequate depth | Bias correction is mandatory because raw Tn5 cut-site profiles carry strong sequence insertion bias that mimics footprint depletion; TOBIAS and HINT-ATAC differ in bias-model construction, so state which was used and do not mix outputs across tools in one figure |
| Cross-sample peak consensus/merging | Per-sample MACS2/3 narrowPeak calling, reproducibility filtering (IDR or overlap-based, e.g., peak present in >=2 replicates), then iterative greedy non-overlapping peak selection ranked by score (ArchR-style) rather than naive `bedtools merge` | Naive union/merge lets peak boundaries drift and inflates or fragments elements across samples; ranked greedy selection preserves fixed peak width and avoids double-counting the same regulatory element under multiple overlapping "peaks" |
| Co-accessibility / cis-regulatory module inference | Cicero (single-cell/metacell co-accessibility, Gaussian-graphical-model-based with a distance penalty) for scATAC; treat resulting modules as candidate cis-regulatory hypotheses requiring orthogonal validation (Hi-C/HiChIP, CRISPRi, eQTL) before causal claims | Cicero's distance-penalized correlation approach is specifically designed for the sparsity of single-cell accessibility data; co-accessibility is correlational evidence of coordinated regulation, not proof of physical looping or causal gene control |
| scATAC doublet removal before clustering | ArchR `addDoubletScores`/`filterDoublets`, or AMULET for fragment-file-based heterotypic doublet detection, applied before iterative LSI | Doublets create hybrid accessibility profiles that can masquerade as intermediate cell states; removing them before dimensionality reduction avoids contaminating the feature-selection step used for all downstream clustering |

See `references/pitfalls.md` before finalizing any of the above — several rows have a matched failure mode documented there (peak-width bias, LSI depth artifacts, motif background bias, footprint depth requirements, naive peak merging, doublet-driven clusters).
