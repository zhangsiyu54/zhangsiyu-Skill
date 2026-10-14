# ATAC-seq Pitfalls

Consult this list before finalizing QC calls, differential accessibility results, motif/footprint interpretation, or scATAC clustering conclusions. Each entry: Trigger (when it happens), Mechanism (why), Symptom (what you would see), Fix (what to do).

## 1. RNA-seq-style DESeq2 defaults applied to ATAC peak counts

- **Trigger**: running `DESeqDataSetFromMatrix` + default `DESeq()` on a peak-by-sample count matrix exactly as if peaks were genes.
- **Mechanism**: DESeq2's median-of-ratios normalization and dispersion model were designed for genes with roughly stable length and no fragment-count-vs-region-size confound; ATAC peaks range from ~150 bp to several kb, and read count scales with both peak width and local GC/mappability, not just "activity." Size factors computed peak-wise can be distorted by width heterogeneity.
- **Symptom**: differential peaks systematically skew toward wider peaks (or GC-rich regions if using bulk WGA-amplified libraries); MA plot shows a length- or GC-dependent trend in log2FC even after normalization; unusually large numbers of "significant" peaks in one width bin.
- **Fix**: fix peak width by using fixed-width (e.g., 250-501 bp) summit-centered peaks before counting (ENCODE/ArchR convention) so length is not a confound; if variable-width peaks must be kept, include peak width as a covariate or use count-per-bp normalization; inspect log2FC vs. peak width and vs. GC content diagnostic plots before trusting the DA list.

## 2. TSS enrichment / FRiP thresholds applied without checking fragment-size distribution

- **Trigger**: accepting a sample because TSS enrichment > 5-7 and FRiP > 0.2-0.3 without ever plotting the fragment-size histogram.
- **Mechanism**: TSS enrichment and FRiP are summary scalars that can pass acceptable thresholds even when the underlying library lacks the expected nucleosomal ladder (sub-nucleosomal ~<100 bp, mono-nucleosome ~180-247 bp, di-nucleosome ~315-473 bp periodicity), which is the classical signature of a successful Tn5 transposition reaction.
- **Symptom**: a sample with borderline-passing FRiP/TSS but a flat or monotonically decaying fragment-size distribution (no ~147 bp periodicity bump) — often indicative of over-digestion, degraded chromatin, or high mitochondrial/free-DNA contamination that got filtered before the metrics were computed, masking the underlying problem.
- **Fix**: always plot and report the fragment-size distribution alongside FRiP/TSS; treat the three metrics as a joint QC signature, not independent gates; flag samples with normal FRiP/TSS but absent nucleosomal periodicity for re-examination of transposition efficiency and chromatin integrity before inclusion.

## 3. Equating "differentially accessible near a promoter" with "differentially transcribed"

- **Trigger**: reporting a promoter-proximal differential peak as if it proves the associated gene changed expression.
- **Mechanism**: promoter accessibility can be necessary but is not sufficient for transcription (poised/paused Pol II, bivalent chromatin, post-transcriptional/translational regulation all decouple the two); distal (enhancer/intronic/intergenic) differential peaks are even further removed from any specific gene's transcription without independent linkage evidence.
- **Symptom**: accessibility-based "differential gene expression" claims that are not corroborated by RNA-seq/scRNA-seq, or bold statements like "peak X activates gene Y" derived only from nearest-gene annotation.
- **Fix**: report promoter peaks as "increased/decreased promoter accessibility," never as "increased/decreased expression," unless paired transcriptomic data confirm the direction; for distal peaks, state nearest-gene annotation as a hypothesis only and require peak-to-gene correlation (Cicero/ArchR/Signac) or orthogonal chromatin-conformation evidence before implying a target gene.

## 4. Motif enrichment reported without a GC/mappability-matched background

- **Trigger**: running HOMER `findMotifsGenome.pl`, chromVAR, or a hypergeometric motif test against a generic genomic background (random genome regions or all peaks genome-wide) instead of a background matched to the foreground peak set's GC content and length distribution.
- **Mechanism**: TF motifs vary widely in GC content (e.g., GC-box/Sp1-family motifs vs. AT-rich homeodomain motifs); if the foreground peak set (e.g., promoter-proximal opening peaks) is systematically more GC-rich than a naive genomic background, GC-rich motifs will appear "enriched" purely from base composition, independent of true TF binding signal.
- **Symptom**: enrichment results dominated by GC-rich motif families (SP/KLF, ZIC, EGR) regardless of biological contrast; enrichment persists even in permuted/shuffled control peak sets.
- **Fix**: use a background of peaks matched for GC content, length, and (ideally) mappability — HOMER's `-bg` with matched peaks, chromVAR's built-in GC-bias-corrected background peak sampling, or explicit dinucleotide-shuffled sequences; always run the same pipeline on a negative-control contrast (e.g., random peak split) to confirm the enrichment signal is not an artifact.

## 5. Independently-called per-sample peaks unioned/merged without controlling multiple testing inflation

- **Trigger**: calling peaks separately in each sample (e.g., MACS2 per BAM) then taking `bedtools merge` or a simple union of all peaks as the "master" peak set for counting.
- **Mechanism**: naive union/merge chains overlapping peaks from different samples end-to-end, so peak boundaries drift and the "same" regulatory element can fragment into multiple overlapping master peaks (or merge distinct elements into one inflated peak); counting reads over this ad hoc set and then running one differential test per row inflates the effective number of independent tests and biases width/summit position sample-to-sample.
- **Symptom**: unusually long or irregular merged peak widths; adjacent "peaks" that are clearly one element split by a few bp; large peak counts that shrink drastically after simple filtering; unstable differential results when peak-calling parameters change slightly.
- **Fix**: build a consensus/reproducible peak set properly — call peaks per sample and per pooled/pseudo-replicate, use IDR (Irreproducible Discovery Rate) or ENCODE-style overlap-reproducibility filtering, then generate fixed-width, non-overlapping consensus peaks by iterative greedy selection ranked by score (highest-scoring summit wins, overlapping lower-scoring peaks removed), as done by ArchR's iterative overlap peak merging or DiffBind's consensus approach; do not simply `bedtools merge -i union.bed`.

## 6. scATAC clustering directly on the binary/raw peak-cell matrix without TF-IDF/LSI

- **Trigger**: running PCA + Louvain/Leiden directly on a raw or binarized peak-by-cell count matrix (as one would for scRNA-seq counts) without a TF-IDF transform or LSI step.
- **Mechanism**: scATAC data is extremely sparse (typically 1-10% of accessible peaks detected per cell) and dominated by per-cell sequencing depth; raw PCA on this matrix is driven overwhelmingly by total fragment count per cell rather than biological state, because depth variance dwarfs biological signal without the log/IDF down-weighting of ubiquitously accessible peaks.
- **Symptom**: the first principal component (or first LSI component in Signac before dropping it) correlates strongly with total fragments/FRiP per cell rather than any biological covariate; clusters that map cleanly onto library-size bins or single-cell coverage instead of cell type.
- **Fix**: always use TF-IDF normalization followed by SVD/LSI (Signac's `RunTFIDF` + `RunSVD`, or ArchR's `IterativeLSI`); check the correlation of each LSI/SVD component with sequencing depth and drop or down-weight components highly correlated with depth (commonly component 1) before clustering/UMAP.

## 7. Batch correction (Harmony/scVI/integration anchors) removing real biological differences

- **Trigger**: applying strong batch-integration methods (Harmony, Seurat CCA/RPCA integration, scVI/PeakVI) across conditions of biological interest (e.g., treatment vs. control, disease vs. healthy) rather than across technical batches (e.g., sequencing runs, chips).
- **Mechanism**: integration algorithms are designed to force overlap between corresponding cell populations across batches; if the "batch" variable is confounded with the biological variable under study, the algorithm can over-correct and erase or dampen the very accessibility shifts that represent the biological effect of interest (differential abundance of cell states, condition-specific accessibility programs).
- **Symptom**: expected condition-specific clusters/cell states disappear after integration; per-cluster differential accessibility between conditions collapses to near-zero after integration despite being present in per-sample analysis pre-integration; integration diagnostics (e.g., kBET, LISI) show suspiciously perfect mixing of conditions that should differ biologically.
- **Fix**: never mix true experimental design factors into a purely technical batch-correction step; run integration diagnostics before and after; perform an unintegrated per-sample/per-condition analysis as a sanity check; prefer methods that explicitly model condition as a covariate to preserve (e.g., MOFA+/GLM-based approaches, or Harmony with condition included as `do_pca=FALSE` reference) rather than blanket removal.

## 8. Over-interpreting TF footprints from low-depth or low-quality data

- **Trigger**: running TOBIAS, HINT-ATAC, or manual Tn5-insertion-profile footprinting on bulk or pseudobulk aggregates with modest depth (well under the ~50-200M unique non-mitochondrial read range typically used in footprinting benchmarks), or on single-cell pseudobulks with few aggregated cells.
- **Mechanism**: footprint depletion is a subtle base-pair-resolution signal (a dip of tens of bp within an open region) riding on top of Tn5 sequence-insertion bias; with low depth, Poisson/sampling noise in per-base cut counts can produce spurious dips or mask real ones, and Tn5 sequence bias itself (not just protein occupancy) generates footprint-like patterns that require an explicit bias-correction model (e.g., k-mer or PDM-based) to distinguish from true occupancy.
- **Symptom**: footprint depth/scores that are unstable across technical replicates or bootstrap subsamples; footprints reported for motifs at only a handful of instances genome-wide; footprint calls that do not reproduce when depth is downsampled.
- **Fix**: state the aggregate read depth used for footprinting and only attempt genome-wide footprinting when depth is adequate (aggregate/pseudobulk on the order of tens of millions of usable reads, per tool documentation); always use a bias-corrected method (TOBIAS `ATACorrect`, HINT-ATAC's bias models) rather than raw cut-site profiles; validate top hits, when possible, against orthogonal TF ChIP-seq/CUT&RUN rather than presenting footprint scores as direct proof of binding.

## 9. Gene activity / gene score treated as a direct expression measurement

- **Trigger**: reporting ArchR gene scores or Signac `GeneActivity` values as if they were normalized RNA expression, or using them directly for expression-style differential testing/marker gene claims without the caveat that they are chromatin proxies.
- **Mechanism**: gene activity models are heuristics that aggregate accessibility signal in a gene body/promoter window (often with a distance-decay weighting and sometimes gene-body-length normalization) as a proxy for potential regulatory input; they do not measure transcript abundance and can be biased by gene length, neighboring gene density, and the specific decay/window parameters chosen by the tool.
- **Symptom**: gene activity ranks used directly to name "top expressed" or "marker" genes for a cluster without cross-referencing RNA data; large discrepancies between gene-activity-based and RNA-based marker lists that go unexplained; longer genes or those with more nearby regulatory peaks systematically ranking higher.
- **Fix**: label gene activity outputs explicitly as "predicted gene activity from accessibility" wherever they are reported; when a matched RNA/multiome modality exists, use it (or cross-modal integration/prediction) for actual expression claims and use gene activity only for label transfer, cluster-to-cell-type mapping, or as a complementary track; state the specific gene-activity method and parameters (ArchR model, Signac window) in the methods text.

## 10. ArchR/Signac version-specific default parameter drift silently changing results

- **Trigger**: re-running an "identical" ArchR or Signac pipeline on a different package version without pinning versions or reviewing changelogs (e.g., default LSI dimensions, default resolution for clustering, default TSS window, default doublet-scoring parameters, default peak-calling backend changing from MACS2 to MACS3-based calling).
- **Mechanism**: both packages have changed function defaults across releases (e.g., ArchR's default `addIterativeLSI` iterations/variable-feature counts, `addClusters` resolution, or Signac's default `FindTopFeatures`/`RunSVD` behavior and peak-calling wrapper defaults); silently relying on "whatever the installed version does by default" makes results non-reproducible across environments and time.
- **Symptom**: re-analysis on a new machine/container produces different cluster counts, different top variable features, or different peak counts from an old report using "the same script"; collaborators using a different package version cannot reproduce figures exactly.
- **Fix**: record exact package/tool versions (R, Bioconductor, ArchR, Signac, Seurat, MACS, TOBIAS) in the methods and `output/run/environment.json`; explicitly set (not rely on defaults for) key parameters — LSI dimensions used, clustering resolution, variable feature count, TSS window, doublet method/threshold, peak caller and version; re-verify key results after any package upgrade before treating them as final.

## 11. ENCODE blacklist regions not filtered before or after peak calling

- **Trigger**: skipping blacklist filtering because "MACS2 already looks reasonable" or because the blacklist file for the exact genome build/version was not readily available.
- **Mechanism**: blacklist regions (e.g., ENCODE/Boyle-Lab blacklist, satellite repeats, centromeric/telomeric regions, high-copy-number segmental duplications) generate anomalously high, structurally-driven read pileups independent of true chromatin state or cell type, because of mapping artifacts and repetitive sequence collapse in the reference assembly.
- **Symptom**: a small number of extremely tall, narrow peaks recurring in essentially every sample/cell type regardless of biology (often at the same few hundred loci); inflated global correlation between unrelated samples driven by these shared artifact peaks; blacklist-region reads inflating FRiP.
- **Fix**: filter reads/peaks against the genome-build-matched ENCODE blacklist (or `excluderanges`-style curated exclusion sets) both before QC-metric computation and before differential/motif analysis; verify blacklist genome build matches the alignment build exactly (hg19 vs. hg38 vs. mm10 blacklist are not interchangeable); re-inspect top peaks by score for any that coincide with un-filtered artifact regions.

## 12. Fixed-width or fragment-based peak sets built with the wrong TSS/gene annotation build

- **Trigger**: computing TSS enrichment, peak annotation, or gene activity using a TSS/gene model (e.g., RefSeq vs. GENCODE, or hg19 TSS coordinates) that does not match the genome build/version the reads were aligned to.
- **Mechanism**: TSS coordinates and gene model boundaries differ meaningfully between annotation sources and genome patch versions; a build mismatch silently shifts every distance-to-TSS, promoter/intron/intergenic classification, and TSS-enrichment aggregate profile without throwing an error.
- **Symptom**: unexpectedly low or noisy TSS enrichment scores despite otherwise good libraries; peak annotation categories that look implausible (e.g., very few promoter peaks); systematic off-center TSS-enrichment aggregate plot (peak not centered on position 0).
- **Fix**: always confirm and record which genome build and annotation source/version (GENCODE vs. RefSeq, release number) were used for TSS enrichment, peak annotation, and gene activity, and confirm they match the alignment build; when in doubt, regenerate the TSS reference from the same TxDb/annotation package used elsewhere in the same pipeline.

## 13. scATAC doublets and mitochondrial/organelle contamination inflating apparent clusters

- **Trigger**: proceeding directly to clustering after basic per-cell fragment/TSS filtering without dedicated doublet detection (e.g., AMULET, ArchR `addDoubletScores`, scDblFinder-style heterotypic doublet estimation) or without checking mitochondrial/chloroplast fragment fraction.
- **Mechanism**: heterotypic doublets (two different cell types profiled as one barcode) produce hybrid accessibility profiles that can form spurious "intermediate" clusters or spuriously connect otherwise distinct populations; high mitochondrial-fragment fraction usually reflects damaged/apoptotic nuclei with degraded, non-representative chromatin.
- **Symptom**: a small cluster expressing marker peaks/genes from two otherwise well-separated cell types simultaneously; a "bridge" of cells connecting two clusters in UMAP space that disappears after doublet removal; unexpectedly high mitochondrial fragment percentage in certain clusters.
- **Fix**: run a dedicated doublet-detection step appropriate to the platform and record the score/threshold used; filter or flag high-mitochondrial-fraction cells with an explicit threshold; re-cluster after removing predicted doublets/damaged cells and compare cluster stability before and after.

## 14. Confusing chromVAR/motif deviation z-scores with absolute per-cell TF activity levels

- **Trigger**: reporting chromVAR deviation z-scores (or bias-corrected deviations) across clusters/conditions as if they were calibrated, comparable-in-magnitude measures of "how much TF X is active," rather than relative, background-set-dependent statistics.
- **Mechanism**: chromVAR z-scores are computed relative to a background peak set matched for GC content and accessibility, and are inherently relative to the specific peak set, motif database, and background sampling used in that run; they are not on an absolute or cross-study-comparable scale, and z-score magnitude conflates motif prevalence in the peak set with true TF occupancy change.
- **Symptom**: cross-study or cross-motif-database comparisons of raw chromVAR z-score magnitudes presented as if directly comparable; claims of "TF X is N-fold more active" derived from z-score ratios.
- **Fix**: treat chromVAR deviations as relative, within-study rankings/effect directions, not absolute activity units; always report the motif database/version and background strategy used; when comparing across studies, re-run chromVAR consistently on a common peak/motif set rather than comparing previously published z-scores directly.
