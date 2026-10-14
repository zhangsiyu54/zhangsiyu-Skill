# CNS Method Guide

Use this guide before running or reporting single-cell RNA-seq analysis.

## Minimal Inputs

- Accept Seurat object, h5ad, 10x matrix, expression matrix plus metadata, or exported embeddings/markers.
- Require cell barcodes, sample IDs, condition/group columns, cluster labels if precomputed, species, and marker genes expected by the user.
- For condition comparisons, require biological sample IDs. Cell-level tests without sample-aware design must be labeled exploratory.

## QC and Preprocessing

- Inspect nCount/UMI, nFeature/gene count, mitochondrial fraction, ribosomal fraction when available, cell complexity, and sample-level cell numbers.
- Use dataset-specific QC thresholds; do not apply universal cutoffs without checking distributions.
- Review doublet and ambient RNA risk. Use DoubletFinder/scDblFinder/Scrublet/SoupX/CellBender-style notes when those outputs exist; otherwise report limitation.
- Normalize with Seurat LogNormalize/SCTransform or Scanpy normalization. Record method and parameters.
- Regress covariates only when biologically justified; do not regress out the variable of interest.

## Integration, Clustering, and Annotation

- Use integration/batch correction only when batches/samples create technical structure; keep an unintegrated assay/layer for expression testing when appropriate.
- Select highly variable genes, run PCA, choose PCs from elbow/variance/JackStraw-style evidence, construct neighbor graph, and cluster with documented resolution.
- UMAP/t-SNE are visualization tools, not proof of lineage distance.
- Annotate cell types with multiple canonical markers, absence markers, and literature consistency. Do not invent labels from one marker.
- Keep original cluster IDs, final annotation labels, and marker evidence table.

## Differential and Advanced Modules

- For condition DE within cell types, prefer pseudobulk by sample/condition using edgeR/DESeq2/limma-voom when replicates exist.
- Cell-level Wilcoxon/MAST can be used for marker discovery or exploratory cell-type-local comparisons, but state sample-replication limitations.
- Trajectory requires relevant continuous biological process and enough intermediate states; do not run pseudotime for disconnected unrelated cell types.
- Cell-cell communication requires reliable cell-type labels and ligand/receptor database version; report it as inference from expression, not direct physical proof.

## Required Figures

- QC violin/histogram, UMAP by cluster/sample/condition/cell type, marker DotPlot/FeaturePlot, cell composition bar plot, top marker heatmap.
- For advanced work: pseudobulk volcano/heatmap, trajectory plot, module score plot, communication bubble/network.

## Methods Text Must Include

- Object/input type, filtering thresholds, normalization method, integration method, number of HVGs/PCs, clustering algorithm/resolution, marker test, annotation evidence, DE strategy, database versions, package versions, and known limitations.

## Algorithm and Design Decision Tree

Use this table to pick the concrete method/tool/parameter for common scenarios. Treat it as a starting point, not a substitute for checking the data (distributions, confounding, replicate counts) before committing to a route. Cross-check any decision here against `references/pitfalls.md`.

| Scenario | Recommended method/tool/parameter | Rationale |
|---|---|---|
| Multi-sample, multi-condition DE within one cell type, with biological replicates | Aggregate counts to pseudobulk per sample x cell type; test with DESeq2/edgeR-LRT or limma-voom, design including condition and known covariates (e.g., muscat-style workflow) | Preserves the sample as the unit of replication; avoids pseudoreplication inflation from cell-level tests (Squair et al. 2021; Crowell et al. 2020) |
| Condition DE with only 1-2 biological replicates per group, or exploratory marker screening within a cluster | Cell-level Wilcoxon/MAST, explicitly labeled exploratory; do not report as confirmatory | Pseudobulk needs enough samples per group for a real error term; without replicates, only descriptive/exploratory comparison is honest |
| Strong batch/technical structure across samples (different runs/chips/days), same expected cell types | Harmony (fast, scales to large data) or Seurat RPCA/CCA integration; verify with before/after UMAP split by sample and by cell type | Harmony offers strong mixing at low computational cost for shared-population data (Korsunsky et al. 2019; benchmarked well by Tran et al. 2020) |
| Cross-platform integration (e.g., 10x v2 vs v3, or different technologies: droplet vs plate-based) | scVI/scANVI (deep generative, models technical dispersion explicitly) or Seurat RPCA with `reference`-based anchors; benchmark with scIB metrics if feasible | Larger technical gaps benefit from a generative batch model rather than purely geometric correction (Lopez et al. 2018; Luecken et al. 2022 scIB benchmark) |
| Batch/condition are confounded by experimental design (e.g., all treated cells run in one batch) | Do not integrate away the confound; report the confound as a limitation, analyze descriptively, and recommend a better-designed follow-up experiment | Integration cannot distinguish technical batch signal from the confounded biological signal; forcing mixing risks erasing the real effect (see pitfalls.md #10) |
| Shallow sequencing / low counts per cell (e.g., high-throughput low-depth droplet runs) | Use looser but justified QC thresholds based on the dataset's own distribution; prefer SCTransform (better behavior at low counts) or careful log-normalization; avoid aggressive HVG/PC counts that fit noise | Low-depth data have more zero-inflation and dropout; universal high-depth thresholds will exclude most real cells |
| Rare cell population detection (expected <1-2% of cells) | Increase clustering resolution locally (subcluster a parent population) or use marker-informed subsetting before reclustering; confirm with a targeted marker panel and cross-sample presence, not a single UMAP island | Global default resolution frequently merges rare populations into a larger neighbor; subclustering with dedicated resolution recovers them without over-splitting common types |
| Trajectory/pseudotime across a continuous differentiation process | Slingshot (robust, tree-based, works well with existing cluster/UMAP structure) or monocle3 (large atlases, complex graphs); use tradeSeq for downstream gene-level dynamics along pseudotime | Both are well-benchmarked for real branching trajectories (Saelens et al. 2019 benchmark; Street et al. 2018; Cao et al. 2019) |
| Reference-based automated annotation available for the tissue/species | Azimuth or SingleR for well-characterized human/mouse references; CellTypist for immune-focused atlases; scANVI when a matched annotated scRNA-seq reference exists for label transfer | Reference mapping should corroborate, not replace, marker-based annotation; always report confidence/mapping scores and mark low-confidence transfers |
| No suitable reference exists (novel tissue, disease state, or non-model species) | Marker-panel-based manual annotation using canonical positive + absence markers per cluster; label ambiguous clusters `needs_annotation` | Forcing a reference-based label onto biology the reference never saw produces confidently wrong annotations |
| Cell-cell communication with a paired ligand-receptor database | CellChat for pathway-level, network-style overview with statistical testing of over-representation; CellPhoneDB or LIANA for multi-subunit complex-aware, resource-comparable analysis; NicheNet when the goal is linking ligands to downstream target gene expression changes | Different tools emphasize different outputs (interaction strength vs ligand-to-target regulatory inference); LIANA additionally provides a consensus across resources/methods (Jin et al. 2021; Efremova et al. 2020; Dimitrov et al. 2022; Browaeys et al. 2020) |
| Doublet detection needed, single-sample 10x droplet data | scDblFinder or DoubletFinder, set expected doublet rate from the loading concentration, run before merging/filtering | Both perform well in independent benchmarks (Xi & Li 2021); running pre-merge avoids cross-sample doublet miscalls (see pitfalls.md #3) |
| Ambient RNA suspected (dominant high-RNA cell type present, e.g., hepatocytes, erythrocytes, plasma cells) | SoupX or DecontX/CellBender correction before finalizing marker calls | Ambient contamination scales with cell size/total UMI and can mimic broad low-level "marker" expression in unrelated clusters (Young & Behjati 2020) |
| Multiple independent datasets need a single joint atlas-scale integration (tens of thousands to millions of cells, many labs/conditions) | scVI/scANVI or Harmony at scale; validate quantitatively with scIB metrics (batch mixing + biological conservation) rather than UMAP alone | Atlas-scale integration needs an integration method that scales and a benchmark-grounded way to check both over- and under-correction (Luecken et al. 2022) |
