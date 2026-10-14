# CNS Method Guide

Use this guide before running spatial transcriptomics analysis or reporting spatial figures.

## Minimal Inputs

- Require spatial coordinates, expression matrix/object, sample/section ID, spot/cell barcode, tissue image if available, group/domain labels, target genes/pathways, and platform.
- For 10x Visium, keep spot coordinates, scale factors, histology image, filtered matrix, and metadata together.
- For deconvolution, require a matched or appropriate single-cell reference and cell-type annotations.

## QC and Preprocessing

- Inspect spot/cell counts, detected genes, mitochondrial/ribosomal fraction, tissue coverage, section quality, and sample balance.
- Do not compare sections without tracking section/sample effects.
- Normalize within platform and sample. Keep raw counts for model-based methods when required.
- Align coordinates and image orientation before plotting; verify that spatial maps are not mirrored or shifted.

## Spatial Analysis

- Spatial feature maps should distinguish expression intensity from detection rate and sequencing depth.
- Spatially variable gene methods require enough spatial units and should report method/statistic/FDR.
- Domain detection/clustering should be interpreted with marker evidence and histology context, not only algorithmic clusters.
- Neighborhood enrichment and ligand-receptor spatial context require cell-type/domain labels and should be framed as spatial co-localization/inference.
- Deconvolution estimates are model-dependent; report reference source, method, and cell-type resolution.

## Required Figures

- Tissue/spatial coordinate QC, spatial domain plot, feature maps for target genes, marker DotPlot/heatmap, deconvolution composition map, spatially variable gene summary, region comparison plot.

## Interpretation Boundaries

- Do not claim cell-cell physical interaction from spot-level co-localization alone.
- Do not infer directionality of signaling without ligand/receptor and spatial evidence.
- State when spot resolution mixes multiple cells or when no histology image is available.

## Methods Text Must Include

- Platform, sample count, coordinate system, normalization, QC thresholds, clustering/domain method, SVG method, deconvolution reference and algorithm, spatial database versions, thresholds, and package versions.

## Algorithm and Design Decision Tree

Use this table to justify method choice before running an analysis; report the chosen method and parameters in the methods text.

| Scenario | Recommended tool/method | Key parameter(s) | Rationale |
|---|---|---|---|
| Unsupervised spatial domain detection, Visium, histology image available | SpaGCN or STAGATE (graph/attention over expression + spatial graph + histology pixel intensity) | neighbor graph radius/k, number of domains, (SpaGCN) histology weight `s` | Explicitly fuses spatial adjacency and image texture into the neighbor graph so domain boundaries track tissue architecture instead of expression noise alone; avoids salt-and-pepper clusters from expression-only Leiden/Louvain. |
| Unsupervised spatial domain detection, no usable image or need for statistical uncertainty on cluster assignment | BayesSpace | Potts spatial smoothing parameter, number of clusters `q`, optional subspot `enhance()` resolution factor | Bayesian spatial prior (Potts model) gives a principled, tunable trade-off between spatial coherence and expression fidelity, plus native subspot resolution enhancement without requiring an image. |
| Spatially variable gene (SVG) detection, large dataset (>10k spots) or need for fast genome-wide screen | SPARK-X (non-parametric) | kernel choice, FDR threshold | Avoids convergence failures of parametric count models at scale; orders-of-magnitude faster than SPARK/SpatialDE while retaining type-I error control, at the cost of slightly less power on small datasets. |
| Spatially variable gene detection, smaller dataset, need explicit spatial pattern classification (periodic/linear/general) | SpatialDE or SPARK (Gaussian process / generalized linear spatial model) | length-scale parameter, FDR threshold | Gaussian-process/GLSM formulation directly tests spatial vs. non-spatial variance and can classify pattern type; more statistically expressive than a simple Moran's I screen when dataset size is manageable. |
| Quick exploratory SVG screen or per-gene spatial autocorrelation score without a full generative model | Moran's I / Geary's C via `squidpy.gr.spatial_autocorr` | neighbor graph (k or Delaunay), permutation or analytic p-value, FDR | Fast, dependency-light first pass to rank candidate genes before committing to a full parametric SVG method; must still control FDR across all tested genes. |
| Cell-type deconvolution of Visium spots, need probabilistic absolute abundance and robustness to reference batch effects | cell2location | reference cell signature model, expected cells-per-spot prior, hardware (GPU recommended) | Full Bayesian model of both reference and spatial data jointly estimates absolute cell abundance (not just proportion) and explicitly models technical variation between reference and spatial platforms; heavier compute cost. |
| Cell-type deconvolution of Visium spots, need doublet-aware proportions with a well-matched scRNA-seq reference and moderate compute | RCTD (Robust Cell Type Decomposition) | mode (`doublet`/`full`/`multi`), platform effect normalization | Explicitly models spot as a controlled mixture of a small number of cell types (doublet mode) and corrects for platform-specific technical effects between reference and spatial assay; faster than cell2location, still count-based. |
| Cell-type deconvolution of Visium spots, small or shallow reference, need fast approximate proportions | SPOTlight (seeded NMF regression) | marker gene set per cell type, number of NMF components | Seeded NMF initialized from marker genes is robust to small or lower-depth references and is computationally light, at some cost in accuracy versus RCTD/cell2location on well-powered references. |
| Ligand-receptor spatial interaction / ROI signaling context | Squidpy `ligrec` (CellPhoneDB-style permutation test) combined with `nhood_enrichment` for co-localization | number of permutations, neighbor graph, FDR/p-value threshold | Co-localization (nhood_enrichment) establishes that two labels are spatially closer than chance; L-R permutation test adds specificity to which ligand-receptor pairs are jointly detected in adjacent/co-located spots — but both remain correlative, not proof of direct signaling (see pitfalls.md). |
| Multi-sample / multi-slide integration before joint clustering or comparison | Harmony or Seurat anchor-based integration (RPCA/CCA) on the PCA embedding, or STAGATE's multi-section 3D mode | batch/section covariate, number of anchors or Harmony theta | Standard batch-correction on the low-dimensional embedding removes systematic capture-area/slide effects while preserving biological structure; STAGATE's 3D mode additionally uses inter-section spatial continuity for methods that support it. |
| Region-of-interest differential expression across domains or conditions | Pseudobulk aggregation per sample/section per region, then edgeR/DESeq2 | aggregation unit (sample x region), FDR correction | Spots within a region are spatially autocorrelated and not independent; pseudobulk aggregation restores a valid statistical unit (the sample) and avoids pseudoreplication-driven false positives (see pitfalls.md #12). |
| Single-cell reference-guided label transfer to spots | Seurat/Signac anchor-based `TransferData`, or cell2location/Tangram for probabilistic spot-to-cell mapping | anchor features, k.weight, reference cell-type resolution | Anchor-based transfer projects reference cell-type scores onto spots using shared variable-feature space; probabilistic mapping methods (Tangram, cell2location) additionally model the many-cells-per-spot mixture explicitly rather than assigning a single nearest label. |
