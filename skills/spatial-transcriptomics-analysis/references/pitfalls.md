# Spatial Transcriptomics Pitfalls

Check this list before finalizing conclusions or handing off a spatial transcriptomics analysis. Each entry gives the trigger condition, why it happens, how it shows up in the output, and the fix.

## 1. Treating a Visium spot as a single cell

- **Trigger**: Assigning one cell type label per spot (e.g. from marker-gene scoring or Louvain cluster identity) and reporting it as "this spot is a macrophage."
- **Mechanism**: A standard 55 µm Visium spot typically captures transcripts from roughly 1-10 cells depending on tissue density (dense tissues such as spleen or tumor can exceed 10-20 nuclei per spot); the measured profile is a mixture, not a single transcriptome.
- **Symptom**: Cell-type maps show biologically implausible sharp single-spot "islands" of a rare cell type, or every spot gets a confident singular label with no mixed/ambiguous category.
- **Fix**: Run reference-based deconvolution (cell2location, RCTD, SPOTlight, or an equivalent method) and report per-spot cell-type *proportions*, not a single call. Only assign a dominant-type label when one proportion clearly exceeds the rest (e.g. >0.5) and say so explicitly as "dominant estimated type," not "cell type."

## 2. Domain clustering that ignores spatial coordinates entirely

- **Trigger**: Running standard scRNA-seq-style clustering (PCA + Louvain/Leiden on expression only) on spot data and calling the result "spatial domains."
- **Mechanism**: Expression-only clustering treats spots as exchangeable, unordered observations; it captures cell-composition similarity but not tissue architecture, so anatomically contiguous regions can be split and non-adjacent regions with similar composition can be merged.
- **Symptom**: Domain maps look "salt-and-pepper" — individual spots of a different color scattered inside an otherwise uniform region, with boundaries that do not track histological structure.
- **Fix**: Use a method that explicitly incorporates spatial neighbor structure — BayesSpace (Potts-model spatial prior), SpaGCN or STAGATE (spatial+histology graph, attention-weighted smoothing), or at minimum a spatial smoothing/neighbor-voting post-process on Leiden clusters. Report the neighbor graph definition (k or radius) used.

## 3. Over-smoothing that erases real fine-grained biology

- **Trigger**: Increasing the spatial smoothing strength, neighbor radius, or Potts/graph regularization weight to "clean up" a noisy domain map.
- **Mechanism**: Spatial-prior methods trade cluster purity for spatial coherence; pushed too far, thin real structures (a single-spot-wide vascular wall, a small immune infiltrate, a tumor-stroma interface only 1-2 spots wide) get absorbed into the majority label of their neighborhood.
- **Symptom**: Boundaries become unnaturally straight or overly large domains that erase known thin anatomical layers (e.g. cortical layers merge, glomeruli disappear in kidney sections) that are visible in the H&E image but not the domain map.
- **Fix**: Cross-check the domain map against the histology image and known tissue anatomy before reporting. Run the algorithm at more than one smoothing/neighbor setting and compare; prefer the setting that recovers known thin structures while still reducing salt-and-pepper noise, and state the parameter chosen and why.

## 4. Spatially variable gene calls confounded by a depth/count gradient

- **Trigger**: Reporting a gene as "spatially variable" from a raw or simply log-normalized count matrix, especially near tissue edges or damaged regions.
- **Mechanism**: Sequencing depth, tissue thickness, permeabilization efficiency, and RNA quality vary smoothly across a section (often lower at edges, higher in dense/nucleus-rich regions). Any gene correlated with total counts will show a spurious spatial pattern that reflects the depth gradient, not true biological spatial expression.
- **Symptom**: The top "spatially variable genes" list is dominated by high-abundance housekeeping/ribosomal genes, and the SVG spatial pattern visually mirrors the total-count or nFeature QC map.
- **Fix**: Run SVG detection on depth-corrected data (SPARK/SPARK-X model total counts as an offset; SpatialDE/Moran's I should be run on normalized, ideally variance-stabilized data) and always overlay the total-count spatial map alongside candidate SVG maps to visually rule out a depth-driven artifact before reporting a gene as spatially variable.

## 5. Tissue detection and spot filtering errors

- **Trigger**: Using Space Ranger's automatic tissue detection or a fixed count/gene threshold without visually reviewing the spot-under-tissue mask.
- **Mechanism**: Automatic tissue segmentation can miss faint or unstained tissue edges (false negatives, dropping real tissue spots) or include background/debris/bubbles as "under tissue" (false positives), and folded or torn tissue regions produce spots with anomalous, non-representative expression profiles.
- **Symptom**: Spatial maps show holes at tissue edges where real tissue was excluded, or isolated bright/dark artifact spots off the visible tissue outline; QC metrics show a bimodal count distribution that doesn't match the histology.
- **Fix**: Visually overlay the tissue detection mask and QC metrics (counts, genes, mito%) on the H&E image before analysis; manually correct the tissue/spot mask when it disagrees with the image, and explicitly exclude/flag spots over folds, bubbles, or tears identified from the image rather than by count thresholds alone.

## 6. Image-coordinate misalignment (scale factor / fiducial mismatch)

- **Trigger**: Re-processing exported coordinates and a downsampled or re-cropped image separately from the original Space Ranger `scalefactors_json`, or manually re-registering an image without carrying over the correct `tissue_hires_scalef`/`tissue_lowres_scalef`.
- **Mechanism**: Visium pixel coordinates in `tissue_positions.csv` are defined relative to the full-resolution image; every downstream image (hires/lowres) requires its own scale factor to map spot coordinates correctly. Using the wrong scale factor, a flipped y-axis, or a re-cropped image without adjusting the offset silently shifts every spot's apparent position.
- **Symptom**: Feature/domain maps look plausible in isolation but are offset, mirrored, or rotated relative to the true tissue structure when overlaid on the H&E image — e.g. a marker known to be epithelial-restricted appears in stroma.
- **Fix**: Always plot spots overlaid directly on the corresponding image using the matched scale factor before trusting any spatial map; treat any expression-histology mismatch (marker expression falling outside its expected tissue compartment) as a red flag for alignment error, not biology.

## 7. Neighborhood enrichment / ligand-receptor co-localization over-interpreted as direct signaling proof

- **Trigger**: Reporting a significant `squidpy.gr.nhood_enrichment` or ligand-receptor co-expression result (e.g. CellPhoneDB-style permutation test) as evidence that "cell type A signals to cell type B via ligand-receptor pair X."
- **Mechanism**: Spot-level co-localization only shows that two labels/domains are spatially closer than a random permutation baseline, and ligand/receptor transcript detection in nearby spots does not establish that the encoded proteins are co-present, secreted, and bound — diffusible ligands can act over multiple spot-widths, and Visium spots average signal from multiple cells so apparent "co-expression" can arise from within-spot mixing rather than adjacent-spot interaction.
- **Symptom**: Report language states signaling direction or causality ("A activates B through ligand X") based solely on spatial statistics, with no protein-level, receptor-downstream, or perturbation evidence.
- **Fix**: Frame results as "spatial co-localization consistent with potential signaling," state the neighbor radius/k used, note the ligand's known diffusion range if known, and recommend orthogonal validation (protein staining, receptor pathway activity score, or perturbation) before claiming a functional interaction.

## 8. Cross-slide/cross-capture-area batch effects ignored before comparison

- **Trigger**: Directly comparing gene expression, domain proportions, or SVG results between spots from different Visium capture areas, slides, or processing batches without any batch correction or explicit batch modeling.
- **Mechanism**: Each capture area is processed independently (permeabilization time, RNA quality, staining batch, sequencing run), introducing systematic technical differences on top of biological differences; because sections are the unit of replication, capture-area effects are fully confounded with any between-condition comparison unless corrected or blocked for.
- **Symptom**: Samples cluster by slide/batch rather than by biological condition in a joint embedding; "differentially expressed" genes between conditions are dominated by known batch markers (e.g. globally elevated mito or ribosomal genes in one batch).
- **Fix**: Integrate with a batch-aware method (Harmony, Seurat CCA/RPCA integration, or a spatial-aware integration method) before joint clustering, and for group comparisons use sample/section as the unit of replication (pseudobulk per sample) rather than pooling spots across sections into one test.

## 9. Deconvolution reference mismatch

- **Trigger**: Using a publicly available single-cell reference atlas that does not match the tissue's disease state, developmental stage, or perturbation condition (e.g. a healthy adult atlas as reference for tumor or fibrotic tissue).
- **Mechanism**: Deconvolution methods assume the reference cell-type signatures are representative of the cell states actually present in the spatial sample; disease-associated or context-specific cell states (tumor-associated macrophages, disease-specific fibroblast subtypes, malignant cell transcriptional states) are absent from a mismatched reference and get incorrectly reassigned to the nearest available "normal" type.
- **Symptom**: Deconvolution proportions show implausible patterns (e.g. no distinct "tumor cell" compartment estimated in an obvious tumor region, or a healthy cell type strongly enriched exactly where histology shows pathology).
- **Fix**: Use a reference derived from matched or closely related tissue, disease state, and, where possible, the same or a similar patient cohort; when only a mismatched reference is available, state this limitation explicitly and treat proportions as approximate/qualitative rather than precise quantitative estimates.

## 10. Post-hoc region comparisons without multiple-testing correction

- **Trigger**: Defining regions of interest by eye (or from an exploratory clustering result) and then running differential expression or statistical tests to "confirm" the region split is meaningful, especially across many candidate regions or many genes.
- **Mechanism**: When the grouping itself was chosen using the same expression data being tested (double-dipping), or when many region pairs/genes are tested, nominal p-values no longer reflect the true false-positive rate; post-hoc-selected boundaries are optimized to maximize apparent differences.
- **Symptom**: Extremely small p-values for genes/regions with boundaries that were manually adjusted after looking at the expression map; no FDR correction reported despite testing many gene-region combinations.
- **Fix**: Apply FDR (Benjamini-Hochberg or similar) correction across all genes/comparisons tested, treat exploratory-boundary-derived findings as hypothesis-generating rather than confirmatory, and where possible validate with an independent sample/section or an orthogonal readout (in situ, IHC) rather than re-testing on the same data used to draw the boundary.

## 11. Normalization choice that ignores spatial depth heterogeneity

- **Trigger**: Applying simple log(CP10K+1) normalization uniformly, or applying SCTransform with default parameters, without checking whether sequencing depth itself varies spatially (e.g. systematically lower at tissue edges or in necrotic regions).
- **Mechanism**: Standard log-normalization assumes a single global size factor is adequate; SCTransform's regularized negative-binomial model handles overdispersion and depth better but was originally tuned on scRNA-seq and can over- or under-correct when Visium depth gradients are large and spatially structured (rather than random per-cell noise).
- **Symptom**: Post-normalization feature maps still show a spatial gradient matching the raw total-count map, or SCTransform residuals show unexpected spatial banding.
- **Fix**: Inspect the total-count and detected-gene spatial maps first; if a strong spatial depth gradient exists, prefer a method with an explicit spatial or count offset (e.g. SCTransform with vst.flavor set appropriately, or SVG methods that model total counts as a covariate) and report the normalization method and any depth-covariate handling in the methods text.

## 12. Pseudoreplication in spot-level statistical tests

- **Trigger**: Running a standard single-cell-style Wilcoxon/t-test for "differential expression between region A and region B" using individual spots as the statistical unit, across one or a few tissue sections.
- **Mechanism**: Neighboring spots are spatially autocorrelated (near-identical due to shared cellular neighborhood and tissue continuity), so they are not independent observations; treating hundreds of non-independent spots as independent replicates inflates the effective sample size and produces artificially small p-values.
- **Symptom**: Very small p-values and huge numbers of "significant" genes between two spatial regions from a single section, with results that don't replicate when tested on an independent section from the same condition.
- **Fix**: Aggregate spots into pseudobulk per region per biological sample/section and test with a count-based model (edgeR/DESeq2) using sample as the replicate unit, or use a method explicitly designed for spatial autocorrelation-aware testing; do not treat spot count as sample size for between-condition inference.

## 13. Species/platform mismatch in ligand-receptor and marker databases

- **Trigger**: Running CellPhoneDB/CellChat/Squidpy ligand-receptor analysis or marker-based annotation on mouse (or other non-human) data using the default human-curated interaction database without ortholog mapping, or applying a Visium-derived marker panel to a different spatial platform (Xenium, Slide-seq, Stereo-seq) with different gene panel coverage or resolution.
- **Mechanism**: Ligand-receptor databases are typically curated primarily from human literature; direct gene-symbol matching against mouse data misses species-specific paralogs/orthologs and can silently drop valid interactions or include invalid ones. Cross-platform panel differences (targeted vs. whole-transcriptome, spot vs. subcellular resolution) change what is even detectable.
- **Symptom**: Unexpectedly short list of detected ligand-receptor pairs, or pairs reported using human gene symbols that don't exist in the mouse dataset's gene list.
- **Fix**: Explicitly convert databases via verified ortholog mapping before running L-R analysis on non-human data, and confirm the interaction database and platform gene panel are compatible before interpreting absence of a pair as biological absence rather than a database/panel gap.

## 14. RNA diffusion / lateral "bleeding" between adjacent spots

- **Trigger**: Interpreting a sharp expression boundary at single-spot resolution as a precise anatomical boundary, particularly for highly abundant transcripts near a strong expression source (e.g. next to a large gland or highly expressing tumor nest).
- **Mechanism**: During tissue permeabilization, mRNA can diffuse laterally a short distance before capture, and highly expressed transcripts from one region can appear at low levels in immediately adjacent spots that don't biologically express the gene; this is a known technical artifact distinct from a true expression gradient.
- **Symptom**: A cell-type-specific marker with a hard biological on/off boundary instead shows a smooth 1-2 spot "halo" of intermediate signal at the border in every section, consistent across unrelated genes with the same high-expression neighbor.
- **Fix**: Treat single-spot-resolution boundaries as approximate (± 1 spot diameter), avoid drawing strong conclusions from marker presence/absence in a single edge spot, and use domain/cluster-level rather than single-spot-level boundaries for anatomical claims when possible.

## 15. Applying single-cell trajectory/pseudotime directly to raw spots

- **Trigger**: Running a standard pseudotime/trajectory method (e.g. Monocle, Slingshot) on spot-level PCA embeddings and reporting a "differentiation trajectory" across the tissue without accounting for spot cell-mixture.
- **Mechanism**: Pseudotime methods assume each unit represents one cell's position along a continuous process; a spot's mixed-cell profile can produce an ordering that reflects changing cell-type composition across space (e.g. an anatomical gradient in cell-type mix) rather than a true single-cell differentiation trajectory.
- **Symptom**: A "trajectory" that aligns suspiciously well with a known anatomical axis (e.g. crypt-to-surface, cortex depth) and correlates strongly with deconvolved cell-type proportions changing together, rather than with a single cell type's internal state changing.
- **Fix**: Check deconvolution proportions along the inferred trajectory first; if the trajectory is driven by cell-type composition change rather than within-type state change, report it as a compositional/spatial gradient, and reserve "differentiation trajectory" language for cases supported by single-cell (or subspot-enhanced) resolution data or explicit compositional deconfounding.
