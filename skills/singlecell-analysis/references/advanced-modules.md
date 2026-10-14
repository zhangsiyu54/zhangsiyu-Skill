# Single-cell Advanced Modules

Load this file when the request goes beyond basic QC/UMAP/marker plots.

## QC and Filtering

- Report cells, genes, median UMI/counts, mitochondrial percentage and detected gene distributions.
- Flag likely issues: very high UMI cells, extreme mitochondrial percentage, low detected genes, batch-specific QC shifts.
- Optional tools by context: DoubletFinder/scDblFinder/Scrublet for doublet review; SoupX/DecontX for ambient RNA review.
- Do not remove cells silently. Record thresholds and before/after counts.

## Normalization and Integration

- Basic single-dataset route: log-normalization, highly variable genes, scaling, PCA, neighbors, UMAP, clustering.
- Multi-sample route: SCT or log-normalization plus integration method such as Seurat integration, Harmony, scVI/scANVI, or user-specified method.
- Always compare integration before/after: UMAP by sample/batch, UMAP by cell type/cluster, batch mixing summary.
- Do not over-integrate away biological condition differences; mention this risk in reports.

## Annotation

- Start from marker evidence and known biology.
- If reference mapping is requested, use Azimuth/SingleR/CellTypist/scANVI only when species and tissue reference are appropriate.
- Deliver `results/cluster_markers.tsv`, `results/annotation_evidence.tsv`, and `figures/umap_celltype.png/.pdf`.
- Label weak clusters as `needs_annotation` instead of inventing confident names.

## Pseudobulk Differential Expression

- Prefer pseudobulk for condition comparisons within cell types when biological replicates exist.
- Aggregate counts by `sample + celltype`, then use DESeq2/edgeR/limma-voom with design including condition and known covariates.
- Stop and ask when there are no biological replicate columns.
- Deliver per-celltype DEG tables, volcano plots, summary table, and method notes.

## Trajectory and Pseudotime

- Require a biologically plausible lineage subset and starting population/root.
- Suitable tools: Slingshot, monocle3, tradeSeq. Avoid trajectory claims on unrelated cell types.
- Deliver trajectory embedding, pseudotime distribution by group, branch markers and limitations.

## Cell-cell Communication

- Require cell type labels and condition/group labels.
- Suitable tools: CellChat, NicheNet, LIANA. State database and species.
- Deliver ligand-receptor table, pathway-level communication, network figures and group comparison if conditions exist.

## Module Scores and Functional Interpretation

- Use curated gene sets or user-provided markers.
- Deliver per-cell module score distributions, UMAP score maps, group/celltype summaries.
- Do not infer pathway activation from one marker gene alone.
