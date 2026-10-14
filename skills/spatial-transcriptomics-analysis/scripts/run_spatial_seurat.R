#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage: Rscript run_spatial_seurat.R <spatial_input> <target_genes> <output_dir> [cluster_resolution]\n",
    "  spatial_input:      path to a 10x Visium output directory (loaded with Load10X_Spatial) or a saved\n",
    "                       Seurat spatial object (.rds)\n",
    "  target_genes:       comma-separated gene symbol(s) for spatial feature maps, e.g. 'Gfap,Pecam1'\n",
    "  output_dir:         output directory\n",
    "  cluster_resolution: optional Louvain resolution for FindClusters (default 0.8)"
  )
}
input <- args[[1]]
target_genes <- trimws(strsplit(args[[2]], ",")[[1]])
out_dir <- args[[3]]
resolution <- if (length(args) >= 4) as.numeric(args[[4]]) else 0.8

required <- c("Seurat", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "), "\nInstall Bioconductor/CRAN dependencies before running formal analysis.")
}

dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)

save_gg <- function(plot, name, width = 7.2, height = 6.4) {
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(name, ".png")), plot, width = width, height = height, dpi = 320, bg = "white")
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(name, ".pdf")), plot, width = width, height = height, bg = "white")
}

## ---- 1. Load input: Visium directory (Space Ranger output) or a saved Seurat spatial object ----
if (dir.exists(input)) {
  message("Loading 10x Visium directory via Seurat::Load10X_Spatial(): ", input)
  obj <- Seurat::Load10X_Spatial(
    data.dir = input,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "slice1",
    filter.matrix = TRUE
  )
} else if (file.exists(input)) {
  message("Loading saved Seurat spatial object: ", input)
  obj <- readRDS(input)
} else {
  stop("Input not found: ", input)
}
assay_name <- Seurat::DefaultAssay(obj)
count_col <- paste0("nCount_", assay_name)
feature_col <- paste0("nFeature_", assay_name)

## ---- 2. QC (spot counts/features, mitochondrial fraction, drop empty spots) ----
n_total_spots <- ncol(obj)
mt_genes <- grep("^MT-|^mt-", rownames(obj), value = TRUE)
if (length(mt_genes) > 0) {
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-|^mt-", assay = assay_name)
}
if (count_col %in% colnames(obj@meta.data)) {
  keep_cells <- colnames(obj)[obj@meta.data[[count_col]] > 0]
  if (length(keep_cells) == 0) {
    stop("All spots have zero ", count_col, "; nothing passes QC. Check the input matrix/filtering.")
  }
  obj <- subset(obj, cells = keep_cells)
}
n_qc_pass <- ncol(obj)
qc_metric <- if (count_col %in% colnames(obj@meta.data)) count_col else feature_col
p_qc <- Seurat::SpatialFeaturePlot(obj, features = qc_metric) + ggplot2::labs(title = "Spot count QC")
save_gg(p_qc, "spot_qc_map")

## ---- 3. Normalization (SCTransform, variance-stabilizing negative-binomial regression) ----
obj <- Seurat::SCTransform(obj, assay = assay_name, verbose = FALSE)
Seurat::DefaultAssay(obj) <- "SCT"

## ---- 4. Spatial feature maps for requested target genes ----
present_genes <- target_genes[target_genes %in% rownames(obj)]
absent_genes <- setdiff(target_genes, present_genes)
if (length(absent_genes)) {
  message("Warning: requested gene(s) not found in object and skipped: ", paste(absent_genes, collapse = ", "))
}
if (length(present_genes) == 0) {
  stop("None of the requested target genes were found in the object: ", paste(target_genes, collapse = ", "))
}
p_feat <- Seurat::SpatialFeaturePlot(obj, features = present_genes)
save_gg(
  p_feat, "spatial_feature_map",
  width = 7.2 * min(length(present_genes), 3),
  height = 6.4 * ceiling(length(present_genes) / 3)
)

## ---- 5. PCA, neighbor graph, clustering, UMAP ----
obj <- Seurat::RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- Seurat::FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- Seurat::FindClusters(obj, resolution = resolution, verbose = FALSE)
obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = 1:30, verbose = FALSE)

p_domain <- Seurat::SpatialDimPlot(obj, label = TRUE, label.size = 3) + ggplot2::labs(title = "Spatial domains")
save_gg(p_domain, "spatial_domain_map")

p_umap <- Seurat::DimPlot(obj, reduction = "umap", label = TRUE) + ggplot2::labs(title = "UMAP of spatial clusters")
save_gg(p_umap, "umap_clusters")

## ---- 6. Spatially variable genes (Moran's I over the neighbor spatial graph) ----
svg_table <- NULL
svg_error <- NULL
tryCatch({
  n_test <- min(2000, length(Seurat::VariableFeatures(obj)))
  obj <- Seurat::FindSpatiallyVariableFeatures(
    obj,
    assay = "SCT",
    features = Seurat::VariableFeatures(obj)[seq_len(n_test)],
    selection.method = "moransi"
  )
  ranked <- Seurat::SpatiallyVariableFeatures(obj, assay = "SCT", selection.method = "moransi")
  svg_table <- data.frame(rank = seq_along(ranked), gene = ranked)
  write.table(svg_table, file.path(out_dir, "results", "spatially_variable_genes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}, error = function(e) {
  svg_error <<- conditionMessage(e)
  message("Spatially variable gene detection failed (continuing without it): ", svg_error)
})

## ---- 7. Result tables ----
qc_df <- obj@meta.data
qc_df <- cbind(spot = rownames(qc_df), qc_df)
write.table(qc_df, file.path(out_dir, "results", "spot_qc.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

domain_summary <- as.data.frame(table(obj@meta.data$seurat_clusters))
colnames(domain_summary) <- c("cluster", "n_spots")
agg_cols <- intersect(c(count_col, feature_col), colnames(obj@meta.data))
for (col in agg_cols) {
  means <- tapply(obj@meta.data[[col]], obj@meta.data$seurat_clusters, mean)
  domain_summary[[paste0("mean_", col)]] <- round(as.numeric(means[as.character(domain_summary$cluster)]), 2)
}
write.table(domain_summary, file.path(out_dir, "results", "domain_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

n_clusters <- nrow(domain_summary)
small_cluster_threshold <- max(5, round(0.03 * n_qc_pass))
small_clusters <- domain_summary[domain_summary$n_spots < small_cluster_threshold, , drop = FALSE]

## ---- 8. Summary + report ----
write.table(
  data.frame(
    item = c("spots", "clusters", "target_genes_found", "svg_computed"),
    value = c(
      ncol(obj),
      length(unique(obj@meta.data$seurat_clusters)),
      paste(present_genes, collapse = ";"),
      !is.null(svg_table)
    )
  ),
  file.path(out_dir, "summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE
)

top_n_svg <- if (!is.null(svg_table)) min(5, nrow(svg_table)) else 0
top_svg_genes <- if (top_n_svg > 0) svg_table$gene[seq_len(top_n_svg)] else character(0)

## ---- automatic caveats tied to references/pitfalls.md, only printed when they apply ----
caveats <- c()
if (nrow(small_clusters) > 0) {
  detail <- paste0(small_clusters$cluster, " (n=", small_clusters$n_spots, ")", collapse = ", ")
  caveats <- c(caveats, paste0(
    "- Cluster(s) with very few spots (< ", small_cluster_threshold, "): ", detail,
    "; per-cluster means and spatial patterns for these clusters are unstable and low-confidence (references/pitfalls.md #2-3)."
  ))
}
if (n_clusters <= 1) {
  caveats <- c(caveats, "- Only one cluster was found at this resolution; the domain map and any domain-level comparison are uninformative for spatial heterogeneity. Consider raising `cluster_resolution` and re-checking against histology.")
}
if (!is.null(svg_table)) {
  caveats <- c(caveats, paste0(
    "- Spatially variable gene ranking (Moran's I on SCT-normalized data) is only depth-corrected to the extent SCTransform's regularized model absorbs the sequencing-depth gradient; it is not a formal depth-offset SVG test (SPARK-X) or FDR-controlled p-value. ",
    "Cross-check the top hits above against `spot_qc_map` (", qc_metric, ") before reporting them as true spatially variable genes, not depth artifacts (references/pitfalls.md #4)."
  ))
}
caveats <- c(caveats, "- Spatial clusters here are expression-driven (SCTransform + PCA + Louvain on the k-NN graph via FindNeighbors), not a spatial-prior method (BayesSpace/SpaGCN/STAGATE); cross-check against histology before treating clusters as anatomical domains, as boundaries can show salt-and-pepper noise (references/pitfalls.md #2-3).")
caveats <- c(caveats, "- A Visium spot mixes multiple cells; do not report a single cell-type identity per spot without reference-based deconvolution (references/pitfalls.md #1).")
caveats <- c(caveats, "- Coordinate/image alignment was not re-verified here beyond what Load10X_Spatial/the saved object already encodes; treat any marker-vs-histology mismatch as a possible alignment error (references/pitfalls.md #6).")

## ---- plain-language conclusion ----
conclusion <- paste0(
  "Of ", n_total_spots, " spots in the input object, ", n_qc_pass, " passed QC (", qc_metric, " > 0) and were used in this analysis, ",
  "forming ", n_clusters, " spatial cluster(s)/domain(s) at resolution ", resolution, ". ",
  if (length(present_genes) > 0) paste0("Spatial feature maps were generated for ", paste(present_genes, collapse = ", "), ". ") else "",
  if (!is.null(svg_table)) {
    paste0("The top candidate spatially variable gene by Moran's I was ", top_svg_genes[1], " (see results/spatially_variable_genes.tsv for the full ranking). ")
  } else {
    "Spatially variable gene detection did not complete in this run. "
  },
  "These clusters and gene rankings are exploratory and expression/graph-driven; confirm them against histology and the caveats below before drawing anatomical or biological conclusions."
)

report_lines <- c(
  "# Spatial Seurat report",
  "",
  "## Input summary",
  paste0("- Input: `", input, "`"),
  paste0("- Assay: ", assay_name, ", total spots: ", n_total_spots, ", spots retained after QC (", qc_metric, " > 0): ", n_qc_pass),
  "",
  "## QC observations",
  paste0("- QC metric plotted (`spot_qc_map`): ", qc_metric),
  if (length(mt_genes) > 0) {
    paste0("- Mitochondrial genes detected and `percent.mt` computed (", length(mt_genes), " genes).")
  } else {
    "- No mitochondrial gene pattern (^MT-/^mt-) detected in this feature panel; percent.mt not computed."
  },
  "",
  "## Target gene and domain maps",
  paste0("- Target genes requested: ", paste(target_genes, collapse = ", ")),
  if (length(absent_genes) > 0) {
    paste0("- Genes not found and skipped: ", paste(absent_genes, collapse = ", "))
  } else {
    "- All requested genes found in the object."
  },
  paste0("- Spatial clusters (domains) at resolution ", resolution, ": ", n_clusters, " (sizes: ",
         paste0(domain_summary$cluster, "=", domain_summary$n_spots, collapse = ", "), ")"),
  "",
  "## Spatially variable genes",
  if (!is.null(svg_table)) {
    paste0("- Top ", top_n_svg, " spatially variable gene(s) (Moran's I), ranked over ", nrow(svg_table),
           " tested variable genes: ", paste(top_svg_genes, collapse = ", "), "; see results/spatially_variable_genes.tsv.")
  } else {
    paste0("- Spatially variable gene detection did not complete: ", svg_error)
  },
  "",
  "## Deconvolution / neighborhood",
  "- Not computed in this run (no single-cell reference or deconvolution table supplied).",
  "",
  "## Conclusion",
  conclusion,
  "",
  "## Caveats specific to this run",
  caveats
)
writeLines(report_lines, file.path(out_dir, "report.md"))

message("Wrote spatial Seurat outputs to ", out_dir)
