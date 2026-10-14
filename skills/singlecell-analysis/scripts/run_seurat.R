#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop(
    "Usage: Rscript run_seurat.R input.rds output_dir [group_col] [celltype_col] [genes_csv] ",
    "[resolution] [npcs] [min_features] [max_mt_pct]"
  )
}
input <- args[[1]]
out_dir <- args[[2]]
group_col <- ifelse(length(args) >= 3 && nzchar(args[[3]]) && args[[3]] != "NA", args[[3]], NA)
celltype_col <- ifelse(length(args) >= 4 && nzchar(args[[4]]) && args[[4]] != "NA", args[[4]], NA)
genes_csv <- ifelse(length(args) >= 5 && nzchar(args[[5]]) && args[[5]] != "NA", args[[5]], NA)
resolution <- ifelse(length(args) >= 6, as.numeric(args[[6]]), 0.8)
npcs <- ifelse(length(args) >= 7, as.integer(args[[7]]), 20L)
min_features <- ifelse(length(args) >= 8, as.numeric(args[[8]]), 200)
max_mt_pct <- ifelse(length(args) >= 9, as.numeric(args[[9]]), 15)

required <- c("Seurat", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

if (!file.exists(input)) stop("Input file not found: ", input)

dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "run"), recursive = TRUE, showWarnings = FALSE)

# Color palette matches assets/palette.json (clusters / feature_low / feature_high).
CLUSTER_COLORS <- c("#3B82C4", "#D84A4A", "#3FA37A", "#D99A2B", "#7C65C1", "#4DA6A6")
FEATURE_LOW <- "#E5E7EB"
FEATURE_HIGH <- "#2563A8"
TEXT_COL <- "#1F2937"

cluster_colors <- function(n) rep(CLUSTER_COLORS, length.out = n)

save_gg <- function(plot, name, subdir = "figures", width = 8, height = 6) {
  dir.create(file.path(out_dir, subdir), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(out_dir, subdir, paste0(name, ".png")), plot, width = width, height = height, dpi = 320, bg = "white")
  ggplot2::ggsave(file.path(out_dir, subdir, paste0(name, ".pdf")), plot, width = width, height = height, bg = "white")
}

publication_theme <- function() {
  ggplot2::theme_classic(base_size = 15) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20, color = TEXT_COL),
      axis.title = ggplot2::element_text(size = 17, color = TEXT_COL),
      axis.text = ggplot2::element_text(size = 13, color = TEXT_COL),
      legend.text = ggplot2::element_text(size = 12, color = TEXT_COL)
    )
}

write_json_simple <- function(x, path) {
  # Minimal dependency-free JSON writer for flat named lists / vectors (no jsonlite requirement).
  kv <- function(name, value) {
    if (is.numeric(value)) {
      paste0("\"", name, "\": ", value)
    } else if (is.logical(value)) {
      paste0("\"", name, "\": ", tolower(as.character(value)))
    } else {
      paste0("\"", name, "\": \"", value, "\"")
    }
  }
  lines <- mapply(kv, names(x), x)
  writeLines(c("{", paste0("  ", paste(lines, collapse = ",\n  ")), "}"), path)
}

obj <- readRDS(input)
if (!methods::is(obj, "Seurat")) stop("Input object is not a Seurat object: ", input)
if (!is.na(group_col) && !(group_col %in% colnames(obj@meta.data))) {
  stop("group_col '", group_col, "' not found in object metadata.")
}
if (!is.na(celltype_col) && !(celltype_col %in% colnames(obj@meta.data))) {
  stop("celltype_col '", celltype_col, "' not found in object metadata.")
}

target_genes <- if (!is.na(genes_csv)) trimws(strsplit(genes_csv, ",")[[1]]) else character(0)

## ---- QC ----
obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
qc_plot <- Seurat::VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
  publication_theme()
save_gg(qc_plot, "qc_violin", width = 10.2, height = 5.2)

cells_before <- ncol(obj)
# Do not remove cells silently (references/pitfalls.md): record thresholds and before/after counts.
keep <- obj$nFeature_RNA >= min_features & obj$percent.mt <= max_mt_pct
obj <- subset(obj, cells = colnames(obj)[keep])
cells_after <- ncol(obj)

## ---- Normalization, dimensionality reduction, clustering ----
obj <- Seurat::NormalizeData(obj)
obj <- Seurat::FindVariableFeatures(obj)
obj <- Seurat::ScaleData(obj)
obj <- Seurat::RunPCA(obj, npcs = max(npcs, 2))
obj <- Seurat::FindNeighbors(obj, dims = 1:npcs)
obj <- Seurat::FindClusters(obj, resolution = resolution)
obj <- Seurat::RunUMAP(obj, dims = 1:npcs)

umap_plot <- Seurat::DimPlot(obj, reduction = "umap", label = TRUE, cols = cluster_colors(nlevels(Seurat::Idents(obj)))) +
  ggplot2::labs(title = "UMAP clusters", x = "UMAP 1", y = "UMAP 2") +
  publication_theme() +
  ggplot2::theme(legend.position = "right")
save_gg(umap_plot, "umap_clusters", width = 7.2, height = 6.4)

if (!is.na(celltype_col)) {
  umap_celltype <- Seurat::DimPlot(obj, reduction = "umap", group.by = celltype_col, label = TRUE) +
    ggplot2::labs(title = "UMAP cell type", x = "UMAP 1", y = "UMAP 2") +
    publication_theme() +
    ggplot2::theme(legend.position = "right")
  save_gg(umap_celltype, "umap_celltype", width = 7.6, height = 6.4)
}

if (!is.na(group_col)) {
  umap_condition <- Seurat::DimPlot(obj, reduction = "umap", group.by = group_col) +
    ggplot2::labs(title = "UMAP by condition", x = "UMAP 1", y = "UMAP 2") +
    publication_theme() +
    ggplot2::theme(legend.position = "right")
  save_gg(umap_condition, "umap_condition", width = 7.6, height = 6.4)
}

## ---- Markers and top-marker heatmap ----
markers <- Seurat::FindAllMarkers(obj, only.pos = TRUE)
write.table(markers, file.path(out_dir, "results", "marker_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

if (nrow(markers) > 0) {
  top_markers <- do.call(rbind, lapply(split(markers, markers$cluster), function(df) {
    df <- df[order(-df$avg_log2FC), ]
    head(df, 10)
  }))
  heatmap_plot <- Seurat::DoHeatmap(
    obj,
    features = unique(top_markers$gene),
    group.colors = cluster_colors(nlevels(Seurat::Idents(obj)))
  ) +
    ggplot2::labs(title = "Top marker heatmap") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20, color = TEXT_COL),
      legend.text = ggplot2::element_text(size = 11)
    )
  save_gg(heatmap_plot, "heatmap_markers", width = 10, height = max(6, 0.18 * length(unique(top_markers$gene)) + 3))
}

## ---- Cell composition (proportion per cluster, split by condition when available) ----
composition_group_col <- if (!is.na(group_col)) group_col else NULL
composition_celltype_col <- if (!is.na(celltype_col)) celltype_col else "seurat_clusters"
meta <- obj@meta.data
if (!is.null(composition_group_col)) {
  prop_tbl <- as.data.frame(table(meta[[composition_group_col]], meta[[composition_celltype_col]]))
  colnames(prop_tbl) <- c("group", "celltype", "count")
} else {
  prop_tbl <- as.data.frame(table(meta[[composition_celltype_col]]))
  colnames(prop_tbl) <- c("celltype", "count")
  prop_tbl$group <- "all"
  prop_tbl <- prop_tbl[, c("group", "celltype", "count")]
}
write.table(prop_tbl, file.path(out_dir, "results", "celltype_proportion.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

composition_plot <- ggplot2::ggplot(prop_tbl, ggplot2::aes(x = group, y = count, fill = celltype)) +
  ggplot2::geom_bar(stat = "identity", position = "fill", width = 0.62) +
  ggplot2::scale_fill_manual(values = cluster_colors(length(unique(prop_tbl$celltype)))) +
  ggplot2::labs(title = "Cell type proportion", x = "Group", y = "Proportion", fill = NULL) +
  publication_theme()
save_gg(composition_plot, "celltype_proportion", width = 7.4, height = 6.0)

## ---- FeaturePlot for requested target genes ----
present_genes <- intersect(target_genes, rownames(obj))
missing_genes <- setdiff(target_genes, present_genes)
if (length(missing_genes)) {
  message("Requested genes not found in object and skipped for FeaturePlot: ", paste(missing_genes, collapse = ", "))
}
if (length(present_genes)) {
  feature_plot <- Seurat::FeaturePlot(obj, features = present_genes) &
    ggplot2::scale_color_gradient(low = FEATURE_LOW, high = FEATURE_HIGH) &
    publication_theme()
  save_gg(feature_plot, "featureplot", width = max(7, 3.6 * min(length(present_genes), 3)), height = 6.2 * ceiling(length(present_genes) / 3))
}

## ---- Module score (only when a gene list is provided via genes_csv) ----
if (length(target_genes) >= 1) {
  obj <- Seurat::AddModuleScore(obj, features = list(target_genes), name = "ModuleScore")
  score_col <- grep("^ModuleScore", colnames(obj@meta.data), value = TRUE)[1]
  module_umap <- Seurat::FeaturePlot(obj, features = score_col) +
    ggplot2::scale_color_gradient(low = FEATURE_LOW, high = FEATURE_HIGH) +
    ggplot2::labs(title = "Module score") +
    publication_theme()
  save_gg(module_umap, "module_score_umap", width = 7.6, height = 6.4)

  module_group <- if (!is.na(celltype_col)) celltype_col else "seurat_clusters"
  module_group_n <- length(unique(obj@meta.data[[module_group]]))
  module_violin <- Seurat::VlnPlot(obj, features = score_col, group.by = module_group, cols = cluster_colors(module_group_n)) +
    ggplot2::labs(title = "Module score by group", y = "Module score") +
    publication_theme()
  save_gg(module_violin, "module_score_violin", width = 7.6, height = 6.0)

  module_summary <- stats::aggregate(obj@meta.data[[score_col]], by = list(group = obj@meta.data[[module_group]]), FUN = mean)
  colnames(module_summary) <- c("group", "mean_module_score")
  write.table(module_summary, file.path(out_dir, "results", "module_score_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}

## ---- Reporting ----
write_json_simple(
  list(
    cells_before_qc = cells_before,
    cells_after_qc = cells_after,
    min_features = min_features,
    max_mt_pct = max_mt_pct,
    resolution = resolution,
    npcs = npcs
  ),
  file.path(out_dir, "run", "qc.json")
)

n_clusters <- length(unique(Seurat::Idents(obj)))

# ---- Composition breakdown: top groups + low-n groups (references/pitfalls.md #11) ----
group_col_for_counts <- meta[[composition_celltype_col]]
group_counts_tbl <- sort(table(group_col_for_counts), decreasing = TRUE)
total_grouped_cells <- sum(group_counts_tbl)
composition_report_lines <- paste0(
  "  - ", names(group_counts_tbl), ": ", as.integer(group_counts_tbl), " cells (",
  sprintf("%.1f%%", 100 * as.numeric(group_counts_tbl) / max(1, total_grouped_cells)), ")"
)
small_group_counts <- group_counts_tbl[group_counts_tbl < 20]

# ---- Top marker genes per cluster (by avg_log2FC) for the report text ----
marker_report_lines <- character(0)
if (nrow(markers) > 0) {
  top3_markers <- do.call(rbind, lapply(split(markers, markers$cluster), function(df) {
    df <- df[order(-df$avg_log2FC), ]
    head(df, 3)
  }))
  for (cl in unique(as.character(top3_markers$cluster))) {
    genes_cl <- top3_markers$gene[as.character(top3_markers$cluster) == cl]
    marker_report_lines <- c(marker_report_lines, paste0("  - Cluster ", cl, ": ", paste(genes_cl, collapse = ", ")))
  }
}

# ---- Caveats tied to references/pitfalls.md, only emitted when they actually apply ----
caveats <- character(0)
if (length(small_group_counts)) {
  caveats <- c(caveats, paste0(
    "Low cell count in '", composition_celltype_col, "': ",
    paste0(names(small_group_counts), " (n=", as.integer(small_group_counts), ")", collapse = "; "),
    " -- annotation/statistics for these groups should be treated as unstable, not confirmatory (references/pitfalls.md #11)."
  ))
}
if (length(target_genes) == 1) {
  caveats <- c(caveats, paste0(
    "Module score / FeaturePlot interpretation relies on a single gene ('", target_genes, "'); a single marker ",
    "cannot confirm cell identity or pathway activity on its own -- treat as exploratory (references/pitfalls.md #6 and #14)."
  ))
}
if (!is.na(group_col)) {
  caveats <- c(caveats, paste0(
    "Cell type composition is pooled across all cells for '", group_col, "'; a statistically valid per-condition ",
    "comparison requires per-sample proportions and a sample-aware test, not a pooled cell-level count ",
    "(references/pitfalls.md #12). This script was not given a sample ID column, so only pooled counts are reported."
  ))
}

# ---- Plain-language conclusion synthesizing composition + markers + caveats ----
top_group_names <- head(names(group_counts_tbl), 3)
conclusion <- paste0(
  "In summary, ", cells_after, " of ", cells_before, " input cells passed QC (nFeature_RNA >= ", min_features,
  ", percent.mt <= ", max_mt_pct, ") and were organized into ", n_clusters, " cluster(s) at resolution ", resolution,
  " using ", npcs, " PCs",
  if (!is.na(celltype_col)) paste0(", annotated via '", celltype_col, "'") else " (clusters only; no cell type annotation column was provided)",
  ". The most abundant group(s) were ", paste(top_group_names, collapse = ", "), ".",
  if (length(marker_report_lines)) " Cluster identity is supported by distinct top marker genes per cluster (see below)." else " No positive marker genes were found for any cluster.",
  if (length(caveats)) " Review the caveats below before treating these results as confirmatory." else " No automatic caveats were triggered for this input."
)

report_lines <- c(
  "# Single-cell Seurat report",
  "",
  "## Results",
  "",
  paste0("- Cells before QC filtering: ", cells_before),
  paste0("- Cells after QC filtering (nFeature_RNA >= ", min_features, ", percent.mt <= ", max_mt_pct, "): ", cells_after),
  paste0("- Clusters (resolution ", resolution, ", ", npcs, " PCs): ", n_clusters),
  paste0("- Marker genes reported (FindAllMarkers, only.pos=TRUE): ", nrow(markers)),
  if (!is.na(celltype_col)) paste0("- Cell type column used: ", celltype_col) else "- Cell type column: not provided (clusters only, needs_annotation)",
  if (!is.na(group_col)) paste0("- Condition/group column used: ", group_col) else "- Condition/group column: not provided",
  if (length(present_genes)) paste0("- FeaturePlot / module score genes: ", paste(present_genes, collapse = ", ")) else "- No target genes requested",
  "",
  paste0("### Cell type/cluster composition (", composition_celltype_col, ")"),
  "",
  composition_report_lines
)
if (length(marker_report_lines)) {
  report_lines <- c(report_lines, "", "### Top marker genes per cluster (by avg_log2FC)", "", marker_report_lines)
}
if (length(caveats)) {
  report_lines <- c(report_lines, "", "### Caveats", "", paste0("- ", caveats))
}
report_lines <- c(report_lines, "", "### Conclusion", "", conclusion)
writeLines(report_lines, file.path(out_dir, "report.md"))
write.table(
  data.frame(item = c("cells_before_qc", "cells_after_qc", "clusters"), value = c(cells_before, cells_after, n_clusters)),
  file.path(out_dir, "summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
saveRDS(obj, file.path(out_dir, "results", "seurat_object_processed.rds"))
