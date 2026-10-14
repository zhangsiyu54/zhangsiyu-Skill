#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript run_trajectory_monocle3.R input_seurat.rds celltype_col output_dir root_celltype")
}
input <- args[[1]]
celltype_col <- args[[2]]
out_dir <- args[[3]]
root_celltype <- args[[4]]

# NOTE: as.cell_data_set() is exported by the SeuratWrappers package (satijalab/seurat-wrappers),
# not by core Seurat. Both packages are required to convert a Seurat object into a monocle3
# cell_data_set.
required <- c("Seurat", "SeuratWrappers", "monocle3", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

if (!file.exists(input)) stop("Input file not found: ", input)

dir.create(file.path(out_dir, "results", "trajectory"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures", "trajectory"), recursive = TRUE, showWarnings = FALSE)

TEXT_COL <- "#173B5F"

obj <- readRDS(input)
if (!methods::is(obj, "Seurat")) stop("Input object is not a Seurat object: ", input)
meta <- obj@meta.data
if (!celltype_col %in% colnames(meta)) stop("celltype_col '", celltype_col, "' not found in object metadata.")
if (!"umap" %in% tolower(names(obj@reductions))) {
  stop("Object has no UMAP reduction. Run Seurat::RunUMAP() (e.g. via run_seurat.R) before trajectory inference.")
}
root_cells <- rownames(meta)[meta[[celltype_col]] == root_celltype]
if (!length(root_cells)) stop("Root celltype '", root_celltype, "' not found in metadata column '", celltype_col, "'.")

## ---- Convert to monocle3 cell_data_set and learn the trajectory graph ----
cds <- SeuratWrappers::as.cell_data_set(obj)
cds <- monocle3::cluster_cells(cds)
cds <- monocle3::learn_graph(cds, use_partition = TRUE)
cds <- monocle3::order_cells(cds, root_cells = root_cells)
pt <- monocle3::pseudotime(cds)

pt_tbl <- data.frame(cell = names(pt), pseudotime = as.numeric(pt), celltype = meta[names(pt), celltype_col])
# Cells assigned Inf pseudotime are unreachable from the root in the principal graph (disconnected
# partition); keep them in the table but flag them rather than silently dropping (pitfalls.md #8).
pt_tbl$reachable_from_root <- is.finite(pt_tbl$pseudotime)
write.table(pt_tbl, file.path(out_dir, "results", "trajectory", "pseudotime.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## ---- Required figure: trajectory embedding colored by pseudotime ----
p <- monocle3::plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = TRUE) +
  ggplot2::labs(title = "Trajectory pseudotime") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20, color = TEXT_COL))
ggplot2::ggsave(file.path(out_dir, "figures", "trajectory", "pseudotime.png"), p, width = 8, height = 6.4, dpi = 320, bg = "white")
ggplot2::ggsave(file.path(out_dir, "figures", "trajectory", "pseudotime.pdf"), p, width = 8, height = 6.4, bg = "white")

## ---- Required (advanced-modules.md): pseudotime distribution by cell type/group ----
finite_tbl <- pt_tbl[pt_tbl$reachable_from_root, ]
dist_plot <- ggplot2::ggplot(finite_tbl, ggplot2::aes(x = celltype, y = pseudotime, fill = celltype)) +
  ggplot2::geom_violin(alpha = 0.65, trim = TRUE) +
  ggplot2::geom_jitter(width = 0.12, size = 0.6, alpha = 0.35, color = TEXT_COL) +
  ggplot2::labs(title = "Pseudotime by cell type", x = NULL, y = "Pseudotime") +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20, color = TEXT_COL),
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
  )
ggplot2::ggsave(file.path(out_dir, "figures", "trajectory", "pseudotime_by_celltype.png"), dist_plot, width = 7.6, height = 6.0, dpi = 320, bg = "white")
ggplot2::ggsave(file.path(out_dir, "figures", "trajectory", "pseudotime_by_celltype.pdf"), dist_plot, width = 7.6, height = 6.0, bg = "white")

## ---- Required (advanced-modules.md): branch/pseudotime-associated genes ----
graph_res <- tryCatch(
  monocle3::graph_test(cds, neighbor_graph = "principal_graph", cores = 1),
  error = function(e) {
    message("monocle3::graph_test failed, skipping pseudotime-gene table: ", conditionMessage(e))
    NULL
  }
)
if (!is.null(graph_res)) {
  graph_res <- graph_res[order(graph_res$q_value), ]
  graph_res$gene <- rownames(graph_res)
  write.table(
    graph_res[, c("gene", "status", "p_value", "q_value", "morans_I")],
    file.path(out_dir, "results", "trajectory", "pseudotime_genes.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

saveRDS(cds, file.path(out_dir, "results", "trajectory", "monocle3_cds.rds"))

n_unreachable <- sum(!pt_tbl$reachable_from_root)
n_sig_genes <- if (!is.null(graph_res)) sum(graph_res$q_value < 0.05, na.rm = TRUE) else NA
top_pt_genes <- if (!is.null(graph_res)) head(graph_res$gene, 5) else character(0)

# ---- Cell type composition along pseudotime + low-n caveat (references/pitfalls.md #11) ----
celltype_counts_tbl <- table(meta[[celltype_col]])
small_celltype_groups <- celltype_counts_tbl[celltype_counts_tbl < 20]

caveats <- character(0)
if (length(root_cells) < 20) {
  caveats <- c(caveats, paste0(
    "Root population '", root_celltype, "' has only ", length(root_cells), " cells; the pseudotime root and ",
    "downstream ordering may be unstable with this few cells (references/pitfalls.md #11)."
  ))
}
if (length(small_celltype_groups)) {
  caveats <- c(caveats, paste0(
    "The following cell type(s) along the trajectory have fewer than 20 cells: ",
    paste0(names(small_celltype_groups), " (n=", as.integer(small_celltype_groups), ")", collapse = "; "),
    " -- their pseudotime distribution should be treated as exploratory (references/pitfalls.md #11)."
  ))
}
caveats <- c(caveats, paste0(
  "Pseudotime is only meaningful for a biologically plausible, continuous lineage; verify the root and subset ",
  "before interpreting (references/pitfalls.md #8)."
))

conclusion <- paste0(
  "Conclusion: starting from ", length(root_cells), " '", root_celltype, "' root cell(s), ",
  sum(pt_tbl$reachable_from_root), " of ", nrow(pt_tbl), " cells were ordered along the trajectory",
  if (n_unreachable > 0) paste0(" (", n_unreachable, " cell(s) were unreachable/disconnected and excluded)") else "",
  ". ",
  if (!is.na(n_sig_genes)) paste0(n_sig_genes, " gene(s) were significantly associated with the trajectory (q < 0.05)",
    if (length(top_pt_genes)) paste0(", led by ", paste(top_pt_genes, collapse = ", "), ". ") else ". ") else "Pseudotime-associated gene testing was not available. ",
  "Treat this trajectory as descriptive of the ordering within this dataset; see caveats below before drawing developmental conclusions."
)

report_lines <- c(
  "# Trajectory report",
  "",
  "## Results",
  "",
  paste0("- Root cell type: ", root_celltype, " (", length(root_cells), " root cells)"),
  paste0("- Cells with finite pseudotime: ", sum(pt_tbl$reachable_from_root), " / ", nrow(pt_tbl)),
  if (n_unreachable > 0) paste0("- WARNING: ", n_unreachable, " cells are unreachable from the root (disconnected partition) and were excluded from the distribution plot.") else "- All cells reachable from the root.",
  if (!is.null(graph_res)) paste0("- Pseudotime-associated genes (q < 0.05): ", n_sig_genes) else "- Pseudotime-associated gene test not available."
)
if (length(top_pt_genes)) {
  report_lines <- c(report_lines, paste0("- Top pseudotime-associated genes (lowest q-value): ", paste(top_pt_genes, collapse = ", ")))
}
report_lines <- c(report_lines, "", "### Caveats", "", paste0("- ", caveats))
report_lines <- c(report_lines, "", "### Conclusion", "", conclusion)
writeLines(report_lines, file.path(out_dir, "report_trajectory.md"))
