#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript run_cellchat.R input_seurat.rds celltype_col output_dir species[human|mouse] [group_col]")
}
input <- args[[1]]
celltype_col <- args[[2]]
out_dir <- args[[3]]
species <- args[[4]]
group_col <- ifelse(length(args) >= 5 && nzchar(args[[5]]) && args[[5]] != "NA", args[[5]], NA)

required <- c("Seurat", "CellChat", "patchwork", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

if (!file.exists(input)) stop("Input file not found: ", input)
if (!tolower(species) %in% c("human", "mouse")) stop("species must be 'human' or 'mouse', got: ", species)

dir.create(file.path(out_dir, "results", "cellchat"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures", "cellchat"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "run"), recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input)
if (!methods::is(obj, "Seurat")) stop("Input object is not a Seurat object: ", input)
meta <- obj@meta.data
if (!celltype_col %in% colnames(meta)) stop("celltype_col '", celltype_col, "' not found in object metadata.")
if (!is.na(group_col) && !(group_col %in% colnames(meta))) {
  stop("group_col '", group_col, "' not found in object metadata.")
}
if (length(unique(meta[[celltype_col]])) < 2) {
  stop("Need at least 2 cell groups in '", celltype_col, "' for cell-cell communication analysis.")
}

run_cellchat_pipeline <- function(seurat_obj, celltype_col, species) {
  # SeuratObject >= 5.0.0 made GetAssayData()'s `slot=` argument defunct (hard error) in favor
  # of `layer=`; older SeuratObject only understands `slot=`. Try the modern name first, fall
  # back to the old one, so this runs across both the Seurat 4.x and 5.x install base.
  data.input <- tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = Seurat::DefaultAssay(seurat_obj), layer = "data"),
    error = function(e) Seurat::GetAssayData(seurat_obj, assay = Seurat::DefaultAssay(seurat_obj), slot = "data")
  )
  cellchat <- CellChat::createCellChat(object = data.input, meta = seurat_obj@meta.data, group.by = celltype_col)
  cellchat@DB <- if (tolower(species) == "mouse") CellChat::CellChatDB.mouse else CellChat::CellChatDB.human
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  cellchat
}

## ---- Main (pooled) CellChat run ----
cellchat <- run_cellchat_pipeline(obj, celltype_col, species)
comm <- CellChat::subsetCommunication(cellchat)
write.table(comm, file.path(out_dir, "results", "cellchat", "ligand_receptor_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

if (length(cellchat@netP$pathways)) {
  pathway_summary <- data.frame(pathway = cellchat@netP$pathways)
  write.table(pathway_summary, file.path(out_dir, "results", "cellchat", "pathway_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}

# Required figure: interaction network (circle plot).
png(file.path(out_dir, "figures", "cellchat", "interaction_circle.png"), width = 1800, height = 1600, res = 250)
CellChat::netVisual_circle(cellchat@net$count, vertex.weight = as.numeric(table(meta[[celltype_col]])), weight.scale = TRUE, label.edge = FALSE)
dev.off()
pdf(file.path(out_dir, "figures", "cellchat", "interaction_circle.pdf"), width = 8, height = 7)
CellChat::netVisual_circle(cellchat@net$count, vertex.weight = as.numeric(table(meta[[celltype_col]])), weight.scale = TRUE, label.edge = FALSE)
dev.off()

# Required figure: ligand-receptor bubble plot (communication "bubble" per cns-method-guide).
if (nrow(comm) > 0) {
  bubble_plot <- CellChat::netVisual_bubble(cellchat, remove.isolate = FALSE)
  ggplot2::ggsave(file.path(out_dir, "figures", "cellchat", "interaction_bubble.png"), bubble_plot, width = max(8, 0.28 * nrow(comm) + 4), height = max(6, 0.22 * length(unique(paste(comm$source, comm$target))) + 3), dpi = 320, bg = "white")
  ggplot2::ggsave(file.path(out_dir, "figures", "cellchat", "interaction_bubble.pdf"), bubble_plot, width = max(8, 0.28 * nrow(comm) + 4), height = max(6, 0.22 * length(unique(paste(comm$source, comm$target))) + 3), bg = "white")
}

saveRDS(cellchat, file.path(out_dir, "results", "cellchat", "cellchat_object.rds"))

# ---- Composition / low-n caveat (references/pitfalls.md #11): CellChat's own filterCommunication
# already drops interactions below min.cells=10, but even 10-30 cells gives an unstable per-group
# communication probability estimate, so flag small groups explicitly in the report. ----
celltype_counts_tbl <- table(meta[[celltype_col]])
small_celltype_groups <- celltype_counts_tbl[celltype_counts_tbl < 20]

# ---- Top senders/receivers/pathways for the report text ----
top_senders <- if (nrow(comm) > 0) sort(table(comm$source), decreasing = TRUE) else integer(0)
top_receivers <- if (nrow(comm) > 0) sort(table(comm$target), decreasing = TRUE) else integer(0)
pathway_col <- if ("pathway_name" %in% colnames(comm)) "pathway_name" else NA
top_pathways <- if (!is.na(pathway_col) && nrow(comm) > 0) sort(table(comm[[pathway_col]]), decreasing = TRUE) else integer(0)

report_lines <- c(
  "# CellChat communication report",
  "",
  "## Results",
  "",
  paste0("- Species: ", species),
  paste0("- Cell groups: ", length(unique(meta[[celltype_col]]))),
  paste0("- Interactions: ", nrow(comm)),
  paste0("- Pathways: ", length(cellchat@netP$pathways))
)
if (nrow(comm) > 0) {
  report_lines <- c(
    report_lines,
    paste0(
      "- Top sending cell type(s): ",
      paste(paste0(names(head(top_senders, 3)), " (", as.integer(head(top_senders, 3)), ")"), collapse = ", ")
    ),
    paste0(
      "- Top receiving cell type(s): ",
      paste(paste0(names(head(top_receivers, 3)), " (", as.integer(head(top_receivers, 3)), ")"), collapse = ", ")
    )
  )
  if (length(top_pathways)) {
    report_lines <- c(
      report_lines,
      paste0(
        "- Top signaling pathway(s): ",
        paste(paste0(names(head(top_pathways, 5)), " (", as.integer(head(top_pathways, 5)), ")"), collapse = ", ")
      )
    )
  }
} else {
  report_lines <- c(report_lines, "- No significant ligand-receptor interactions passed CellChat's default filtering (min.cells = 10).")
}

caveats <- character(0)
if (length(small_celltype_groups)) {
  caveats <- c(caveats, paste0(
    "The following cell group(s) have fewer than 20 cells: ",
    paste0(names(small_celltype_groups), " (n=", as.integer(small_celltype_groups), ")", collapse = "; "),
    " -- communication probabilities for these groups are statistically unstable (references/pitfalls.md #11)."
  ))
}
caveats <- c(caveats, "Communication is inferred from expression (CellChatDB over-representation test), not directly observed physical signaling; report as inference, not proof.")
report_lines <- c(report_lines, "", "### Caveats", "", paste0("- ", caveats))

## ---- Optional group comparison (only when group_col has >=2 usable levels) ----
if (!is.na(group_col)) {
  group_levels <- names(table(meta[[group_col]]))[table(meta[[group_col]]) >= 30]
  if (length(group_levels) < 2) {
    message(
      "Skipping group comparison: fewer than 2 levels of '", group_col,
      "' have at least 30 cells (CellChat needs enough cells per group to estimate communication probabilities reliably)."
    )
    report_lines <- c(report_lines, "", "### Group comparison", "", "Skipped: insufficient cells per group level.")
  } else {
    group_levels <- group_levels[1:2]
    cc_list <- list()
    for (lvl in group_levels) {
      sub_obj <- subset(obj, cells = rownames(meta)[meta[[group_col]] == lvl])
      cc_list[[lvl]] <- run_cellchat_pipeline(sub_obj, celltype_col, species)
    }
    merged <- CellChat::mergeCellChat(cc_list, add.names = names(cc_list))

    compare_plot <- CellChat::compareInteractions(merged, show.legend = FALSE, group = c(1, 2))
    ggplot2::ggsave(file.path(out_dir, "figures", "cellchat", "interaction_comparison.png"), compare_plot, width = 6, height = 5.4, dpi = 320, bg = "white")
    ggplot2::ggsave(file.path(out_dir, "figures", "cellchat", "interaction_comparison.pdf"), compare_plot, width = 6, height = 5.4, bg = "white")

    rank_plot <- CellChat::rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
    ggplot2::ggsave(file.path(out_dir, "figures", "cellchat", "pathway_rank_comparison.png"), rank_plot, width = 7.5, height = 7, dpi = 320, bg = "white")
    ggplot2::ggsave(file.path(out_dir, "figures", "cellchat", "pathway_rank_comparison.pdf"), rank_plot, width = 7.5, height = 7, bg = "white")

    png(file.path(out_dir, "figures", "cellchat", "interaction_diff_circle.png"), width = 1800, height = 1600, res = 250)
    CellChat::netVisual_diffInteraction(merged, weight.scale = TRUE)
    dev.off()
    pdf(file.path(out_dir, "figures", "cellchat", "interaction_diff_circle.pdf"), width = 8, height = 7)
    CellChat::netVisual_diffInteraction(merged, weight.scale = TRUE)
    dev.off()

    comparison_counts <- data.frame(
      group = names(cc_list),
      interactions = vapply(cc_list, function(cc) nrow(CellChat::subsetCommunication(cc)), numeric(1)),
      cells = vapply(group_levels, function(lvl) sum(meta[[group_col]] == lvl), numeric(1))
    )
    write.table(comparison_counts, file.path(out_dir, "results", "cellchat", "group_comparison_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    saveRDS(merged, file.path(out_dir, "results", "cellchat", "cellchat_merged_object.rds"))

    comparison_lines <- paste0(
      "  - ", comparison_counts$group, ": ", comparison_counts$interactions, " interaction(s), ",
      comparison_counts$cells, " cell(s)"
    )
    report_lines <- c(
      report_lines, "",
      "### Group comparison", "",
      paste0("Compared `", group_col, "`: ", paste(group_levels, collapse = " vs ")),
      comparison_lines
    )
  }
}

conclusion <- paste0(
  "Conclusion: across ", length(unique(meta[[celltype_col]])), " cell group(s), CellChat inferred ", nrow(comm),
  " significant ligand-receptor interaction(s) spanning ", length(cellchat@netP$pathways), " pathway(s).",
  if (nrow(comm) > 0) paste0(
    " Signaling is most concentrated from ", names(head(top_senders, 1)), " towards ", names(head(top_receivers, 1)), "."
  ) else "",
  if (length(small_celltype_groups)) " Some cell group(s) had low cell counts (see caveats); interpret their signaling with caution." else "",
  " These results describe expression-based interaction inference, not confirmed physical ligand-receptor signaling; validate key pairs experimentally before drawing mechanistic conclusions."
)
report_lines <- c(report_lines, "", "### Conclusion", "", conclusion)

writeLines(report_lines, file.path(out_dir, "report_cellchat.md"))
