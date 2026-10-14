#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(
    "Usage: Rscript run_pseudobulk_de.R input_seurat.rds sample_col celltype_col condition_col output_dir ",
    "[case] [control] [padj_threshold=0.05] [log2fc_threshold=1]"
  )
}
input <- args[[1]]
sample_col <- args[[2]]
celltype_col <- args[[3]]
condition_col <- args[[4]]
out_dir <- args[[5]]
case_label <- ifelse(length(args) >= 6 && nzchar(args[[6]]) && args[[6]] != "NA", args[[6]], NA)
control_label <- ifelse(length(args) >= 7 && nzchar(args[[7]]) && args[[7]] != "NA", args[[7]], NA)
padj_threshold <- ifelse(length(args) >= 8 && nzchar(args[[8]]) && args[[8]] != "NA", as.numeric(args[[8]]), 0.05)
log2fc_threshold <- ifelse(length(args) >= 9 && nzchar(args[[9]]) && args[[9]] != "NA", as.numeric(args[[9]]), 1)
if (is.na(padj_threshold) || padj_threshold <= 0 || padj_threshold > 1) {
  stop("padj_threshold must be a number in (0, 1], got: ", if (length(args) >= 8) args[[8]] else NA)
}
if (is.na(log2fc_threshold) || log2fc_threshold < 0) {
  stop("log2fc_threshold must be a non-negative number, got: ", if (length(args) >= 9) args[[9]] else NA)
}

required <- c("Seurat", "Matrix", "DESeq2", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

dir.create(file.path(out_dir, "results", "pseudobulk"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures", "pseudobulk"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "run"), recursive = TRUE, showWarnings = FALSE)

save_gg <- function(plot, name, width = 8, height = 6) {
  ggplot2::ggsave(file.path(out_dir, "figures", "pseudobulk", paste0(name, ".png")), plot, width = width, height = height, dpi = 320, bg = "white")
  ggplot2::ggsave(file.path(out_dir, "figures", "pseudobulk", paste0(name, ".pdf")), plot, width = width, height = height, bg = "white")
}

obj <- readRDS(input)
meta <- obj@meta.data
stopifnot(sample_col %in% colnames(meta), celltype_col %in% colnames(meta), condition_col %in% colnames(meta))
## SeuratObject >= 5.0.0 made the `slot=` argument of GetAssayData() defunct (hard error, not
## just a deprecation warning) in favor of `layer=`; older SeuratObject (<5.0.0) does not accept
## `layer=` at all. Try the modern argument first and fall back to the old one so this script
## works across the Seurat 4.x/5.x split rather than hard-erroring on whichever version is
## actually installed.
counts <- tryCatch(
  Seurat::GetAssayData(obj, assay = Seurat::DefaultAssay(obj), layer = "counts"),
  error = function(e) Seurat::GetAssayData(obj, assay = Seurat::DefaultAssay(obj), slot = "counts")
)
meta$pb_group <- paste(meta[[sample_col]], meta[[celltype_col]], sep = "||")
groups <- split(seq_len(nrow(meta)), meta$pb_group)
pb_counts <- do.call(cbind, lapply(groups, function(idx) Matrix::rowSums(counts[, rownames(meta)[idx], drop = FALSE])))
pb_meta <- do.call(rbind, strsplit(colnames(pb_counts), "\\|\\|"))
pb_meta <- data.frame(pb_id = colnames(pb_counts), sample = pb_meta[, 1], celltype = pb_meta[, 2], stringsAsFactors = FALSE)
sample_meta <- unique(meta[, c(sample_col, condition_col), drop = FALSE])
colnames(sample_meta) <- c("sample", "condition")
pb_meta <- merge(pb_meta, sample_meta, by = "sample", all.x = TRUE)
rownames(pb_meta) <- pb_meta$pb_id
pb_counts <- pb_counts[, pb_meta$pb_id, drop = FALSE]

summary_rows <- list()
top_gene_lines <- character(0)
low_replicate_lines <- character(0)
skipped_celltypes <- character(0)
for (ct in sort(unique(pb_meta$celltype))) {
  keep <- pb_meta$celltype == ct
  sub_counts <- pb_counts[, keep, drop = FALSE]
  sub_meta <- pb_meta[keep, , drop = FALSE]
  if (length(unique(sub_meta$condition)) < 2) {
    skipped_celltypes <- c(skipped_celltypes, paste0(ct, " (only one condition level present)"))
    next
  }
  if (!is.na(case_label) && !is.na(control_label)) {
    sub_meta$condition <- factor(sub_meta$condition, levels = c(control_label, case_label))
  } else {
    sub_meta$condition <- factor(sub_meta$condition)
  }
  cond_counts <- table(sub_meta$condition)
  dds <- DESeq2::DESeqDataSetFromMatrix(round(as.matrix(sub_counts)), sub_meta, ~ condition)
  dds <- DESeq2::DESeq(dds)
  res <- as.data.frame(DESeq2::results(dds))
  res$gene <- rownames(res)
  res$celltype <- ct
  res$change <- ifelse(!is.na(res$padj) & res$padj < padj_threshold & res$log2FoldChange >= log2fc_threshold, "up",
    ifelse(!is.na(res$padj) & res$padj < padj_threshold & res$log2FoldChange <= -log2fc_threshold, "down", "ns"))
  res <- res[order(res$padj), c("celltype", "gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "change")]
  safe_ct <- gsub("[^A-Za-z0-9_-]+", "_", ct)
  write.table(res, file.path(out_dir, "results", "pseudobulk", paste0(safe_ct, "_DEG.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  volcano <- ggplot2::ggplot(res, ggplot2::aes(log2FoldChange, -log10(padj), color = change)) +
    ggplot2::geom_point(alpha = 0.82, size = 1.8, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(up = "#E45B5B", down = "#4EA3F1", ns = "#B8C2CC")) +
    ggplot2::labs(title = paste0("Pseudobulk DE: ", ct), x = "log2 fold change", y = "-log10 adjusted p value", color = NULL) +
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20, color = "#173B5F"))
  save_gg(volcano, paste0(safe_ct, "_volcano"), width = 8.2, height = 6.2)

  # Top-gene heatmap of pseudobulk normalized expression across samples (required alongside volcano).
  top_n <- min(30, sum(!is.na(res$padj) & res$change != "ns"))
  if (top_n > 1) {
    top_genes <- res$gene[order(res$padj)][seq_len(top_n)]
    norm_counts <- DESeq2::counts(dds, normalized = TRUE)[top_genes, , drop = FALSE]
    log_counts <- log2(norm_counts + 1)
    zscored <- t(scale(t(log_counts)))
    heat_df <- as.data.frame(as.table(zscored))
    colnames(heat_df) <- c("gene", "pb_sample", "zscore")
    heat_df$condition <- sub_meta[as.character(heat_df$pb_sample), "condition"]
    sample_order <- rownames(sub_meta)[order(sub_meta$condition)]
    heat_df$pb_sample <- factor(heat_df$pb_sample, levels = sample_order)
    heat_df$gene <- factor(heat_df$gene, levels = rev(top_genes))
    heatmap_plot <- ggplot2::ggplot(heat_df, ggplot2::aes(x = pb_sample, y = gene, fill = zscore)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.3) +
      ggplot2::scale_fill_gradient2(low = "#4EA3F1", mid = "white", high = "#E45B5B", midpoint = 0, na.value = "grey80") +
      ggplot2::labs(title = paste0("Pseudobulk top-gene heatmap: ", ct), x = "Pseudobulk sample", y = NULL, fill = "z-score") +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 18, color = "#173B5F"),
        axis.text.x = ggplot2::element_text(angle = 40, hjust = 1)
      )
    save_gg(heatmap_plot, paste0(safe_ct, "_heatmap"), width = max(6.5, 0.5 * ncol(norm_counts) + 4), height = max(5.5, 0.25 * top_n + 2))
  }

  n_per_condition <- paste(paste0(names(cond_counts), "=", as.integer(cond_counts)), collapse = "; ")
  min_replicates <- min(cond_counts)

  # Top DE genes by significance for the report text (res is already ordered by padj above).
  top_hits <- head(res[!is.na(res$padj), c("gene", "padj", "log2FoldChange")], 5)
  if (nrow(top_hits) > 0) {
    gene_desc <- paste0(top_hits$gene, " (padj=", signif(top_hits$padj, 3), ", log2FC=", round(top_hits$log2FoldChange, 2), ")")
    n_sig <- sum(res$change != "ns", na.rm = TRUE)
    sig_note <- if (n_sig == 0) paste0(" [none reached padj<", padj_threshold, "]") else ""
    top_gene_lines <- c(top_gene_lines, paste0("  - ", ct, ": ", paste(gene_desc, collapse = "; "), sig_note))
  } else {
    top_gene_lines <- c(top_gene_lines, paste0("  - ", ct, ": no gene had a non-missing adjusted p-value"))
  }

  # Pseudoreplication / small-n caveat (references/pitfalls.md #2): pseudobulk still needs enough
  # biological replicates per condition for the error term to be meaningful.
  if (min_replicates < 3) {
    low_replicate_lines <- c(low_replicate_lines, paste0("  - ", ct, ": ", n_per_condition, " sample(s) per condition"))
  }

  summary_rows[[ct]] <- data.frame(
    celltype = ct, genes = nrow(res), up = sum(res$change == "up", na.rm = TRUE), down = sum(res$change == "down", na.rm = TRUE),
    n_per_condition = n_per_condition, min_replicates = min_replicates
  )
}
summary <- if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame()
write.table(summary, file.path(out_dir, "results", "pseudobulk_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

total_up <- if (nrow(summary)) sum(summary$up) else 0
total_down <- if (nrow(summary)) sum(summary$down) else 0
report_lines <- c(
  "# Pseudobulk differential expression",
  "",
  "## Results",
  "",
  paste0("- Cell types tested: ", nrow(summary)),
  paste0("- Significance thresholds used: padj < ", padj_threshold, ", |log2FC| >= ", log2fc_threshold),
  paste0("- Total DE genes (padj<", padj_threshold, ", |log2FC|>=", log2fc_threshold, "): ", total_up + total_down, " (", total_up, " up, ", total_down, " down)"),
  if (length(skipped_celltypes)) paste0("- Skipped (only one condition level present): ", paste(skipped_celltypes, collapse = "; ")) else NULL
)
if (length(top_gene_lines)) {
  report_lines <- c(report_lines, "", "### Top DE genes by significance per cell type", "", top_gene_lines)
}
caveats <- character(0)
if (length(low_replicate_lines)) {
  caveats <- c(
    caveats,
    paste0(
      "The following cell type(s) were tested with fewer than 3 biological replicates in at least one condition; ",
      "treat these results as exploratory, not confirmatory, since pseudobulk still needs a real per-sample error ",
      "term (references/pitfalls.md #2):"
    )
  )
  report_lines <- c(report_lines, "", "### Caveats", "", paste0("- ", caveats), low_replicate_lines)
}
conclusion <- paste0(
  "Conclusion: pseudobulk DESeq2 testing was performed per sample x cell type to avoid cell-level pseudoreplication. ",
  if (nrow(summary)) paste0(total_up + total_down, " gene(s) were called differentially expressed across ", nrow(summary), " cell type(s) tested. ") else "No cell type had at least 2 condition levels with data to test. ",
  if (length(low_replicate_lines)) "Some cell type(s) had low replicate depth (see caveats above) and should not be treated as confirmatory. " else "",
  "Review the per-cell-type DEG tables and volcano/heatmap figures for full detail."
)
report_lines <- c(report_lines, "", "### Conclusion", "", conclusion)
writeLines(report_lines, file.path(out_dir, "report_pseudobulk.md"))
