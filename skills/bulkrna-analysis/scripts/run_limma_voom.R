#!/usr/bin/env Rscript
# Formal route for RAW COUNTS with multi-factor/batch/paired/continuous-covariate designs
# (edgeR::filterByExpr + calcNormFactors + limma::voom + lmFit + eBayes). For a simple
# two-group raw-count contrast, prefer run_deseq2.R. For already-normalized/log-scale
# expression (not raw counts), use run_limma_trend.R instead -- voom's precision weights need
# real per-sample read depth, which only raw counts carry.
args_raw <- commandArgs(trailingOnly = TRUE)

# Flag-style optional args (--name=value) can appear anywhere on the command
# line; everything else is a positional arg, so existing invocations with
# only the 5 positional args keep working unchanged.
flag_args <- grep("^--", args_raw, value = TRUE)
args <- args_raw[!grepl("^--", args_raw)]
if (length(args) < 5) {
  stop("Usage: Rscript run_limma_voom.R count_matrix.csv metadata.csv design_formula contrast output_dir [--padj-threshold=0.05] [--log2fc-threshold=1] [--top-n-heatmap=50]")
}
counts_file <- args[[1]]
metadata_file <- args[[2]]
design_formula <- args[[3]]
contrast_expr <- args[[4]]
out_dir <- args[[5]]

parse_flag <- function(flags, name, default) {
  pattern <- paste0("^--", name, "=")
  hit <- grep(pattern, flags, value = TRUE)
  if (length(hit) == 0) return(default)
  suppressWarnings(as.numeric(sub(pattern, "", hit[[length(hit)]])))
}
padj_threshold <- parse_flag(flag_args, "padj-threshold", 0.05)
log2fc_threshold <- parse_flag(flag_args, "log2fc-threshold", 1.0)
top_n_heatmap <- as.integer(parse_flag(flag_args, "top-n-heatmap", 50))
if (is.na(padj_threshold) || is.na(log2fc_threshold) || is.na(top_n_heatmap)) {
  stop("--padj-threshold, --log2fc-threshold, and --top-n-heatmap must be numeric, e.g. --padj-threshold=0.01")
}

required <- c("edgeR", "limma", "ggplot2", "pheatmap")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "qc"), recursive = TRUE, showWarnings = FALSE)

save_gg <- function(plot, name, width = 8, height = 6) {
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(name, ".png")), plot, width = width, height = height, dpi = 320, bg = "white")
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(name, ".pdf")), plot, width = width, height = height, bg = "white")
}

if (!file.exists(counts_file)) stop(sprintf("count matrix '%s' does not exist. Check the path.", counts_file))
if (!file.exists(metadata_file)) stop(sprintf("metadata file '%s' does not exist. Check the path.", metadata_file))

counts <- read.csv(counts_file, check.names = FALSE)
metadata <- read.csv(metadata_file, check.names = FALSE)
if (!("gene" %in% colnames(counts))) {
  stop(sprintf("count matrix '%s' is missing the required 'gene' column. Found columns: %s.", counts_file, paste(colnames(counts), collapse = ", ")))
}
if (!("sample" %in% colnames(metadata))) {
  stop(sprintf("metadata file '%s' is missing the required 'sample' column. Found columns: %s.", metadata_file, paste(colnames(metadata), collapse = ", ")))
}
rownames(counts) <- counts$gene
counts$gene <- NULL
matrix_only <- setdiff(colnames(counts), metadata$sample)
metadata_only <- setdiff(metadata$sample, colnames(counts))
if (length(matrix_only) || length(metadata_only)) {
  detail <- character(0)
  if (length(matrix_only)) detail <- c(detail, sprintf("in count matrix but not metadata: %s", paste(matrix_only, collapse = ", ")))
  if (length(metadata_only)) detail <- c(detail, sprintf("in metadata but not count matrix: %s", paste(metadata_only, collapse = ", ")))
  stop(sprintf("Sample ID mismatch between counts and metadata (%s).", paste(detail, collapse = "; ")))
}
metadata <- metadata[match(colnames(counts), metadata$sample), ]
if (any(is.na(metadata$sample))) stop("Sample IDs in count matrix do not match metadata.")

design_terms <- all.vars(as.formula(design_formula))
extra_covariate_cols <- setdiff(design_terms, "sample")
group_col <- design_terms[length(design_terms)]  # last term in the formula is treated as the primary group of interest for QC/caveat purposes

y <- edgeR::DGEList(counts = as.matrix(counts))
design <- model.matrix(as.formula(design_formula), metadata)
keep <- edgeR::filterByExpr(y, design = design)
if (!any(keep)) stop("No genes pass edgeR::filterByExpr(). Check that the count matrix and design formula are correct.")
genes_input <- nrow(y)
y <- y[keep, , keep.lib.sizes = FALSE]
genes_filtered_out <- genes_input - nrow(y)
y <- edgeR::calcNormFactors(y)
v <- limma::voom(y, design, plot = FALSE)
fit <- limma::lmFit(v, design)
contrast <- limma::makeContrasts(contrasts = contrast_expr, levels = design)
fit2 <- limma::contrasts.fit(fit, contrast)
fit2 <- limma::eBayes(fit2)
res <- limma::topTable(fit2, number = Inf, sort.by = "P")
res$gene <- rownames(res)
res$change <- ifelse(res$adj.P.Val < padj_threshold & res$logFC >= log2fc_threshold, "up",
  ifelse(res$adj.P.Val < padj_threshold & res$logFC <= -log2fc_threshold, "down", "ns"))
res <- res[, c("gene", setdiff(colnames(res), "gene"))]
write.table(res, file.path(out_dir, "results", "limma_voom_differential_expression.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(sample = metadata$sample, lib_size = colSums(counts)), file.path(out_dir, "qc", "library_size.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(design, file.path(out_dir, "qc", "design_matrix.tsv"), sep = "\t", quote = FALSE)

sample_annotation <- if (group_col %in% colnames(metadata)) data.frame(group = metadata[[group_col]], row.names = metadata$sample) else NULL

pca <- prcomp(t(v$E), scale. = TRUE)
percent_var <- round(100 * (pca$sdev^2) / sum(pca$sdev^2))
pca_dat <- data.frame(sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
if (!is.null(sample_annotation)) pca_dat$group <- sample_annotation[pca_dat$sample, "group"]
p <- ggplot2::ggplot(pca_dat, ggplot2::aes(PC1, PC2, label = sample, color = if (!is.null(sample_annotation)) group else NULL)) +
  ggplot2::geom_point(size = 3.2) +
  ggplot2::geom_text(vjust = -0.8, size = 3.4, show.legend = FALSE, check_overlap = TRUE) +
  ggplot2::labs(title = "PCA of voom expression",
                x = paste0("PC1 (", percent_var[1], "% variance)"),
                y = paste0("PC2 (", percent_var[2], "% variance)"), color = NULL) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20, color = "#173B5F"), legend.position = "top")
save_gg(p, "pca", width = 7.8, height = 6.2)

volcano <- ggplot2::ggplot(res, ggplot2::aes(logFC, -log10(adj.P.Val), color = change)) +
  ggplot2::geom_point(alpha = 0.82, size = 2, na.rm = TRUE) +
  ggplot2::geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "#A7C6DD", linewidth = 0.45) +
  ggplot2::geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "#A7C6DD", linewidth = 0.45) +
  ggplot2::scale_color_manual(values = c(up = "#E45B5B", down = "#4EA3F1", ns = "#B8C2CC")) +
  ggplot2::labs(
    title = "limma-voom volcano",
    subtitle = sprintf("padj < %g, |log2FC| >= %g", padj_threshold, log2fc_threshold),
    x = "log2 fold change", y = "-log10 adjusted p value", color = NULL
  ) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20, color = "#173B5F"),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 13, color = "#173B5F")
  )
save_gg(volcano, "limma_voom_volcano", width = 8.4, height = 6.4)

top <- head(na.omit(res$gene[order(res$adj.P.Val)]), top_n_heatmap)
pheatmap::pheatmap(v$E[top, , drop = FALSE], scale = "row",
                   annotation_col = sample_annotation,
                   color = colorRampPalette(c("#4EA3F1", "white", "#E45B5B"))(101),
                   filename = file.path(out_dir, "figures", "limma_voom_heatmap.png"),
                   width = 8.4, height = 8.0)
pheatmap::pheatmap(v$E[top, , drop = FALSE], scale = "row",
                   annotation_col = sample_annotation,
                   color = colorRampPalette(c("#4EA3F1", "white", "#E45B5B"))(101),
                   filename = file.path(out_dir, "figures", "limma_voom_heatmap.pdf"),
                   width = 8.4, height = 8.0)

sample_cor <- cor(v$E, method = "pearson")
if (!is.null(sample_annotation)) {
  pheatmap::pheatmap(sample_cor, annotation_col = sample_annotation, annotation_row = sample_annotation,
                     color = colorRampPalette(c("#4EA3F1", "white", "#E45B5B"))(101),
                     display_numbers = TRUE, number_color = "#173B5F", number_format = "%.2f",
                     main = "Sample correlation heatmap",
                     filename = file.path(out_dir, "figures", "sample_correlation_heatmap.png"), width = 7.6, height = 7.0)
  pheatmap::pheatmap(sample_cor, annotation_col = sample_annotation, annotation_row = sample_annotation,
                     color = colorRampPalette(c("#4EA3F1", "white", "#E45B5B"))(101),
                     display_numbers = TRUE, number_color = "#173B5F", number_format = "%.2f",
                     main = "Sample correlation heatmap",
                     filename = file.path(out_dir, "figures", "sample_correlation_heatmap.pdf"), width = 7.6, height = 7.0)
}

# ---------------------------------------------------------------------------
# QC verdict + caveats (same pattern as run_deseq2.R/run_limma_trend.R): flag
# outlier samples via same-group correlation gap, and any extra covariate
# fully confounded with the primary group term.
# ---------------------------------------------------------------------------
outlier_notes <- character(0)
if (group_col %in% colnames(metadata)) {
  for (g in unique(metadata[[group_col]])) {
    members <- metadata$sample[metadata[[group_col]] == g]
    if (length(members) < 2) {
      outlier_notes <- c(outlier_notes, sprintf("Group '%s' has only %d sample(s); within-group correlation consistency cannot be assessed.", g, length(members)))
      next
    }
    own_group_corr <- vapply(members, function(s) mean(sample_cor[s, setdiff(members, s)]), numeric(1))
    group_avg <- mean(own_group_corr)
    for (s in members) {
      gap <- group_avg - own_group_corr[[s]]
      if (gap > 0.1) {
        outlier_notes <- c(outlier_notes, sprintf(
          "Sample '%s' (group '%s') correlates with its own group at r=%.2f, %.2f below the group average (r=%.2f) -- possible outlier, mislabel, or batch effect.",
          s, g, own_group_corr[[s]], gap, group_avg
        ))
      }
    }
  }
}

n_samples <- nrow(metadata)
min_group_n <- if (group_col %in% colnames(metadata)) min(table(metadata[[group_col]])) else n_samples
caveats <- character(0)
caveat_titles <- character(0)
if (min_group_n < 3) {
  caveat_titles <- c(caveat_titles, "underpowered design")
  caveats <- c(caveats, sprintf(
    "**Underpowered design** (pitfalls.md #10): the smallest group has only %d sample(s), below the recommended minimum of 3 biological replicates per group. Treat DEG calls as exploratory.",
    min_group_n
  ))
}
if (min_group_n > 30) {
  caveat_titles <- c(caveat_titles, "effect-size floor at large n")
  caveats <- c(caveats, sprintf(
    "**Effect-size floor at large n** (pitfalls.md #12): with a smallest group size of %d, statistical significance can be reached for biologically trivial fold-changes; the abs(log2FC) >= %g filter is kept as a hard floor.",
    min_group_n, log2fc_threshold
  ))
}
if (length(extra_covariate_cols) > 1 && group_col %in% colnames(metadata)) {
  for (col in setdiff(extra_covariate_cols, group_col)) {
    if (!(col %in% colnames(metadata))) next
    tab <- table(metadata[[col]], metadata[[group_col]])
    if (nrow(tab) < 2) next
    levels_per_group <- apply(tab > 0, 1, sum)
    if (all(levels_per_group == 1)) {
      caveat_titles <- c(caveat_titles, sprintf("'%s' confounded with '%s'", col, group_col))
      caveats <- c(caveats, sprintf(
        "**Confounded design** (pitfalls.md #13): metadata column '%s' is perfectly confounded with '%s' -- the design matrix cannot separate these two effects.",
        col, group_col
      ))
    }
  }
}

up_count <- sum(res$change == "up", na.rm = TRUE)
down_count <- sum(res$change == "down", na.rm = TRUE)
# "Top genes" means significant (up/down) genes, not just small-padj genes -- see the same fix
# applied in run_deseq2.R/run_limma_trend.R/quickplot_bulk.py.
sig_res <- res[!is.na(res$adj.P.Val) & res$change != "ns", ]
top_n <- min(5, nrow(sig_res))
top_genes <- head(sig_res[order(sig_res$adj.P.Val), ], top_n)

report_lines <- c(
  "# limma-voom report", "",
  sprintf("Formal limma-voom analysis. Design: `%s`. Contrast: `%s`.", design_formula, contrast_expr),
  "",
  "## Input summary", "",
  sprintf("- Samples: %d", n_samples),
  sprintf("- Genes in input matrix: %d; genes tested after edgeR::filterByExpr: %d", genes_input, nrow(res)),
  sprintf("- Thresholds: padj < %g and abs(log2FC) >= %g", padj_threshold, log2fc_threshold),
  "",
  "## Key findings", "",
  sprintf("- Significant up: %d gene(s)", up_count),
  sprintf("- Significant down: %d gene(s)", down_count)
)
if (top_n > 0) {
  report_lines <- c(report_lines, sprintf("- Top %d gene(s) by adjusted p value:", top_n))
  for (i in seq_len(top_n)) {
    report_lines <- c(report_lines, sprintf("  %d. `%s`: log2FC = %.2f, padj = %.3g (%s)", i, top_genes$gene[i], top_genes$logFC[i], top_genes$adj.P.Val[i], top_genes$change[i]))
  }
}

report_lines <- c(report_lines, "", "## QC verdict", "")
if (length(outlier_notes)) {
  report_lines <- c(report_lines, paste0("- ", outlier_notes))
} else if (group_col %in% colnames(metadata)) {
  report_lines <- c(report_lines, "- No sample fell outside the same-group correlation heuristic; no outlier flagged.")
} else {
  report_lines <- c(report_lines, sprintf("- Design formula's last term ('%s') is not a metadata column; per-group QC skipped. Sample correlation heatmap still reflects overall sample similarity.", group_col))
}

report_lines <- c(report_lines, "", "## Caveats (automatically checked against references/pitfalls.md)", "")
if (length(caveats)) {
  report_lines <- c(report_lines, paste0("- ", caveats))
} else {
  report_lines <- c(report_lines, "- No data-triggered caveats from the automatic checks. This does not replace manual review of `references/pitfalls.md`.")
}

qc_sentence <- if (length(outlier_notes) == 0) "sample correlation shows no evidence of outlier samples" else "sample correlation flagged one or more issues noted above"
caveat_clause <- if (length(caveat_titles)) paste0(", with the following flagged concerns: ", paste(caveat_titles, collapse = ", ")) else ""
report_lines <- c(report_lines, "", "## Conclusion", "",
  sprintf(
    "%d gene(s) were significantly upregulated and %d downregulated for contrast `%s` (design `%s`) at padj<%g and abs(log2FC)>=%g. Voom-weighted %s. Results are based on a formal limma-voom empirical-Bayes moderated t-test with read-depth-aware precision weighting%s.",
    up_count, down_count, contrast_expr, design_formula, padj_threshold, log2fc_threshold, qc_sentence, caveat_clause
  )
)

report_lines <- c(report_lines, "", "## Files", "",
  "Figures: `limma_voom_volcano`, `limma_voom_heatmap`, `pca`, `sample_correlation_heatmap` (when a primary group column is available), each exported as PNG and PDF. Full gene table: `results/limma_voom_differential_expression.tsv`. QC: `qc/library_size.tsv`, `qc/design_matrix.tsv`."
)
writeLines(report_lines, file.path(out_dir, "report_limma_voom.md"))
write.table(data.frame(item = c("samples", "genes_input", "genes_filtered_out", "genes_tested", "up_genes", "down_genes", "padj_threshold", "log2fc_threshold", "top_n_heatmap"),
                       value = c(n_samples, genes_input, genes_filtered_out, nrow(res), up_count, down_count, padj_threshold, log2fc_threshold, top_n_heatmap)),
            file.path(out_dir, "summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
