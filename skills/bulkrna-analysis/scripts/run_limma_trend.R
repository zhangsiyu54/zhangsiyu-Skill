#!/usr/bin/env Rscript
# Formal route for ALREADY-NORMALIZED / LOG-SCALE bulk RNA-seq expression input (TPM/FPKM/CPM/
# logCPM/microarray-style/log-ratio matrices exported from a public dataset, a collaborator's
# pipeline, or an Excel workbook) -- NOT raw integer counts. Neither `run_deseq2.R` nor
# `run_limma_voom.R` is correct for this input (both assume raw counts and re-normalize
# internally); this script fills that gap with limma's `eBayes(trend=TRUE)` empirical Bayes
# model, which is the standard limma route for already-normalized continuous expression
# (Ritchie et al. 2015 Nucleic Acids Res; Smyth 2004 Stat Appl Genet Mol Biol -- the same
# empirical-Bayes moderated t-test used by run_deseq2.R/run_limma_voom.R, applied directly to
# provided expression values instead of raw counts run through voom's read-depth weighting).
#
# Usage:
#   Rscript run_limma_trend.R expression_matrix.csv metadata.csv group_col case control [output_dir] \
#       [--padj-threshold=0.05] [--log2fc-threshold=1] [--top-n-heatmap=50]
#
# expression_matrix.csv: first column `gene`, remaining columns sample IDs, already log-scale
#   continuous expression (may be negative, e.g. z-scored or already-centered values).
# metadata.csv: `sample` column plus a group column and any extra covariates to check for
#   confounding.

args_raw <- commandArgs(trailingOnly = TRUE)
flag_args <- grep("^--", args_raw, value = TRUE)
args <- args_raw[!grepl("^--", args_raw)]
if (length(args) < 5) {
  stop("Usage: Rscript run_limma_trend.R expression_matrix.csv metadata.csv group_col case control [output_dir] [--padj-threshold=0.05] [--log2fc-threshold=1] [--top-n-heatmap=50]")
}
expr_file <- args[[1]]
metadata_file <- args[[2]]
group_col <- args[[3]]
case_label <- args[[4]]
control_label <- args[[5]]
out_dir <- ifelse(length(args) >= 6, args[[6]], "output")

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

required <- c("limma", "ggplot2", "pheatmap")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "qc"), recursive = TRUE, showWarnings = FALSE)

save_gg <- function(plot, name, width = 8, height = 6) {
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(name, ".png")), plot, width = width, height = height, dpi = 320, bg = "white")
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(name, ".pdf")), plot, width = width, height = height, bg = "white")
}

# ---------------------------------------------------------------------------
# Input validation: fail fast with actionable messages, never a raw stack
# trace / cryptic subscript error (references/pitfalls.md).
# ---------------------------------------------------------------------------
if (!file.exists(expr_file)) stop(sprintf("expression matrix '%s' does not exist. Check the path.", expr_file))
if (!file.exists(metadata_file)) stop(sprintf("metadata file '%s' does not exist. Check the path.", metadata_file))

read_table_auto <- function(path) {
  first_line <- readLines(path, n = 1, warn = FALSE)
  sep <- if (grepl("\t", first_line)) "\t" else ","
  read.delim(path, sep = sep, check.names = FALSE, stringsAsFactors = FALSE)
}
expr_tbl <- read_table_auto(expr_file)
metadata <- read_table_auto(metadata_file)

if (!("gene" %in% colnames(expr_tbl))) {
  stop(sprintf(
    "expression matrix '%s' is missing the required 'gene' column (first column). Found columns: %s.",
    expr_file, paste(colnames(expr_tbl), collapse = ", ")
  ))
}
if (!("sample" %in% colnames(metadata))) {
  stop(sprintf(
    "metadata file '%s' is missing the required 'sample' column. Found columns: %s.",
    metadata_file, paste(colnames(metadata), collapse = ", ")
  ))
}
if (!(group_col %in% colnames(metadata))) {
  other_cols <- setdiff(colnames(metadata), "sample")
  stop(sprintf(
    "metadata file '%s' has no column named '%s'. Available non-sample columns: %s.",
    metadata_file, group_col, paste(other_cols, collapse = ", ")
  ))
}
if (nrow(expr_tbl) == 0) stop(sprintf("expression matrix '%s' has a header but zero gene rows.", expr_file))
if (nrow(metadata) == 0) stop(sprintf("metadata file '%s' has a header but zero sample rows.", metadata_file))

dup_genes <- unique(expr_tbl$gene[duplicated(expr_tbl$gene)])
if (length(dup_genes)) {
  stop(sprintf(
    "expression matrix '%s' has duplicate gene IDs: %s. Aggregate duplicates before re-running.",
    expr_file, paste(head(dup_genes, 10), collapse = ", ")
  ))
}
dup_samples <- unique(metadata$sample[duplicated(metadata$sample)])
if (length(dup_samples)) {
  stop(sprintf("metadata file '%s' has duplicate 'sample' values: %s.", metadata_file, paste(dup_samples, collapse = ", ")))
}

rownames(expr_tbl) <- expr_tbl$gene
expr_tbl$gene <- NULL
expr_all <- as.matrix(expr_tbl)
storage.mode(expr_all) <- "double"

matrix_only <- setdiff(colnames(expr_all), metadata$sample)
metadata_only <- setdiff(metadata$sample, colnames(expr_all))
if (length(matrix_only) || length(metadata_only)) {
  detail <- character(0)
  if (length(matrix_only)) detail <- c(detail, sprintf("in expression matrix but not metadata: %s", paste(matrix_only, collapse = ", ")))
  if (length(metadata_only)) detail <- c(detail, sprintf("in metadata but not expression matrix: %s", paste(metadata_only, collapse = ", ")))
  stop(sprintf("Sample ID mismatch between expression matrix and metadata (%s). Sample IDs must match exactly.", paste(detail, collapse = "; ")))
}
metadata <- metadata[match(colnames(expr_all), metadata$sample), ]
if (any(is.na(metadata$sample))) stop("Sample IDs in expression matrix do not match metadata.")

if (!(case_label %in% metadata[[group_col]])) {
  stop(sprintf("No samples have %s == '%s' (case label). Observed values: %s.", group_col, case_label, paste(sort(unique(metadata[[group_col]])), collapse = ", ")))
}
if (!(control_label %in% metadata[[group_col]])) {
  stop(sprintf("No samples have %s == '%s' (control label). Observed values: %s.", group_col, control_label, paste(sort(unique(metadata[[group_col]])), collapse = ", ")))
}
n_case <- sum(metadata[[group_col]] == case_label)
n_control <- sum(metadata[[group_col]] == control_label)
metadata[[group_col]] <- factor(metadata[[group_col]], levels = c(control_label, case_label))
extra_covariate_cols <- setdiff(colnames(metadata), c("sample", group_col))

# Genes with any non-finite (NA/Inf) value across samples cannot be modeled by lmFit; drop them
# with an explicit, reported count rather than letting lmFit silently mis-handle them.
genes_input <- nrow(expr_all)
keep <- rowSums(is.finite(expr_all)) == ncol(expr_all)
if (!any(keep)) stop("No genes have finite values across all samples. Check the expression matrix for formatting issues.")
expr <- expr_all[keep, , drop = FALSE]
genes_filtered_out <- genes_input - nrow(expr)

# A near-all-non-negative-integer matrix is a strong signal the caller actually has raw counts
# and should be using run_deseq2.R/run_limma_voom.R instead of this script (pitfalls.md-style
# scale-detection check, mirrored from quickplot_bulk.py's expression-scale handling).
looks_like_raw_counts <- all(expr_all >= 0, na.rm = TRUE) && mean(abs(expr_all - round(expr_all)) < 1e-6, na.rm = TRUE) > 0.98

# ---------------------------------------------------------------------------
# limma trend model (no voom: there is no raw read-depth information to build
# voom's precision weights from, since the input is already normalized/log-scale).
# ---------------------------------------------------------------------------
design <- model.matrix(~0 + metadata[[group_col]])
colnames(design) <- levels(metadata[[group_col]])
contrast <- limma::makeContrasts(contrasts = paste0(case_label, " - ", control_label), levels = design)
fit <- limma::lmFit(expr, design)
fit2 <- limma::contrasts.fit(fit, contrast)
fit2 <- limma::eBayes(fit2, trend = TRUE)

res <- limma::topTable(fit2, number = Inf, sort.by = "P")
res$gene <- rownames(res)
res <- res[, c("gene", setdiff(colnames(res), "gene"))]
colnames(res)[colnames(res) == "logFC"] <- "log2FC"
colnames(res)[colnames(res) == "AveExpr"] <- "average_expression"
colnames(res)[colnames(res) == "P.Value"] <- "pvalue"
colnames(res)[colnames(res) == "adj.P.Val"] <- "padj"
res$change <- ifelse(res$padj < padj_threshold & res$log2FC >= log2fc_threshold, "up",
  ifelse(res$padj < padj_threshold & res$log2FC <= -log2fc_threshold, "down", "ns"))
res <- res[order(res$padj, -abs(res$log2FC)), ]
write.table(res, file.path(out_dir, "results", "differential_expression.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(gene = rownames(expr), expr, check.names = FALSE),
            file.path(out_dir, "results", "normalized_expression_used.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(res$gene[res$change == "up"], file.path(out_dir, "results", "up_genes.txt"))
writeLines(res$gene[res$change == "down"], file.path(out_dir, "results", "down_genes.txt"))
write.table(data.frame(gene = res$gene, rank_score = sign(res$log2FC) * -log10(pmax(res$pvalue, 1e-300))),
            file.path(out_dir, "results", "gsea_ranked_list.rnk"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

sample_annotation <- data.frame(group = metadata[[group_col]], row.names = metadata$sample)

volcano_plot <- ggplot2::ggplot(res, ggplot2::aes(log2FC, -log10(padj), color = change)) +
  ggplot2::geom_point(alpha = 0.82, size = 2.1, na.rm = TRUE) +
  ggplot2::geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "#A7C6DD", linewidth = 0.45) +
  ggplot2::geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "#A7C6DD", linewidth = 0.45) +
  ggplot2::scale_color_manual(values = c(up = "#E45B5B", down = "#4EA3F1", ns = "#B8C2CC")) +
  ggplot2::labs(
    title = "Volcano plot",
    subtitle = sprintf("limma trend, padj < %g, |log2FC| >= %g", padj_threshold, log2fc_threshold),
    x = "log2 fold change", y = "-log10 adjusted p value", color = NULL
  ) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20, color = "#173B5F"),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 13, color = "#173B5F"),
    axis.title = ggplot2::element_text(size = 17, color = "#173B5F"),
    axis.text = ggplot2::element_text(size = 13, color = "#173B5F"),
    legend.position = "top"
  )
save_gg(volcano_plot, "volcano", width = 8.8, height = 6.6)

top <- head(na.omit(res$gene[order(res$padj)]), top_n_heatmap)
pheatmap::pheatmap(expr[top, , drop = FALSE], scale = "row",
                   annotation_col = sample_annotation,
                   color = colorRampPalette(c("#4EA3F1", "white", "#E45B5B"))(101),
                   fontsize = 12, fontsize_row = 8, fontsize_col = 11,
                   filename = file.path(out_dir, "figures", "heatmap_top_degs.png"), width = 8.4, height = 8.0)
pheatmap::pheatmap(expr[top, , drop = FALSE], scale = "row",
                   annotation_col = sample_annotation,
                   color = colorRampPalette(c("#4EA3F1", "white", "#E45B5B"))(101),
                   fontsize = 12, fontsize_row = 8, fontsize_col = 11,
                   filename = file.path(out_dir, "figures", "heatmap_top_degs.pdf"), width = 8.4, height = 8.0)

# PCA directly on the provided expression values -- no vst/rlog needed since the input is
# already normalized/log-scale.
pca <- prcomp(t(expr), scale. = FALSE)
percent_var <- round(100 * (pca$sdev^2) / sum(pca$sdev^2))
pca_data <- data.frame(sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = metadata[[group_col]])
pca_plot <- ggplot2::ggplot(pca_data, ggplot2::aes(PC1, PC2, color = group, label = sample)) +
  ggplot2::geom_point(size = 3.4, alpha = 0.9) +
  ggplot2::geom_text(vjust = -0.9, size = 3.4, show.legend = FALSE, check_overlap = TRUE) +
  ggplot2::scale_color_manual(values = stats::setNames(c("#4EA3F1", "#E45B5B"), c(control_label, case_label))) +
  ggplot2::labs(
    title = "PCA of samples",
    x = paste0("PC1 (", percent_var[1], "% variance)"),
    y = paste0("PC2 (", percent_var[2], "% variance)"),
    color = NULL
  ) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20, color = "#173B5F"),
    axis.title = ggplot2::element_text(size = 17, color = "#173B5F"),
    axis.text = ggplot2::element_text(size = 13, color = "#173B5F"),
    legend.position = "top"
  )
save_gg(pca_plot, "pca", width = 7.8, height = 6.4)

sample_cor <- cor(expr, method = "pearson")
pheatmap::pheatmap(sample_cor,
                   annotation_col = sample_annotation, annotation_row = sample_annotation,
                   color = colorRampPalette(c("#4EA3F1", "white", "#E45B5B"))(101),
                   display_numbers = TRUE, number_color = "#173B5F", number_format = "%.2f",
                   fontsize = 12, fontsize_row = 11, fontsize_col = 11,
                   main = "Sample correlation heatmap",
                   filename = file.path(out_dir, "figures", "sample_correlation_heatmap.png"), width = 7.6, height = 7.0)
pheatmap::pheatmap(sample_cor,
                   annotation_col = sample_annotation, annotation_row = sample_annotation,
                   color = colorRampPalette(c("#4EA3F1", "white", "#E45B5B"))(101),
                   display_numbers = TRUE, number_color = "#173B5F", number_format = "%.2f",
                   fontsize = 12, fontsize_row = 11, fontsize_col = 11,
                   main = "Sample correlation heatmap",
                   filename = file.path(out_dir, "figures", "sample_correlation_heatmap.pdf"), width = 7.6, height = 7.0)

# No raw read count exists for this route, so "library size"/"detected genes" are approximated
# from the provided expression values (total signal, and features above the matrix minimum) --
# a QC proxy, not a literal sequencing-depth QC metric.
total_signal <- colSums(expr)
detected_features <- colSums(expr > min(expr, na.rm = TRUE))
qc_long <- data.frame(
  sample = factor(rep(metadata$sample, 2), levels = metadata$sample),
  group = rep(metadata[[group_col]], 2),
  metric = rep(c("Total signal (proxy)", "Detected features (proxy)"), each = nrow(metadata)),
  value = c(total_signal[metadata$sample], detected_features[metadata$sample])
)
qc_plot <- ggplot2::ggplot(qc_long, ggplot2::aes(sample, value, fill = group)) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::facet_wrap(~metric, scales = "free_y") +
  ggplot2::scale_fill_manual(values = stats::setNames(c("#4EA3F1", "#E45B5B"), c(control_label, case_label))) +
  ggplot2::labs(title = "Sample QC summary", x = NULL, y = NULL, fill = NULL) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20, color = "#173B5F"),
    axis.text.x = ggplot2::element_text(angle = 38, hjust = 1, size = 11, color = "#173B5F"),
    axis.text.y = ggplot2::element_text(size = 13, color = "#173B5F"),
    strip.text = ggplot2::element_text(size = 15, face = "bold", color = "#173B5F"),
    legend.position = "top"
  )
save_gg(qc_plot, "qc_summary", width = 10.5, height = 6.0)

write.table(data.frame(sample = metadata$sample, group = metadata[[group_col]],
                       total_signal = as.numeric(total_signal[metadata$sample]),
                       detected_features = as.numeric(detected_features[metadata$sample])),
            file.path(out_dir, "qc", "sample_qc.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(design, file.path(out_dir, "qc", "design_matrix.tsv"), sep = "\t", quote = FALSE)

# ---------------------------------------------------------------------------
# QC verdict (same pattern as run_deseq2.R): flag outlier samples via
# same-group correlation gap, and any extra covariate fully confounded with
# the group of interest.
# ---------------------------------------------------------------------------
outlier_notes <- character(0)
for (g in levels(metadata[[group_col]])) {
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

confounded_covariates <- character(0)
if (length(extra_covariate_cols)) {
  for (col in extra_covariate_cols) {
    tab <- table(metadata[[col]], metadata[[group_col]])
    if (nrow(tab) < 2) next
    levels_per_group <- apply(tab > 0, 1, sum)
    if (all(levels_per_group == 1)) confounded_covariates <- c(confounded_covariates, col)
  }
}

# ---------------------------------------------------------------------------
# Caveats: only the ones this run's data actually triggers (references/pitfalls.md).
# ---------------------------------------------------------------------------
caveats <- character(0)
caveat_titles <- character(0)
if (n_case < 3 || n_control < 3) {
  caveat_titles <- c(caveat_titles, "underpowered design")
  caveats <- c(caveats, sprintf(
    "**Underpowered design** (pitfalls.md #10): n=%d (%s) vs n=%d (%s) is below the recommended minimum of 3 biological replicates per group. Variance estimates are unstable at this replication level even with empirical-Bayes moderation; treat results as exploratory.",
    n_case, case_label, n_control, control_label
  ))
}
if (n_case > 30 && n_control > 30) {
  caveat_titles <- c(caveat_titles, "effect-size floor at large n")
  caveats <- c(caveats, sprintf(
    "**Effect-size floor at large n** (pitfalls.md #12): with n=%d/%d per group, statistical significance can be reached for biologically trivial fold-changes; the abs(log2FC) >= %g filter is kept as a hard floor.",
    n_case, n_control, log2fc_threshold
  ))
}
for (col in confounded_covariates) {
  caveat_titles <- c(caveat_titles, sprintf("'%s' confounded with '%s'", col, group_col))
  caveats <- c(caveats, sprintf(
    "**Confounded design** (pitfalls.md #13): metadata column '%s' is perfectly confounded with '%s' -- the design cannot separate these two effects.",
    col, group_col
  ))
}
if (looks_like_raw_counts) {
  caveat_titles <- c(caveat_titles, "input may actually be raw counts")
  caveats <- c(caveats, sprintf(
    "**Input may actually be raw counts** (pitfalls.md #11): %.0f%% of values in this matrix are non-negative integers, which is also consistent with raw/unnormalized counts. If this is raw count data, `run_deseq2.R` or `run_limma_voom.R` (which model count-specific mean-variance structure and normalize internally) are more appropriate than this already-normalized-input route; re-run there if in doubt.",
    100 * mean(expr_all >= 0 & abs(expr_all - round(expr_all)) < 1e-6, na.rm = TRUE)
  ))
}
caveat_titles <- c(caveat_titles, "trend-based, not voom-weighted")
caveats <- c(caveats, "**No read-depth-based precision weighting** (Ritchie et al. 2015): this route uses `eBayes(trend=TRUE)`, which moderates variance using the mean-expression trend across genes, not `voom`'s per-observation weights derived from sequencing depth -- because no raw read-depth information exists for already-normalized input. This is the standard and correct limma approach for this input type, but is a different (not worse, just different-inputs) model from `run_limma_voom.R`.")

up_count <- sum(res$change == "up", na.rm = TRUE)
down_count <- sum(res$change == "down", na.rm = TRUE)
# "Top genes" means significant (up/down) genes, not just small-padj genes -- a gene can have a
# tiny padj but sit below the log2FC floor (change == "ns") and would otherwise misleadingly
# rank #1 in a "top gene" list despite not clearing the effect-size threshold.
sig_res <- res[!is.na(res$padj) & res$change != "ns", ]
top_n <- min(5, nrow(sig_res))
top_genes <- head(sig_res, top_n)

report_lines <- c(
  "# bulkRNA limma (trend) report",
  "",
  sprintf(
    "Formal limma analysis for **%s** vs **%s** (%s = reference/control) on already-normalized/log-scale expression input, contrast `%s - %s`.",
    case_label, control_label, control_label, case_label, control_label
  ),
  "",
  "## Input summary",
  "",
  sprintf("- Samples: %d %s vs %d %s", n_case, case_label, n_control, control_label),
  sprintf("- Genes in input matrix: %d; genes tested after removing %d gene(s) with non-finite values: %d", genes_input, genes_filtered_out, nrow(res)),
  sprintf("- Contrast: %s vs %s; thresholds: padj < %g and abs(log2FC) >= %g", case_label, control_label, padj_threshold, log2fc_threshold),
  "- Model: limma `lmFit` + `contrasts.fit` + `eBayes(trend=TRUE)` on provided expression values (no voom, no DESeq2/edgeR re-normalization)."
)
if (length(extra_covariate_cols)) {
  report_lines <- c(report_lines, sprintf("- Additional metadata columns present (checked for confounding): %s", paste(extra_covariate_cols, collapse = ", ")))
}
report_lines <- c(report_lines, "", "## Key findings", "",
  sprintf("- Significant up in %s: %d gene(s)", case_label, up_count),
  sprintf("- Significant down in %s: %d gene(s)", case_label, down_count)
)
if (top_n > 0) {
  report_lines <- c(report_lines, sprintf("- Top %d gene(s) by adjusted p value:", top_n))
  for (i in seq_len(top_n)) {
    report_lines <- c(report_lines, sprintf("  %d. `%s`: log2FC = %.2f, padj = %.3g (%s)", i, top_genes$gene[i], top_genes$log2FC[i], top_genes$padj[i], top_genes$change[i]))
  }
}

report_lines <- c(report_lines, "", "## QC verdict", "")
if (length(outlier_notes)) {
  report_lines <- c(report_lines, paste0("- ", outlier_notes))
} else {
  report_lines <- c(report_lines, "- No sample fell outside the same-group correlation heuristic; no outlier flagged.")
}
if (length(confounded_covariates)) {
  report_lines <- c(report_lines, sprintf("- Covariate '%s' is fully confounded with the group of interest (see caveats below).", confounded_covariates))
} else if (length(extra_covariate_cols)) {
  report_lines <- c(report_lines, "- No additional metadata covariate is fully confounded with the group of interest.")
}

report_lines <- c(report_lines, "", "## Caveats (automatically checked against references/pitfalls.md)", "", paste0("- ", caveats))

qc_sentence <- if (length(outlier_notes) == 0 && length(confounded_covariates) == 0) {
  "PCA and sample correlation show no evidence of outlier samples or covariate confounding"
} else {
  "PCA/sample correlation flagged one or more issues noted above"
}
caveat_clause <- if (length(caveat_titles) > 1) paste0(", with the following flagged concerns: ", paste(caveat_titles, collapse = ", ")) else ""
report_lines <- c(report_lines, "", "## Conclusion", "",
  sprintf(
    "%d gene(s) were significantly upregulated and %d downregulated in %s vs %s (n=%d/%d) at padj<%g and abs(log2FC)>=%g. %s. Results are based on a formal limma empirical-Bayes moderated t-test on already-normalized expression%s.",
    up_count, down_count, case_label, control_label, n_case, n_control, padj_threshold, log2fc_threshold, qc_sentence, caveat_clause
  )
)

report_lines <- c(report_lines, "", "## Files", "",
  "Figures: `volcano`, `heatmap_top_degs`, `pca`, `sample_correlation_heatmap`, `qc_summary`, each exported as PNG and PDF. Full gene table: `results/differential_expression.tsv`. Enrichment-ready lists: `results/up_genes.txt`, `results/down_genes.txt`, `results/gsea_ranked_list.rnk`. Per-sample QC: `qc/sample_qc.tsv`."
)

writeLines(report_lines, file.path(out_dir, "report.md"))
write.table(data.frame(item = c("case_samples", "control_samples", "genes_input", "genes_filtered_out", "genes_tested", "up_genes", "down_genes", "padj_threshold", "log2fc_threshold", "top_n_heatmap"),
                       value = c(n_case, n_control, genes_input, genes_filtered_out, nrow(res), up_count, down_count, padj_threshold, log2fc_threshold, top_n_heatmap)),
            file.path(out_dir, "summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
