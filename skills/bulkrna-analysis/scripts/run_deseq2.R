#!/usr/bin/env Rscript
args_raw <- commandArgs(trailingOnly = TRUE)

# Flag-style optional args (--name=value) can appear anywhere on the command
# line; everything else is treated as a positional arg, so existing
# invocations with only the 5-6 positional args keep working unchanged.
flag_args <- grep("^--", args_raw, value = TRUE)
args <- args_raw[!grepl("^--", args_raw)]
if (length(args) < 5) {
  stop("Usage: Rscript run_deseq2.R count_matrix.csv metadata.csv group_col case control [output_dir] [--padj-threshold=0.05] [--log2fc-threshold=1] [--top-n-heatmap=50]")
}
counts_file <- args[[1]]
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

required <- c("DESeq2", "ggplot2", "pheatmap")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "), "\nInstall Bioconductor/CRAN dependencies before running formal analysis.")
}

dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "qc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "run"), recursive = TRUE, showWarnings = FALSE)

save_gg <- function(plot, name, width = 8, height = 6) {
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(name, ".png")), plot, width = width, height = height, dpi = 320, bg = "white")
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(name, ".pdf")), plot, width = width, height = height, bg = "white")
}

# ---------------------------------------------------------------------------
# Input validation: fail fast with actionable messages (never a raw stack
# trace / cryptic subscript error) for the malformed-input cases in
# references/pitfalls.md and references/input-contract.md.
# ---------------------------------------------------------------------------
if (!file.exists(counts_file)) stop(sprintf("count matrix '%s' does not exist. Check the path.", counts_file))
if (!file.exists(metadata_file)) stop(sprintf("metadata file '%s' does not exist. Check the path.", metadata_file))

counts <- read.csv(counts_file, check.names = FALSE)
metadata <- read.csv(metadata_file, check.names = FALSE)

if (!("gene" %in% colnames(counts))) {
  stop(sprintf(
    "count matrix '%s' is missing the required 'gene' column (first column). Found columns: %s.",
    counts_file, paste(colnames(counts), collapse = ", ")
  ))
}
if (!("sample" %in% colnames(metadata))) {
  stop(sprintf(
    "metadata file '%s' is missing the required 'sample' column. Found columns: %s. Add a 'sample' column matching the count matrix header.",
    metadata_file, paste(colnames(metadata), collapse = ", ")
  ))
}
if (!(group_col %in% colnames(metadata))) {
  other_cols <- setdiff(colnames(metadata), "sample")
  stop(sprintf(
    "metadata file '%s' has no column named '%s'. Available non-sample columns: %s. Pass the correct group column argument.",
    metadata_file, group_col, paste(other_cols, collapse = ", ")
  ))
}
if (nrow(counts) == 0) stop(sprintf("count matrix '%s' has a header but zero gene rows.", counts_file))
if (nrow(metadata) == 0) stop(sprintf("metadata file '%s' has a header but zero sample rows.", metadata_file))

dup_genes <- unique(counts$gene[duplicated(counts$gene)])
if (length(dup_genes)) {
  stop(sprintf(
    "count matrix '%s' has duplicate gene IDs: %s. Aggregate duplicates (e.g. sum counts per gene ID) before re-running so each gene has exactly one row.",
    counts_file, paste(head(dup_genes, 10), collapse = ", ")
  ))
}
dup_samples <- unique(metadata$sample[duplicated(metadata$sample)])
if (length(dup_samples)) {
  stop(sprintf(
    "metadata file '%s' has duplicate 'sample' values: %s. Each sample must appear exactly once; deduplicate or rename before re-running.",
    metadata_file, paste(dup_samples, collapse = ", ")
  ))
}

rownames(counts) <- counts$gene
counts$gene <- NULL

matrix_only <- setdiff(colnames(counts), metadata$sample)
metadata_only <- setdiff(metadata$sample, colnames(counts))
if (length(matrix_only) || length(metadata_only)) {
  detail <- character(0)
  if (length(matrix_only)) detail <- c(detail, sprintf("in count matrix but not metadata: %s", paste(matrix_only, collapse = ", ")))
  if (length(metadata_only)) detail <- c(detail, sprintf("in metadata but not count matrix: %s", paste(metadata_only, collapse = ", ")))
  stop(sprintf(
    "Sample ID mismatch between counts and metadata (%s). Sample IDs must match exactly, including case, between the count matrix header and metadata$sample.",
    paste(detail, collapse = "; ")
  ))
}

metadata <- metadata[match(colnames(counts), metadata$sample), ]
if (any(is.na(metadata$sample))) stop("Sample IDs in count matrix do not match metadata.")

count_values <- as.matrix(counts)
if (any(is.na(count_values)) || any(!is.finite(count_values))) {
  stop(sprintf("count matrix '%s' has non-numeric or missing values. Fix or remove the offending cells before re-running.", counts_file))
}
if (any(count_values < 0)) {
  stop(sprintf(
    "count matrix '%s' contains negative values. Counts must be non-negative; check upstream quantification/export before re-running.",
    counts_file
  ))
}
has_non_integer <- any(abs(count_values - round(count_values)) > 1e-6)

if (!(case_label %in% metadata[[group_col]])) {
  stop(sprintf(
    "No samples have %s == '%s' (case label). Observed %s values: %s.",
    group_col, case_label, group_col, paste(sort(unique(metadata[[group_col]])), collapse = ", ")
  ))
}
if (!(control_label %in% metadata[[group_col]])) {
  stop(sprintf(
    "No samples have %s == '%s' (control label). Observed %s values: %s.",
    group_col, control_label, group_col, paste(sort(unique(metadata[[group_col]])), collapse = ", ")
  ))
}

n_case <- sum(metadata[[group_col]] == case_label)
n_control <- sum(metadata[[group_col]] == control_label)
metadata[[group_col]] <- factor(metadata[[group_col]], levels = c(control_label, case_label))

# Extra covariate columns (beyond sample/group) are checked later for
# confounding with the group of interest (pitfalls.md #13).
extra_covariate_cols <- setdiff(colnames(metadata), c("sample", group_col))

# ---------------------------------------------------------------------------
# Low-count gene filtering (references/cns-method-guide.md: "Filter low-count
# genes before modeling"), tied to library depth and the smallest group size.
# ---------------------------------------------------------------------------
genes_input <- nrow(counts)
min_group_n <- min(n_case, n_control)
keep <- rowSums(count_values >= 10) >= max(1, min_group_n)
if (!any(keep)) {
  stop("No genes pass the minimal detection filter (count >= 10 in at least one full group's worth of samples). Check that the count matrix is correct.")
}
counts <- counts[keep, , drop = FALSE]
genes_filtered_out <- genes_input - nrow(counts)

dds <- DESeq2::DESeqDataSetFromMatrix(round(as.matrix(counts)), metadata, as.formula(paste("~", group_col)))
dds <- DESeq2::DESeq(dds)
res <- as.data.frame(DESeq2::results(dds, contrast = c(group_col, case_label, control_label)))
res$gene <- rownames(res)
res$change <- ifelse(!is.na(res$padj) & res$padj < padj_threshold & res$log2FoldChange >= log2fc_threshold, "up",
  ifelse(!is.na(res$padj) & res$padj < padj_threshold & res$log2FoldChange <= -log2fc_threshold, "down", "ns"))
res <- res[order(res$padj), c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "change")]
write.table(res, file.path(out_dir, "results", "differential_expression.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# vst()'s default nsub=1000 subsamples genes to fit the mean-dispersion trend;
# with fewer genes than that (small panels, demo/pilot data) it errors out.
# Capping nsub at the available gene count keeps vst() usable at any matrix size.
vsd <- DESeq2::vst(dds, blind = FALSE, nsub = min(1000, nrow(dds)))
norm <- as.data.frame(SummarizedExperiment::assay(vsd))
norm$gene <- rownames(norm)
write.table(norm, file.path(out_dir, "results", "normalized_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

sample_annotation <- data.frame(group = metadata[[group_col]], row.names = metadata$sample)

volcano_plot <- ggplot2::ggplot(res, ggplot2::aes(log2FoldChange, -log10(padj), color = change)) +
  ggplot2::geom_point(alpha = 0.82, size = 2.1, na.rm = TRUE) +
  ggplot2::geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "#A7C6DD", linewidth = 0.45) +
  ggplot2::geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "#A7C6DD", linewidth = 0.45) +
  ggplot2::scale_color_manual(values = c(up = "#E45B5B", down = "#4EA3F1", ns = "#B8C2CC")) +
  ggplot2::labs(
    title = "Volcano plot",
    subtitle = sprintf("padj < %g, |log2FC| >= %g", padj_threshold, log2fc_threshold),
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
pheatmap::pheatmap(SummarizedExperiment::assay(vsd)[top, , drop = FALSE], scale = "row",
                   annotation_col = data.frame(group = metadata[[group_col]], row.names = metadata$sample),
                   color = colorRampPalette(c("#4EA3F1", "white", "#E45B5B"))(101),
                   fontsize = 12, fontsize_row = 8, fontsize_col = 11,
                   filename = file.path(out_dir, "figures", "heatmap_top_degs.png"),
                   width = 8.4, height = 8.0)
pheatmap::pheatmap(SummarizedExperiment::assay(vsd)[top, , drop = FALSE], scale = "row",
                   annotation_col = data.frame(group = metadata[[group_col]], row.names = metadata$sample),
                   color = colorRampPalette(c("#4EA3F1", "white", "#E45B5B"))(101),
                   fontsize = 12, fontsize_row = 8, fontsize_col = 11,
                   filename = file.path(out_dir, "figures", "heatmap_top_degs.pdf"),
                   width = 8.4, height = 8.0)

# QC summary: PCA, sample correlation, library size / detected genes
# (see references/cns-method-guide.md Required Figures / advanced-modules.md).
pca_data <- DESeq2::plotPCA(vsd, intgroup = group_col, returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
pca_plot <- ggplot2::ggplot(pca_data, ggplot2::aes(PC1, PC2, color = .data[[group_col]], label = name)) +
  ggplot2::geom_point(size = 3.4, alpha = 0.9) +
  ggplot2::geom_text(vjust = -0.9, size = 3.4, show.legend = FALSE) +
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

sample_cor <- cor(SummarizedExperiment::assay(vsd), method = "pearson")
pheatmap::pheatmap(sample_cor,
                   annotation_col = sample_annotation,
                   annotation_row = sample_annotation,
                   color = colorRampPalette(c("#4EA3F1", "white", "#E45B5B"))(101),
                   display_numbers = TRUE, number_color = "#173B5F", number_format = "%.2f",
                   fontsize = 12, fontsize_row = 11, fontsize_col = 11,
                   main = "Sample correlation heatmap",
                   filename = file.path(out_dir, "figures", "sample_correlation_heatmap.png"),
                   width = 7.6, height = 7.0)
pheatmap::pheatmap(sample_cor,
                   annotation_col = sample_annotation,
                   annotation_row = sample_annotation,
                   color = colorRampPalette(c("#4EA3F1", "white", "#E45B5B"))(101),
                   display_numbers = TRUE, number_color = "#173B5F", number_format = "%.2f",
                   fontsize = 12, fontsize_row = 11, fontsize_col = 11,
                   main = "Sample correlation heatmap",
                   filename = file.path(out_dir, "figures", "sample_correlation_heatmap.pdf"),
                   width = 7.6, height = 7.0)

lib_size <- colSums(count_values)
detected_genes <- colSums(count_values > 0)
qc_long <- data.frame(
  sample = factor(rep(metadata$sample, 2), levels = metadata$sample),
  group = rep(metadata[[group_col]], 2),
  metric = rep(c("Library size", "Detected genes"), each = nrow(metadata)),
  value = c(lib_size[metadata$sample], detected_genes[metadata$sample])
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
                       library_size = as.numeric(lib_size[metadata$sample]),
                       detected_genes = as.numeric(detected_genes[metadata$sample])),
            file.path(out_dir, "qc", "sample_qc.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# QC verdict: reuse the same sample-correlation matrix rendered above to flag
# (1) any sample whose mean same-group correlation trails its group average
# (possible outlier/mislabel/batch effect), and (2) any extra metadata
# covariate that is perfectly confounded with the group of interest
# (contingency check, pitfalls.md #13).
# ---------------------------------------------------------------------------
outlier_notes <- character(0)
for (g in levels(metadata[[group_col]])) {
  members <- metadata$sample[metadata[[group_col]] == g]
  if (length(members) < 2) {
    outlier_notes <- c(outlier_notes, sprintf(
      "Group '%s' has only %d sample(s); within-group correlation consistency cannot be assessed.", g, length(members)
    ))
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
    if (nrow(tab) < 2) next  # constant covariate (zero variance): no information to check
    levels_per_group <- apply(tab > 0, 1, sum)
    if (all(levels_per_group == 1)) {
      confounded_covariates <- c(confounded_covariates, col)
    }
  }
}

# ---------------------------------------------------------------------------
# Caveats: check the run's actual data against references/pitfalls.md and
# report only the ones triggered -- never a static list of all pitfalls.
# ---------------------------------------------------------------------------
caveats <- character(0)
caveat_titles <- character(0)
if (n_case < 3 || n_control < 3) {
  caveat_titles <- c(caveat_titles, "underpowered design")
  caveats <- c(caveats, sprintf(
    "**Underpowered design** (pitfalls.md #10): n=%d (%s) vs n=%d (%s) is below the recommended minimum of 3 biological replicates per group. Gene-wise dispersion cannot be reliably estimated at this replication level; treat DEG calls as exploratory, not confirmed differential expression.",
    n_case, case_label, n_control, control_label
  ))
}
if (n_case > 30 && n_control > 30) {
  caveat_titles <- c(caveat_titles, "effect-size floor at large n")
  caveats <- c(caveats, sprintf(
    "**Effect-size floor at large n** (pitfalls.md #12): with n=%d/%d per group, statistical significance can be reached for biologically trivial fold-changes. The abs(log2FC) >= %g filter is kept as a hard floor alongside padj < %g; consider ranking by effect size rather than padj alone.",
    n_case, n_control, log2fc_threshold, padj_threshold
  ))
}
for (col in confounded_covariates) {
  caveat_titles <- c(caveat_titles, sprintf("'%s' confounded with '%s'", col, group_col))
  caveats <- c(caveats, sprintf(
    "**Confounded design** (pitfalls.md #13): metadata column '%s' is perfectly confounded with '%s' -- every '%s' value maps to exactly one group. The design matrix cannot separate these two effects; do not attribute observed differences to '%s' alone.",
    col, group_col, col, group_col
  ))
}
if (has_non_integer) {
  caveat_titles <- c(caveat_titles, "non-integer count values")
  caveats <- c(caveats, "**Non-integer count values** (pitfalls.md #11): the input matrix contained non-integer values that were rounded before DESeq2 modeling. If this matrix was already normalized (TPM/FPKM/CPM) rather than raw counts, DESeq2's negative-binomial model is not appropriate -- use limma on log-scale values instead.")
}

up_count <- sum(res$change == "up", na.rm = TRUE)
down_count <- sum(res$change == "down", na.rm = TRUE)
# "Top genes" means significant (up/down) genes, not just small-padj genes -- a gene can have a
# tiny padj but sit below the log2FC floor (change == "ns") and would otherwise misleadingly
# rank #1 in a "top gene" list despite not clearing the effect-size threshold.
sig_res <- res[!is.na(res$padj) & res$change != "ns", ]
top_n <- min(5, nrow(sig_res))
top_genes <- head(sig_res, top_n)

report_lines <- c(
  "# bulkRNA DESeq2 report",
  "",
  sprintf(
    "Formal DESeq2 analysis for **%s** vs **%s** (%s = reference/control), design `~ %s`.",
    case_label, control_label, control_label, group_col
  ),
  "",
  "## Input summary",
  "",
  sprintf("- Samples: %d %s vs %d %s", n_case, case_label, n_control, control_label),
  sprintf(
    "- Genes in input matrix: %d; genes tested after the count>=10-in->=%d-samples filter: %d",
    genes_input, min_group_n, nrow(res)
  ),
  sprintf("- Contrast: %s vs %s; thresholds: padj < %g and abs(log2FC) >= %g (apeglm/shrinkage not applied here)", case_label, control_label, padj_threshold, log2fc_threshold)
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
    report_lines <- c(report_lines, sprintf(
      "  %d. `%s`: log2FC = %.2f, padj = %.3g (%s)",
      i, top_genes$gene[i], top_genes$log2FoldChange[i], top_genes$padj[i], top_genes$change[i]
    ))
  }
}
na_padj_genes <- sum(is.na(res$padj))
if (na_padj_genes > 0) {
  report_lines <- c(report_lines, sprintf(
    "- %d gene(s) have padj = NA (DESeq2 independent filtering by low mean count, or a Cook's-distance outlier flag; see pitfalls.md #4) -- not the same as non-significant, and excluded from the counts above.",
    na_padj_genes
  ))
}

report_lines <- c(report_lines, "", "## QC verdict", "")
if (length(outlier_notes)) {
  report_lines <- c(report_lines, paste0("- ", outlier_notes))
} else {
  report_lines <- c(report_lines, "- No sample fell outside the same-group correlation heuristic; no outlier flagged.")
}
if (length(confounded_covariates)) {
  report_lines <- c(report_lines, paste0(sprintf("- Covariate '%s' is fully confounded with the group of interest (see caveats below).", confounded_covariates)))
} else if (length(extra_covariate_cols)) {
  report_lines <- c(report_lines, "- No additional metadata covariate is fully confounded with the group of interest.")
}

report_lines <- c(report_lines, "", "## Caveats (automatically checked against references/pitfalls.md)", "")
if (length(caveats)) {
  report_lines <- c(report_lines, paste0("- ", caveats))
} else {
  report_lines <- c(report_lines, "- No data-triggered caveats from the automatic checks (design >= 3/group, no full covariate confounding detected, counts are integer-valued). This does not replace manual review of `references/pitfalls.md`.")
}

qc_sentence <- if (length(outlier_notes) == 0 && length(confounded_covariates) == 0) {
  "PCA and sample correlation show no evidence of outlier samples or covariate confounding"
} else {
  "PCA/sample correlation flagged one or more issues noted above"
}
robustness <- "based on a formal DESeq2 negative-binomial model with explicit contrast and BH FDR correction"
caveat_clause <- if (length(caveat_titles)) paste0(", with the following flagged concerns: ", paste(caveat_titles, collapse = ", ")) else ""
report_lines <- c(report_lines, "", "## Conclusion", "",
  sprintf(
    "%d gene(s) were significantly upregulated and %d downregulated in %s vs %s (n=%d/%d) at padj<%g and abs(log2FC)>=%g. %s. Results are %s%s.",
    up_count, down_count, case_label, control_label, n_case, n_control, padj_threshold, log2fc_threshold, qc_sentence, robustness, caveat_clause
  )
)

report_lines <- c(report_lines, "", "## Files", "",
  "Figures: `volcano`, `heatmap_top_degs`, `pca`, `sample_correlation_heatmap`, `qc_summary` (library size and detected genes), each exported as PNG and PDF. Full gene table: `results/differential_expression.tsv`. Normalized (VST) expression: `results/normalized_counts.tsv`. Per-sample QC: `qc/sample_qc.tsv`."
)

writeLines(report_lines, file.path(out_dir, "report.md"))
write.table(data.frame(item = c("case_samples", "control_samples", "genes_input", "genes_filtered_out", "genes_tested", "up_genes", "down_genes", "padj_threshold", "log2fc_threshold", "top_n_heatmap"),
                       value = c(n_case, n_control, genes_input, genes_filtered_out, nrow(res), up_count, down_count, padj_threshold, log2fc_threshold, top_n_heatmap)),
            file.path(out_dir, "summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
