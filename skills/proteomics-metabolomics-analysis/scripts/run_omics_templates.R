#!/usr/bin/env Rscript
# Formal differential proteomics/metabolomics pipeline (limma moderated t-test route).
#
# SCOPE NOTE (read before extending): this script covers the single-omics,
# two-group differential workflow that is realistically portable as one
# script: missingness QC, group-aware feature filtering, limma lmFit/eBayes
# differential testing, and the core CNS-required figures/tables. Per
# references/pitfalls.md and references/cns-method-guide.md, the following
# are OUT OF SCOPE here and must not be faked with a shortcut:
#   - TMT/iTRAQ ratio-compression correction (needs MS2/MS3 co-isolation
#     metadata and per-plex bridge-channel design; see cns-method-guide.md
#     and pitfalls.md #2, #4).
#   - MOFA2/mixOmics(DIABLO) multi-omics integration (needs matched sample
#     IDs across omics layers and a much larger interactive workflow; see
#     references/advanced-modules.md).
#   - Pathway / compound-ID mapping and enrichment (needs a verified
#     name-to-InChIKey-to-KEGG/HMDB crosswalk; see pitfalls.md #10 and
#     references/advanced-modules.md).
# If the user needs any of the above, read references/advanced-modules.md
# and references/cns-method-guide.md and use the dedicated DEP/MSstats/
# MetaboAnalyst/mixOmics packages/objects instead of extending this script.
#
# DATA CONVENTION (see references/input-contract.md): the abundance matrix
# is assumed to be ALREADY NORMALIZED, roughly log2-scale continuous
# intensity (e.g. MaxLFQ log2 intensity, TMT log2 reporter ratio to a
# bridge channel, or a normalized/log-transformed metabolomics matrix).
# This script therefore does NOT re-normalize (no size-factor/VST/TIC/PQN
# step) — applying a count-based normalization on top of already-normalized
# log-scale values would distort them. If the input is raw/non-normalized
# intensity, normalize it first per references/cns-method-guide.md.

args <- commandArgs(trailingOnly = TRUE)

# ---- Optional --flag value overrides (parsed out of `args` before positional
# parsing below, so existing positional usage keeps working unchanged) ----
parse_flag <- function(args, flag, default, as_numeric = TRUE) {
  hit <- which(args == flag)
  value <- default
  if (length(hit) > 0) {
    i <- hit[[1]]
    if (i == length(args)) stop("Flag ", flag, " requires a value.")
    value <- args[[i + 1]]
    if (as_numeric) value <- suppressWarnings(as.numeric(value))
    if (as_numeric && is.na(value)) stop("Flag ", flag, " requires a numeric value.")
    args <<- args[-c(i, i + 1)]
  }
  value
}
opt_padj_threshold <- parse_flag(args, "--padj-threshold", 0.05)
opt_log2fc_threshold <- parse_flag(args, "--log2fc-threshold", 1)
opt_min_group_completeness <- parse_flag(args, "--min-group-completeness", 0.5)
if (!(opt_padj_threshold > 0 && opt_padj_threshold <= 1)) stop("--padj-threshold must be in (0, 1].")
if (opt_log2fc_threshold < 0) stop("--log2fc-threshold must be >= 0.")
if (!(opt_min_group_completeness >= 0 && opt_min_group_completeness <= 1)) stop("--min-group-completeness must be in [0, 1].")

if (length(args) < 5) {
  stop("Usage: Rscript run_omics_templates.R matrix.csv metadata.csv group_col case control [output_dir] ",
       "[--padj-threshold 0.05] [--log2fc-threshold 1] [--min-group-completeness 0.5]")
}
matrix_file <- args[[1]]
metadata_file <- args[[2]]
group_col <- args[[3]]
case_label <- args[[4]]
control_label <- args[[5]]
out_dir <- ifelse(length(args) >= 6, args[[6]], "output")

required <- c("limma", "ggplot2", "pheatmap")
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

# ---- Adjustable thresholds (documented; override via CLI flags shown in the
# Usage message above, or edit the defaults here; the value actually used is
# always reported in report.md/methods.md) ----
MIN_GROUP_COMPLETENESS <- opt_min_group_completeness   # keep a feature only if >=MIN_GROUP_COMPLETENESS
                                 # fraction of samples are non-missing WITHIN
                                 # EACH group (case AND control), per
                                 # cns-method-guide.md's group-aware
                                 # missingness filtering guidance.
                                 # Default 0.5; override with --min-group-completeness.
SIG_PADJ <- opt_padj_threshold   # default 0.05; override with --padj-threshold
SIG_LOG2FC <- opt_log2fc_threshold  # default 1; override with --log2fc-threshold
TOP_N_HEATMAP <- 40              # cap heatmap rows; visual-style.md asks to
                                 # avoid crowded heatmaps unless the user
                                 # explicitly wants all features.
HIGH_MISSING_FLAG <- 0.3         # reference line only (visual flag), not a
                                 # hard cutoff, on the missingness bar chart.
SMALL_GROUP_N <- 3               # below this per-group n, flag pitfalls.md #8
                                 # (multiple-testing burden at small sample size)

## ---- Load and validate inputs ----
mat_raw <- read.csv(matrix_file, check.names = FALSE, na.strings = c("NA", ""))
metadata <- read.csv(metadata_file, check.names = FALSE, na.strings = c("NA", ""))
stopifnot(
  "feature" %in% colnames(mat_raw),
  "sample" %in% colnames(metadata),
  group_col %in% colnames(metadata)
)
if (anyDuplicated(mat_raw$feature)) stop("Duplicate feature IDs in matrix file.")

sample_cols <- setdiff(colnames(mat_raw), "feature")
rownames(mat_raw) <- mat_raw$feature
mat_all <- as.matrix(mat_raw[, sample_cols, drop = FALSE])
storage.mode(mat_all) <- "double"

metadata <- metadata[match(sample_cols, metadata$sample), ]
if (any(is.na(metadata$sample))) {
  missing_ids <- sample_cols[!(sample_cols %in% metadata$sample)]
  stop("Sample IDs in matrix are missing from metadata: ", paste(missing_ids, collapse = ", "))
}
group_raw <- as.character(metadata[[group_col]])
if (!all(c(case_label, control_label) %in% group_raw)) {
  stop("case/control labels not found in metadata column '", group_col, "'.")
}
group <- factor(group_raw, levels = c(control_label, case_label))
if (any(table(group) < 2)) {
  warning("Fewer than 2 samples in one group; differential statistics will be unstable.")
}
samples <- sample_cols
n_control <- sum(group == control_label)
n_case <- sum(group == case_label)

## ---- Missingness QC (per-feature and per-sample) ----
is_na_mat <- is.na(mat_all)
feat_n_missing_control <- rowSums(is_na_mat[, group == control_label, drop = FALSE])
feat_n_missing_case <- rowSums(is_na_mat[, group == case_label, drop = FALSE])
feat_n_missing_total <- rowSums(is_na_mat)
feat_frac_missing_total <- feat_n_missing_total / ncol(mat_all)

# Group-aware completeness rule (see MIN_GROUP_COMPLETENESS above): require
# at least MIN_GROUP_COMPLETENESS fraction of NON-missing values in EACH
# group. This protects against a feature that looks "complete enough"
# overall but is actually all-or-nothing in one condition (pitfalls.md #11).
pass_filter <- ((n_control - feat_n_missing_control) / max(n_control, 1) >= MIN_GROUP_COMPLETENESS) &
  ((n_case - feat_n_missing_case) / max(n_case, 1) >= MIN_GROUP_COMPLETENESS)

feature_missingness <- data.frame(
  scope = "feature",
  id = rownames(mat_all),
  group = NA_character_,
  n_total = ncol(mat_all),
  n_missing = feat_n_missing_total,
  frac_missing = round(feat_frac_missing_total, 4),
  n_missing_control = feat_n_missing_control,
  frac_missing_control = round(feat_n_missing_control / max(n_control, 1), 4),
  n_missing_case = feat_n_missing_case,
  frac_missing_case = round(feat_n_missing_case / max(n_case, 1), 4),
  pass_filter = pass_filter,
  stringsAsFactors = FALSE
)

sample_n_missing <- colSums(is_na_mat)
sample_missingness <- data.frame(
  scope = "sample",
  id = samples,
  group = as.character(group),
  n_total = nrow(mat_all),
  n_missing = sample_n_missing,
  frac_missing = round(sample_n_missing / nrow(mat_all), 4),
  n_missing_control = NA_integer_,
  frac_missing_control = NA_real_,
  n_missing_case = NA_integer_,
  frac_missing_case = NA_real_,
  pass_filter = NA,
  stringsAsFactors = FALSE
)

missingness_summary <- rbind(feature_missingness, sample_missingness)
write.table(missingness_summary, file.path(out_dir, "results", "missingness_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

filtered_mat <- mat_all[pass_filter, , drop = FALSE]
n_dropped <- sum(!pass_filter)
message(sprintf(
  "Missingness filter (>= %.0f%% complete per group): kept %d / %d features, dropped %d.",
  MIN_GROUP_COMPLETENESS * 100, nrow(filtered_mat), nrow(mat_all), n_dropped
))
if (nrow(filtered_mat) < 2) stop("Fewer than 2 features pass the missingness filter; relax MIN_GROUP_COMPLETENESS or inspect input quality.")

write.table(
  data.frame(sample = samples, group = as.character(group), n_missing = sample_n_missing,
             frac_missing = round(sample_n_missing / nrow(mat_all), 4)),
  file.path(out_dir, "qc", "sample_qc.tsv"), sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  data.frame(item = c("features_total", "features_kept", "features_dropped_missingness"),
             value = c(nrow(mat_all), nrow(filtered_mat), n_dropped)),
  file.path(out_dir, "qc", "feature_filtering.tsv"), sep = "\t", row.names = FALSE, quote = FALSE
)

## ---- limma differential analysis ----
# lmFit handles per-row missing values natively (each feature is fit with
# only its own non-missing samples), so no imputation is performed before
# statistical testing — only the group-completeness filter above is applied.
safe_control <- make.names(control_label)
safe_case <- make.names(case_label)
design <- stats::model.matrix(~0 + group)
colnames(design) <- c(safe_control, safe_case)

fit <- limma::lmFit(filtered_mat, design)
contrast_expr <- paste0(safe_case, " - ", safe_control)
cont.matrix <- limma::makeContrasts(contrasts = contrast_expr, levels = design)
fit2 <- limma::contrasts.fit(fit, cont.matrix)
fit2 <- limma::eBayes(fit2)

tt <- limma::topTable(fit2, coef = 1, number = Inf, adjust.method = "BH", sort.by = "P")
tt$feature <- rownames(tt)
tt$change <- ifelse(!is.na(tt$adj.P.Val) & tt$adj.P.Val < SIG_PADJ & tt$logFC >= SIG_LOG2FC, "up",
  ifelse(!is.na(tt$adj.P.Val) & tt$adj.P.Val < SIG_PADJ & tt$logFC <= -SIG_LOG2FC, "down", "ns"))
res <- tt[, c("feature", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "change")]
colnames(res) <- c("feature", "log2FC", "AveExpr", "t", "pvalue", "padj", "B", "change")
res <- res[order(res$padj), ]
write.table(res, file.path(out_dir, "results", "differential_features.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

n_up <- sum(res$change == "up")
n_down <- sum(res$change == "down")

## ---- Derived stats for the report (missingness overview, top named
## features, small-n flag) — computed once here so the report block below
## stays declarative. ----
n_high_miss_features <- sum(feature_missingness$frac_missing > HIGH_MISSING_FLAG)
high_miss_samples <- sample_missingness$id[sample_missingness$frac_missing > HIGH_MISSING_FLAG]
small_groups <- c(control_label, case_label)[c(n_control, n_case) < SMALL_GROUP_N]
top_report <- head(res[order(res$padj), ], 10)

## ---- Figures ----
# Palette matches assets/palette.json and references/visual-style.md.
COL_UP <- "#E45B5B"
COL_DOWN <- "#4EA3F1"
COL_NS <- "#B8C2CC"
COL_TEXT <- "#173B5F"
COL_GRID <- "#DCEAF5"
COL_AMBER <- "#E2A13B"

base_theme <- ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 24, color = COL_TEXT),
    axis.title = ggplot2::element_text(size = 21, color = COL_TEXT),
    axis.text = ggplot2::element_text(size = 16, color = COL_TEXT),
    legend.position = "top"
  )

# 1) PCA — center-only (no per-feature scaling) since the input is already a
# common log2-scale unit across features; use only features with ZERO
# missing values so the PCA is not distorted by any imputation choice.
complete_idx <- stats::complete.cases(filtered_mat)
if (sum(complete_idx) < 3) {
  warning("Fewer than 3 fully-observed features; skipping PCA plot.")
} else {
  pca_mat <- filtered_mat[complete_idx, , drop = FALSE]
  pca <- stats::prcomp(t(pca_mat), center = TRUE, scale. = FALSE)
  var_pct <- round(100 * (pca$sdev^2) / sum(pca$sdev^2), 1)
  pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], sample = samples, group = as.character(group))
  pca_plot <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, color = group, label = sample)) +
    ggplot2::geom_point(size = 3.2, alpha = 0.88) +
    ggplot2::geom_text(size = 3.6, vjust = -0.9, show.legend = FALSE, check_overlap = TRUE, color = COL_TEXT) +
    ggplot2::scale_color_manual(values = setNames(c(COL_DOWN, COL_UP), c(control_label, case_label))) +
    ggplot2::labs(title = "PCA", x = paste0("PC1 (", var_pct[1], "%)"), y = paste0("PC2 (", var_pct[2], "%)"), color = NULL) +
    base_theme
  save_gg(pca_plot, "pca", width = 7.8, height = 6.4)
}

# 2) Missingness bar chart (per sample), grouped/colored by condition.
miss_plot_df <- data.frame(sample = samples, group = as.character(group),
                            frac_missing = round(sample_n_missing / nrow(mat_all), 4))
miss_plot_df$sample <- factor(miss_plot_df$sample, levels = samples[order(group)])
missingness_plot <- ggplot2::ggplot(miss_plot_df, ggplot2::aes(sample, frac_missing, fill = group)) +
  ggplot2::geom_col(width = 0.68, alpha = 0.92) +
  ggplot2::geom_hline(yintercept = HIGH_MISSING_FLAG, linetype = "dashed", color = COL_AMBER, linewidth = 0.5) +
  ggplot2::scale_fill_manual(values = setNames(c(COL_DOWN, COL_UP), c(control_label, case_label))) +
  ggplot2::labs(title = "Per-sample missingness", x = NULL, y = "Fraction missing features", fill = NULL) +
  base_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
save_gg(missingness_plot, "missingness", width = max(7.4, length(samples) * 0.42 + 3), height = 5.8)

# 3) Sample correlation heatmap (pairwise-complete Pearson correlation;
# no imputation needed for a correlation matrix).
cor_mat <- stats::cor(filtered_mat, use = "pairwise.complete.obs", method = "pearson")
ann_col <- data.frame(group = group, row.names = samples)
# Match the annotation bar colors to the same case/control palette used in
# PCA and the volcano plot (COL_DOWN = control, COL_UP = case), so group
# color-coding is consistent across every figure in this report.
ann_colors <- list(group = setNames(c(COL_DOWN, COL_UP), c(control_label, case_label)))
pheatmap::pheatmap(cor_mat,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  color = grDevices::colorRampPalette(c("white", COL_DOWN, COL_TEXT))(101),
  display_numbers = ncol(cor_mat) <= 20, number_color = "black", fontsize_number = 8,
  fontsize = 12, fontsize_row = 10, fontsize_col = 10,
  main = "Sample correlation",
  filename = file.path(out_dir, "figures", "sample_correlation.png"),
  width = 7.6, height = 6.8
)
pheatmap::pheatmap(cor_mat,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  color = grDevices::colorRampPalette(c("white", COL_DOWN, COL_TEXT))(101),
  display_numbers = ncol(cor_mat) <= 20, number_color = "black", fontsize_number = 8,
  fontsize = 12, fontsize_row = 10, fontsize_col = 10,
  main = "Sample correlation",
  filename = file.path(out_dir, "figures", "sample_correlation.pdf"),
  width = 7.6, height = 6.8
)

# 4) Volcano plot (all tested features).
volcano_plot <- ggplot2::ggplot(res, ggplot2::aes(log2FC, -log10(pmax(padj, 1e-12)), color = change)) +
  ggplot2::geom_point(alpha = 0.82, size = 2.1, na.rm = TRUE) +
  ggplot2::geom_vline(xintercept = c(-SIG_LOG2FC, SIG_LOG2FC), linetype = "dashed", color = "#A7C6DD", linewidth = 0.45) +
  ggplot2::geom_hline(yintercept = -log10(SIG_PADJ), linetype = "dashed", color = "#A7C6DD", linewidth = 0.45) +
  ggplot2::scale_color_manual(values = c(up = COL_UP, down = COL_DOWN, ns = COL_NS)) +
  ggplot2::labs(title = sprintf("Differential features (padj<%g, |log2FC|>=%g)", SIG_PADJ, SIG_LOG2FC),
                x = "log2 fold change", y = "-log10 adjusted p value", color = NULL) +
  base_theme
save_gg(volcano_plot, "volcano", width = 8.8, height = 6.6)

# 5) Top-feature heatmap (capped at TOP_N_HEATMAP; visualization-only
# row-mean fill for any residual NA cell within the selected subset — NEVER
# used for statistics above; flagged if it happens).
top_features <- head(res$feature[order(res$padj)], TOP_N_HEATMAP)
top_mat <- filtered_mat[top_features, , drop = FALSE]
resid_na <- sum(is.na(top_mat))
if (resid_na > 0) {
  message(sprintf("Note: %d residual missing cell(s) among top %d heatmap features filled with the per-feature row mean for display only (not used in statistics).", resid_na, length(top_features)))
  row_means <- rowMeans(top_mat, na.rm = TRUE)
  na_idx <- which(is.na(top_mat), arr.ind = TRUE)
  top_mat[na_idx] <- row_means[na_idx[, 1]]
}
z <- t(scale(t(top_mat)))
pheatmap::pheatmap(z,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  color = grDevices::colorRampPalette(c(COL_DOWN, "white", COL_UP))(101),
  breaks = seq(-2.2, 2.2, length.out = 102),
  fontsize = 12, fontsize_row = 8, fontsize_col = 10,
  main = "Top feature heatmap",
  filename = file.path(out_dir, "figures", "top_feature_heatmap.png"),
  width = 8.4, height = max(6.2, length(top_features) * 0.18 + 3)
)
pheatmap::pheatmap(z,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  color = grDevices::colorRampPalette(c(COL_DOWN, "white", COL_UP))(101),
  breaks = seq(-2.2, 2.2, length.out = 102),
  fontsize = 12, fontsize_row = 8, fontsize_col = 10,
  main = "Top feature heatmap",
  filename = file.path(out_dir, "figures", "top_feature_heatmap.pdf"),
  width = 8.4, height = max(6.2, length(top_features) * 0.18 + 3)
)

## ---- Report and summary ----

top_lines <- if (nrow(top_report) > 0) {
  sprintf("- %s: log2FC=%.2f, padj=%.2e (%s)", top_report$feature, top_report$log2FC, top_report$padj, top_report$change)
} else {
  "- No features available (empty differential table)."
}

caveats <- character(0)
if (length(high_miss_samples) > 0) {
  caveats <- c(caveats, sprintf(
    "- %d sample(s) exceed %d%% missing features (%s); inspect these samples for run quality before trusting features that are largely observed only in the other samples (references/pitfalls.md #4, #11).",
    length(high_miss_samples), HIGH_MISSING_FLAG * 100, paste(high_miss_samples, collapse = ", ")
  ))
}
if (n_high_miss_features > 0) {
  caveats <- c(caveats, sprintf(
    "- %d feature(s) exceed %d%% missing values; check whether their missingness is concentrated in one group (possible MNAR / detection-limit signal, references/pitfalls.md #1) before ranking them by fold change.",
    n_high_miss_features, HIGH_MISSING_FLAG * 100
  ))
}
caveats <- c(caveats,
  "- Missingness handling here is filter-only (group-aware completeness filter); no MNAR/MAR imputation was performed before statistical testing (references/pitfalls.md #1, references/cns-method-guide.md).",
  "- TMT/iTRAQ ratio-compression correction (references/pitfalls.md #2) and multi-omics integration (MOFA2/mixOmics) require dedicated, non-portable workflows not run here — see references/advanced-modules.md."
)
if (length(small_groups) > 0) {
  group_n_txt <- paste(sprintf("%s=%d", c(control_label, case_label), c(n_control, n_case))[c(control_label, case_label) %in% small_groups], collapse = ", ")
  caveats <- c(caveats, sprintf(
    "- Group(s) with n < %d (%s): per-feature variance estimates are unstable at this sample size; moderated eBayes shrinkage (already used here) partially compensates, but treat results as preliminary (references/pitfalls.md #8).",
    SMALL_GROUP_N, group_n_txt
  ))
}

if (n_up + n_down == 0) {
  conclusion <- sprintf(
    "Of %d tested features across %d samples (%s vs %s), none met the log2FC >= %s / padj < %s threshold in either direction. This may reflect a true absence of a strong group effect, insufficient statistical power at the current sample size, or a threshold that is too strict for this dataset — inspect the volcano plot and consider whether the study is adequately powered before concluding there is no effect.",
    nrow(res), length(samples), case_label, control_label, SIG_LOG2FC, SIG_PADJ
  )
} else {
  top1 <- top_report[1, ]
  direction_txt <- if (identical(top1$change, "up")) paste0("higher in ", case_label) else if (identical(top1$change, "down")) paste0("lower in ", case_label) else "inconclusive direction"
  conclusion <- sprintf(
    "Of %d tested features across %d samples (%s vs %s), %d were up and %d were down (log2FC >= %s, padj < %s); the strongest signal was %s (log2FC=%.2f, padj=%.2e, %s). These results describe association/direction of change only under the filter-only missingness handling above; see the caveats before treating any single feature's fold-change magnitude as a precisely measured, imputation-independent effect, and confirm the annotation/identification confidence of any feature central to the conclusion.",
    nrow(res), length(samples), case_label, control_label, n_up, n_down, SIG_LOG2FC, SIG_PADJ, top1$feature, top1$log2FC, top1$padj, direction_txt
  )
}

writeLines(c(
  "# Proteomics/Metabolomics limma differential report",
  "",
  "## Input matrix summary",
  paste0("- Features in input matrix: ", nrow(mat_all)),
  paste0("- Samples: ", length(samples), " (", control_label, " = ", n_control, ", ", case_label, " = ", n_case, ")"),
  paste0("- Data convention: input treated as already-normalized, log2-scale intensity (no re-normalization performed)."),
  "",
  "## Normalization and missingness",
  paste0("- Missingness filter: keep features with >= ", MIN_GROUP_COMPLETENESS * 100, "% non-missing values in EACH group."),
  paste0("- Features kept: ", nrow(filtered_mat), " / ", nrow(mat_all), " (dropped ", n_dropped, " for excess missingness)."),
  paste0("- No imputation was used for statistical testing; limma::lmFit fits each feature on its own non-missing samples."),
  paste0("- A per-feature row-mean fill was used ONLY for display in the top-feature heatmap (", resid_na, " residual cell(s)); never for statistics."),
  "",
  "## Missingness overview",
  paste0("- Features flagged high-missingness (>", HIGH_MISSING_FLAG * 100, "% missing across samples): ", n_high_miss_features, " / ", nrow(mat_all),
         " (", sprintf("%.1f", 100 * n_high_miss_features / max(nrow(mat_all), 1)), "%)"),
  paste0("- Samples flagged high-missingness (>", HIGH_MISSING_FLAG * 100, "% missing features): ", length(high_miss_samples), " / ", length(samples),
         if (length(high_miss_samples) > 0) paste0(" (", paste(high_miss_samples, collapse = ", "), ")") else ""),
  "",
  "## Differential results",
  paste0("- Contrast: ", case_label, " vs ", control_label, " (limma lmFit + eBayes moderated t-test, BH-FDR)."),
  paste0("- Features tested: ", nrow(res)),
  paste0("- Significant up (log2FC >= ", SIG_LOG2FC, ", padj < ", SIG_PADJ, "): ", n_up),
  paste0("- Significant down (log2FC <= -", SIG_LOG2FC, ", padj < ", SIG_PADJ, "): ", n_down),
  "",
  "### Top differential features (by padj)",
  top_lines,
  "",
  "## Pathway/multi-omics interpretation",
  "- Not run by this script. Pathway/compound-ID mapping, TMT ratio-compression correction, and MOFA2/mixOmics multi-omics integration require dedicated, non-portable workflows — see references/advanced-modules.md and references/cns-method-guide.md before attempting these.",
  "",
  "## Caveats",
  caveats,
  "",
  "## Conclusion",
  conclusion,
  "",
  "## Limitations",
  "- Two-group contrast only; no batch/covariate term is modeled in this run (see references/pitfalls.md #4 if samples span multiple batches/plexes/runs).",
  "- Missingness handling here is filter-only (no MNAR/MAR imputation); see references/cns-method-guide.md's imputation decision tree and references/advanced-modules.md if imputation is required for your downstream use case.",
  "- Feature identity/annotation confidence (e.g., MSI level for metabolites, razor-peptide protein grouping) is not evaluated by this script; carry that check through separately before biological interpretation."
), file.path(out_dir, "report.md"))

write.table(
  data.frame(item = c("features_total", "features_kept", "features_dropped_missingness", "samples", "up_features", "down_features"),
             value = c(nrow(mat_all), nrow(filtered_mat), n_dropped, length(samples), n_up, n_down)),
  file.path(out_dir, "summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE
)

message("Done. Wrote outputs to ", out_dir)
