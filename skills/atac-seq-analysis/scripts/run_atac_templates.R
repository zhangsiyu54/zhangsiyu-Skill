#!/usr/bin/env Rscript
# Bulk ATAC-seq differential accessibility pipeline (DESeq2-based).
#
# Scope (see references/input-contract.md, references/cns-method-guide.md,
# references/pitfalls.md):
# - Input is a peak-by-sample COUNT MATRIX (first column `peak`, remaining
#   columns sample IDs) plus SAMPLE METADATA (`sample` + group column), i.e.
#   the "existing outputs" / bulk-count route of the input contract.
# - Genome build, peak-calling, consensus/fixed-width peak construction, and
#   ENCODE blacklist filtering are the CALLER'S responsibility upstream of
#   this script (pitfalls #1, #11, #12) -- this script does not fabricate
#   genome-specific annotation logic and does not require TxDb/BSgenome
#   packages. Peak-to-gene annotation (ChIPseeker/HOMER annotatePeak) is out
#   of scope here; if requested, run it as a separate step with an explicit
#   genome build/annotation database and report it as a nearest-gene
#   hypothesis, not a proven target (pitfall #3).
# - scATAC workflows (Signac/ArchR: TF-IDF/LSI, UMAP, clustering, gene
#   activity, motif deviation with chromVAR) are separate, much larger,
#   non-portable, interactive workflows that depend on genome-specific
#   objects (fragment files, TxDb/BSgenome, Signac/ArchR project state) and
#   are intentionally NOT reproduced in this self-contained script; run them
#   interactively in their native Signac/ArchR session instead.
#
# Usage:
#   Rscript run_atac_templates.R peak_counts.csv metadata.csv group_col case control [output_dir] [motif_table] [options...]
#
# Options (optional, order-independent, may appear after the positional args):
#   --padj-threshold=<num>   Adjusted p-value cutoff for calling a peak significant (default: 0.05)
#   --log2fc-threshold=<num> Minimum |log2FoldChange| for calling a peak significant (default: 1)
#   --top-n-heatmap=<int>    Number of top peaks shown in the top-peak accessibility heatmap (default: 50)
#   --top-n-motifs=<int>     Number of top motifs shown in the motif enrichment barplot (default: 12)
#
# Example:
#   Rscript run_atac_templates.R peak_counts.csv metadata.csv group treatment control output motif_enrichment.csv --padj-threshold=0.01 --log2fc-threshold=1.5

raw_args <- commandArgs(trailingOnly = TRUE)
is_flag <- grepl("^--[A-Za-z0-9_-]+=", raw_args)
opt_args <- raw_args[is_flag]
args <- raw_args[!is_flag]

parse_opt <- function(opt_args, name, default, type = "numeric") {
  pattern <- paste0("^--", name, "=")
  hit <- opt_args[grepl(pattern, opt_args)]
  if (length(hit) == 0) {
    return(default)
  }
  value <- sub(pattern, "", hit[length(hit)])
  if (type == "numeric") as.numeric(value) else as.integer(value)
}

if (length(args) < 5) {
  stop("Usage: Rscript run_atac_templates.R peak_counts.csv metadata.csv group_col case control [output_dir] [motif_table.csv] [--padj-threshold=0.05] [--log2fc-threshold=1] [--top-n-heatmap=50] [--top-n-motifs=12]")
}
counts_file <- args[[1]]
metadata_file <- args[[2]]
group_col <- args[[3]]
case_label <- args[[4]]
control_label <- args[[5]]
out_dir <- ifelse(length(args) >= 6, args[[6]], "output")
motif_file <- ifelse(length(args) >= 7 && nzchar(args[[7]]), args[[7]], NA)

padj_threshold <- parse_opt(opt_args, "padj-threshold", 0.05)
log2fc_threshold <- parse_opt(opt_args, "log2fc-threshold", 1)
top_n_heatmap <- parse_opt(opt_args, "top-n-heatmap", 50, type = "integer")
top_n_motifs <- parse_opt(opt_args, "top-n-motifs", 12, type = "integer")

required <- c("DESeq2", "ggplot2", "pheatmap")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "), "\nInstall Bioconductor/CRAN dependencies before running formal analysis.")
}
shrink_type <- if (requireNamespace("apeglm", quietly = TRUE)) {
  "apeglm"
} else if (requireNamespace("ashr", quietly = TRUE)) {
  "ashr"
} else {
  NA
}
if (is.na(shrink_type)) {
  message("Neither apeglm nor ashr is installed; log2FoldChange will NOT be shrunk. ",
          "Install one of them for CNS-standard shrinkage (see references/cns-method-guide.md decision tree).")
}

dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "qc"), recursive = TRUE, showWarnings = FALSE)

# Palette from assets/palette.json (opening/closing/neutral/text/grid).
OPEN <- "#E45B5B"
CLOSE <- "#4EA3F1"
NS <- "#B8C2CC"
TEXT <- "#173B5F"
GRID <- "#DCEAF5"

atac_theme <- ggplot2::theme_classic(base_size = 16) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 24, color = TEXT),
    axis.title = ggplot2::element_text(size = 21, color = TEXT),
    axis.text = ggplot2::element_text(size = 16, color = TEXT),
    legend.text = ggplot2::element_text(size = 15),
    legend.title = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = GRID),
    plot.margin = ggplot2::margin(10, 14, 10, 10)
  )

save_gg <- function(plot, name, width = 8.4, height = 6.4) {
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(name, ".png")), plot, width = width, height = height, dpi = 320, bg = "white")
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(name, ".pdf")), plot, width = width, height = height, bg = "white")
}

message("NOTE: genome build, peak-calling parameters and blacklist filtering status ",
        "are assumed to already be validated upstream (see references/input-contract.md). ",
        "This script only models the supplied peak count matrix.")

## ---- Load and validate inputs -------------------------------------------

counts <- read.csv(counts_file, check.names = FALSE)
metadata <- read.csv(metadata_file, check.names = FALSE)
stopifnot(
  "peak" %in% colnames(counts),
  "sample" %in% colnames(metadata),
  group_col %in% colnames(metadata)
)
rownames(counts) <- counts$peak
counts$peak <- NULL
metadata <- metadata[match(colnames(counts), metadata$sample), ]
if (any(is.na(metadata$sample))) {
  stop("Sample IDs in peak count matrix do not match sample IDs in metadata.")
}
if (!all(c(case_label, control_label) %in% unique(as.character(metadata[[group_col]])))) {
  stop("Case/control labels not found in metadata column '", group_col, "'.")
}
metadata[[group_col]] <- factor(metadata[[group_col]], levels = c(control_label, case_label))

## ---- Sample QC (matrix-derived only) -------------------------------------
# Library depth, FRiP, duplicate rate, insert-size periodicity and TSS
# enrichment require raw BAM/fragment-level data and cannot be recomputed
# from a peak count matrix alone (see references/cns-method-guide.md QC
# list and pitfall #2); report what the matrix can support and flag the rest
# as caller-supplied upstream QC.
count_mat <- as.matrix(counts)
mode(count_mat) <- "numeric"
sample_qc <- data.frame(
  sample = colnames(count_mat),
  group = as.character(metadata[[group_col]]),
  total_count = colSums(count_mat),
  peaks_detected = colSums(count_mat > 0),
  frac_peaks_detected = round(colSums(count_mat > 0) / nrow(count_mat), 4)
)
write.table(sample_qc, file.path(out_dir, "qc", "sample_qc.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## ---- Peak-width bias diagnostic (pitfall #1) ------------------------------
# If peak IDs look like "chr:start-end", check whether width correlates with
# discovered log2FC later; this flags the classic ATAC-on-DESeq2 pitfall of
# wider peaks systematically ranking as more "significant".
peak_ids <- rownames(count_mat)
width_match <- regmatches(peak_ids, regexec("^[^:]+:([0-9]+)-([0-9]+)$", peak_ids))
peak_width <- vapply(width_match, function(m) {
  if (length(m) == 3) as.numeric(m[3]) - as.numeric(m[2]) else NA_real_
}, numeric(1))
names(peak_width) <- peak_ids

## ---- DESeq2 differential accessibility ------------------------------------

dds <- DESeq2::DESeqDataSetFromMatrix(round(count_mat), metadata, as.formula(paste("~", group_col)))
dds <- tryCatch(
  DESeq2::DESeq(dds),
  error = function(e) {
    message("Default parametric dispersion fit failed (", conditionMessage(e),
            "); retrying with fitType = 'local' (common with a small number of peaks).")
    DESeq2::DESeq(dds, fitType = "local")
  }
)

raw_res <- DESeq2::results(dds, contrast = c(group_col, case_label, control_label))
res_df <- as.data.frame(raw_res)
res_df$peak <- rownames(res_df)

if (!is.na(shrink_type)) {
  result_names <- DESeq2::resultsNames(dds)
  coef_name <- result_names[result_names != "Intercept"][1]
  shrunk <- tryCatch({
    if (shrink_type == "apeglm") {
      DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm", res = raw_res)
    } else {
      DESeq2::lfcShrink(dds, contrast = c(group_col, case_label, control_label), type = "ashr", res = raw_res)
    }
  }, error = function(e) {
    message("lfcShrink failed (", conditionMessage(e), "); falling back to unshrunken log2FoldChange.")
    NULL
  })
  if (!is.null(shrunk)) {
    shrunk_df <- as.data.frame(shrunk)
    res_df$log2FoldChange <- shrunk_df$log2FoldChange
    res_df$lfcSE <- shrunk_df$lfcSE
    shrink_used <- shrink_type
  } else {
    shrink_used <- "none"
  }
} else {
  shrink_used <- "none"
}

res_df$change <- ifelse(!is.na(res_df$padj) & res_df$padj < padj_threshold & res_df$log2FoldChange >= log2fc_threshold, "opening",
  ifelse(!is.na(res_df$padj) & res_df$padj < padj_threshold & res_df$log2FoldChange <= -log2fc_threshold, "closing", "ns"))
res_df <- res_df[order(res_df$padj), c("peak", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "change")]
write.table(res_df, file.path(out_dir, "results", "differential_peaks.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

width_note <- "Peak IDs were not parseable as chr:start-end; peak-width bias diagnostic skipped."
if (sum(!is.na(peak_width)) >= 10) {
  merged_width <- peak_width[res_df$peak]
  keep <- !is.na(merged_width) & !is.na(res_df$log2FoldChange)
  if (sum(keep) >= 10) {
    width_cor <- suppressWarnings(cor(merged_width[keep], res_df$log2FoldChange[keep], method = "spearman"))
    if (is.na(width_cor)) {
      width_note <- "Peak-width vs. log2FC check skipped: peak widths or log2FC values had no variance (e.g. all peaks the same fixed width)."
    } else if (abs(width_cor) > 0.3) {
      width_note <- sprintf(
        "Peak-width vs. log2FoldChange Spearman r = %.3f. WARNING: peak width correlates with log2FC beyond the |r|>0.3 rule of thumb -- consider fixed-width summit-centered peaks or a peak-width covariate (see references/pitfalls.md #1).",
        width_cor
      )
    } else {
      width_note <- sprintf(
        "Peak-width vs. log2FoldChange Spearman r = %.3f. No strong peak-width vs. log2FC trend detected (|r|<=0.3).",
        width_cor
      )
    }
  }
}

## ---- VST for QC transforms and figures ------------------------------------
# vst() needs enough peaks to fit its fast dispersion trend (nsub, default
# 1000); fall back to the exact varianceStabilizingTransformation for small
# peak sets (e.g. a demo/subset run) as DESeq2 itself recommends.
vsd <- if (nrow(count_mat) >= 1000) {
  DESeq2::vst(dds, blind = FALSE)
} else {
  DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
}
norm_mat <- SummarizedExperiment::assay(vsd)
norm_out <- as.data.frame(norm_mat)
norm_out$peak <- rownames(norm_out)
norm_out <- norm_out[, c("peak", setdiff(colnames(norm_out), "peak"))]
write.table(norm_out, file.path(out_dir, "results", "normalized_accessibility.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## ---- Figure 1: differential peak volcano ----------------------------------

volcano_plot <- ggplot2::ggplot(res_df, ggplot2::aes(log2FoldChange, -log10(padj), color = change)) +
  ggplot2::geom_point(alpha = 0.82, size = 2.1, na.rm = TRUE) +
  ggplot2::geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "#A7C6DD", linewidth = 0.45) +
  ggplot2::geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "#A7C6DD", linewidth = 0.45) +
  ggplot2::scale_color_manual(values = c(opening = OPEN, closing = CLOSE, ns = NS)) +
  ggplot2::labs(
    title = sprintf("Differential accessibility (padj<%s, |log2FC|>=%s)", padj_threshold, log2fc_threshold),
    x = "log2 fold change", y = "-log10 adjusted p value"
  ) +
  atac_theme
save_gg(volcano_plot, "differential_peak_volcano", width = 8.8, height = 6.6)

## ---- Figure 2: PCA ---------------------------------------------------------

pca_data <- DESeq2::plotPCA(vsd, intgroup = group_col, returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
pca_colors <- stats::setNames(c(CLOSE, OPEN), c(control_label, case_label))
pca_plot <- ggplot2::ggplot(pca_data, ggplot2::aes(PC1, PC2, color = .data[[group_col]], label = name)) +
  ggplot2::geom_point(size = 3.6, alpha = 0.9) +
  ggplot2::geom_text(size = 3.2, vjust = -1.1, show.legend = FALSE, color = TEXT) +
  ggplot2::scale_color_manual(values = pca_colors) +
  ggplot2::labs(
    title = "Sample PCA (VST accessibility)",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  atac_theme
save_gg(pca_plot, "pca", width = 8.6, height = 6.4)

## ---- Figure 3: sample correlation heatmap ---------------------------------

cor_mat <- cor(norm_mat, method = "pearson")
group_annotation <- data.frame(group = metadata[[group_col]], row.names = metadata$sample)
pheatmap::pheatmap(
  cor_mat,
  color = colorRampPalette(c("white", CLOSE))(100),
  annotation_col = group_annotation,
  display_numbers = ncol(cor_mat) <= 12,
  number_color = TEXT,
  fontsize = 12, fontsize_row = 11, fontsize_col = 11,
  main = "Sample correlation (VST accessibility)",
  filename = file.path(out_dir, "figures", "sample_correlation_heatmap.png"),
  width = 7.8, height = 6.8
)
pheatmap::pheatmap(
  cor_mat,
  color = colorRampPalette(c("white", CLOSE))(100),
  annotation_col = group_annotation,
  display_numbers = ncol(cor_mat) <= 12,
  number_color = TEXT,
  fontsize = 12, fontsize_row = 11, fontsize_col = 11,
  main = "Sample correlation (VST accessibility)",
  filename = file.path(out_dir, "figures", "sample_correlation_heatmap.pdf"),
  width = 7.8, height = 6.8
)

## ---- Figure 4: top-peak accessibility heatmap -----------------------------

top_peaks <- head(na.omit(res_df$peak[order(res_df$padj)]), top_n_heatmap)
pheatmap::pheatmap(
  norm_mat[top_peaks, , drop = FALSE], scale = "row",
  annotation_col = group_annotation,
  color = colorRampPalette(c(CLOSE, "white", OPEN))(101),
  fontsize = 12, fontsize_row = 8, fontsize_col = 11,
  main = "Top differential peak accessibility (row z-score)",
  filename = file.path(out_dir, "figures", "peak_heatmap.png"),
  width = 8.4, height = 8.4
)
pheatmap::pheatmap(
  norm_mat[top_peaks, , drop = FALSE], scale = "row",
  annotation_col = group_annotation,
  color = colorRampPalette(c(CLOSE, "white", OPEN))(101),
  fontsize = 12, fontsize_row = 8, fontsize_col = 11,
  main = "Top differential peak accessibility (row z-score)",
  filename = file.path(out_dir, "figures", "peak_heatmap.pdf"),
  width = 8.4, height = 8.4
)

## ---- Optional: motif enrichment pass-through + barplot --------------------

motif_note <- "No motif enrichment table supplied; motif figures/tables skipped."
motif_n <- 0
top_motifs_report <- NULL
motif_sort_col_report <- NA
motif_bg_matched <- FALSE
if (!is.na(motif_file)) {
  if (!file.exists(motif_file)) {
    stop("Motif table not found: ", motif_file)
  }
  motif_sep <- if (grepl("\\.tsv$", motif_file, ignore.case = TRUE)) "\t" else ","
  motifs <- read.table(motif_file, header = TRUE, sep = motif_sep, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(motifs) == 0) {
    stop("Motif table '", motif_file, "' has no data rows.")
  }
  stopifnot("motif" %in% colnames(motifs))
  sort_col <- if ("padj" %in% colnames(motifs)) "padj" else if ("pvalue" %in% colnames(motifs)) "pvalue" else {
    stop("Motif table must contain a 'padj' or 'pvalue' column (see references/input-contract.md).")
  }
  motifs <- motifs[order(motifs[[sort_col]]), ]
  write.table(motifs, file.path(out_dir, "results", "motif_enrichment.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  top_motifs <- head(motifs, top_n_motifs)
  label_col <- if ("tf" %in% colnames(top_motifs)) "tf" else "motif"
  top_motifs$neglog10 <- -log10(pmax(top_motifs[[sort_col]], 1e-12))
  top_motifs[[label_col]] <- factor(top_motifs[[label_col]], levels = rev(top_motifs[[label_col]]))
  motif_plot <- ggplot2::ggplot(top_motifs, ggplot2::aes(.data[[label_col]], neglog10)) +
    ggplot2::geom_col(fill = CLOSE, alpha = 0.86, width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Motif enrichment", x = NULL, y = paste0("-log10 adjusted ", sort_col)) +
    atac_theme
  save_gg(motif_plot, "motif_enrichment_barplot", width = 8.2, height = max(4.8, nrow(top_motifs) * 0.42 + 2))
  motif_note <- sprintf("Motif enrichment table with %d motifs written; top %d plotted by %s.", nrow(motifs), nrow(top_motifs), sort_col)
  motif_n <- nrow(motifs)
  # Keep the ranking (by raw motif table order, i.e. sorted ascending by
  # sort_col) for the plain-language top-motif list in report.md.
  top_motifs_report <- motifs[seq_len(min(5, nrow(motifs))), c(label_col, sort_col)]
  colnames(top_motifs_report) <- c("label", "sort_value")
  motif_sort_col_report <- sort_col
  bg_col <- intersect(c("background", "background_matched"), colnames(motifs))
  motif_bg_matched <- length(bg_col) > 0 && all(tolower(as.character(motifs[[bg_col[1]]])) %in% c("1", "true", "yes", "matched"))
}

## ---- Summary, report, methods text ----------------------------------------

n_opening <- sum(res_df$change == "opening", na.rm = TRUE)
n_closing <- sum(res_df$change == "closing", na.rm = TRUE)

write.table(
  data.frame(
    item = c("peaks_tested", "opening_peaks", "closing_peaks", "samples", "shrinkage_method", "motif_count", "padj_threshold", "log2fc_threshold", "top_n_heatmap", "top_n_motifs"),
    value = c(nrow(res_df), n_opening, n_closing, ncol(count_mat), shrink_used, motif_n, padj_threshold, log2fc_threshold, top_n_heatmap, top_n_motifs)
  ),
  file.path(out_dir, "summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE
)

deseq2_version <- as.character(utils::packageVersion("DESeq2"))
r_version <- paste(R.version$major, R.version$minor, sep = ".")

## ---- Plain-language top-peak / top-motif lists ----------------------------

n_sig <- n_opening + n_closing
sig_frac <- if (nrow(res_df) > 0) n_sig / nrow(res_df) else 0

top_peaks_report <- head(res_df[order(res_df$padj), c("peak", "log2FoldChange", "padj", "change")], 5)
top_peak_lines <- if (nrow(top_peaks_report) == 0) {
  "No peaks were tested."
} else {
  vapply(seq_len(nrow(top_peaks_report)), function(i) {
    row <- top_peaks_report[i, ]
    sprintf("%d. `%s` -- log2FC = %.3g, padj = %.3g (%s)", i, row$peak, row$log2FoldChange, row$padj, row$change)
  }, character(1))
}

top_motif_lines <- if (is.null(top_motifs_report)) {
  "No motif enrichment table supplied."
} else {
  vapply(seq_len(nrow(top_motifs_report)), function(i) {
    row <- top_motifs_report[i, ]
    sprintf("%d. **%s** -- %s = %.3g", i, row$label, motif_sort_col_report, row$sort_value)
  }, character(1))
}

## ---- Automatic caveats (references/pitfalls.md) ---------------------------

caveats <- character(0)
if (grepl("WARNING", width_note)) {
  caveats <- c(caveats, paste0("**Peak-width bias (pitfalls.md #1)**: ", width_note))
} else {
  caveats <- c(caveats, width_note)
}

if (nrow(res_df) > 0 && n_sig == 0) {
  caveats <- c(caveats, paste0(
    "**No significant differential peaks were detected.** This may reflect a true null result, an ",
    "underpowered design (too few replicates/peaks), or overly strict thresholds -- do not report a ",
    "positive differential accessibility finding from this run."
  ))
} else if (nrow(res_df) >= 20 && n_sig > 0 && n_sig < 3) {
  caveats <- c(caveats, sprintf(
    "**Possible underpowered result**: only %d/%d peaks reached significance. Treat the differential peak list as preliminary and check replicate number/sequencing depth.",
    n_sig, nrow(res_df)
  ))
}
if (nrow(res_df) >= 20 && sig_frac > 0.5) {
  caveats <- c(caveats, sprintf(
    "**Implausibly large differential peak fraction (pitfalls.md #1)**: %d/%d peaks (%.0f%%) are called significant. This can indicate a normalization or majority-shift problem (e.g., global accessibility shift, size-factor distortion) rather than genuine focal differential accessibility -- inspect normalization and diagnostic plots before trusting this list.",
    n_sig, nrow(res_df), 100 * sig_frac
  ))
}
if (!is.null(top_motifs_report) && !motif_bg_matched) {
  caveats <- c(caveats, paste0(
    "**Motif background not confirmed (pitfalls.md #4)**: this table does not state that a GC/length-matched ",
    "background peak set was used. If the motif test was run against a generic genomic background, GC-rich ",
    "motifs may appear enriched purely from base composition -- confirm the background strategy before ",
    "interpreting the ranked list."
  ))
}

## ---- Plain-language conclusion --------------------------------------------

conclusion <- if (nrow(res_df) == 0) {
  "No peaks were available to test; no conclusion can be drawn from this run."
} else if (n_sig == 0) {
  sprintf(
    "Across %d tested peaks, none met the significance thresholds (padj < %s, |log2FC| >= %s), so this run does not support a differential accessibility conclusion for %s vs %s. ",
    nrow(res_df), padj_threshold, log2fc_threshold, case_label, control_label
  )
} else {
  lead <- top_peaks_report[1, ]
  sprintf(
    "Of %d tested peaks, %d (%.0f%%) show significant differential accessibility (%d opening, %d closing) between %s and %s, led by `%s` (log2FC = %.3g, padj = %.3g). ",
    nrow(res_df), n_sig, 100 * sig_frac, n_opening, n_closing, case_label, control_label, lead$peak, lead$log2FoldChange, lead$padj
  )
}
conclusion <- paste0(conclusion, if (!is.null(top_motifs_report)) {
  sprintf(
    "Motif enrichment highlights %s as the top candidate transcription-factor motifs overlapping the tested peak set; this indicates motif accessibility, not confirmed TF binding or activity (pitfalls.md #4, #14), and any nearest-gene or promoter-accessibility interpretation should be reported as a hypothesis, not proven expression change (pitfalls.md #3), unless corroborated by independent data.",
    paste(head(top_motifs_report$label, 3), collapse = ", ")
  )
} else {
  "No motif enrichment results were available to interpret."
})

writeLines(c(
  "# Bulk ATAC differential accessibility report",
  "",
  "## Overview",
  "",
  paste0("- Contrast: ", case_label, " vs ", control_label, " (column `", group_col, "`)"),
  paste0("- Peaks tested: ", nrow(res_df)),
  paste0("- Samples: ", ncol(count_mat)),
  paste0("- Motifs tested: ", motif_n),
  paste0("- Significant differential peaks (padj<", padj_threshold, ", |log2FC|>=", log2fc_threshold, "): ", n_sig,
         " total -- ", n_opening, " opening, ", n_closing, " closing, ",
         nrow(res_df) - n_sig, " not significant"),
  paste0("- Log2FC shrinkage: ", shrink_used),
  "",
  "## Top differential peaks",
  "",
  top_peak_lines,
  "",
  "## Top enriched motifs",
  "",
  top_motif_lines,
  "",
  "## Caveats",
  "",
  paste0("- ", caveats),
  "",
  "## Conclusion",
  "",
  conclusion,
  "",
  "## Methods",
  "",
  paste0("DESeq2 (v", deseq2_version, ", R ", r_version, ") was run on a peak-by-sample count matrix ",
         "treating peaks as count features analogous to genes, using a single-factor design (~", group_col, "). ",
         "Size factors and dispersions used DESeq2 defaults (median-of-ratios normalization; parametric dispersion ",
         "fit with a local fallback). Differential peaks were called from the Wald test in `results()`; ",
         "log2 fold changes reported in `differential_peaks.tsv` were shrunk with `", shrink_used, "` where available ",
         "(p-values/adjusted p-values are from the unshrunk Wald test, per DESeq2 convention). ",
         "Peaks with padj < ", padj_threshold, " and |shrunk log2FC| >= ", log2fc_threshold, " were called opening/closing; all others 'ns'. ",
         "VST (or exact variance-stabilizing transform for small peak sets) was used for PCA, sample correlation, ",
         "and heatmap figures."),
  "",
  "## Limitations",
  "",
  "- Genome build, peak-calling parameters, consensus/fixed-width peak construction, and ENCODE blacklist ",
  "  filtering were assumed to be handled upstream and are not verified by this script.",
  "- Sample QC here (`qc/sample_qc.tsv`) is limited to library depth and detected-peak counts derivable from ",
  "  the count matrix; FRiP, duplicate rate, fragment-size periodicity and TSS enrichment require raw ",
  "  alignment/fragment data and must be supplied/checked separately.",
  "- Peak-to-gene annotation was not run; any gene names in the input are the caller's own annotation and are ",
  "  not validated or re-derived here.",
  "- Motif enrichment (when supplied) is passed through and plotted as-is; this script does not verify that a ",
  "  GC/length-matched background was used upstream (see references/pitfalls.md #4).",
  "- scATAC analysis (Signac/ArchR: TF-IDF/LSI, UMAP, clustering, gene activity, chromVAR) is out of scope for ",
  "  this script; run those workflows interactively in their native object/session."
), file.path(out_dir, "report.md"))

message("Done. Wrote differential_peaks.tsv, normalized_accessibility.tsv, sample_qc.tsv, summary.tsv, report.md, ",
        "and figures (differential_peak_volcano, pca, sample_correlation_heatmap, peak_heatmap",
        if (!is.na(motif_file)) ", motif_enrichment_barplot" else "", ") to ", out_dir)
