#!/usr/bin/env Rscript
# GSEA on a full ranked gene list via fgsea::fgsea() (adaptive multilevel).
#
# Usage:
#   Rscript run_gsea_fgsea.R <ranked_table.tsv> <gene_col> <score_col> <gene_sets.gmt> [output_dir] [padj_threshold]
#
# <ranked_table.tsv> must contain every tested gene with a signed ranking
# statistic (Wald/t statistic or shrunken log2FC) in <score_col> — do not rank
# by raw/unshrunken log2FC (see references/pitfalls.md #2).
#
# [padj_threshold] (arg 6, default 0.05) controls which tested pathways are
# flagged "significant" in the results table, NES dotplot title, and report_gsea.md.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript run_gsea_fgsea.R ranked_table.tsv gene_col score_col gene_sets.gmt [output_dir] [padj_threshold]")
}
ranked_file <- args[[1]]
gene_col <- args[[2]]
score_col <- args[[3]]
gmt_file <- args[[4]]
out_dir <- ifelse(length(args) >= 5, args[[5]], "output")
padj_threshold <- if (length(args) >= 6 && nzchar(args[[6]])) as.numeric(args[[6]]) else 0.05

required <- c("fgsea", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

TEXT_COLOR <- "#173B5F"
set.seed(1)   # reproducibility of fgsea's adaptive multilevel estimation (see pitfalls.md #13)
EPS <- 1e-50

wrap_labels <- function(x, width = 45) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

## Simple Jaccard-overlap heuristic (pitfalls.md #3-style, applied to gene-set
## membership rather than GO terms): flag pairs of significant pathways whose full
## member-gene sets overlap heavily (e.g. overlapping MSigDB/Hallmark collections).
jaccard_redundancy_flags <- function(pathway_names, pathways, jaccard_threshold = 0.5, top_n = 10) {
  names_use <- head(pathway_names, top_n)
  flags <- list()
  n <- length(names_use)
  if (n < 2) return(flags)
  for (i in seq_len(n - 1)) {
    gi <- pathways[[names_use[i]]]
    if (is.null(gi) || length(gi) == 0) next
    for (j in seq(i + 1, n)) {
      gj <- pathways[[names_use[j]]]
      if (is.null(gj) || length(gj) == 0) next
      inter <- length(intersect(gi, gj))
      uni <- length(union(gi, gj))
      jac <- if (uni > 0) inter / uni else 0
      if (jac >= jaccard_threshold) {
        flags[[length(flags) + 1]] <- list(a = names_use[i], b = names_use[j], jac = jac, n = inter)
      }
    }
  }
  if (length(flags) == 0) return(flags)
  flags[order(-vapply(flags, function(f) f$jac, numeric(1)))]
}

## ---- 1. Load ranked list -----------------------------------------------------

ranked <- read.delim(ranked_file, check.names = FALSE)
stopifnot(gene_col %in% colnames(ranked), score_col %in% colnames(ranked))
stats_raw <- ranked[[score_col]]
names(stats_raw) <- ranked[[gene_col]]
n_total <- length(stats_raw)
finite <- is.finite(stats_raw)
dup <- duplicated(names(stats_raw))
stats <- sort(stats_raw[finite & !dup], decreasing = TRUE)
n_dropped_nonfinite <- sum(!finite)
n_dropped_dup <- sum(finite & dup)

## ---- 2. Run fgsea -------------------------------------------------------------

pathways <- fgsea::gmtPathways(gmt_file)
res <- fgsea::fgsea(pathways = pathways, stats = stats, minSize = 10, maxSize = 500, eps = EPS)
res <- as.data.frame(res)
res <- res[order(res$padj), ]

## Long-format leading-edge table, kept before collapsing the list column.
leading_edge_long <- do.call(rbind, lapply(seq_len(nrow(res)), function(i) {
  genes <- res$leadingEdge[[i]]
  if (length(genes) == 0) return(NULL)
  data.frame(pathway = res$pathway[i], gene = genes, NES = res$NES[i], padj = res$padj[i])
}))
if (!is.null(leading_edge_long)) {
  write.table(leading_edge_long, file.path(out_dir, "results", "gsea_leading_edge.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}

res$n_leading_edge <- vapply(res$leadingEdge, length, integer(1))
res$leadingEdge <- vapply(res$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
write.table(res, file.path(out_dir, "results", "fgsea_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## ---- 3. NES summary dotplot (top pathways) ------------------------------------

top <- head(res[order(res$padj), ], 20)
top$pathway_wrapped <- factor(wrap_labels(top$pathway, 45), levels = rev(wrap_labels(top$pathway, 45)))
p_dot <- ggplot2::ggplot(top, ggplot2::aes(NES, pathway_wrapped, color = padj, size = size)) +
  ggplot2::geom_point(alpha = 0.9) +
  ggplot2::scale_color_gradient(low = "#E45B5B", high = "#4EA3F1", trans = "reverse", name = "adjusted p value") +
  ggplot2::scale_size_continuous(name = "Gene set size") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "#A7C6DD") +
  ggplot2::labs(title = paste0("GSEA top pathways (significant: padj < ", padj_threshold, ")"), x = "Normalized enrichment score", y = "Pathway") +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20, color = TEXT_COLOR),
    axis.title = ggplot2::element_text(size = 17, color = TEXT_COLOR),
    axis.text = ggplot2::element_text(size = 12, color = TEXT_COLOR)
  )
ggplot2::ggsave(file.path(out_dir, "figures", "gsea_dotplot.png"), p_dot, width = 9.4, height = max(6, 0.35 * nrow(top) + 2.5), dpi = 320, bg = "white")
ggplot2::ggsave(file.path(out_dir, "figures", "gsea_dotplot.pdf"), p_dot, width = 9.4, height = max(6, 0.35 * nrow(top) + 2.5), bg = "white")

## ---- 4. GSEA running-score plot for the top pathway ---------------------------

res_ranked_by_sig <- res[order(res$padj, -abs(res$NES)), ]
top_pathway <- res_ranked_by_sig$pathway[1]
p_run <- fgsea::plotEnrichment(pathways[[top_pathway]], stats) +
  ggplot2::labs(
    title = paste0("GSEA running score: ", wrap_labels(top_pathway, 60)),
    x = "Rank in ordered gene list", y = "Running enrichment score"
  ) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 18, color = TEXT_COLOR),
    axis.title = ggplot2::element_text(size = 17, color = TEXT_COLOR),
    axis.text = ggplot2::element_text(size = 12, color = TEXT_COLOR)
  )
ggplot2::ggsave(file.path(out_dir, "figures", "gsea_running_score_top.png"), p_run, width = 9, height = 6, dpi = 320, bg = "white")
ggplot2::ggsave(file.path(out_dir, "figures", "gsea_running_score_top.pdf"), p_run, width = 9, height = 6, bg = "white")

## ---- 5. Reports ---------------------------------------------------------------

sig <- res[!is.na(res$padj) & res$padj < padj_threshold, , drop = FALSE]
sig <- sig[order(sig$padj), , drop = FALSE]
n_up <- sum(sig$NES >= 0, na.rm = TRUE)
n_down <- sum(sig$NES < 0, na.rm = TRUE)
top_source <- if (nrow(sig) > 0) sig else res[order(res$padj), , drop = FALSE]
top_label <- if (nrow(sig) > 0) {
  "Top pathways by adjusted p value:"
} else {
  sprintf("No pathway reached padj < %s; top pathways by padj are listed for reference only (not claimed as significant):", padj_threshold)
}
top_n_report <- min(5, nrow(top_source))
top_lines <- if (top_n_report > 0) {
  vapply(seq_len(top_n_report), function(i) {
    r <- top_source[i, ]
    sprintf("  - %s - NES=%.2f, padj=%.3g, leading-edge genes=%d, gene set size=%d", r$pathway, r$NES, r$padj, r$n_leading_edge, r$size)
  }, character(1))
} else character(0)

drop_fraction <- (n_dropped_dup + n_dropped_nonfinite) / max(n_total, 1)

caveats <- character(0)
if (drop_fraction > 0.05) {
  caveats <- c(caveats, sprintf(
    "%.1f%% of input rows (%d of %d) were dropped before ranking (%d duplicate gene ID(s), %d non-finite score(s)); a high drop rate can indicate upstream ID mapping or ranking-metric problems (references/pitfalls.md #4).",
    100 * drop_fraction, n_dropped_dup + n_dropped_nonfinite, n_total, n_dropped_dup, n_dropped_nonfinite
  ))
}
small_edge <- sig[!is.na(sig$n_leading_edge) & sig$n_leading_edge < 3, , drop = FALSE]
if (nrow(small_edge) > 0) {
  caveats <- c(caveats, sprintf(
    "%d significant pathway(s) have fewer than 3 leading-edge genes, which is statistically nominal rather than coordinated pathway regulation (references/pitfalls.md #6): %s.",
    nrow(small_edge), paste0(small_edge$pathway, " (n=", small_edge$n_leading_edge, ")", collapse = ", ")
  ))
}
redundancy_flags <- jaccard_redundancy_flags(sig$pathway, pathways)
if (length(redundancy_flags) > 0) {
  examples <- vapply(head(redundancy_flags, 3), function(f) sprintf("'%s' vs '%s' share %d genes (Jaccard=%.2f)", f$a, f$b, f$n, f$jac), character(1))
  caveats <- c(caveats, sprintf(
    "%d pair(s) of significant pathways share most of their member genes, suggesting overlapping/redundant gene-set definitions rather than independent pathway hits (references/pitfalls.md #3): %s.",
    length(redundancy_flags), paste(examples, collapse = "; ")
  ))
}

conclusion <- if (nrow(sig) == 0) {
  sprintf(
    "None of the %d tested gene sets reached padj < %s for this ranked list; the pathways listed above are the closest to significance but should not be reported as enriched findings without independent support.",
    length(pathways), padj_threshold
  )
} else {
  top1 <- top_source[1, ]
  direction <- if (top1$NES >= 0) "up" else "down"
  sprintf(
    "%d of %d tested gene sets were significant at padj < %s (%d skewed up, %d skewed down), led by '%s' (%s, NES=%.2f, padj=%.3g, %d leading-edge genes). %s",
    nrow(sig), length(pathways), padj_threshold, n_up, n_down, top1$pathway, direction, top1$NES, top1$padj, top1$n_leading_edge,
    if (length(caveats) > 0) "See caveats above before treating every listed pathway as an independent biological finding." else "No redundancy or small-leading-edge issues were flagged among the top pathways."
  )
}

report_lines <- c(
  "# GSEA report", "",
  "## Input", "",
  paste0("- Ranked genes used: ", length(stats), " of ", n_total, " input rows"),
  paste0("- Dropped (non-finite score): ", n_dropped_nonfinite),
  paste0("- Dropped (duplicate gene ID, largest |score| occurrence kept): ", n_dropped_dup),
  paste0("- Gene sets tested (FDR correction scope): ", length(pathways)),
  paste0("- Enrichment significance threshold (which tested pathways are flagged significant here): padj < ", padj_threshold),
  "",
  "## Results", "",
  paste0("- Significant pathways at padj < ", padj_threshold, ": ", nrow(sig), " (", n_up, " up / ", n_down, " down by NES direction)"),
  paste0("- ", top_label),
  top_lines,
  paste0("- Top pathway shown in running-score plot: ", top_pathway),
  paste0("- fgsea eps: ", EPS, "; seed: 1 (see references/pitfalls.md #13 on permutation resolution)"),
  if (length(caveats) > 0) c("", "## Caveats", "", paste0("- ", caveats)) else character(0),
  "",
  "## Conclusion", "",
  conclusion
)
writeLines(report_lines, file.path(out_dir, "report_gsea.md"))
