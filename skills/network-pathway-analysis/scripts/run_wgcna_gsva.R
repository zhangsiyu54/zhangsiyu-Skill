#!/usr/bin/env Rscript
# Formal WGCNA (co-expression network / module-trait / hub gene) + GSVA/ssGSEA
# (pathway activity) pipeline. See references/cns-method-guide.md,
# references/pitfalls.md, and references/deliverables.md for the required
# figures, statistics, and output files this script must produce.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(
    "Usage: Rscript run_wgcna_gsva.R expression_matrix.csv metadata.csv trait_cols gene_sets.gmt output_dir [scoring_method] [min_samples] [padj_threshold] [min_module_size]\n",
    "  trait_cols: comma-separated phenotype column name(s) in metadata.csv, e.g. \"group,severity\"\n",
    "  scoring_method: \"gsva\" (default) or \"ssgsea\"\n",
    "  min_samples: below this sample count, run exploratory correlation mode instead of full WGCNA module detection (default 15, see pitfalls.md #2)\n",
    "  padj_threshold: BH-FDR adjusted-p significance cutoff for module-trait/pathway-group significance stars and report.md text (default 0.05)\n",
    "  min_module_size: minimum WGCNA module size passed to blockwiseModules() (default 30, see pitfalls.md #8; ignored in exploratory mode)"
  )
}
expr_file <- args[[1]]
metadata_file <- args[[2]]
trait_cols <- strsplit(args[[3]], ",")[[1]]
trait_cols <- trimws(trait_cols)
gmt_file <- args[[4]]
out_dir <- args[[5]]
scoring_method <- ifelse(length(args) >= 6, args[[6]], "gsva")
min_samples <- ifelse(length(args) >= 7, as.integer(args[[7]]), 15L)
padj_threshold <- ifelse(length(args) >= 8, as.numeric(args[[8]]), 0.05)
min_module_size <- ifelse(length(args) >= 9, as.integer(args[[9]]), 30L)

if (!scoring_method %in% c("gsva", "ssgsea")) stop("scoring_method must be 'gsva' or 'ssgsea'")
if (is.na(padj_threshold) || padj_threshold <= 0 || padj_threshold >= 1) stop("padj_threshold must be a number between 0 and 1")
if (is.na(min_module_size) || min_module_size < 2) stop("min_module_size must be an integer >= 2")

required <- c("WGCNA", "GSVA", "GSEABase", "ggplot2", "pheatmap")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

## WGCNA must be attached (not just namespace-loaded) because blockwiseModules() and
## related internals resolve the correlation function by name (match.fun("cor")) against
## the search path at call time; without library(WGCNA) that resolves to base stats::cor,
## which lacks the extra arguments (weights.x/weights.y/cosine) WGCNA's internals pass,
## crashing with "unused arguments" deep inside blockwiseModules. All other packages here
## are only ever called via `::` and do not need to be attached.
suppressPackageStartupMessages(library(WGCNA))

dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "run"), recursive = TRUE, showWarnings = FALSE)

## ---- Fixed analysis thresholds (reported in methods for transparency) ----
## MIN_MODULE_SIZE is CLI-configurable (min_module_size arg, default 30 above); the
## rest remain fixed constants since they are conventional defaults (pitfalls.md #8, #12)
## rather than knobs users routinely need to change per run.
MIN_MODULE_SIZE <- min_module_size
MERGE_CUT_HEIGHT <- 0.25       # pitfalls.md #8
NETWORK_TYPE <- "signed"       # pitfalls.md #12
KME_THRESHOLD <- 0.8           # pitfalls.md #4
GS_THRESHOLD <- 0.2            # pitfalls.md #4
MIN_GENE_SET_SIZE <- 10        # pitfalls.md #7
MAX_GENE_SET_SIZE <- 500       # pitfalls.md #7
TOP_EXPLORATORY_GENES <- 500   # cap for exploratory small-n mode

save_gg <- function(plot, name, width = 8, height = 6) {
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(name, ".png")), plot, width = width, height = height, dpi = 320, bg = "white")
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(name, ".pdf")), plot, width = width, height = height, bg = "white")
}

save_base_plot <- function(draw_fn, name, width = 9, height = 6, res = 320) {
  png(file.path(out_dir, "figures", paste0(name, ".png")), width = width, height = height, units = "in", res = res)
  draw_fn()
  dev.off()
  pdf(file.path(out_dir, "figures", paste0(name, ".pdf")), width = width, height = height)
  draw_fn()
  dev.off()
}

## pheatmap's `main` title is centered on the plot device at a fixed font size and does
## not wrap or shrink to fit; a long title (e.g. the exploratory-mode heatmap caption)
## silently runs off both edges of a narrow device instead of erroring. Widen the device
## to fit the longest title line rather than truncating/wrapping the caption text.
title_safe_width <- function(title, base_width) {
  lines <- strsplit(title, "\n", fixed = TRUE)[[1]]
  needed <- max(nchar(lines)) * 0.13 + 1
  max(base_width, needed)
}

corr_palette <- colorRampPalette(c("#4EA3F1", "white", "#E45B5B"))(101)

## ---- Load and align inputs ----
expr <- read.csv(expr_file, check.names = FALSE)
metadata <- read.csv(metadata_file, check.names = FALSE)
stopifnot("gene" %in% colnames(expr), "sample" %in% colnames(metadata))
missing_traits <- setdiff(trait_cols, colnames(metadata))
if (length(missing_traits)) stop("Trait column(s) not found in metadata: ", paste(missing_traits, collapse = ", "))
rownames(expr) <- expr$gene
expr$gene <- NULL
metadata <- metadata[match(colnames(expr), metadata$sample), ]
if (any(is.na(metadata$sample))) stop("Sample IDs do not match metadata.")

datExpr <- t(as.matrix(expr))          # samples (rows) x genes (columns), as required by WGCNA
storage.mode(datExpr) <- "numeric"

## ---- QC: goodSamplesGenes (remove zero-variance/low-detection features, pitfalls.md #13 + cns-method-guide "Matrix QC") ----
gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  message("goodSamplesGenes flagged ", sum(!gsg$goodGenes), " genes and ", sum(!gsg$goodSamples), " samples for removal.")
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  metadata <- metadata[gsg$goodSamples, , drop = FALSE]
}
nSamples <- nrow(datExpr)
nGenes <- ncol(datExpr)

## ---- Required figure: sample clustering QC (pitfalls.md #13) ----
save_base_plot(function() {
  plot(hclust(dist(datExpr), method = "average"), main = "Sample clustering (Euclidean, average linkage)",
       xlab = "", sub = "", cex = 0.8)
}, "sample_clustering", width = 8.4, height = 5.4)

## ---- Trait table (numeric-encoded for correlation) ----
traitData <- metadata[, trait_cols, drop = FALSE]
for (col in trait_cols) {
  if (!is.numeric(traitData[[col]])) traitData[[col]] <- as.numeric(factor(traitData[[col]]))
}

exploratory <- nSamples < min_samples
softPower <- NA_integer_
achieved_R2 <- NA_real_
net <- NULL

if (exploratory) {
  ## pitfalls.md #2 / cns-method-guide decision tree row 1: below min_samples, do not run
  ## formal WGCNA module detection. Fall back to a correlation-based screen of the most
  ## variable genes against traits and label everything exploratory.
  message("Sample size (", nSamples, ") is below min_samples (", min_samples, "); running EXPLORATORY mode (no WGCNA module detection).")
  gene_var <- apply(datExpr, 2, stats::var, na.rm = TRUE)
  keep_n <- min(TOP_EXPLORATORY_GENES, nGenes)
  top_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(keep_n)]
  MEs <- datExpr[, top_genes, drop = FALSE]
  colnames(MEs) <- top_genes
  moduleColors <- setNames(rep("exploratory", nGenes), colnames(datExpr))
} else {
  ## ---- Soft-thresholding power selection with explicit scale-free-fit check (pitfalls.md #1) ----
  powers <- c(seq(1, 10, 1), seq(12, 20, 2))
  sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers, networkType = NETWORK_TYPE, verbose = 0)
  fit <- sft$fitIndices
  fit$signedR2 <- -sign(fit$slope) * fit$SFT.R.sq
  candidates <- fit$Power[fit$signedR2 >= 0.80]
  softPower <- if (length(candidates)) min(candidates) else fit$Power[which.max(fit$signedR2)]
  achieved_R2 <- fit$signedR2[fit$Power == softPower][1]
  if (is.na(sft$powerEstimate) || achieved_R2 < 0.80) {
    message("Scale-free fit never reached R^2 >= 0.80 (achieved ", round(achieved_R2, 3), " at power ", softPower, "); treat modules as more exploratory (pitfalls.md #1).")
  }

  save_base_plot(function() {
    par(mfrow = c(1, 2))
    plot(fit$Power, fit$signedR2, xlab = "Soft threshold (power)", ylab = "Scale-free topology fit (signed R^2)",
         type = "b", pch = 19, col = "#4EA3F1", main = "Scale-free topology fit")
    abline(h = 0.80, col = "#E45B5B", lty = 2)
    text(fit$Power, fit$signedR2, labels = fit$Power, pos = 3, cex = 0.7)
    plot(fit$Power, fit$mean.k., xlab = "Soft threshold (power)", ylab = "Mean connectivity",
         type = "b", pch = 19, col = "#4CB984", main = "Mean connectivity")
    text(fit$Power, fit$mean.k., labels = fit$Power, pos = 3, cex = 0.7)
  }, "soft_threshold_selection", width = 10, height = 5)

  ## ---- Network construction and module detection ----
  net <- WGCNA::blockwiseModules(
    datExpr,
    power = softPower,
    networkType = NETWORK_TYPE,
    TOMType = NETWORK_TYPE,
    minModuleSize = MIN_MODULE_SIZE,
    mergeCutHeight = MERGE_CUT_HEIGHT,
    numericLabels = TRUE,
    saveTOMs = FALSE,
    verbose = 0
  )
  moduleColors <- WGCNA::labels2colors(net$colors)
  names(moduleColors) <- colnames(datExpr)

  ## ---- Required figure: module dendrogram ----
  save_base_plot(function() {
    WGCNA::plotDendroAndColors(
      net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
      groupLabels = "Module", dendroLabels = FALSE, addGuide = TRUE,
      hang = 0.03, guideHang = 0.05, main = "Gene clustering and module assignment"
    )
  }, "module_dendrogram", width = 9, height = 6)

  MEs <- WGCNA::moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
  MEs <- WGCNA::orderMEs(MEs)
}

## ---- Module (or exploratory gene)-trait correlation matrix with full-matrix FDR (pitfalls.md #3) ----
moduleTraitCor <- stats::cor(MEs, traitData, use = "pairwise.complete.obs")
moduleTraitP <- WGCNA::corPvalueStudent(moduleTraitCor, nSamples)
module_trait <- data.frame(
  module = rep(rownames(moduleTraitCor), times = ncol(moduleTraitCor)),
  trait = rep(colnames(moduleTraitCor), each = nrow(moduleTraitCor)),
  correlation = as.vector(moduleTraitCor),
  pvalue = as.vector(moduleTraitP),
  stringsAsFactors = FALSE
)
module_trait$padj <- stats::p.adjust(module_trait$pvalue, method = "BH")
module_trait$n_samples <- nSamples
module_trait$mode <- if (exploratory) "exploratory_top_variable_gene" else "wgcna_module"
write.table(module_trait, file.path(out_dir, "results", "module_trait_correlation.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## ---- Required figure: module-trait heatmap (blue-white-red, adjusted-p stars only, pitfalls.md #3) ----
## Star tiers scale off padj_threshold (default 0.05) so a custom cutoff still nests
## ***/**/* consistently: loosest tier is the threshold itself, tighter tiers at
## threshold/5 and threshold/50 (matches the conventional 0.05/0.01/0.001 spacing).
padj_mat <- matrix(stats::p.adjust(as.vector(moduleTraitP), method = "BH"), nrow = nrow(moduleTraitCor), dimnames = dimnames(moduleTraitCor))
label_mat <- matrix(
  sprintf("%.2f%s", moduleTraitCor, ifelse(padj_mat < padj_threshold / 50, "***", ifelse(padj_mat < padj_threshold / 5, "**", ifelse(padj_mat < padj_threshold, "*", "")))),
  nrow = nrow(moduleTraitCor), dimnames = dimnames(moduleTraitCor)
)

## Cap heatmap rows for readability (exploratory mode can screen hundreds of individual
## genes, which would otherwise render an unusably tall figure); the full, uncapped
## matrix is still written to results/module_trait_correlation.tsv above. Rank by the
## strongest |correlation| to any trait.
MAX_HEATMAP_ROWS <- 40
plot_moduleTraitCor <- moduleTraitCor
plot_label_mat <- label_mat
heatmap_capped <- nrow(moduleTraitCor) > MAX_HEATMAP_ROWS
if (heatmap_capped) {
  row_rank <- order(-apply(abs(moduleTraitCor), 1, max, na.rm = TRUE))[seq_len(MAX_HEATMAP_ROWS)]
  plot_moduleTraitCor <- moduleTraitCor[row_rank, , drop = FALSE]
  plot_label_mat <- label_mat[row_rank, , drop = FALSE]
}
heatmap_title <- if (exploratory) {
  paste0("Top-variable-gene x trait correlation (exploratory, n < min_samples)", if (heatmap_capped) paste0("\n(top ", MAX_HEATMAP_ROWS, " of ", nrow(moduleTraitCor), " genes by |r|; full table in results/)") else "")
} else "Module-trait correlation"
## Row clustering needs >= 2 rows (hclust fails on a single module/gene); guard rather
## than crash when only one module survives (e.g. tiny cohorts, weak module structure).
can_cluster_rows <- !exploratory && nrow(plot_moduleTraitCor) > 1
mt_width <- title_safe_width(heatmap_title, 7)
pheatmap::pheatmap(plot_moduleTraitCor, color = corr_palette, breaks = seq(-1, 1, length.out = 102),
                    display_numbers = plot_label_mat, number_color = "#173B5F", fontsize_number = 10,
                    cluster_rows = can_cluster_rows, cluster_cols = FALSE,
                    main = heatmap_title,
                    filename = file.path(out_dir, "figures", "module_trait_heatmap.png"), width = mt_width, height = max(5, 0.35 * nrow(plot_moduleTraitCor) + 3))
pheatmap::pheatmap(plot_moduleTraitCor, color = corr_palette, breaks = seq(-1, 1, length.out = 102),
                    display_numbers = plot_label_mat, number_color = "#173B5F", fontsize_number = 10,
                    cluster_rows = can_cluster_rows, cluster_cols = FALSE,
                    main = heatmap_title,
                    filename = file.path(out_dir, "figures", "module_trait_heatmap.pdf"), width = mt_width, height = max(5, 0.35 * nrow(plot_moduleTraitCor) + 3))

## ---- Module assignment table ----
write.table(data.frame(gene = colnames(datExpr), module = moduleColors[colnames(datExpr)]),
            file.path(out_dir, "results", "gene_modules.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## ---- Hub genes: intramodular connectivity (kME) + gene significance (GS) (pitfalls.md #4) ----
primary_trait <- traitData[[trait_cols[1]]]
GS <- as.numeric(stats::cor(datExpr, primary_trait, use = "pairwise.complete.obs"))
GS_p <- as.numeric(WGCNA::corPvalueStudent(GS, nSamples))

if (exploratory) {
  hub_table <- data.frame(
    gene = colnames(datExpr),
    module = "exploratory",
    kME_own_module = NA_real_,
    GS = GS,
    GS_pvalue = GS_p,
    is_hub = FALSE,
    note = "exploratory mode (n < min_samples): ranked by |correlation to primary trait|, not a WGCNA hub gene",
    stringsAsFactors = FALSE
  )
  hub_table <- hub_table[order(-abs(hub_table$GS)), ]
  hub_table <- head(hub_table, TOP_EXPLORATORY_GENES)
} else {
  kME_table <- WGCNA::signedKME(datExpr, MEs, outputColumnName = "kME")
  kME_own <- rep(NA_real_, nGenes)
  for (i in seq_len(nGenes)) {
    col <- paste0("kME", moduleColors[i])
    if (col %in% colnames(kME_table)) kME_own[i] <- kME_table[i, col]
  }
  hub_table <- data.frame(
    gene = colnames(datExpr),
    module = moduleColors[colnames(datExpr)],
    kME_own_module = kME_own,
    GS = GS,
    GS_pvalue = GS_p,
    stringsAsFactors = FALSE
  )
  hub_table <- hub_table[hub_table$module != "grey", ]  # grey = unassigned genes, not a real module
  hub_table$is_hub <- !is.na(hub_table$kME_own_module) & hub_table$kME_own_module >= KME_THRESHOLD & abs(hub_table$GS) >= GS_THRESHOLD
  # rep(..., length.out=nrow()) instead of a bare scalar assignment: when every gene lands in the
  # "grey" (unassigned) module -- e.g. weak/no real co-expression structure -- hub_table has 0 rows,
  # and `df$col <- scalar` errors ("replacement has 1 row, data has 0") on a 0-row data.frame.
  hub_table$note <- rep(sprintf("hub thresholds: kME >= %.2f and |GS| >= %.2f (topologically central candidates, not causal drivers; pitfalls.md #4)", KME_THRESHOLD, GS_THRESHOLD), length.out = nrow(hub_table))
  hub_table <- hub_table[order(hub_table$module, -hub_table$kME_own_module), ]
}
write.table(hub_table, file.path(out_dir, "results", "hub_genes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## ---- Required figure: hub gene / network summary ----
if (exploratory) {
  top_hub_plot <- head(hub_table, 15)
  top_hub_plot$gene <- factor(top_hub_plot$gene, levels = rev(top_hub_plot$gene))
  p_hub <- ggplot2::ggplot(top_hub_plot, ggplot2::aes(gene, GS)) +
    ggplot2::geom_col(fill = "#8B72D6", alpha = 0.85) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Top trait-correlated genes (exploratory)", x = NULL, y = "Correlation to primary trait") +
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20, color = "#173B5F"),
                   axis.title = ggplot2::element_text(size = 17, color = "#173B5F"))
} else if (nrow(hub_table) == 0) {
  # Every gene landed in the unassigned "grey" module -- no real co-expression modules were
  # detected at this soft power/cut height (weak signal, tiny toy dataset, or genuinely no
  # structure). Skip the hub-gene plot with a clear placeholder rather than crashing on an
  # empty top_per_module (split() on 0 rows -> rbind(list()) -> NULL downstream).
  message("No non-grey modules detected: skipping hub-gene network summary figure (see report caveats).")
  top_per_module <- hub_table[0, ]
  p_hub <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = "No non-grey modules detected\n(all genes assigned to the unassigned 'grey' module)", size = 5, color = "#173B5F") +
    ggplot2::theme_void()
} else {
  top_per_module <- do.call(rbind, lapply(split(hub_table, hub_table$module), function(d) d[which.max(d$kME_own_module), ]))
  top_per_module <- top_per_module[order(-top_per_module$kME_own_module), ]
  top_per_module$gene <- factor(top_per_module$gene, levels = rev(top_per_module$gene))
  p_hub <- ggplot2::ggplot(top_per_module, ggplot2::aes(gene, kME_own_module, fill = module)) +
    ggplot2::geom_col(alpha = 0.9, show.legend = FALSE) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Top hub gene per module (intramodular connectivity)", x = NULL, y = "kME (module membership)") +
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20, color = "#173B5F"),
                   axis.title = ggplot2::element_text(size = 17, color = "#173B5F"))
}
save_gg(p_hub, "hub_gene_network_summary", width = 8, height = max(4.5, 0.35 * nrow(if (exploratory) head(hub_table, 15) else top_per_module) + 2))

## ---- Pathway activity scoring: GSVA/ssGSEA on the same expression matrix (genes x samples) ----
## Filter gene sets by size after matching to the matrix (pitfalls.md #7) and require
## normalized/log expression, not raw counts (pitfalls.md #6) -- caller is responsible
## for supplying normalized/log expression as documented in references/input-contract.md.
exprMatrix <- as.matrix(expr)
storage.mode(exprMatrix) <- "numeric"
gene_sets <- GSEABase::getGmt(gmt_file)
set_sizes <- vapply(GSEABase::geneIds(gene_sets), function(g) length(intersect(g, rownames(exprMatrix))), integer(1))
keep <- set_sizes >= MIN_GENE_SET_SIZE & set_sizes <= MAX_GENE_SET_SIZE
n_dropped <- sum(!keep)
gene_sets <- gene_sets[keep]
if (length(gene_sets) == 0) stop("No gene sets remain after size filtering (min=", MIN_GENE_SET_SIZE, ", max=", MAX_GENE_SET_SIZE, "); check gene ID matching.")
if (n_dropped > 0) message(n_dropped, " gene set(s) dropped for having <", MIN_GENE_SET_SIZE, " or >", MAX_GENE_SET_SIZE, " matched genes.")

score_pathways <- function(exprMatrix, gene_sets, method) {
  ## GSVA >= 1.50 (Bioc 3.18+) uses parameter-object constructors; older releases used
  ## gsva(expr, gset.idx.list, method=...). Try the current API first, fall back for
  ## older installs so this script works across GSVA versions without guessing internals.
  result <- tryCatch({
    param <- if (method == "ssgsea") GSVA::ssgseaParam(exprMatrix, gene_sets) else GSVA::gsvaParam(exprMatrix, gene_sets, kcdf = "Gaussian")
    GSVA::gsva(param, verbose = FALSE)
  }, error = function(e) NULL)
  if (is.null(result)) {
    result <- GSVA::gsva(exprMatrix, gene_sets, method = method, kcdf = "Gaussian", verbose = FALSE)
  }
  result
}
gsva_scores <- score_pathways(exprMatrix, gene_sets, scoring_method)
pathway_out <- as.data.frame(gsva_scores)
pathway_out <- cbind(pathway = rownames(pathway_out), pathway_out)
write.table(pathway_out, file.path(out_dir, "results", "pathway_scores.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## ---- Required figure: pathway activity heatmap ----
score_mat <- as.matrix(gsva_scores)
lim <- max(abs(score_mat), na.rm = TRUE)
if (!is.finite(lim) || lim <= 0) lim <- 1e-6
pathway_heatmap_title <- paste0("Pathway activity (", scoring_method, ")")
pathway_heatmap_width <- title_safe_width(pathway_heatmap_title, max(7, 0.4 * ncol(score_mat) + 3))
pheatmap::pheatmap(score_mat, color = corr_palette, breaks = seq(-lim, lim, length.out = 102),
                    main = pathway_heatmap_title,
                    cluster_rows = nrow(score_mat) > 1, cluster_cols = ncol(score_mat) > 1,
                    fontsize_row = 10, fontsize_col = 10,
                    filename = file.path(out_dir, "figures", "pathway_activity_heatmap.png"), width = pathway_heatmap_width, height = max(5, 0.35 * nrow(score_mat) + 2.5))
pheatmap::pheatmap(score_mat, color = corr_palette, breaks = seq(-lim, lim, length.out = 102),
                    main = pathway_heatmap_title,
                    cluster_rows = nrow(score_mat) > 1, cluster_cols = ncol(score_mat) > 1,
                    fontsize_row = 10, fontsize_col = 10,
                    filename = file.path(out_dir, "figures", "pathway_activity_heatmap.pdf"), width = pathway_heatmap_width, height = max(5, 0.35 * nrow(score_mat) + 2.5))

## ---- Pathway boxplot + group summary when the primary trait is categorical (Required Figures: pathway boxplot) ----
primary_raw <- metadata[[trait_cols[1]]]
is_categorical <- !is.numeric(primary_raw) || length(unique(primary_raw)) <= 6
group_summary <- NULL
if (is_categorical) {
  groups <- as.character(primary_raw)
  long <- do.call(rbind, lapply(rownames(score_mat), function(pw) {
    data.frame(pathway = pw, sample = colnames(score_mat), score = as.numeric(score_mat[pw, ]), group = groups, stringsAsFactors = FALSE)
  }))
  p_box <- ggplot2::ggplot(long, ggplot2::aes(group, score, fill = group)) +
    ggplot2::geom_boxplot(outlier.size = 0.8, alpha = 0.85) +
    ggplot2::facet_wrap(~pathway, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c("#4EA3F1", "#E45B5B", "#4CB984", "#E2A13B", "#8B72D6")) +
    ggplot2::labs(title = paste0("Pathway activity by ", trait_cols[1]), x = NULL, y = "Pathway score") +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20, color = "#173B5F"),
                   axis.title = ggplot2::element_text(size = 17, color = "#173B5F"),
                   legend.position = "top")
  save_gg(p_box, "pathway_activity_boxplot", width = max(8, 3 * length(unique(long$pathway)) ^ 0.5 + 4), height = max(6, 3 * ceiling(length(unique(long$pathway)) / 3)))

  group_levels <- unique(groups)
  summary_rows <- lapply(rownames(score_mat), function(pw) {
    vals <- score_mat[pw, ]
    means <- tapply(vals, groups, mean, na.rm = TRUE)
    if (length(group_levels) == 2) {
      test <- "wilcoxon"
      fit <- tryCatch(stats::wilcox.test(vals ~ groups), error = function(e) NULL)
    } else {
      test <- "kruskal_wallis"
      fit <- tryCatch(stats::kruskal.test(vals ~ groups), error = function(e) NULL)
    }
    stat <- if (is.null(fit)) NA_real_ else unname(fit$statistic)
    pv <- if (is.null(fit)) NA_real_ else fit$p.value
    c(list(pathway = pw), as.list(means), list(test = test, statistic = stat, pvalue = pv))
  })
  group_summary <- do.call(rbind, lapply(summary_rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  group_summary$padj <- stats::p.adjust(group_summary$pvalue, method = "BH")
  write.table(group_summary, file.path(out_dir, "results", "pathway_group_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  message("Primary trait '", trait_cols[1], "' is continuous; skipping pathway boxplot/group summary (no discrete groups).")
}

## ---- Methods / report text ----
pkg_versions <- vapply(c("WGCNA", "GSVA", "GSEABase", "ggplot2", "pheatmap"), function(p) as.character(utils::packageVersion(p)), character(1))

n_modules_detected <- if (exploratory) NA_integer_ else length(setdiff(unique(moduleColors), "grey"))

## Top module/gene-trait associations, ranked by |correlation| (pitfalls.md #10: association
## language only, no causal wording).
mt_ranked <- module_trait[order(-abs(module_trait$correlation)), , drop = FALSE]
mt_top <- head(mt_ranked, min(3, nrow(mt_ranked)))
fdr_tested <- sum(!is.na(module_trait$padj))
fdr_survived <- sum(!is.na(module_trait$padj) & module_trait$padj < padj_threshold)

## Top pathway group differences, ranked by padj (skip rows with no test result).
pw_top <- NULL
if (!is.null(group_summary)) {
  pw_ranked <- group_summary[!is.na(group_summary$padj), , drop = FALSE]
  pw_ranked <- pw_ranked[order(pw_ranked$padj), , drop = FALSE]
  pw_top <- head(pw_ranked, min(3, nrow(pw_ranked)))
}

## Hub gene names for the top-associated module (formal WGCNA mode only; exploratory mode
## has no discrete modules, so report the top individually-correlated gene instead).
top_hub_genes <- character(0)
top_module_label <- mt_top$module[1]
if (!exploratory && nrow(mt_top) > 0) {
  top_module_color <- sub("^ME", "", top_module_label)
  module_hubs <- hub_table[hub_table$module == top_module_color & hub_table$is_hub, , drop = FALSE]
  if (nrow(module_hubs) == 0) module_hubs <- hub_table[hub_table$module == top_module_color, , drop = FALSE]
  module_hubs <- module_hubs[order(-module_hubs$kME_own_module), , drop = FALSE]
  top_hub_genes <- head(module_hubs$gene, 5)
}

fmt_p <- function(p) if (is.na(p)) "NA" else trimws(formatC(p, format = "g", digits = 3))
stars <- function(p) if (is.na(p)) "" else if (p < padj_threshold / 50) "***" else if (p < padj_threshold / 5) "**" else if (p < padj_threshold) "*" else ""

mt_lines <- if (nrow(mt_top) > 0) {
  vapply(seq_len(nrow(mt_top)), function(i) {
    r <- mt_top[i, ]
    sprintf("- **%s** vs **%s**: r=%.2f, p=%s, padj=%s%s", r$module, r$trait, r$correlation, fmt_p(r$pvalue), fmt_p(r$padj), stars(r$padj))
  }, character(1))
} else "- No module/gene-trait correlations available."

pw_lines <- if (!is.null(pw_top) && nrow(pw_top) > 0) {
  c(paste0("Group comparison across: ", paste(group_levels, collapse = ", "), "."),
    vapply(seq_len(nrow(pw_top)), function(i) {
      r <- pw_top[i, ]
      sprintf("- **%s**: %s statistic=%s, p=%s, padj=%s%s", r$pathway, r$test, fmt_p(r$statistic), fmt_p(r$pvalue), fmt_p(r$padj), stars(r$padj))
    }, character(1)))
} else if (is_categorical) {
  "Group summary computed but no pathway had a valid test result."
} else {
  paste0("Primary trait '", trait_cols[1], "' is continuous; no discrete-group pathway comparison was run.")
}

hub_lines <- if (length(top_hub_genes) > 0) {
  c("", "## Hub genes", sprintf("Top module **%s** hub gene candidate(s) (kME/GS-based, topological centrality only): %s.", top_module_label, paste(top_hub_genes, collapse = ", ")))
} else if (exploratory && nrow(mt_top) > 0) {
  c("", "## Hub genes", sprintf("Exploratory mode has no discrete modules; top individually trait-correlated gene: **%s** (r=%.2f vs %s).", mt_top$module[1], mt_top$correlation[1], mt_top$trait[1]))
} else character(0)

## Caveats: only print what applies (references/pitfalls.md).
caveats <- character(0)
if (nSamples < 15) {
  caveats <- c(caveats, sprintf("Sample size (n=%d) is well below the ~15-20 samples recommended for stable WGCNA module detection; treat all module-trait associations as strongly exploratory (pitfalls.md #2).", nSamples))
} else if (nSamples < 20) {
  caveats <- c(caveats, sprintf("Sample size (n=%d) is below the conventional 15-20 threshold for stable WGCNA modules; module composition and hub genes may not be reproducible (pitfalls.md #2).", nSamples))
}
if (fdr_tested > 0) {
  caveats <- c(caveats, sprintf("BH-FDR multiple-testing correction was applied across the full module x trait matrix (%d cell(s) tested); %d/%d remained significant at padj<%g (pitfalls.md #3).", fdr_tested, fdr_survived, fdr_tested, padj_threshold))
} else {
  caveats <- c(caveats, "No module-trait p-values could be computed; correlations are descriptive only.")
}
if (!is.null(pw_top) || !is.null(group_summary)) {
  pw_fdr_tested <- if (!is.null(group_summary)) sum(!is.na(group_summary$padj)) else 0
  pw_fdr_survived <- if (!is.null(group_summary)) sum(!is.na(group_summary$padj) & group_summary$padj < padj_threshold) else 0
  if (pw_fdr_tested > 0) {
    caveats <- c(caveats, sprintf("Pathway group-comparison p-values were BH-FDR corrected across all scored pathways (%d tested); %d/%d remained significant at padj<%g.", pw_fdr_tested, pw_fdr_survived, pw_fdr_tested, padj_threshold))
  }
}
caveats <- c(caveats, sprintf("Pathway scores (%s) are only valid within this single joint scoring run; do not directly compare these values to scores computed in a separate GSVA/ssGSEA invocation on a different sample set (pitfalls.md #5).", scoring_method))
if (heatmap_capped) {
  caveats <- c(caveats, sprintf("The module-trait heatmap figure shows only the top %d of %d rows ranked by |correlation| for readability; the full matrix is in results/module_trait_correlation.tsv.", MAX_HEATMAP_ROWS, nrow(moduleTraitCor)))
}
if (length(top_hub_genes) > 0) {
  caveats <- c(caveats, "Hub gene status reflects intramodular connectivity (kME) and gene significance (GS) only -- descriptive network centrality, not evidence of a causal or regulatory role (pitfalls.md #4).")
}

## Plain-language conclusion.
conclusion_parts <- character(0)
if (nrow(mt_top) > 0) {
  r <- mt_top[1, ]
  direction <- if (r$correlation > 0) "positively" else "negatively"
  conclusion_parts <- c(conclusion_parts, sprintf("%s %s is most strongly associated with %s (r=%.2f, %s correlated, padj=%s)",
                                                    if (exploratory) "Gene" else "Module", r$module, r$trait, r$correlation, direction, fmt_p(r$padj)))
}
if (!is.null(pw_top) && nrow(pw_top) > 0) {
  r <- pw_top[1, ]
  conclusion_parts <- c(conclusion_parts, sprintf("pathway activity differs most between groups for %s (%s, p=%s)", r$pathway, r$test, fmt_p(r$pvalue)))
}
if (length(top_hub_genes) > 0) {
  conclusion_parts <- c(conclusion_parts, sprintf("candidate hub gene(s) for that module include %s", paste(head(top_hub_genes, 3), collapse = ", ")))
}
conclusion <- if (length(conclusion_parts) > 0) {
  paste0(paste(conclusion_parts, collapse = "; "), ".")
} else {
  "No module-trait or pathway-group signal could be summarized from this run."
}
conclusion <- paste(conclusion,
  "These results are correlational and exploratory: they describe co-expression and pathway-activity patterns in this dataset, not causal or mechanistic relationships. Validate hub genes and pathway differences in an independent cohort before drawing biological conclusions (pitfalls.md #4, #10).")

report_lines <- c(
  "# WGCNA / GSVA report",
  "",
  paste0("Samples analyzed: ", nSamples, if (exploratory) " (EXPLORATORY MODE: below min_samples threshold, see pitfalls.md #2)" else ""),
  paste0("Genes analyzed: ", nGenes),
  paste0("Network type: ", NETWORK_TYPE, "; min module size: ", MIN_MODULE_SIZE, "; merge cut height: ", MERGE_CUT_HEIGHT),
  if (!exploratory) paste0("Soft power: ", softPower, " (achieved scale-free fit R^2 = ", round(achieved_R2, 3), ")") else "Soft power: not applicable (exploratory mode)",
  if (!exploratory) paste0("Modules detected (excl. grey): ", n_modules_detected) else paste0("Top variable genes screened: ", ncol(MEs)),
  paste0("Module/gene-trait correlation: Pearson; traits: ", paste(trait_cols, collapse = ", "), "; multiple-testing correction: BH-FDR across full matrix, significance threshold padj<", padj_threshold, " (pitfalls.md #3)"),
  paste0("Hub gene criteria: kME >= ", KME_THRESHOLD, " and |GS| >= ", GS_THRESHOLD, " (descriptive network centrality, not causal; pitfalls.md #4)"),
  paste0("Pathway scoring: ", scoring_method, " on normalized/log expression; gene sets filtered to size [", MIN_GENE_SET_SIZE, ", ", MAX_GENE_SET_SIZE, "] (", n_dropped, " dropped); scores are only comparable within this joint run (pitfalls.md #5)"),
  paste0("Gene set file: ", basename(gmt_file)),
  "Package versions:",
  paste0("  ", names(pkg_versions), " ", pkg_versions),
  "",
  "## Module-trait results",
  mt_lines,
  "",
  "## Pathway activity results",
  pw_lines,
  hub_lines,
  "",
  "## Caveats",
  paste0("- ", caveats),
  "",
  "## Conclusion",
  conclusion
)
writeLines(report_lines, file.path(out_dir, "report.md"))

writeLines(
  c(paste0("item\tvalue"),
    paste0("samples\t", nSamples),
    paste0("genes\t", nGenes),
    paste0("exploratory_mode\t", exploratory),
    paste0("soft_power\t", ifelse(exploratory, NA, softPower)),
    paste0("modules\t", ifelse(exploratory, NA, length(setdiff(unique(moduleColors), "grey")))),
    paste0("pathways\t", nrow(score_mat)),
    paste0("padj_threshold\t", padj_threshold),
    paste0("min_module_size\t", MIN_MODULE_SIZE)),
  file.path(out_dir, "summary.tsv")
)
