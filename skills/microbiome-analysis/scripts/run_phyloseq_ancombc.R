#!/usr/bin/env Rscript
# Formal phyloseq / ANCOM-BC2 microbiome workflow.
# Usage: Rscript run_phyloseq_ancombc.R abundance.csv taxonomy.csv metadata.csv group_col output_dir [tax_level] [permutations]
#   tax_level    optional taxonomic rank to agglomerate to before ANCOM-BC2 (e.g. "genus"); default = ASV/OTU level
#   permutations optional PERMANOVA/betadisper permutation count; default = 999
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript run_phyloseq_ancombc.R abundance.csv taxonomy.csv metadata.csv group_col output_dir [tax_level] [permutations]")
}
abundance_file <- args[[1]]
taxonomy_file  <- args[[2]]
metadata_file  <- args[[3]]
group_col      <- args[[4]]
out_dir        <- args[[5]]
tax_level      <- if (length(args) >= 6 && nzchar(args[[6]])) args[[6]] else NA
n_perm         <- if (length(args) >= 7 && nzchar(args[[7]])) as.integer(args[[7]]) else 999L

required <- c("phyloseq", "ANCOMBC", "ggplot2", "vegan", "mia")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "),
       "\nInstall Bioconductor/CRAN dependencies before running formal analysis.")
}

dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "run"), recursive = TRUE, showWarnings = FALSE)

PALETTE <- c("#4EA3F1", "#E45B5B", "#4CB984", "#E2A13B", "#8B72D6", "#56BFC0", "#F08AB4", "#8AA1B4", "#B8C2CC", "#D9C25C")
TEXT <- "#173B5F"

theme_micro <- function() {
  ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 24, color = TEXT),
      axis.title = ggplot2::element_text(size = 21, color = TEXT),
      axis.text = ggplot2::element_text(size = 16, color = TEXT),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 13)
    )
}

save_gg <- function(plot, name, width = 8, height = 6) {
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(name, ".png")), plot, width = width, height = height, dpi = 320, bg = "white")
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(name, ".pdf")), plot, width = width, height = height, bg = "white")
}

write_json <- function(obj, path) {
  kv <- vapply(names(obj), function(n) {
    v <- obj[[n]]
    if (is.numeric(v) || is.logical(v)) sprintf('"%s": %s', n, jsonlite_scalar(v)) else sprintf('"%s": "%s"', n, v)
  }, character(1))
  writeLines(paste0("{\n  ", paste(kv, collapse = ",\n  "), "\n}"), path)
}
jsonlite_scalar <- function(v) if (is.logical(v)) tolower(as.character(v)) else format(v, scientific = FALSE)

## ---- 1. Load tables and build phyloseq object ------------------------------
otu  <- read.csv(abundance_file, check.names = FALSE)
tax  <- read.csv(taxonomy_file, check.names = FALSE)
meta <- read.csv(metadata_file, check.names = FALSE)
if (nrow(otu) == 0) stop("Abundance table '", abundance_file, "' has no data rows.")
if (nrow(tax) == 0) stop("Taxonomy table '", taxonomy_file, "' has no data rows.")
if (nrow(meta) == 0) stop("Metadata table '", metadata_file, "' has no data rows.")

id_col_otu <- if ("taxon" %in% colnames(otu)) "taxon" else if ("feature" %in% colnames(otu)) "feature" else
  stop("Abundance table must have a 'taxon' or 'feature' id column; found: ", paste(colnames(otu), collapse = ", "))
id_col_tax <- if ("taxon" %in% colnames(tax)) "taxon" else if ("feature" %in% colnames(tax)) "feature" else
  stop("Taxonomy table must have a 'taxon' or 'feature' id column; found: ", paste(colnames(tax), collapse = ", "))
stopifnot("sample" %in% colnames(meta), group_col %in% colnames(meta))

if (anyDuplicated(otu[[id_col_otu]])) stop("Abundance table has duplicate taxon/feature IDs; make IDs unique before running.")
if (anyDuplicated(tax[[id_col_tax]])) stop("Taxonomy table has duplicate taxon/feature IDs; make IDs unique before running.")
if (anyDuplicated(meta$sample)) stop("Metadata table has duplicate sample IDs; make sample IDs unique before running.")

rownames(otu) <- otu[[id_col_otu]]
otu[[id_col_otu]] <- NULL
rownames(tax) <- tax[[id_col_tax]]
tax[[id_col_tax]] <- NULL
rownames(meta) <- meta$sample

if (!all(vapply(otu, is.numeric, logical(1)))) {
  stop("Abundance table contains non-numeric values in sample columns; check for blanks/text.")
}

common <- intersect(colnames(otu), rownames(meta))
dropped_samples <- setdiff(colnames(otu), rownames(meta))
if (length(common) == 0) stop("No overlapping sample IDs between abundance table and metadata.")
if (length(dropped_samples)) {
  warning(sprintf("Dropping %d sample(s) present in the abundance table but missing from metadata: %s",
                   length(dropped_samples), paste(dropped_samples, collapse = ", ")))
}
otu  <- otu[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]
meta[[group_col]] <- factor(meta[[group_col]])
n_groups <- nlevels(meta[[group_col]])
singleton_groups <- names(which(table(meta[[group_col]]) == 1))

# Drop all-zero taxa (housekeeping; ANCOM-BC2's prv_cut handles the formal prevalence filter below).
zero_taxa <- rowSums(as.matrix(otu)) == 0
if (any(zero_taxa)) {
  message(sprintf("Removing %d all-zero taxa before analysis", sum(zero_taxa)))
  otu <- otu[!zero_taxa, , drop = FALSE]
}
if (nrow(otu) == 0) stop("All taxa are all-zero after filtering; nothing left to analyze.")

# Fill taxonomy gaps: taxa present in the abundance table but missing from the taxonomy table
# would otherwise create duplicate/NA row names when subsetting and crash phyloseq::tax_table().
missing_taxa <- setdiff(rownames(otu), rownames(tax))
if (length(missing_taxa)) {
  message(sprintf("Taxonomy table is missing %d taxon/taxa present in the abundance table; filling with 'Unclassified': %s",
                   length(missing_taxa), paste(utils::head(missing_taxa, 5), collapse = ", ")))
  filler <- as.data.frame(matrix("Unclassified", nrow = length(missing_taxa), ncol = ncol(tax),
                                  dimnames = list(missing_taxa, colnames(tax))))
  tax <- rbind(tax, filler)
}
tax <- tax[rownames(otu), , drop = FALSE]

low_depth <- colSums(otu) < 1000
low_depth_samples <- names(which(low_depth))
if (any(low_depth)) {
  message("Low-depth samples detected (<1000 reads): ", paste(low_depth_samples, collapse = ", "),
          " -- review before treating alpha/beta diversity as comparable across samples (references/pitfalls.md #3).")
}
if (length(singleton_groups)) {
  message("Group(s) with only 1 sample: ", paste(singleton_groups, collapse = ", "),
          " -- group-level statistics for these groups have essentially no within-group variance estimate.")
}

ps <- phyloseq::phyloseq(
  phyloseq::otu_table(as.matrix(otu), taxa_are_rows = TRUE),
  phyloseq::tax_table(as.matrix(tax)),
  phyloseq::sample_data(meta)
)
sd_df <- as(phyloseq::sample_data(ps), "data.frame")

## ---- 2. Read-depth / QC -----------------------------------------------------
qc <- data.frame(
  sample = phyloseq::sample_names(ps),
  library_size = as.numeric(phyloseq::sample_sums(ps)),
  observed_taxa = as.numeric(colSums(as(phyloseq::otu_table(ps), "matrix") > 0))
)
qc$group <- sd_df[qc$sample, group_col]
write.table(qc, file.path(out_dir, "run", "qc.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write_json(list(samples = phyloseq::nsamples(ps), taxa = phyloseq::ntaxa(ps),
                min_library_size = min(qc$library_size), max_library_size = max(qc$library_size)),
           file.path(out_dir, "run", "qc.json"))

p_depth <- ggplot2::ggplot(qc, ggplot2::aes(x = stats::reorder(sample, -library_size), y = library_size, fill = group)) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::scale_fill_manual(values = PALETTE) +
  ggplot2::labs(title = "Read depth per sample", x = "Sample", y = "Library size (reads)", fill = "Group") +
  theme_micro() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
save_gg(p_depth, "read_depth", width = max(8.5, nrow(qc) * 0.5 + 3), height = 5.8)

## ---- 3. Alpha diversity: Observed / Shannon / Simpson -----------------------
alpha <- phyloseq::estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson"))
alpha$sample <- rownames(alpha)
alpha$group <- sd_df[alpha$sample, group_col]
alpha <- alpha[, c("sample", "group", "Observed", "Shannon", "Simpson")]
write.table(alpha, file.path(out_dir, "results", "alpha_diversity.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

alpha_long <- do.call(rbind, lapply(c("Observed", "Shannon", "Simpson"), function(m) {
  data.frame(sample = alpha$sample, group = alpha$group, metric = m, value = alpha[[m]])
}))
alpha_long$metric <- factor(alpha_long$metric, levels = c("Observed", "Shannon", "Simpson"))

alpha_stats <- do.call(rbind, lapply(c("Observed", "Shannon", "Simpson"), function(m) {
  vals <- split(alpha[[m]], alpha$group)
  means <- vapply(vals, mean, numeric(1))
  higher_group <- if (length(means) >= 2) names(means)[which.max(means)] else NA_character_
  means_txt <- paste(sprintf("%s=%.3g", names(means), means), collapse = "; ")
  if (length(vals) < 2) {
    return(data.frame(metric = m, test = NA_character_, p_value = NA_real_,
                       higher_group = higher_group, group_means = means_txt))
  }
  test <- tryCatch(
    if (length(vals) == 2) stats::wilcox.test(vals[[1]], vals[[2]]) else stats::kruskal.test(vals),
    error = function(e) NULL
  )
  data.frame(metric = m, test = if (is.null(test)) NA_character_ else if (length(vals) == 2) "wilcox.test" else "kruskal.test",
             p_value = if (is.null(test)) NA_real_ else test$p.value,
             higher_group = higher_group, group_means = means_txt)
}))
alpha_stats$p_adj <- stats::p.adjust(alpha_stats$p_value, method = "BH")
write.table(alpha_stats, file.path(out_dir, "results", "alpha_diversity_stats.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

p_alpha <- ggplot2::ggplot(alpha_long, ggplot2::aes(group, value, fill = group)) +
  ggplot2::geom_boxplot(alpha = 0.55, outlier.shape = NA, width = 0.55) +
  ggplot2::geom_jitter(ggplot2::aes(color = group), width = 0.12, size = 1.6, alpha = 0.75, show.legend = FALSE) +
  ggplot2::facet_wrap(~metric, scales = "free_y") +
  ggplot2::scale_fill_manual(values = PALETTE) +
  ggplot2::scale_color_manual(values = PALETTE) +
  ggplot2::labs(title = "Alpha diversity", x = "Group", y = "Diversity value", fill = "Group") +
  theme_micro()
save_gg(p_alpha, "alpha_diversity", width = 10.5, height = 5.6)

## ---- 4. Beta diversity: Bray-Curtis + PCoA + PERMANOVA + betadisper ---------
bray <- phyloseq::distance(ps, method = "bray")
bray_labels <- attr(bray, "Labels")
grp_vec <- factor(sd_df[bray_labels, group_col])

pcoa <- stats::cmdscale(bray, k = 2, eig = TRUE)
var_explained <- pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100
ord_df <- data.frame(sample = rownames(pcoa$points), Axis1 = pcoa$points[, 1], Axis2 = pcoa$points[, 2])
ord_df$group <- sd_df[ord_df$sample, group_col]

p_ord <- ggplot2::ggplot(ord_df, ggplot2::aes(Axis1, Axis2, color = group)) +
  ggplot2::geom_point(size = 3.2, alpha = 0.86) +
  ggplot2::scale_color_manual(values = PALETTE) +
  ggplot2::labs(
    title = "PCoA (Bray-Curtis)",
    x = sprintf("PCoA1 (%.1f%%)", var_explained[1]),
    y = sprintf("PCoA2 (%.1f%%)", var_explained[2]),
    color = "Group"
  ) +
  theme_micro()
save_gg(p_ord, "ordination", width = 8.4, height = 6.4)

perm_data <- data.frame(group = grp_vec)
if (nlevels(grp_vec) >= 2) {
  permanova <- tryCatch(vegan::adonis2(bray ~ group, data = perm_data, permutations = n_perm),
                         error = function(e) { message("PERMANOVA failed: ", conditionMessage(e)); NULL })
  disp <- tryCatch(vegan::betadisper(bray, grp_vec), error = function(e) { message("betadisper failed: ", conditionMessage(e)); NULL })
  disp_test <- if (!is.null(disp)) tryCatch(vegan::permutest(disp, permutations = n_perm),
                                             error = function(e) { message("betadisper permutest failed: ", conditionMessage(e)); NULL }) else NULL
} else {
  message("Only one group present; PERMANOVA/betadisper group comparison skipped.")
  permanova <- NULL
  disp_test <- NULL
}

permanova_tbl <- if (!is.null(permanova)) {
  data.frame(
    term = rownames(permanova),
    Df = permanova$Df,
    SumOfSqs = permanova$SumOfSqs,
    R2 = permanova$R2,
    F = permanova$F,
    p_value = permanova[["Pr(>F)"]]
  )
} else {
  data.frame(term = "group", Df = NA_real_, SumOfSqs = NA_real_, R2 = NA_real_, F = NA_real_, p_value = NA_real_)
}
dispersion_tbl <- if (!is.null(disp_test)) {
  data.frame(
    term = "betadisper (group dispersion)",
    Df = disp_test$tab$Df[1],
    SumOfSqs = NA_real_,
    R2 = NA_real_,
    F = disp_test$tab$F[1],
    p_value = disp_test$tab$`Pr(>F)`[1]
  )
} else {
  data.frame(term = "betadisper (group dispersion)", Df = NA_real_, SumOfSqs = NA_real_, R2 = NA_real_, F = NA_real_, p_value = NA_real_)
}
write.table(rbind(permanova_tbl, dispersion_tbl), file.path(out_dir, "results", "permanova.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

## ---- 5. Taxa composition: top taxa + Other ----------------------------------
otu_mat <- as(phyloseq::otu_table(ps), "matrix")
rel_mat <- sweep(otu_mat, 2, colSums(otu_mat), FUN = "/")
mean_rel <- sort(rowMeans(rel_mat), decreasing = TRUE)
top_n <- min(9, length(mean_rel))
top_taxa <- names(mean_rel)[seq_len(top_n)]
other_taxa <- setdiff(rownames(rel_mat), top_taxa)
taxa_labels <- c(top_taxa, if (length(other_taxa)) "Other")

comp_long <- do.call(rbind, lapply(taxa_labels, function(t) {
  vals <- if (identical(t, "Other")) colSums(rel_mat[other_taxa, , drop = FALSE]) else rel_mat[t, ]
  data.frame(sample = colnames(rel_mat), taxon = t, abundance = as.numeric(vals))
}))
comp_long$taxon <- factor(comp_long$taxon, levels = taxa_labels)
comp_long$group <- sd_df[comp_long$sample, group_col]
sample_order <- unique(comp_long$sample[order(comp_long$group, comp_long$sample)])
comp_long$sample <- factor(comp_long$sample, levels = sample_order)

p_comp <- ggplot2::ggplot(comp_long, ggplot2::aes(sample, abundance * 100, fill = taxon)) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::scale_fill_manual(values = PALETTE[seq_along(taxa_labels)]) +
  ggplot2::labs(title = "Taxa composition", x = "Sample", y = "Relative abundance (%)", fill = "Taxon") +
  theme_micro() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
save_gg(p_comp, "taxa_composition", width = max(8.5, length(sample_order) * 0.55 + 3), height = 6.2)

rel_out <- data.frame(taxon = rownames(rel_mat), rel_mat, check.names = FALSE)
write.table(rel_out, file.path(out_dir, "results", "taxa_relative_abundance.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

## ---- 6. Differential abundance: ANCOM-BC2 (compositionality-aware) ---------
n_significant <- NA_integer_
top_taxa_txt <- character(0)
res_df <- NULL
if (n_groups < 2) {
  message("Only one group present; ANCOM-BC2 differential abundance skipped.")
} else {
  # Convert to TreeSummarizedExperiment as recommended by the ANCOM-BC2 vignette
  # (phyloseq -> mia::convertFromPhyloseq -> ANCOMBC::ancombc2). makeTreeSEFromPhyloseq()
  # is deprecated in newer mia releases in favor of convertFromPhyloseq(); fall back to the
  # old name on older mia installs where the new function does not exist yet.
  tse <- if (exists("convertFromPhyloseq", where = asNamespace("mia"))) {
    mia::convertFromPhyloseq(ps)
  } else {
    mia::makeTreeSEFromPhyloseq(ps)
  }
  ancombc_out <- tryCatch(
    ANCOMBC::ancombc2(
      data = tse,
      assay_name = "counts",
      tax_level = if (is.na(tax_level)) NULL else tax_level,
      fix_formula = group_col,
      p_adj_method = "BH",
      prv_cut = 0.10,       # taxa present in >=10% of samples, per references/cns-method-guide.md decision tree
      lib_cut = 0,
      group = group_col,    # enables structural-zero detection across group_col levels
      struc_zero = TRUE,
      neg_lb = TRUE,
      alpha = 0.05,
      n_cl = 1
    ),
    error = function(e) {
      message("ANCOM-BC2 failed (often caused by too few samples per group for a stable model): ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(ancombc_out)) res_df <- ancombc_out$res
}
if (!is.null(res_df)) {
  write.table(res_df, file.path(out_dir, "results", "differential_taxa.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  message("Skipping results/differential_taxa.tsv and differential_taxa figure (ANCOM-BC2 did not run or produced no result).")
}

lfc_cols <- if (!is.null(res_df)) {
  cols <- grep("^lfc_", colnames(res_df), value = TRUE)
  cols[!grepl("Intercept", cols)]
} else character(0)
if (length(lfc_cols) >= 1) {
  contrast <- lfc_cols[1]
  suffix <- sub("^lfc_", "", contrast)
  q_col <- paste0("q_", suffix)
  diff_col <- paste0("diff_", suffix)
  plot_df <- data.frame(
    taxon = res_df$taxon,
    lfc = res_df[[contrast]],
    q_val = res_df[[q_col]]
  )
  plot_df$significant <- if (diff_col %in% colnames(res_df)) res_df[[diff_col]] else (!is.na(plot_df$q_val) & plot_df$q_val < 0.05)
  plot_df$neglog10q <- -log10(pmax(plot_df$q_val, 1e-300))
  plot_df$call <- ifelse(plot_df$significant & plot_df$lfc > 0, "enriched",
                   ifelse(plot_df$significant & plot_df$lfc < 0, "depleted", "ns"))
  n_significant <- sum(plot_df$significant, na.rm = TRUE)
  sig_df <- plot_df[!is.na(plot_df$significant) & plot_df$significant, ]
  sig_df <- sig_df[order(sig_df$q_val), ]
  top_sig <- utils::head(sig_df, 5)
  if (nrow(top_sig)) {
    top_taxa_txt <- sprintf("%s: %s in %s (log fold change = %.2f, q = %.3g)",
                             top_sig$taxon, top_sig$call, suffix, top_sig$lfc, top_sig$q_val)
  }

  p_diff <- ggplot2::ggplot(plot_df, ggplot2::aes(lfc, neglog10q, color = call)) +
    ggplot2::geom_point(size = 2.6, alpha = 0.85) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#A7C6DD") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "#A7C6DD") +
    ggplot2::scale_color_manual(values = c(enriched = "#E45B5B", depleted = "#4EA3F1", ns = "#B8C2CC")) +
    ggplot2::labs(
      title = paste0("Differential abundance (ANCOM-BC2): ", suffix),
      x = "Log fold change (bias-corrected)", y = "-log10(q value)", color = NULL
    ) +
    theme_micro()
  save_gg(p_diff, "differential_taxa", width = 8.6, height = 6.6)
} else {
  message("No fixed-effect contrast column (lfc_*) found in ANCOM-BC2 output; skipping differential taxa plot.")
}

## ---- 7. Summary, report, artifacts -----------------------------------------
summary_df <- data.frame(
  item = c("samples", "taxa", "significant_taxa_q0.05", "permanova_R2", "permanova_p", "betadisper_p"),
  value = c(phyloseq::nsamples(ps), phyloseq::ntaxa(ps), n_significant,
            round(permanova_tbl$R2[1], 4), permanova_tbl$p_value[1], dispersion_tbl$p_value[1])
)
write.table(summary_df, file.path(out_dir, "summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

group_counts <- table(meta[[group_col]])
group_counts_txt <- paste(sprintf("%s: n=%d", names(group_counts), as.integer(group_counts)), collapse = ", ")

## Alpha diversity: per-metric direction + significance -----------------------
alpha_lines <- unlist(lapply(seq_len(nrow(alpha_stats)), function(i) {
  row <- alpha_stats[i, ]
  if (is.na(row$p_value)) {
    sprintf("- %s: %s; statistical test skipped or unavailable (insufficient replicates in at least one group).",
            row$metric, row$group_means)
  } else {
    sig <- if (!is.na(row$p_adj) && row$p_adj < 0.05) "significant" else "not significant"
    sprintf("- %s: %s; highest in **%s**; %s p = %.4g, BH-adjusted p = %.4g (%s).",
            row$metric, row$group_means, row$higher_group, row$test, row$p_value, row$p_adj, sig)
  }
}))

## Beta diversity: PERMANOVA + conditional dispersion caveat ------------------
beta_lines <- if (is.na(permanova_tbl$p_value[1])) {
  "- Only one group present, or PERMANOVA failed to run; beta diversity group comparison skipped."
} else {
  perm_sig <- if (permanova_tbl$p_value[1] < 0.05) "significant group separation" else "no significant group separation"
  disp_line <- if (is.na(dispersion_tbl$p_value[1])) {
    "- betadisper dispersion check could not be computed for this data."
  } else if (dispersion_tbl$p_value[1] < 0.05) {
    sprintf(paste0("- betadisper dispersion p = %.4g (significant) -- within-group dispersion differs between ",
                   "groups, so part of the PERMANOVA result may reflect unequal spread rather than a pure ",
                   "centroid shift (references/pitfalls.md #5); interpret the PERMANOVA result with this caveat."),
            dispersion_tbl$p_value[1])
  } else {
    sprintf(paste0("- betadisper dispersion p = %.4g (not significant) -- dispersion is comparable between ",
                   "groups, supporting a centroid-shift interpretation of the PERMANOVA result."),
            dispersion_tbl$p_value[1])
  }
  c(sprintf("- PERMANOVA (%d permutations): R2 = %.4g, F = %.4g, p = %.4g -> %s",
            n_perm, permanova_tbl$R2[1], permanova_tbl$F[1], permanova_tbl$p_value[1], perm_sig),
    disp_line)
}

## Differential abundance: top taxa named with direction -----------------------
diff_lines <- if (n_groups < 2) {
  "- Skipped: only one group present in metadata."
} else if (is.null(res_df)) {
  "- ANCOM-BC2 did not produce a result for this dataset (see run log); no differential taxa reported."
} else {
  c("- Prevalence filter: taxa present in >=10% of samples (prv_cut = 0.10)",
    "- Structural zero detection: enabled (struc_zero = TRUE, neg_lb = TRUE)",
    sprintf("- FDR method: BH; significant taxa (q < 0.05): %s", n_significant),
    if (length(top_taxa_txt)) c("- Top differentially abundant taxa:", paste0("  - ", top_taxa_txt))
    else "- No taxa passed the q < 0.05 threshold.")
}

## Automatic caveats tied to references/pitfalls.md, only when they apply -----
caveats <- character(0)
if (length(low_depth_samples)) {
  caveats <- c(caveats, sprintf(
    "Low read-depth sample(s) (<1000 reads): %s -- review before treating diversity/composition as comparable across samples (references/pitfalls.md #3).",
    paste(low_depth_samples, collapse = ", ")))
}
if (length(dropped_samples)) {
  caveats <- c(caveats, sprintf(
    "%d sample(s) in the abundance table had no metadata match and were dropped: %s -- verify sample IDs.",
    length(dropped_samples), paste(dropped_samples, collapse = ", ")))
}
if (length(missing_taxa)) {
  caveats <- c(caveats, sprintf(
    "%d taxon/taxa were missing from the taxonomy table and filled as 'Unclassified': %s -- verify the taxonomy table is complete.",
    length(missing_taxa), paste(utils::head(missing_taxa, 5), collapse = ", ")))
}
if (length(singleton_groups)) {
  caveats <- c(caveats, sprintf(
    "Group(s) with only 1 sample: %s -- statistics for these groups have essentially no within-group variance estimate.",
    paste(singleton_groups, collapse = ", ")))
}
if (!is.na(dispersion_tbl$p_value[1]) && dispersion_tbl$p_value[1] < 0.05) {
  caveats <- c(caveats, "betadisper flagged a significant dispersion difference between groups (see Beta diversity section above).")
}

## Plain-language conclusion ---------------------------------------------------
conclusion <- character(0)
if (n_groups < 2) {
  conclusion <- "Only one group is present in the metadata, so no between-group diversity or composition comparison could be made."
} else {
  shannon_row <- alpha_stats[alpha_stats$metric == "Shannon", ]
  if (nrow(shannon_row) && !is.na(shannon_row$p_adj[1])) {
    if (shannon_row$p_adj[1] < 0.05) {
      conclusion <- c(conclusion, sprintf("Shannon diversity was significantly higher in %s (BH-adjusted p = %.4g).",
                                           shannon_row$higher_group[1], shannon_row$p_adj[1]))
    } else {
      conclusion <- c(conclusion, sprintf(
        "Shannon diversity did not differ significantly between groups (BH-adjusted p = %.4g), though the mean was nominally higher in %s.",
        shannon_row$p_adj[1], shannon_row$higher_group[1]))
    }
  } else {
    conclusion <- c(conclusion, "Alpha diversity could not be tested statistically due to insufficient replicates in at least one group.")
  }
  if (!is.na(permanova_tbl$p_value[1])) {
    if (permanova_tbl$p_value[1] < 0.05) {
      disp_note <- if (!is.na(dispersion_tbl$p_value[1]) && dispersion_tbl$p_value[1] < 0.05) {
        ", though this may partly reflect a dispersion difference rather than a pure centroid shift"
      } else ""
      conclusion <- c(conclusion, sprintf(
        "Community composition (Bray-Curtis PERMANOVA) differed significantly between groups (R2 = %.4g, p = %.4g)%s.",
        permanova_tbl$R2[1], permanova_tbl$p_value[1], disp_note))
    } else {
      conclusion <- c(conclusion, sprintf(
        "Community composition did not differ significantly between groups by PERMANOVA (R2 = %.4g, p = %.4g).",
        permanova_tbl$R2[1], permanova_tbl$p_value[1]))
    }
  }
  if (length(top_taxa_txt)) {
    conclusion <- c(conclusion, sprintf(
      "ANCOM-BC2 identified %d taxon/taxa as differentially abundant at q < 0.05; the strongest hit was %s.",
      n_significant, top_sig$taxon[1]))
  } else if (!is.null(res_df)) {
    conclusion <- c(conclusion, "No taxa were called differentially abundant by ANCOM-BC2 at q < 0.05.")
  }
}

report_lines <- c(
  "# Microbiome phyloseq / ANCOM-BC2 report", "",
  paste0("- Samples: ", phyloseq::nsamples(ps), " (", group_counts_txt, ")"),
  paste0("- Taxa (features): ", phyloseq::ntaxa(ps)),
  paste0("- Group column: ", group_col),
  "",
  "## Alpha diversity", "",
  alpha_lines, "",
  "## Beta diversity (Bray-Curtis PERMANOVA)", "",
  beta_lines, "",
  "## Differential abundance (ANCOM-BC2)", "",
  diff_lines, ""
)
if (length(caveats)) {
  report_lines <- c(report_lines, "## Caveats", "", paste0("- ", caveats), "")
}
report_lines <- c(report_lines,
  "## Compositional limitations", "",
  "- Relative abundances are constrained; ANCOM-BC2 corrects for sampling-fraction bias but reported",
  "  fold changes remain relative, not absolute biomass (references/pitfalls.md #2, #13).",
  "- See references/pitfalls.md for rarefaction, dispersion, batch/kitome and database-version caveats.",
  "",
  "## Conclusion", "",
  paste(conclusion, collapse = " ")
)
writeLines(report_lines, file.path(out_dir, "report.md"))

artifacts <- list.files(out_dir, recursive = TRUE)
writeLines(paste0("[\n  \"", paste(artifacts, collapse = "\",\n  \""), "\"\n]"),
           file.path(out_dir, "run", "artifacts.json"))

message("Done. Outputs written under: ", out_dir)
