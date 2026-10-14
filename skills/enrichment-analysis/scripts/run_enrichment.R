#!/usr/bin/env Rscript
# ORA (over-representation analysis) via clusterProfiler::enrichGO / enrichKEGG.
#
# Usage:
#   Rscript run_enrichment.R <input_file> <organism_db> <keytype> <database> <output_dir> \
#       [universe_file] [padj_cutoff] [kegg_organism] [enrichment_padj_threshold]
#
# <input_file> is one of:
#   - a plain gene list, one ID per line (foreground only). In this case
#     <universe_file> is REQUIRED (a background gene list, one ID per line) —
#     see references/pitfalls.md #1: this script refuses to silently fall back
#     to a whole-genome background.
#   - a differential-expression table (.tsv/.csv) with a header containing at
#     least "gene" and "padj" columns (see references/input-contract.md). In
#     this case the foreground is genes with padj < <padj_cutoff> (default
#     0.05) and the background universe is ALL genes in the table (i.e. every
#     gene that was actually tested), unless <universe_file> overrides it.
#
# <database> is one of "BP", "CC", "MF" (GO ontology) or "KEGG".
# <kegg_organism> is the KEGG organism code (e.g. "hsa"); if omitted, it is
# looked up from <organism_db> for a handful of common organisms.
#
# [padj_cutoff] (arg 7) and [enrichment_padj_threshold] (arg 9) are deliberately
# separate knobs:
#   - padj_cutoff: which DE-table genes become the foreground gene list tested.
#   - enrichment_padj_threshold: which *tested terms* in the ORA output are
#     flagged "significant" (results table, dotplot/barplot title, report.md).
# They commonly differ (e.g. a lenient padj_cutoff to keep a workable foreground,
# but a strict enrichment_padj_threshold for reporting).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(
    "Usage: Rscript run_enrichment.R <input_file> <organism_db> <keytype> <database:BP|CC|MF|KEGG> ",
    "<output_dir> [universe_file] [padj_cutoff] [kegg_organism] [enrichment_padj_threshold]"
  )
}
input_file     <- args[[1]]
organism_db    <- args[[2]]
keytype        <- args[[3]]
database       <- toupper(args[[4]])
out_dir        <- args[[5]]
universe_file  <- if (length(args) >= 6 && nzchar(args[[6]])) args[[6]] else NA_character_
padj_cutoff    <- if (length(args) >= 7 && nzchar(args[[7]])) as.numeric(args[[7]]) else 0.05
kegg_organism  <- if (length(args) >= 8 && nzchar(args[[8]])) args[[8]] else NA_character_
enrichment_padj_threshold <- if (length(args) >= 9 && nzchar(args[[9]])) as.numeric(args[[9]]) else 0.05

if (!database %in% c("BP", "CC", "MF", "KEGG")) {
  stop("database must be one of BP, CC, MF, KEGG (got: ", database, ")")
}

required <- c("clusterProfiler", organism_db, "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "run"), recursive = TRUE, showWarnings = FALSE)

PALETTE <- c(BP = "#3FA37A", CC = "#D99A2B", MF = "#7C65C1", KEGG = "#3B82C4")
TEXT_COLOR <- "#1F2937"

save_gg <- function(plot, name, width = 9, height = 7) {
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(name, ".png")), plot, width = width, height = height, dpi = 320, bg = "white")
  ggplot2::ggsave(file.path(out_dir, "figures", paste0(name, ".pdf")), plot, width = width, height = height, bg = "white")
}

wrap_labels <- function(x, width = 45) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

split_genes <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character(0))
  trimws(unlist(strsplit(x, "/")))
}

## Simple Jaccard-overlap heuristic for pitfalls.md #3 (redundant GO/KEGG terms):
## flag pairs of top significant terms whose member-gene sets overlap heavily.
jaccard_redundancy_flags <- function(df, padj_threshold = 0.05, jaccard_threshold = 0.5, top_n = 10) {
  if (!"geneID" %in% colnames(df)) return(list())
  sig <- df[!is.na(df$p.adjust) & df$p.adjust < padj_threshold, ]
  if (nrow(sig) < 2) return(list())
  sig <- head(sig[order(sig$p.adjust), ], top_n)
  gene_sets <- lapply(sig$geneID, split_genes)
  flags <- list()
  n <- nrow(sig)
  for (i in seq_len(n - 1)) {
    gi <- gene_sets[[i]]
    if (length(gi) == 0) next
    for (j in seq(i + 1, n)) {
      gj <- gene_sets[[j]]
      if (length(gj) == 0) next
      inter <- length(intersect(gi, gj))
      uni <- length(union(gi, gj))
      jac <- if (uni > 0) inter / uni else 0
      if (jac >= jaccard_threshold) {
        flags[[length(flags) + 1]] <- list(a = sig$Description[i], b = sig$Description[j], jac = jac, n = inter)
      }
    }
  }
  flags[order(-vapply(flags, function(f) f$jac, numeric(1)))]
}

## ---- 1. Load input and build foreground / background universe -------------

first_line <- readLines(input_file, n = 1, warn = FALSE)
looks_like_table <- grepl("gene", first_line, ignore.case = TRUE) &&
  (grepl("\t", first_line) || grepl(",", first_line))

if (looks_like_table) {
  sep <- if (grepl("\t", first_line)) "\t" else ","
  de_table <- read.delim(input_file, sep = sep, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("gene", "padj") %in% colnames(de_table))) {
    stop("Input looks like a DE table but is missing required columns 'gene'/'padj'. ",
         "See references/input-contract.md.")
  }
  de_table <- de_table[!is.na(de_table$gene) & nzchar(de_table$gene), ]
  universe_ids <- unique(de_table$gene)
  foreground_ids <- unique(de_table$gene[!is.na(de_table$padj) & de_table$padj < padj_cutoff])
  input_mode <- "differential_expression_table"
} else {
  foreground_ids <- unique(trimws(readLines(input_file, warn = FALSE)))
  foreground_ids <- foreground_ids[nzchar(foreground_ids)]
  universe_ids <- NULL
  input_mode <- "gene_list"
}

if (!is.na(universe_file)) {
  universe_ids <- unique(trimws(readLines(universe_file, warn = FALSE)))
  universe_ids <- universe_ids[nzchar(universe_ids)]
  universe_source <- paste0("explicit universe_file (", basename(universe_file), ")")
} else if (input_mode == "differential_expression_table") {
  universe_source <- "all genes present in the input differential-expression table (i.e. every tested gene, not a whole-genome default)"
} else {
  universe_source <- NA_character_
}

if (is.null(universe_ids)) {
  stop(
    "No background/universe available. A plain gene-list input requires an explicit ",
    "universe_file argument (the full set of genes tested/detected in the experiment). ",
    "Refusing to silently default to a whole-genome background — see references/pitfalls.md #1."
  )
}
if (length(foreground_ids) == 0) stop("Foreground gene list is empty after filtering — nothing to test.")
if (!all(foreground_ids %in% universe_ids)) {
  n_missing <- sum(!foreground_ids %in% universe_ids)
  message(n_missing, " foreground gene(s) not present in the supplied universe; adding them so the universe always contains the foreground.")
  universe_ids <- unique(c(universe_ids, foreground_ids))
}

OrgDb <- getExportedValue(organism_db, organism_db)

## ---- 2. ID mapping (KEGG requires ENTREZID) --------------------------------

unmapped <- data.frame(gene = character(0), stage = character(0))
mapping_rates <- list()
map_to_entrez <- function(ids, label) {
  m <- clusterProfiler::bitr(ids, fromType = keytype, toType = "ENTREZID", OrgDb = OrgDb, drop = TRUE)
  mapped_rate <- nrow(m) / length(ids)
  mapping_rates[[label]] <<- list(mapped = nrow(m), total = length(ids), rate = mapped_rate)
  message(sprintf("%s: %d/%d (%.1f%%) mapped to ENTREZID", label, nrow(m), length(ids), 100 * mapped_rate))
  unmapped_ids <- setdiff(ids, m[[keytype]])
  if (length(unmapped_ids)) {
    unmapped <<- rbind(unmapped, data.frame(gene = unmapped_ids, stage = label))
  }
  unique(m$ENTREZID)
}

if (database == "KEGG" && keytype != "ENTREZID") {
  foreground_use <- map_to_entrez(foreground_ids, "foreground")
  universe_use <- map_to_entrez(universe_ids, "universe")
} else {
  foreground_use <- foreground_ids
  universe_use <- universe_ids
}

if (nrow(unmapped) > 0) {
  write.table(unmapped, file.path(out_dir, "results", "unmapped_genes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}

## ---- 3. Run ORA -------------------------------------------------------------

if (database == "KEGG") {
  if (is.na(kegg_organism)) {
    kegg_lookup <- c(
      "org.Hs.eg.db" = "hsa", "org.Mm.eg.db" = "mmu", "org.Rn.eg.db" = "rno",
      "org.Dr.eg.db" = "dre", "org.Dm.eg.db" = "dme", "org.Ce.eg.db" = "cel",
      "org.Sc.sgd.db" = "sce"
    )
    kegg_organism <- kegg_lookup[[organism_db]]
    if (is.null(kegg_organism) || is.na(kegg_organism)) {
      stop("Cannot infer KEGG organism code from '", organism_db, "'. Pass it explicitly as the 8th argument.")
    }
  }
  res <- clusterProfiler::enrichKEGG(
    gene = foreground_use, organism = kegg_organism, keyType = "kegg",
    universe = universe_use, pAdjustMethod = "BH",
    pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500
  )
  out_name <- "kegg_enrichment"
} else {
  res <- clusterProfiler::enrichGO(
    gene = foreground_use, OrgDb = OrgDb, keyType = keytype, ont = database,
    universe = universe_use, pAdjustMethod = "BH",
    pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500
  )
  out_name <- "go_enrichment"
}

if (is.null(res) || nrow(as.data.frame(res)) == 0) {
  stop("No terms were tested — check organism_db/keytype/database and gene ID formatting.")
}

tab <- as.data.frame(res)
tab$significant <- !is.na(tab$p.adjust) & tab$p.adjust < enrichment_padj_threshold
write.table(tab, file.path(out_dir, "results", paste0(out_name, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

## Optional: reduce redundant GO terms (references/pitfalls.md #3), best-effort.
simplified_written <- FALSE
if (database %in% c("BP", "CC", "MF") && requireNamespace("GOSemSim", quietly = TRUE)) {
  simplified <- tryCatch(clusterProfiler::simplify(res, cutoff = 0.7, by = "p.adjust", select_fun = min), error = function(e) NULL)
  if (!is.null(simplified) && nrow(as.data.frame(simplified)) > 0) {
    stab <- as.data.frame(simplified)
    stab$significant <- !is.na(stab$p.adjust) & stab$p.adjust < enrichment_padj_threshold
    write.table(stab, file.path(out_dir, "results", paste0(out_name, "_simplified.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
    simplified_written <- TRUE
  }
}

## ---- 4. Plots (manual ggplot for full control over house style) -----------

plot_df <- tab[order(tab$p.adjust), ]
plot_df <- head(plot_df, 20)
plot_df$term_wrapped <- wrap_labels(plot_df$Description, 45)
plot_df$term_wrapped <- factor(plot_df$term_wrapped, levels = rev(plot_df$term_wrapped))
plot_df$GeneRatioNum <- vapply(strsplit(plot_df$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))

dot_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(GeneRatioNum, term_wrapped, size = Count, color = -log10(pmax(p.adjust, 1e-300)))) +
  ggplot2::geom_point(alpha = 0.9) +
  ggplot2::scale_color_gradient(low = "#4EA3F1", high = "#E45B5B", name = "-log10 adjusted p value") +
  ggplot2::scale_size_continuous(name = "Count", range = c(2, 9)) +
  ggplot2::labs(title = paste0(out_name, " dotplot (significant: padj < ", enrichment_padj_threshold, ")"), x = "Gene ratio", y = "Enriched term") +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 22, color = TEXT_COLOR),
    axis.title = ggplot2::element_text(size = 17, color = TEXT_COLOR),
    axis.text = ggplot2::element_text(size = 12, color = TEXT_COLOR)
  )
save_gg(dot_plot, paste0(out_name, "_dotplot"), height = max(6, 0.35 * nrow(plot_df) + 2.5))

bar_df <- head(plot_df, 15)
bar_df$term_wrapped <- factor(wrap_labels(bar_df$Description, 45), levels = rev(wrap_labels(bar_df$Description, 45)))
bar_color <- PALETTE[[database]]
bar_plot <- ggplot2::ggplot(bar_df, ggplot2::aes(term_wrapped, Count)) +
  ggplot2::geom_col(fill = bar_color, alpha = 0.88, width = 0.68) +
  ggplot2::coord_flip() +
  ggplot2::labs(title = paste0(out_name, " top terms (significant: padj < ", enrichment_padj_threshold, ")"), x = "Enriched term", y = "Gene count") +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 22, color = TEXT_COLOR),
    axis.title = ggplot2::element_text(size = 17, color = TEXT_COLOR),
    axis.text = ggplot2::element_text(size = 12, color = TEXT_COLOR)
  )
save_gg(bar_plot, paste0(out_name, "_barplot"), height = max(5.5, 0.4 * nrow(bar_df) + 2.2))

## ---- 5. Reports --------------------------------------------------------------

sig_tab <- tab[tab$significant, , drop = FALSE]
sig_tab <- sig_tab[order(sig_tab$p.adjust), , drop = FALSE]
top_source <- if (nrow(sig_tab) > 0) sig_tab else tab[order(tab$p.adjust), , drop = FALSE]
top_n <- min(5, nrow(top_source))
top_label <- if (nrow(sig_tab) > 0) {
  "Top terms by adjusted p value:"
} else {
  sprintf("No term reached padj < %s; top terms by padj are listed for reference only (not claimed as significant):", enrichment_padj_threshold)
}
top_lines <- if (top_n > 0) {
  vapply(seq_len(top_n), function(i) {
    r <- top_source[i, ]
    sprintf("  - %s - padj=%.3g, count=%d, GeneRatio=%s", r$Description, r$p.adjust, r$Count, r$GeneRatio)
  }, character(1))
} else character(0)

small_overlap <- sig_tab[!is.na(sig_tab$Count) & sig_tab$Count < 3, , drop = FALSE]
redundancy_flags <- jaccard_redundancy_flags(tab, padj_threshold = enrichment_padj_threshold)

caveats <- character(0)
for (label in names(mapping_rates)) {
  info <- mapping_rates[[label]]
  if (info$rate < 0.9) {
    caveats <- c(caveats, sprintf(
      "ID mapping rate for %s genes is low: %d/%d (%.1f%%) mapped to ENTREZID. A low mapping rate can silently shrink both foreground and background and bias results (references/pitfalls.md #4); see results/unmapped_genes.tsv.",
      label, info$mapped, info$total, 100 * info$rate
    ))
  }
}
if (nrow(small_overlap) > 0) {
  caveats <- c(caveats, sprintf(
    "%d significant term(s) have a gene overlap count below 3, which is statistically nominal rather than a coherent biological signal (references/pitfalls.md #6): %s.",
    nrow(small_overlap), paste0(small_overlap$Description, " (n=", small_overlap$Count, ")", collapse = ", ")
  ))
}
if (length(redundancy_flags) > 0) {
  examples <- vapply(head(redundancy_flags, 3), function(f) sprintf("'%s' vs '%s' share %d genes (Jaccard=%.2f)", f$a, f$b, f$n, f$jac), character(1))
  caveats <- c(caveats, sprintf(
    "%d pair(s) of top significant terms show high gene-set overlap, suggesting redundant (parent/child or near-duplicate) biology rather than independent findings (references/pitfalls.md #3): %s.%s",
    length(redundancy_flags), paste(examples, collapse = "; "),
    if (simplified_written) " A simplify()-reduced table was also written; prefer it for reporting non-redundant terms." else " Consider clusterProfiler::simplify()/REVIGO-style reduction (GOSemSim was unavailable here) before reporting each term as an independent hit."
  ))
}

conclusion <- if (nrow(sig_tab) == 0) {
  sprintf(
    "None of the %d tested terms reached padj < %s; the terms listed above are the closest to significance but should not be reported as enriched findings without independent support.",
    nrow(tab), enrichment_padj_threshold
  )
} else {
  top1 <- top_source[1, ]
  sprintf(
    "%d of %d tested terms were significant at padj < %s, led by '%s' (padj=%.3g, %d overlapping genes out of a background universe of %d). %s",
    nrow(sig_tab), nrow(tab), enrichment_padj_threshold, top1$Description, top1$p.adjust, top1$Count, length(universe_ids),
    if (length(caveats) > 0) "See caveats above before treating every listed term as an independent biological finding." else "No redundancy or small-overlap issues were flagged among the top terms."
  )
}

report_lines <- c(
  "# Enrichment report", "",
  "## Input", "",
  paste0("- Route: ORA (", ifelse(database == "KEGG", "clusterProfiler::enrichKEGG", paste0("clusterProfiler::enrichGO, ontology=", database)), ")"),
  paste0("- Input mode: ", input_mode),
  paste0("- Foreground genes tested: ", length(foreground_ids)),
  paste0("- Background universe size: ", length(universe_ids), " (source: ", universe_source, ")"),
  paste0("- padj cutoff used to define foreground from the DE table (table mode only): ", padj_cutoff),
  paste0("- Enrichment significance threshold (which tested terms are flagged significant here): padj < ", enrichment_padj_threshold),
  "",
  "## Results", "",
  paste0("- Terms tested (full table, before significance filtering): ", nrow(tab)),
  paste0("- Terms significant at padj < ", enrichment_padj_threshold, ": ", sum(tab$significant)),
  paste0("- GO redundancy reduction (simplify) applied: ", simplified_written),
  paste0("- Unmapped genes during ID conversion: ", if (nrow(unmapped) > 0) paste0(nrow(unmapped), " (see results/unmapped_genes.tsv)") else "0 (no ID mapping step was needed for this keytype/database)"),
  paste0("- ", top_label),
  top_lines,
  if (length(caveats) > 0) c("", "## Caveats", "", paste0("- ", caveats)) else character(0),
  "",
  "## Conclusion", "",
  conclusion
)
writeLines(report_lines, file.path(out_dir, "report.md"))

write.table(
  data.frame(item = c("foreground_genes", "universe_size", "terms_tested", "terms_significant", "enrichment_padj_threshold"),
             value = c(length(foreground_ids), length(universe_ids), nrow(tab), sum(tab$significant), enrichment_padj_threshold)),
  file.path(out_dir, "summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE
)
