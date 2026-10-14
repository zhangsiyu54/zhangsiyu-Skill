# Enrichment Analysis Pitfalls

Read this file before finalizing any enrichment conclusion. Each entry: Trigger (when it happens), Mechanism (why it happens), Symptom (how it shows up in output), Fix (what to do instead).

## 1. Whole-genome background in ORA

- Trigger: running `enrichGO`/`enrichKEGG`-style hypergeometric tests with the default "all annotated genes in the organism package" as background instead of the genes actually tested/detected in the experiment.
- Mechanism: the hypergeometric/Fisher test compares the DEG list against a background universe; if the universe includes tens of thousands of genes never expressed or never testable in the assay (e.g., filtered out by `filterByExpr` or low counts), the apparent enrichment odds ratio is inflated because the "non-DEG" pool is artificially large and unrepresentative.
- Symptom: dozens of GO/KEGG terms reach significance with implausibly small p values; terms related to ubiquitously expressed housekeeping processes appear enriched; results do not replicate when background is corrected.
- Fix: set `universe` to the full list of genes that passed expression filtering / were tested for differential expression (i.e., all rows in the DE results table, not only significant ones), not `NULL`/whole-genome default. Report the universe size in methods text.

## 2. GSEA ranked on a noise-dominated metric

- Trigger: ranking genes for `fgsea`/GSEA by raw (unshrunken) log2FC, or by log2FC alone without incorporating variance/expression-level confidence.
- Mechanism: unshrunken log2FC estimates for low-count genes are extremely noisy (large fold changes from small, high-variance counts) and can dominate the extremes of the ranking, which fgsea's running-sum statistic weights most heavily; this pulls "significant" gene sets toward noise-driven genes rather than true signal.
- Symptom: leading-edge genes for top pathways are unfamiliar low-expression genes with huge but statistically uncertain fold changes; NES/padj is unstable across re-runs or small pipeline changes; results don't align with the DEG table's adjusted-p ranking.
- Fix: rank by a signed test statistic (e.g., DESeq2 Wald stat, limma `t` or `-log10(p)*sign(logFC)`) or shrunken log2FC (`lfcShrink`/apeglm), not raw log2FC. Document the exact ranking formula in methods.

## 3. Redundant GO terms reported as independent findings

- Trigger: reporting a long list of "significant" GO BP terms without checking for parent-child or highly overlapping-gene-set relationships (e.g., "immune response" and "regulation of immune response" both appear, driven by the same 15 genes).
- Mechanism: GO's hierarchical, non-independent structure means many terms share most of their annotated genes; standard ORA/GSEA tests each term independently and does not account for this structural non-independence, so semantically redundant terms all pass significance together.
- Symptom: top-20 term list reads like near-duplicates (parent/child pairs, synonyms); a bubble plot shows many bubbles clustering around the same handful of driver genes.
- Fix: apply semantic-similarity-based reduction before reporting/plotting — `clusterProfiler::simplify()` (GOSemSim-based) or an external REVIGO pass — and state the similarity cutoff used. Do not present un-reduced GO BP lists as if each term were an independent biological finding.

## 4. Silent gene loss during ID mapping

- Trigger: converting between ENTREZ/ENSEMBL/SYMBOL (e.g., `bitr()`, `AnnotationDbi::select()`) without checking the mapping's success rate.
- Mechanism: ID mapping is not 1:1 — one Entrez ID can map to multiple symbols (aliases) and one symbol can map to multiple Entrez IDs (paralogs, retired/merged IDs); default mapping functions silently drop unmapped IDs and/or arbitrarily collapse many-to-one duplicates by keeping the first match.
- Symptom: input gene count and post-mapping gene count differ but the drop is never reported; duplicate rows for the same gene after mapping, sometimes with conflicting expression values; a gene the collaborator expects to see is missing with no explanation.
- Fix: always compute and report the mapped/unmapped counts and the mapping rate (e.g., "18,412 of 19,850 input genes mapped, 93%"); write an `unmapped_genes.tsv`; for many-to-one collisions, use an explicit resolution rule (e.g., keep the row with the largest absolute statistic) and state it in methods rather than relying on function defaults.

## 5. Multiple-testing correction applied per-database instead of across the full test space

- Trigger: running BH/FDR correction separately within GO BP, then separately within KEGG, then separately within Reactome, and reporting each `padj` as if it controls a single overall false discovery rate.
- Mechanism: GO alone can involve testing several thousand terms simultaneously; correcting each database's p values independently understates the true multiplicity of the combined search (test all databases together, or at minimum, always disclose that corrections are per-database) and can let borderline terms across multiple small databases collectively slip through.
- Symptom: a term is "significant" in one small custom/database collection with few tested terms but would not survive correction if pooled with the thousands of GO terms tested in the same analysis; reviewers ask "corrected across how many tests?" and the report has no clear answer.
- Fix: report, for every enrichment table, the exact number of terms/gene sets tested that FDR correction was applied over, and keep per-database correction only when explicitly comparing databases separately (label this clearly). Never imply a single global FDR when correction was actually done per collection.

## 6. Leading-edge overlap too small to be biologically meaningful

- Trigger: highlighting a GSEA/ORA result as a key finding when the leading-edge/core-enrichment or overlapping-gene count is 2-3 genes, often out of a gene set with 100+ members.
- Mechanism: with small overlaps, the statistical signal is driven almost entirely by 1-2 highly ranked genes; NES/p value can still be nominally significant (especially with permutation-based p values at low N) but does not reflect coordinated regulation of the pathway.
- Symptom: a "significant" pathway with count = 2-3 genes and gene ratio well under 5%; the same 2-3 genes recur as leading-edge across many unrelated-sounding pathways (a hub-gene artifact).
- Fix: apply and report a minimum leading-edge/overlap size (e.g., >=5 genes, or a gene ratio floor) before calling a term biologically interpretable; flag terms below threshold as "statistically nominal, not emphasized" rather than omitting the number silently.

## 7. Species mismatch and naive ortholog matching

- Trigger: applying human MSigDB/Hallmark gene sets to mouse (or other non-human) data by simply uppercasing mouse gene symbols and joining on string match, instead of using a curated ortholog mapping.
- Mechanism: uppercase-symbol matching works for many 1:1 orthologs by coincidence (mouse symbols are typically the human symbol capitalized) but silently fails or mismatches for genes with non-trivial ortholog relationships (many-to-many paralog families, species-specific duplications/losses, differently named orthologs), and it cannot detect when no true ortholog exists.
- Symptom: gene set sizes after "conversion" look suspiciously close to the original human set size (near 100% apparent match) with no reported failures; a handful of pathway hits trace back to a false-positive symbol coincidence rather than a true ortholog.
- Fix: use a maintained ortholog mapping resource (e.g., `orthogene`, `babelgene`, MGI/HGNC homology tables, biomaRt) and report conversion success rate and how many-to-many cases were resolved. State the organism and source-organism-of-gene-sets explicitly in methods; never assume symbol case conversion is sufficient.

## 8. Treating KEGG "disease" or generic reference maps as mechanistic evidence

- Trigger: reporting a hit on a KEGG map such as "Pathways in cancer," "Human papillomavirus infection," or another broad disease/reference map as if it demonstrates a specific mechanism relevant to the study.
- Mechanism: KEGG's disease and several reference maps are large curated aggregates spanning many unrelated sub-mechanisms and cell types; a handful of overlapping genes (e.g., common signaling genes like PIK3CA, TP53) will hit many such broad maps regardless of the specific biology, because those genes participate in dozens of curated diagrams.
- Symptom: a generic hub map dominates the top results across many unrelated projects/datasets; the gene overlap driving the hit is a small set of promiscuous, broadly-annotated genes rather than a disease- or tissue-specific signature.
- Fix: inspect which genes drive the overlap on broad/disease maps before interpreting; prefer more specific KEGG/Reactome sub-pathways or Hallmark gene sets for mechanistic claims, and explicitly caveat disease-category hits as descriptive/associative rather than mechanistic in the report text.

## 9. Conflating ORA and GSEA as answering the same question

- Trigger: presenting ORA and GSEA results for the same comparison in a single merged table/figure without distinguishing what each tests, or using them interchangeably to "confirm" each other.
- Mechanism: ORA tests whether a gene set is over-represented among a fixed, thresholded subset of genes (binary in/out, hypergeometric test) — it discards magnitude/rank information for genes near but outside the cutoff. GSEA tests whether a gene set is skewed toward one end of a full ranked continuum, using every gene's rank/magnitude, and does not require a significance threshold. Agreement or disagreement between the two does not have a single, simple statistical interpretation.
- Symptom: report claims "ORA and GSEA both confirm pathway X is enriched" as though this doubles the evidence, or claims a discrepancy invalidates one method, without acknowledging they test different hypotheses on different information.
- Fix: always label ORA and GSEA results separately with their respective hypotheses stated in one sentence each; when both are run, describe convergent results as "consistent under two different enrichment frameworks" rather than as independent replication, and explain discordance as expected given the different null hypotheses rather than as an error.

## 10. Circular enrichment (testing a pathway-derived gene list against its source database)

- Trigger: deriving a gene list from a pathway/database-driven step (e.g., genes annotated to "immune response" pulled from GO itself, or a gene panel curated from KEGG pathway membership) and then re-testing that list for GO/KEGG enrichment against the same database.
- Mechanism: if the input gene list's construction already used the annotation database as a filter, re-testing against that same database is tautological — the enrichment result reflects the selection criterion used to build the list, not an independent biological signal from the experiment.
- Symptom: the single most "significant" term is the exact category used to build/filter the input gene list, or a set of near-synonymous parent/child terms of it, with an extreme enrichment p value (often the most significant term in the whole table).
- Fix: trace how the input gene list was constructed before enrichment; if it was filtered/curated using GO/KEGG/Reactome membership, either exclude that database from re-testing or explicitly disclose the circularity and interpret only the remaining, independently-selected gene sets.

## 11. Gene-length / detection-power bias in ORA

- Trigger: RNA-seq-derived DEG lists tested with standard ORA without acknowledging that statistical power to call a gene "differentially expressed" correlates with gene length and read count, not biology.
- Mechanism: longer transcripts and more highly expressed genes accumulate more reads and thus more statistical power to reach significance at a given effect size; GO/KEGG categories enriched for long genes (e.g., signaling receptors, extracellular matrix, kinases) can appear enriched purely from this detection-power artifact (documented by Timmons et al. 2015; Young et al. goseq method addresses this for GO).
- Symptom: recurring enrichment for "signal transduction," "cell adhesion," or other long-gene-rich categories across unrelated datasets/tissues; enrichment persists even in permuted/label-shuffled control DEG lists.
- Fix: when gene-length bias is a concern (RNA-seq ORA on GO terms), consider length-aware methods (e.g., `goseq`) or at minimum flag the possibility in limitations text; prefer GSEA on a full ranked list (less sensitive to a hard length-correlated threshold) when this bias is a known concern for the dataset.

## 12. Ignoring the "background problem" — background silently narrowed by the tool

- Trigger: some ORA tools/web servers automatically restrict the background universe to only genes annotated in at least one category of the tested database, silently dropping unannotated genes from both the target list and the background.
- Mechanism: if annotation coverage differs between the significant list and the background (e.g., a batch of newly discovered/annotated genes concentrated in the DEG list), auto-restricting the background changes the effective denominator asymmetrically, biasing the test (described as the "background problem" in Ziemann et al. 2024).
- Symptom: enrichment results change substantially when the same DEG list and background are re-run in a different tool defaulting to a different background-restriction behavior; total gene counts reported by the tool don't match the counts supplied.
- Fix: confirm what background the specific tool/function actually uses internally (read the documentation, not just the argument name), and report the effective universe size used in the test, not just the size supplied.

## 13. Unstable or over-precise p values from too few GSEA permutations

- Trigger: reporting fgsea/GSEA p values with many digits of precision (e.g., p = 1e-12) from a default `nperm`/permutation-based run, or comparing padj across pathways where the nominal p is already at the resolution floor of the permutation count.
- Mechanism: permutation-based p values have a resolution floor set by the number of permutations (e.g., 1,000 permutations cannot resolve p below ~1e-3 reliably without extrapolation); `fgsea`'s adaptive multi-level splitting (`fgseaMultilevel`) approximates smaller p values but very small p values still carry wide relative uncertainty unless `eps` and sample depth are set deliberately.
- Symptom: two runs of the same GSEA with different random seeds give visibly different p values/rankings for top pathways; extremely small p values reported without a corresponding permutation count in methods.
- Fix: use `fgseaMultilevel` (or increase `nperm`) for stable estimates of small p values, set and report a fixed seed for reproducibility, and state the eps/permutation settings in the methods text; treat p values near the resolution floor as approximate.

## 14. Database/version drift breaking reproducibility

- Trigger: re-running or comparing an enrichment analysis months/years apart, or across collaborators, without pinning GO/KEGG/Reactome/MSigDB release versions.
- Mechanism: GO, KEGG, Reactome, and MSigDB are updated on independent schedules (GO terms are added/obsoleted, KEGG pathways are restructured, MSigDB adds/removes gene set members between releases); an analysis re-run against a newer database version can yield different gene set memberships and therefore different enrichment results for the identical input gene list.
- Symptom: a "replication" analysis on the same DEG list months later shows different top terms with no change to the input data; a collaborator using a different `org.Hs.eg.db`/`msigdbr` version cannot reproduce the reported numbers.
- Fix: record and report the exact database/package version (e.g., `org.Hs.eg.db` build date, MSigDB release number such as "v2023.2.Hs", KEGG release date) in methods; when reproducibility across time matters, archive the GMT/annotation files used alongside the results.
