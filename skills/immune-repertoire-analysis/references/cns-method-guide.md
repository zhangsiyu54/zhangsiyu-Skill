# CNS Method Guide

Use this guide before TCR/BCR repertoire, clonotype, or scVDJ-coupled analysis.

## Minimal Inputs

- Require clonotype/contig table, sample IDs, chain type, clonotype ID, CDR3 sequence, V/D/J gene columns, count/frequency, group labels, and cell metadata for scVDJ.
- Confirm whether clonotype is defined by exact CDR3 nucleotide, CDR3 amino acid, paired alpha-beta/heavy-light chain, or tool-specific clonotype ID.
- Keep productive/nonproductive chains and multi-chain cells explicit.

## QC and Harmonization

- Summarize cells/contigs per sample, productive contig fraction, chain recovery, clonotype counts, singleton fraction, and sample-level sequencing depth.
- Normalize diversity and expansion comparisons for sampling depth when possible.
- Do not compare raw clonotype counts across samples with very different cell numbers without frequency or rarefaction-aware summaries.

## Core Analyses

- Clonal expansion: report clone size classes, top clonotypes, frequency, and group/cell-type distribution.
- Diversity: report Shannon, Simpson/Gini, Chao/observed richness when appropriate and state whether downsampling was used.
- V/J usage: compare usage by group with sample-level awareness.
- Overlap: use Morisita/Jaccard-style metrics and clarify whether overlap is sequence-level or clonotype-ID-level.
- scTCR/scBCR coupling: link clonotypes to cell types, clusters, pseudotime, phenotype, or marker states only after barcode matching.

## Interpretation Boundaries

- Clonotype expansion suggests antigen-driven or proliferative processes but does not identify antigen specificity without validation.
- Public clones and convergent motifs require careful cohort/sample context.
- BCR lineage/SHM analysis requires immunoglobulin-specific preprocessing and germline references.

## Algorithm and Design Decision Tree

Read `references/pitfalls.md` alongside this table before finalizing method choice; the "Rationale" column names the pitfall each choice avoids.

| Scenario | Tool/Method/Parameter | Rationale |
|---|---|---|
| Diversity comparison across uneven cell numbers/depth | Downsample all samples to common cell/template count first (repeat over multiple random subsamples); report rarefaction curves; use Chao1 (bias-corrected) or Recon/DivE for extrapolated richness with error bars | Raw richness and clonality scale with sampling effort; unequal depth alone can flip rank order (pitfall 1, 4) |
| BCR clonal lineage / somatic hypermutation analysis | Do not exact-match CDR3; cluster by same V+J gene, same CDR3 length, sequence similarity threshold (Change-O `DefineClones`/SCOPer/partis), then build lineage trees with IgPhyML or dowser using TIgGER-corrected germline alleles | SHM causes descendants of one clone to diverge in sequence; exact-match clonotyping undercounts BCR clonal expansion (pitfall 7) |
| Public clonotype overlap testing | Compute observed overlap (Morisita-Horn or Jaccard) and compare to a permutation-based or generation-probability (OLGA/IGoR pgen) null model before calling overlap biologically meaningful | Convergent recombination produces nonzero background overlap between unrelated repertoires; raw shared-clone counts alone cannot distinguish signal from chance (pitfall 6) |
| Bulk vs. single-cell repertoire data handling | For bulk (amplicon or RNA-seq derived, e.g. MiXCR/TRUST4 assemblies): treat frequency as read/UMI-fraction, disclose whether UMI-corrected; for single-cell (10x VDJ): treat frequency as cell-fraction, apply doublet/multiplet filtering, and do not merge bulk read-fraction and single-cell cell-fraction into one comparable axis | The two data types estimate clone size from fundamentally different units (molecules vs. cells); direct numeric comparison without unit conversion is invalid |
| Convergent recombination-aware clonotype definition | Default to CDR3-amino-acid + V/J gene as the "clonotype" unit for cross-sample comparison; when expansion claims matter, verify top clones share consistent V/J and low nucleotide diversity; use CDR3-nucleotide only when studying receptor generation/repertoire-formation questions | aa-only clonotyping can merge independently-generated convergent sequences into a false "clone"; nt-only clonotyping can split one true clone into artificial singletons (pitfall 5) |
| scTCR/scBCR multiplet suspicion (3+ chains per barcode) | Apply pipeline doublet test (e.g., Cell Ranger VDJ doublet filter) before analysis; report excluded fraction; for retained multi-chain cells, collapse to the dominant chain pair by UMI count and state this rule explicitly | Multiplets/ambient contamination fabricate spurious multi-chain clonotypes; not all multi-chain cells are artifacts (allelic inclusion is real, ~10-20% baseline), so blanket removal without reporting the rule is also wrong (pitfall 3) |
| V/J usage comparison across studies/cohorts | Restrict comparison to samples processed with the same library prep chemistry (same multiplex primer set, or same 5'-RACE/UMI protocol); flag cross-protocol comparisons as descriptive only | Multiplex PCR primer sets and RACE/UMI protocols have unequal, gene-specific amplification efficiency that is indistinguishable from biological usage bias unless protocol is held constant (pitfall 8) |
| Small-N samples (few hundred cells or fewer) | Compute bootstrap or jackknife confidence intervals for every diversity/overlap metric; avoid ranking or hypothesis-testing samples whose CIs overlap; state cell count next to every number | Point diversity estimates from small samples carry high sampling variance and give false precision if reported without uncertainty (pitfall 9) |
| Clonally expanded cells in per-cell differential tests (scRNA+VDJ) | Aggregate by clonotype (pseudobulk-by-clone) or use mixed-effects models with clonotype as a random effect rather than naive per-cell tests; report clone size distribution alongside any significance claim | Clonally related cells are not independent observations; a single large clone can dominate a per-cell test and masquerade as a group effect (pitfall 12) |
| Longitudinal/multi-tissue clone tracking (alluvial/Sankey) | Limit to top clonotypes with adequate cell support per node; test shared-clone counts against a permutation null before using causal language ("migrated," "differentiated") | Co-occurrence of a clonotype label across timepoints/tissues does not by itself establish a causal or temporally ordered relationship (pitfall 10) |

## Required Figures

- Clonal expansion stacked/bar plot, top clonotype frequency plot, diversity boxplot, V/J usage heatmap/bar, overlap heatmap, alluvial/sankey clone tracking, UMAP colored by clone size for scVDJ.

## Methods Text Must Include

- Tool/source, clonotype definition, chain handling, productive filtering, sample normalization, diversity metrics, overlap metric, statistical tests, thresholds, and package versions.
