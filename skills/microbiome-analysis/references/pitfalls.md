# Microbiome Analysis Pitfalls

Read this before finalizing any conclusion, figure caption, or methods paragraph. Each entry lists the trigger condition, the underlying mechanism, how it typically shows up in results, and the fix. When a pitfall has genuine scientific debate, both sides are noted so conclusions are not overstated in either direction.

## 1. Rarefaction discards real data and can be statistically inefficient

- **Trigger**: Subsampling every sample down to the minimum library size before alpha/beta diversity or differential abundance.
- **Mechanism**: Rarefying throws away sequences from every sample larger than the smallest one, discarding real counts and adding subsampling noise. McMurdie & Holmes (2014, PLOS Computational Biology) argue this is inefficient compared with mixture-model-based normalization (e.g., variance-stabilizing or negative-binomial approaches from RNA-seq).
- **Symptom**: Reduced power to detect true differences, especially when library sizes vary a lot across samples; results that shift when the rarefaction depth or random seed changes.
- **Fix**: For differential abundance, prefer count-aware compositional methods (ANCOM-BC, ALDEx2, DESeq2/edgeR-style models) that use the full count matrix with a library-size offset instead of subsampling. Reserve rarefaction for alpha/beta diversity comparisons where equal-depth sampling is the simplest way to remove a depth artifact, and always report the depth chosen and how many samples/reads were dropped.
- **Counterpoint to keep in mind**: A 2023 reanalysis by Schloss (mSphere) found methodological issues in the original McMurdie & Holmes simulations and concluded rarefaction is a robust way to control uneven sequencing effort for many practical cases. Treat "never rarefy" as a strong default, not an absolute rule — state the choice and justification explicitly rather than picking a side silently.

## 2. Applying RNA-seq-style normalization without addressing compositionality

- **Trigger**: Running DESeq2 or edgeR directly on an ASV/OTU count table the same way it would be run on RNA-seq gene counts.
- **Mechanism**: Microbiome count data are compositional — each sample's total is an arbitrary artifact of sequencing depth, not absolute biomass (Gloor et al. 2017, *Frontiers in Microbiology*). A real increase in one taxon forces the apparent relative abundance of every other taxon down even if their absolute abundance is unchanged, so naive size-factor normalization can flip the sign of a fold-change for taxa that did not actually change.
- **Symptom**: A taxon reported as "significantly decreased" when only a different, dominant taxon increased; differential-abundance direction that is not reproducible when a different reference/denominator taxon is used.
- **Fix**: Use a log-ratio approach (CLR/ALR) or a compositionally aware method (ANCOM, ANCOM-BC, ALDEx2) that estimates a sampling-fraction/bias-correction term rather than treating the count sum as informative. State explicitly in methods that abundances are relative, not absolute.

## 3. Alpha diversity compared across very different sequencing depths

- **Trigger**: Comparing Shannon, Simpson, Chao1, or Observed Features across samples/groups without checking read-depth distribution.
- **Mechanism**: Richness-sensitive metrics (Observed Features, Chao1) and, to a lesser extent, Shannon/Simpson are systematically higher in deeper-sequenced samples because more reads reveal more rare taxa. If depth correlates with a covariate (e.g., cases sequenced on an older run with fewer reads), the depth artifact is confounded with the biological comparison.
- **Symptom**: A group difference in alpha diversity that disappears (or reverses) after rarefying to equal depth or after plotting diversity against read count and seeing a strong trend.
- **Fix**: Rarefy to a common, justified depth for alpha diversity comparisons (or use depth-corrected estimators / rarefaction curves extrapolated with iNEXT-style methods), and always plot/report diversity vs. sequencing depth per sample to rule out a depth confound before reporting group differences.

## 4. Beta diversity distance choice silently changes conclusions

- **Trigger**: Reporting only one beta diversity distance (usually Bray-Curtis) without justification, or switching metrics until one gives a significant PERMANOVA.
- **Mechanism**: Bray-Curtis and unweighted Jaccard emphasize compositional (often rare-taxon presence/absence) differences; weighted UniFrac emphasizes phylogenetically close, abundant lineages; unweighted UniFrac emphasizes phylogenetic membership regardless of abundance (Lozupone & Knight 2005). Each answers a different biological question, so they can legitimately disagree.
- **Symptom**: Ordination and PERMANOVA results that look convincingly separated with one metric and null with another, with no explanation of why the metric was chosen.
- **Fix**: Choose the distance based on the hypothesis first (phylogenetic conservatism vs. raw compositional/abundance shift vs. rare-taxon membership), state the choice in methods before looking at results, and report more than one metric when the question is exploratory rather than confirmatory.

## 5. PERMANOVA significance driven by dispersion, not centroid location

- **Trigger**: Reporting a significant PERMANOVA (adonis2) result as evidence that group centroids differ without checking group dispersion.
- **Mechanism**: PERMANOVA is sensitive to both differences in multivariate location (centroid) and differences in within-group spread (dispersion). A group with more heterogeneous community composition can produce a significant PERMANOVA even if the average community composition is the same.
- **Symptom**: Significant PERMANOVA paired with an ordination plot where one group's points are simply more scattered rather than shifted; result that does not hold up when checked with betadisper (PERMDISP2).
- **Fix**: Always run `betadisper`/PERMDISP2 alongside PERMANOVA and report both. If dispersion differs significantly, qualify the PERMANOVA interpretation (or use a distance/test less sensitive to dispersion) instead of claiming a compositional shift.

## 6. Batch effects from sequencing run, extraction kit, or low-biomass contamination ("kitome")

- **Trigger**: Samples processed across multiple extraction kit lots, sequencing runs, or plates, especially with low-biomass material (skin, lung, environmental swabs, low-input clinical samples).
- **Mechanism**: DNA extraction kits and PCR reagents carry their own background microbial DNA. In low-biomass samples this contaminating signal can rival or exceed the true biological signal (Salter et al. 2014, *BMC Biology*). Run/kit batch can also be confounded with the group of interest if cases and controls were processed on different days or kits.
- **Symptom**: Clustering by extraction batch or sequencing run rather than by biological group in ordination plots; taxa known as common reagent contaminants (e.g., *Ralstonia*, *Burkholderia*, *Delftia*, *Pseudomonas* in some kits) appearing as "differentially abundant."
- **Fix**: Always sequence extraction blanks/negative controls and use a contaminant-removal tool (e.g., decontam, prevalence or frequency method) before analysis; check whether batch/kit is confounded with the group variable and include it as a covariate or blocking factor if not fully confounded; flag low-biomass studies explicitly as higher risk in the limitations section.

## 7. Taxonomic classification database and version mismatches break cross-study comparability

- **Trigger**: Comparing taxa composition or differential taxa lists across studies, or across time within the same lab, without checking the reference database and version used for classification.
- **Mechanism**: SILVA, Greengenes/Greengenes2, and GTDB use different taxonomic frameworks, naming conventions, and even different phylogenies (GTDB is genome-phylogeny-based and has substantially revised bacterial taxonomy relative to legacy schemes; Parks et al. 2018, *Nature Biotechnology*). The same ASV sequence can be assigned different genus/family names, or a different resolution, depending on the database and version.
- **Symptom**: A taxon named as "enriched" in one study cannot be matched to a taxon in another study; apparent disagreement between two datasets that is actually a naming/database artifact rather than a biological one.
- **Fix**: Always report the classifier, database, and exact version/release in methods; when comparing across studies, re-classify with a common database and version where feasible, or compare at the ASV/sequence level (or a coarser, more stable rank) instead of relying on genus/species labels alone.

## 8. Spurious correlations from unaddressed compositional constraints

- **Trigger**: Computing pairwise Pearson/Spearman correlations between taxa (or between a taxon and a metadata variable) directly on relative abundances.
- **Mechanism**: Because relative abundances sum to a constant, an increase in one taxon mechanically forces others down, inducing negative correlations among unrelated taxa (a well-known compositional data analysis artifact). Standard correlation networks built on raw relative abundance can therefore report co-occurrence/exclusion patterns that are pure normalization artifacts.
- **Symptom**: Dense networks of "negative correlations" that do not replicate across datasets or subsets; correlation strength that scales with the number of taxa included.
- **Fix**: Use log-ratio transforms (CLR, or pairwise log-ratio methods such as SparCC, SPIEC-EASI, or proportionality metrics like phi/rho) designed for compositional data instead of raw Pearson/Spearman on relative abundance.

## 9. ASV vs. OTU clustering choice changes resolution and reproducibility

- **Trigger**: Mixing ASV-based (DADA2/Deblur) and OTU-based (97% similarity clustering) results across a project, or comparing an ASV table against a historical OTU table.
- **Mechanism**: ASV methods resolve single-nucleotide differences and are deterministic/reusable across studies (the same true sequence gets the same ASV in any dataset processed the same way), while OTU clustering groups sequences at a similarity threshold that depends on the reference/de novo algorithm and the specific sample set, so cluster membership is not stable across runs or studies.
- **Symptom**: Apparent loss or gain of "species" richness when switching pipelines that is actually a resolution/method artifact; inability to directly merge or compare ASV and OTU tables.
- **Fix**: Pick ASVs (DADA2/Deblur) as the default for new analyses given their reproducibility, state the method explicitly, and avoid directly pooling ASV and OTU tables from different pipelines without a common re-clustering or sequence-matching step.

## 10. Treating a single differential abundance method's result as definitive

- **Trigger**: Reporting the output of one differential abundance tool (e.g., a single DESeq2 or LEfSe run) as the final list of significant taxa.
- **Mechanism**: A large benchmarking study across 38 real 16S datasets found that different differential abundance methods (edgeR, DESeq2, LEfSe, limma-voom, ALDEx2, ANCOM, corncob, etc.) can produce drastically different sets and numbers of significant taxa on the same data, with false-positive/false-discovery behavior varying widely by method and by dataset characteristics such as sample size and effect size (Nearing et al. 2022, *Nature Communications*). ALDEx2 and ANCOM tended to be the most conservative/consistent in that benchmark.
- **Symptom**: A dramatically different "top hits" list if the differential abundance tool is swapped, with no sensitivity analysis reported.
- **Fix**: Run at least two methodologically distinct approaches (e.g., a compositional method such as ANCOM-BC/ALDEx2 plus a count-based model such as MaAsLin2) and report taxa consistent across methods with more confidence than single-method hits; always report the specific tool, version, and default vs. custom parameters used.

## 11. Structural zeros and pseudocount choice distort log-ratio methods

- **Trigger**: Adding an arbitrary pseudocount (e.g., +1) to a sparse ASV table before log or CLR transformation, or ignoring taxa that are truly absent (structural zero) in one group vs. taxa that are just under the detection limit (sampling zero).
- **Mechanism**: Log-ratio-based methods are undefined at zero, so a pseudocount must be added; the size of that pseudocount can materially change CLR values for low-abundance taxa, especially in sparse, low-depth samples. Structural zeros (a taxon that biologically cannot occur in one group) need different handling than sampling zeros (a taxon that is simply below detection in a particular sample).
- **Symptom**: Significance or fold-change magnitude for rare taxa that changes noticeably with the pseudocount value (e.g., 0.5 vs. 1 vs. a multiplicative replacement); a taxon flagged as differentially abundant purely because it is present in one group and completely absent in the other, which some methods treat differently than a partial change.
- **Fix**: Prefer methods with principled zero handling (ANCOM-BC's structural zero test, ALDEx2's Bayesian Dirichlet-multinomial replacement) over a hand-picked additive pseudocount; report the prevalence filter and zero-handling approach in methods, and treat "present in one group, absent in the other" taxa as a distinct reporting category rather than folding them into a continuous fold-change.

## 12. Longitudinal/repeated-measures data analyzed as if samples were independent

- **Trigger**: Running standard PERMANOVA, alpha diversity ANOVA, or differential abundance on a dataset with multiple timepoints per subject without accounting for subject identity.
- **Mechanism**: Repeated samples from the same subject are correlated (pseudoreplication), inflating apparent degrees of freedom and often making tests anti-conservative (p-values too small) because the model treats correlated observations as independent evidence.
- **Symptom**: Very small p-values that shrink further simply by adding more timepoints per subject rather than more subjects; effect sizes that do not replicate when the same number of independent subjects is resampled.
- **Fix**: Use subject as a blocking/strata term in PERMANOVA (e.g., `vegan::adonis2` with `strata=`), or use mixed-effects models with a random intercept/slope for subject (MaAsLin2 with a random effect, or `lme4`); state the repeated-measures structure explicitly in methods.

## 13. Relative abundance mistaken for absolute abundance or microbial load

- **Trigger**: Concluding a taxon or the whole community "increased" or "bloomed" based only on relative abundance/sequencing proportions.
- **Mechanism**: Sequencing yields relative, not absolute, abundance; a taxon's relative abundance can rise purely because other taxa declined (or total biomass dropped), with no real increase in that taxon's absolute count, and vice versa.
- **Symptom**: A claimed "bloom" of a taxon that is not corroborated by orthogonal quantification (qPCR, flow cytometry, spike-in), or that disappears when a total-abundance-aware analysis is applied.
- **Fix**: Where the question concerns absolute load (e.g., "did total *Bacteroides* increase"), pair sequencing with an absolute quantification method (16S qPCR, spike-in standards, cell counting) or explicitly caveat that all reported changes are relative-abundance shifts, not absolute abundance changes.

## 14. Primer/region choice and cross-region merging bias community profiles

- **Trigger**: Merging or directly comparing ASV/OTU tables generated from different 16S hypervariable regions (e.g., V3-V4 vs. V4 vs. V1-V2), or comparing results across studies without checking the primer set.
- **Mechanism**: Different primer pairs have different amplification efficiency and taxonomic coverage biases (some primers under-amplify or miss certain phyla/genera entirely), and DADA2 error models and ASV sequences are region-specific, so ASVs from different regions are not directly comparable sequences even for the same organism.
- **Symptom**: A taxon appearing "absent" in one dataset and "present" in another purely because of primer coverage differences, not true biological absence; failed ASV merging or inflated apparent diversity when tables from different regions are concatenated.
- **Fix**: Always record and report the primer set/region in methods; do not merge ASV tables across regions — instead, compare at a coarser taxonomic rank after independent classification, or restrict cross-region comparisons to closed-reference OTU/exact-match approaches designed for that purpose.
