# Microbiome Advanced Modules

## Raw Read Processing

- DADA2/QIIME2 route: primer trimming, quality filtering, denoising, chimera removal, ASV table, taxonomy assignment.
- Record truncation lengths, database, classifier, retained reads and chimera rate.

## Diversity

- Alpha diversity: Observed features, Shannon, Simpson or user-specified metrics.
- Beta diversity: Bray-Curtis, Jaccard, UniFrac when phylogeny is available.
- PERMANOVA: report distance, permutations, formula and R2.

## Differential Abundance

- Prefer ANCOM-BC/ALDEx2/DESeq2-style count-aware methods over simple relative-abundance t tests.
- Record prevalence and abundance filtering.
- State compositional bias and multiple testing correction.

## Functional Profiles

- HUMAnN/KO/MetaCyc profiles can be plotted like pathway abundance.
- Do not mix taxonomic and functional conclusions without stating the input type.
