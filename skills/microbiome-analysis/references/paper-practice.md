# Paper Practice

Follow QIIME 2/DADA2/ANCOM-BC-style microbiome papers: compositional constraints and metadata design must be explicit.

## Reference Sources

- QIIME 2: https://qiime2.org/
- QIIME 2 Nature Biotechnology: https://www.nature.com/articles/s41587-019-0209-9
- DADA2: https://benjjneb.github.io/dada2/
- DADA2 Nature Methods: https://www.nature.com/articles/nmeth.3869
- ANCOM-BC Bioconductor: https://bioconductor.org/packages/ANCOMBC/

## Paper-style Flow

1. Identify raw reads, ASV/OTU table, taxonomy, metadata, or existing diversity outputs.
2. Confirm sample matching, taxonomy rank, sequencing region, and whether counts or relative abundances are provided.
3. For raw reads, follow DADA2/QIIME2 denoise, taxonomy, diversity, and differential abundance routes.
4. For tables, report alpha diversity, beta diversity/ordination, composition, PERMANOVA when possible, and differential abundance with compositional caveats.
5. Document filtering, rarefaction/normalization, distance metric, and model formula.

## Visualization Pattern

- PCoA/PCA: legend outside, axis labels include variance when available.
- Composition: top taxa plus Other; avoid more than 10 taxa.
- Differential taxa: horizontal bars or volcano-style plots with readable taxonomy labels.
