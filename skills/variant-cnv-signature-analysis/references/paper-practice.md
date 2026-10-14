# Paper Practice

Follow maftools/SigProfiler/inferCNV-style cancer genomics papers: summarize cohort-level mutation structure before interpretation.

## Reference Sources

- maftools GitHub: https://github.com/PoisonAlien/maftools
- maftools Genome Research: https://genome.cshlp.org/content/28/11/1747
- SigProfilerExtractor GitHub: https://github.com/AlexandrovLab/SigProfilerExtractor
- inferCNV GitHub: https://github.com/broadinstitute/infercnv

## Paper-style Flow

1. Confirm MAF-like columns, genome build, cohort, sample IDs, CNV score definition, and signature source.
2. Summarize mutation burden, top mutated genes, variant classifications, and sample-level event counts.
3. Generate mutation frequency/oncoplot-style summaries before gene-specific interpretation.
4. Add CNV heatmaps and signature exposure only when inputs are available and definitions are clear.
5. Avoid claiming driver status from frequency alone.

## Visualization Pattern

- Mutation summary: genes sorted by frequency and samples ordered by burden/group.
- CNV heatmap: blue-white-red centered at neutral copy number.
- Signature exposure: stacked bars with legend outside and limited signature count.
