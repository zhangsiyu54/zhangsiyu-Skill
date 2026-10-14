# Enrichment Advanced Modules

Load this file for pathway interpretation beyond a simple bubble plot.

## ORA vs GSEA

- ORA requires a thresholded gene list and a background universe.
- GSEA/fgsea requires a ranked named vector, not only significant genes.
- Do not mix ORA and GSEA interpretation without labeling the method.

## Databases

- GO BP/CC/MF: broad ontology terms, often redundant.
- KEGG/Reactome: pathway mechanisms, better for interpretable pathway names.
- MSigDB Hallmark: compact gene sets for high-level biological programs.
- Custom GMT: use when user has project-specific signatures.

## ID Mapping

- Record organism, input ID type, mapped ID type, mapped genes and unmapped genes.
- Do not guess human vs mouse or SYMBOL vs ENSEMBL.
- Keep an unmapped gene table.

## Interpretation

- Report leading-edge genes for GSEA.
- Collapse redundant GO terms only when a semantic similarity tool is available.
- Avoid overclaiming: enrichment means association of a gene set with the input list/rank, not direct pathway activation.

## Deliverables

- ORA: enrichment table, up/down separate results, bubble/bar plots.
- GSEA: ranked metric file, NES/padj table, leading-edge genes, running score plots when available.
- Summary: top terms grouped by database and direction.
