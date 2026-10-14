# Enrichment Input Contract

Accepted inputs:
- `differential_expression.tsv` with `gene`, `log2FC`, `padj`.
- `up_genes.txt` and/or `down_genes.txt`.
- ranked table with `gene` and ranking metric.
- existing enrichment table with `term`, `category`, `gene_ratio`, `padj`, `count`.

Required parameters for new enrichment:
- organism.
- gene ID type.
- databases: GO BP/CC/MF, KEGG, Reactome, GSEA.
- separate up/down or combined analysis.

Stop and ask if organism or ID type is missing.
