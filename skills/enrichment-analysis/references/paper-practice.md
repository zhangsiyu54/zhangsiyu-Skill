# Paper Practice

Follow GSEA/fgsea/clusterProfiler papers: separate gene selection, ranked-list testing, database choice, leading-edge interpretation, and visualization.

## Reference Sources

- GSEA PNAS: https://www.pnas.org/doi/10.1073/pnas.0506580102
- fgsea Bioconductor: https://bioconductor.org/packages/fgsea/
- clusterProfiler Bioconductor: https://bioconductor.org/packages/clusterProfiler/
- enrichplot Bioconductor: https://bioconductor.org/packages/enrichplot/
- MSigDB/GSEA: https://www.gsea-msigdb.org/gsea/index.jsp

## Paper-style Flow

1. Confirm species, gene ID type, DEG thresholds, ranked metric, and database.
2. For ORA, split up/down genes when biologically meaningful.
3. For GSEA, use a full ranked list rather than only significant genes.
4. Report table fields: ID, description, gene ratio, count, p value, adjusted p value, leading-edge/core genes when applicable.
5. Provide database-specific summaries and avoid overinterpreting redundant terms.

## Visualization Pattern

- Include bubble/bar plots for overview and GSEA running score for selected pathways.
- Keep at most 20 terms per panel; wrap long terms.
- Use enrichment map/cnet-style plots only for a small curated subset.
