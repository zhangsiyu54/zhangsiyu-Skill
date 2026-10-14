# bulkRNA Input Contract

Required files:
- `count_matrix.csv` or `counts.tsv`: first column `gene`, remaining columns sample IDs.
- `metadata.csv`: one row per sample. Must include `sample` and a group column.

Required parameters:
- group column name.
- case label and control label.
- organism and gene ID type for downstream annotation.

Stop and ask if:
- matrix columns do not match metadata samples.
- values are not non-negative counts but the user asks for DESeq2.
- contrast, group labels, organism, or ID type are missing.

Expression scale handling for `scripts/quickplot_bulk.py`:
- Default `--expression-scale auto`.
- Integer non-negative matrices are treated as `raw-count` and transformed with `log2(count + 1)` for exploratory statistics and plots.
- Small non-integer matrices are treated as `log-expression` and used as-is; do not apply a second `log2(x + 1)`.
- Large non-integer abundance/intensity matrices can be treated as `linear-expression` and transformed with `log2(value + 1)`.
- When the upstream scale is known, pass `--expression-scale raw-count`, `--expression-scale log-expression`, or `--expression-scale linear-expression` explicitly and record it in `methods.md`.

Default exploratory thresholds:
- `abs(log2FC) >= 1`
- `padj < 0.05`

Formal-route equivalent for already-normalized/log-expression input: `scripts/run_limma_trend.R` (limma `eBayes(trend=TRUE)`, no voom). Use this instead of hand-writing a custom script when R/Bioconductor is available and the input is not raw counts -- this includes matrices exported from Excel workbooks or public/collaborator pipelines, as long as they're first extracted into the `gene` + sample-column shape this contract describes.
