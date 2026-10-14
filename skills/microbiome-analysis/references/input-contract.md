# Microbiome Input Contract

Accepted inputs:
- ASV/OTU table: first column `taxon` or `feature`, remaining columns sample IDs.
- taxonomy table: `taxon`, `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species` when available.
- metadata table: one row per sample with `sample` and group columns.
- existing alpha diversity, beta distance, PCoA or differential abundance tables.

Required parameters:
- group/condition column.
- target taxonomy rank for composition plots.
- whether counts are raw counts, rarefied counts, or relative abundance.

Stop and ask if:
- sample IDs do not match.
- taxonomy rank is requested but taxonomy table is missing.
- differential abundance is requested from compositional relative abundance only.
