# Immune Repertoire Input Contract

Accepted inputs:
- clonotype table with `sample`, `group`, `clonotype`, `count`, optional `celltype`, `v_gene`, `j_gene`, `cdr3`.
- VDJ contig annotation table from Cell Ranger or exported scTCR/scBCR workflow.

Required parameters:
- chain type: TCR alpha/beta, BCR heavy/light, or paired.
- group/sample column and clonotype count/frequency column.

Stop and ask if clonotype ID or sample grouping is missing.
