# Variant CNV Signature Input Contract

Accepted inputs:
- mutation frequency table with `gene`, `count`, `frequency`.
- CNV matrix with first column `gene`, remaining columns sample IDs.
- signature exposure table with `signature`, sample columns.
- optional per-sample TMB table with `sample`, `assay`, `panel_size_mb`, `mutation_count` (used by `scripts/quickplot_cancer_genomics.py --tmb`; TMB is computed as `mutation_count / panel_size_mb`). Always state the assay per sample and never pool TMB across assays without harmonization (see `references/pitfalls.md` #1).
- MAF-like mutation table (tab-delimited, standard MAF columns: `Hugo_Symbol`, `Chromosome`, `Start_Position`, `End_Position`, `Reference_Allele`, `Tumor_Seq_Allele2`, `Variant_Classification`, `Variant_Type`, `Tumor_Sample_Barcode`) for the formal R/maftools route (`scripts/run_maftools_signature.R`), with an optional clinical/metadata table (`Tumor_Sample_Barcode` plus a grouping column) for oncoplot annotation.

Required parameters:
- cohort/sample definition.
- genome/build and gene identifier style when raw variant calls are used.
- grouping column when comparing subgroups.
- for TMB: the assay/panel name and its callable footprint in Mb (never a generic "30 Mb exome" constant).

Stop and ask if sample IDs or mutation/CNV units are ambiguous.
