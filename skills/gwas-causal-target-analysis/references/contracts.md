# Input and output contracts

## Minimum study manifest

`study_id, source, release, trait, trait_definition, role, ancestry, genome_build, n_total, n_cases, n_controls, unit, sample_overlap, access_date`

## Harmonized variant table

`study_id, chr, pos, rsid, effect_allele, other_allele, beta, se, p, eaf, n, info, qc_status, qc_reason`

Represent each biallelic SNP or indel as one row. Alleles must be explicit A/C/G/T sequences. Split multi-allelic records upstream; do not encode several alleles with comma, slash, or pipe separators. The generic schema validator also accepts QTL summary statistics, but platform, molecular trait, tissue/cell context, and cis/trans definitions require study-specific manual validation.

## Target evidence table

`target, locus, evidence_class, dataset, tissue, cell_type, direction, estimate, uncertainty, significance, independence_group, qc_status, source`

## Deliverables

- study and phenotype manifest
- harmonization/QC report
- locus and replication table
- fine-mapping credible sets
- MR sensitivity report
- colocalization prior-sensitivity report
- target evidence matrix and ranked list
- commands, versions, seeds, logs, and null/failed analyses
