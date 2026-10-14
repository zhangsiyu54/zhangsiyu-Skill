# Data contracts

## Target evidence

`target, disease, population, tissue, cell_type, spatial_region, desired_direction, evidence_class, source, estimate, uncertainty, direction, independence_group, qc_status`

## Drug evidence

`compound, target, mechanism, modulation_direction, approval_status, indication, potency, selectivity, dose, exposure_margin, route, tissue_penetration, target_organ_risk, population_constraint, renal_hepatic_constraint, interaction_flag, contraindication, trial_status, failure_reason, safety_flag, source`

## Candidate output

`candidate, target, tier, direction_match, genetics_support, mechanism_support, disease_context, pharmacology, clinical_evidence, safety, tractability, uncertainty, decisive_next_experiment`

## Deliverables

- target and compound identity map
- evidence matrix with provenance and independence groups
- exclusion and contradiction ledger
- tiered candidate list
- safety and clinical-failure summary
- validation experiment plan
