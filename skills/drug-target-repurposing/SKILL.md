---
name: drug-target-repurposing
description: Prioritize therapeutic targets and repurposing candidates by integrating human genetics, molecular QTL, bulk and single-cell expression, spatial localization, perturbation signatures, target-drug relationships, clinical evidence, pharmacology, safety, and tractability. Use for old-drug repurposing, target validation, signature reversal, drug ranking, mechanism-aligned intervention direction, combination screening, or experimental candidate selection.
---

# Drug target repurposing

Rank hypotheses, not prescriptions. Require the proposed intervention direction to agree with genetic and molecular evidence whenever direction can be inferred.

## Define the therapeutic hypothesis

1. Specify disease, stage, target population, tissue, cell type, spatial compartment, target, desired modulation direction, endpoint, and acceptable safety constraints.
2. Separate target-first, drug-first, and disease-signature-reversal analyses.
3. Freeze primary evidence sources and scoring rules before viewing final ranks.

## Assemble independent evidence

- Human genetics: locus, fine-mapping, MR, colocalization, loss/gain-of-function, phenome-wide effects.
- Molecular mechanism: eQTL/pQTL/mQTL, pathway, regulatory network, perturbation, rescue experiments.
- Disease context: bulk expression, single-cell subtype, spatial localization, stage specificity.
- Pharmacology: binding, potency, selectivity, agonist/antagonist/degrader mechanism, exposure, tissue penetration.
- Translation: approved indication, clinical trials, failures, contraindications, interactions, toxicity, and patent/formulation constraints.

Do not double-count mirrored databases or publications reporting the same experiment.

## Rank candidates

1. Apply hard exclusions for direction conflict, unacceptable safety, absent relevant exposure, or invalid target identity.
2. Score evidence classes separately and display missing evidence; do not convert missing data into positive evidence.
3. Penalize correlated evidence, contradictory direction, broad promiscuity, and non-independent validation.
4. Report uncertainty and produce tiers rather than a falsely precise single rank.
5. Propose the cheapest decisive validation experiment for each top candidate.

For expression-signature reversal, require disease and drug signatures to be compatible in cell type, dose, duration, platform, and direction. For combinations, report Bliss/Loewe/HSA/ZIP only when the dose matrix supports them.

Read [references/evidence-model.md](references/evidence-model.md) for scoring and [references/contracts.md](references/contracts.md) for required tables.

## Executable candidate scoring

Use the dependency-light [scripts/score_candidates.py](scripts/score_candidates.py) to strictly validate every target/drug evidence field in [references/contracts.md](references/contracts.md) and create a deterministic audit table. Target evidence may include an optional `score` from 0 to 1; absent scores mean presence of qualifying evidence, not efficacy. Evidence duplicated within the same `independence_group` contributes only its maximum value. Target identity passes only when at least one QC-passing evidence row remains; blank or unknown modulation direction fails the mechanism-direction gate.

```bash
python3 scripts/score_candidates.py --target-evidence target_evidence.csv --drug-evidence drug_evidence.csv --output candidate_scores.csv
```

Default evidence-class weights are declared in the script. Freeze overrides before ranking in a JSON file and pass `--weights weights.json`. Output includes the exact weights and class scores, hard-gate failures, direction contradictions, missing classes, independent-group counts, and ordinal tiers. `Excluded` means a hard gate failed; Tier 1 requires aligned human genetics, molecular or disease-context support, acceptable safety, and no unresolved contradiction. Computational tiers are prioritization aids, not clinical recommendations.

## Hand off

- Obtain genetic evidence from `gwas-causal-target-analysis`.
- Obtain disease expression from `public-omics-meta-analysis`, `bulkrna-analysis`, `singlecell-analysis`, and `spatial-transcriptomics-analysis`.
- Use `crispr-screen-analysis` for perturbational target support and `network-pathway-analysis` for pathway context.

## Completion gate

Require target identity, intervention direction, evidence provenance, independence groups, safety review, clinical-failure review, uncertainty, alternative targets, validation plan, and explicit statement that computational ranking is not clinical efficacy evidence.
