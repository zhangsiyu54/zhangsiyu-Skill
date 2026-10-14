---
name: large-scale-singlecell-engineering
description: Engineer and operate disk-backed, checkpointed, parallel single-cell workflows for roughly 500,000 cells or more across many samples, studies, tissues, and platforms. Use for large Seurat, AnnData, h5ad, 10x, BPCells, HDF5, Zarr, scVI/scANVI, HPC, Slurm, memory planning, chunking, atlas integration, donor-level pseudobulk, or recovery from failed large single-cell runs.
---

# Large-scale single-cell engineering

Own the computational architecture; hand biological interpretation to `singlecell-analysis`. Treat donors, not cells, as biological replicates.

## Plan before loading data

1. Inventory studies, donors, samples, tissues, conditions, platforms, chemistry, genome build, matrix type, and raw-data availability.
2. Estimate cells, genes, nonzero entries, disk footprint, peak memory, CPUs, GPUs, wall time, and scratch requirements.
3. Select a disk-backed representation: BPCells/Seurat v5, HDF5-backed AnnData, Zarr, or another validated sparse backend.
4. Define checkpoints and immutable stage outputs before execution.
5. Preserve raw counts; never overwrite them with normalized or integrated values.

## Execute in stages

1. Ingest and validate each sample independently.
2. Run per-sample QC, ambient-RNA review, doublet detection, and metadata validation. Store QC decisions rather than silently filtering.
3. Write a stable per-sample disk-backed object and checksum.
4. Build donor/sample summaries before merging.
5. Integrate with a method appropriate to scale and biology. Use scVI/scANVI, reciprocal mapping, Harmony, or hierarchical integration only after testing that biological signals are retained.
6. Use stratified sketches for method selection and visualization, then project or fit the selected method at full scale.
7. Perform condition testing at donor/sample level with pseudobulk or justified mixed models. Do not use cells as independent replicates.
8. Run expensive communication, trajectory, and regulatory analyses within justified compartments or representative subsets, then validate conclusions at donor level.

Read [references/architecture.md](references/architecture.md) for the 1 TB shared-server design, [references/checkpoints.md](references/checkpoints.md) for required stage gates, and [references/operations.md](references/operations.md) before generating scheduler jobs or writing shared disk-backed objects.

## Executable planning tools

These dependency-light Python tools use only the standard library and do not load expression matrices or submit jobs:

1. Validate a CSV/TSV manifest containing `study_id`, `donor_id`, `sample_id`, and `input_path`; optionally include `cells`, `nonzero_entries`, and `file_size_bytes`:
   `python3 scripts/validate_manifest.py manifest.tsv --json-out manifest_validation.json`
2. Produce a transparent heuristic resource plan:
   `python3 scripts/estimate_resources.py manifest.tsv --out resource_plan.json`
3. Generate a non-overwriting Slurm array/checkpoint layout without submitting it:
   `python3 scripts/generate_slurm_layout.py manifest.tsv resource_plan.json --out-dir planned_jobs`

The generator binds the resource plan to the exact manifest, creates collision-resistant per-row output identifiers, and emits the full 01-09 checkpoint DAG. Only stage `01_ingested` receives an executable template; stages 02-09 remain explicit placeholders until reviewed stage-specific drivers exist. The ingestion wrapper uses unique temporary directories, quarantines failed/old partials, validates non-empty output, atomically renames successful output, and records a structured completion marker.

Treat the resource plan as a pre-pilot estimate only. Before production, benchmark small, median, and large samples and replace the generated envelope with measured upper-quantile requirements plus headroom. Review scheduler account, partition, QOS, concurrency, paths, and the analysis driver before manually calling `sbatch`.

## Safety rules

- Keep matrices sparse; audit every conversion that could densify data.
- Do not concatenate studies before validating gene identifiers and metadata.
- Do not remove study effects until confirming that study is not confounded with the biology of interest.
- Do not claim increased power solely from more cells when donor count is unchanged.
- Use scheduler resource limits on shared servers; never request all RAM by default.

## Deliverables

Produce an input manifest, resource plan, stage DAG, scheduler scripts, environment lock, QC ledger, checkpoint manifest, integrated object, donor-level analysis tables, figures, logs, checksums, and a rerun guide.
