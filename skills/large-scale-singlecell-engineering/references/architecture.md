# Architecture for a shared 1 TB server

## Storage tiers

- Immutable raw: downloaded matrices, FASTQ, fragments, images, and source metadata.
- Working: disk-backed objects, indexes, models, and temporary shards on fast scratch.
- Results: small tables, figures, manifests, and versioned final objects.

## Scheduler policy

Use Slurm/PBS job arrays for per-sample stages. Request measured memory plus headroom. Cap a normal large job below total server memory and preserve operating-system headroom. Separate CPU-heavy, RAM-heavy, GPU, and I/O-heavy queues when possible.

## Representation choices

- Use BPCells with Seurat v5 for R-centric workflows that need disk-backed sparse matrices.
- Use backed AnnData/HDF5 or Zarr for Python/Scanpy workflows.
- Use scVI/scANVI with GPU for latent integration when the design supports it.
- Export compact interoperable objects at milestones; do not repeatedly convert the full atlas.

## Statistical unit

Maintain `study -> donor -> sample -> library -> cell` identifiers. Use donor/sample units for group comparisons. Include study/batch covariates or random effects when supported by the design.
