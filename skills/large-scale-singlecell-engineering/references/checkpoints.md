# Required checkpoints

1. `00_manifest`: source, checksum, study, donor, sample, tissue, condition, platform.
2. `01_ingested`: one validated sparse object per sample.
3. `02_qc`: retained cells, thresholds, doublets, ambient-RNA decisions, QC plots.
4. `03_harmonized`: stable gene IDs and controlled metadata vocabulary.
5. `04_sketch`: representative stratified subset and integration benchmark.
6. `05_integrated`: full latent space, neighbors, embeddings, model files.
7. `06_annotated`: labels, confidence, evidence, unresolved cells.
8. `07_pseudobulk`: donor-by-cell-type counts and design matrices.
9. `08_mechanism`: compartment-specific communication, trajectory, and regulatory results.
10. `09_release`: objects, tables, figures, environment, logs, and checksums.

Each checkpoint must be restartable without recomputing earlier valid stages.
