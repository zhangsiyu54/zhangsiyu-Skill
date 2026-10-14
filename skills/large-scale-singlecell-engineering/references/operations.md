# Production operations contract

## Resource profiling

Benchmark one small, one median, and one large sample. Record peak RSS, output bytes, nonzero entries, wall time, and CPU efficiency. Derive Slurm resource classes from measured upper quantiles plus headroom; cap concurrent jobs by RAM, scratch throughput, and file-descriptor limits rather than CPU count alone. Generate sbatch files only after this pilot.

## Atomic and restartable stages

Write each output to a unique temporary path, validate schema/counts/checksum, atomically rename, then create a completion marker containing input hashes, command, versions, seed, resources, and exit status. A retry must ignore a valid completed stage, quarantine corrupt partial output, and use a bounded retry count. Clean temporary files only after validation.

## Concurrent I/O

Assume HDF5/AnnData is unsafe for uncontrolled concurrent writes. Use one writer per object or independent shards followed by a validated merge. Never let job arrays write to the same file. Prefer immutable stage outputs.

## Sketch and integration acceptance

Stratify by study, donor, condition, tissue, and preliminary cell state; set minimum quotas to protect rare populations. Compare methods on biological conservation, batch mixing, label transfer confidence, donor representation, neighborhood stability, and rare-cell retention. Pre-register acceptance thresholds after pilot benchmarking.

## Reproducibility and downloads

Pin container digest or lockfile, R/Python packages, reference versions, random seeds, and GPU/CUDA stack. Record source license, release, URL, checksum, download date, retry/resume behavior, and cache path. Never treat a mutable URL as a version identifier.
