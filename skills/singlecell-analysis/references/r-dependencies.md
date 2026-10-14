# R Package Requirements

Different scripts need different package sets — check per-script before running:

| Script | Packages |
| --- | --- |
| `run_seurat.R` | `Seurat`, `ggplot2` |
| `run_pseudobulk_de.R` | `Seurat`, `DESeq2`, `Matrix`, `ggplot2` |
| `run_cellchat.R` | `Seurat`, `CellChat` (GitHub: `jinworks/CellChat`) |
| `run_trajectory_monocle3.R` | `Seurat`, `monocle3` (GitHub: `cole-trapnell-lab/monocle3`), `SeuratWrappers` (GitHub: `satijalab/seurat-wrappers`) |

```r
BiocManager::install(c("Seurat","DESeq2","Matrix","patchwork","ggplot2"))
remotes::install_github(c("cole-trapnell-lab/monocle3","satijalab/seurat-wrappers","jinworks/CellChat"))
```

```bash
python3 scripts/check_environment.py --out output --r-packages Seurat,DESeq2,Matrix,ggplot2,monocle3,SeuratWrappers,CellChat
```

Note: `monocle3`'s dependency `BPCells` needs a working system C++17 toolchain and `libhdf5`
(`brew install hdf5` on macOS) to compile from source — if `remotes::install_github` fails with
"Unable to compile program with C++17 std::filesystem" even after installing hdf5, the actual
root cause is usually an incomplete Xcode Command Line Tools install (missing libc++ headers),
not a missing R/system package: `sudo rm -rf /Library/Developer/CommandLineTools && xcode-select --install`.

Verified (2026-07): `run_seurat.R` and `run_pseudobulk_de.R` ran end to end with real Seurat +
DESeq2 (clustering correctly recovered a synthetic 3-cell-type structure). `run_cellchat.R` and
`run_trajectory_monocle3.R` were fixed for a `GetAssayData(slot=...)` incompatibility with
SeuratObject >= 5.0 (now uses `layer=` with a fallback) but could not be executed end to end in
the authoring environment due to the monocle3/BPCells toolchain issue above — validate once in
an environment with a working C++ toolchain.
