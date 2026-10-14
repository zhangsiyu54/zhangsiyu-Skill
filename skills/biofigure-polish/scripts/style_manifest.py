#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def guess_type(name):
    lower = name.lower()
    for key in ["volcano", "heatmap", "umap", "dotplot", "boxplot", "survival", "km", "bubble", "barplot", "pca"]:
        if key in lower:
            return key
    return "figure"


def main():
    parser = argparse.ArgumentParser(description="Create a figure polishing manifest.")
    parser.add_argument("--figures", default="figures")
    parser.add_argument("--out", default="output/figure_manifest.tsv")
    args = parser.parse_args()
    fig_dir = Path(args.figures)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if not fig_dir.exists():
        print(f"warning: figures directory not found: {fig_dir} -- writing an empty manifest")
    elif not fig_dir.is_dir():
        raise SystemExit(f"error: --figures path is not a directory: {fig_dir}")
    candidates = sorted([p for p in fig_dir.glob("**/*") if p.suffix.lower() in {".png", ".pdf"}]) if fig_dir.exists() else []
    chosen = {}
    for path in candidates:
        # Key by path relative to fig_dir (not just the filename stem) so that two
        # different figures in different subfolders that happen to share a basename
        # (e.g. "panel1.png" under two experiment folders) are inventoried as
        # separate figures instead of silently colliding into one manifest row.
        # Same-name PNG/PDF pairs *within the same folder* still dedupe, preferring
        # the PNG as the primary source.
        key = str(path.relative_to(fig_dir).with_suffix(""))
        if key not in chosen or path.suffix.lower() == ".png":
            chosen[key] = path
    files = sorted(chosen.values())
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["source_file", "output_stem", "panel_label", "figure_type", "changes", "status"])
        for idx, path in enumerate(files):
            label = chr(ord("A") + idx)
            writer.writerow([str(path), f"output/figures/{label}_{path.stem}", label, guess_type(path.name), "standardize panel label/export png+pdf", "planned"])
    if fig_dir.exists() and not files:
        print(f"warning: no .png/.pdf files found under {fig_dir}")
    print(f"manifest written to {out}; figures={len(files)}")


if __name__ == "__main__":
    main()
