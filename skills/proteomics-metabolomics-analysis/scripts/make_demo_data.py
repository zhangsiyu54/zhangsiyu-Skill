#!/usr/bin/env python3
import argparse
import csv
import random
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Create demo proteomics/metabolomics matrix and DEG-like table.")
    parser.add_argument("--out", default="demo/omics")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    random.seed(67)
    samples = [f"C{i}" for i in range(1, 6)] + [f"T{i}" for i in range(1, 6)]
    with (out / "metadata.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", "group"])
        for s in samples:
            writer.writerow([s, "Treat" if s.startswith("T") else "Control"])
    features = [f"Feature_{i:03d}" for i in range(1, 81)]
    with (out / "matrix.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["feature", *samples])
        for i, feature in enumerate(features, 1):
            row = [feature]
            for s in samples:
                shift = 1.1 if s.startswith("T") and i <= 12 else -0.9 if s.startswith("T") and 13 <= i <= 22 else 0
                row.append(round(random.gauss(24 + shift, 0.55), 3))
            writer.writerow(row)
    with (out / "differential_features.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["feature", "log2FC", "padj"])
        for i, feature in enumerate(features, 1):
            fc = random.gauss(0, 0.35) + (1.2 if i <= 12 else -1.0 if 13 <= i <= 22 else 0)
            padj = min(1, 10 ** (-abs(fc) * 2.0))
            writer.writerow([feature, round(fc, 3), round(padj, 6)])
    print(f"demo written to {out}")


if __name__ == "__main__":
    main()
