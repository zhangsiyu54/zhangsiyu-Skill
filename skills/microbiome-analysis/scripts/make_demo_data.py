#!/usr/bin/env python3
import argparse
import csv
import random
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Create demo microbiome abundance, taxonomy and metadata tables.")
    parser.add_argument("--out", default="demo/microbiome")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    random.seed(37)
    samples = [f"C{i}" for i in range(1, 7)] + [f"T{i}" for i in range(1, 7)]
    taxa = ["Bacteroides", "Faecalibacterium", "Akkermansia", "Lactobacillus", "Prevotella", "Bifidobacterium", "Escherichia", "Roseburia"]
    with (out / "metadata.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", "group"])
        for s in samples:
            writer.writerow([s, "Treat" if s.startswith("T") else "Control"])
    with (out / "taxonomy.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["taxon", "genus"])
        for t in taxa:
            writer.writerow([t, t])
    with (out / "abundance.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["taxon", *samples])
        for taxon in taxa:
            row = [taxon]
            for sample in samples:
                shift = 55 if sample.startswith("T") and taxon in {"Akkermansia", "Lactobacillus"} else 0
                row.append(max(0, int(random.gauss(120 + shift, 35))))
            writer.writerow(row)
    print(f"demo written to {out}")


if __name__ == "__main__":
    main()
