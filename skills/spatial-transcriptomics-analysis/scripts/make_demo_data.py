#!/usr/bin/env python3
import argparse
import csv
import math
import random
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Create demo spatial spot and expression tables.")
    parser.add_argument("--out", default="demo/spatial")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    random.seed(47)
    with (out / "spots.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["spot", "x", "y", "domain", "nCount"])
        for i in range(18):
            for j in range(14):
                if (i - 8.5) ** 2 + (j - 6.5) ** 2 > 70:
                    continue
                domain = "Tumor" if i > 9 and j > 5 else "Stroma" if i < 7 else "Immune"
                ncount = int(random.gauss(6200, 900))
                writer.writerow([f"spot_{i}_{j}", i, j, domain, ncount])
    with (out / "expression.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["spot", "GeneA", "GeneB"])
        for i in range(18):
            for j in range(14):
                if (i - 8.5) ** 2 + (j - 6.5) ** 2 > 70:
                    continue
                gene_a = max(0, random.gauss(1.0 + (i / 18) * 2.4 + math.sin(j / 2), 0.35))
                gene_b = max(0, random.gauss(2.0 + (j / 14) * 1.8, 0.35))
                writer.writerow([f"spot_{i}_{j}", round(gene_a, 3), round(gene_b, 3)])
    print(f"demo written to {out}")


if __name__ == "__main__":
    main()
