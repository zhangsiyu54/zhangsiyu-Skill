#!/usr/bin/env python3
import argparse
import csv
import random
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Create demo bulk RNA-seq count matrix and metadata.")
    parser.add_argument("--out", default="demo/bulkrna", help="Output directory")
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    random.seed(7)

    samples = [f"C{i}" for i in range(1, 4)] + [f"T{i}" for i in range(1, 4)]
    with (out / "metadata.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", "group"])
        for sample in samples:
            writer.writerow([sample, "Control" if sample.startswith("C") else "Treat"])

    with (out / "count_matrix.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["gene", *samples])
        for i in range(1, 121):
            base = random.randint(40, 900)
            effect = 1.0
            if i <= 12:
                effect = 2.6
            elif 13 <= i <= 24:
                effect = 0.38
            row = [f"Gene{i:03d}"]
            for sample in samples:
                value = base * (effect if sample.startswith("T") else 1.0)
                row.append(max(1, int(random.gauss(value, max(3, value * 0.12)))))
            writer.writerow(row)

    print(f"demo written to {out}")


if __name__ == "__main__":
    main()
