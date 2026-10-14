#!/usr/bin/env python3
import argparse
import csv
import random
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Create demo ATAC differential peak and motif tables.")
    parser.add_argument("--out", default="demo/atac")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    random.seed(57)
    samples = [f"C{i}" for i in range(1, 5)] + [f"T{i}" for i in range(1, 5)]
    with (out / "peak_counts.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["peak", *samples])
        for i in range(1, 51):
            peak = f"chr1:{100000+i*900}-{100420+i*900}"
            row = [peak]
            for s in samples:
                shift = 50 if s.startswith("T") and i <= 12 else 0
                row.append(max(0, int(random.gauss(140 + shift, 30))))
            writer.writerow(row)
    with (out / "differential_peaks.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["peak", "log2FC", "padj"])
        for i in range(1, 51):
            fc = random.gauss(0, 0.7) + (1.4 if i <= 10 else -1.2 if 11 <= i <= 18 else 0)
            padj = min(1, 10 ** (-abs(fc) * 1.5))
            writer.writerow([f"chr1:{100000+i*900}-{100420+i*900}", round(fc, 3), round(padj, 6)])
    motifs = ["IRF1", "STAT1", "JUN", "FOS", "CTCF", "RUNX1", "NFkB"]
    with (out / "motif_enrichment.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["motif", "tf", "padj", "enrichment"])
        for motif in motifs:
            enr = random.uniform(1.3, 4.5)
            writer.writerow([motif, motif, round(10 ** (-enr), 6), round(enr, 3)])
    print(f"demo written to {out}")


if __name__ == "__main__":
    main()
