#!/usr/bin/env python3
import argparse
import csv
import random
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Create demo expression-clinical survival table.")
    parser.add_argument("--out", default="demo/public-db-survival")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    random.seed(17)
    with (out / "expression_clinical.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", "type", "GeneA", "survival_month", "status"])
        for i in range(1, 61):
            tumor = i <= 42
            expr = random.gauss(7.8 if tumor else 4.8, 0.9)
            high_risk = expr > 7.4
            time = max(3, random.gauss(28 if high_risk else 46, 10))
            status = 1 if random.random() < (0.68 if high_risk else 0.36) else 0
            writer.writerow([f"S{i:03d}", "Tumor" if tumor else "Normal", round(expr, 3), round(time, 1), status])
    print(f"demo written to {out}")


if __name__ == "__main__":
    main()
