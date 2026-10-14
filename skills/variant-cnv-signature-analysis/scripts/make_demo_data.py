#!/usr/bin/env python3
import argparse
import csv
import random
from pathlib import Path

# Variant classifications used for the demo MAF, biased toward the functional
# classes maftools treats as non-synonymous by default (see read.maf vc_nonSyn),
# with a minority of silent variants included for realism (pitfall #6).
VC_CHOICES = (
    ["Missense_Mutation"] * 6
    + ["Nonsense_Mutation"]
    + ["Frame_Shift_Del"]
    + ["Frame_Shift_Ins"]
    + ["In_Frame_Del"]
    + ["Splice_Site"]
    + ["Silent"]
)
CHROM_MAP = {
    "TP53": "17", "KRAS": "12", "PIK3CA": "3", "EGFR": "7", "PTEN": "10",
    "APC": "5", "BRAF": "7", "IDH1": "2", "SMAD4": "18", "ARID1A": "1",
}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="demo/cancer-genomics")
    args = p.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    random.seed(77)
    genes = ["TP53", "KRAS", "PIK3CA", "EGFR", "PTEN", "APC", "BRAF", "IDH1", "SMAD4", "ARID1A"]
    samples = [f"S{i:02d}" for i in range(1, 13)]

    counts = {}
    with (out / "mutation_frequency.csv").open("w", newline="", encoding="utf-8") as h:
        w = csv.writer(h)
        w.writerow(["gene", "count", "frequency"])
        for i, g in enumerate(genes):
            count = max(1, 10 - i + random.randint(-1, 2))
            count = min(count, len(samples))
            counts[g] = count
            w.writerow([g, count, round(count / len(samples), 3)])

    with (out / "cnv_matrix.csv").open("w", newline="", encoding="utf-8") as h:
        w = csv.writer(h)
        w.writerow(["gene", *samples])
        for g in genes:
            w.writerow([g, *[random.choice([-2, -1, 0, 0, 0, 1, 2]) for _ in samples]])

    with (out / "signature_exposure.csv").open("w", newline="", encoding="utf-8") as h:
        w = csv.writer(h)
        w.writerow(["signature", *samples])
        for sig in ["SBS1", "SBS5", "SBS7", "SBS13"]:
            w.writerow([sig, *[round(random.random(), 3) for _ in samples]])

    # --- MAF-like mutation table (demo.maf) + clinical metadata for the optional
    # R/maftools route (scripts/run_maftools_signature.R). Mutation counts per
    # gene are drawn from `counts` above so the demo MAF and mutation_frequency.csv
    # stay internally consistent.
    maf_rows = []
    for g in genes:
        mutated_samples = random.sample(samples, counts[g])
        for s in mutated_samples:
            ref = random.choice("ACGT")
            alt = random.choice([b for b in "ACGT" if b != ref])
            maf_rows.append(
                {
                    "Hugo_Symbol": g,
                    "Chromosome": CHROM_MAP.get(g, "1"),
                    "Start_Position": random.randint(1_000_000, 100_000_000),
                    "End_Position": None,
                    "Reference_Allele": ref,
                    "Tumor_Seq_Allele2": alt,
                    "Variant_Classification": random.choice(VC_CHOICES),
                    "Variant_Type": "SNP",
                    "Tumor_Sample_Barcode": s,
                    "NCBI_Build": "GRCh37",
                }
            )
    for row in maf_rows:
        row["End_Position"] = row["Start_Position"]
    random.shuffle(maf_rows)
    maf_fields = [
        "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
        "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification",
        "Variant_Type", "Tumor_Sample_Barcode", "NCBI_Build",
    ]
    with (out / "demo.maf").open("w", newline="", encoding="utf-8") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(maf_fields)
        for row in maf_rows:
            w.writerow([row[f] for f in maf_fields])

    with (out / "demo_clinical.tsv").open("w", newline="", encoding="utf-8") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["Tumor_Sample_Barcode", "cohort_group"])
        for i, s in enumerate(samples):
            w.writerow([s, "GroupA" if i % 2 == 0 else "GroupB"])

    # --- Per-sample TMB table (tmb_summary.csv) for the Python quickplot TMB
    # boxplot. Mixes a "targeted_panel" (small, hotspot-enriched footprint) and
    # a "WES" assay on purpose so the boxplot demonstrates why TMB must be
    # reported/compared per-assay rather than pooled (see pitfalls.md #1).
    assay_size_mb = {"targeted_panel": 1.5, "WES": 34.0}
    with (out / "tmb_summary.csv").open("w", newline="", encoding="utf-8") as h:
        w = csv.writer(h)
        w.writerow(["sample", "assay", "panel_size_mb", "mutation_count"])
        for i, s in enumerate(samples):
            assay = "WES" if i % 4 == 0 else "targeted_panel"
            panel_gene_hits = sum(1 for row in maf_rows if row["Tumor_Sample_Barcode"] == s)
            background = random.randint(30, 160) if assay == "WES" else random.randint(0, 4)
            mutation_count = panel_gene_hits + background
            w.writerow([s, assay, assay_size_mb[assay], mutation_count])

    print(f"demo written to {out}")


if __name__ == "__main__":
    main()
