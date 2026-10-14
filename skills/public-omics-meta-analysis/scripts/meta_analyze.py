#!/usr/bin/env python3
"""Validate cohort/effect tables and run transparent inverse-variance meta-analysis."""

import argparse
import csv
import math
import os
import sys
from collections import defaultdict

REQUIRED_EFFECT = {"cohort_id", "feature_id", "contrast", "effect", "se", "p", "n_case", "n_control",
                   "direction", "model", "covariates", "qc_status"}
REQUIRED_COHORT = {"cohort_id", "accession", "repository", "publication", "country", "ancestry", "tissue",
                   "disease_definition", "comparator", "assay", "platform", "n_donors", "n_samples", "paired",
                   "role", "overlap_group", "inclusion_status", "exclusion_reason"}


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        return reader.fieldnames, list(reader)


def require(fields, required, label):
    missing = sorted(required - set(fields))
    if missing:
        raise ValueError(f"{label}: missing columns: {', '.join(missing)}")


def finite(value, label):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{label}: expected numeric value, got {value!r}")
    if not math.isfinite(number):
        raise ValueError(f"{label}: value must be finite")
    return number


def normal_p(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def synthesize(studies):
    k = len(studies)
    weights = [1.0 / (se * se) for _, _, se in studies]
    sw = sum(weights)
    fixed = sum(w * effect for w, (_, effect, _) in zip(weights, studies)) / sw
    q = sum(w * (effect - fixed) ** 2 for w, (_, effect, _) in zip(weights, studies))
    df = k - 1
    c = sw - sum(w * w for w in weights) / sw
    tau2 = max(0.0, (q - df) / c) if k > 1 and c > 0 else 0.0
    random_w = [1.0 / (se * se + tau2) for _, _, se in studies]
    srw = sum(random_w)
    random = sum(w * effect for w, (_, effect, _) in zip(random_w, studies)) / srw
    fixed_se = math.sqrt(1.0 / sw)
    random_se = math.sqrt(1.0 / srw)
    i2 = max(0.0, (q - df) / q) * 100.0 if q > 0 and k > 1 else 0.0
    pred_half = 1.96 * math.sqrt(tau2 + random_se * random_se) if k >= 3 else None
    return {
        "k": k, "fixed_effect": fixed, "fixed_se": fixed_se,
        "fixed_ci_low": fixed - 1.96 * fixed_se, "fixed_ci_high": fixed + 1.96 * fixed_se,
        "fixed_p": normal_p(fixed / fixed_se), "random_effect": random, "random_se": random_se,
        "random_ci_low": random - 1.96 * random_se, "random_ci_high": random + 1.96 * random_se,
        "random_p": normal_p(random / random_se), "q": q, "q_df": df,
        "i2_percent": i2, "tau2_dl": tau2,
        "prediction_low": random - pred_half if pred_half is not None else "",
        "prediction_high": random + pred_half if pred_half is not None else "",
    }


def write_rows(path, rows, fields):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--effects", required=True)
    parser.add_argument("--cohorts")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--min-cohorts", type=int, default=2)
    args = parser.parse_args()

    errors = []
    cohort_ids = None
    included_cohorts = None
    if args.cohorts:
        fields, cohorts = read_csv(args.cohorts)
        require(fields, REQUIRED_COHORT, "cohort table")
        cohort_ids, included_cohorts = set(), set()
        for line, row in enumerate(cohorts, 2):
            cid = row["cohort_id"].strip()
            if not cid or cid in cohort_ids:
                errors.append(f"cohort row {line}: blank or duplicate cohort_id {cid!r}")
            cohort_ids.add(cid)
            try:
                n_donors = int(row["n_donors"])
                if n_donors <= 0:
                    raise ValueError
            except ValueError:
                errors.append(f"cohort row {line}: n_donors must be a positive integer")
            if row["inclusion_status"].strip().lower() not in {"included", "excluded", "pending"}:
                errors.append(f"cohort row {line}: inclusion_status must be included, excluded, or pending")
            elif row["inclusion_status"].strip().lower() == "included":
                included_cohorts.add(cid)
            try:
                n_samples = int(row["n_samples"])
                if n_samples <= 0:
                    raise ValueError
            except ValueError:
                errors.append(f"cohort row {line}: n_samples must be a positive integer")

    fields, effects = read_csv(args.effects)
    require(fields, REQUIRED_EFFECT, "effect table")
    groups = defaultdict(list)
    seen = set()
    for line, row in enumerate(effects, 2):
        key = (row["cohort_id"].strip(), row["feature_id"].strip(), row["contrast"].strip())
        if not all(key):
            errors.append(f"effect row {line}: cohort_id, feature_id, and contrast are required")
            continue
        if key in seen:
            errors.append(f"effect row {line}: duplicate cohort/feature/contrast {key}")
            continue
        seen.add(key)
        if cohort_ids is not None and key[0] not in cohort_ids:
            errors.append(f"effect row {line}: unknown cohort_id {key[0]!r}")
        elif included_cohorts is not None and key[0] not in included_cohorts:
            errors.append(f"effect row {line}: cohort_id {key[0]!r} is not included")
        try:
            effect = finite(row["effect"], f"effect row {line} effect")
            se = finite(row["se"], f"effect row {line} se")
            if se <= 0:
                raise ValueError(f"effect row {line}: se must be > 0")
            p_value = finite(row["p"], f"effect row {line} p")
            if not 0 <= p_value <= 1:
                raise ValueError(f"effect row {line}: p must be between 0 and 1")
            for field in ("n_case", "n_control"):
                if int(row[field]) <= 0:
                    raise ValueError(f"effect row {line}: {field} must be a positive integer")
            if row["qc_status"].strip().lower() not in {"pass", "included", "ok"}:
                raise ValueError(f"effect row {line}: qc_status must pass before synthesis")
            for field in ("direction", "model", "covariates"):
                if not row[field].strip():
                    raise ValueError(f"effect row {line}: {field} is required")
            groups[(key[1], key[2])].append((key[0], effect, se))
        except (ValueError, TypeError) as exc:
            errors.append(str(exc))
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 2

    summary, loso = [], []
    for (feature, contrast), studies in sorted(groups.items()):
        if len(studies) < args.min_cohorts:
            continue
        result = synthesize(studies)
        summary.append({"feature_id": feature, "contrast": contrast, **result})
        if len(studies) >= 3:
            for omitted, _, _ in studies:
                reduced = [x for x in studies if x[0] != omitted]
                rr = synthesize(reduced)
                loso.append({"feature_id": feature, "contrast": contrast, "omitted_cohort": omitted,
                             "k": rr["k"], "random_effect": rr["random_effect"],
                             "random_se": rr["random_se"], "random_p": rr["random_p"],
                             "i2_percent": rr["i2_percent"], "tau2_dl": rr["tau2_dl"]})
    if not summary:
        print("ERROR: no feature/contrast group met --min-cohorts", file=sys.stderr)
        return 2
    summary_fields = ["feature_id", "contrast", "k", "fixed_effect", "fixed_se", "fixed_ci_low",
                      "fixed_ci_high", "fixed_p", "random_effect", "random_se", "random_ci_low",
                      "random_ci_high", "random_p", "q", "q_df", "i2_percent", "tau2_dl",
                      "prediction_low", "prediction_high"]
    loso_fields = ["feature_id", "contrast", "omitted_cohort", "k", "random_effect", "random_se",
                   "random_p", "i2_percent", "tau2_dl"]
    write_rows(args.out_prefix + ".meta.csv", summary, summary_fields)
    write_rows(args.out_prefix + ".loso.csv", loso, loso_fields)
    print(f"validated {len(effects)} effects; wrote {len(summary)} meta-analysis rows and {len(loso)} LOSO rows")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
