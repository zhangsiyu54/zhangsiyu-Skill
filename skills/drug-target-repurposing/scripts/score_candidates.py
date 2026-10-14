#!/usr/bin/env python3
"""Validate drug/target evidence and produce deterministic, auditable tiers."""

import argparse
import csv
import json
import os
import sys
from collections import defaultdict

TARGET_REQUIRED = {"target", "disease", "population", "tissue", "cell_type", "spatial_region", "desired_direction",
                   "evidence_class", "source", "estimate", "uncertainty", "direction", "independence_group", "qc_status"}
DRUG_REQUIRED = {"compound", "target", "mechanism", "modulation_direction", "approval_status", "indication",
                 "potency", "selectivity", "dose", "exposure_margin", "route", "tissue_penetration",
                 "target_organ_risk", "population_constraint", "renal_hepatic_constraint", "interaction_flag",
                 "contraindication", "trial_status", "failure_reason", "safety_flag", "source"}
DEFAULT_WEIGHTS = {"human_genetics": 3.0, "molecular_causality": 2.0, "disease_context": 1.5,
                   "perturbation_rescue": 2.0, "pharmacology": 2.0, "clinical_evidence": 1.5,
                   "safety": 2.0, "tractability": 1.0}


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


def truthy(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y", "pass", "acceptable", "plausible"}


def score_value(row):
    raw = row.get("score", "")
    if raw == "":
        return 1.0
    value = float(raw)
    if not 0 <= value <= 1:
        raise ValueError(f"score must be between 0 and 1, got {raw!r}")
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-evidence", required=True)
    parser.add_argument("--drug-evidence", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--weights", help="JSON file overriding predeclared class weights")
    args = parser.parse_args()
    weights = dict(DEFAULT_WEIGHTS)
    if args.weights:
        with open(args.weights, encoding="utf-8") as handle:
            custom = json.load(handle)
        unknown = sorted(set(custom) - set(weights))
        if unknown:
            raise ValueError(f"unknown weight classes: {', '.join(unknown)}")
        for key, value in custom.items():
            if float(value) < 0:
                raise ValueError(f"weight {key} must be non-negative")
            weights[key] = float(value)

    tf, targets = read_csv(args.target_evidence)
    df, drugs = read_csv(args.drug_evidence)
    require(tf, TARGET_REQUIRED, "target evidence")
    require(df, DRUG_REQUIRED, "drug evidence")
    errors, target_map = [], defaultdict(list)
    for line, row in enumerate(targets, 2):
        if not row["target"].strip() or not row["source"].strip() or not row["independence_group"].strip():
            errors.append(f"target row {line}: target, source, and independence_group are required")
        if row["evidence_class"] not in weights:
            errors.append(f"target row {line}: unknown evidence_class {row['evidence_class']!r}")
        try:
            row["_score"] = score_value(row)
        except ValueError as exc:
            errors.append(f"target row {line}: {exc}")
        target_map[row["target"].strip()].append(row)
    for line, row in enumerate(drugs, 2):
        if not row["compound"].strip() or not row["target"].strip() or not row["source"].strip():
            errors.append(f"drug row {line}: compound, target, and source are required")
        for field in ("target_organ_risk", "population_constraint", "renal_hepatic_constraint",
                      "interaction_flag", "contraindication", "safety_flag"):
            if not row[field].strip():
                errors.append(f"drug row {line}: safety field {field} is required; use an explicit none/unknown value")
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 2

    output = []
    for drug in sorted(drugs, key=lambda r: (r["compound"], r["target"])):
        target = drug["target"].strip()
        rows = [r for r in target_map.get(target, []) if r["qc_status"].strip().lower() in {"pass", "ok", "included"}]
        desired = sorted({r["desired_direction"].strip().lower() for r in rows if r["desired_direction"].strip()})
        modulation = drug["modulation_direction"].strip().lower()
        contradictions = []
        if len(desired) > 1:
            contradictions.append("conflicting desired directions")
        if desired and modulation not in desired and modulation != "unknown":
            contradictions.append(f"drug direction {modulation} conflicts with desired {desired[0]}")
        evidence_directions = {r["direction"].strip().lower() for r in rows if r["direction"].strip().lower() not in {"", "unknown", "neutral"}}
        if len(evidence_directions) > 1:
            contradictions.append("target evidence has contradictory directions")

        gates = []
        if not rows:
            gates.append("unresolved target identity")
        if not drug["mechanism"].strip() or modulation in {"", "unknown"}:
            gates.append("mechanism direction unknown")
        if not truthy(drug["tissue_penetration"]) or not truthy(drug["exposure_margin"]):
            gates.append("relevant exposure not plausible")
        if drug["safety_flag"].strip().lower() in {"fail", "unsafe", "disqualifying", "red"} or truthy(drug["contraindication"]):
            gates.append("disqualifying safety or contraindication")
        if contradictions and any("conflicts with desired" in x for x in contradictions):
            gates.append("direction conflict")

        by_class = defaultdict(dict)
        for row in rows:
            group = row["independence_group"].strip()
            by_class[row["evidence_class"]][group] = max(by_class[row["evidence_class"]].get(group, 0), row["_score"])
        class_scores = {key: (sum(groups.values()) / len(groups) if groups else None) for key, groups in by_class.items()}
        missing = [key for key in weights if key not in class_scores]
        numerator = sum(weights[key] * value for key, value in class_scores.items())
        denominator = sum(weights[key] for key in class_scores) or 1.0
        score = numerator / denominator
        has_genetics = class_scores.get("human_genetics", 0) > 0
        context = class_scores.get("molecular_causality", 0) > 0 or class_scores.get("disease_context", 0) > 0
        acceptable_safety = not any("safety" in gate for gate in gates)
        if gates:
            tier = "Excluded"
        elif has_genetics and context and acceptable_safety and not contradictions:
            tier = "Tier 1"
        elif score >= 0.5 and len(missing) <= 4:
            tier = "Tier 2"
        else:
            tier = "Tier 3"
        output.append({
            "candidate": drug["compound"], "target": target, "tier": tier,
            "normalized_score": f"{score:.6f}", "hard_gate_failures": "; ".join(gates),
            "contradictions": "; ".join(contradictions), "missing_evidence": "; ".join(missing),
            "class_scores_json": json.dumps(class_scores, sort_keys=True),
            "weights_json": json.dumps(weights, sort_keys=True),
            "evidence_rows_used": len(rows), "independent_groups": len({r['independence_group'] for r in rows}),
        })
    fields = ["candidate", "target", "tier", "normalized_score", "hard_gate_failures", "contradictions",
              "missing_evidence", "class_scores_json", "weights_json", "evidence_rows_used", "independent_groups"]
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(output)
    print(f"validated {len(targets)} target-evidence and {len(drugs)} drug-evidence rows; wrote {len(output)} candidates")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
