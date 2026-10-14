#!/usr/bin/env python3
"""Generate (but never submit) a generic Slurm array and checkpoint layout."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import stat
from pathlib import Path


def delimiter(path: Path) -> str:
    return "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","


def output_id(index: int, row: dict[str, str]) -> str:
    sample = (row.get("sample_id") or "").strip()
    slug = re.sub(r"[^A-Za-z0-9._-]", "_", sample)[:80] or "sample"
    identity = "\0".join((str(index), row.get("study_id", ""), row.get("donor_id", ""), sample, (row.get("input_path") or "").strip()))
    return f"{index:06d}_{slug}_{hashlib.sha256(identity.encode()).hexdigest()[:12]}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("resource_plan", type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--max-concurrent", type=int, default=4)
    args = parser.parse_args()
    if args.max_concurrent < 1:
        parser.error("--max-concurrent must be positive")
    with args.manifest.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter=delimiter(args.manifest)))
    if not rows:
        parser.error("manifest has no rows")
    plan = json.loads(args.resource_plan.read_text(encoding="utf-8"))
    manifest_hash = hashlib.sha256(args.manifest.read_bytes()).hexdigest()
    if plan.get("manifest_sha256") != manifest_hash:
        parser.error("resource plan was not created from this exact manifest")
    if plan.get("manifest_row_count") != len(rows) or len(plan.get("samples", [])) != len(rows):
        parser.error("resource plan row count does not match manifest")
    manifest_keys = [(r.get("study_id", ""), r.get("donor_id", ""), r.get("sample_id", "")) for r in rows]
    plan_keys = [(r.get("study_id", ""), r.get("donor_id", ""), r.get("sample_id", "")) for r in plan["samples"]]
    if manifest_keys != plan_keys:
        parser.error("resource plan sample order/set does not match manifest")
    output_ids = [output_id(i, row) for i, row in enumerate(rows)]
    if len(set(output_ids)) != len(output_ids):
        parser.error("generated output/checkpoint identifiers are not unique")
    env = plan["array_envelope"]
    out = args.out_dir.resolve()
    if out.exists() and any(out.iterdir()):
        parser.error("--out-dir must be absent or empty; refusing to overwrite")
    stages = [
        ("01_ingested", "00_manifest", True), ("02_qc", "01_ingested", False),
        ("03_harmonized", "02_qc", False), ("04_sketch", "03_harmonized", False),
        ("05_integrated", "04_sketch", False), ("06_annotated", "05_integrated", False),
        ("07_pseudobulk", "06_annotated", False), ("08_mechanism", "07_pseudobulk", False),
        ("09_release", "08_mechanism", False),
    ]
    for subdir in ("bin", "slurm", "logs", "partials"):
        (out / subdir).mkdir(parents=True, exist_ok=True)
    for stage, _, _ in stages:
        (out / "checkpoints" / stage).mkdir(parents=True, exist_ok=True)
        (out / "outputs" / stage).mkdir(parents=True, exist_ok=True)
    snapshot = out / "manifest.tsv"
    fields = list(rows[0].keys())
    with snapshot.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    with (out / "output_ids.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["row_index", "sample_id", "output_id"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({"row_index": i, "sample_id": row.get("sample_id", ""), "output_id": output_ids[i]} for i, row in enumerate(rows))
    shutil.copyfile(args.resource_plan, out / "resource_plan.json")
    runner = out / "bin/run_sample.sh"
    runner.write_text('''#!/usr/bin/env bash
set -euo pipefail
row_index="${1:?row index required}"
root="${2:?layout root required}"
manifest="$root/manifest.tsv"
IFS=$'\\t' read -r sample_id input_path safe_id < <(python3 - "$manifest" "$row_index" <<'PY'
import csv, hashlib, pathlib, re, sys
with pathlib.Path(sys.argv[1]).open(encoding="utf-8", newline="") as h:
    rows=list(csv.DictReader(h, delimiter="\\t"))
i=int(sys.argv[2])
if i < 0 or i >= len(rows): raise SystemExit("array index outside manifest")
r=rows[i]
values=[]
for key in ("sample_id", "input_path"):
    value=(r.get(key) or "").strip()
    if not value or any(x in value for x in ("\\n", "\\r", "\\t")): raise SystemExit(f"invalid {key}")
    values.append(value)
slug=re.sub(r"[^A-Za-z0-9._-]", "_", values[0])[:80] or "sample"
identity="\\0".join((str(i), r.get("study_id", ""), r.get("donor_id", ""), values[0], values[1]))
values.append(f"{i:06d}_{slug}_{hashlib.sha256(identity.encode()).hexdigest()[:12]}")
print("\\t".join(values))
PY
)
checkpoint="$root/checkpoints/01_ingested/${safe_id}.done"
output_dir="$root/outputs/01_ingested/$safe_id"
if [[ -s "$checkpoint" ]]; then echo "checkpoint exists: $sample_id"; exit 0; fi
: "${PIPELINE_DRIVER:?Set PIPELINE_DRIVER to a reviewed executable; this template does not run analysis by itself}"
if [[ ! -x "$PIPELINE_DRIVER" ]]; then echo "PIPELINE_DRIVER is not executable" >&2; exit 64; fi
run_tag="${SLURM_JOB_ID:-local}.${SLURM_ARRAY_TASK_ID:-$row_index}.$$"
temp_dir="$root/outputs/01_ingested/.${safe_id}.tmp.${run_tag}"
quarantine="$root/partials/${safe_id}.${run_tag}"
if [[ -e "$output_dir" ]]; then mv "$output_dir" "$root/partials/${safe_id}.prior.${run_tag}"; fi
mkdir -p "$temp_dir"
cleanup() {
  status=$?
  if [[ $status -ne 0 && -e "$temp_dir" ]]; then mv "$temp_dir" "$quarantine"; fi
  exit "$status"
}
trap cleanup EXIT
"$PIPELINE_DRIVER" --sample-id "$sample_id" --input "$input_path" --output "$temp_dir"
if [[ -z "$(find "$temp_dir" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "driver succeeded but produced no output" >&2; exit 65
fi
mv "$temp_dir" "$output_dir"
python3 - "$checkpoint.tmp.${run_tag}" "$sample_id" "$safe_id" "$input_path" "$PIPELINE_DRIVER" "$row_index" <<'PY'
import hashlib, json, os, pathlib, platform, sys
marker, sample, safe_id, input_path, driver, row_index=sys.argv[1:]
p=pathlib.Path(input_path)
if p.is_file():
    h=hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024), b""): h.update(block)
    input_hash="sha256:"+h.hexdigest()
else: input_hash="unavailable:path:"+input_path
payload={"sample_id":sample,"output_id":safe_id,"row_index":int(row_index),"status":"complete","exit_status":0,
 "input_path":input_path,"input_hash":input_hash,"driver":str(pathlib.Path(driver).resolve()),
 "driver_version":os.environ.get("PIPELINE_DRIVER_VERSION","unrecorded"),
 "command":[driver,"--sample-id",sample,"--input",input_path,"--output","<atomic-temp-dir>"],
 "python_version":platform.python_version(),"seed":os.environ.get("ANALYSIS_SEED","1"),
 "resources":{"slurm_job_id":os.environ.get("SLURM_JOB_ID","local"),"cpus":os.environ.get("SLURM_CPUS_PER_TASK","unrecorded"),"memory_per_node":os.environ.get("SLURM_MEM_PER_NODE","unrecorded")}}
pathlib.Path(marker).write_text(json.dumps(payload,indent=2)+"\\n",encoding="utf-8")
PY
mv "$checkpoint.tmp.${run_tag}" "$checkpoint"
trap - EXIT
''', encoding="utf-8")
    runner.chmod(runner.stat().st_mode | stat.S_IXUSR)
    sbatch = out / "slurm/01_ingest_array.sbatch"
    sbatch.write_text(f'''#!/usr/bin/env bash
#SBATCH --job-name=sc_ingest
#SBATCH --array=0-{len(rows)-1}%{args.max_concurrent}
#SBATCH --cpus-per-task={int(env['cpus'])}
#SBATCH --mem={int(env['memory_gb'])}G
#SBATCH --time={int(env['wall_hours']):02d}:00:00
#SBATCH --output={out}/logs/%x_%A_%a.out
#SBATCH --error={out}/logs/%x_%A_%a.err
set -euo pipefail
# REVIEW partition/account/QOS and pilot-derived resources before manual sbatch.
"{runner}" "$SLURM_ARRAY_TASK_ID" "{out}"
''', encoding="utf-8")
    dag = {"schema_version": 1, "stages": [
        {"id": stage, "depends_on": [dependency], "status": "executable_template" if executable else "placeholder_requires_reviewed_driver"}
        for stage, dependency, executable in stages
    ]}
    (out / "stage_dag.json").write_text(json.dumps(dag, indent=2) + "\n", encoding="utf-8")
    (out / "README.md").write_text(f'''# Generated single-cell job-array layout

Planning output only: nothing has been submitted. Review `resource_plan.json`, set a tested `PIPELINE_DRIVER`, and edit Slurm account/partition/QOS before manually submitting.

- Samples: {len(rows)}
- Maximum concurrent tasks: {args.max_concurrent}
- Each sample writes only to its own output and checkpoint path.
- Output identifiers include row index and a stable identity hash, preventing sanitized-name collisions.
- Failed or pre-existing partial outputs are moved under `partials/`; successful outputs are atomically renamed.
- A non-empty `.done` marker makes retries skip completed samples.
- The envelope is heuristic until small/median/large pilot benchmarks replace it.
- **Only `01_ingested` has an executable template.** Directories and DAG entries for stages 02-09 are placeholders and must not be run until stage-specific drivers are implemented and reviewed.
''', encoding="utf-8")
    print(str(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
