#!/usr/bin/env python3
"""Validate the minimum GWAS study manifest contract."""
import argparse, csv, json, re, sys
from collections import Counter
from datetime import date, datetime, timezone
from pathlib import Path

REQ="study_id source release trait trait_definition role ancestry genome_build n_total n_cases n_controls unit sample_overlap access_date".split()
BUILDS={"GRCH37","HG19","GRCH38","HG38"}; ROLES={"discovery","replication","exposure","outcome","qtl","validation"}
def missing(x): return not x.strip() or x.strip().lower() in {"na","n/a","null"}
def num(x):
    try:return float(x)
    except ValueError:return None
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("input",type=Path); ap.add_argument("--report",type=Path,required=True); a=ap.parse_args()
    try:
        with open(a.input,encoding="utf-8-sig",newline="") as f:
            sample=f.read(32768); f.seek(0); d=csv.Sniffer().sniff(sample,delimiters="\t,").delimiter; r=csv.DictReader(f,delimiter=d)
            cols=r.fieldnames or []; errors=[]; warnings=[]; examples=[]; counts=Counter(); ids=set()
            miss=[x for x in REQ if x not in cols]
            if miss: errors.append({"code":"missing_required_columns","fields":miss})
            for line,row in enumerate(r,start=2):
                counts["rows"]+=1; er=[]; wa=[]; sid=row.get("study_id","").strip()
                for k in REQ:
                    if k in cols and missing(row.get(k,"")): er.append(f"missing:{k}")
                if sid in ids: er.append("duplicate:study_id")
                ids.add(sid)
                build=row.get("genome_build","").strip().upper()
                if build and build not in BUILDS: er.append("invalid:genome_build")
                role=row.get("role","").strip().lower()
                if role and role not in ROLES: wa.append("unrecognized:role")
                nt=num(row.get("n_total","")); nc=num(row.get("n_cases","")); nn=num(row.get("n_controls",""))
                if any(v is not None and v<0 for v in (nt,nc,nn)): er.append("invalid:sample_size")
                if None not in (nt,nc,nn) and nc+nn>nt: er.append("inconsistent:case_control_total")
                try: date.fromisoformat(row.get("access_date",""))
                except ValueError: er.append("invalid:access_date")
                overlap=row.get("sample_overlap","").strip().lower()
                if overlap in {"unknown","unclear","not assessed"}: wa.append("unresolved:sample_overlap")
                counts["error_rows" if er else "valid_rows"]+=1; counts["warning_rows"]+=bool(wa)
                if er or wa: examples.append({"line":line,"study_id":sid,"errors":er,"warnings":wa})
            if not counts["rows"]: errors.append({"code":"no_data_rows"})
        report={"schema_version":"1.0","generated_at":datetime.now(timezone.utc).isoformat(),"input":str(a.input),"columns":cols,"counts":dict(counts),"errors":errors,"warnings":warnings,"row_findings":examples}
        a.report.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
        print(f"rows={counts['rows']} valid={counts['valid_rows']} error_rows={counts['error_rows']} warning_rows={counts['warning_rows']} schema_errors={len(errors)}")
        return 1 if errors or counts["error_rows"] else 0
    except Exception as e: print(f"ERROR: {e}",file=sys.stderr); return 2
if __name__=="__main__": raise SystemExit(main())
