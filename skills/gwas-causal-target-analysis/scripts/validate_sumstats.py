#!/usr/bin/env python3
"""Validate and normalize a delimited GWAS summary-statistics table."""
import argparse, csv, gzip, json, math, re, sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ALIASES = {
 "chr":["chr","chrom","chromosome"], "pos":["pos","position","bp","base_pair_location"],
 "rsid":["rsid","snp","variant_id","markername"], "effect_allele":["effect_allele","ea","a1","alt"],
 "other_allele":["other_allele","oa","nea","a2","ref"], "beta":["beta","effect","estimate","b"],
 "or":["or","odds_ratio"], "se":["se","standard_error","stderr"], "p":["p","pval","p_value","pvalue"],
 "eaf":["eaf","effect_allele_frequency","af","a1freq"], "n":["n","n_total","samplesize"], "info":["info","imputation_info","rsq"]}
REQUIRED = ["chr","pos","effect_allele","other_allele","se","p"]
ALLELE_RE = re.compile(r"^[ACGT]+$")
MULTIALLELIC_RE = re.compile(r"[,/|]")

def op(path): return gzip.open(path, "rt", encoding="utf-8-sig", newline="") if str(path).endswith(".gz") else open(path, encoding="utf-8-sig", newline="")
def clean(x): return re.sub(r"[^a-z0-9]+", "_", x.strip().lower()).strip("_")
def number(x):
    try:
        v=float(x); return v if math.isfinite(v) else None
    except (ValueError, TypeError): return None
def allele(x): return x.strip().upper()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("input", type=Path); ap.add_argument("--report", type=Path, required=True); ap.add_argument("--normalized", type=Path); ap.add_argument("--delimiter", choices=["tab","comma"]); ap.add_argument("--max-examples", type=int, default=20); a=ap.parse_args()
    try:
        with op(a.input) as fh:
            sample=fh.read(65536); fh.seek(0)
            delim="\t" if a.delimiter=="tab" else "," if a.delimiter=="comma" else csv.Sniffer().sniff(sample, delimiters="\t,").delimiter
            reader=csv.DictReader(fh, delimiter=delim)
            if not reader.fieldnames: raise ValueError("missing header")
            lookup={clean(x):x for x in reader.fieldnames}; mapping={}
            for canon,names in ALIASES.items():
                hits=[lookup[n] for n in names if n in lookup]
                if hits: mapping[canon]=hits[0]
            errors=[]; warnings=[]
            missing=[x for x in REQUIRED if x not in mapping]
            if "beta" not in mapping and "or" not in mapping: missing.append("beta_or_or")
            if missing: errors.append({"code":"missing_required_columns","fields":missing})
            counts=Counter(); examples=[]; seen=set(); rows=[]
            for i,row in enumerate(reader, start=2):
                counts["rows"]+=1; reasons=[]
                get=lambda k: row.get(mapping.get(k,""),"").strip()
                ch=get("chr").upper().removeprefix("CHR"); pv=number(get("p")); pos=number(get("pos")); se=number(get("se")); ea=allele(get("effect_allele")); oa=allele(get("other_allele")); ev=number(get("beta")); ov=number(get("or")); eaf=number(get("eaf")); info=number(get("info")); n=number(get("n"))
                if mapping and (ch not in {str(x) for x in range(1,23)}|{"X","Y","MT","M"}): reasons.append("invalid_chr")
                if pos is None or pos<=0 or int(pos)!=pos: reasons.append("invalid_pos")
                if MULTIALLELIC_RE.search(ea) or MULTIALLELIC_RE.search(oa):
                    reasons.append("unsplit_multiallelic")
                elif not ALLELE_RE.fullmatch(ea) or not ALLELE_RE.fullmatch(oa) or ea==oa:
                    reasons.append("invalid_alleles")
                if pv is None or not (0<=pv<=1): reasons.append("invalid_p")
                if se is None or se<=0: reasons.append("invalid_se")
                if "beta" in mapping and ev is None: reasons.append("invalid_beta")
                if "or" in mapping and (ov is None or ov<=0): reasons.append("invalid_or")
                if "eaf" in mapping and (eaf is None or not 0<=eaf<=1): reasons.append("invalid_eaf")
                if "info" in mapping and (info is None or not 0<=info<=1): reasons.append("invalid_info")
                if "n" in mapping and (n is None or n<=0): reasons.append("invalid_n")
                key=(ch, get("pos"), ea, oa)
                if key in seen: reasons.append("duplicate_variant")
                seen.add(key)
                if {ea,oa} in ({"A","T"},{"C","G"}): reasons.append("palindromic")
                status="FAIL" if any(x.startswith("invalid") or x in {"duplicate_variant","unsplit_multiallelic"} for x in reasons) else "WARN" if reasons else "PASS"
                counts[status.lower()]+=1
                for reason in reasons: counts[reason]+=1
                if reasons and len(examples)<a.max_examples: examples.append({"line":i,"reasons":reasons})
                if a.normalized:
                    rows.append({k:get(k) for k in ["chr","pos","rsid","effect_allele","other_allele","beta","or","se","p","eaf","n","info"]}|{"qc_status":status,"qc_reason":";".join(reasons)})
            if counts["rows"]==0: errors.append({"code":"no_data_rows"})
        report={"schema_version":"1.0","generated_at":datetime.now(timezone.utc).isoformat(),"input":str(a.input),"delimiter":"tab" if delim=="\t" else "comma","columns":reader.fieldnames,"column_mapping":mapping,"counts":dict(counts),"errors":errors,"warnings":warnings,"examples":examples,"interpretation":"Schema and row-level QC only; not allele harmonization or genome-build verification."}
        a.report.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
        if a.normalized:
            fields=["chr","pos","rsid","effect_allele","other_allele","beta","or","se","p","eaf","n","info","qc_status","qc_reason"]
            with open(a.normalized,"w",encoding="utf-8",newline="") as out:
                w=csv.DictWriter(out,fieldnames=fields,delimiter="\t"); w.writeheader(); w.writerows(rows)
        print(f"rows={counts['rows']} pass={counts['pass']} warn={counts['warn']} fail={counts['fail']} errors={len(errors)}")
        return 1 if errors or counts["fail"] else 0
    except Exception as e:
        print(f"ERROR: {e}",file=sys.stderr); return 2
if __name__=="__main__": raise SystemExit(main())
