#!/usr/bin/env python3
"""Normalize aggregate estimates, evaluate Bayes factors, and collect paper outputs.

Only aggregate analysis exports enter the web report. Participant RDS files are
never read here. Native manuscript tables and plots retain their source layout.
"""
from __future__ import annotations
import argparse
import csv
import html
import json
import math
import os
from pathlib import Path
import re
import shutil
import tempfile

BASE = Path(__file__).resolve().parents[1]

def rows(path):
    path = Path(path)
    if not path.is_file(): return []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="," if path.suffix.lower()==".csv" else "\t"))

def write_rows(path, records, columns=None):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    if columns is None:
        columns=list(dict.fromkeys(k for row in records for k in row))
    with tempfile.NamedTemporaryFile("w",encoding="utf-8",newline="",dir=path.parent,delete=False) as handle:
        tmp=Path(handle.name)
        writer=csv.DictWriter(handle,fieldnames=columns,delimiter="\t",extrasaction="ignore")
        writer.writeheader()
        for row in records:
            writer.writerow({k: ("NA" if v is None or isinstance(v,float) and not math.isfinite(v) else v) for k,v in row.items()})
    os.replace(tmp,path)

def choose(row,*names,default=""):
    for key in names:
        if row.get(key) not in (None,"","NA","NaN","nan"): return row[key]
    return default

def number(value):
    try:
        x=float(value)
        return x if math.isfinite(x) else None
    except (ValueError,TypeError): return None

def normalise(row,study,cohort,source):
    estimate=number(choose(row,"estimate","HR","OR","Effect"))
    lo=number(choose(row,"lower","conf.low","LCI","LCI95","CI_low"))
    hi=number(choose(row,"upper","conf.high","HCI","UCI","UCI95","CI_high"))
    beta=number(choose(row,"beta","logHR","logOR"))
    se=number(choose(row,"se","std.error","SE"))
    if beta is None and estimate is not None and estimate>0: beta=math.log(estimate)
    # Exact coefficient SE is required for NHANES. Its t-based CIs cannot
    # universally be inverted with 1.96. Omit BF if the native SE is absent.
    se_source="coefficient" if se is not None else "unavailable"
    if se is None and cohort!="nhanes" and lo and hi and 0<lo<hi:
        se=(math.log(hi)-math.log(lo))/(2*1.959963984540054)
        se_source="normal_CI_reconstruction"
    if estimate is None and beta is not None and abs(beta)<700: estimate=math.exp(beta)
    if lo is None and beta is not None and se is not None and abs(beta-1.959964*se)<700: lo=math.exp(beta-1.959964*se)
    if hi is None and beta is not None and se is not None and abs(beta+1.959964*se)<700: hi=math.exp(beta+1.959964*se)
    status=choose(row,"status","error",default="ok")
    if estimate is None: status="not_estimable" if status=="ok" else status
    return dict(study=study,cohort=cohort,family=choose(row,"family",default="manuscript_primary"),
      outcome=choose(row,"outcome","Outcome","Outcome_label",default="death" if cohort in ("nhanes","chns") and study=="maha" else ""),
      diet=choose(row,"diet","Diet","Score","score","Exposure"),
      model=choose(row,"model","Model","Model_Type",default="primary"),
      scale=choose(row,"scale","Scale","Contrast","contrast",default="SD" if study=="maha" else ""),
      effect_type=choose(row,"effect_type","Effect_measure","effect",default="OR" if study=="dose" and cohort=="nhanes" else "HR"),
      beta=beta,se=se,estimate=estimate,lower=lo,upper=hi,
      p=number(choose(row,"p","P","p.value","p_value")),
      fdr=number(choose(row,"fdr","FDR","P_FDR")),
      n=number(choose(row,"n","N","N_total")),events=number(choose(row,"events","Events","N_event")),
      df=number(choose(row,"df","df_residual")),increment=number(choose(row,"increment","Increment")),
      ci_method=choose(row,"ci_method",default="normal" if cohort!="nhanes" else "survey"),se_source=se_source,
      adjustment=choose(row,"adjustment","Covariates"),status=status,source=source)

def collect_effects(out):
    records=[]
    for path in sorted((out/"assoc").glob("*.tsv")):
        if path.name in {"effects.tsv","charls.failures.tsv"} or ".failures." in path.name: continue
        parts=path.stem.split(".")
        if len(parts)<2: continue
        cohort,study=parts[0],parts[1]
        if cohort not in {"ukb","nhanes","chns","charls"}: continue
        for row in rows(path): records.append(normalise(row,study,cohort,str(path.relative_to(out))))
    for cohort in ("ukb","nhanes"):
        path=out/"dose"/cohort/"effects.tsv"
        # A normalized adapter may already have copied to assoc/<cohort>.dose.tsv.
        if (out/"assoc"/f"{cohort}.dose.tsv").exists(): continue
        for row in rows(path): records.append(normalise(row,"dose",cohort,str(path.relative_to(out))))
    write_rows(out/"assoc/effects.tsv",records,columns=list(normalise({},"","","").keys()))
    return records

def bf_row(beta,se,prior95=1.10,margin=1.05):
    """Normal-estimate likelihood, normal H1 prior centred at zero; BF01=P(D|H0)/P(D|H1)."""
    if not (math.isfinite(beta) and math.isfinite(se) and se>0 and prior95>1 and margin>1):
        raise ValueError("Bayes requires finite beta, positive SE, prior95>1 and margin>1")
    tau=math.log(prior95)/1.96; v=se*se; w=tau*tau
    log_bf01=.5*math.log1p(w/v)-.5*(beta/se)**2*(w/(v+w))
    def bounded_exp(x):
        return math.exp(x) if -745<x<709 else (0.0 if x<=-745 else float("inf"))
    m=math.log(margin)
    # max of two one-sided tests; 90% CI must lie strictly inside +/- log(margin).
    p1=.5*math.erfc((beta+m)/se/math.sqrt(2))
    p2=.5*math.erfc((m-beta)/se/math.sqrt(2))
    low=beta-1.6448536269514722*se; high=beta+1.6448536269514722*se
    return dict(prior_95_ratio=prior95,prior_sd_log_effect=tau,margin_ratio=margin,
      log_BF01=log_bf01,BF01=bounded_exp(log_bf01),BF10=bounded_exp(-log_bf01),
      evidence="point_null" if log_bf01>=math.log(3) else "alternative" if log_bf01<=-math.log(3) else "inconclusive",
      TOST_p=max(p1,p2),equivalent=low>-m and high<m,
      CI90_lower=bounded_exp(low),CI90_upper=bounded_exp(high))

def bayes(out,cohorts,studies,priors,margins):
    effects=collect_effects(out); result=[]
    for row in effects:
        if row["cohort"] not in cohorts or row["study"] not in studies: continue
        b,se=row["beta"],row["se"]
        if b is None or se is None or se<=0 or row["status"]!="ok": continue
        for prior in priors:
            for margin in margins:
                result.append({**row,**bf_row(b,se,prior,margin),
                    "method":"Normal likelihood approximation; zero-centred normal prior",
                    "note":"Survey coefficient normal approximation" if row["cohort"]=="nhanes" else ""})
    path=out/"bayes/effects.tsv"
    # Preserve completed results for cohorts/studies not selected in this call.
    old=[r for r in rows(path) if r.get("cohort") not in cohorts or r.get("study") not in studies]
    write_rows(path,old+result,columns=list(normalise({},"","","").keys())+list(bf_row(0,.1).keys())+["method","note"])
    return result

def is_aggregate_table(path):
    """Exclude participant/cache/model payloads; accept manuscript result tables only."""
    if path.suffix.lower() not in {".csv",".tsv",".txt",".xlsx"}: return False
    if any(p.startswith(".") for p in path.parts): return False
    if re.search(r"(participant|individual|dat[._-]|cohort[._-]data|train|prediction|^scores?$)",path.stem,re.I): return False
    if path.suffix.lower()==".xlsx": return True  # Publication exports only, from a scoped source root.
    try:
        # Inspect only the header and one row. A participant export can be many
        # gigabytes and must never be loaded just to decide whether to omit it.
        with path.open(encoding="utf-8-sig",newline="") as handle:
            reader=csv.DictReader(handle,delimiter="," if path.suffix.lower()==".csv" else "\t")
            fields=reader.fieldnames or []
            forbidden={"eid","id","iid","seqn","idind","participant_id","patient_id"}
            if any(k.lower().strip() in forbidden for k in fields): return False
            return next(reader,None) is not None
    except (UnicodeError,csv.Error): return False

def collect_final(out,studies,registry=None):
    registry=registry or BASE/"config/papers.json"
    defs=json.loads(Path(registry).read_text())
    all_entries=[]
    for study in studies:
        if study not in {"maha","dose"}: continue
        dest=out/"final"/study; dest.mkdir(parents=True,exist_ok=True)
        source_roots=[out/study]
        # MAHA's renderer writes directly to final/maha; preserve those native composites.
        candidates=[]
        for source in source_roots:
            if source.is_dir(): candidates += sorted(source.rglob("*"))
        old_index=rows(dest/"index.tsv")
        emitted=set(); entries=[]
        for path in candidates:
            if not path.is_file() or any(p.startswith(".") for p in path.relative_to(out).parts): continue
            suffix=path.suffix.lower()
            if suffix not in {".png",".pdf",".svg"} and not is_aggregate_table(path): continue
            rel=path.relative_to(out/study)
            # A flat article directory; cohort/stage prefixes prevent basename collisions.
            name="__".join(rel.parts)
            target=dest/name
            if len(name)>180:
                import hashlib
                name=path.stem[:110]+"_"+hashlib.sha256(str(rel).encode()).hexdigest()[:10]+suffix
                target=dest/name
            shutil.copy2(path,target); emitted.add(name)
            entries.append(dict(study=study,kind="figure" if suffix in {".png",".pdf",".svg"} else "table",
                title=path.stem,file=str(target.relative_to(out)),source=str(path.relative_to(out)),status="available"))
        # Keep files made by the native publication renderer, avoiding our previous copies.
        previous={Path(x.get("file","")).name for x in old_index if x.get("source","").startswith(study+"/")}
        for path in sorted(dest.rglob("*")):
            if not path.is_file() or path.name in emitted or path.name in previous or path.name=="index.tsv": continue
            if path.suffix.lower() not in {".png",".pdf",".svg",".xlsx",".csv",".tsv"}: continue
            if path.suffix.lower() in {".csv",".tsv"} and not is_aggregate_table(path): continue
            entries.append(dict(study=study,kind="figure" if path.suffix.lower() in {".png",".pdf",".svg"} else "table",
                title=path.stem,file=str(path.relative_to(out)),source="native_publication_renderer",status="available"))
        # Remove only files owned by our preceding index which are no longer sourced.
        for name in previous-emitted:
            target=dest/name
            if target.is_file(): target.unlink()
        publication={r.get("item"):r for r in rows(dest/"publication.status.tsv")}
        composites={r.get("figure"):r for r in rows(dest/"composite_status.tsv")}
        for item in defs.get(study,[]):
            # Manuscript figures have exact registered names. A source forest,
            # sensitivity plot or MRI panel is never mistaken for a main figure.
            path=dest/item["file"]
            exists=path.is_file() and path.stat().st_size>0
            state="available" if exists else "not_generated"
            detail=item.get("detail","")
            if study=="dose":
                record=composites.get(item["id"],{})
                if record.get("status")=="partial" and exists: state="partial"
                elif record.get("status")=="input_missing": state="not_generated";exists=False
                elif exists and record.get("status")!="complete": state="partial"
                if record.get("missing_panels"): detail+=" Missing panels: "+record["missing_panels"]
            else:
                record=publication.get(item.get("cohort"),{})
                if exists and record.get("status")!="complete":
                    state="partial"
                    detail+=" Cohort publication build is incomplete; this figure is available for inspection."
                if record.get("reason"):detail+=" Build note: "+record["reason"]
            entries.append(dict(study=study,kind="manuscript_item",title=item["title"],
                file=str(path.relative_to(out)) if exists else "",source=item["stage"],
                status=state,manuscript_id=item["id"],detail=detail.strip()))
        write_rows(dest/"index.tsv",entries)
        all_entries.extend(entries)
    # Consolidated index includes previously generated studies too.
    merged=[]
    for study in ("maha","dose"): merged += rows(out/"final"/study/"index.tsv")
    write_rows(out/"final/index.tsv",merged)
    return all_entries

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("action",choices=["collect","bayes","final"])
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--cohorts",default="ukb,nhanes,chns,charls")
    p.add_argument("--studies",default="maha,dose,charls")
    p.add_argument("--priors",default="1.05,1.10,1.20")
    p.add_argument("--margins",default="1.05,1.075,1.10")
    a=p.parse_args(argv); out=a.output.resolve()
    if a.action=="collect": result=collect_effects(out)
    elif a.action=="bayes": result=bayes(out,a.cohorts.split(","),a.studies.split(","),
            [float(x) for x in a.priors.split(",")],[float(x) for x in a.margins.split(",")])
    else: result=collect_final(out,a.studies.split(","))
    print(f"{a.action}: {len(result)} aggregate records")
if __name__=="__main__": main()
