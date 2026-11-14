#!/usr/bin/env python3
"""Current C1 reliability and PGS aggregate reports; no fitting or old-format adapters."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from file_lock import flock as shared_flock

RELIABILITY_TABLES = ("reliability_status", "reliability_population_performance", "risk_coverage", "aurc_contrasts",
    "reliability_increment_50", "risk_standardized_selection", "support_trend", "reliability_prevalent_support_audit")


def reliability_tables(scope):
    """Use verified aggregates or this run's explicit unavailable-model status."""
    index=Path(scope)/"qc/C1_FINAL_INDEX.json"
    if index.is_file():
        latest=json.loads(index.read_text())
        unavailable={**latest.get("unavailable_methods",{}),**latest.get("failed_methods",{})}
        if "glmnet" in unavailable:
            # Preserve old fitted files, but never report their metrics as a
            # current result when this cohort cannot support the risk model.
            item=unavailable["glmnet"]
            status=pd.DataFrame([dict(head=head,status="unavailable_glmnet",
                reason=item.get("reason",item.get("status","unavailable"))) for head in
                ("P0_risk_only","P1_incident_support","P2_incident_plus_prevalent")])
            return {"reliability_status":status},dict(status="skipped_unavailable_glmnet")
    directory=Path(scope)/"qc/reliability"
    marker=directory/"COMPLETE.json"
    if not marker.is_file():return {},{}
    try:
        receipt=json.loads(marker.read_text());files=receipt.get("files",{})
        contract=directory/"reliability_contract.json"
        requested=[contract,*[directory/(name+".tsv") for name in RELIABILITY_TABLES if (directory/(name+".tsv")).is_file()]]
        for file in requested:
            if files.get(file.name)!=hashlib.sha256(file.read_bytes()).hexdigest():
                raise ValueError(f"Reliability aggregate integrity failed: {file}")
        meta=json.loads(contract.read_text())
        if meta.get("version")!="LE8-C1-reliability-2026-10-26-v2":return {},{}
        tables={name:pd.read_csv(directory/(name+".tsv"),sep="\t") for name in RELIABILITY_TABLES if (directory/(name+".tsv")).is_file()}
        return tables,meta
    except (OSError,json.JSONDecodeError,pd.errors.EmptyDataError) as error:
        raise ValueError(f"Cannot read reliability results: {directory}") from error



def reliability_summary_row(trait, tables, contract):
    row=dict(trait=trait,discovery_or_replication=contract.get("discovery_or_replication","unavailable"))
    status=tables.get("reliability_status",pd.DataFrame())
    for prefix,name in (("P0","P0_risk_only"),("P1","P1_incident_support"),("P2","P2_incident_plus_prevalent")):
        values=status[status["head"].eq(name)] if "head" in status else pd.DataFrame()
        row[prefix+"_status"]=values.iloc[0]["status"] if len(values) else "unavailable"
    primary=[]
    def extract(table,model,reference,metric,prefix):
        frame=tables.get(table,pd.DataFrame());result=None
        if {"model","reference","metric"}<=set(frame):
            z=frame[frame.model.eq(model)&frame.reference.eq(reference)&frame.metric.eq(metric)]
            if len(z):result=z.iloc[0]
        row[prefix]=float(result["difference"]) if result is not None else np.nan
        row[prefix+"_CI"]=f"[{result['lower']:.6g}, {result['upper']:.6g}]" if result is not None and np.isfinite(result.get("lower",np.nan)) and np.isfinite(result.get("upper",np.nan)) else "unavailable"
        return result
    for metric in ("Brier","LogLoss"):
        primary.append(extract("aurc_contrasts","P1_incident_support","P0_risk_only","AURC_"+metric,"P1_minus_P0_AURC_"+metric))
        primary.append(extract("reliability_increment_50","P1_incident_support","P0_risk_only",metric,"P1_minus_P0_50pct_"+metric))
    extract("aurc_contrasts","P2_incident_plus_prevalent","P1_incident_support","AURC_Brier","P2_minus_P1_AURC_Brier")
    trend=tables.get("support_trend",pd.DataFrame())
    for support,prefix in (("eligible_match_count","twin_count_Brier_slope"),("family_ess","family_ESS_Brier_slope")):
        z=trend[trend.support_metric.eq(support)&trend.error_metric.eq("Brier")] if {"support_metric","error_metric"}<=set(trend) else pd.DataFrame()
        row[prefix]=float(z.iloc[0].slope) if len(z) else np.nan
        row[prefix+"_CI"]=f"[{z.iloc[0].lower:.6g}, {z.iloc[0].upper:.6g}]" if len(z) and np.isfinite(z.iloc[0].lower) and np.isfinite(z.iloc[0].upper) else "unavailable"
    if contract.get("status")!="completed_primary_P1_vs_P0":row["primary_conclusion"]="primary_unavailable"
    elif any(r is None or r.get("inference_status")!="evaluated" for r in primary):row["primary_conclusion"]="insufficient_primary_inference"
    elif all(r["upper"]<0 for r in primary):row["primary_conclusion"]="lower_error_supported_across_primary_contrasts"
    else:row["primary_conclusion"]="primary_support_increment_not_consistently_established"
    return row



def write_pgs_cross_trait_summary(root, layer):
    """Summarize saved C1 PGS scans, including traits outside the current batch.

    The run/skip switch controls fitting only. Missing/empty/all-NA scans are
    availability records, never tested features or negative findings.
    """
    import os, tempfile, fcntl
    from final import PGS_SCANS, PGS_SUMMARY_COLUMNS, saved_pgs_tables
    if layer not in {"prot", "met"}:
        raise ValueError("PGS summary layer must be prot or met")
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    from outputs import report_directory, summary_lock_directory
    output = report_directory(root)
    output.mkdir(parents=True, exist_ok=True)
    lockdir = summary_lock_directory(root)
    lockdir.mkdir(parents=True, exist_ok=True)
    def read(path):
        if not path.is_file() or not path.stat().st_size:
            return pd.DataFrame()
        try:
            return pd.read_csv(path)
        except pd.errors.EmptyDataError:
            return pd.DataFrame()
    def publish(path, write):
        with tempfile.NamedTemporaryFile(dir=output, prefix=".pgs-summary-", delete=False) as handle:
            temp = Path(handle.name)
        try:
            write(temp)
            os.replace(temp, path)
        finally:
            temp.unlink(missing_ok=True)
    with (lockdir / ("pgs_summary_" + layer + ".lock")).open("a") as lock:
        shared_flock(lock, fcntl.LOCK_EX)
        collected = {name: [] for name in ("summary", "status", "associations")}
        prefix = "pwas" if layer == "prot" else "mwas"
        for scope in sorted(root.glob("*/" + layer + "/c1_correlate")):
            scans = {}
            for role in PGS_SCANS:
                path = scope / (prefix + "_" + role + ".csv")
                d = read(path)
                if len(d):
                    d = d.assign(source_file=str(path.relative_to(root)))
                scans[role] = d
            saved = saved_pgs_tables(scans, read(scope / "c1.pgs_status.csv"))
            for name, d in saved.items():
                collected[name].append(d.assign(Y=scope.parent.parent.name, layer=layer))
        empty_columns = dict(summary=PGS_SUMMARY_COLUMNS,
            status=["status", "declared_status", "analyses_summarized", "detail"],
            associations=["term", "analysis", "cohort_scope", "p.value", "FDR", "source_file"])
        result = {}
        for name, frames in collected.items():
            nonempty = [d for d in frames if len(d)]
            d = pd.concat(nonempty, ignore_index=True) if nonempty else pd.DataFrame(columns=[*empty_columns[name], "Y", "layer"])
            result[name] = d
            target = output / (f"c1_pgs_cross_trait_{name}.{layer}.tsv")
            publish(target, lambda path, frame=d: frame.to_csv(path, sep="\t", index=False))
        body = '<!doctype html><meta charset="utf-8"><title>C1 saved PGS summary</title><h1>C1 saved PGS analyses</h1>'
        body += '<p>Saved estimates and original within-scan FDR; no refitting. Missing PGS is omitted, not a negative result. Original cohort and model definitions apply.</p>'
        body += result["status"].to_html(index=False, escape=True)
        if len(result["summary"]):
            body += '<h2>Completed analyses</h2>' + result["summary"].to_html(index=False, escape=True)
        target = output / f"c1_pgs_cross_trait_summary.{layer}.html"
        publish(target, lambda path: path.write_text(body, encoding="utf-8"))
        print(f"[LE8 C1] PGS summary: {len(result['summary'])} saved analyses retained; missing PGS skipped; {target}", flush=True)
        return result



def write_reliability_cross_trait_summary(root, layer):
    import os,tempfile,fcntl
    from outputs import report_directory, summary_lock_directory
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    output=report_directory(root);output.mkdir(parents=True,exist_ok=True)
    lockdir=summary_lock_directory(root);lockdir.mkdir(parents=True,exist_ok=True)
    with (lockdir/("reliability_summary_"+layer+".lock")).open("a") as lock:
        shared_flock(lock,fcntl.LOCK_EX)
        rows=[]
        for scope in sorted(root.glob("*/"+layer+"/c1_correlate")):
            tables,contract=reliability_tables(scope)
            if tables:rows.append(reliability_summary_row(scope.parent.parent.name,tables,contract))
        frame=pd.DataFrame(rows)
        if not rows:return frame
        target=output/("c1_reliability_cross_trait_summary."+layer+".tsv")
        with tempfile.NamedTemporaryFile(dir=output,prefix=".reliability-summary-",delete=False) as handle:temp=Path(handle.name)
        try:frame.to_csv(temp,sep="\t",index=False);os.replace(temp,target)
        finally:temp.unlink(missing_ok=True)
        return frame



def write_reliability_report(root, traits, layers, out):
    """Primary report for current frozen-model reliability."""
    import html
    out=Path(out);root=Path(root);collected={};summaries=[]
    for trait in traits:
        for layer in layers:
            tables,contract=reliability_tables(root/trait/layer/"c1_correlate")
            if not tables:continue
            summaries.append(dict(reliability_summary_row(trait,tables,contract),layer=layer))
            for name,frame in tables.items():collected.setdefault(name,[]).append(frame.assign(trait=trait,layer=layer))
    if not summaries:return False
    out.mkdir(parents=True,exist_ok=True)
    merged={k:pd.concat(v,ignore_index=True) for k,v in collected.items()}
    merged["cross_trait_summary"]=pd.DataFrame(summaries)
    with pd.ExcelWriter(out/"C1_FINAL_reliability_tables.xlsx",engine="openpyxl") as writer:
        for name,frame in merged.items():
            frame.to_excel(writer,sheet_name=name[:31],index=False)
            frame.to_csv(out/("c1."+name+".tsv"),sep="\t",index=False)
    body=['<!doctype html><html lang="en"><meta charset="utf-8"><title>C1 FINAL reliability</title>',
        '<style>body{font:15px system-ui;margin:32px;line-height:1.5}table{border-collapse:collapse;font-size:12px}td,th{padding:6px;border:1px solid #ddd}.table{overflow:auto}</style>',
        '<h1>C1 FINAL: reliability of frozen glmnet risk</h1>',
        '<p>Primary evidence: direct P1−P0 selected-set losses at 50% coverage and AURC on a common target grid. Negative differences favor molecular support. Intervals use paired family bootstrap conditional on frozen models and heads. Selected-vs-nonselected tables are supplementary.</p>']
    def section(title,name,head=None):
        frame=merged.get(name,pd.DataFrame())
        if head is not None and "head" in frame:frame=frame[frame["head"].eq(head)]
        if head is not None and "model" in frame:frame=frame[frame.model.eq(head)]
        body.extend(["<h2>"+html.escape(title)+"</h2>",'<div class="table">'+(frame.to_html(index=False,escape=True) if len(frame) else '<p>Unavailable; see head status.</p>')+'</div>'])
    section("1. All-person frozen glmnet performance","reliability_population_performance")
    section("2. P0 risk-only reliability","risk_coverage","P0_risk_only")
    section("3. P1 risk and incident molecular support","risk_coverage","P1_incident_support")
    section("4. Primary P1−P0: direct 50% contrast","reliability_increment_50","P1_incident_support")
    section("Primary and secondary target-grid AURC","aurc_contrasts")
    section("5. P2 risk, incident and prevalent support","risk_coverage","P2_incident_plus_prevalent")
    section("6. Secondary P2−P1: direct 50% contrast","reliability_increment_50","P2_incident_plus_prevalent")
    section("7. Continuous support-error trend, adjusted for frozen risk strata","support_trend")
    section("8. Discovery / replication and primary conclusions","cross_trait_summary")
    section("Availability of every head","reliability_status")
    body.append('<p>HTN is discovery. COPD is locked replication only with documented prior non-inspection; otherwise it remains an unverified replication candidate. No K/PCA/radius/architecture/coverage tuning using test outcomes.</p></html>')
    (out/"C1_FINAL_reliability.html").write_text("\n".join(body))
    (out/"C1_FINAL_REPORT.json").write_text(json.dumps(dict(primary_report="C1_FINAL_reliability.html",source="verified FINAL aggregate tables",risk_model_refit=False,scopes=summaries),indent=2,default=str)+"\n")
    return True
