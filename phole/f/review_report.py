#!/usr/bin/env python3
"""Build the offline cross-check view from checksum-matched audit supplements.

The authoritative cross-check JSON is never changed by this rendering step.
This preserves the input hashes used by specificity and feasibility analyses.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def validate_files(folder,hashes):
    for name,expected in hashes.items():
        if Path(name).name!=name:raise ValueError("Unsafe input filename: "+name)
        path=folder/name
        if not path.is_file() or sha(path)!=expected:raise ValueError("Stale or missing audit input: "+str(path))


PAIR_DISPLAY_FIELDS = set((
    "PAIR_ID HYPOTHESIS_ID PROT PRIMARY PRIMARY_ASSAY ALTERNATIVE_PROT REFERENCE REFERENCE_ASSAY "
    "TIER PRIOR_PRIORITY_CASE SCREEN_STRATUM SCREEN_STRATA N_ASSAY_HYPOTHESES N_PRIMARY_ASSAYS "
    "N_REFERENCE_COHORTS N_REFERENCE_ASSAYS REFERENCE_COHORTS PROFILE_METRIC_SCOPE "
    "FULL_SHARED_REGIONS FULL_EXACT_RSID FULL_JACCARD FULL_INTENDED_MAX_JACCARD "
    "REGION_HELDOUT_SHARED_REGIONS REGION_HELDOUT_NONHOTSPOT_BOTH "
    "REGION_HELDOUT_PREFERRED_TO_ALL_INTENDED EXACT_ANCHOR_MULTILOCUS"
).split())


def fingerprint_display(folder,cross_dir,specificity_dir):
    """Verify supplements, then embed all pair summaries and a bounded detail set.

    Selection limits browser payload size only; it never changes scientific tiers
    or the complete TSV outputs. Every exact-anchor multi-region case is included.
    """
    folder=Path(folder)
    summary=json.loads((folder/"fingerprint-dossiers.summary.json").read_text())
    validate_files(folder,summary["output_sha256"])
    data=json.loads((folder/"fingerprint-dossiers.data.json").read_text())
    if data.get("summary")!={k:v for k,v in summary.items() if k!="output_sha256"}:
        raise ValueError("Fingerprint data and summary disagree")
    for name,expected in summary["input_sha256"].items():
        if name=="protein_map":
            provenance=json.loads((cross_dir/"provenance.json").read_text())
            if provenance.get("mapping_sha256")!=expected:
                raise ValueError("Fingerprint protein map does not match this cross-check")
            continue
        parts=Path(name).parts
        if len(parts)==2 and parts[0]=="specificity":path=specificity_dir/parts[1]
        elif len(parts)==1:path=cross_dir/name
        else:raise ValueError("Unsafe fingerprint input filename: "+name)
        if not path.is_file() or sha(path)!=expected:
            raise ValueError("Stale or missing fingerprint input: "+str(path))
    selected=[h for h in data.get("hypotheses",[]) if h.get("PRIOR_PRIORITY_CASE")==1
              or h.get("EXACT_ANCHOR_MULTILOCUS")==1 or h.get("TIER")=="T1_distinctive_heldout"]
    ids={h["HYPOTHESIS_ID"] for h in selected}
    display={"summary":summary,"pairs":[{k:v for k,v in r.items() if k in PAIR_DISPLAY_FIELDS}
             for r in data.get("pairs",[])],"hypotheses":selected}
    for key in ("anchors","matched_loci","probe_comparisons","same_platform_context",
                "context_matches","profile_signals","sensitivity"):
        display[key]=[r for r in data.get(key,[]) if r.get("HYPOTHESIS_ID") in ids]
    for key in ("same_cohort_screen","same_cohort_screen_matches","same_cohort_profile_signals"):
        display[key]=data.get(key,[])
    display["display_scope"]={"pair_summaries":"all directed primary/protein/alternative pairs",
        "detailed_hypotheses":len(selected),"detail_rule":"every prior priority, T1, or exact-anchor multi-region hypothesis",
        "complete_tables":"fingerprints/assay-evidence.tsv and fingerprints/pair-summary.tsv",
        "same_cohort_rows":len(display["same_cohort_screen"])}
    return display


def inline_assets(template,payload):
    body=json.dumps(payload,ensure_ascii=False,separators=(",",":"),allow_nan=False).replace("<","\\u003c")
    return (template.replace("__PAYLOAD__",body)
            .replace("__FINGERPRINT_CSS__",Path(__file__).with_name("fingerprint_view.css").read_text())
            .replace("__FINGERPRINT_SCRIPT__",Path(__file__).with_name("fingerprint_view.js").read_text()))


def render(cross_dir,specificity_dir,readiness_dir,out=None,fingerprints_dir=None):
    cross_dir,specificity_dir,readiness_dir=map(Path,(cross_dir,specificity_dir,readiness_dir))
    base=json.loads((cross_dir/"cross-check.data.json").read_text())
    spec=json.loads((specificity_dir/"assay-specificity.data.json").read_text())
    spec_summary=json.loads((specificity_dir/"assay-specificity.summary.json").read_text())
    ready=json.loads((readiness_dir/"ml-readiness.json").read_text())
    if spec.get("summary")!=spec_summary:raise ValueError("Specificity data and summary disagree")
    validate_files(specificity_dir,spec_summary.get("output_sha256",{}))
    validate_files(cross_dir,ready["input_hashes"])
    # The specificity receipt records its cross-check input hashes. Verify the
    # exact shared JSON, even if other receipt keys are added in a later version.
    receipt=spec_summary.get("input_provenance",{})
    def contains_hash(value,expected):
        if isinstance(value,dict):return any(contains_hash(v,expected) for v in value.values())
        if isinstance(value,list):return any(contains_hash(v,expected) for v in value)
        return value==expected
    if not contains_hash(receipt,sha(cross_dir/"cross-check.data.json")):
        raise ValueError("Specificity analysis does not match this cross-check JSON")
    display={k:spec.get(k,[]) for k in ("discordance","canonical")}
    alternatives=spec.get("alternatives",[])
    display["alternatives"]=[r for r in alternatives if r["PRIMARY"]=="ukb" or r["FOLLOW_UP_CLASS"]!="hotspot_or_ambiguous_anchor"]
    display["alternative_counts"]=dict(Counter(r["PRIMARY"]+"::"+r["PROT"] for r in alternatives))
    display["alternative_display_scope"]="All UKB hypotheses and other-cohort hypotheses outside hotspot_or_ambiguous_anchor; complete alternatives remain in alternative-targets.tsv"
    payload={**base,"specificity":display,"ml_readiness":ready}
    payload["specificity"]["summary"]=spec_summary
    if fingerprints_dir:
        payload["fingerprints"]=fingerprint_display(fingerprints_dir,cross_dir,specificity_dir)
        payload["fingerprints"]["axes"]=base.get("axes",[])
    template=Path(__file__).with_name("cross_check.html").read_text()
    out=Path(out) if out else cross_dir/"cross-check.html"
    out.write_text(inline_assets(template,payload),encoding="utf-8")
    receipt={"schema":"phole-review-render-v2","source_json_sha256":sha(cross_dir/"cross-check.data.json"),
        "specificity_json_sha256":sha(specificity_dir/"assay-specificity.data.json"),"readiness_sha256":sha(readiness_dir/"ml-readiness.json"),
        "template_sha256":sha(Path(__file__).with_name("cross_check.html")),"html_sha256":sha(out)}
    if fingerprints_dir:
        focused=out.parent/"fingerprint-casebook.html"
        focused_template='''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>PHOLE · 检测靶标遗传指纹</title>
<style>body{margin:0;background:#f2f6f8;color:#244052;font-family:system-ui,sans-serif}header{background:#16384a;color:white;padding:32px max(22px,calc((100vw - 1420px)/2))}header h1{font-size:28px;margin:6px 0}header p{color:#c3d7e2;font-size:13px}main{max-width:1420px;margin:24px auto;padding:24px;background:white;border:1px solid #d4e1e8;border-radius:9px}__FINGERPRINT_CSS__</style></head><body>
<header><span>PHOLE / CROSS-CHECK</span><h1>检测靶标遗传指纹 · 完整候选复核</h1><p>全部定向蛋白对摘要 · exact-anchor 多区域病例 · 同队列高相似性补充</p></header><main><div id="fingerprintReview"></div></main>
<script id="payload" type="application/json">__PAYLOAD__</script><script>__FINGERPRINT_SCRIPT__</script>
<script>window.fingerprintReview=PHOLEFingerprintReview.mount(document.getElementById('fingerprintReview'),JSON.parse(document.getElementById('payload').textContent));</script></body></html>'''
        focused.write_text(inline_assets(focused_template,payload["fingerprints"]),encoding="utf-8")
        receipt.update(fingerprint_json_sha256=sha(Path(fingerprints_dir)/"fingerprint-dossiers.data.json"),
            fingerprint_css_sha256=sha(Path(__file__).with_name("fingerprint_view.css")),
            fingerprint_js_sha256=sha(Path(__file__).with_name("fingerprint_view.js")),
            focused_html_sha256=sha(focused),fingerprint_display_scope=payload["fingerprints"]["display_scope"])
    (out.parent/"review-render.provenance.json").write_text(json.dumps(receipt,indent=2)+"\n")
    return receipt


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ("cross-dir","specificity-dir","readiness-dir"):ap.add_argument("--"+key,type=Path,required=True)
    ap.add_argument("--out",type=Path);ap.add_argument("--fingerprints-dir",type=Path);a=ap.parse_args()
    print(json.dumps(render(a.cross_dir,a.specificity_dir,a.readiness_dir,a.out,a.fingerprints_dir),indent=2))
if __name__=="__main__":main()
