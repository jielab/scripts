#!/usr/bin/env python3
"""Check exact-variant signed pQTL replication across independent cohorts.

Primary candidates and the external assay-test family are frozen before any
external P-value is compared with a threshold. Multiple assays never count as
multiple cohort confirmations. This is separate from latent-model discovery.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import tempfile
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

import numpy as np

DATASETS=("ukb","ckb","decode","scot")
PLATFORMS={"ukb":"Olink","ckb":"Olink","decode":"SomaScan","scot":"MS"}
SCHEMA="phole-edge-evidence-v2"


def digest_file(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for data in iter(lambda:stream.read(1024*1024),b""):h.update(data)
    return h.hexdigest()


def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as stream:
        return list(csv.DictReader(stream,delimiter="\t"))


def atomic_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fd,name=tempfile.mkstemp(prefix="."+path.name,dir=path.parent)
    try:
        with os.fdopen(fd,"w") as stream:
            json.dump(value,stream,ensure_ascii=False,indent=2,allow_nan=False);stream.write("\n")
        os.replace(name,path)
    finally:Path(name).unlink(missing_ok=True)


def atomic_table(path,rows,fields):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fd,name=tempfile.mkstemp(prefix="."+path.name,dir=path.parent)
    try:
        with os.fdopen(fd,"wb") as raw:
            output=gzip.GzipFile(filename="",fileobj=raw,mode="wb",mtime=0) if path.suffix==".gz" else raw
            with io.TextIOWrapper(output,encoding="utf-8",newline="") as stream:
                writer=csv.DictWriter(stream,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
                writer.writeheader();writer.writerows(rows)
        os.replace(name,path)
    finally:Path(name).unlink(missing_ok=True)


def source_identity(rows):
    required={"SOURCE_INDEX","LOCUS","CHR","POS","EA","NEA","BUILD"}
    if not rows or not required.issubset(rows[0]):raise ValueError("source.tsv lacks complete, allele-oriented anchor identity")
    identity=[]
    for i,row in enumerate(rows):
        if int(row["SOURCE_INDEX"])!=i:raise ValueError("SOURCE_INDEX does not match matrix row order")
        if not row["EA"] or not row["NEA"] or row["EA"]==row["NEA"]:raise ValueError("Invalid source EA/NEA")
        identity.append((row["LOCUS"],row["CHR"],int(row["POS"]),row["EA"],row["NEA"],row["BUILD"]))
    if len(set(identity))!=len(identity):raise ValueError("Duplicate source anchors")
    return identity


def load_dataset(root,dataset,reference=None):
    folder=Path(root)/dataset
    if not folder.exists():return None
    meta_path=folder/"manifest.json"
    if not meta_path.is_file():raise ValueError(f"Incomplete matrix directory: {folder}")
    meta=json.loads(meta_path.read_text())
    if (meta.get("schema")!="phole-signed-matrix-v2" or meta.get("directional") is not True or
            meta.get("input_scope")!="full" or meta.get("dataset")!=dataset):
        raise ValueError(f"{dataset}: signed complete-summary matrix receipt required; sig/MAGMA matrices are invalid")
    source=read_tsv(folder/"source.tsv");targets=read_tsv(folder/"target.tsv")
    source_hash=digest_file(folder/"source.tsv");target_hash=digest_file(folder/"target.tsv")
    if source_hash!=meta.get("source_hash"):raise ValueError(f"{dataset}: source.tsv checksum mismatch")
    if target_hash!=meta.get("target_hash"):raise ValueError(f"{dataset}: target.tsv checksum mismatch")
    identity=source_identity(source)
    if reference is not None and (source_hash!=reference["source_hash"] or identity!=reference["identity"]):
        raise ValueError(f"{dataset}: source order, coordinate, EA/NEA or provenance differs from the primary panel")
    shape=(len(source),len(targets))
    if list(shape)!=meta.get("shape"):raise ValueError(f"{dataset}: matrix shape receipt mismatch")
    arrays={};checksums={}
    for name in ("beta","se","p","observed","cis"):
        path=folder/(name+".npy")
        actual=digest_file(path)
        if actual!=meta.get("array_checksums",{}).get(name):raise ValueError(f"{dataset}: {name}.npy checksum mismatch")
        array=np.load(path,mmap_mode="r",allow_pickle=False)
        if array.shape!=shape:raise ValueError(f"{dataset}: {name} shape mismatch")
        if name in {"observed","cis"} and array.dtype!=np.dtype(bool):raise ValueError(f"{dataset}: {name} must be boolean")
        arrays[name]=array;checksums[name]=actual
    known=np.load(folder/"target_cis_known.npy",allow_pickle=False)
    if known.shape!=(len(targets),) or known.dtype!=np.dtype(bool):raise ValueError(f"{dataset}: invalid target_cis_known")
    groups=defaultdict(list);seen=set()
    for j,t in enumerate(targets):
        if int(t["TARGET_INDEX"])!=j or not t.get("PROT") or t["PROT"]==".":raise ValueError(f"{dataset}: invalid target identity/order")
        if t["ASSAY"] in seen:raise ValueError(f"{dataset}: duplicate assay identifier")
        if t.get("CIS_KNOWN") not in {"0","1"} or bool(int(t["CIS_KNOWN"]))!=bool(known[j]):
            raise ValueError(f"{dataset}: cis-annotation receipt mismatch")
        seen.add(t["ASSAY"]);groups[t["PROT"]].append(j)
    valid=(arrays["observed"] & np.isfinite(arrays["beta"]) & np.isfinite(arrays["se"]) & (arrays["se"]>0) &
           np.isfinite(arrays["p"]) & (arrays["p"]>=0) & (arrays["p"]<=1))
    return {"dataset":dataset,"folder":folder,"meta":meta,"source":source,"source_hash":source_hash,"identity":identity,
            "target":targets,"groups":dict(groups),"known":known,"valid":valid,**arrays,
            "provenance":{"manifest_sha256":digest_file(meta_path),"source_sha256":source_hash,"target_sha256":target_hash,
                          "array_sha256":checksums,"target_cis_known_sha256":digest_file(folder/"target_cis_known.npy")},
            "invalid_observed_cells":int(np.sum(arrays["observed"] & ~valid))}


def primary_candidates(primary,threshold):
    candidates=[];screen=[]
    for prot,columns in sorted(primary["groups"].items()):
        hits=np.flatnonzero(np.any(primary["valid"][:,columns] & (primary["p"][:,columns]<=threshold),axis=1))
        for i0 in hits:
            i=int(i0);observed=[j for j in columns if primary["valid"][i,j]]
            signs={int(np.sign(primary["beta"][i,j])) for j in observed}
            eligible=len(signs)==1 and 0 not in signs
            reason="eligible" if eligible else "primary_assays_disagree_or_include_zero"
            results=[{"assay":primary["target"][j]["ASSAY"],"observed":bool(primary["valid"][i,j]),
                      "beta":float(primary["beta"][i,j]) if primary["valid"][i,j] else None,
                      "p":float(primary["p"][i,j]) if primary["valid"][i,j] else None} for j in columns]
            row={"SOURCE_INDEX":i,"LOCUS":primary["source"][i]["LOCUS"],"PROT":prot,"PRIMARY_ELIGIBLE":int(eligible),
                 "N_EXPECTED":len(columns),"N_OBSERVED":len(observed),"PRIMARY_THRESHOLD":threshold,
                 "MIN_P":min(float(primary["p"][i,j]) for j in observed),"STATUS":reason,
                 "ASSAY_RESULTS":json.dumps(results,separators=(",",":"),allow_nan=False)}
            screen.append(row)
            if eligible:candidates.append({"i":i,"prot":prot,"sign":next(iter(signs)),"primary":row})
    candidates.sort(key=lambda r:(r["i"],r["prot"]))
    screen.sort(key=lambda r:(r["SOURCE_INDEX"],r["PROT"]))
    return candidates,screen


def edge_class(candidate,datasets):
    states=[];unknown=[]
    for name,data in datasets.items():
        if data is None:continue
        for j in data["groups"].get(candidate["prot"],[]):
            if not data["known"][j]:unknown.append(name+":"+data["target"][j]["ASSAY"])
            else:states.append(bool(data["cis"][candidate["i"],j]))
    if unknown:return "unknown",["unknown_target_annotation:"+";".join(unknown)]
    if not states or len(set(states))!=1:return "unknown",["target_cis_classification_conflict"]
    return ("cis" if states[0] else "trans"),[]


def cohort_evidence(candidate,data,threshold,is_primary=False):
    if data is None:
        return {"N_EXPECTED":0,"N_OBSERVED":0,"N_SAME_SIGN":0,"N_OPPOSITE_SIGN":0,"N_OPPOSITE_SIGNIFICANT":0,
                "MIN_P":".","MAX_P":".","MISSING_ASSAYS":".","ASSAY_RESULTS":"[]","STATUS":"cohort_unavailable","PASS":0,"FLAGS":"cohort_unavailable"}
    columns=data["groups"].get(candidate["prot"],[]);i=candidate["i"]
    observed=[j for j in columns if data["valid"][i,j]];missing=[data["target"][j]["ASSAY"] for j in columns if not data["valid"][i,j]]
    signs=[int(np.sign(data["beta"][i,j])) for j in observed]
    same=sum(s==candidate["sign"] for s in signs)
    opposite=sum(s==-candidate["sign"] for s in signs)
    opposite_sig=int(sum(int(np.sign(data["beta"][i,j]))==-candidate["sign"] and threshold is not None and data["p"][i,j]<=threshold for j in observed))
    passed=bool(columns) and len(observed)==len(columns) and same==len(columns) and threshold is not None and all(data["p"][i,j]<=threshold for j in observed)
    if is_primary:passed=True
    flags=[]
    if missing:flags.append("incomplete_assay_coverage")
    if len(set(signs))>1:flags.append("assay_direction_heterogeneity")
    if opposite_sig:flags.append("opposite_significant_assay")
    if any(data["p"][i,j]==0 for j in observed):flags.append("p_underflow_zero")
    status=("primary_strict" if is_primary else "not_measured" if not columns else "discordant" if opposite_sig else
            "unavailable" if not observed else "incomplete" if missing else "replicated" if passed else "not_confirmed")
    results=[{"assay":data["target"][j]["ASSAY"],"observed":bool(data["valid"][i,j]),
              "beta":float(data["beta"][i,j]) if data["valid"][i,j] else None,
              "se":float(data["se"][i,j]) if data["valid"][i,j] else None,
              "p":float(data["p"][i,j]) if data["valid"][i,j] else None,
              "same_primary_sign":bool(int(np.sign(data["beta"][i,j]))==candidate["sign"]) if data["valid"][i,j] else None} for j in columns]
    return {"N_EXPECTED":len(columns),"N_OBSERVED":len(observed),"N_SAME_SIGN":same,"N_OPPOSITE_SIGN":opposite,
            "N_OPPOSITE_SIGNIFICANT":opposite_sig,
            "MIN_P":min(float(data["p"][i,j]) for j in observed) if observed else ".",
            "MAX_P":max(float(data["p"][i,j]) for j in observed) if observed else ".",
            "MISSING_ASSAYS":";".join(missing) or ".","ASSAY_RESULTS":json.dumps(results,separators=(",",":"),allow_nan=False),
            "STATUS":status,"PASS":int(passed),"FLAGS":";".join(flags) or "."}


BASE_FIELDS=["EDGE","SOURCE_INDEX","LOCUS","SNP","CHR","POS","EA","NEA","PROT","CLASS","PRIMARY_DATASET","PRIMARY_SIGN",
             "PRIMARY_THRESHOLD","PRIMARY_MIN_P","REPLICATION_THRESHOLD","N_EXTERNAL_ASSAY_TESTS","N_REPLICATED_COHORTS",
             "REPLICATED_COHORTS","ORTHOGONAL_REPLICATED_COHORTS","N_OPPOSITE_SIGNIFICANT","STABLE","FLAGS"]
COHORT_FIELDS=["N_EXPECTED","N_OBSERVED","N_SAME_SIGN","N_OPPOSITE_SIGN","N_OPPOSITE_SIGNIFICANT","MIN_P","MAX_P",
               "MISSING_ASSAYS","ASSAY_RESULTS","STATUS","PASS","FLAGS"]
SCREEN_FIELDS=["SOURCE_INDEX","LOCUS","PROT","PRIMARY_ELIGIBLE","N_EXPECTED","N_OBSERVED","PRIMARY_THRESHOLD","MIN_P","STATUS","ASSAY_RESULTS"]


def run(matrix_root,out,primary="ukb",external=None,alpha=0.05,primary_p=5e-8):
    matrix_root=Path(matrix_root).resolve();out=Path(out).resolve()
    external=list(external if external is not None else [d for d in DATASETS if d!=primary])
    if primary not in DATASETS or len(set(external))!=len(external) or any(d not in DATASETS or d==primary for d in external):
        raise ValueError("External datasets must be distinct known cohorts, excluding primary")
    if not 0<alpha<1 or not 0<primary_p<1:raise ValueError("Thresholds must be in (0,1)")
    base=load_dataset(matrix_root,primary)
    if base is None:raise ValueError("Primary matrix is unavailable: "+primary)
    if not base["target"]:raise ValueError("Primary matrix has no assays")
    loaded={primary:base}
    for d in external:loaded[d]=load_dataset(matrix_root,d,base)
    if len({bool(d["meta"].get("synthetic",False)) for d in loaded.values() if d is not None})!=1:
        raise ValueError("Synthetic and real matrices must not be combined")
    primary_threshold=primary_p/len(base["target"])
    candidates,screen=primary_candidates(base,primary_threshold)
    # Freeze the external assay-test family using observation availability only.
    # No external P-value cutoff, sign or assay result changes its denominator.
    counts={d:sum(int(np.sum(loaded[d]["valid"][c["i"],loaded[d]["groups"].get(c["prot"],[])])) for c in candidates)
            if loaded[d] is not None else 0 for d in external}
    total_tests=sum(counts.values());threshold=alpha/total_tests if total_tests else None
    cohort_status={d:Counter() for d in external};edges=[];stable=[];class_counts=Counter()
    for number,c in enumerate(candidates,1):
        source=base["source"][c["i"]];kind,flags=edge_class(c,loaded);class_counts[kind]+=1
        evidence={d:cohort_evidence(c,data,primary_threshold if d==primary else threshold,d==primary) for d,data in loaded.items()}
        replicated=[d for d in external if evidence[d]["PASS"]]
        orthogonal=[d for d in replicated if PLATFORMS[d]!=PLATFORMS[primary]]
        opposite=sum(evidence[d]["N_OPPOSITE_SIGNIFICANT"] for d in external)
        is_stable=len(replicated)>=2 and bool(orthogonal) and opposite==0 and kind!="unknown"
        if evidence[primary]["N_OBSERVED"]<evidence[primary]["N_EXPECTED"]:flags.append("primary_assay_missing")
        if opposite:flags.append("significant_opposite_external_assay")
        if len(replicated)>=2 and not orthogonal:flags.append("no_orthogonal_platform_replication")
        row={"EDGE":f"EDGE{number:07d}","SOURCE_INDEX":c["i"],"LOCUS":source["LOCUS"],"SNP":source.get("SNP","."),
             "CHR":source["CHR"],"POS":source["POS"],"EA":source["EA"],"NEA":source["NEA"],"PROT":c["prot"],"CLASS":kind,
             "PRIMARY_DATASET":primary,"PRIMARY_SIGN":c["sign"],"PRIMARY_THRESHOLD":primary_threshold,"PRIMARY_MIN_P":c["primary"]["MIN_P"],
             "REPLICATION_THRESHOLD":threshold if threshold is not None else ".","N_EXTERNAL_ASSAY_TESTS":sum(evidence[d]["N_OBSERVED"] for d in external),
             "N_REPLICATED_COHORTS":len(replicated),"REPLICATED_COHORTS":";".join(replicated) or ".",
             "ORTHOGONAL_REPLICATED_COHORTS":";".join(orthogonal) or ".","N_OPPOSITE_SIGNIFICANT":opposite,
             "STABLE":int(is_stable),"FLAGS":";".join(flags) or "."}
        for d in loaded:
            row.update({d.upper()+"_"+k:v for k,v in evidence[d].items()})
            if d!=primary:cohort_status[d][evidence[d]["STATUS"]]+=1
        edges.append(row)
        if is_stable:stable.append(row)
    fields=BASE_FIELDS+[d.upper()+"_"+field for d in loaded for field in COHORT_FIELDS]
    out.mkdir(parents=True,exist_ok=True)
    atomic_table(out/"primary-screen.tsv.gz",screen,SCREEN_FIELDS)
    atomic_table(out/"edge-replication.tsv.gz",edges,fields)
    atomic_table(out/"stable-edges.tsv",stable,fields)
    atomic_table(out/"stable-cis.tsv",[r for r in stable if r["CLASS"]=="cis"],fields)
    atomic_table(out/"stable-trans.tsv",[r for r in stable if r["CLASS"]=="trans"],fields)
    summary={"schema":SCHEMA,"generated":datetime.now(timezone.utc).isoformat(timespec="seconds"),
             "primary":primary,"external_datasets":external,"platforms":PLATFORMS,"matrix_root":str(matrix_root),
             "synthetic":bool(base["meta"].get("synthetic",False)),"n_primary_assays":len(base["target"]),
             "n_primary_significant_edges_screened":len(screen),"n_primary_sign_heterogeneous_excluded":sum(not r["PRIMARY_ELIGIBLE"] for r in screen),
             "n_primary_candidates":len(candidates),"primary_threshold":primary_threshold,
             "primary_selection":f"all observed assays have one identical nonzero sign; at least one P <= {primary_p:g} / all primary assays; missing primary assays retained and flagged",
             "external_family_alpha":alpha,"external_assay_tests":total_tests,"external_assay_tests_by_cohort":counts,
             "replication_threshold":threshold,"external_test_family":"all observed valid external assay cells for the frozen primary candidates, counted before comparing any external P-value",
             "cohort_replication":"all expected assays must be observed, have the primary sign and pass the single Bonferroni threshold; one cohort contributes at most one confirmation",
             "n_stable":len(stable),"n_stable_cis":sum(r["CLASS"]=="cis" for r in stable),"n_stable_trans":sum(r["CLASS"]=="trans" for r in stable),
             "n_edges_with_opposite_significant_assay":sum(r["N_OPPOSITE_SIGNIFICANT"]>0 for r in edges),
             "candidate_class_counts":dict(class_counts),"external_cohort_status_counts":{d:dict(c) for d,c in cohort_status.items()},
             "stable_rule":"primary eligibility + >=2 external cohorts + >=1 platform different from primary + no significant opposite assay + known consistent cis/trans annotation",
             "missing_semantics":"unavailable cohort, unmeasured protein, unobserved assay and incomplete assay coverage are distinct; absence is never zero or a negative association",
             "evidence_scope":"exact same allele-oriented variant and effect direction; no LD equivalence, colocalization, mediation or causal-gene claim",
             "selection_separation":"external stability lists are independent evidence/positive controls/sensitivity strata; never used to select targets or candidates for formal held-out latent discovery",
             "panel_completeness":{d:data["meta"].get("panel_completeness","provided_unverified") for d,data in loaded.items() if data is not None},
             "panel_limit":"Results are conditional on the frozen source panel; even a complete candidate export is reduced to chosen physical-bin anchors",
             "invalid_observed_cells":{d:data["invalid_observed_cells"] for d,data in loaded.items() if data is not None},
             "source_sha256":base["source_hash"],"source_identity_fields":["LOCUS","CHR","POS","EA","NEA","BUILD"],
             "provenance":{d:data["provenance"] for d,data in loaded.items() if data is not None},
             "implementation_sha256":digest_file(__file__),
             "output_sha256":{name:digest_file(out/name) for name in ["primary-screen.tsv.gz","edge-replication.tsv.gz","stable-edges.tsv","stable-cis.tsv","stable-trans.tsv"]}}
    atomic_json(out/"edge-evidence.summary.json",summary)
    print(json.dumps({k:summary[k] for k in ["schema","synthetic","n_primary_candidates","external_assay_tests","replication_threshold","n_stable","n_stable_cis","n_stable_trans"]},indent=2))
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix-root",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    p.add_argument("--primary",choices=DATASETS,default="ukb");p.add_argument("--datasets",help="External cohorts, comma separated [all others]")
    p.add_argument("--alpha",type=float,default=0.05);p.add_argument("--primary-p",type=float,default=5e-8)
    a=p.parse_args();external=None if a.datasets is None else [d.strip() for d in a.datasets.split(",") if d.strip()]
    run(a.matrix_root,a.out,a.primary,external,a.alpha,a.primary_p)


if __name__=="__main__":
    try:main()
    except (ValueError,OSError,KeyError) as error:raise SystemExit("ERROR: "+str(error)) from error
