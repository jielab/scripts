#!/usr/bin/env python3
"""Count evidence and effective study units without pretending summaries form a matrix.

This is a descriptive feasibility audit. Physical bins/blocks are not independent
LD loci, multi-cohort evidence is not an untouched validation set, and selected
significant leads cannot supply the nonsignificant cells needed by PHOLE.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import itertools
import json
from collections import Counter, defaultdict
from decimal import Decimal
from pathlib import Path

DATASETS=("ukb","ckb","decode","scot")
BAD={"identical_alleles_in_id","invalid_effect_or_se","reported_nonsignificant"}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def write_tsv(path,rows,fields):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore");w.writeheader();w.writerows(rows)
def save_json(path,obj):Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
def valid_target(p):
    t=p.get("target")
    return bool(t) and t.get("status")=="single_gene" and t.get("identity_verified",True) not in (False,0,"0","false")
def qpass(s):return not (set(s.get("flags",[]))&BAD)
def assay_metadata_usable(v):
    return all(a.get("status")=="complete" and not any(str(f).startswith("metadata_") or f in
        {"file_dataset_build_conflict","target_or_build_unresolved","cis_effect_invalid"} for f in a.get("flags",[]))
        for a in v.get("assays",[]))
def strict(s,d,thresholds):return Decimal(str(s["p"]))<=Decimal(str(thresholds[d])) and qpass(s)
def bin_for(ch,pos,size):return str(ch),(int(pos)-1)//size
def block_for(ch,pos):
    # Keep the extended MHC together; retain X/Y/MT as real chromosome classes.
    if str(ch)=="6" and 25000000<=int(pos)<=35000000:return "6:MHC25-35Mb"
    return str(ch)+":"+str((int(pos)-1)//5000000)


def audit(data,regions,primary="ukb"):
    proteins=data["proteins"];thresholds=data["config"]["strong_thresholds"]
    measured={d:{p["prot"] for p in proteins if p["data"][d].get("aliases")} for d in DATASETS}
    valid={p["prot"] for p in proteins if valid_target(p)}
    priority={p["prot"] for p in proteins if p.get("priority") and valid_target(p)}
    own={p["prot"] for p in proteins if valid_target(p) and p["data"][primary].get("status")=="complete"
         and any(s.get("kind")=="cis" and strict(s,primary,thresholds) for s in p["data"][primary].get("signals",[]))
         and not p["data"][primary].get("assay_heterogeneity")
         and not p["data"][primary].get("assay_support_partial") and assay_metadata_usable(p["data"][primary])}
    cohorts=[];all_primary=[]
    for d in DATASETS:
        edges=set();strict_edges=set();positions=set();blocks=set();targets=set();cis=set()
        for p in proteins:
            for s in p["data"][d].get("signals",[]):
                if s.get("kind") not in {"cis","trans"}:continue
                key=(p["prot"],str(s["chr"]),int(s["pos"]),str(s["snp"]),s["kind"])
                edges.add(key)
                if not strict(s,d,thresholds) or not assay_metadata_usable(p["data"][d]):continue
                strict_edges.add(key);positions.add((str(s["chr"]),int(s["pos"])));blocks.add(block_for(s["chr"],s["pos"]))
                if s["kind"]=="cis":cis.add(p["prot"])
                else:
                    targets.add(p["prot"])
                    if d==primary:all_primary.append((p["prot"],s))
        cohorts.append({"DATASET":d,"MAPPED_PROTEINS":len(measured[d]),"READABLE_PROTEINS":sum(p["data"][d].get("status") in {"complete","partial"} for p in proteins),
            "REPORTED_EDGES":len(edges),"STRICT_QC_EDGES":len(strict_edges),"STRICT_QC_CIS_PROTEINS":len(cis),
            "STRICT_QC_TRANS_PROTEINS":len(targets),"STRICT_QC_VARIANT_POSITIONS":len(positions),"PHYSICAL_BLOCKS_5MB":len(blocks),
            "COMPLETE_EXPORTS":data["summary"]["datasets"][d].get("complete_full_summaries",False)})
    overlap=[{"DATASETS":";".join(ds),"MAPPED_PROTEINS":len(set.intersection(*(measured[d] for d in ds))),
              "PRIORITY_CIS_PROTEINS":len(priority&set.intersection(*(measured[d] for d in ds)))}
             for n in (2,3,4) for ds in itertools.combinations(DATASETS,n)]
    region_rows=[r for r in regions if str(r.get("REGION_PRIORITY")) in {"1","True","true"}]
    # A protein-region record can have several lead positions. Assign its centre
    # to a fixed physical bin only for a transparent, reproducible size audit.
    bybin=defaultdict(set);byblock=defaultdict(set);bytarget=defaultdict(set);region_bins=set()
    for r in region_rows:
        middle=(int(r["START"])+int(r["END"]))//2;key=bin_for(r["CHR"],middle,1000000)
        bybin[key].add(r["prot"]);byblock[block_for(r["CHR"],middle)].add(r["prot"])
        bytarget[r["prot"]].add(key);region_bins.add(key)
    hot=sorted([{"CHR":k[0],"START":k[1]*1000000+1,"END":(k[1]+1)*1000000,"TARGETS":len(v),"PROTEINS":";".join(sorted(v))}
                for k,v in bybin.items()],key=lambda r:(-r["TARGETS"],str(r["CHR"]),r["START"]))
    top5keys={(r["CHR"],(r["START"]-1)//1000000) for r in hot[:5]}
    top5edges=sum(1 for r in region_rows if bin_for(r["CHR"],(int(r["START"])+int(r["END"]))//2,1000000) in top5keys)
    cohort_qc={r["DATASET"]:r for r in cohorts}
    panels=[]
    for name,ps,scope in [("primary_strict_cis",own,"primary_only"),("cross_cohort_cis_priority",priority,"external_screened")]:
        trans={(pr,str(s["chr"]),int(s["pos"]),s["snp"]) for pr,s in all_primary if pr in ps}
        pos={(r[1],r[2]) for r in trans};bins={bin_for(c,b,1000000) for c,b in pos};bs={block_for(c,b) for c,b in pos}
        ct=Counter(r[0] for r in trans)
        panels.append({"PANEL":name,"SELECTION_SCOPE":scope,"PROTEINS":len(ps),"PRIMARY_MAPPED":len(ps&measured[primary]),
            "PRIMARY_STRICT_TRANS_EDGES":len(trans),"PRIMARY_TRANS_POSITIONS":len(pos),"PRIMARY_TRANS_BINS_1MB":len(bins),
            "PRIMARY_TRANS_BLOCKS_5MB":len(bs),"TARGETS_WITH_TRANS":len(ct),"TARGETS_WITH_3PLUS_TRANS":sum(v>=3 for v in ct.values()),
            "MAPPED_IN_DECODE":len(ps&measured["decode"]),"MAPPED_IN_CKB":len(ps&measured["ckb"]),"MAPPED_IN_SCOT":len(ps&measured["scot"]),
            "DENSE_PANEL_CELL_ILLUSTRATION":len(bins)*len(ps&measured[primary]),
            "CELL_ILLUSTRATION_SCOPE":"1 anchor per occupied physical trans bin x mapped targets; hypothetical, not observed cells"})
    genes={p["target"].get("gene_id") or p["target"].get("gene") for p in proteins if p["prot"] in priority}
    result={"schema":"phole-ml-readiness-v1","primary":primary,"summary_only":True,"signed_matrix_available":False,
        "conclusion":"full-GWAS extraction is required; these counts describe candidate pilot scope, not a model-power or latent-node discovery guarantee",
        "priority_cis_proteins":len(priority),"priority_cis_unique_coding_ids":len(genes),"priority_cis_with_ms":sum(p["prot"] in priority and p.get("ms_support",False) for p in proteins),
        "primary_strict_cis_screen_proteins":len(own),"cross_platform_strict_trans_regions":len(region_rows),
        "cross_platform_trans_targets":len(bytarget),"cross_platform_trans_bins_1mb":len(region_bins),"cross_platform_trans_blocks_5mb":len(byblock),
        "cross_platform_trans_targets_3plus_regions":sum(len(v)>=3 for v in bytarget.values()),
        "cross_platform_trans_top5_bin_edges":top5edges,"cross_platform_trans_top5_bin_fraction":top5edges/len(region_rows) if region_rows else None,
        "cis_and_trans_priority_intersection":len(priority&set(bytarget)),
        "signed_replicated_edges":None,"independent_ld_loci":None,"observed_matrix_cells":None,"observed_matrix_density":None,
        "cohorts":cohorts,"overlap":overlap,"panels":panels,"top_trans_bins":hot[:20],
        "limitations":["Selected lead summaries do not contain nonsignificant effects, variant coverage or a complete signed matrix.",
            "Physical bins and 5 Mb blocks count positional diversity, not LD-independent mechanisms.",
            "The cross-cohort priority panel used external cohorts for selection; those cohorts are not untouched validation for that panel.",
            "The primary-cis panel is a conservative calibration set. Trans-only targets may be biological and are not declared wrong proteins.",
            "No universal edge-count cutoff establishes machine-learning power; evaluate held-out losses, block uncertainty, coverage and null controls.",
            "A low-rank latent factor is not uniquely identified with one causal gene or protein."]}
    return result,own,priority,hot


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--cross-dir",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True);ap.add_argument("--primary",choices=DATASETS,default="ukb")
    args=ap.parse_args(argv);folder=args.cross_dir
    data=json.loads((folder/"cross-check.data.json").read_text());summary=json.loads((folder/"summary.json").read_text())
    provenance=json.loads((folder/"provenance.json").read_text())
    for name,expected in provenance.get("output_sha256",{}).items():
        candidate=folder/name
        if Path(name).name!=name:raise ValueError("Unsafe filename in cross-check provenance")
        if not candidate.is_file() or sha(candidate)!=expected:raise ValueError("Cross-check output checksum mismatch: "+name)
    if data.get("summary")!=summary:raise ValueError("Cross-check JSON and summary disagree; rerun cross-check")
    if len({p["prot"] for p in data["proteins"]})!=len(data["proteins"]):raise ValueError("Duplicate canonical proteins")
    if summary["counts"]["proteins"]!=len(data["proteins"]):raise ValueError("Cross-check row count mismatch")
    regions=read_tsv(folder/"trans-region-priority.tsv")
    if len(regions)!=summary["counts"]["trans_regions_strict3plus"]:raise ValueError("Trans region table and summary count disagree")
    result,own,priority,hot=audit(data,regions,args.primary)
    hashes={name:sha(folder/name) for name in ("cross-check.data.json","summary.json","trans-region-priority.tsv","provenance.json")}
    result["input_hashes"]=hashes
    args.out.mkdir(parents=True,exist_ok=True)
    save_json(args.out/"ml-readiness.json",result)
    for name,rows in (("cohort-size.tsv",result["cohorts"]),("cohort-overlap.tsv",result["overlap"]),("ml-panels.tsv",result["panels"]),("priority-trans-hotspots.tsv",hot)):
        fields=list(rows[0]) if rows else ["CHR","START","END","TARGETS","PROTEINS"]
        write_tsv(args.out/name,rows,fields)
    for name,ps in (("ml-targets.primary-cis.txt",own),("ml-targets.cross-cohort-priority.txt",priority)):
        (args.out/name).write_text("\n".join(sorted(ps))+"\n")
    fields=["METRIC","VALUE","SCOPE"]
    rows=[{"METRIC":k,"VALUE":"not_estimated" if v is None else v,"SCOPE":"reported-summary feasibility audit"}
          for k,v in result.items() if isinstance(v,(int,float,bool)) or v is None]
    write_tsv(args.out/"ml-readiness.tsv",rows,fields)
    print(json.dumps({k:v for k,v in result.items() if not isinstance(v,(dict,list))},ensure_ascii=False,indent=2))


if __name__=="__main__":main()
