#!/usr/bin/env python3
"""Build signed locus-by-assay matrices from complete allele-aware GWAS files.

The discovery panel is fixed before other cohorts are read. Non-significance is
retained; absent, invalid and unresolved variants remain masked. MAGMA gene
statistics and mplot's selected .sig.txt rows are never signed matrix inputs.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from summary_metadata import read_metadata, metadata_path

import numpy as np

from signed_input import (
    DATASETS, SCHEMA, aligned_at_anchor, atomic_json, atomic_tsv, chrom_order,
    digest, discovery_record, extract_positions, file_hash, file_identity,
    normal_build, normal_chr, palindromic, read_tsv, safe_assay, summary_identity, validate_manifest,
)

IMPLEMENTATION={"prepare":file_hash(Path(__file__)),"signed_input":file_hash(Path(__file__).with_name("signed_input.py"))}


def log(message):
    print(message, file=sys.stderr, flush=True)


def clean_json(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {k:clean_json(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)):
        return [clean_json(v) for v in value]
    return value


def arguments():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--assay-manifest", type=Path)
    p.add_argument("--candidate-loci", type=Path)
    p.add_argument("--protein-annotation", type=Path)
    p.add_argument("--cross-check-dir", type=Path)
    p.add_argument("--panel-completeness", choices=["full_summary_candidates","provided_unverified"],
                   default="provided_unverified",help="Scope of candidate-locus enumeration, distinct from complete per-assay summary inputs")
    p.add_argument("--panel-receipt",type=Path,help="Cross-check summary certifying candidate-enumeration scope")
    p.add_argument("--out", type=Path)
    p.add_argument("--state-dir", type=Path, help="Reusable preparation checkpoints")
    p.add_argument("--discovery", choices=DATASETS, default="ukb")
    p.add_argument("--discovery-protein-list", type=Path,
                   help="Optional prespecified training-only PROT list; only these proteins select panel anchors")
    p.add_argument("--target-protein-list", type=Path,
                   help="Optional explicit canonical PROT set retained in every cohort; selection scope is recorded")
    p.add_argument("--target-selection-scope", choices=["prespecified","primary_only","external_screened","unknown"], default="unknown",
                   help="Evidence used to choose --target-protein-list; external_screened cannot supply untouched cohort validation")
    p.add_argument("--p-threshold", type=float, default=5e-8)
    p.add_argument("--cis-flank", type=int, default=500000)
    p.add_argument("--bin-size", type=int, default=1000000)
    p.add_argument("--block-size", type=int, default=5000000)
    p.add_argument("--locus-bed", type=Path)
    p.add_argument("--locus-build", choices=["37","38"])
    p.add_argument("--locus-bed-coordinates", choices=["bed0","one-based"],default="bed0")
    p.add_argument("--palindrome-maf",type=float,default=0.42)
    p.add_argument("--palindrome-tolerance",type=float,default=0.10)
    p.add_argument("--allow-scan",action="store_true")
    p.add_argument("--tabix",default="tabix")
    p.add_argument("--replace",default="TRUE",help="TRUE discards per-assay extraction checkpoints")
    p.add_argument("--gwas-root",type=Path,help="Used only for exact-path manifest discovery")
    p.add_argument("--category",default="common",help="Current GWAS layout: ROOT/dataset/CATEGORY/assay/gwas/assay.gz")
    p.add_argument("--protein-map",type=Path,help="prot.names.tsv for --export-manifest")
    p.add_argument("--export-manifest",type=Path,help="Discover unique full-summary paths, write manifest and stop")
    p.add_argument("--dataset-builds",default="",help="Explicit receipts, e.g. ukb=38,ckb=38,decode=38,scot=38")
    p.add_argument("--dataset-ancestries",default="",help="E.g. ukb=EUR,ckb=EAS,decode=EUR,scot=EUR; unknown is kept unknown")
    p.add_argument("--build-receipt",type=Path,help="TSV with DATASET and BUILD, optionally ASSAY/PATH, for auto-manifest discovery")
    return p.parse_args()


def key_values(value):
    result={}
    for item in str(value).split(","):
        if item.strip():
            key,val=item.split("=",1)
            result[key.strip().lower()]=val.strip()
    return result


def local_build_receipts(path, folder, category, assay):
    """Read fresh receipts tied to the selected full summary, never its raw input.

    Upstream .grch caches can predate a liftOver output. Their freshness check is
    therefore material, as is the association_input path in a plot sidecar.
    """
    receipts=[];warnings=[]
    summary_mtime=path.stat().st_mtime_ns
    trait=folder/category/assay
    caches=set()
    if path.parent==(trait/"gwas").resolve():caches.add(trait/"qc"/(assay+".grch"))
    if path.parent.name in {"gwas","clean"}:
        caches.add(path.parent.parent/"qc"/(assay+".grch"))
    for receipt in sorted(caches):
        if not receipt.is_file():continue
        if receipt.stat().st_mtime_ns<summary_mtime:
            warnings.append(("stale_build_receipt",str(receipt)));continue
        text=receipt.read_text().strip().split()
        if len(text)!=1:
            warnings.append(("invalid_build_receipt",str(receipt)));continue
        try:build=normal_build(text[0])
        except ValueError:
            warnings.append(("invalid_build_receipt",str(receipt)));continue
        receipts.append((build,str(receipt.resolve())))
    meta_path=metadata_path(folder)
    if meta_path.is_file():
        meta=read_metadata(folder).get(assay,{})
        raw_path=meta.get("association_input","")
        # Exact absolute path identity only: never infer correspondence from an
        # assay basename, nor guess a Windows-to-WSL volume remapping.
        same_input=bool(raw_path) and Path(raw_path).is_absolute() and Path(raw_path).resolve()==path
        if not same_input or meta.get("association_scope","")!="full_gwas":
            warnings.append(("unmatched_build_metadata",str(meta_path)))
        elif meta_path.stat().st_mtime_ns<summary_mtime:
            warnings.append(("stale_build_receipt",str(meta_path)))
        else:
            try:build=normal_build(meta.get("build",""))
            except ValueError:warnings.append(("invalid_build_receipt",str(meta_path)))
            else:receipts.append((build,str(meta_path.resolve())))
    return receipts,warnings


def export_manifest(args):
    if args.gwas_root is None or args.protein_map is None:
        raise ValueError("--export-manifest requires --gwas-root and --protein-map")
    category=safe_assay(args.category)
    builds={};assay_builds={}
    for dataset,value in key_values(args.dataset_builds).items():
        d="ukb" if dataset=="prot" else dataset
        if d not in DATASETS:raise ValueError(f"Unknown dataset build receipt: {dataset}")
        builds[d]=(normal_build(value),"explicit --dataset-builds")
    if args.build_receipt:
        for raw in read_tsv(args.build_receipt):
            r={k.upper():v for k,v in raw.items()}
            d=r.get("DATASET","").lower()
            if d=="prot":d="ukb"
            if d in DATASETS:
                value=r.get("BUILD") or r.get("INFERRED_BUILD") or r.get("BEST_BUILD")
                if value:
                    assay=r.get("ASSAY","").strip()
                    key=(d,safe_assay(assay)) if assay and assay!="." else d
                    bank=assay_builds if isinstance(key,tuple) else builds
                    build=normal_build(value)
                    if key in bank and bank[key][0]!=build:
                        raise ValueError(f"Conflicting build receipts for {d}")
                    evidence=str(args.build_receipt.resolve())+(f"#ASSAY={assay}" if isinstance(key,tuple) else "#DATASET="+d)
                    receipt_path=r.get("PATH","")
                    if receipt_path and receipt_path!=".":
                        if not isinstance(key,tuple):raise ValueError("A PATH-specific build receipt also needs ASSAY")
                        receipt_file=Path(receipt_path).expanduser()
                        if not receipt_file.is_absolute():receipt_file=args.build_receipt.parent/receipt_file
                        receipt_path=str(receipt_file.resolve())
                    bank[key]=(build,evidence,receipt_path) if isinstance(key,tuple) else (build,evidence)
    ancestry=key_values(args.dataset_ancestries)
    mapping=read_tsv(args.protein_map)
    if not mapping or not {"prot",*DATASETS}.issubset(mapping[0]):
        raise ValueError("--protein-map must contain prot/ukb/ckb/decode/scot")
    known={}
    for row in mapping:
        for d in DATASETS:
            for alias in row[d].split(";"):
                if alias.strip() not in {"","."}:
                    key=(d,safe_assay(alias.strip()))
                    if key in known and known[key]!=row["prot"]:
                        raise ValueError(f"Assay maps to two proteins: {key}")
                    known[key]=row["prot"]
    # Cross-check may explicitly resolve panel-specific UKB assays. Only exact
    # discovered names are added; no wildcard prefix is guessed here.
    if args.candidate_loci:
        for r in read_tsv(args.candidate_loci):
            d=r["DATASET"].lower();d="ukb" if d=="prot" else d
            if d in DATASETS:
                key=(d,safe_assay(r["ASSAY"]))
                if key in known and known[key]!=r["PROT"]:
                    raise ValueError(f"Candidate/mapping disagreement: {key}")
                known[key]=r["PROT"]
    rows=[];issues=[]
    extensions=(".gz",".tsv.gz",".txt.gz",".bgz",".tsv",".txt")
    for (d,assay),prot in sorted(known.items()):
        folder=args.gwas_root/DATASETS[d]
        directories=(folder/category/assay/"gwas",folder/"clean",folder)
        choices=sorted({(directory/(assay+ext)).resolve() for directory in directories
                        for ext in extensions if (directory/(assay+ext)).is_file()})
        if len(choices)!=1:
            issues.append({"DATASET":d,"ASSAY":assay,"PROT":prot,"STATUS":"missing" if not choices else "ambiguous_paths",
                           "PATHS":";".join(map(str,choices)) or "."})
            continue
        path=choices[0]
        receipts,warnings=local_build_receipts(path,folder,category,assay)
        for status,warning_path in warnings:
            issues.append({"DATASET":d,"ASSAY":assay,"PROT":prot,"STATUS":status,"PATHS":warning_path})
        if (d,assay) in assay_builds:
            b,e,receipt_path=assay_builds[d,assay]
            if receipt_path and receipt_path!=str(path):
                raise ValueError(f"Build receipt PATH does not match selected full summary: {d}/{assay}")
            receipts.append((b,e))
        if d in builds:receipts.append(builds[d])
        if not receipts:
            raise ValueError(f"Full summary found for {d}/{assay}, but no fresh per-file or explicit build receipt; pass --dataset-builds or --build-receipt")
        if len({r[0] for r in receipts})!=1:
            raise ValueError(f"Conflicting build receipts for {d}/{assay}: {receipts}")
        rows.append({"DATASET":d,"ASSAY":assay,"PROT":prot,"PATH":str(path),"BUILD":receipts[0][0],
                     "ANCESTRY":ancestry.get(d,"UNKNOWN").upper(),"EFFECT_SCALE":"unknown","INPUT_SCOPE":"full",
                     "SCOPE_EVIDENCE":"standard_full_summary_path; verify upstream file is complete",
                     "BUILD_EVIDENCE":";".join(e for _,e in receipts)})
    fields=["DATASET","ASSAY","PROT","PATH","BUILD","ANCESTRY","EFFECT_SCALE","INPUT_SCOPE","SCOPE_EVIDENCE","BUILD_EVIDENCE"]
    atomic_tsv(args.export_manifest,rows,fields)
    atomic_tsv(Path(str(args.export_manifest)+".issues.tsv"),issues,["DATASET","ASSAY","PROT","STATUS","PATHS"])
    log(f"Manifest exported: {args.export_manifest}; unique full paths={len(rows)}, issues={len(issues)}. No GWAS matrix has been built.")


def read_candidates(path, assays, discovery, protein_list=None):
    rows=read_tsv(path)
    required={"DATASET","ASSAY","PROT","CHR","POS","SNP","P"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Candidate loci need {sorted(required)}")
    allowed=None
    if protein_list:
        allowed={line.strip().split()[0] for line in protein_list.read_text().splitlines()
                 if line.strip() and not line.lstrip().startswith("#") and line.strip().upper() not in {"PROT","PROTEIN"}}
        if not allowed:
            raise ValueError("Empty discovery protein list")
    known={(a["DATASET"],a["ASSAY"]):a for a in assays}
    by_assay=defaultdict(list);qc=Counter()
    for r in rows:
        d=r["DATASET"].lower();d="ukb" if d=="prot" else d
        if d!=discovery or (allowed is not None and r["PROT"] not in allowed):
            continue
        key=(d,r["ASSAY"])
        if key not in known:
            qc["candidate_without_full_summary_manifest"]+=1
            continue
        if known[key]["PROT"]!=r["PROT"]:
            raise ValueError(f"Candidate/manifest PROT mismatch: {key}")
        try:
            p=float(r["P"]);pos=int(r["POS"])
            if not math.isfinite(p) or not 0<=p<=1 or pos<1:raise ValueError()
        except ValueError:
            qc["invalid_candidate_row"]+=1
            continue
        by_assay[key].append({**r,"CHR":normal_chr(r["CHR"]),"POS":pos})
    if not by_assay:
        raise ValueError("No discovery candidate rows have a matching complete-summary manifest entry")
    return by_assay,qc,allowed


def read_locus_bed(args, build):
    if not args.locus_bed:return []
    lines=[line for line in args.locus_bed.read_text().splitlines() if line.strip() and not line.startswith("#")]
    if not lines:raise ValueError("Empty locus BED")
    first=lines[0].split()
    if first[0].upper() in {"CHR","CHROM","CHROMOSOME"}:
        header=[x.upper() for x in first];data=[dict(zip(header,line.split())) for line in lines[1:]]
    else:
        data=[dict(zip(["CHR","START","END","LOCUS","LD_BLOCK"],line.split())) for line in lines]
    out=[]
    for row in data:
        b=row.get("BUILD") or args.locus_build
        if b is None or normal_build(b)!=build:
            raise ValueError("Locus BED needs a matching explicit BUILD column or --locus-build")
        start,end=int(row["START"]),int(row["END"])
        if args.locus_bed_coordinates=="bed0":start+=1
        if start<1 or end<start:raise ValueError("Invalid locus BED interval")
        chrom=normal_chr(row["CHR"])
        locus=row.get("LOCUS") or f"{chrom}:{start}-{end}"
        out.append({"CHR":chrom,"START":start,"END":end,"LOCUS":locus,
                    "LD_BLOCK":row.get("LD_BLOCK") or locus,"BLOCK_METHOD":"provided_intervals; LD provenance must be documented"})
    out.sort(key=lambda r:(chrom_order(r["CHR"]),r["START"]))
    for a,b in zip(out,out[1:]):
        if a["CHR"]==b["CHR"] and b["START"]<=a["END"]:
            raise ValueError("Overlapping locus BED intervals are ambiguous")
    if len({r['LOCUS'] for r in out})!=len(out):raise ValueError("Duplicate LOCUS in BED")
    return out


def locus_for(record,args,bed):
    c,p=record["chr"],record["pos"]
    if bed:
        return next((r for r in bed if r["CHR"]==c and r["START"]<=p<=r["END"]),None)
    start=((p-1)//args.bin_size)*args.bin_size+1
    group=((p-1)//args.block_size)*args.block_size+1
    block=f"{c}:{group}-{group+args.block_size-1}"
    if c=="6" and 25000000<=p<=35000000:block="6:MHC_25-35Mb"
    return {"CHR":c,"START":start,"END":start+args.bin_size-1,"LOCUS":f"{c}:{start}-{start+args.bin_size-1}",
            "LD_BLOCK":block,"BLOCK_METHOD":"physical_5Mb_groups_with_MHC; not LD independence" if args.block_size==5000000 else "physical_groups_with_MHC; not LD independence"}


def discovery_panel(args,assays,candidates,build):
    bed=read_locus_bed(args,build)
    by_key={(a["DATASET"],a["ASSAY"]):a for a in assays}
    pool=[];qc=[]
    cache_dir=(getattr(args, "state_dir", None) or args.out/"state")/"discovery";cache_dir.mkdir(parents=True,exist_ok=True)
    for k,(key,rr) in enumerate(sorted(candidates.items()),1):
        a=by_key[key]
        positions={(r["CHR"],r["POS"]) for r in rr}
        signature=digest({"schema":SCHEMA,"implementation":IMPLEMENTATION,"stage":"anchors","assay":a,"input":summary_identity(a),"candidates":rr,
                          "p":args.p_threshold,"pal_maf":args.palindrome_maf})
        cached=cache_dir/(digest(key)[:20]+".json")
        obj=json.loads(cached.read_text()) if cached.exists() and args.replace.upper()!="TRUE" else {}
        if obj.get("signature")==signature:
            selected=obj["selected"]
            for s in selected:
                for field in ("eaf","n"):
                    if s[field] is None:s[field]=math.nan
            stat=obj["qc"];stat["checkpoint_reused"]=True
        else:
            records,extraction=extract_positions(a,positions,args.allow_scan,args.tabix)
            counts=Counter();selected=[];seen=set()
            for r in rr:
                hit,why=discovery_record(records.get((r["CHR"],r["POS"]),[]),r["SNP"])
                if hit is None:counts[why]+=1;continue
                if hit["p"]>args.p_threshold:counts["full_summary_not_discovery_significant"]+=1;continue
                if palindromic(hit["ea"],hit["nea"]):
                    f=hit["eaf"]
                    if not math.isfinite(f) or min(f,1-f)>args.palindrome_maf:
                        counts["palindrome_missing_or_near_half_frequency"]+=1;continue
                identity=(hit["chr"],hit["pos"],hit["ea"],hit["nea"])
                if identity in seen:continue
                seen.add(identity)
                selected.append({**hit,"prot":a["PROT"],"assay":a["ASSAY"],"ancestry":a["ANCESTRY"]})
            stat={"DATASET":key[0],"ASSAY":key[1],"PROT":a["PROT"],"candidate_rows":len(rr),"validated_anchors":len(selected),
                  "extraction":extraction,"rejections":dict(counts),"checkpoint_reused":False}
            atomic_json(cached,{"signature":signature,"selected":clean_json(selected),"qc":stat})
        pool.extend(selected);qc.append(stat)
        if k%100==0 or k==len(candidates):log(f"Discovery anchor extraction {k}/{len(candidates)}; validated candidate rows={len(pool)}")
    bins=defaultdict(list);outside=0
    for s in pool:
        loc=locus_for(s,args,bed)
        if loc is None:outside+=1;continue
        bins[loc["LOCUS"]].append((s,loc))
    panel=[]
    for name,members in bins.items():
        hit,loc=min(members,key=lambda v:(v[0]["p"],-abs(v[0]["z"]),v[0]["pos"],v[0]["ea"],v[0]["nea"],v[0]["assay"]))
        exact=[v[0] for v in members if (v[0]["pos"],v[0]["ea"],v[0]["nea"])==(hit["pos"],hit["ea"],hit["nea"])]
        panel.append({"SOURCE_INDEX":0,"LOCUS":name,"CHR":hit["chr"],"POS":hit["pos"],"EA":hit["ea"],"NEA":hit["nea"],
                      "EAF":hit["eaf"],"SNP":hit["snp"],"BUILD":build,"LD_BLOCK":loc["LD_BLOCK"],"BLOCK_METHOD":loc["BLOCK_METHOD"],
                      "SOURCE_PROTS":";".join(sorted({s["prot"] for s in exact})),"SOURCE_ASSAYS":";".join(sorted({s["assay"] for s in exact})),
                      "CANDIDATE_PROTS":";".join(sorted({s["prot"] for s,_ in members})),"PANEL_SOURCE":args.discovery,
                      "ANCESTRY":hit["ancestry"],"DISCOVERY_P":hit["p_text"],"DISCOVERY_Z":hit["z"],
                      "LOCUS_START":loc["START"],"LOCUS_END":loc["END"]})
    panel.sort(key=lambda s:(chrom_order(s["CHR"]),s["POS"],s["EA"],s["NEA"]))
    for i,s in enumerate(panel):s["SOURCE_INDEX"]=i
    if not panel:raise ValueError("No validated discovery anchors; inspect full-summary paths, alleles and source build")
    return panel,{"assays":qc,"candidate_rows_outside_locus_bed":outside,"n_validated_rows":len(pool),"n_fixed_anchors":len(panel)}


def annotation_lookup(path,build):
    out=defaultdict(list)
    for row in read_tsv(path):
        r={k.upper():v for k,v in row.items()}
        prot=r.get("PROT") or r.get("PROTEIN")
        if not prot:raise ValueError("Protein annotation needs prot/PROT")
        if r.get("BUILD") not in {None,"","."} and normal_build(r["BUILD"])!=build:
            raise ValueError(f"Protein annotation build mismatch for {prot}")
        try:
            start,end=int(r["START"]),int(r["END"])
            if start<1 or end<start:continue
            c=normal_chr(r["CHR"])
            if c not in {*(str(x) for x in range(1,23)),"X","Y","MT"}:continue
        except (KeyError,ValueError):continue
        out[prot].append({"CHR":c,"START":start,"END":end,"GENE":r.get("GENE_SYMBOL") or r.get("GENE") or prot,
                          "GENE_ID":r.get("GENE_ID","."),
                          "ANNOTATION_STATUS":r.get("ANNOTATION_STATUS","provided_coordinate"),
                          "COORDINATE_SOURCE":r.get("COORDINATE_SOURCE",r.get("SOURCE","provided")),
                          "COORDINATE_STATUS":r.get("COORDINATE_STATUS","provided_annotation"),
                          "IDENTITY_VERIFIED":r.get("IDENTITY_VERIFIED",""),
                          "SHARED_ANCHOR_TARGETS":r.get("SHARED_ANCHOR_TARGETS",".")})
    return out


def cross_tiers(folder):
    if folder is None:return {}
    result={}
    for name in ("cross-check.tsv","protein_evidence.tsv","protein-evidence.tsv"):
        p=folder/name
        if p.exists():
            for raw in read_tsv(p):
                r={k.upper():v for k,v in raw.items()}
                prot=r.get("PROT") or r.get("PROTEIN")
                if prot:result[prot]=r.get("TIER") or r.get("EVIDENCE_TIER") or r.get("CONFIDENCE_TIER") or "."
            break
    return result


def target_record(a,index,annotation,tiers):
    ann=annotation.get(a["PROT"],[])
    cis_known=bool(ann) and not any(r["ANNOTATION_STATUS"] in {"complex_partial","unresolved","point_unverified"}
        or str(r.get("IDENTITY_VERIFIED","")).lower() in {"0","false"} for r in ann)
    genes={r["GENE"] for r in ann}
    # Keep assays/proteoforms distinct as observations, while aliases and products
    # of the same coding gene stay in one held-out fold and bootstrap group.
    group="GENE:"+next(iter(genes)) if cis_known and len(genes)==1 and all(r["ANNOTATION_STATUS"]=="single_gene" for r in ann) else "PROT:"+a["PROT"]
    return {"TARGET_INDEX":index,"PROT":a["PROT"],"PROTEIN":a["PROT"],"ASSAY":a["ASSAY"],"DATASET":a["DATASET"],
            "PROT_GROUP":group,
            "PATH":a["PATH"],"BUILD":a["BUILD"],"ANCESTRY":a["ANCESTRY"],"EFFECT_SCALE":a.get("EFFECT_SCALE","unknown"),
            "CHR":";".join(dict.fromkeys(r["CHR"] for r in ann)) or ".",
            "START":";".join(str(r["START"]) for r in ann) or ".","END":";".join(str(r["END"]) for r in ann) or ".",
            "GENE":";".join(dict.fromkeys(r["GENE"] for r in ann)) or ".",
            "GENE_ID":";".join(dict.fromkeys(r.get("GENE_ID",".") for r in ann)) or ".",
            "TARGET_INTERVALS":";".join(f'{r["CHR"]}:{r["START"]}-{r["END"]}' for r in ann) or ".",
            "ANNOTATION_STATUS":";".join(dict.fromkeys(r["ANNOTATION_STATUS"] for r in ann)) or "unresolved",
            "COORDINATE_SOURCE":";".join(dict.fromkeys(r.get("COORDINATE_SOURCE","provided") for r in ann)) or ".",
            "COORDINATE_STATUS":";".join(dict.fromkeys(r.get("COORDINATE_STATUS","provided_annotation") for r in ann)) or "unresolved",
            "SHARED_ANCHOR_TARGETS":";".join(dict.fromkeys(r.get("SHARED_ANCHOR_TARGETS",".") for r in ann)) or ".",
            "CIS_KNOWN":int(cis_known),"CROSS_TIER":tiers.get(a["PROT"],"."),"INPUT_SCOPE":"full"}


def build_dataset(args,d,assays,panel,source_hash,annotation,tiers):
    directory=args.out/d;directory.mkdir(parents=True,exist_ok=True)
    (directory/"prepare.done").unlink(missing_ok=True)
    cache_dir=(getattr(args, "state_dir", None) or args.out/"state")/d;cache_dir.mkdir(parents=True,exist_ok=True)
    positions={(s["CHR"],s["POS"]) for s in panel}
    shape=(len(panel),len(assays));targets=[];qc=[];all_aligned=0
    array_names=("beta","se","z","p","eaf","n","observed","cis")
    with tempfile.TemporaryDirectory(prefix=".matrix-",dir=directory) as temporary:
        temporary=Path(temporary);arrays={}
        for name in array_names:
            boolean=name in {"observed","cis"}
            dtype="bool" if boolean else "float64" if name=="p" else "float32"
            arrays[name]=np.lib.format.open_memmap(temporary/(name+".npy"),mode="w+",dtype=dtype,shape=shape)
            arrays[name][:]=False if boolean else np.nan
        for j,a in enumerate(assays):
            signature=digest({"schema":SCHEMA,"implementation":IMPLEMENTATION,"stage":"matrix","assay":a,"input":summary_identity(a),"source_hash":source_hash,
                              "palindrome_maf":args.palindrome_maf,"palindrome_tolerance":args.palindrome_tolerance})
            cache=cache_dir/(digest((d,a["ASSAY"]))[:20]+".npz")
            reused=False
            if cache.exists() and args.replace.upper()!="TRUE":
                with np.load(cache,allow_pickle=False) as saved:
                    if str(saved["signature"].item())==signature:
                        vectors={name:saved[name].copy() for name in ("beta","se","z","p","eaf","n","observed")}
                        stat=json.loads(str(saved["qc"].item()));reused=True
            if not reused:
                records,extraction=extract_positions(a,positions,args.allow_scan,args.tabix)
                vectors={name:np.full(len(panel),np.nan,dtype=np.float64 if name=="p" else np.float32) for name in ("beta","se","z","p","eaf","n")}
                vectors["observed"]=np.zeros(len(panel),dtype=bool);counts=Counter()
                for i,s in enumerate(panel):
                    hit,why,_=aligned_at_anchor(records.get((s["CHR"],s["POS"]),[]),s,a["ANCESTRY"],
                                               palindrome_maf=args.palindrome_maf,palindrome_tolerance=args.palindrome_tolerance)
                    counts[why]+=1
                    if hit is None:continue
                    for name in ("beta","se","p","eaf","n"):vectors[name][i]=hit[name]
                    vectors["z"][i]=hit["beta"]/hit["se"];vectors["observed"][i]=True
                stat={"DATASET":d,"ASSAY":a["ASSAY"],"PROT":a["PROT"],"N_ANCHORS":len(panel),"N_OBSERVED":int(vectors["observed"].sum()),
                      "alignment_counts":dict(counts),"extraction":extraction}
                fd,name=tempfile.mkstemp(prefix=".checkpoint-",suffix=".npz",dir=cache_dir);os.close(fd)
                try:
                    np.savez_compressed(name,signature=np.array(signature),qc=np.array(json.dumps(stat)),**vectors)
                    os.replace(name,cache)
                finally:Path(name).unlink(missing_ok=True)
            for name,vector in vectors.items():arrays[name][:,j]=vector
            ann=annotation.get(a["PROT"],[])
            for i,s in enumerate(panel):
                arrays["cis"][i,j]=any(s["CHR"]==r["CHR"] and r["START"]-args.cis_flank<=s["POS"]<=r["END"]+args.cis_flank for r in ann)
            targets.append(target_record(a,j,annotation,tiers))
            stat["checkpoint_reused"]=reused;qc.append(stat);all_aligned+=stat["N_OBSERVED"]
            if (j+1)%100==0 or j+1==len(assays):log(f"{d}: aligned {j+1}/{len(assays)} assays; observed edges={all_aligned:,}")
        for arr in arrays.values():arr.flush()
        arrays.clear()
        for name in array_names:os.replace(temporary/(name+".npy"),directory/(name+".npy"))
        np.save(temporary/"target_cis_known.npy",np.array([bool(t["CIS_KNOWN"]) for t in targets]))
        os.replace(temporary/"target_cis_known.npy",directory/"target_cis_known.npy")
    atomic_tsv(directory/"source.tsv",panel)
    atomic_tsv(directory/"target.tsv",targets)
    atomic_json(directory/"input_qc.json",qc)
    atomic_tsv(directory/"assay_qc.tsv",[{"DATASET":q["DATASET"],"ASSAY":q["ASSAY"],"PROT":q["PROT"],"N_ANCHORS":q["N_ANCHORS"],
                                        "N_OBSERVED":q["N_OBSERVED"],"OBSERVED_FRACTION":q["N_OBSERVED"]/len(panel),"CHECKPOINT_REUSED":q["checkpoint_reused"],
                                        "ALIGNMENT_COUNTS":json.dumps(q["alignment_counts"],sort_keys=True)} for q in qc])
    manifest={"schema":"phole-signed-matrix-v2","format":"phole-signed-matrix-v2","dataset":d,"directional":True,"input_scope":"full",
              "synthetic":all(a.get("SYNTHETIC")=="TRUE" for a in assays),
              "panel_completeness":args.panel_completeness,
              "panel_receipt_sha256":file_hash(args.panel_receipt) if args.panel_receipt else None,
              "shape":list(shape),"n_source":len(panel),"n_target":len(assays),"n_observed":all_aligned,"source_hash":source_hash,
              "target_hash":file_hash(directory/"target.tsv"),"build":assays[0]["BUILD"],"discovery":args.discovery,
              "arrays":{name:name+".npy" for name in array_names},"target_cis_known":"target_cis_known.npy",
              "missing_semantics":"unobserved values are NaN with observed=False; absence is never zero",
              "cis_semantics":f"any supplied effective annotation interval +/- {args.cis_flank} bp; unknown or identity-unverified targets require target_cis_known=False",
              "effect_semantics":"BETA and signed BETA/SE oriented to source EA; beta units remain assay specific",
              "full_summary_evidence":"INPUT_SCOPE=full in the explicit manifest; headers and allele/numeric rows validated; completeness cannot be proven from queried sites alone",
              "cache_fingerprint":"summary/index path/size/mtime_ns plus implementation, parameters and source SHA256; not a full GWAS content hash",
              "target_selection":"explicit canonical PROT list applied to every cohort" if getattr(args,"target_protein_list",None) else "all provided full-summary assays; cross-cohort evidence tier is annotation only",
              "target_selection_scope":getattr(args,"target_selection_scope","unknown") if getattr(args,"target_protein_list",None) else "all_provided",
              "protein_annotation_sha256":file_hash(args.protein_annotation),
              "array_checksums":{name:file_hash(directory/(name+".npy")) for name in array_names}}
    atomic_json(directory/"manifest.json",manifest)
    (directory/"prepare.done").write_text(source_hash+"\n")
    return manifest


def main():
    args=arguments()
    if args.export_manifest:
        export_manifest(args);return
    required={"assay_manifest":args.assay_manifest,"candidate_loci":args.candidate_loci,"protein_annotation":args.protein_annotation,"out":args.out}
    if any(v is None for v in required.values()):raise ValueError("Signed prepare requires "+", ".join("--"+k.replace("_","-") for k,v in required.items() if v is None))
    if args.replace.upper() not in {"TRUE","FALSE"}:raise ValueError("--replace must be TRUE or FALSE")
    if not 0<args.p_threshold<=1 or args.cis_flank<0 or args.bin_size<1 or args.block_size<args.bin_size:
        raise ValueError("Invalid P/window/bin/group settings")
    if not 0<args.palindrome_maf<0.5 or not 0<args.palindrome_tolerance<0.5:
        raise ValueError("Invalid palindrome frequency settings")
    assays=validate_manifest(args.assay_manifest)
    target_set=None
    if args.target_protein_list:
        target_set={line.strip() for line in args.target_protein_list.read_text().splitlines()
                    if line.strip() and not line.lstrip().startswith("#")}
        if not target_set or any("\t" in name or name=="." for name in target_set):
            raise ValueError("--target-protein-list requires one canonical PROT per line, without a header")
        known={a["PROT"] for a in assays}
        if target_set-known:raise ValueError("Target PROTs absent from the assay manifest: "+", ".join(sorted(target_set-known)[:20]))
        assays=[a for a in assays if a["PROT"] in target_set]
        if not any(a["DATASET"]==args.discovery for a in assays):raise ValueError("Target list retains no primary assays")
    if len({a["SYNTHETIC"] for a in assays})!=1:
        raise ValueError("Synthetic fixture and real summaries must not be mixed in one matrix")
    builds={a["BUILD"] for a in assays}
    if len(builds)!=1:raise ValueError(f"All summaries and panel must use the same explicit build, got {builds}; perform allele-aware liftover upstream")
    build=next(iter(builds))
    effective=args.cross_check_dir/"coordinates-used.tsv" if args.cross_check_dir else None
    if effective is not None and effective.is_file() and file_hash(effective)!=file_hash(args.protein_annotation):
        raise ValueError("Protein annotation differs from cross-check/coordinates-used.tsv; use the exact effective coordinates from this cross-check")
    if effective is not None and effective.is_file():
        receipt_path=args.cross_check_dir/"provenance.json"
        if not receipt_path.is_file():raise ValueError("Effective coordinates lack cross-check provenance")
        receipt=json.loads(receipt_path.read_text())
        expected=receipt.get("coordinates_used_sha256") or receipt.get("output_sha256",{}).get("coordinates-used.tsv")
        if not expected or expected!=file_hash(effective):raise ValueError("Effective coordinate checksum differs from cross-check provenance")
        if receipt.get("cis_flank")!=args.cis_flank:raise ValueError("Cis flank differs from cross-check; rerun cross-check with the intended window")
    args.out.mkdir(parents=True,exist_ok=True);(args.out/"prepare.done").unlink(missing_ok=True)
    candidates,candidate_qc,allowed=read_candidates(args.candidate_loci,assays,args.discovery,args.discovery_protein_list)
    panel,discovery_qc=discovery_panel(args,assays,candidates,build)
    atomic_tsv(args.out/"source.tsv",panel);source_hash=file_hash(args.out/"source.tsv")
    atomic_json(args.out/"discovery_qc.json",discovery_qc)
    annotation=annotation_lookup(args.protein_annotation,build);tiers=cross_tiers(args.cross_check_dir)
    output={}
    for d in DATASETS:
        selected=[a for a in assays if a["DATASET"]==d]
        if selected:output[d]=build_dataset(args,d,selected,panel,source_hash,annotation,tiers)
    if args.discovery not in output or output[args.discovery]["n_observed"]==0:
        raise ValueError("No aligned discovery observations; signed model cannot run")
    manifest={"schema":"phole-signed-panel-v2","format":"phole-signed-panel-v2","directional":True,"input_scope":"full","build":build,
              "synthetic":all(a.get("SYNTHETIC")=="TRUE" for a in assays),
              "panel_completeness":args.panel_completeness,
              "panel_receipt":str(args.panel_receipt.resolve()) if args.panel_receipt else None,
              "panel_receipt_sha256":file_hash(args.panel_receipt) if args.panel_receipt else None,
              "panel_limit":"Candidate enumeration is conditional and incomplete; this is not the complete significant-locus graph" if args.panel_completeness!="full_summary_candidates" else "Candidate positions came from full-summary exports; physical bin selection still defines a conditional locus panel",
              "source_hash":source_hash,"n_source":len(panel),"discovery":args.discovery,
              "datasets":{d:{"path":str((args.out/d).resolve()),"n_target":m["n_target"],"n_observed":m["n_observed"]} for d,m in output.items()},
              "panel_selection":"strongest validated discovery-only variant per provided locus" if args.locus_bed else "strongest validated discovery-only variant per physical bin; bins are not LD-independent",
              "internal_evaluation_scope":"training-protein-list-selected panel; fit must keep other PROTs held out" if allowed is not None else "conditional on a primary-panel selected using all discovery proteins; not unbiased protein cold-start validation",
              "discovery_protein_list":sorted(allowed) if allowed is not None else None,
              "target_protein_list":sorted(target_set) if target_set is not None else None,
              "target_protein_list_sha256":file_hash(args.target_protein_list) if args.target_protein_list else None,
              "target_selection_scope":args.target_selection_scope if target_set is not None else "all_provided",
              "external_validation_scope":"external cohorts participated in target screening; replication is conditional, not untouched cohort validation" if target_set is not None and args.target_selection_scope in {"external_screened","unknown"}
                  else "independent cohort replication at frozen primary anchors, conditional on declared target selection",
              "cross_check_used_for_selection":target_set is not None and args.target_selection_scope in {"external_screened","unknown"},"cis_flank":args.cis_flank,"bin_size":args.bin_size,"block_size":args.block_size,
              "assay_manifest":str(args.assay_manifest.resolve()),"assay_manifest_sha256":file_hash(args.assay_manifest),
              "candidate_loci":str(args.candidate_loci.resolve()),"candidate_loci_sha256":file_hash(args.candidate_loci),
              "protein_annotation":str(args.protein_annotation.resolve()),"protein_annotation_sha256":file_hash(args.protein_annotation),
              "candidate_qc":dict(candidate_qc),"palindromes":"same known ancestry + unambiguous frequency only; cross-ancestry and unknown ancestry excluded",
              "ancestry":{d:sorted({a["ANCESTRY"] for a in assays if a["DATASET"]==d}) for d in output}}
    atomic_json(args.out/"manifest.json",manifest)
    (args.out/"prepare.done").write_text(source_hash+"\n")
    log(f"Signed matrices ready: {args.out}; {len(panel)} frozen anchors; datasets={','.join(output)}")


if __name__=="__main__":
    try:
        main()
    except (ValueError,FileNotFoundError,RuntimeError) as exc:
        raise SystemExit("ERROR: "+str(exc)) from exc
