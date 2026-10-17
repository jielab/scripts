#!/usr/bin/env python3
from __future__ import annotations
import argparse, fcntl, hashlib, json, os, subprocess, tempfile
from pathlib import Path
import pandas as pd

RELEASE="2026-10-17.c1-lockpaper-v1"
PRIMARY_SELECTOR = "risk_resid_error"
COVERAGES = (.2, .4, .5, .6, .8, 1.0)


def eid_hash(ids):
    """Language-independent SHA256 of sorted unique IDs, one per line."""
    return hashlib.sha256("".join(x + "\n" for x in sorted(set(map(str, ids)))).encode()).hexdigest()


def roster_record(path):
    path = Path(path)
    meta = path.with_suffix(path.suffix + ".protocol.json")
    record = json.loads(meta.read_text())
    if record.get("release") != RELEASE or record.get("roster_sha256") != sha(path):
        raise RuntimeError("Frozen outer roster integrity/protocol mismatch: " + str(path))
    return record


def configure_roster(args, *, settings=None, dry_run=False):
    """Resolve shared CLI/environment settings before inspecting omics or Y."""
    if settings is None:
        from le8_common import shared_analysis_settings
        settings, _, _ = shared_analysis_settings(args)
    for key, value in settings.items():
        setattr(args, key, value)
    roster = Path(args.outer_roster) if args.outer_roster else (
        Path(args.analysis_root) / "_protocol" / f"c1_outer_roster_seed{args.seed}.tsv"
    )
    if not dry_run:
        ensure_roster(ukb_phe=args.ukb_phe or os.getenv("UKB_PHE", "/mnt/d/data/ukb/phe"),
                      output=roster, seed=args.seed, group_file=args.group_file,
                      group_col=args.group_col, r_bin=args.r_bin)
        print(f"[LE8 C1 LOCK] shared outer roster: {roster}; SHA256={sha(roster)}", flush=True)
    else:
        print(f"[LE8 C1 LOCK] planned shared outer roster: {roster}", flush=True)
    args.outer_roster = str(roster)
    os.environ["LE8_OUTER_ROSTER"] = str(roster)
    return roster

def sha(path, block=1<<20):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(block),b""): h.update(b)
    return h.hexdigest()

def stable_hash(x,seed):
    return int.from_bytes(hashlib.sha256(f"{seed}:{x}".encode()).digest()[:8],"big")

def _rds_export(rds,out,r_bin,group_col=""):
    code = r'''
args<-commandArgs(trailingOnly=TRUE)
x<-readRDS(args[1])
if(!is.data.frame(x)||!"eid"%in%names(x))stop("all.rds must be a data.frame with eid")
g<-args[3]
z<-data.frame(eid=as.character(x$eid),stringsAsFactors=FALSE)
if(nzchar(g)){
 if(!g%in%names(x))stop("Group column missing: ",g)
 z$group<-as.character(x[[g]])
}
data.table::fwrite(z,args[2],sep="\t",quote=FALSE,na="NA")
'''
    subprocess.run([r_bin,"--vanilla","-e",code,str(rds),str(out),group_col],check=True)

def _read_group_file(path,r_bin):
    path=Path(path)
    if path.suffix.lower()==".rds":
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/"group.tsv"
            code = r'''
args<-commandArgs(trailingOnly=TRUE)
x<-readRDS(args[1]);if(!is.data.frame(x))stop("group RDS must be data.frame")
if(!all(c("eid","group")%in%names(x)))stop("group file requires eid,group")
data.table::fwrite(data.frame(eid=as.character(x$eid),group=as.character(x$group)),
                   args[2],sep="\t",quote=FALSE,na="NA")
'''
            subprocess.run([r_bin,"--vanilla","-e",code,str(path),str(out)],check=True)
            return pd.read_csv(out,sep="\t",dtype=str)
    sep="\t" if path.suffix.lower() in {".tsv",".txt"} else ","
    return pd.read_csv(path,sep=sep,dtype=str)

def ensure_roster(*,ukb_phe,output,seed=2026,test_fraction=.20,group_file="",group_col="",r_bin="Rscript",force=False):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(str(output) + ".lock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        return _ensure_roster(ukb_phe=ukb_phe, output=output, seed=seed,
                              test_fraction=test_fraction, group_file=group_file,
                              group_col=group_col, r_bin=r_bin, force=force)


def _ensure_roster(*,ukb_phe,output,seed,test_fraction,group_file,group_col,r_bin,force):
    if test_fraction != .20:
        raise ValueError("C1 LOCKPAPER fixes the outer test fraction at 0.20")
    ukb_phe=Path(ukb_phe); output=Path(output)
    allrds=ukb_phe/"Rdata"/"all.rds"
    if not allrds.is_file(): raise FileNotFoundError(allrds)
    meta=output.with_suffix(output.suffix+".protocol.json")
    requested=dict(release=RELEASE,seed=int(seed),test_fraction=float(test_fraction),
                   source=str(allrds.resolve()),source_sha256=sha(allrds),
                   group_file=str(Path(group_file).resolve()) if group_file else "",
                   group_file_sha256=sha(group_file) if group_file else "",
                   group_col=str(group_col or ""))
    if output.exists() != meta.exists() and not force:
        raise RuntimeError("Frozen outer roster or its protocol sidecar is missing")
    if output.is_file() and meta.is_file() and not force:
        old=roster_record(output)
        for k,v in requested.items():
            if old.get(k)!=v: raise RuntimeError(f"Frozen outer roster protocol mismatch for {k}: {old.get(k)!r} != {v!r}")
        return output
    output.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as td:
        td=Path(td); ph=td/"phenotype.tsv"
        _rds_export(allrds,ph,r_bin,group_col)
        d=pd.read_csv(ph,sep="\t",dtype=str)
        if d.eid.isna().any() or d.eid.duplicated().any() or d.eid.str.strip().eq("").any(): raise ValueError("Invalid/duplicate eid in all.rds")
        if group_file:
            g=_read_group_file(group_file,r_bin)
            if not {"eid","group"}<=set(g): raise ValueError("Group file requires eid,group")
            if g.eid.isna().any() or g.eid.duplicated().any() or g.group.isna().any(): raise ValueError("Invalid group file")
            if "group" in d:
                supplied=d.eid.map(g.set_index("eid").group)
                if not d.group.eq(supplied).all(): raise ValueError("Conflicting family definitions in file and phenotype column")
                d=d.drop(columns="group")
            d=d.merge(g[["eid","group"]],on="eid",how="left",validate="one_to_one")
            if d.group.isna().any(): raise ValueError("Group file does not cover all phenotype IDs")
        elif "group" not in d:
            d["group"]=d.eid
        if d.group.isna().any() or d.group.astype(str).str.strip().eq("").any(): raise ValueError("Incomplete group mapping")
        groups=sorted(d.group.astype(str).unique(),key=lambda x:(stable_hash(x,seed),x))
        if len(groups)<2: raise ValueError("At least two independent groups are required")
        ntest=max(1,min(len(groups)-1,round(len(groups)*test_fraction)))
        test=set(groups[:ntest])
        d["role"]=d.group.astype(str).map(lambda g:"test" if g in test else "training")
        if d.groupby("group").role.nunique().max()!=1: raise AssertionError("Family crosses outer roles")
        outtmp=td/"outer.tsv"; d.sort_values("eid")[["eid","role"]].to_csv(outtmp,sep="\t",index=False); os.replace(outtmp,output)
    record={**requested,"N":int(len(d)),"groups":int(d.group.nunique()),
            "test_N":int((d.role=="test").sum()),"training_N":int((d.role=="training").sum()),
            "test_groups":int(len(test)),"roster_sha256":sha(output),
            "group_definition":group_col or ("external eid,group mapping" if group_file else "individual ID; family mapping unavailable"),
            "test_eid_hash":eid_hash(d.loc[d.role=="test","eid"]),
            "rule":"sort groups by SHA256(seed:group); exact prespecified group count assigned to test; outcome/omics blind"}
    tmp=Path(str(meta)+".tmp");tmp.write_text(json.dumps(record,indent=2),encoding="utf-8");os.replace(tmp,meta)
    return output

def main(argv=None):
    p=argparse.ArgumentParser()
    p.add_argument("--ukb-phe",required=True);p.add_argument("--output",required=True)
    p.add_argument("--seed",type=int,default=2026);p.add_argument("--test-fraction",type=float,default=.20)
    p.add_argument("--group-file",default="");p.add_argument("--group-col",default="");p.add_argument("--r-bin",default="Rscript")
    p.add_argument("--force",action="store_true");a=p.parse_args(argv)
    print(ensure_roster(ukb_phe=a.ukb_phe,output=a.output,seed=a.seed,test_fraction=a.test_fraction,
                        group_file=a.group_file,group_col=a.group_col,r_bin=a.r_bin,force=a.force))
if __name__=="__main__": main()
