"""Compare GRCh37 IBDmix exports and genotype recurrent, reference-supported markers."""
from __future__ import annotations

import argparse
import gzip
import os
import re
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from ibdmix_paper import AFR5, AUTOSOMES, CHROMOSOMES, CORE, LENGTHS, ROOT, inputs, union, configured_path

UK_REFS = ("Altai", "Chagyr", "Vindija", "Denisova3", "Denisova25")
KG_REFS = ("Altai.2013", "Denisova.2013")
X_SPAN = (2699520, 154931043)
POPCOUNT = np.array([int(x).bit_count() for x in range(256)], dtype=np.uint8)
TABLES = ("audit", "burden", "benchmark", "control", "windows", "batches", "loci", "concordance", "spatial", "x_ratio", "size_matched", "male_tiles")
BASE = ROOT / "config"


def save_table(frame, path):
    frame.to_csv(path, sep="\t", index=False, na_rep="NA", float_format="%.8g")


def read_bed(path):
    by_chr = defaultdict(list)
    opener = gzip.open if path.suffix in {".gz",".bgz"} else open
    with opener(path,"rt") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            c, a, b, *_ = line.split()
            c = c.removeprefix("chr")
            if c not in CHROMOSOMES:
                continue
            a, b = int(a), int(b)
            if not 0 <= a < b <= LENGTHS["37"][c]:
                raise ValueError(f"Invalid GRCh37 BED interval in {path}: {line.strip()}")
            by_chr[c].append((a, b))
    return {c: np.asarray(list(union(sorted(data))), dtype=np.int64) for c, data in by_chr.items()}


def intersect_mask(lo, hi, regions):
    if len(regions) == 0:
        return np.zeros(len(lo), bool)
    j = np.searchsorted(regions[:, 1], lo, side="right")
    return (j < len(regions)) & (regions[np.minimum(j, len(regions)-1), 0] < hi)


def accessible_bp(lo,hi,regions):
    if len(regions)==0:
        return np.zeros(len(lo),dtype=np.int64)
    starts,ends=regions.T
    cumulative=np.r_[0,np.cumsum(ends-starts)]
    def before(x):
        j=np.searchsorted(ends,x,side="right")
        return cumulative[j]+np.where(j<len(ends),np.maximum(0,x-starts[np.minimum(j,len(starts)-1)]),0)
    return before(hi)-before(lo)


def read_calls(path, chrom, refs):
    from ibdmix_raw import frame
    df = frame(path)
    if len(df) and (not df.chrom.str.removeprefix("chr").replace("23","X").eq(chrom).all()
                    or not df.anc.isin(refs).all() or df.ID.eq("").any()
                    or not np.isfinite(df.slod).all() or not df.start.ge(0).all()
                    or not df.end.le(LENGTHS["37"][chrom]).all() or not df.length.gt(0).all()
                    or not df.length.eq(df.end-df.start).all()):
        raise ValueError(f"Unexpected reference, coordinate, field or build: {path}")
    return df


def batches(spec):
    names = []
    for item in spec.split(","):
        m = re.fullmatch(r"ygb(\d{3})(?:-(?:ygb)?(\d{3}))?", item)
        if not m or int(m[1]) > int(m[2] or m[1]):
            raise ValueError("Use batch names such as ygb001 or ygb001-049")
        names.extend(f"ygb{x:03d}" for x in range(int(m[1]), int(m[2] or m[1])+1))
    if len(names)!=len(set(names)):
        raise ValueError("Duplicate batch in --batches")
    return names


def roster_for(root, name, paths, override):
    if name == "1kg":
        p = pd.read_csv(BASE/"1kg.phase3.panel.tsv", sep="\t", dtype=str)
        return p, p.loc[p.gender.eq("male"),"sample"].to_numpy(), "bundled Phase 3 panel; assumed tested"
    expected = {}
    for c in CHROMOSOMES:
        logs = list((root/"logs").glob(f"chr{c}.*prepare.log"))
        sizes = {int(m) for log in logs for m in re.findall(r'modern_passes=\d+[^\n]*?\bsamples=(\d+)',log.read_text())}
        if len(sizes)!=1:
            raise ValueError(f"{name} chr{c}: missing or inconsistent preparation sample count")
        expected[c] = sizes.pop()
        meta = root/"jobs"/f"chr{c}.run.meta.tsv"
        cfg = dict(row.split("\t",1) for row in meta.read_text().splitlines() if "\t" in row)
        if (cfg.get("genome_build")!="b37" or cfg.get("coordinate_system")!="0-based-half-open"
                or set(cfg.get("refs","").split())!=set(UK_REFS)
                or (c=="X" and cfg.get("x_mode")!="male_haploid_nonpar")):
            raise ValueError(f"{name} chr{c}: incompatible build, refs, mask profile or X mode")
    if override:
        p = pd.read_csv(override/f"{name}.tsv", sep="\t", dtype=str)
        if not {"sample","x_tested"}.issubset(p) or not p.x_tested.isin(["0","1"]).all():
            raise ValueError("Roster needs sample and x_tested (actual inclusion 0/1)")
        xids = p.loc[p.x_tested.eq("1"),"sample"].to_numpy()
        source = "provided tested roster and per-chromosome log counts"
    else:
        def observed(c):
            from ibdmix_raw import frames
            return np.asarray(sorted({iid for d in frames(paths[c]) for iid in d.ID}), dtype=str)
        # The positive IDs exhaust an independently logged denominator. Otherwise require the true roster.
        p = pd.DataFrame({"sample":observed("1")})
        xids = observed("X")
        source = "all IDs recovered by equality with independent preparation counts"
    if (p["sample"].duplicated().any() or not set(xids)<=set(p["sample"])
            or any(expected[c]!=len(p) for c in AUTOSOMES) or len(xids)!=expected["X"]):
        raise ValueError(f"{name}: tested roster unresolved. Supply --roster-dir with actual batch rosters.")
    return p, xids, source


def person_union(idx, a, b, n, clip=None):
    if clip:
        a,b = np.maximum(a,clip[0]), np.minimum(b,clip[1])
        keep = b>a
        idx,a,b = idx[keep],a[keep],b[keep]
    if not len(idx):
        return np.zeros(n, dtype=np.int64)
    order = np.lexsort((a,idx))
    idx,a,b = idx[order],a[order],b[order]
    # Offset each person's coordinates into its own nonoverlapping block.
    shift = idx.astype(np.int64) * 300_000_000
    a,b = a+shift,b+shift
    last = np.maximum.accumulate(b)
    added = np.maximum(0,b-np.maximum(a,np.r_[0,last[:-1]]))
    return np.bincount(idx, weights=added, minlength=n).astype(np.int64)


def tile_hits(matrix, indices, lo, hi, bit, width):
    # A single exported tract must cover at least 80% of a complete tile.
    needed = (4*width+4)//5
    long_enough = hi-lo >= needed
    indices, lo, hi = indices[long_enough], lo[long_enough], hi[long_enough]
    first = np.maximum(0, -np.floor_divide(-(lo+needed-width), width))
    stop = np.minimum(len(matrix), (hi-needed)//width+1)
    sizes = np.maximum(0,stop-first)
    if sizes.sum():
        offsets = np.repeat(first - (np.cumsum(sizes)-sizes),sizes)
        tile = offsets + np.arange(sizes.sum())
        matrix[tile, np.repeat(indices,sizes)] |= bit


def masks(bits):
    nea, den = bits & 7, bits & 24
    return {"N": (POPCOUNT[nea]>=2) & (den==0), "D": (den==24) & (nea==0)}


def run_batch(name, args, stage, gaps, callable_bed, tables, seen):
    root = args.dir_ibdmix/name
    paths = inputs(root/"raw")
    if set(paths)!=set(CHROMOSOMES):
        raise ValueError(f"{name}: need all 22 autosomes and X")
    kg = name=="1kg"
    refs = KG_REFS if kg else UK_REFS
    roster, xids, provenance = roster_for(root,name,paths,args.roster_dir)
    raw_ids = roster["sample"].to_numpy(str)
    raw_xids = np.asarray(xids,dtype=str)
    sentinels = set(roster.loc[roster["sample"].str.fullmatch(r"-\d+"),"sample"]) if not kg else set()
    roster = roster[~roster["sample"].isin(sentinels)].reset_index(drop=True)
    xids = [i for i in xids if i not in sentinels]
    ids = roster["sample"].to_numpy(str)
    pp = pd.Index(ids)
    if not kg:
        if seen.intersection(ids):
            raise ValueError(f"{name}: overlapping individuals across YGB batches")
        seen.update(ids)
    xids = np.asarray(xids,dtype=str)
    xmask = np.isin(ids,xids)
    if not kg:
        private = stage/"genotypes"/name
        private.mkdir(parents=True, mode=0o700)
        np.save(private/"autosome.samples.npy", ids, allow_pickle=False)
        np.save(private/"X.samples.npy", xids, allow_pickle=False)
    groups = {"ALL":np.ones(len(ids),bool),"Xtested":xmask}
    if kg:
        groups.update({s:roster.super_pop.eq(s).to_numpy() for s in sorted(set(roster.super_pop))})
        groups["AFR5"] = roster["pop"].isin(AFR5).to_numpy()
        groups["GBR"] = roster["pop"].eq("GBR").to_numpy()
    totals = defaultdict(lambda:np.zeros(len(ids),np.int64))
    global_union = 0
    manifest = {"name":name, "N":len(ids),"N_X":len(xids),"roster":provenance,"prepared_N":len(raw_ids),
                "prepared_N_X":len(raw_xids),"excluded_negative_IDs":len(sentinels)}
    for chrom in CHROMOSOMES:
        cx = chrom=="X"
        current = xids if cx else ids
        index = pd.Index(current)
        local_to_full = pp.get_indexer(current)
        n = len(current)
        df = read_calls(paths[chrom],chrom,refs)
        raw_rows = df.anc.value_counts()
        if not kg and args.roster_dir is None and df.ID.nunique()!=len(raw_xids if cx else raw_ids):
            raise ValueError(f"{name} chr{chrom}: unresolved sample identities; supply --roster-dir")
        df = df[~df.ID.isin(sentinels)]
        df["i"] = index.get_indexer(df.ID)
        if (df.i<0).any():
            raise ValueError(f"{name} chr{chrom}: segment ID outside tested roster")
        if not kg and args.roster_dir is None and df.i.nunique()!=n:
            raise ValueError(f"{name} chr{chrom}: cannot verify roster from positives; supply --roster-dir")
        w = args.window_bp
        left = np.arange(LENGTHS["37"][chrom]//w)*w
        right = left+w
        regions = gaps.get(chrom,[])
        eligible = ~intersect_mask(left,right,regions)
        if cx:
            eligible &= (left>=X_SPAN[0]) & (right<=X_SPAN[1])
        if callable_bed is not None:
            eligible &= 5*accessible_bp(left,right,callable_bed.get(chrom,[])) >= 4*w
        bits = np.zeros((len(left),n),np.uint8)
        strict = np.zeros_like(bits)
        scope = "X" if cx else "A"
        for k,ref in enumerate(refs):
            part = df[df.anc==ref]
            qualified = part.length.ge(50000) & part.slod.ge(4)
            r = part[qualified]
            ix = r.i.to_numpy()
            lo,hi = r.start.to_numpy(), r.end.to_numpy()
            length = hi-lo
            bridges = intersect_mask(lo,hi,regions)
            kept = ~bridges
            hi_lod = r.slod.to_numpy()>=6
            if cx and (np.any(lo<X_SPAN[0]) or np.any(hi>X_SPAN[1])):
                raise ValueError(f"{name}: X calls outside common central non-PAR analysis span")
            tables["audit"].append(dict(batch=name,chrom=chrom,ref=ref,n_tested=n,
                n_rows=int(raw_rows.get(ref,0)),non_sentinel_rows=len(part),pass_length_lod=len(r),observed_carriers=len(np.unique(ix)),
                gap_touching_rows=int(bridges.sum()),longest_bp=int(length.max(initial=0)),
                longest_start0=int(lo[length.argmax()]) if len(lo) else None,
                longest_end0=int(hi[length.argmax()]) if len(hi) else None,
                gap_free_max_bp=int(length[kept].max(initial=0))))
            for label,select in [("export",np.ones(len(r),bool)),("gap_free",kept)]:
                bp = person_union(ix[select],lo[select],hi[select],n)
                totals[scope,ref,label][local_to_full] += bp
            if kg and not cx and ref=="Altai.2013":
                global_union += int(person_union(np.zeros(len(r),int),lo,hi,1)[0])
            if chrom=="3":
                covered = person_union(ix[kept],lo[kept],hi[kept],n,CORE)
                for group,members in groups.items():
                    hit = covered[members]
                    tables["control"].append(dict(batch=name,group=group,ref=ref,n=len(hit),
                        any_bp=int((hit>0).sum()),core80=int((5*hit>=4*(CORE[1]-CORE[0])).sum())))
            tile_hits(bits,ix[kept],lo[kept],hi[kept],1<<k,w)
            s = kept & hi_lod
            tile_hits(strict,ix[s],lo[s],hi[s],1<<k,w)
        win = pd.DataFrame({"batch":name,"chrom":chrom,"tile":np.arange(len(left)),
                            "eligible":eligible,"n":n})
        if kg:
            for group,member in groups.items():
                m = member[local_to_full]
                win[f"n_{group}"] = int(m.sum())
                for k,key in enumerate(("alt","den")):
                    win[f"{key}_{group}"] = np.sum((bits[:,m] & (1<<k))!=0, axis=1)
        else:
            primary = masks(bits)
            strict_hits = masks(strict)
            packed = {}
            for k,key in enumerate(UK_REFS):
                win[key] = ((bits & (1<<k))!=0).sum(axis=1)
            for label,field in [("N",7),("D",24)]:
                any_hit = (bits & field)!=0
                yes = primary[label]
                unknown = any_hit & ~yes
                # Stronger LOD required for supporting refs; keep competing-reference veto at LOD4.
                strict_yes = strict_hits[label] & yes
                win[label] = yes.sum(axis=1)
                win[label+"_LOD6"] = strict_yes.sum(axis=1)
                win[label+"_uncertain"] = unknown.sum(axis=1)
                for suffix,state in [("carrier",yes),("uncertain",unknown)]:
                    packed[label+"_"+suffix] = np.packbits(state,axis=1,bitorder="little")
            win["mixed"] = (((bits & 7)!=0)&((bits & 24)!=0)).sum(axis=1)
            win["any_ref"] = (bits!=0).sum(axis=1)
            np.savez_compressed(private/f"chr{chrom}.npz", **packed)
        men = xmask[local_to_full]
        for k,ref in enumerate(refs):
            hits = ((bits[eligible][:,men] & (1<<k))!=0).sum()
            tables["male_tiles"].append(dict(batch=name,scope=scope,chrom=chrom,ref=ref,
                carrier_tiles=int(hits),tested_person_tiles=int(eligible.sum())*int(men.sum())))
        tables["windows"].append(win)
        print(f"{name} chr{chrom}: {len(df):,} rows", flush=True)
    for (scope,ref,kind),bp in totals.items():
        # Physical footprint is intentionally separate from callable-normalized tile comparisons.
        span = X_SPAN[1]-X_SPAN[0] if scope=="X" else sum(LENGTHS["37"][c] for c in AUTOSOMES)
        for group,which in groups.items():
            members = which & (xmask if scope=="X" else True)
            values = bp[members]
            tables["burden"].append(dict(batch=name,scope=scope,group=group,ref=ref,filter=kind,
                n=len(values),total_bp=int(values.sum()),span_bp=span,
                mean_Mb=np.mean(values)/1e6,median_Mb=np.median(values)/1e6,physical_fraction=np.mean(values)/span))
    if kg:
        total = totals["A","Altai.2013","export"]
        for measure,value,published in [
            ("Person-summed Gb",total.sum()/1e9,110.98),("Union Gb",global_union/1e9,1.29),
            *[(group+" Mb/person",total[groups[group]].mean()/1e6,orig)
              for group,orig in [("AFR5",17),("EUR",51),("EAS",55),("SAS",55)]]]:
            tables["benchmark"].append(dict(measure=measure,observed=value,Chen_2020=published,ratio=value/published))
    return manifest


def main():
    p = argparse.ArgumentParser(prog="gu.paper.sh ibdmix", description=__doc__,
        epilog="For a single-cohort segment report use --target NAME --dir-ibdmix RAW.")
    p.add_argument("--dir-ibdmix",type=Path,default=Path(configured_path("ANALYSIS_ROOT","gu","ibdmix")))
    p.add_argument("--dir-out",type=Path,default=Path(configured_path("ANALYSIS_ROOT","gu","ibdmix","paper")),
                   help="Report directory [ANALYSIS_ROOT/gu/ibdmix/paper]")
    p.add_argument("--batches",default="ygb001",help="YGB batches, e.g. ygb001-049 [ygb001]")
    p.add_argument("--roster-dir",type=Path)
    p.add_argument("--callable-bed",type=Path,help="GRCh37 common callable BED shared by all selected datasets/references")
    p.add_argument("--window-bp",type=int,default=50000)
    p.add_argument("--min-carriers",type=int,default=5)
    p.add_argument("--source",default="",help="Input release or commit identifiers for provenance")
    args = p.parse_args()
    if not 10000<=args.window_bp<=200000 or args.min_carriers<1:
        p.error("window must be 10–200 kb and min-carriers positive")
    names = batches(args.batches)
    for name in ["1kg",*names]:
        if not (args.dir_ibdmix/name/"raw").is_dir():
            p.error(f"missing {name}/raw in {args.dir_ibdmix}")
    args.dir_out.mkdir(parents=True, exist_ok=True)
    dest = args.dir_out/"asv"
    # Serialize replacement, build privately, then publish a real directory.
    import fcntl
    from shutil import rmtree
    with (args.dir_out/".asv.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if dest.is_symlink():
            p.error(f"Migrate legacy ASV symlink before running: {dest}")
        stage = args.dir_out/".asv.build"
        previous = args.dir_out/".asv.previous"
        if previous.exists():
            if not dest.exists(): os.replace(previous,dest)
            else: rmtree(previous)
        if stage.exists(): rmtree(stage)
        stage.mkdir()
        try:
            tables = {k:[] for k in TABLES}
            gaps = read_bed(BASE/"asv.hg19.major_gaps.bed")
            callable_bed = read_bed(args.callable_bed) if args.callable_bed else None
            seen = set()
            manifests = [run_batch(name,args,stage,gaps,callable_bed,tables,seen) for name in ["1kg",*names]]
            from asv_report import finish
            finish(stage,args,manifests,tables)
            if dest.exists(): os.replace(dest,previous)
            try: os.replace(stage,dest)
            except BaseException:
                if previous.exists(): os.replace(previous,dest)
                raise
            if previous.exists(): rmtree(previous)
            (args.dir_out/"REBUILD_REQUIRED.json").unlink(missing_ok=True)
            print(dest/"asv.html")
        except BaseException:
            if stage.exists(): rmtree(stage)
            raise



if __name__=="__main__":
    main()
