#!/usr/bin/env python3
"""GRID ABM LD-block preparation.

Supported boundary layouts
--------------------------
1. Four ancestry-specific GRCh37 maps:
   AFR.37.bed,EAS.37.bed,EUR.37.bed,SAS.37.bed
   Each population's PRS-CSx posterior weights are aggregated within its own map.
   Output tensor: N x sum(blocks) x 1.

2. One shared map:
   ALL.37.bed
   All four PRS-CSx population-specific posterior weights are aggregated within the
   same blocks. Output tensor: N x blocks x 4, channels AFR/EAS/EUR/SAS.

The shared-map implementation therefore preserves ancestry-specific PRS information
without forcing four different block boundaries into one token sequence.

Current GRID genotypes/PRS-CSx weights are GRCh37. GRCh38 maps are rejected here.
"""
from __future__ import annotations
import argparse, hashlib, json, os, re, warnings, shutil, time
from pathlib import Path
import numpy as np
import pandas as pd
import block_scoring

POPS = ("AFR","EAS","EUR","SAS")
PRODUCER_VERSION = "csx-assigned-dosage-blocks.1"
BINDING_NAMESPACE = "producer-bound"
# Keep one digest per file for the whole run. A bounded LRU evicted the 132 GiB
# genotype digest while hashing the 44 PVAR/PSAM metadata files for each trait.
_DIGESTS = {}



def file_identity(path):
    stat = Path(path).stat()
    return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)


def sha256(path):
    h = hashlib.sha256()
    size = Path(path).stat().st_size
    large = size >= 1024**3
    started = last = time.monotonic()
    done = 0
    if large:
        print(f"[GRID ABM] START SHA256 {path}: {size/2**30:.2f} GiB (disk/CPU)", flush=True)
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(8*1024*1024), b""):
            h.update(b)
            done += len(b)
            now = time.monotonic()
            if large and now-last >= 60:
                print(f"[GRID ABM] PROGRESS SHA256 {Path(path).name}: {done/2**30:.2f}/{size/2**30:.2f} GiB ({100*done/size:.1f}%)", flush=True)
                last = now
    if large:
        print(f"[GRID ABM] DONE SHA256 {Path(path).name}: {time.monotonic()-started:.1f}s", flush=True)
    return h.hexdigest()


def sha256_cached(path):
    """Reuse a recorded legacy digest while resolved path, size and mtime agree."""
    p = Path(path).resolve()
    identity = file_identity(p)
    from runtime_state import file_identity as disk_identity
    stamp = disk_identity(p)
    marker = p.with_name(p.name + ".verification.json")
    if marker.is_file():
        saved = json.loads(marker.read_text())
        if saved.get("identity") == stamp and saved.get("file_generation") == list(identity):
            return saved["legacy_digest"]
    # A legacy completed supplement already binds its shape, exact cohort and
    # SNP order; adopt its recorded digest without rescanning dosage bytes.
    complete = p.parent / "COMPLETE.json"
    if p.name == "genotypes.npy" and complete.is_file() and not marker.exists():
        receipt = json.loads(complete.read_text())
        if receipt.get("genotype_sha256") and p.stat().st_mtime_ns <= complete.stat().st_mtime_ns:
            array = np.load(p, mmap_mode="r", allow_pickle=False)
            if list(array.shape) != receipt.get("shape"):
                raise ValueError("Completed dosage supplement shape changed")
            array._mmap.close()
            digest = receipt["genotype_sha256"]
            json_write(marker, dict(identity=stamp, file_generation=list(identity), legacy_digest=digest, basis="completed_producer_shape_and_metadata"))
            return digest
    cached = _DIGESTS.get(p)
    if cached is not None and cached[0] == identity:
        if identity[2] >= 1024**3:
            print(f"[GRID ABM] REUSE SHA256 {p.name}: {identity[2]/2**30:.2f} GiB; unchanged in this run", flush=True)
        return cached[1]
    digest = sha256(p)
    if file_identity(p) != identity:
        raise RuntimeError(f"File changed during SHA256 verification: {p}")
    _DIGESTS[p] = (identity, digest)
    if identity[2] >= 64 * 1024**2:
        json_write(marker, dict(identity=stamp, file_generation=list(identity), legacy_digest=digest))
    return digest


def receipt_matches(receipt, payload, signature):
    """Reuse only an exact current source/input signature."""
    return receipt.get("signature") == signature and receipt.get("signature_payload") == payload


def json_write(path,value):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_name(p.name+".tmp")
    tmp.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False,default=lambda x: str(x) if isinstance(x,Path) else x.item())+"\n")
    os.replace(tmp,p)


def read_config(path):
    c=json.loads(Path(path).read_text())
    if int(c.get("genome_build",37)) != 37:
        raise ValueError("Current GRID ABM supports the GRCh37 analysis only")
    if tuple(c["populations"])!=POPS:
        raise ValueError("Population order must be AFR,EAS,EUR,SAS")
    return c


def resolve_block_files(block_dir, block_files, *, lightweight=False):
    root=Path(block_dir)
    names=[x.strip() for x in str(block_files).split(",") if x.strip()]
    if len(names) not in (1,4):
        raise ValueError("--ld-block-files must contain either one ALL file or four ancestry-specific files")
    allowed = ({"ALL.37.bed"} if len(names)==1 else {"AFR.37.bed","EAS.37.bed","EUR.37.bed","SAS.37.bed"})
    if any(name not in allowed for name in names):
        raise ValueError("Use ALL.37.bed or the four supported GRCh37 population maps")
    resolved=[]
    for name in names:
        p=root/name
        if not p.is_file(): raise FileNotFoundError(p)
        resolved.append(p.resolve())
    if len(names)==4:
        pops=[]
        for p in resolved:
            base=p.name.upper()
            pop=next((x for x in POPS if base.startswith(x+".")),None)
            if pop is None: raise ValueError(f"Cannot infer population from {p.name}")
            pops.append(pop)
        if tuple(pops)!=POPS:
            raise ValueError("Four-file mode must be ordered AFR,EAS,EUR,SAS")
        mode="population_multiview"
    else:
        mode="shared_universal"
    # Current analysis must use GRCh37 boundary files.
    for p in resolved:
        if p.name not in {"AFR.37.bed","EAS.37.bed","EUR.37.bed","SAS.37.bed","ALL.37.bed"}:
            raise ValueError(f"Current GRID input is GRCh37; refusing non-.37 boundary file: {p}")
    if lightweight:
        identities={p.name:dict(resolved_path=str(p),size=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in resolved}
        return dict(mode=mode,files=resolved,requested=names,sha256={},
                    identities=identities,layout_id=mode + "_37_lightweight")
    hashes={p.name:sha256(p) for p in resolved}
    sig=hashlib.sha256(json.dumps({"mode":mode,"sha256":hashes},
                                  sort_keys=True).encode()).hexdigest()
    label=("multi37_" if mode=="population_multiview" else "all37_")+sig[:12]
    return dict(mode=mode,files=resolved,requested=names,sha256=hashes,
                signature=sig,layout_id=label)


def read_bed(path, source_population="ALL"):
    d=pd.read_csv(path,sep=r"\s+",header=None,comment="#",usecols=[0,1,2],
                  names=["chrom","start","end"])
    d["chrom"]=d.chrom.astype(str).str.replace("^chr","",regex=True)
    if not d.chrom.isin([str(i) for i in range(1,23)]).all():
        raise ValueError(f"{path}: non-autosomal chromosome")
    d["CHR"]=d.chrom.astype(int)
    d["start"]=pd.to_numeric(d.start,errors="raise").astype(np.int64)
    d["end"]=pd.to_numeric(d.end,errors="raise").astype(np.int64)
    if (d.start<0).any() or (d.end<=d.start).any():raise ValueError(f"{path}: invalid interval")
    d=d.sort_values(["CHR","start","end"],kind="stable").reset_index(drop=True)
    for c,q in d.groupby("CHR",sort=True):
        if len(q)>1 and np.any(q.start.to_numpy()[1:] < q.end.to_numpy()[:-1]):
            raise ValueError(f"{path}: overlapping LD blocks on chr{c}")
    d["source_population"]=source_population
    d["block_local"]=np.arange(len(d),dtype=np.int32)
    d["block_id"]=[f"{source_population}.chr{c}.b{i+1:04d}" for i,c in enumerate(d.CHR)]
    d["midpoint"]=(d.start+d.end)/2
    d["length"]=d.end-d.start
    return d[["source_population","block_local","block_id","CHR","start","end","midpoint","length"]]


def build_manifest(layout, config):
    frames=[]
    if layout["mode"]=="population_multiview":
        for pop,path in zip(POPS,layout["files"]):
            frames.append(read_bed(path,pop))
        manifest=pd.concat(frames,ignore_index=True)
        manifest["token_index"]=np.arange(len(manifest),dtype=np.int32)
        channels=["scalar"]  # ancestry is recorded on each token
        counts={p:int((manifest.source_population==p).sum()) for p in POPS}
        expected={"AFR":2581,"EAS":1445,"EUR":1703,"SAS":1445}
        mismatch={p:(counts[p],expected[p]) for p in POPS if counts[p]!=expected[p]}
        if mismatch:
            # Do not fail: local files may be updated maps, but make the difference explicit.
            status="noncanonical_2016_counts"
        else: status="matches_expected_local_2016_counts"
    else:
        manifest=read_bed(layout["files"][0],"ALL")
        manifest["token_index"]=np.arange(len(manifest),dtype=np.int32)
        channels=list(POPS)
        counts={"ALL":len(manifest)}
    audit=dict(mode=layout["mode"],layout_id=layout["layout_id"],block_files=[str(p) for p in layout["files"]],
               block_file_sha256=layout["sha256"],counts=counts,
               total_tokens=int(len(manifest)),channels=channels,channel_count=len(channels),
               genome_build=37,layout_mode="shared" if layout["mode"]=="shared_universal" else layout["mode"],
               coordinate_semantics="BED 0-based half-open; PRS BP converted to BP-1")
    audit['files']=[dict(resolved_path=str(path),filename=path.name,
        block_count=counts[pop],sha256=layout['sha256'].get(path.name))
        for pop,path in zip(POPS if layout['mode']=='population_multiview' else ('ALL',),layout['files'])]
    if layout['mode']=='population_multiview':audit['population_count_status']=status
    print('[GRID ABM] LD audit '+json.dumps(audit),flush=True)
    return manifest,channels,audit


def read_weights(path):
    d=pd.read_csv(path,sep=r"\s+",compression="infer",low_memory=False)
    up={str(c).upper():c for c in d.columns}
    def c(*names):return next((up[x] for x in names if x in up),None)
    snp,a1,beta,chrcol,bp=c("SNP","RSID","ID"),c("A1","EA","EFFECT_ALLELE"),c("BETA","EFFECT"),c("CHR","CHROM"),c("BP","POS")
    a2=c("A2","NEA")
    if None in (snp,a1,beta,chrcol,bp):raise ValueError(f"Invalid PRS-CSx weights: {path}")
    z=pd.DataFrame(dict(SNP=d[snp].astype(str),A1=d[a1].astype(str).str.upper(),
                        BETA=pd.to_numeric(d[beta],errors="coerce"),
                        CHR=pd.to_numeric(d[chrcol],errors="coerce"),
                        BP=pd.to_numeric(d[bp],errors="coerce"),
                        A2=d[a2].astype(str).str.upper() if a2 is not None else ""))
    z=z[np.isfinite(z.BETA)&z.CHR.between(1,22)&np.isfinite(z.BP)].copy()
    z["CHR"]=z.CHR.astype(int);z["BP"]=z.BP.astype(np.int64)
    if z.SNP.duplicated().any():raise ValueError(f"Duplicate SNPs in {path}")
    return z.reset_index(drop=True)


def assign_blocks(weights, bed):
    out=np.full(len(weights),-1,dtype=np.int32)
    for chrom in sorted(weights.CHR.unique()):
        wi=np.flatnonzero(weights.CHR.to_numpy()==chrom)
        q=bed[bed.CHR==chrom]
        if q.empty:continue
        starts=q.start.to_numpy();ends=q.end.to_numpy();blocks=q.block_local.to_numpy()
        pos=weights.BP.to_numpy()[wi]-1
        j=np.searchsorted(starts,pos,side="right")-1
        jj=np.clip(j,0,len(q)-1)
        good=(j>=0)&(j<len(q))&(pos<ends[jj])
        out[wi[good]]=blocks[jj[good]]
    return out


def locate_weights(config,trait, *, lightweight=False):
    prov=Path(config["csx_analysis_root"])/trait/"provenance.json"
    if not prov.is_file():raise FileNotFoundError(prov)
    meta=json.loads(prov.read_text());result={}
    for pop in POPS:
        rec=meta.get("models",{}).get(f"csx.{pop}")
        if not rec or not rec.get("weights_file"):raise ValueError(f"{prov}: no csx.{pop}")
        p=Path(rec["weights_file"])
        if not p.is_file():raise FileNotFoundError(p)
        if not lightweight and rec.get("weights_sha256") and sha256(p)!=rec["weights_sha256"]:
            raise ValueError(f"Published CSx weights changed: {p}")
        result[pop]=p
    return result,meta


def read_rds(path):
    # Native R handles data.table RAW/ALTREP attributes unsupported by rdata.
    # Reuse GRID's reader, which also preserves character IDs and numeric precision.
    from prsformer import grid_data_helpers
    return grid_data_helpers().read_table(path)


def identified(d):
    d=d.copy();c=next((x for x in ("eid","IID","#IID","ID") if x in d),None)
    if c is None:raise ValueError("Missing ID")
    if c!="eid":d=d.rename(columns={c:"eid"})
    d["eid"]=d.eid.astype(str)
    if d.eid.duplicated().any():raise ValueError("Duplicate IDs")
    return d


def prsformer_cache(config):
    root=Path(config["prsformer_cache"])
    required=[root/"prepare.json",root/"data.tsv.gz",root/"variants.tsv.gz",root/"genotypes.npy"]
    for p in required:
        if not p.is_file():raise FileNotFoundError(f"PRSformer preparation cache missing: {p}. Run ./grid.sh abm --module-only --through prepare")
    meta=json.loads((root/"prepare.json").read_text())
    cohort=pd.read_csv(root/"data.tsv.gz",sep="\t",compression="gzip",
                       dtype={"eid":str, "group":str, "family_id":str, "group_id":str},low_memory=False)
    variants=pd.read_csv(root/"variants.tsv.gz",sep="\t",compression="gzip",low_memory=False)
    variants["SNP"]=variants.SNP.astype(str);variants["REF"]=variants.REF.astype(str).str.upper();variants["ALT"]=variants.ALT.astype(str).str.upper()
    if variants.SNP.duplicated().any():raise ValueError("Duplicate PRSformer-cache SNP")
    geno=np.load(root/"genotypes.npy",mmap_mode="r",allow_pickle=False)
    if geno.shape!=(len(cohort),len(variants)):raise ValueError("PRSformer genotype cache shape mismatch")
    return root,meta,cohort,variants,geno


def weight_map_for_population(pop,weights,bed,variants,channel,token_offset=0):
    blocks=assign_blocks(weights,bed)
    assigned=blocks>=0
    mass=weights.BETA.abs().to_numpy(float)
    total=float(mass.sum());covered=float(mass[assigned].sum())
    missing=weights.loc[~assigned,["CHR","SNP","BP","A1","A2","BETA"]]
    qc=dict(weighted_snps_total=len(weights),boundary_assigned_snps=int(assigned.sum()),
        boundary_unassigned_snps=int((~assigned).sum()),
        boundary_assigned_snp_fraction=float(assigned.mean()) if len(weights) else 0.,
        boundary_unassigned_snp_fraction=float((~assigned).mean()) if len(weights) else 0.,
        total_abs_beta_mass=total,boundary_assigned_abs_beta_mass=covered,
        boundary_unassigned_abs_beta_mass=total-covered,
        boundary_assigned_abs_beta_fraction=covered/total if total else 1.,
        boundary_unassigned_abs_beta_fraction=(total-covered)/total if total else 0.,
        boundary_unassigned_by_chr=json.dumps({str(k):int(v) for k,v in missing.groupby('CHR').size().items()}),
        first_unassigned_examples=missing.head(20).to_json(orient='records'))
    lookup=pd.Series(np.arange(len(variants),dtype=np.int32),index=variants.SNP).to_dict()
    w=weights.loc[assigned].reset_index(drop=True);blocks=blocks[assigned]
    cols=np.array([lookup.get(x,-1) for x in w.SNP],dtype=np.int32)
    present=cols>=0
    fraction=float(present.mean()) if len(cols) else 0.
    qc.update(assigned_weighted_snps=len(w),assigned_weighted_snps_in_cache=int(present.sum()),
        assigned_cache_fraction=fraction,
        assigned_abs_beta_mass_cache_fraction=float(w.loc[present,'BETA'].abs().sum())/covered if covered else 1.)
    w=w.loc[present].reset_index(drop=True);cols=cols[present];blocks=blocks[present]
    v=variants.iloc[cols].reset_index(drop=True)
    if {"CHR","BP"} <= set(v) and (not np.array_equal(w.CHR.to_numpy(),v.CHR.to_numpy()) or not np.array_equal(w.BP.to_numpy(),v.BP.to_numpy())):
        raise ValueError(f"{pop}: dosage cache CHR/BP differs from weights")
    a1=w.A1.to_numpy(str);ref=v.REF.to_numpy(str);alt=v.ALT.to_numpy(str)
    alt_effect=a1==alt;ref_effect=a1==ref
    if not np.all(alt_effect|ref_effect):
        bad=w.loc[~(alt_effect|ref_effect),["SNP","A1"]].head()
        raise ValueError(f"{pop}: effect allele mismatch: {bad.to_dict('records')}")
    a2=w.A2.to_numpy(str)
    if np.any((a2!="") & (a2!=np.where(alt_effect,ref,alt))):
        raise ValueError(f"{pop}: non-effect allele mismatch")
    beta=w.BETA.to_numpy(np.float32)
    coef=np.where(alt_effect,beta,-beta).astype(np.float32)
    offset=np.where(ref_effect,2*beta,0).astype(np.float32)
    return dict(pop=pop,cache_col=cols,block_local=blocks,token_index=blocks+token_offset,
                channel=int(channel),coef=coef,offset=offset,weighted_snps=len(weights),
                cached_snps=int(present.sum()),cached_fraction=fraction,**qc)


def build_maps(config,layout,manifest,variants):
    paths,_=locate_weights(config,config["_trait"])
    maps=[];qc=[]
    if layout["mode"]=="population_multiview":
        for channel,(pop,path) in enumerate(zip(POPS,layout["files"])):
            bed=manifest[manifest.source_population==pop]
            weights=read_weights(paths[pop])
            offset=int(bed.token_index.min())
            m=weight_map_for_population(pop,weights,bed,variants,0,offset)
            maps.append(m);qc.append({k:v for k,v in m.items() if k not in {"cache_col","block_local","token_index","coef","offset"}})
    else:
        bed=manifest
        for channel,pop in enumerate(POPS):
            weights=read_weights(paths[pop])
            m=weight_map_for_population(pop,weights,bed,variants,channel,0)
            maps.append(m);qc.append({k:v for k,v in m.items() if k not in {"cache_col","block_local","token_index","coef","offset"}})
    for q in qc:
        q.update(trait=config['_trait'],ld_layout=layout['layout_id'])
        for field in ('boundary_assigned_snp_fraction','boundary_assigned_abs_beta_fraction'):
            if q[field]<float(config.get('min_'+field,0.95)):
                raise ValueError(f"{q['pop']}: boundary coverage sanity threshold failed: {field}={q[field]}; check map/build")
    minfrac=float(config["min_weight_snp_fraction_in_prsformer_cache"])
    bad=[q for q in qc if q["cached_fraction"]<minfrac]
    if bad:
        warnings.warn(
            f"PRSformer dosage cache coverage is below {minfrac:.2%} for PRS-CSx block reconstruction: "
            + json.dumps(bad) +
            ". Continuing with available cached SNPs.",
            RuntimeWarning,
            stacklevel=2,
        )
    return maps,pd.DataFrame(qc)


def complete_weighted_dosage(config,layout,manifest,prep_meta,cohort,variants,geno):
    """Supplement every assigned missing SNP, even when coverage exceeds 99%."""
    paths,_=locate_weights(config,config['_trait'])
    missing=[]
    for pop in POPS:
        w=read_weights(paths[pop])
        bed=manifest[manifest.source_population==pop] if layout['mode']=='population_multiview' else manifest
        missing.append(w.loc[(assign_blocks(w,bed)>=0)&~w.SNP.isin(variants.SNP)])
    missing=pd.concat(missing,ignore_index=True)
    if missing.empty:return variants,geno,None
    from block_dosage import supplement
    return supplement(config,prep_meta,cohort,variants,geno,missing,sha256_cached)


def whole_csx(config,trait):
    path=Path(config["score_root"])/trait/"1csx.scores.rds"
    if not path.is_file():raise FileNotFoundError(path)
    return identified(read_rds(path)),path


def lincheck(x,y):
    ok=np.isfinite(x)&np.isfinite(y);x=np.asarray(x)[ok];y=np.asarray(y)[ok]
    if len(x)<10 or np.std(x)==0 or np.std(y)==0:return dict(N=len(x),correlation=np.nan,slope=np.nan,intercept=np.nan,nrmse=np.nan)
    slope=float(np.cov(x,y,ddof=0)[0,1]/np.var(x));intercept=float(y.mean()-slope*x.mean())
    corr=float(np.corrcoef(x,y)[0,1]);pred=intercept+slope*x
    return dict(N=len(x),correlation=corr,slope=slope,intercept=intercept,nrmse=float(np.sqrt(np.mean((y-pred)**2))/np.std(y)))


def direct_assigned_score(geno,rows,weights,bed,variants):
    w=weights.loc[assign_blocks(weights,bed)>=0]
    lookup={s:i for i,s in enumerate(variants.SNP)}
    w=w[w.SNP.isin(lookup)]
    result=np.zeros(len(rows),dtype=np.float64)
    # Bound memory even for full UKB sample QC.
    for start in range(0,len(w),1024):
        chunk=w.iloc[start:start+1024];cols=np.array([lookup[x] for x in chunk.SNP])
        g=np.asarray(geno[np.ix_(rows,cols)],dtype=np.float64)
        effect_alt=chunk.A1.to_numpy()==variants.iloc[cols].ALT.to_numpy()
        dosage=np.where(effect_alt,g,2-g)
        result+=np.where(g>=0,dosage*chunk.BETA.to_numpy(),0).sum(axis=1)
    return result


def assigned_reconstruction(block_score,direct,config):
    x=np.asarray(block_score,float);y=np.asarray(direct,float)
    if not len(x) or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Assigned-subset reconstruction has missing/nonfinite scores')
    scale=max(float(np.std(y)),float(np.max(np.abs(y)))*1e-12,1e-12)
    nrmse=float(np.sqrt(np.mean((x-y)**2))/scale)
    if np.std(y)==0:
        corr=1. if np.allclose(x,y,rtol=1e-6,atol=1e-8) else 0.;slope=1. if corr==1. else 0.
    else:
        corr=float(np.corrcoef(x,y)[0,1]) if np.std(x)>0 else 0.
        slope=float(np.cov(y,x,ddof=0)[0,1]/np.var(y))
    q=dict(assigned_reconstruction_correlation=corr,assigned_reconstruction_slope=slope,
           assigned_reconstruction_nrmse=nrmse)
    if (corr<float(config.get('reconstruction_corr_min',.999999)) or
        not float(config.get('reconstruction_slope_min',.99999))<=slope<=float(config.get('reconstruction_slope_max',1.00001)) or
        nrmse>float(config.get('reconstruction_nrmse_max',1e-5))):
        raise ValueError(f'Assigned-subset reconstruction failed: {q}')
    return q


def cache_paths(config,layout,trait):
    # Isolate producer-bound outputs; never stamp an unbound existing array.
    base=Path(config["cache_root"])/BINDING_NAMESPACE/layout["layout_id"]/trait/"_LOCAL_ONLY"
    return base,base/"block_scores.npy",base/"tokens.tsv.gz",base/"cohort.tsv.gz",base/"BLOCKS_COMPLETE.json"


def prepare_trait(config_path,trait,block_dir=None,block_files=None,replace=False,device="auto"):
    from dosage_storage import dosage_lock
    config = read_config(config_path)
    layout = resolve_block_files(block_dir or config["ld_block_dir"],block_files or config["ld_block_files"])
    paths = cache_paths(config,layout,trait)
    if not replace and paths[-1].is_file() and not (Path(config["prsformer_cache"])/"genotypes.npy").exists():
        from abm_core.grid_inputs import locate_blocks
        arr, *_ = locate_blocks(dict(block_config=str(config_path),execution={"device":device}),
                                trait,block_files,block_dir,allow_build=False)
        arr._mmap.close()
        return layout,paths[0],json.loads(paths[-1].read_text())
    with dosage_lock(config, exclusive=True):
        return _prepare_trait(config_path,trait,block_dir,block_files,replace,device)


def _prepare_trait(config_path,trait,block_dir=None,block_files=None,replace=False,device="auto"):
    from abm_core.contracts import array_record, verify_array_binding, binding_metadata
    config=read_config(config_path);config["_trait"]=trait
    layout=resolve_block_files(block_dir or config["ld_block_dir"],block_files or config["ld_block_files"])
    manifest,channels,boundary_audit=build_manifest(layout,config)
    cache,prep_meta,cohort,variants,geno=prsformer_cache(config)
    variants,geno,dosage_supplement=complete_weighted_dosage(config,layout,manifest,prep_meta,cohort,variants,geno)
    maps,weight_qc=build_maps(config,layout,manifest,variants)
    base,score_path,token_path,cohort_path,receipt_path=cache_paths(config,layout,trait)
    base.mkdir(parents=True,exist_ok=True)
    signature_payload=dict(version=config["version"],trait=trait,layout=layout["signature"],
        prsformer_prepare=sha256(cache/"prepare.json"),variants=sha256(cache/"variants.tsv.gz"),
        dosage_supplement=dosage_supplement,
        weights={q.pop:sha256(locate_weights(config,trait)[0][q.pop]) for q in weight_qc.itertuples()},
        cohort_sha256=sha256(cache/"data.tsv.gz"),genotype_sha256=sha256_cached(cache/"genotypes.npy"),
        whole_score_sha256=sha256(Path(config["score_root"])/trait/"1csx.scores.rds"),
        source_sha256=sha256(Path(__file__)),scoring_source_sha256=sha256(Path(block_scoring.__file__)),qc_config={k:config[k] for k in config if k.startswith(("reconstruction_","min_weight_","min_boundary_"))},
        cohort_N=len(cohort),tokens=len(manifest),channels=len(channels))
    signature=hashlib.sha256(json.dumps(signature_payload,sort_keys=True).encode()).hexdigest()
    if receipt_path.is_file() and score_path.is_file() and not replace:
        old=json.loads(receipt_path.read_text())
        if receipt_matches(old,signature_payload,signature) and old.get("status")=="complete":
            arr=np.load(score_path,mmap_mode="r",allow_pickle=False)
            expected=[len(cohort),len(manifest),len(channels)]
            if (list(arr.shape)==old["shape"]==expected and arr.dtype==np.float32
                    and token_path.is_file() and cohort_path.is_file()):
                print(f"[GRID ABM] REUSE block cache {trait} {layout['layout_id']}",flush=True)
                verify_array_binding(score_path, binding_metadata(cohort),
                                     old.get("array_binding", {}), tokens=manifest, channels=channels)
                return layout,base,old
        raise ValueError(f"Existing incompatible block cache: {base}; use a new boundary file/content or --replace after review")

    if replace and base.exists() and any(base.iterdir()):
        import uuid
        archive = base.with_name("_LOCAL_ONLY_ARCHIVE_" + uuid.uuid4().hex[:12])
        base.rename(archive)
        base.mkdir(parents=True)
    # No stale completion marker may survive a failed validation.
    receipt_path.unlink(missing_ok=True)
    try:
        import torch
    except ImportError as exc:raise RuntimeError("PyTorch required") from exc
    from abm_core.support import _backend
    selected, _ = _backend(device)
    device=torch.device(selected)
    n=len(cohort);t=len(manifest);c=len(channels)
    tmp=score_path.with_name(score_path.name+".tmp")
    tmp.unlink(missing_ok=True)
    out=None
    try:
        out=np.lib.format.open_memmap(tmp,mode="w+",dtype=np.float32,shape=(n,t,c))
        batch=int(config["score_batch_size"])
        for st,en,values in block_scoring.score_batches(geno,maps,t,c,batch,device):
            out[st:en]=values
            if st==0 or en==n or (st//batch)%100==0:print(f"[GRID ABM] BLOCKS {trait} {en:,}/{n:,}",flush=True)
        out.flush()
        out._mmap.close()
        if device.type=="cuda":torch.cuda.synchronize()
        os.replace(tmp,score_path)
    finally:
        if out is not None:out._mmap.close()
        tmp.unlink(missing_ok=True)
    manifest.to_csv(token_path,sep="\t",index=False,compression="gzip")
    cohort.to_csv(cohort_path,sep="\t",index=False,compression="gzip")

    saved,saved_path=whole_csx(config,trait)
    joined=cohort[["eid"]].merge(saved,on="eid",how="left",validate="one_to_one")
    arr=np.load(score_path,mmap_mode="r")
    rng=np.random.default_rng(int(config["seed"])+17)
    rr=np.sort(rng.choice(n,min(n,int(config["reconstruction_sample_n"])),replace=False))
    reconstruction=[]
    for pi,pop in enumerate(POPS):
        if layout["mode"]=="population_multiview":
            tok=np.flatnonzero(manifest.source_population.to_numpy()==pop)
            x=np.asarray(arr[rr][:,tok,0],np.float64).sum(axis=1)
        else:
            x=np.asarray(arr[rr][:,:,pi],np.float64).sum(axis=1)
        col=next((z for z in (f"csx.{pop}",f"CSX_{pop}",f"csx_{pop}") if z in joined),None)
        if col is None:raise ValueError(f"{saved_path}: missing csx.{pop}")
        y=pd.to_numeric(joined.iloc[rr][col],errors="coerce").to_numpy(float)
        # Independent direct computation uses original alleles/betas, not token indices or coefficients.
        weights=read_weights(locate_weights(config,trait)[0][pop])
        bed=manifest[manifest.source_population==pop] if layout['mode']=='population_multiview' else manifest
        direct=direct_assigned_score(geno,rr,weights,bed,variants)
        q=assigned_reconstruction(x,direct,config)
        full=lincheck(x,y)
        q.update({f'full_csx_{k}':float(v) if np.isfinite(v) else None for k,v in full.items()})
        q.update(population=pop,comparison='covered_block_score_vs_full_csx')
        reconstruction.append(q)

    analysis=Path(config["output_root"])/layout["layout_id"]/trait/"blocks";analysis.mkdir(parents=True,exist_ok=True)
    manifest.to_csv(analysis/"block_manifest.tsv.gz",sep="\t",index=False,compression="gzip")
    weight_qc.to_csv(analysis/"block_weight_qc.tsv",sep="\t",index=False)
    pd.DataFrame(reconstruction).to_csv(analysis/"block_reconstruction.tsv",sep="\t",index=False)
    receipt=dict(status="complete",signature=signature,signature_payload=signature_payload,trait=trait,layout=layout,boundary_audit=boundary_audit,
        shape=[n,t,c],channels=channels,score_file=str(score_path),tokens=str(token_path),cohort=str(cohort_path),
        weight_qc=weight_qc.to_dict("records"),reconstruction=reconstruction,
        prsformer_training_required=False,prsformer_preparation_cache_only=dosage_supplement is None,
        dosage_supplement=dosage_supplement,
        shared_map_semantics=("one genomic token per shared block with four ancestry-specific PRS channels"
                              if layout["mode"]=="shared_universal" else
                              "one scalar token per ancestry-specific LD block; four maps retained separately"))
    identity = hashlib.sha256(json.dumps(dict(layout=layout["signature"],
        weights=signature_payload["weights"]), sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    receipt["array_binding"] = array_record(score_path,
        binding_metadata(cohort), geometry_sha256=identity,
        tokens=manifest, channels=channels, producer_version=PRODUCER_VERSION)
    receipt["producer_version"] = PRODUCER_VERSION
    receipt["executed_device"] = selected
    # Completion is published only after flush, all scientific checks and full binding.
    json_write(analysis/"BLOCKS_COMPLETE.json",receipt)
    json_write(receipt_path,receipt)
    print(json.dumps({"trait":trait,"layout_id":layout["layout_id"],"shape":[n,t,c],"status":"prepared"}),flush=True)
    return layout,base,receipt


def main(argv=None):
    a=argparse.ArgumentParser()
    a.add_argument("--config",default=str(Path(__file__).resolve().parents[1]/"config/abm.blocks.json"))
    a.add_argument("--trait",required=True,choices=("height","ldl","t2dm"))
    a.add_argument("--ld-block-dir")
    a.add_argument("--ld-block-files")
    a.add_argument("--replace",action="store_true")
    x=a.parse_args(argv);prepare_trait(x.config,x.trait,x.ld_block_dir,x.ld_block_files,x.replace)

if __name__=="__main__":main()
