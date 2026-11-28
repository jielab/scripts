"""Aggregate IBDmix carrier sets across independent batches and render the offline report."""
import base64
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, hypergeom

from asv_scan import CHROMOSOMES, POPCOUNT, save_table


def rate(a,b):
    return a/b if b else np.nan


def rank_corr(a,b):
    if len(a)<3 or np.ptp(a)==0 or np.ptp(b)==0:
        return np.nan
    return spearmanr(a,b).statistic


def comparison(win, tables):
    rng = np.random.default_rng(20261126)
    for scope,subset in [("A",win.chrom!="X"),("X",win.chrom=="X")]:
        v = win[subset & win.eligible]
        ch_sizes = v.groupby("chrom",sort=False).size().to_numpy()
        cuts = np.cumsum(ch_sizes)[:-1]
        for ref,kgref in [("Altai","alt"),("Denisova3","den")]:
            for group in ("ALL","EUR","GBR","AFR5","EAS","SAS","AMR"):
                uf = (v[ref]/v.n).to_numpy()
                kf = (v[f"{kgref}_{group}"]/v[f"n_{group}"]).to_numpy()
                pos = (uf+kf)>0
                rho_full = rank_corr(kf,uf)
                rho_pos = rank_corr(kf[pos],uf[pos])
                if group=="EUR":
                    # Exact marginal probability under uniform sampling without replacement.
                    sample_n = int(v.n_EUR.iat[0])
                    if v.n.nunique()==1 and sample_n <= int(v.n.iat[0]):
                        counts = v[ref].to_numpy()
                        vals,inv = np.unique(counts,return_inverse=True)
                        prob = hypergeom.sf(4,int(v.n.iat[0]),vals,sample_n)[inv]
                        known = v[f"{kgref}_EUR"].to_numpy()>=5
                        both = prob[known].sum()
                        either = known.sum()+prob.sum()-both
                        tables["size_matched"].append(dict(scope=scope,ref=ref,draw_n=sample_n,
                            kg_ge5=int(known.sum()),expected_ygb_ge5=prob.sum(),expected_shared=both,
                            expected_kg_seen=rate(both,known.sum()),expected_ygb_seen=rate(both,prob.sum()),
                            ratio_expected_intersection_union=rate(both,either)))
                for rule in ["5 carriers","20 carriers","1%","5%"]:
                    if "%" in rule:
                        fr = int(rule.strip("%"))/100
                        a,b = kf>=fr,uf>=fr
                    else:
                        n = int(rule.split()[0])
                        a,b = (v[f"{kgref}_{group}"]>=n).to_numpy(),(v[ref]>=n).to_numpy()
                    shared = int((a & b).sum())
                    either = int((a | b).sum())
                    j = rate(shared,either)
                    tables["concordance"].append(dict(scope=scope,ref=ref,kg_group=group,threshold=rule,
                        tiles=len(v),kg=int(a.sum()),ygb=int(b.sum()),shared=shared,union=either,
                        kg_seen_in_ygb=rate(shared,a.sum()),ygb_seen_in_kg=rate(shared,b.sum()),
                        Jaccard=j,rho_all=rho_full,rho_any=rho_pos,any_tiles=int(pos.sum()),
                        median_abs_freq_diff_any=np.median(np.abs(kf[pos]-uf[pos])) if pos.any() else np.nan))
                    if group=="EUR" and rule=="1%" and np.isfinite(j):
                        # Spatial calibration: wrap around within each chromosome's ordered eligible tiles.
                        parts = np.split(b,cuts)
                        null = []
                        for _ in range(499):
                            shifted = np.concatenate([np.roll(p,rng.integers(len(p))) for p in parts])
                            shared_null = (a & shifted).sum()
                            null.append(rate(shared_null,a.sum()+b.sum()-shared_null))
                        tables["spatial"].append(dict(scope=scope,ref=ref,kg_group=group,threshold=rule,
                            observed=j,null_median=np.median(null),null_q025=np.quantile(null,.025),
                            null_q975=np.quantile(null,.975),p_shift=(1+np.sum(np.asarray(null)>=j))/500,
                            permutations=499))


def locus_catalogue(stage,args,names,win):
    rows = []
    w = args.window_bp
    for chrom in CHROMOSOMES:
        part = win[win.chrom==chrom].reset_index(drop=True)
        for label,analogue in [("N","alt"),("D","den")]:
            bits = []
            for name in names:
                with np.load(stage/"genotypes"/name/f"chr{chrom}.npz") as data:
                    bits.append(data[label+"_carrier"])
            # At 49 batches the largest chromosome takes ~305 MB per packed track.
            # Compare every neighbor to the original sentinel, preventing single-linkage chains.
            anchor = previous = -2
            for t in np.flatnonzero(part.eligible & (part[label]>=args.min_carriers)):
                if t==previous+1:
                    inter = sum(int(POPCOUNT[m[anchor]&m[t]].sum()) for m in bits)
                    size = part.loc[anchor,label]+part.loc[t,label]-inter
                    if 5*inter>=4*size:
                        row = rows[-1]
                        row["block_end0"] = (t+1)*w
                        row["n_tiles"] += 1
                        row["min_carrier_Jaccard"] = min(row["min_carrier_Jaccard"],inter/size)
                        previous = t
                        continue
                anchor = previous = int(t)
                s = part.iloc[t]
                n = int(s.n)
                ca = int(s[label])
                unknown = int(s[label+"_uncertain"])
                rows.append(dict(id=f"ASV37-{label}-{chrom}-{t*w+1}-{(t+1)*w}",
                    lineage=label,chrom=chrom,start0=t*w,end0=(t+1)*w,block_end0=(t+1)*w,n_tiles=1,
                    min_carrier_Jaccard=1.,YGB_n=n,YGB_carriers=ca,YGB_uncertain=unknown,
                    YGB_no_match=n-ca-unknown,YGB_detected_freq=ca/n,YGB_resolved_freq=rate(ca,n-unknown),
                    YGB_LOD6=int(s[label+"_LOD6"]),batches_with_carrier=int(s[label+"_batches"]),
                    batches_ge5=int(s[label+"_batches5"]),kg_ALL=int(s[analogue+"_ALL"]),
                    kg_EUR=int(s[analogue+"_EUR"]),kg_n_EUR=int(s.n_EUR),
                    kg_SAS=int(s[analogue+"_SAS"]),kg_EAS=int(s[analogue+"_EAS"]),
                    kg_AFR5=int(s[analogue+"_AFR5"]),
                    support="EUR ≥5" if s[analogue+"_EUR"]>=5 else (
                        "1KG ALL ≥5" if s[analogue+"_ALL"]>=5 else "YGB only / <5 in 1KG")))
    fields = ("id lineage chrom start0 end0 block_end0 n_tiles min_carrier_Jaccard YGB_n YGB_carriers "
              "YGB_uncertain YGB_no_match YGB_detected_freq YGB_resolved_freq YGB_LOD6 batches_with_carrier "
              "batches_ge5 kg_ALL kg_EUR kg_n_EUR kg_SAS kg_EAS kg_AFR5 support").split()
    return pd.DataFrame(rows,columns=fields)


def combined(tables, args):
    frames = tables.pop("windows")
    kg = pd.concat([v for v in frames if v.batch.iat[0]=="1kg"],ignore_index=True)
    uk = pd.concat([v for v in frames if v.batch.iat[0]!="1kg"],ignore_index=True)
    for label in ("N","D"):
        uk[label+"_batches"] = (uk[label]>0).astype(int)
        uk[label+"_batches5"] = (uk[label]>=5).astype(int)
    idx = ["chrom","tile"]
    if uk.groupby(idx).eligible.nunique().max()!=1:
        raise ValueError("Unequal batch analysis masks")
    sparse = uk[uk.eligible].copy()
    sparse["start0"] = sparse.tile*args.window_bp
    save_table(sparse.drop(columns="eligible"),args.stage/"batch_windows.tsv.gz")
    uk = uk.drop(columns=["eligible","batch"]).groupby(idx,sort=False).sum().reset_index()
    win = kg.drop(columns=["batch","n"]).merge(uk,on=idx,validate="one_to_one",sort=False)
    win["start0"] = win.tile * args.window_bp
    win["end0"] = win.start0 + args.window_bp
    return win


def plot_figures(stage, frames, win):
    plt.rcParams.update({"svg.fonttype":"none","pdf.fonttype":42,"font.size":9,
                         "axes.spines.top":False,"axes.spines.right":False,"axes.titlesize":11})
    teal, rust, blue = "#00796b","#ac4831","#345886"
    out = stage/"figures"
    out.mkdir()
    images = []
    def finish(fig,name):
        fig.savefig(out/(name+".pdf"),bbox_inches="tight")
        fig.savefig(out/(name+".svg"),bbox_inches="tight")
        fig.savefig(out/(name+".png"),dpi=220,bbox_inches="tight")
        plt.close(fig)
        images.append(base64.b64encode((out/(name+".png")).read_bytes()).decode())
    fig,axes = plt.subplots(2,2,figsize=(8.2,6.7),constrained_layout=True)
    for i, (ref,key) in enumerate([("Altai","alt"),("Denisova3","den")]):
        for j, (scope,select) in enumerate([("Autosomes",win.chrom!="X"),("Male X",win.chrom=="X")]):
            ax=axes[j,i]
            v=win[win.eligible & select]
            a=v[key+"_EUR"]/v.n_EUR*100
            b=v[ref]/v.n*100
            h=ax.hexbin(np.log10(1+a),np.log10(1+b),gridsize=45,mincnt=1,cmap="viridis",bins="log")
            ax.plot([0,2.01],[0,2.01],color="#888",ls="--",lw=.8)
            ticks=np.log10(1+np.array([0,1,5,20,100]))
            ax.set(xticks=ticks,yticks=ticks,xticklabels=["0","1","5","20","100"],
                   yticklabels=["0","1","5","20","100"],title=f"{scope}  ·  {ref}")
            if j==1: ax.set_xlabel("1KG EUR carriers (%)")
            if i==0: ax.set_ylabel("YGB carriers (%)")
            fig.colorbar(h,ax=ax,label="Tiles",shrink=.78)
    finish(fig,"F1_concordance")
    fig,axes=plt.subplots(2,1,figsize=(9.1,5.3),sharex=True,constrained_layout=True)
    x=win[win.chrom=="X"]
    for ax,ref,k in zip(axes,["Altai","Denisova3"],["alt","den"]):
        for y,label,c in [(100*x[ref]/x.n,"YGB",rust),(100*x[k+"_EUR"]/x.n_EUR,"1KG EUR",teal)]:
            ax.plot((x.start0+x.end0)/2e6,y.where(x.eligible),lw=.9,alpha=.9,color=c,label=label)
        ax.axvspan(58.632012,61.632012,color="0.82",label="Centromere")
        ax.set(ylabel="Carriers (%)",title=f"Male non-PAR X  ·  {ref}")
        ax.legend(frameon=False,ncol=3,loc="upper right")
    axes[-1].set_xlabel("GRCh37 chrX position (Mb)")
    finish(fig,"F2_chrX")
    fig,ax=plt.subplots(1,2,figsize=(9,3.8),constrained_layout=True)
    b=frames["benchmark"].iloc[2:]
    ix=np.arange(len(b))
    ax[0].bar(ix-.18,b.Chen_2020,.36,color="#a3b2bb",label="Cell 2020")
    ax[0].bar(ix+.18,b.observed,.36,color=teal,label="Current 1KG")
    ax[0].set(xticks=ix,xticklabels=[s.split()[0] for s in b.measure],ylabel="Autosomal Altai Mb / person",
              title="Phase 3 benchmark")
    ax[0].legend(frameon=False)
    cr=frames["control"]
    picks=[("1kg",g,"Altai.2013") for g in ["EUR","SAS","EAS","AFR5"]]
    labels=["EUR","SAS","EAS","AFR5"]
    for ref in ("Altai","Vindija"):
        s=cr[(cr.batch!="1kg")&(cr.group=="ALL")&(cr.ref==ref)]
        picks.append(("YGB","ALL",ref));labels.append("YGB\n"+ref)
    values=[]
    for batch,g,ref in picks:
        v=cr[(cr.group==g)&(cr.ref==ref)&((cr.batch!="1kg") if batch=="YGB" else (cr.batch==batch))]
        values.append(100*v.core80.sum()/v.n.sum())
    ax[1].bar(range(len(values)),values,color=[teal]*4+[rust]*2)
    ax[1].set(xticks=range(len(labels)),xticklabels=labels,ylabel="Carriers covering ≥80% of core (%)",
              title="3p21.31  ·  49.4 kb published core")
    finish(fig,"F3_controls")
    fig,ax=plt.subplots(1,2,figsize=(9,3.8),constrained_layout=True)
    cat=frames["loci"]
    for scope,style in [("Autosomes","-"),("Male X","--")]:
        for label,c in [("N",teal),("D",rust)]:
            s=cat[(cat.lineage==label)&((cat.chrom!="X") if scope=="Autosomes" else (cat.chrom=="X"))]
            if len(s):
                x=np.sort(s.YGB_carriers)
                ax[0].step(x,np.arange(len(s),0,-1),where="pre",lw=1.5,color=c,ls=style,label=f"{scope} · {label}")
    ax[0].set(xscale="log",yscale="log",xlabel="Minimum YGB carriers",ylabel="Sentinel loci",title="Carrier spectrum")
    ax[0].legend(frameon=False,fontsize=8)
    ratio=frames["x_ratio"]
    for j,(batch,ref) in enumerate([("1kg","Altai.2013"),("YGB","Altai"),("YGB","Vindija"),("YGB","Denisova3")]):
        s=ratio[(ratio.batch==batch)&(ratio.ref==ref)]
        for k,(f,c) in enumerate([("export","#a3b2bb"),("gap_free",blue)]):
            ax[1].bar(j+k*.35,s.loc[s["filter"]==f,"ratio_X_A"].iloc[0],.32,color=c,
                      label=(f if j==0 else None))
    ax[1].axhline(1,color="0.4",ls="--",lw=.8)
    ax[1].set(xticks=np.arange(4)+.175,xticklabels=["1KG\nAltai","YGB\nAltai","YGB\nVindija","YGB\nDenisova3"],
              ylabel="X / autosomal physical footprint",title="Same X-tested men")
    ax[1].legend(frameon=False)
    finish(fig,"F4_catalogue")
    amb=frames["reference_support"]
    fig,ax=plt.subplots(figsize=(9,3.9),constrained_layout=True)
    bottom=np.zeros(len(amb))
    for key,title,col in [("N","N: ≥2 N only",teal),("D","D: both D only",blue),
                          ("mixed","Both lineages",rust),("single_ref","One reference only","#b1beb8")]:
        heights=100*amb[key].to_numpy()/amb.any_ref.replace(0,np.nan).to_numpy()
        ax.bar(range(len(amb)),heights,bottom=bottom,color=col,width=.82,label=title)
        bottom+=heights
    ax.set(xticks=range(len(amb)),xticklabels=amb.chrom.tolist(),ylim=(0,100),xlabel="Chromosome",
           ylabel="Positive person–tiles (%)",title="Reference support within a person and a tile")
    ax.legend(frameon=False,ncol=4,loc="upper center",bbox_to_anchor=(.5,-.21),fontsize=8)
    finish(fig,"F5_reference_support")
    return images


def finish(stage,args,manifests,tables):
    args.stage=stage
    win=combined(tables,args)
    if not win.loc[win.eligible,"chrom"].eq("X").any() or not win.loc[win.eligible,"chrom"].ne("X").any():
        raise ValueError("No comparable autosomal or central non-PAR X tiles; check GRCh37 callable BED coverage.")
    comparison(win,tables)
    names=[m["name"] for m in manifests[1:]]
    tables["loci"]=locus_catalogue(stage,args,names,win)
    frames={k:pd.DataFrame(v) for k,v in tables.items() if k not in ("batches","x_ratio")}
    loc=frames["loci"]
    loc["priority_EUR_LOD6"]=(loc.kg_EUR>=5) & (loc.YGB_LOD6>=5)
    save_table(loc[loc.priority_EUR_LOD6],stage/"priority_EUR_LOD6.tsv.gz")
    counts=win[win.eligible].groupby("chrom",sort=False)[["any_ref","mixed","N","D"]].sum().reset_index()
    counts["single_ref"]=counts.any_ref-counts.mixed-counts.N-counts.D
    frames["reference_support"]=counts
    b=frames["burden"]
    xp=b[b.group.eq("Xtested")].copy()
    xp["batch"]=np.where(xp.batch=="1kg","1kg","YGB")
    idx=["batch","ref","filter"]
    xa=xp.groupby(idx+["scope"]).agg(total_bp=("total_bp","sum"),n=("n","sum"),span=("span_bp","first"))
    xa["frac"]=xa.total_bp/(xa.n*xa.span)
    xa=xa.frac.unstack("scope").reset_index()
    xa["ratio_X_A"]=xa.X/xa.A
    frames["x_ratio"]=xa
    meta={"build":"GRCh37","coordinate_system":"0-based half-open","window_bp":args.window_bp,
          "min_carriers":args.min_carriers,"gap_min_bp":100000,
          "common_callable":str(args.callable_bed or ""),"minimum_callable_tile_fraction":0.8,"source":args.source,
          "callable_sha256":hashlib.sha256(args.callable_bed.read_bytes()).hexdigest() if args.callable_bed else None,
          "cohorts":manifests}
    meta["raw_interpretation"]={"coordinate_conversion":"native start/end minus 1", "african_denisovan_background":"disabled", "thresholds":"per-chromosome run.meta.tsv lod_cut and len_cut"}
    (stage/"release.json").write_text(json.dumps(meta,indent=2)+"\n")
    for key,v in frames.items():
        save_table(v,stage/(key+".tsv"+(".gz" if key=="loci" else "")))
    save_table(win,stage/"windows.tsv.gz")
    fig=plot_figures(stage,frames,win)
    from asv_html import write_html
    write_html(stage,meta,frames,win,fig)
