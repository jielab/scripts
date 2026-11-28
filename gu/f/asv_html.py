"""Render a self-contained, searchable report; no individual identifiers in HTML."""
import json
from html import escape
from pathlib import Path

import numpy as np
import pandas as pd


def write_html(stage,meta,frames,win,images):
    source=Path(__file__).with_name("asv_report.html").read_text()
    def records(df):
        return json.loads(df.to_json(orient="values",double_precision=6))
    def table(df, digits=3):
        return df.to_html(index=False, border=0, na_rep="—",float_format=lambda x:f"{x:,.{digits}f}")
    cat=frames["loci"]
    kg=meta["cohorts"][0]
    uk=meta["cohorts"][1:]
    def num(x):return format(int(x),",")
    summary=f'1KG {num(kg["N"])} 人 · YGB {num(sum(c["N"] for c in uk))} 人 · {len(uk)} 个 YGB 批次 · GRCh37'
    audit=frames["audit"].groupby("batch",sort=False).agg(rows=("n_rows","sum"),included_rows=("non_sentinel_rows","sum"),gap_rows=("gap_touching_rows","sum")).reset_index()
    counts=cat.assign(scope=np.where(cat.chrom.eq("X"),"X","Autosomes"),
                      kg5=cat.kg_ALL.ge(5),eur5=cat.kg_EUR.ge(5),strong=cat.YGB_LOD6.ge(5),priority=cat.priority_EUR_LOD6).groupby(["scope","lineage"],sort=False).agg(
                          sentinels=("id","size"),tiles=("n_tiles","sum"),kg5=("kg5","sum"),eur5=("eur5","sum"),LOD6_ge5=("strong","sum"),EUR5_and_LOD6=("priority","sum")).reset_index()
    comp=frames["concordance"]
    rows=comp[(comp.kg_group=="EUR") & (comp.threshold=="1%")]
    a=rows[(rows.scope=="A")&(rows.ref=="Altai")].iloc[0]
    xr=frames["x_ratio"]
    male=frames["male_tiles"].copy()
    male["cohort"]=np.where(male.batch.eq("1kg"),"1kg","YGB")
    s=male.groupby(["cohort","ref","scope"])[["carrier_tiles","tested_person_tiles"]].sum()
    s["carrier_fraction"]=s.carrier_tiles/s.tested_person_tiles
    mz=s.carrier_fraction.unstack("scope").reset_index()
    mz["ratio_X_A"]=mz.X/mz.A
    ctrl=frames["control"].copy()
    ctrl["cohort"]=np.where(ctrl.batch.eq("1kg"),"1kg","YGB")
    cs=ctrl.groupby(["cohort","group","ref"],as_index=False)[["n","any_bp","core80"]].sum()
    cs["core80_pct"]=100*cs.core80/cs.n
    cs=cs[(cs["group"].isin(["EUR","SAS","ALL"])) & cs.ref.isin(["Altai.2013","Altai","Vindija","Chagyr"])]
    fields=["id","lineage","chrom","start0","end0","block_end0","n_tiles","YGB_carriers","YGB_n",
            "YGB_uncertain","YGB_LOD6","kg_ALL","kg_EUR","kg_n_EUR","batches_ge5"]
    data={"w":meta["window_bp"],"m":meta["min_carriers"],"l":records(cat[fields]),
          "c":list(comp.columns),"comp":records(comp),"tracks":{}}
    for ch,v in win.groupby("chrom",sort=False):
        data["tracks"][ch]=records(pd.DataFrame({"e":v.eligible.astype(int),"N":v.Altai/v.n,"n":v.alt_EUR/v.n_EUR,
                                                "D":v.Denisova3/v.n,"d":v.den_EUR/v.n_EUR}))
    substitutions={"SUMMARY":escape(summary),"W":num(meta["window_bp"]),"MIN":str(meta["min_carriers"]),
           "NLOCI":num(len(cat)),"NX":num(sum(c["N_X"] for c in uk)),"RECOVERY":f"{a.kg_seen_in_ygb:.1%}",
           "RHO":f"{a.rho_any:.3f}","NRECORDS":num(audit.rows.sum()),"COUNTS":table(counts),
           "BENCH":table(frames["benchmark"],4),"SIZE":table(frames["size_matched"]),"NULL":table(frames["spatial"]),
           "X":table(xr),"XWIN":table(mz),"CONTROL":table(cs),"AUDIT":table(audit,0),
           "ROSTERS":table(pd.DataFrame(meta["cohorts"])),
           "PROV":escape(meta["source"] or "Unspecified"),"MASK":escape(meta["common_callable"] or "未提供共同 callable BED。分析区域为物理坐标代理。"),
           "DATA":json.dumps(data,ensure_ascii=False,separators=(",",":"),allow_nan=False).replace("<","\\u003c"),
           **{f"FIG{i+1}":im for i,im in enumerate(images)}}
    import re
    source=re.sub(r"@@([A-Z0-9]+)@@",lambda m:substitutions[m[1]],source)
    (stage/"asv.html").write_text(source,encoding="utf-8")
