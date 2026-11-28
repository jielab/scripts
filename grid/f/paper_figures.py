"""Build GRID manuscript exhibits from released aggregate TSVs. No prediction fitting."""
import argparse
import fcntl
import html
import json
import os
import shutil
import tempfile
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator, ScalarFormatter
import numpy as np
import pandas as pd


TRAITS = ("height", "ldl", "t2dm")
NAMES = {"height": "Height", "ldl": "LDL cholesterol", "t2dm": "Type 2 diabetes",
         "ldl.high": "High LDL cholesterol"}
CSX, DISCO, STACK = "csx.population", "disco", "block.stack"
LABEL = {CSX: "PRS-CSx", DISCO: "DiscoDivas", STACK: "Four-channel stack", "covariate": "Covariates",
         "abm": "GRID ABM", "direct.linear": "Global ridge", "direct.global_mlp": "Global MLP",
         "direct.blocknet": "Direct BlockNet", "direct.regional_ridge": "Regional ridge",
         "block.glmnet.disco.offset": "Global offset", "block.glmnet.disco": "Disco regional elastic net",
         "blocknet.disco.global": "Global Disco network", "blocknet.disco": "Disco regional BlockNet"}
POLICY = {"all": "Full test set", "raw_global": "Global count", "regional_count": "Regional count",
          "latent_count": "Latent count", "risk_only": "Risk only", "risk_count": "Risk + count",
          "risk_local": "Risk + local features"}
INK, MUTED, LIGHT, TEAL, BLUE, GOLD, PURPLE = (
    "#182e3a", "#506775", "#d9e2e6", "#007867", "#326798", "#a6541d", "#765486")
LABEL.update({"universal.auto":"Universal GRID","universal.fixed":"Fixed-φ GRID","universal.global.auto":"Auto-φ global sum","universal.global.fixed":"Fixed-φ global sum"})
STYLE = {CSX: (BLUE, "s"), DISCO: (GOLD, "^"), STACK: (PURPLE, "o"), "universal.auto":(TEAL,"D")}
POPS = ("AFR", "EAS", "EUR", "SAS")
WIDTH = 7.2
RC = {"font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
      "axes.titleweight": "bold", "axes.labelcolor": INK, "text.color": INK,
      "xtick.color": MUTED, "ytick.color": INK, "axes.edgecolor": MUTED,
      "axes.spines.right": False, "axes.spines.top": False, "axes.linewidth": .55,
      "xtick.labelsize": 7.5, "ytick.labelsize": 8, "pdf.fonttype": 42,
      "ps.fonttype": 42, "svg.fonttype": "none", "savefig.facecolor": "white"}
SUPP = {
    "universal_atlas": ("S8", "TRAIN regional score atlas for the auto and fixed META weights"),
    "universal_coefficients": ("S9", "Tuned regional correction coefficients in original block-score units"),
    "discovery": ("S1a", "Discovery GWAS sample sizes and overlap metadata"),
    "ld_maps": ("S1b", "LD block maps and SNP assignment"),
    "score_coverage": ("S1c", "Score variant and chromosome coverage"),
    "cohort": ("S1d", "Cohort roles and genetic ancestry"),
    "test_demographics": ("S1e", "Test-roster baseline age and source-coded sex distribution"),
    "stack_weights": ("S2e", "Frozen calibration weights in the convex stack"),
    "availability": ("S2a", "Available and omitted methods"),
    "model_settings": ("S2b", "Model and validation settings"),
    "anchor_training": ("S2c", "DiscoDivas anchor training"),
    "csx_weights": ("S2d", "PRS-CSx and single-ancestry score coefficients"),
    "performance": ("S3", "Absolute predictive performance by endpoint and ancestry"),
    "contrasts": ("S4", "Paired method comparisons and family-bootstrap intervals"),
    "screening": ("S5", "ABM screening, coverage and risk-standardized loss"),
    "ldl_twins": ("S6", "LDL twin policies and paired error contrasts"),
}


def subset(frame, **criteria):
    mask = pd.Series(True, index=frame.index)
    for key, value in criteria.items():
        mask &= frame[key].eq(value)
    return frame.loc[mask]


def fmt(value, digits=0):
    return f"{value:,.{digits}f}" if pd.notna(value) and np.isfinite(value) else "—"


def axis(ax, xlabel="", zero=False):
    ax.set_axisbelow(True)
    ax.grid(axis="x", color=LIGHT, lw=.55)
    ax.tick_params(length=2.5, width=.5)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel(xlabel)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=2))
    f = ScalarFormatter(useOffset=False)
    f.set_scientific(False)
    ax.xaxis.set_major_formatter(f)
    if zero:
        ax.axvline(0, color=MUTED, lw=.9, ls=(0, (3, 2)), zorder=1)


def panel(ax, letter, title):
    ax.set_title(f"{letter}   {title}", loc="left", pad=10)


def blank(ax, label="Not available in this run"):
    ax.set_axis_off()
    ax.text(.5, .5, label, ha="center", va="center", color=MUTED, transform=ax.transAxes)


def interval(ax, x, low, high, y, color=TEAL, marker="o"):
    if not np.isfinite(x):
        return
    if np.isfinite([low, high]).all():
        ax.plot([low, high], [y, y], color=color, lw=1.6, solid_capstyle="round", zorder=2)
    ax.scatter(x, y, s=22, marker=marker, facecolor=color, edgecolor="white", lw=.35, zorder=3)


def ordered_rows(frame, column, order):
    return frame.set_index(column).reindex(order)


def finite_rows(frame, column):
    return frame[np.isfinite(frame[column])]


def make_table1(perf, cohort, demographics):
    # One aggregate row per stratum is enough; outcomes never depend on the model.
    test = perf[perf.method == CSX].set_index(["endpoint", "ancestry"])
    pops = list(dict.fromkeys(["ALL", *[p for p in POPS if p in cohort.ancestry.values],
                               *sorted(set(cohort.ancestry) - set(POPS) - {"ALL"})]))
    roster = cohort.pivot(index="role", columns="ancestry", values="n").fillna(0)
    roster["ALL"] = roster.sum(axis=1)
    rows = []
    def add(section, name, values):
        rows.append([section, name, *values])
    for role, label in [("train", "Training"), ("tune", "Tuning"), ("calibration", "Calibration"), ("test", "Test")]:
        if role in roster.index:
            add("Cohort allocation", label + ", n", [fmt(roster.at[role,p]) for p in pops])
    add("Cohort allocation", "Total roster, n", [fmt(roster[p].sum()) for p in pops])
    if not demographics.empty:
        demo=demographics.set_index("ancestry").reindex(pops)
        if (demo.n_roster.notna() & (demo.n_roster != roster.loc["test",pops])).any():
            raise ValueError("Baseline demographic counts disagree with the test roster.")
        if "age_mean" in demo:
            add("Baseline in TEST", "Age, mean (SD), years",
                [f"{fmt(r.age_mean,1)} ({fmt(r.age_sd,1)})" if pd.notna(r.age_mean) else "—" for r in demo.itertuples()])
        codes=[c for c in demo if c.startswith("sex_code_")]
        if codes:
            # Report the stored coding. Custom cohort inputs can reverse the sex coding.
            for code in sorted(codes):
                counts=demo[code].fillna(0)
                add("Baseline in TEST",f"Sex code {code[9:]}, n (% known)",
                    [f"{fmt(counts[p])} ({fmt(100*counts[p]/demo.at[p,'sex_n'],1)}%)"
                     if demo.at[p,'sex_n']>0 else "—" for p in pops])
    for trait in TRAITS:
        if trait not in test.index.get_level_values(0):
            continue
        d = test.loc[trait].reindex(pops)
        add(NAMES[trait] + " in TEST", "Evaluable, n (% of test roster)",
            [f"{fmt(d.at[p,'n'])} ({fmt(100*d.at[p,'n']/roster.at['test',p],1)}%)"
             if pd.notna(d.at[p,'n']) and roster.at['test',p] else "—" for p in pops])
        if trait == "t2dm":
            add(NAMES[trait] + " in TEST", "Prevalent cases, n (%)",
                [f"{fmt(r.n_cases)} ({fmt(100*r.event_fraction,2)}%)" if pd.notna(r.n_cases) else "—"
                 for r in d.itertuples()])
        else:
            unit = "cm" if trait == "height" else "mmol/L"
            add(NAMES[trait] + " in TEST", f"Mean (SD), {unit}",
                [f"{fmt(r.outcome_mean,2)} ({fmt(np.sqrt(r.outcome_variance*r.n/(r.n-1)),2)})"
                 if r.n > 1 else "—" for r in d.itertuples()])
    return pd.DataFrame(rows, columns=["Section", "Characteristic", *pops])


def figure_table(frame):
    fig, ax = plt.subplots(figsize=(9.35, 5.5))
    fig.subplots_adjust(left=.026, right=.984, bottom=.05, top=.96)
    ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis("off")
    ax.text(0, .975, "Table 1", fontsize=13, weight="bold", va="top", color=TEAL)
    ax.text(.105, .975, "Study population and test-set phenotypes", fontsize=13, weight="bold", va="top")
    cols = frame.columns[2:]
    right = np.linspace(.405, .993, len(cols))
    y = .865
    for val, xx in zip(cols, right):
        ax.text(xx, y, "Overall" if val == "ALL" else val, ha="right", weight="bold")
    ax.plot([0, 1], [y-.017]*2, color=INK, lw=.8)
    prev = None
    step = min(.058, .70/(len(frame)+frame.Section.nunique()+1))
    for row in frame.itertuples(index=False, name=None):
        group, label, *values = row
        if group != prev:
            y -= step*1.2
            ax.text(0, y, group, va="center", color=TEAL, weight="bold")
            prev = group
        y -= step
        ax.text(.014, y, label, va="center")
        for val, xx in zip(values, right):
            ax.text(xx, y, val, ha="right", va="center", fontsize=7.7)
    ax.plot([0,1], [y-step*.65]*2, color=INK, lw=.7)
    notes = ("Percentages use the test roster of the corresponding ancestry, except for case prevalence (evaluable participants) and sex, if recorded (% of known). "
             "SD is the sample standard deviation. OTH denotes other/unassigned ancestry. Participants can contribute to multiple traits.")
    ax.text(0, y-step, textwrap.fill(notes, 148), va="top", fontsize=7.4, linespacing=1.45, color=MUTED)
    return fig


def design_figure(d):
    fig, ax = plt.subplots(figsize=(WIDTH,5.05))
    fig.subplots_adjust(left=.042,right=.98,top=.975,bottom=.025)
    ax.set(xlim=(0,1),ylim=(0,1)); ax.axis("off")
    def txt(x,y,s,**style):
        ax.text(x,y,s,**style)
    def link(x1,y1,x2,y2):
        ax.annotate("",xy=(x2,y2),xytext=(x1,y1),
                    arrowprops=dict(arrowstyle="-|>",color=MUTED,lw=.9,shrinkA=0,shrinkB=1))
    def tag(x,y,s):
        txt(x,y,s,fontsize=8.3,color=TEAL,weight="bold")
    maps=d["ld_maps"]
    shared=(not maps.empty and maps.blocks.nunique()==1 and maps.layout.eq("consensus").all())
    n=f"{int(maps.blocks.iloc[0]):,}" if shared else "B"
    txt(0,.95,"A   Ancestry interpolation within each LD block",weight="bold",fontsize=10.1)
    tag(.01,.866,"1   Score four channels")
    tag(.385,.866,"2   Interpolate slopes")
    tag(.727,.866,"3   Weight each block")
    txt(.01,.805,"Regional PRS-CSx scores",fontsize=8)
    # Diagram colours identify the discovery channel, not an observed effect size.
    for i,(name,color) in enumerate(zip(("AFR","EAS","EUR","SAS"),(BLUE,GOLD,TEAL,PURPLE))):
        y=.742-i*.049
        txt(.014,y,name,va="center",fontsize=7.5,color=color,weight="bold")
        for j in range(5):
            ax.add_patch(Rectangle((.076+j*.052,y-.017),.047,.033,facecolor=color,edgecolor="none",alpha=.65))
    for j,name in enumerate(("1","2","3","…","B")):
        txt(.099+j*.052,.534,name,ha="center",fontsize=7.2)
    txt(.20,.479,f"{n} blocks" if shared else "LD block index",ha="center",fontsize=8)
    txt(.39,.805,"Target projected into 1KG PC space",fontsize=7.7)
    anchors=[(.43,.61,BLUE),(.48,.733,GOLD),(.626,.69,TEAL),(.579,.56,PURPLE)]
    target=(.544,.65)
    for x,y,color in anchors:
        ax.plot([x,target[0]],[y,target[1]],color=LIGHT,lw=1,zorder=0)
        ax.scatter(x,y,s=60,c=color,marker="o")
    ax.scatter(*target,s=46,c=INK,marker="*",zorder=3)
    txt(.39,.49,"●  Anchor      ★  Target person",fontsize=7.2)
    link(.665,.657,.706,.657)
    txt(.73,.73,r"$z_{ib}=\sum_{a} w_{ia}\,s_{iba}$",fontsize=13)
    txt(.73,.637,"One weighted score per\nindividual i and block b",fontsize=8,linespacing=1.45)
    txt(.73,.506,"Weights vary by person.",fontsize=7.4,color=MUTED)
    ax.plot([0,1],[.431,.431],color=LIGHT,lw=.8)
    txt(0,.373,"B   Two uses of regional information",weight="bold",fontsize=10.1)
    tag(.01,.295,"Predict the phenotype")
    txt(.01,.245,"DiscoDivas + block elastic net + unblended\nDisco regional elastic net and BlockNet",fontsize=7.7,va="top",linespacing=1.45)
    link(.415,.212,.459,.212)
    txt(.47,.253,"Fit convex stack weights\non calibration families",va="top",fontsize=7.8,linespacing=1.45)
    ax.plot([.704,.704],[.146,.316],lw=.8,color=LIGHT)
    tag(.731,.295,"Screen for lower error")
    txt(.731,.245,"Count TRAIN neighbours.\nFreeze selection rules.\nEvaluate error in TEST.",va="top",fontsize=7.5,linespacing=1.35)
    ax.plot([0,1],[.106,.106],color=LIGHT,lw=.8)
    co=d["cohort"].groupby("role").n.sum()
    txt(.01,.044,"Family splits",weight="bold",fontsize=8)
    s="    ·    ".join(f"{name}  {fmt(co.get(role))}" for role,name in
                       [("train","Train"),("tune","Tune"),("calibration","Calibrate"),("test","Test")])
    txt(.18,.044,s,fontsize=8)
    return fig


def atlas_figure(d, traits):
    atlas=d["universal_atlas"]
    if atlas.empty:
        return None
    atlas=atlas[atlas["mode"].eq("auto")].copy()
    measured=[t for t in traits if t in atlas.trait.values]
    if not measured: return None
    coords=["chrom38","start38","end38","block_id"]
    frames=[]
    for trait in measured:
        fr=atlas[atlas.trait.eq(trait)].sort_values(coords).set_index(coords)
        if fr.index.has_duplicates or (frames and not fr.index.equals(frames[0].index)):
            raise ValueError("Fig. 1 requires the same ordered universal BED for every trait")
        frames.append(fr)
    n=len(frames[0])
    mat=np.stack([f.variance_share.to_numpy() for f in frames])
    if not np.isfinite(mat).all() or np.any(mat<0) or not np.allclose(mat.sum(1),1):
        raise ValueError("Invalid within-trait block variance shares")
    fig=plt.figure(figsize=(WIDTH,6.35))
    gs=fig.add_gridspec(3,2,left=.12,right=.97,bottom=.105,top=.92,
                         width_ratios=[1.12,1],height_ratios=[.58,1,1.1],hspace=.82,wspace=.46)
    design=fig.add_subplot(gs[0,:]);design.axis("off")
    panel(design,"A","From one posterior weight per SNP to a regional predictor")
    stages=[(.00,"SNP weights","One META effect per SNP\nwithin each trait"),
            (.36,"Regional scores","Sum weighted dosages within\neach shared LD interval"),
            (.71,"Regional learning","Fit block coefficients on TRAIN;\nchoose shrinkage on TUNE")]
    for i,(x,title,detail) in enumerate(stages):
        design.add_patch(Rectangle((x,.74),.25,.035,transform=design.transAxes,color=[BLUE,TEAL,GOLD][i]))
        design.text(x,.66,title,transform=design.transAxes,fontsize=8.5,weight="bold",va="top")
        design.text(x,.35,detail,transform=design.transAxes,fontsize=7.1,va="top",linespacing=1.35)
        if i<2:
            design.annotate("",xy=(x+.335,.52),xytext=(x+.276,.52),xycoords="axes fraction",
                            arrowprops={"arrowstyle":"->","color":MUTED,"lw":1})
    a=fig.add_subplot(gs[1,:]);b=fig.add_subplot(gs[2,0]);c=fig.add_subplot(gs[2,1])
    panel(a,"B",f"One shared BED across traits   ·   {n:,} ordered blocks")
    heat=a.imshow(np.log10(1+n*mat),aspect="auto",interpolation="none",
                  cmap="magma",vmin=0,vmax=2,extent=(-.5,n-.5,len(measured)-.5,-.5))
    a.set_yticks(range(len(measured)),[NAMES[t].replace(" cholesterol","").replace("Type 2 diabetes","T2DM") for t in measured])
    ch=np.asarray(frames[0].index.get_level_values("chrom38"))
    ticks=[];labs=[]
    for chrom in pd.unique(ch):
        ix=np.flatnonzero(ch==chrom);ticks.append(ix.mean());labs.append(str(chrom))
        if ix[0]: a.axvline(ix[0]-.5,color="white",alpha=.6,lw=.4)
    # Fewer tick labels on short terminal chromosomes; every boundary remains visible.
    a.set_xticks(ticks,[s if int(s)<=14 or int(s)%2==0 else "" for s in labs],fontsize=6.7)
    a.set_xlabel("Autosome · one column per LD block");a.tick_params(length=0)
    color=fig.colorbar(heat,ax=a,location="right",pad=.025,shrink=.95,ticks=[0,1,2],extend="max")
    color.ax.set_yticklabels(["0","1","2"]);color.set_label("log₁₀(1 + relative variance)",fontsize=7)
    panel(b,"C","Concentration of score variance")
    b.plot([0,100],[0,100],color=MUTED,ls=":",lw=.7)
    pal={"height":BLUE,"ldl":TEAL,"t2dm":GOLD}
    for t,v in zip(measured,mat):
        b.plot(np.linspace(0,100,n+1),100*np.r_[0,np.sort(v)[::-1].cumsum()],color=pal[t],label=t.upper() if t!="height" else "Height")
    b.set(xlim=(0,100),ylim=(0,101),xlabel="Blocks ranked by variance (%)",ylabel="Cumulative variance share (%)")
    b.legend(frameon=False,fontsize=6.8,loc="lower right")
    panel(c,"D","Effective number of blocks")
    for i,(t,v) in enumerate(zip(measured,mat)):
        count=1/np.sum(v*v)
        c.barh(i,count,height=.48,color=pal[t])
        c.annotate(f"{count:,.0f}",(count,i),xytext=(3,0),textcoords="offset points",va="center",fontsize=8)
    c.set_yticks(range(len(measured)),[t.upper() if t!="height" else "Height" for t in measured]);c.invert_yaxis();c.margins(x=.22)
    axis(c,"Inverse sum of squared variance shares")
    return fig


def benchmark_figure(d, traits):
    pts, pairs=d["performance"], d["contrasts"]
    methods=["universal.auto","universal.global.auto",STACK,CSX,DISCO]
    names=["Universal GRID","Matched global sum","Four-channel stack","PRS-CSx","DiscoDivas"]
    if not pts.method.eq("universal.auto").any():
        return None
    colors=[TEAL,MUTED,PURPLE,BLUE,GOLD,"#62475f"]
    markers=["D","o","s","s","^","v"]
    fig,axes=plt.subplots(2,len(traits),squeeze=False,figsize=(WIDTH,4.6))
    fig.subplots_adjust(left=.23,right=.96,top=.90,bottom=.12,wspace=.34,hspace=.67)
    for j,t in enumerate(traits):
        a,b=axes[:,j]
        met,delta,label=("auroc","AUROC","AUROC") if t=="t2dm" else ("r2","predictive_R2","Predictive R²")
        panel(a,chr(65+j),NAMES[t])
        s=ordered_rows(subset(pts,endpoint=t,ancestry="ALL"),"method",methods)
        for i,(m,r) in enumerate(s.iterrows()):
            x=r[met]
            if pd.notna(x):
                a.scatter(x,i,color=colors[i],marker=markers[i],s=23)
                a.annotate(f"{x:.4f}",(x,i),xytext=(0,4),textcoords="offset points",ha="center",fontsize=6.5)
            else:
                a.text(.50,i,"n/a",transform=a.get_yaxis_transform(),ha="center",color=MUTED,fontsize=7)
        a.set_yticks(range(len(methods)),names if j==0 else [""]*len(methods));a.set_ylim(len(methods)-.6,-.8)
        axis(a,label)
        finite=s[met].dropna()
        if not finite.empty:
            if t=="t2dm":
                a.set_xlim(min(.5,float(finite.min())-.015),1.)
            else:
                a.set_xlim(min(-.015,float(finite.min())-.035),max(.25,float(finite.max())+.07))
        panel(b,chr(65+len(traits)+j),"Paired difference")
        q=ordered_rows(subset(pairs,endpoint=t,ancestry="ALL",method="universal.auto",metric=delta),"reference",methods[1:])
        for i,r in enumerate(q.itertuples()):
            interval(b,r.difference,r.ci_low,r.ci_high,i,colors[i+1],markers[i+1])
        b.set_yticks(range(len(q)),names[1:] if j==0 else [""]*len(q));b.set_ylim(len(q)-.55,-.55)
        axis(b,"ΔAUROC" if t=="t2dm" else "ΔR²",True)
    return fig


def ldl_figure(d):
    p=subset(d["performance"],endpoint="ldl",ancestry="ALL")
    methods=["universal.global.auto","universal.auto","universal.global.fixed","universal.fixed","direct.linear","direct.regional_ridge","direct.global_mlp","direct.blocknet",
             "block.glmnet.disco.offset","block.glmnet.disco","blocknet.disco",STACK]
    p=finite_rows(ordered_rows(p,"method",methods),"r2")
    if p.empty:
        return None
    fig,(a,b)=plt.subplots(2,1,figsize=(WIDTH,5.7),gridspec_kw={"height_ratios":[1.8,1]})
    fig.subplots_adjust(left=.32,right=.94,top=.93,bottom=.10,hspace=.58)
    panel(a,"A","Prediction on the full LDL test set")
    for i,(m,r) in enumerate(p.iterrows()):
        color=TEAL if m==STACK or m in ["block.glmnet.disco","blocknet.disco"] else BLUE
        a.scatter(r.r2,i,s=28,marker="D" if m==STACK else "o",color=color,zorder=3)
        a.annotate(f"{r.r2:.5f}",(r.r2,i),xytext=(7,0),textcoords="offset points",va="center",fontsize=7.7)
    a.set_yticks(range(len(p)),[LABEL[m] for m in p.index]); a.invert_yaxis()
    a.margins(x=.20,y=.10);axis(a,"Predictive R²")
    panel(b,"B","Paired regional gain over the matched global model")
    tests=[("universal.auto","universal.global.auto","Auto-φ blocks − matched sum"),
           ("universal.fixed","universal.global.fixed","Fixed-φ blocks − matched sum"),
           ("direct.blocknet","direct.global_mlp","Direct BlockNet − global MLP"),
           ("blocknet.disco","blocknet.disco.global","Regional − global Disco network"),
           ("block.glmnet.disco","block.glmnet.disco.offset","Regional elastic net − offset")]
    labels=[]
    for m,ref,label in tests:
        q=subset(d["contrasts"],endpoint="ldl",ancestry="ALL",method=m,reference=ref,metric="predictive_R2")
        if len(q):
            r=q.iloc[0];interval(b,r.difference,r.ci_low,r.ci_high,len(labels));labels.append(label)
    b.set_yticks(range(len(labels)),labels);b.set_ylim(len(labels)-.45,-.55)
    axis(b,"ΔR² (regional minus global)",True)
    if not labels:blank(b,"Matched contrasts unavailable")
    return fig


def ancestry_figure(d, traits):
    if "universal.auto" not in d["performance"].method.values:
        return None
    pairs=d["contrasts"]
    base=d["performance"]
    fig=plt.figure(figsize=(WIDTH,4.65 if "t2dm" in traits else 2.9))
    gs=fig.add_gridspec(2 if "t2dm" in traits else 1,len(traits),left=.095,right=.975,
                        top=.90,bottom=.12,wspace=.38,hspace=.78)
    for j,t in enumerate(traits):
        a=fig.add_subplot(gs[0,j]);panel(a,chr(65+j),NAMES[t])
        metric="AUROC" if t=="t2dm" else "predictive_R2"
        for k,ref in enumerate([CSX,DISCO,STACK]):
            sl=subset(pairs,endpoint=t,method="universal.auto",reference=ref,metric=metric)
            for i,(pop,r) in enumerate(ordered_rows(sl,"ancestry",POPS).iterrows()):
                interval(a,r.difference,r.ci_low,r.ci_high,i+(k-1)*.23,*STYLE[ref])
        a.set_yticks(range(len(POPS)),POPS if j==0 else [""]*len(POPS));a.set_ylim(3.55,-.6)
        axis(a,"ΔAUROC" if t=="t2dm" else "ΔR²",True)
    if "t2dm" in traits:
        a=fig.add_subplot(gs[1,:]);a.axis("off")
        panel(a,chr(65+len(traits)),"Type 2 diabetes discrimination within each ancestry")
        sl=subset(base,endpoint="t2dm");xs=[.0,.23,.43,.60,.77,.95]
        for x,s in zip(xs,["Ancestry","Cases / evaluable","PRS-CSx","DiscoDivas","Four-channel","Universal"]):
            a.text(x,.84,s,transform=a.transAxes,weight="bold",ha="left" if x==.0 else "center",fontsize=8)
        a.plot([0,1],[.78,.78],transform=a.transAxes,color=LIGHT,lw=.8)
        for i,pop in enumerate(POPS):
            rr=ordered_rows(subset(sl,ancestry=pop),"method",[CSX,DISCO,STACK,"universal.auto"])
            r=rr.loc[CSX]
            vals=[pop,f"{fmt(r.n_cases)} / {fmt(r.n)}",*[fmt(v,4) for v in rr.auroc]]
            for j,(x,s) in enumerate(zip(xs,vals)):
                a.text(x,.60-i*.19,s,transform=a.transAxes,ha="left" if j==0 else "center",
                       color=TEAL if j==5 else INK,fontsize=8)
    leg=[Line2D([],[],ls="-",color=STYLE[m][0],marker=STYLE[m][1],ms=4,lw=1.4,label="vs "+LABEL[m])
         for m in [CSX,DISCO,STACK]]
    fig.legend(handles=leg,loc="lower center",frameon=False,ncol=3,bbox_to_anchor=(.52,-.003),fontsize=8)
    return fig


def twins_figure(d):
    w=d["ldl_twins"]
    w=w[w.trait.eq("ldl")]
    if w.empty or "policy" not in w:
        return None
    full=finite_rows(ordered_rows(subset(w,ancestry="ALL"),"policy",POLICY),"mse")
    if full.empty:
        return None
    fig=plt.figure(figsize=(WIDTH,5.25))
    g=fig.add_gridspec(2,2,width_ratios=[1,1.1],height_ratios=[1.6,1],left=.175,right=.97,
                      top=.92,bottom=.115,wspace=.62,hspace=.60)
    a=fig.add_subplot(g[0,0]);b=fig.add_subplot(g[0,1]);c=fig.add_subplot(g[1,0]);e=fig.add_subplot(g[1,1])
    panel(a,"A","Fixed BlockNet error")
    for i,(k,r) in enumerate(full.iterrows()):
        color=PURPLE if k.startswith("risk") else TEAL
        a.scatter(r.mse,i,color=color,s=22)
        a.annotate(f"{r.mse:.4f}",(r.mse,i),xytext=(4,0),textcoords="offset points",va="center",fontsize=7)
    a.set_yticks(range(len(full)),[POLICY[k] for k in full.index])
    a.set_ylim(len(full)-.5,-.5);a.margins(x=.30);axis(a,"MSE, (mmol/L)²")
    panel(b,"B","Test-roster coverage")
    rules=["raw_global","regional_count","latent_count","risk_only","risk_count","risk_local"]
    pops=[i for i in [*POPS,"OTH"] if i in w.ancestry.values]
    matrix=w.pivot(index="policy",columns="ancestry",values="coverage").reindex(index=rules,columns=pops)*100
    bound=max(25.,float(np.nanmax(matrix.values)) if np.isfinite(matrix.values).any() else 25.)
    cmap=matplotlib.colors.LinearSegmentedColormap.from_list("coverage",["#f4f7f7","#91c9be","#006b60"])
    cmap.set_bad("#e6e6e6")
    b.imshow(matrix,aspect="auto",vmin=0,vmax=bound,cmap=cmap,interpolation="none")
    for i,arr in enumerate(matrix.values):
        for j,v in enumerate(arr):
            b.text(j,i,fmt(v,1),ha="center",va="center",color="white" if v>bound*.7 else INK,fontsize=6.7)
    b.set_xticks(range(len(pops)),pops);b.set_yticks(range(len(rules)),["Global","Regional","Latent","Risk only","+ count","+ local"],fontsize=7)
    b.set_xlabel("Selected within ancestry (%)");b.tick_params(length=0,pad=2);b.spines[:].set_visible(False)
    panel(c,"C","Gain beyond risk only")
    ref=full.reindex(["risk_count","risk_local"])
    for i,r in enumerate(ref.itertuples()):
        interval(c,r.selected_minus_risk_only_mse,r.risk_only_ci_low,r.risk_only_ci_high,i,PURPLE)
    c.set_yticks([0,1],["Risk + count","Risk + local features"]);c.set_ylim(1.55,-.55)
    axis(c,"ΔMSE, (mmol/L)²",True)
    panel(e,"D","Matched TRAIN donors")
    sub=full.reindex(["all","latent_count"]).dropna(subset=["median_latent_twins"])
    for i,(name,row) in enumerate(sub.iterrows()):
        val=row.median_latent_twins
        e.barh(i,val,height=.40,color=TEAL if name=="latent_count" else "#7c919b")
        e.annotate(fmt(val),(val,i),xytext=(4,0),textcoords="offset points",va="center",fontsize=8)
    e.set_yticks(range(len(sub)),["All" if k=="all" else "Latent" for k in sub.index])
    e.invert_yaxis();e.margins(x=.22,y=.32);axis(e,"Median latent twin count")
    return fig


def curves_figure(d, classification, traits):
    curves=d["curves"]
    if curves.empty:
        return None
    endpoints=[t for t in (("ldl.high","t2dm") if classification else TRAITS)
               if t in curves.endpoint.values and (("ldl" if t=="ldl.high" else t) in traits)]
    if not endpoints:return None
    rows=len(endpoints) if classification else 1
    cols=2 if classification else len(endpoints)
    fig,axs=plt.subplots(rows,cols,squeeze=False,figsize=(WIDTH,2.6*rows+.35))
    fig.subplots_adjust(left=.105,right=.975,top=.88 if rows==1 else .93,bottom=.28 if rows==1 else .15,
                        wspace=.46,hspace=.6)
    for j,t in enumerate(endpoints):
        for k,mode in enumerate(["ROC","precision–recall"] if classification else ["calibration"]):
            a=axs[j,k] if classification else axs[0,j]
            letter=chr(65+(j*2+k if classification else j))
            panel(a,letter,NAMES[t]+(" · "+("ROC" if k==0 else "PR") if classification else ""))
            for m,ls in zip(STYLE,[":","--","-.","-"]):
                q=subset(curves,endpoint=t,ancestry="ALL",method=m,curve=mode)
                # ROC/PR rows retain upstream point order, including vertical segments.
                if len(q):
                    a.plot(q.x,q.y,ls=ls,color=STYLE[m][0],lw=1.25,marker=None if classification else STYLE[m][1],ms=2.8)
            if classification:
                a.set(xlim=(0,1),ylim=(0,1),xlabel="Recall" if k else "False-positive rate",ylabel="Precision" if k else "Sensitivity")
                if k==0: a.plot([0,1],[0,1],ls=":",color=MUTED,lw=.65,zorder=0)
                else:
                    q=subset(d["performance"],endpoint=t,ancestry="ALL",method=CSX)
                    if len(q):a.axhline(q.event_fraction.iloc[0],ls=":",color=MUTED,lw=.65)
            else:
                ran=[min(a.get_xlim()[0],a.get_ylim()[0]),max(a.get_xlim()[1],a.get_ylim()[1])]
                a.plot(ran,ran,ls=":",color=MUTED,lw=.7,zorder=0)
                a.set(xlim=ran,ylim=ran,xlabel="Mean prediction",ylabel="Observed mean")
                a.xaxis.set_major_locator(MaxNLocator(4));a.yaxis.set_major_locator(MaxNLocator(4))
    handles=[Line2D([],[],color=STYLE[m][0],ls=ls,lw=1.6,label=LABEL[m]) for m,ls in zip(STYLE,[":","--","-.","-"]) if m in curves.method.values]
    fig.legend(handles=handles,ncol=3,frameon=False,loc="lower center",bbox_to_anchor=(.5,.005))
    return fig


def abm_figure(d, traits):
    s=d["screening"]
    if s.empty:return None
    selections=["selected_twin_count_0p20","selected_support_0p20","selected_risk_error_0p20","selected_support_error_0p20"]
    fig,axs=plt.subplots(1,len(traits),squeeze=False,figsize=(WIDTH,2.6))
    fig.subplots_adjust(left=.21,right=.975,top=.81,bottom=.22,wspace=.38)
    for j,t in enumerate(traits):
        a=axs[0,j];panel(a,chr(65+j),NAMES[t])
        q=ordered_rows(subset(s,trait=t,method="abm"),"selection",selections)
        for i,r in enumerate(q.itertuples()):
            interval(a,r.selected_minus_all_loss,r.selection_low,r.selection_high,i,PURPLE if i>=2 else TEAL)
        axis(a,"ΔBrier score" if t=="t2dm" else "ΔMSE",True)
        a.set_yticks(range(4),["Twin count","Family support","Risk-error model","Support-error model"] if j==0 else [""]*4)
        a.set_ylim(3.6,-.6)
    return fig


def table_workbook(tables, path):
    """TSVs remain available in installations without the optional Excel writer."""
    try:
        import xlsxwriter
    except ImportError:
        print("Excel export skipped (XlsxWriter unavailable). All tables are available as TSVs.")
        return False
    with xlsxwriter.Workbook(str(path), {"strings_to_formulas":False,"strings_to_urls":False}) as book:
        header=book.add_format({"bold":True,"bg_color":"#1f4651","font_color":"white","text_wrap":True,"valign":"vcenter"})
        title=book.add_format({"bold":True,"font_size":12,"font_color":TEAL})
        index=book.add_worksheet("Contents")
        index.set_column(0,0,16);index.set_column(1,1,85)
        index.write(0,0,"GRID manuscript tables",title)
        index.write_row(2,0,["Table","Description"],header)
        for i,(sheet,(name,_)) in enumerate(tables.items(),3):
            index.write_url(i,0,f"internal:'{sheet}'!A1",string=sheet)
            index.write(i,1,name)
        index.freeze_panes(3,0);index.hide_gridlines(2)
        for sheet,(name,frame) in tables.items():
            ws=book.add_worksheet(sheet)
            ws.write(0,0,name,title);ws.hide_gridlines(2)
            for i,col in enumerate(frame.columns):
                ws.write(2,i,col,header)
                width=max(len(str(col)),frame[col].astype(str).str.len().quantile(.9) if len(frame) else 0)
                ws.set_column(i,i,min(40,max(12,width+2)))
            ws.set_row(2,32)
            for i,row in enumerate(frame.itertuples(index=False,name=None),3):
                for j,v in enumerate(row):
                    if pd.notna(v) and not (isinstance(v,float) and not np.isfinite(v)):
                        ws.write(i,j,v)
            ws.freeze_panes(3,2);ws.autofilter(2,0,2+len(frame),len(frame.columns)-1)
    return True


def load_inputs(results):
    data={}
    for name in [*SUPP,"curves"]:
        path=results/(name+".tsv")
        if name in {"performance","contrasts","cohort"} and not path.is_file():
            raise FileNotFoundError(f"Missing required GRID aggregate: {path}")
        data[name]=pd.read_csv(path,sep="\t",low_memory=False) if path.is_file() else pd.DataFrame()
    p,c=data["performance"],data["contrasts"]
    keys=["trait","endpoint","ancestry","method"]
    if p.duplicated(keys).any() or c.duplicated(keys+["reference","metric"]).any():
        raise ValueError("Repeated method or paired contrast in a TEST stratum.")
    if p.groupby(keys[:3]).n.nunique().gt(1).any():
        raise ValueError("Methods have different evaluable TEST counts in the same stratum.")
    # A single-trait export can omit columns for the other endpoint type.
    for col in ["auroc","n_cases","event_fraction","r2","outcome_mean","outcome_variance"]:
        if col not in p: p[col]=np.nan
    if data["cohort"].duplicated(["role","ancestry"]).any():
        raise ValueError("Repeated cohort role and ancestry.")
    for name,keys2 in [("ldl_twins",["ancestry","policy"]),("screening",["trait","method","selection"])]:
        df=data[name]
        if all(k in df for k in keys2) and df.duplicated(keys2).any():
            raise ValueError("Repeated selection row in "+name)
    return data


def validation_note(d):
    info=d["model_settings"]
    if "validation_status" not in info:return "Validation status was not recorded in the supplied aggregates."
    status=info.validation_status.dropna().unique()
    if list(status)==["exploratory_previously_inspected_test"]:
        return "Results are exploratory. This internal UK Biobank test set was previously inspected."
    return "Recorded validation status: "+", ".join(map(str,status))+"."


def build(results, out, threshold=None, raster=True):
    from matplotlib.transforms import Bbox
    d=load_inputs(results)
    traits=[t for t in TRAITS if t in d["performance"].trait.values]
    if not traits:raise ValueError("No supported trait in performance.tsv")
    t1=make_table1(d["performance"],d["cohort"],d["test_demographics"])
    figures=out/"figures";tables=out/"tables"
    figures.mkdir();tables.mkdir()
    books={"Table 1":("Study population and test-set phenotypes",t1)}
    t1.to_csv(tables/"Table_1.tsv",sep="\t",index=False)
    for key,(sheet,title) in SUPP.items():
        if not d[key].empty:
            frame=d[key]
            dest=tables/f"Table_{sheet}_{key}.tsv"
            if key=="ld_maps" and "bed_files" in frame:
                frame=frame.copy()
                frame["bed_files"]=frame.bed_files.map(lambda v: "|".join(Path(z).name for z in str(v).split("|")) if pd.notna(v) else "")
            if key=="model_settings" and "phenotype" in frame:
                frame=frame.copy()
                def describe_phenotype(value):
                    if pd.isna(value): return ""
                    try: record=json.loads(value)
                    except (ValueError, TypeError): return str(value)
                    if not isinstance(record,dict): return str(value)
                    return " / ".join(str(record[k]) for k in ("source_column","transformation","endpoint","administrative_end_date") if k in record)
                frame["phenotype"]=frame.phenotype.map(describe_phenotype)
            if frame is d[key]:
                shutil.copyfile(results/(key+".tsv"),dest)
            else:
                frame.to_csv(dest,sep="\t",index=False)
            books[sheet]=(title,frame)
    if not d["curves"].empty:
        q=d["curves"]
        q=q[q.method.isin(STYLE)&q.ancestry.eq("ALL")&q.curve.isin(["calibration","ROC","precision–recall"])]
        q.to_csv(tables/"Table_S7_plotted_curves.tsv",sep="\t",index=False)
        books["S7"]=("Source points for supplementary calibration and discrimination curves",q)
    table_workbook(books,out/"Manuscript_tables.xlsx")
    status=validation_note(d)
    thresh=(f" at LDL ≥ {threshold:g} mmol/L" if threshold is not None else
            " at the threshold used in the original evaluation (value not supplied)")
    abm=d["model_settings"]
    alpha=(", ".join(f"{NAMES.get(r.trait,r.trait)} {fmt(r.alpha,2)}" for r in abm.itertuples())
           if "alpha" in abm else "unavailable")
    ci="Whiskers show 95% percentile intervals from paired family bootstrap resampling. A missing interval is not estimable."
    docs=[
        ("Table_1","Table 1","Study population and test-set phenotypes",
         "Cohort roles refer to the full roster. Phenotype summaries use participants with an observed outcome in TEST. Outcomes are observed height, measured LDL cholesterol and baseline prevalent type 2 diabetes. "+status,lambda:figure_table(t1)),
        ("Fig_1_atlas","Figure 1","From a genome-wide score to learnable LD-block profiles",
         "A. The same posterior allelic effects are used for every target participant within a trait. B. Each column is one universal LD block in genomic order. Colour shows log10(1 + Vb / mean Vb), where Vb is the across-participant TRAIN variance of that block score within a trait. Empty and constant blocks contribute zero. C. Blocks are ranked separately by trait. D. Effective block count is 1 / sum(pb²), with pb = Vb / sum Vb. These shares describe score dispersion; correlated blocks can also contribute covariance to the genome-wide sum.",lambda:atlas_figure(d,traits)),
        ("Fig_2_benchmarks","Figure 2","One-weight GRID and comparator performance on the same participants",
         "Upper panels report full-model predictive R² or AUROC. Universal GRID learns block coefficients from auto-shrinkage META weights. The matched global model uses the same covered sum, its square, age, sex, squared age and pairwise interactions among age, sex and the sum. Comparator pipelines use their own covariates. The four-channel comparator is block.stack. Lower panels use paired differences. "+ci+" "+status,lambda:benchmark_figure(d,traits)),
        ("Fig_3_ldl","Figure 3","Global nonlinearity and the contribution of regional information to LDL prediction",
         "Panel A shows full-model predictive R². For the universal models, matched pairs compare regional learning to its own covered global sum at auto or fixed shrinkage. Panel B also compares the direct BlockNet versus the global MLP, the Disco regional network versus its global counterpart, and the Disco regional elastic net versus its matched global offset. The offset includes smooth age effects and global score interactions. Each pair uses the same test participants. The stack is displayed as a point estimate in A. "+ci,lambda:ldl_figure(d) if "ldl" in traits else None),
        ("Fig_4_ancestry","Figure 4","Variation in predictive gain across genetic ancestry groups",
         "Points show the one-weight universal model minus each comparator (PRS-CSx, DiscoDivas and the four-channel stack) for each target ancestry group. Each panel uses its own horizontal scale. Positive values favour the one-weight universal model. OTH is retained in the full supplementary tables. When available, the diabetes table gives absolute AUROC and the number of prevalent cases relative to the evaluable sample. "+ci,lambda:ancestry_figure(d,traits)),
        ("Fig_5_twins","Figure 5","Training-neighbour counts and error screening for LDL",
         "A. The LDL predictor is the same frozen four-channel direct BlockNet in every row. Only the test subset changes. B. Coverage uses the entire roster within ancestry, including people without an LDL outcome. The global count uses a pooled cutoff. Regional and latent count policies use ancestry-specific cutoffs. Risk only predicts squared error from predicted LDL, age, sex and PC1–4. Risk + count adds donor and family counts; the local model also uses distance and donor LDL dispersion. C. Differences are selected-group MSE minus risk-only MSE. Family resampling preserves overlap between groups. D. Median latent twin counts use the full or selected roster, including participants without observed LDL.",lambda:twins_figure(d) if "ldl" in traits else None),
        ("Fig_S1_calibration","Figure S1","Calibration of the full-model predictions",
         "Observed means are plotted against mean predictions in test groups defined using calibration-set score deciles. The diagonal indicates agreement. Height is in cm, LDL in mmol/L and diabetes on the probability scale.",lambda:curves_figure(d,False,traits)),
        ("Fig_S2_discrimination","Figure S2","ROC and precision–recall curves for the binary endpoints",
         "The binary endpoints, when available, are high LDL cholesterol"+thresh+" and prevalent type 2 diabetes. Curves use the released display points. Tabulated average precision and AUROC use the full scores. Dotted grey lines show chance ROC performance or the observed event fraction. High LDL is ranked using the continuous LDL prediction.",lambda:curves_figure(d,True,traits)),
        ("Fig_S3_abm","Figure S3","Selective prediction with the primary ABM across traits",
         "Points show the difference in loss between the selected group and the full test set, for an unchanged ABM predictor. Policies target 20% coverage; actual coverage, case retention and risk-standardized comparisons are in Table S5. Lower values indicate lower selected-group error. For a binary outcome, the selected group's prevalence affects Brier score. Fitted ABM residual weights (alpha): "+alpha+". "+ci,lambda:abm_figure(d,traits))]
    docs.append(("Fig_S4_four_channel_design","Figure S4","Architecture of the four-channel comparator",
                 "Each discovery ancestry supplies a channel. The interpolated score and four calibrated components form the previously evaluated block.stack comparator.",lambda:design_figure(d)))
    records=[]
    with plt.rc_context(RC), PdfPages(out/"GRID_figure_proof.pdf",metadata={"Title":"GRID manuscript figures and tables"}) as proof:
        for stem,number,title,caption,draw in docs:
            fig=draw()
            if fig is None:
                records.append(dict(number=number,title=title,status="source unavailable",caption=caption,file=""))
                continue
            folder=tables if number=="Table 1" else figures
            target=folder/(stem+".pdf")
            fig.savefig(target,dpi=600)
            if number!="Table 1" and raster:
                fig.savefig(folder/(stem+".tiff"),dpi=600,pil_kwargs={"compression":"tiff_lzw"})
            w,h=fig.get_size_inches()
            lines=textwrap.wrap(caption,width=int(14.4*w))
            foot=.28+.17*len(lines)
            if number!="Table 1":
                fig.text(.025,1+.10/h,number+"   "+title,fontsize=9.1,weight="bold",va="bottom")
            fig.text(.025,-.10/h,"\n".join(lines),fontsize=7.9,va="top",linespacing=1.45)
            proof.savefig(fig,bbox_inches=Bbox.from_bounds(0,-foot,w,h+foot+(.08 if number=="Table 1" else .39)),pad_inches=0)
            plt.close(fig)
            records.append(dict(number=number,title=title,status="available",caption=caption,file=target.relative_to(out).as_posix()))
    pd.DataFrame(records).to_csv(out/"Figure_legends.tsv",sep="\t",index=False)
    pd.DataFrame([{"table":sheet,"description":name,"rows":len(frame)} for sheet,(name,frame) in books.items()]).to_csv(out/"Table_index.tsv",sep="\t",index=False)
    links="".join(f"<li><b>{html.escape(r['number'])}.</b> "+(f"<a href=\"{r['file']}\">{html.escape(r['title'])}</a>" if r["file"] else html.escape(r["title"])+" (source unavailable)")+f"<p>{html.escape(r['caption'])}</p></li>" for r in records)
    table_links="".join(f'<li><a href="tables/{f.name}">{html.escape(f.name)}</a></li>' for f in sorted(tables.glob("*.tsv")))
    wb='<a href="Manuscript_tables.xlsx">Manuscript tables (Excel)</a> · ' if (out/"Manuscript_tables.xlsx").is_file() else ""
    (out/"index.html").write_text('<!doctype html><meta charset="utf-8"><title>GRID manuscript</title><style>body{font:16px/1.5 system-ui;max-width:54rem;margin:3rem auto;color:#182e3a}li{margin:.5em 0}a{color:#007867}p{margin:.2em 0 .8em}</style><h1>GRID manuscript</h1><p>'+html.escape(status)+'</p><p>'+wb+'<a href="GRID_figure_proof.pdf">Figure proof and legends</a></p><ol>'+links+'</ol><h2>Tables</h2><ul>'+table_links+'</ul>',encoding="utf-8")
    return records


def render(results, out=None, threshold=None, raster=True):
    results=Path(results).resolve()
    out=Path(out).resolve() if out else results.parent/"manuscript"
    if out==results or out in results.parents:
        raise ValueError("Choose an output separate from the aggregate input directory.")
    out.parent.mkdir(parents=True,exist_ok=True)
    with (out.parent/("."+out.name+".lock")).open("a") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        with tempfile.TemporaryDirectory(dir=out.parent,prefix="."+out.name+"-") as tmp:
            stage=Path(tmp)/"new";stage.mkdir()
            def snapshot():
                return {p.name:(p.stat().st_size,p.stat().st_mtime_ns) for p in results.glob("*.tsv")}
            before=snapshot()
            records=build(results,stage,threshold,raster)
            if before!=snapshot():
                raise RuntimeError("Aggregate tables changed during rendering. Repeat the render.")
            # Install the complete build together. Restore the prior edition on a failed rename.
            # Retain untracked files if the user added anything to the managed output.
            marker=".grid-generated.json"
            files=sorted(str(f.relative_to(stage)) for f in stage.rglob("*") if f.is_file())
            (stage/marker).write_text(json.dumps(files,indent=2)+"\n")
            backup=Path(tmp)/"prior"
            if out.exists():
                prior=set(json.loads((out/marker).read_text())) if (out/marker).is_file() else set()
                for f in out.rglob("*"):
                    rel=f.relative_to(out)
                    if f.is_file() and str(rel) not in prior and str(rel)!=marker:
                        destination=stage/rel
                        if destination.exists():
                            raise FileExistsError(f"An existing unmanaged file would be overwritten: {out/rel}")
                        destination.parent.mkdir(parents=True,exist_ok=True)
                        shutil.copy2(f,destination)
                out.rename(backup)
            try:
                stage.rename(out)
            except BaseException:
                if backup.exists(): backup.rename(out)
                raise
    print(f"GRID manuscript: {out} ({sum(bool(r['file']) for r in records)} numbered exhibits)")
    return out


def add_run_context(root, results):
    """Write aggregate context available in a live GRID run to the render inputs."""
    from core2 import atomic_table, read_table
    from output_layout import project_path
    p=results/"test_demographics.tsv"
    participants=project_path(root,"tmp/cohort/participants.tsv.gz")
    if participants.is_file():
        co=read_table(participants,["role","ancestry","age","sex"])
        co=co[co.role.eq("test")].copy()
        rows=[]
        for ancestry,g in [("ALL",co),*co.groupby("ancestry")]:
            r=dict(ancestry=ancestry,n_roster=len(g))
            if "age" in g:
                age=pd.to_numeric(g.age,errors="coerce").dropna()
                r.update(age_n=len(age),age_mean=age.mean(),age_sd=age.std())
            if "sex" in g:
                codes=pd.to_numeric(g.sex,errors="coerce").dropna()
                r["sex_n"]=len(codes)
                for value,n in codes.value_counts().items():
                    r["sex_code_"+format(value,"g")]=int(n)
            rows.append(r)
        atomic_table(pd.DataFrame(rows),p)
    else:
        p.unlink(missing_ok=True)
    performance=pd.read_csv(results/"performance.tsv",sep="\t",usecols=["trait","method"])
    available=pd.read_csv(results/"availability.tsv",sep="\t")
    from methods import marker,published
    for table,source in [("universal_atlas","atlas.tsv"),("universal_coefficients","block_coefficients.tsv")]:
        frames=[]
        for t in performance.trait.unique():
            for mode in ("auto","fixed"):
                name=("universal.prepare." if table=="universal_atlas" else "universal.")+mode
                if not len(subset(available,trait=t,method=name,status="complete")):
                    continue
                rec=marker(root,t,name)
                if published(rec):
                    frames.append(pd.read_csv(rec.parent/source,sep="\t").assign(mode=mode))
        path=results/(table+".tsv")
        if frames: atomic_table(pd.concat(frames,ignore_index=True),path)
        else: path.unlink(missing_ok=True)
    rows=[]
    for t in performance.loc[performance.method.eq(STACK),"trait"].unique():
        f=root/t/STACK/"model.json"
        if f.is_file():
            info=json.loads(f.read_text())
            for method,w in info["weights"].items():
                rows.append(dict(trait=t,component=method,weight=w,calibration_n=info.get("calibration"),
                                 objective=info.get("objective"),test_labels_used=info.get("test_labels_used")))
    p=results/"stack_weights.tsv"
    if rows:
        atomic_table(pd.DataFrame(rows),p)
    else:
        p.unlink(missing_ok=True)


def publish_manuscript(cfg,root):
    from file_lock import file_lock
    from methods_report import publish
    from output_layout import project_path
    publish(cfg,root)
    paper=project_path(root,"grid.paper/main")
    cache=project_path(root,"grid.paper/logs")/"cache"
    with file_lock(cache/"merge.lock"):
        add_run_context(root,paper/"results")
        return render(paper/"results",threshold=cfg["evaluation"].get("ldl_threshold",4.9))


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--results",type=Path,required=True,help="Existing paper/main/results TSV directory")
    cli.add_argument("--out",type=Path,help="Default: manuscript/ beside results/")
    cli.add_argument("--ldl-threshold",type=float,help="Original threshold for the high-LDL caption only; does not change classification")
    cli.add_argument("--pdf-only",action="store_true",help="Skip 600-dpi TIFF export")
    args=cli.parse_args()
    threshold=args.ldl_threshold
    if threshold is None:
        # Look for the evaluation receipt, including named GRID runs.
        root=args.results.resolve().parent
        for parent in (root,*root.parents):
            if parent.name=="main" and parent.parent.name=="paper":
                receipt=parent.parent/"logs"/root.relative_to(parent)/"cache/ldl/inventory.json"
                if receipt.is_file():
                    threshold=json.loads(receipt.read_text()).get("evaluation",{}).get("ldl_threshold")
                break
    render(args.results,args.out,threshold,not args.pdf_only)


if __name__=="__main__":
    main()
