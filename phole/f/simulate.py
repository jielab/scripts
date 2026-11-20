#!/usr/bin/env python3
"""Create explicitly synthetic four-cohort full-GWAS fixtures for PHOLE.

These are controlled numerical experiments, not human GWAS findings. Background
protein loadings are independent random draws and never functions of gene
position. Scenarios differ only in the planted sparse programs, preserving the
background, sampling noise, missingness and allele transformations for a seed.
"""
from __future__ import annotations

import argparse
import gzip
import io
import math
from pathlib import Path

import numpy as np

from signed_input import DATASETS, atomic_json, atomic_tsv, chrom_order

HEADER=["SNP","CHR","POS","EA","NEA","EAF","N","BETA","SE","P"]
COHORTS={"ukb":{"n":54000,"ancestry":"EUR"},"ckb":{"n":4000,"ancestry":"EAS"},
         "decode":{"n":35000,"ancestry":"EUR"},"scot":{"n":2000,"ancestry":"EUR"}}


def arguments():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--scenario",choices=["null","latent","platform_artifact"],default="latent")
    p.add_argument("--seed",type=int,default=20260828)
    p.add_argument("--n-source",type=int,default=120)
    p.add_argument("--n-protein",type=int,default=48)
    p.add_argument("--missing",type=float,default=0.10)
    return p.parse_args()


def write_full_summary(path,rows):
    """Deterministic gzip bytes (no timestamp or path embedded in the header)."""
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="",mode="wb",fileobj=raw,mtime=0) as zipped:
            with io.TextIOWrapper(zipped,encoding="utf-8",newline="") as stream:
                stream.write("\t".join(HEADER)+"\n")
                for row in rows:
                    stream.write("\t".join(format(v,".12g") if isinstance(v,(float,np.floating)) else str(v) for v in row)+"\n")


def generate(args):
    if args.n_source<100 or args.n_source<args.n_protein+8:
        raise ValueError("Need at least 100 source positions and at least n-protein+8 to reserve unmeasured source programs")
    if args.n_protein<16 or not 0<=args.missing<0.5:
        raise ValueError("Need at least 16 proteins and missingness in [0,0.5)")
    out=args.out
    if out.exists() and any(out.iterdir()) and not (out/"SYNTHETIC_FIXTURE.txt").exists():
        raise ValueError("Refusing to overwrite a nonempty directory not marked as a synthetic fixture")
    out.mkdir(parents=True,exist_ok=True)
    (out/"SYNTHETIC_FIXTURE.txt").write_text("SYNTHETIC DATA ONLY. Never report these planted signals as biological discoveries.\n")
    rng=np.random.default_rng(args.seed)
    m,n=args.n_source,args.n_protein
    chromosomes=[str(i) for i in range(1,23)]+["X","Y"]
    chrom=np.array([chromosomes[i%24] for i in range(m)],dtype=object)
    pos=np.array([1000001+(i//24)*6000000 for i in range(m)],dtype=np.int64)
    snp=np.array([f"rs990000{100000+i}" for i in range(m)],dtype=object)
    order=sorted(range(m),key=lambda i:(chrom_order(chrom[i]),pos[i]))
    alleles=[("A","G"),("C","T"),("A","C"),("G","T")]
    ea=np.empty(m,dtype=object);nea=np.empty(m,dtype=object)
    base_eaf=rng.uniform(0.08,0.42,m)
    for i in range(m):
        ea[i],nea[i]=("A","T") if i%23==0 else alleles[i%len(alleles)]
        if i%23==0:base_eaf[i]=0.12
    proteins=np.array([f"SYNTH_PROT_{i:03d}" for i in range(n)],dtype=object)
    coding_sources=rng.choice(m,size=n,replace=False)
    unmeasured=np.setdiff1d(np.arange(m),coding_sources)
    program_sources=[];program_blocks=set()
    for i in rng.permutation(unmeasured):
        block=(chrom[i],"MHC") if chrom[i]=="6" and 25000000<=pos[i]<=35000000 else (chrom[i],int((pos[i]-1)//5000000))
        if block in program_blocks:continue
        program_blocks.add(block);program_sources.append(int(i))
        if len(program_sources)==4:break
    program_sources=np.array(program_sources,dtype=int)
    scot_proteins=set(rng.choice(n,size=max(12,int(math.ceil(n*0.35))),replace=False).tolist())

    # Exactly rank-two random background. A modest minimum row amplitude keeps
    # the toy panel broad enough for block-matched null sampling after the
    # discovery significance filter. It does not introduce a sparse program.
    source_loadings=rng.normal(0,0.04,(m,2))
    protein_loadings=rng.normal(0,1,(2,n))
    background=source_loadings@protein_loadings
    strength=np.max(np.abs(background),axis=1)
    rescale=np.maximum(1.,0.12/np.maximum(strength,1e-8))
    source_loadings*=rescale[:,None]
    background=source_loadings@protein_loadings
    cis_effect=rng.uniform(0.30,0.45,n)
    baseline=background.copy()
    for j,i in enumerate(coding_sources):baseline[i,j]+=cis_effect[j]

    planted=np.zeros((m,n),dtype=float)
    programs=[]
    for h,i in enumerate(program_sources):
        targets=np.sort(rng.choice(n,size=min(16,n),replace=False))
        effect=rng.choice([-1.,1.],len(targets))*rng.uniform(0.18,0.30,len(targets))
        planted[i,targets]=effect
        programs.append({"PROGRAM":f"SPARSE_{h+1}","KIND":"latent_program","SOURCE_ID":f"SYNTH_LOCUS_{i:04d}",
                         "SNP":snp[i],"CHR":chrom[i],"POS":int(pos[i]),"TARGET_PROTS":";".join(proteins[targets]),
                         "TARGET_COUNT":len(targets),"SOURCE_CODING_PROTEIN_MEASURED":False,
                         "ACTIVE_UKB":args.scenario!="null","ACTIVE_CKB":args.scenario=="latent",
                         "ACTIVE_DECODE":args.scenario=="latent","ACTIVE_SCOT":args.scenario=="latent",
                         "SCOT_COVERED_TARGET_COUNT":sum(int(j) in scot_proteins for j in targets),
                         "SCENARIO":args.scenario})

    annotation=[]
    for j,prot in enumerate(proteins):
        i=coding_sources[j]
        annotation.append({"prot":prot,"gene_symbol":prot,"gene_id":f"SYNTH_GENE_{j:05d}","CHR":chrom[i],
                           "START":int(pos[i]-2000),"END":int(pos[i]+2000),"strand":"+","biotype":"protein_coding",
                           "build":"38","match_method":"synthetic_known_mapping","annotation_status":"single_gene",
                           "component":prot,"source":"SYNTHETIC_FIXTURE","notes":"fictional coding gene, not an annotation claim"})
    atomic_tsv(out/"protein_coordinates.GRCh38.tsv",annotation)
    manifest=[];candidates=[];truth_effects=[];coverage={};summary_paths={}
    for d_index,(d,spec) in enumerate(COHORTS.items()):
        cohort_rng=np.random.default_rng(np.random.SeedSequence([args.seed,1000+d_index]))
        represented=sorted(scot_proteins) if d=="scot" else list(range(n))
        f=base_eaf+cohort_rng.normal(0,0.07 if d=="ckb" else 0.01,m)
        f=np.clip(f,0.03,0.45)
        if d=="ckb":f[np.arange(m)%23==0]=0.25
        else:f[np.arange(m)%23==0]=0.12+cohort_rng.normal(0,0.005,np.sum(np.arange(m)%23==0))
        active=args.scenario=="latent" or (args.scenario=="platform_artifact" and d=="ukb")
        true_beta=baseline+(planted if active else 0.)
        assay_list=[]
        for j in represented:
            assay_list.append((j,str(proteins[j])))
            if j in {2,min(17,n-1)}:assay_list.append((j,str(proteins[j])+".1"))
        nonmissing=0;swap_count=0;complement_count=0
        for j,assay in assay_list:
            measurement_scale=float(cohort_rng.lognormal(0,0.08))
            variant_n=np.full(m,spec["n"],dtype=int)
            variant_n[chrom=="Y"]=spec["n"]//2
            variance=2*f*(1-f);variance[chrom=="Y"]=f[chrom=="Y"]*(1-f[chrom=="Y"])
            se=measurement_scale/np.sqrt(variant_n*variance)
            estimate=true_beta[:,j]+cohort_rng.normal(size=m)*se
            missing=cohort_rng.random(m)<args.missing
            transformations=cohort_rng.choice(4,size=m,p=[0.60,0.20,0.10,0.10])
            rows=[]
            for i in order:
                if missing[i]:continue
                beta=float(estimate[i]);x,y=str(ea[i]),str(nea[i]);eaf=float(f[i])
                p=math.erfc(abs(beta/float(se[i]))/math.sqrt(2))
                transform=int(transformations[i])
                if transform in {1,3}:
                    x,y=y,x;beta=-beta;eaf=1-eaf;swap_count+=1
                if transform in {2,3}:
                    x=x.translate(str.maketrans("ACGT","TGCA"));y=y.translate(str.maketrans("ACGT","TGCA"));complement_count+=1
                rows.append([snp[i],chrom[i],int(pos[i]),x,y,eaf,int(variant_n[i]),beta,float(se[i]),p])
                nonmissing+=1
                if d=="ukb" and p<=5e-8:
                    candidates.append({"DATASET":"ukb","ASSAY":assay,"PROT":str(proteins[j]),"CHR":chrom[i],"POS":int(pos[i]),
                                       "SNP":snp[i],"P":format(p,".12g"),"KIND":"cis" if i==coding_sources[j] else "trans",
                                       "CLASS_STATUS":"synthetic_known_coordinates"})
                if i in program_sources:
                    truth_effects.append({"DATASET":d,"ASSAY":assay,"PROT":str(proteins[j]),"SNP":snp[i],"CHR":chrom[i],"POS":int(pos[i]),
                                          "CANONICAL_EA":ea[i],"CANONICAL_NEA":nea[i],"TRUE_BETA":float(true_beta[i,j]),
                                          "PLANTED_PROGRAM_BETA":float(planted[i,j]) if active else 0.,"ACTIVE_PROGRAM":bool(active),
                                          "SYNTHETIC":True})
            path=out/"gwas"/DATASETS[d]/"clean"/(assay+".gz")
            write_full_summary(path,rows)
            manifest.append({"DATASET":d,"ASSAY":assay,"PROT":str(proteins[j]),"PATH":str(path.resolve()),"BUILD":"38",
                             "ANCESTRY":spec["ancestry"],"EFFECT_SCALE":"synthetic_protein_SD_per_effect_allele","INPUT_SCOPE":"full",
                             "SYNTHETIC":"TRUE","SCOPE_EVIDENCE":"complete defined toy variant universe except simulated missingness",
                             "BUILD_EVIDENCE":"synthetic_coordinates_only"})
            summary_paths[d,assay]=str(path.resolve())
        coverage[d]={"n_proteins":len(represented),"n_assays":len(assay_list),"n_rows":nonmissing,
                     "swap_rows":swap_count,"complement_rows":complement_count,"ancestry":spec["ancestry"],"sample_size":spec["n"]}
    atomic_tsv(out/"assay_manifest.tsv",manifest)
    atomic_tsv(out/"candidate_loci.tsv",candidates)
    atomic_tsv(out/"truth.tsv",programs)
    atomic_tsv(out/"truth_effects.tsv",truth_effects)
    mapping=[]
    for j,prot in enumerate(proteins):
        row={"prot":prot}
        for d in DATASETS:
            row[d]=";".join(a['ASSAY'] for a in manifest if a['DATASET']==d and a['PROT']==prot) or "."
        mapping.append(row)
    atomic_tsv(out/"prot.names.tsv",mapping)
    cis_rows=[{"PROT":str(prot),"SNP":snp[coding_sources[j]],"CHR":chrom[coding_sources[j]],"POS":int(pos[coding_sources[j]]),
               "TRUE_CIS_BETA":float(cis_effect[j]),"SYNTHETIC":True} for j,prot in enumerate(proteins)]
    atomic_tsv(out/"cis_controls.tsv",cis_rows)
    atomic_json(out/"synthetic.json",{"synthetic":True,"schema":"phole-synthetic-full-gwas-v2","scenario":args.scenario,"seed":args.seed,
                                     "n_source":m,"n_protein":n,"background_rank":2,"sparse_programs":4 if args.scenario!="null" else 0,
                                     "missing_probability":args.missing,"coverage":coverage,"genome_build":"38",
                                     "positions":"fictional sentinel positions including X/Y; no real gene or rsID claim",
                                     "background":"rank two; protein loadings are independent random draws, unrelated to coding positions",
                                     "noise":"independent cohort and assay errors; SE from N and genotype variance",
                                     "scenario_interpretation":{"null":"background rank two plus cis controls; no added sparse program",
                                                                "latent":"same planted sparse programs in all four cohorts",
                                                                "platform_artifact":"sparse programs in UKB only; external cohorts retain background and cis controls"}[args.scenario],
                                     "validation_limit":"a numerical pipeline check, not evidence that latent factors identify unknown proteins",
                                     "files":{"manifest":"assay_manifest.tsv","candidates":"candidate_loci.tsv","annotation":"protein_coordinates.GRCh38.tsv",
                                              "truth":"truth.tsv","effects":"truth_effects.tsv","cis_controls":"cis_controls.tsv"}})
    print(f"SYNTHETIC {args.scenario}: {m} source positions, {n} proteins, {len(manifest)} assays; {out}",flush=True)


if __name__=="__main__":
    try:generate(arguments())
    except ValueError as exc:raise SystemExit("ERROR: "+str(exc)) from exc
