"""Use the existing 1KG projection/ancestry R routines with GRCh38 imputed PGEN."""
from __future__ import annotations
from collections import Counter
import fcntl
from pathlib import Path
import subprocess
import pandas as pd
from core2 import atomic_table, atomic_json
from genomics import pfile, pvar_rows, bridge_for_chromosome, COMPLEMENT, allele_alignment

def project(cfg, root):
    # Multiple analysis roots can share one permanent projected-PC location.
    # Refuse concurrent writers instead of racing on R's output temporary files.
    folder=Path(cfg["inputs"]["pca_file"]).parent;folder.mkdir(parents=True,exist_ok=True)
    with (folder/"grid_projection.lock").open("a+") as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise RuntimeError("Another PCA projection is writing this permanent output directory")
        return _project(cfg,root)

def _project(cfg, root):
    inputs=cfg["inputs"];opts=cfg["pca"];root=Path(root)
    work=root/"tmp/pca";work.mkdir(parents=True,exist_ok=True)
    qc=root/"qc/pca";qc.mkdir(parents=True,exist_ok=True)
    pcfile=Path(inputs["pca_file"]);pcfile.parent.mkdir(parents=True,exist_ok=True)
    ancestry=Path(inputs["ancestry_file"]);ancestry.parent.mkdir(parents=True,exist_ok=True)
    if cfg.get("chromosomes",list(range(1,23)))!=list(range(1,23)):
        raise ValueError("Reference PCA projection requires all22autosomes; use a supplied projected PCA table for a chromosome-only smoke run")
    weights=pd.read_csv(inputs["pca_weight"],sep=r"\s+",dtype=str)
    ident=int(opts.get("id_column",2))-1;allele=int(opts.get("allele_column",6))-1
    first=int(opts.get("first_pc_column",7))-1;npc=int(opts.get("n_pcs",20))
    if weights.shape[1]<first+npc:raise ValueError("PCA weight file lacks requested PC columns")
    ref=weights.iloc[:,[ident,allele,*range(first,first+npc)]].copy()
    ref.columns=["SNP","A1",*[f"PC{i}" for i in range(1,npc+1)]]
    ref["SNP"]=ref.SNP.str.lower()
    if ref.SNP.duplicated().any():raise ValueError("PCA weights have duplicate IDs")
    pc_cols=[f"PC{i}" for i in range(1,npc+1)]
    names={str(c).lstrip("#").upper():c for c in weights}
    if {"REF","ALT"}<=set(names):pair_columns=(names["REF"],names["ALT"])
    elif {"A1","A2"}<=set(names):pair_columns=(names["A1"],names["A2"])
    elif opts.get("other_allele_column"):
        pair_columns=(weights.columns[allele],weights.columns[int(opts["other_allele_column"])-1])
    else:raise ValueError("PCA weight header needs REF/ALT or A1/A2; otherwise declare pca.other_allele_column explicitly")
    other=[]
    for effect,a,b in zip(ref.A1,weights[pair_columns[0]],weights[pair_columns[1]]):
        effect=str(effect).upper();pair={str(a).upper(),str(b).upper()}
        if effect not in pair:pair={x.translate(COMPLEMENT) for x in pair}
        if effect not in pair or len(pair)!=2:raise ValueError("PCA effect allele disagrees with its declared allele pair")
        other.append(next(x for x in pair if x!=effect))
    ref["A2"]=other
    ref=ref.set_index("SNP");wanted=set(ref.index)
    counts=[];total=0
    for chrom in range(1,23):
        prefix=Path(inputs["genotype_dir"])/f"chr{chrom}";arguments,pvar=pfile(prefix)
        bridge=bridge_for_chromosome(inputs.get("rsid_map"),chrom,wanted)
        found={};bad=set();target_count=Counter()
        for row in pvar_rows(pvar):
            if row["CHR"]!=chrom:raise ValueError("PVAR chromosome differs from its chr filename")
            ids=row["RSIDS"]&wanted
            for other,bp in bridge.get(row["TARGET_ID"],[]):
                if bp and int(bp)!=row["BP38"]:raise ValueError("PCA rsID bridge disagrees with PVAR BP38")
                ids|=other
            if row["TARGET_ID"] in {".", "0"} or row["BP38"] < 1:continue
            pair={row["REF38"],row["ALT38"]}
            if pair in ({"A","T"},{"C","G"}) or len(pair)!=2 or not pair<={"A","C","G","T"}:continue
            for sid in ids:
                orientation=allele_alignment(ref.loc[sid,"A1"],ref.loc[sid,"A2"],row["REF38"],row["ALT38"])
                if orientation is None:continue
                effect=orientation[0]
                if sid in found:bad.add(sid)
                found[sid]=[row["TARGET_ID"],effect,*ref.loc[sid,pc_cols].tolist()]
                target_count[row["TARGET_ID"]]+=1
        rows=[r for s,r in found.items() if s not in bad and target_count[r[0]]==1]
        # --extract selects every row sharing an ID, including incompatible
        # alleles/indels and records without the reference rsID in INFO. Count
        # all occurrences of selected native IDs, not only aligned records.
        # A second streaming pass keeps memory bounded by the PCA weight set.
        selected_ids={r[0] for r in rows}
        occurrences=Counter()
        for row in pvar_rows(pvar):
            if row["TARGET_ID"] in selected_ids:
                occurrences[row["TARGET_ID"]]+=1
        duplicated={ident for ident,n in occurrences.items() if n>1}
        rows=[r for r in rows if r[0] not in duplicated]
        atomic_json({"chromosome":chrom,"duplicate_target_ids":sorted(duplicated),
                     "policy":"exclude every occurrence of a repeated target ID"},
                    qc/f"chr{chrom}.duplicate_ids.json")
        if not rows:raise ValueError(f"No PCA variants on chr{chrom}; checkrsID mapping/reference alleles and duplicate IDs")
        weight_file=work/f"chr{chrom}.weights.tsv"
        atomic_table(pd.DataFrame(rows,columns=["ID","A1",*[f"PC{i}" for i in range(1,npc+1)]]),weight_file)
        extract=work/f"chr{chrom}.extract.txt"
        extract.write_text("\n".join(r[0] for r in rows)+"\n")
        out=work/f"chr{chrom}"
        command=[cfg["runtime"].get("plink","plink2"),*arguments,"--extract",str(extract),"--score",str(weight_file),"1","2","header-read","no-mean-imputation","list-variants","cols=nallele,scoresums",
                 "--score-col-nums",f"3-{npc+2}","--threads",str(cfg["runtime"]["threads"]),"--memory",str(cfg["runtime"].get("plink_memory_mb",4096)),"--out",str(out)]
        with (work/f"chr{chrom}.command.log").open("w") as log:
            log.write(repr(command)+"\n");log.flush()
            try:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
            except subprocess.CalledProcessError as exc:
                raise RuntimeError(f"PCA chr{chrom}: PLINK failed; see {work/f'chr{chrom}.command.log'}") from exc
        used_ids=[line.strip() for line in Path(str(out)+".sscore.vars").read_text().splitlines() if line.strip() and not line.startswith("#")]
        if len(set(used_ids))!=len(used_ids) or set(used_ids)!={r[0] for r in rows}:
            raise ValueError(f"chr{chrom}: PLINK's accepted PCA variant IDs differ from the aligned weights")
        used=len(used_ids)
        counts.append({"chromosome":chrom,"mapped_weights":len(rows),"used_weights":used,"duplicate_target_ids":len(duplicated)});total+=used
    coverage=total/len(ref)
    atomic_table(pd.DataFrame(counts),qc/"variants.tsv")
    atomic_json({"reference_variants":len(ref),"used_variants":total,"coverage":coverage,"genome_build":38,
                 "reference_space":"provided1KGallele-weights and matchingPCcenters; raw dosage SUM convention preserved",
                 "missing_dosage_rule":"no-mean-imputation, matching the original reference projection convention"},qc/"projection.json")
    if coverage<float(opts.get("min_variant_fraction",.9)):
        raise ValueError(f"Only{coverage:.1%}ofreferencePCAweightsmatched; inspectqc/pca/variants.tsv beforeusingincompatiblyscaledPCs")
    worker=Path(__file__).with_name("1.pca.R");rscript=cfg["runtime"].get("rscript","Rscript")
    commands=[
      [rscript,str(worker),"projection","--dir",str(work),"--out",str(pcfile),"--outdir",str(qc),"--n-pc",str(npc),"--distance-pcs",str(cfg["disco"].get("distance_pcs",10)),"--projection-only","TRUE"],
      [rscript,str(worker),"projection","--pca",str(pcfile),"--med",inputs["pca_centers"],"--outdir",str(qc),"--data-outdir",str(pcfile.parent),"--n-pc",str(npc),"--distance-pcs",str(cfg["disco"].get("distance_pcs",10)),"--score-dir",str(work)],
      [rscript,str(worker),"ancestry","--phe",inputs["phenotype_file"],"--pca",str(pcfile),"--med",inputs["pca_centers"],"--out",str(ancestry),"--outdir",str(qc),"--ethnicity-col",str(opts.get("ethnicity_column","auto")),"--n-pc",str(npc),"--distance-pcs",str(cfg["disco"].get("distance_pcs",10)),"--prob-threshold",str(opts.get("ancestry_probability",.90)),"--anchor-max-per-group",str(opts.get("anchor_max_n",10000))]
    ]
    for i,command in enumerate(commands):
        with (work/f"Rstep{i}.log").open("w") as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)

