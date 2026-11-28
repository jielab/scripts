"""Decode a single tag exposure from the matching catalogue release."""
import argparse
import gzip
import json
import os
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd


def export(report,locus,dest):
    root=report.resolve()
    meta=json.loads((root/"release.json").read_text())
    cat=pd.read_csv(root/"loci.tsv.gz",sep="\t",dtype={"chrom":str},usecols=["id","lineage","chrom","start0","YGB_carriers","YGB_uncertain"])
    found=cat[cat.id==locus]
    if len(found)!=1:
        raise ValueError(f"Unknown ASV in this release: {locus}")
    row=found.iloc[0]
    t=int(row.start0)//meta["window_bp"]
    dest.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=dest.parent,delete=False) as tmp:
        temp=Path(tmp.name)
    opener = gzip.open if dest.suffix==".gz" else open
    yes_total=uncertain_total=0
    try:
        with opener(temp,"wt") as out:
            out.write("IID\tbatch\t"+locus+"\n")
            for cohort in meta["cohorts"][1:]:
                name=cohort["name"]
                folder=root/"genotypes"/name
                all_ids=np.load(folder/"autosome.samples.npy",allow_pickle=False)
                current=(np.load(folder/"X.samples.npy",allow_pickle=False) if row.chrom=="X" else all_ids)
                with np.load(folder/f"chr{row.chrom}.npz") as file:
                    yes,unknown=[np.unpackbits(file[row.lineage+"_"+key][t],count=len(current),bitorder="little").astype(bool)
                                 for key in ("carrier","uncertain")]
                if (yes & unknown).any():
                    raise ValueError("Overlapping carrier and uncertain states")
                yes_total+=yes.sum()
                uncertain_total+=unknown.sum()
                values=pd.Series(np.select([yes,unknown],["1","NA"],default="0"),index=current).reindex(all_ids,fill_value="NA")
                pd.DataFrame({"IID":all_ids,"batch":name,locus:values.to_numpy()}).to_csv(out,sep="\t",index=False,header=False)
        if (yes_total,uncertain_total)!=(row.YGB_carriers,row.YGB_uncertain):
            raise ValueError("Genotype counts disagree with this catalogue release")
        os.replace(temp,dest)
    finally:
        temp.unlink(missing_ok=True)


def export_matrix(report,dest):
    """Write tag × person int8 exposures, including X-ineligible people as missing."""
    import h5py

    root=report.resolve()
    meta=json.loads((root/"release.json").read_text())
    cat=pd.read_csv(root/"loci.tsv.gz",sep="\t",dtype={"chrom":str})
    cohorts=[c["name"] for c in meta["cohorts"][1:]]
    folders=[root/"genotypes"/c for c in cohorts]
    ids=[np.load(f/"autosome.samples.npy",allow_pickle=False) for f in folders]
    n=sum(map(len,ids))
    if cat.empty:
        raise ValueError("No catalogue markers for matrix export")
    dest.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=dest.parent,delete=False) as tmp:
        temp=Path(tmp.name)
    checks=np.zeros((len(cat),2),np.int64)
    try:
        with h5py.File(temp,"w") as h:
            h.attrs["definition"]="Rows=ASV sentinel; columns=IID; carrier=1, no qualifying lineage match=0, ambiguous or untested=-1"
            h.attrs["build"]=meta["build"]
            h.attrs["source"]=meta["source"]
            h.create_dataset("ASV",data=np.array(cat.id,dtype="S"))
            h.create_dataset("IID",data=np.concatenate(ids).astype("S"))
            h.create_dataset("batch",data=np.repeat(cohorts,[len(i) for i in ids]).astype("S"))
            matrix=h.create_dataset("exposure",(len(cat),n),dtype="i1",fillvalue=-1,
                                    chunks=(min(64,len(cat)),min(2048,n)),compression="gzip",compression_opts=1)
            offset=0
            for folder,people in zip(folders,ids):
                for (chrom,label),block in cat.groupby(["chrom","lineage"],sort=False):
                    x=chrom=="X"
                    tested=np.load(folder/"X.samples.npy",allow_pickle=False) if x else people
                    positions=pd.Index(people).get_indexer(tested)
                    if (positions<0).any():
                        raise ValueError("Missing sample in chromosome roster")
                    with np.load(folder/f"chr{chrom}.npz") as archive:
                        yes=archive[label+"_carrier"]
                        uncertain=archive[label+"_uncertain"]
                    for start in range(0,len(block),64):
                        chunk=block.iloc[start:start+64]
                        rows=chunk.index.to_numpy()
                        tiles=chunk.start0.to_numpy()//meta["window_bp"]
                        a=np.unpackbits(yes[tiles],axis=1,count=len(tested),bitorder="little").astype(np.int8)
                        u=np.unpackbits(uncertain[tiles],axis=1,count=len(tested),bitorder="little").astype(np.int8)
                        checks[rows,0]+=a.sum(axis=1)
                        checks[rows,1]+=u.sum(axis=1)
                        values=np.full((len(chunk),len(people)),-1,np.int8)
                        values[:,positions]=a-u
                        matrix[rows,offset:offset+len(people)]=values
                offset+=len(people)
            if not np.array_equal(checks,cat[["YGB_carriers","YGB_uncertain"]].to_numpy()):
                raise ValueError("Exposure matrix disagrees with catalogue counts")
        os.replace(temp,dest)
    finally:
        temp.unlink(missing_ok=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--report",required=True,type=Path)
    select=p.add_mutually_exclusive_group(required=True)
    select.add_argument("--locus")
    select.add_argument("--matrix",action="store_true",help="HDF5 ASV × IID matrix, states 1, 0 and -1 (missing)")
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    if a.matrix:
        export_matrix(a.report,a.out)
    else:
        export(a.report,a.locus,a.out)
