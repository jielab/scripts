#!/usr/bin/env python3
"""Build DISCOVERY centres in the same PC coordinates as the target projection.

Accepts projected discovery participants, or precomputed discovery PC centres.
For linear dosage-sum projections (the current 0pca.sh), discovery EAF can
also produce mean PC coordinates without access to discovery genotypes.
"""
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=['individuals','eaf'],required=True)
    p.add_argument('--input',required=True,help='Individuals: POP,PC1...; EAF: manifest POP,FILE,N_GWAS')
    p.add_argument('--output',required=True);p.add_argument('--pcs',type=int,default=10)
    p.add_argument('--pca-space',required=True,help='Coordinate-system identifier shared with Yeval --pca-space')
    p.add_argument('--source',required=True,help='Discovery cohort/projection provenance')
    p.add_argument('--pca-weights',help='Current 0pca dosage-sum weight table: SNP col2, scored allele col6, PCs col7 onward')
    p.add_argument('--projection-snps',help='Exact SNP set used in target projection (one ID per line)')
    p.add_argument('--frequency-allele',required=False,help='EAF mode: name of effect-allele column, typically A1')
    a=p.parse_args();pc=[f'PC{i+1}' for i in range(a.pcs)]
    if a.pcs<2:raise ValueError('Need at least two PCs')
    d=pd.read_csv(a.input,sep='\t',dtype={'POP':str})
    if a.mode=='individuals':
        if not set(['POP']+pc)<=set(d):raise ValueError('Discovery projections need POP and all PCs')
        if not np.isfinite(d[pc].to_numpy(float)).all() or d.POP.isna().any():raise ValueError('Nonfinite PCs or missing POP')
        out=d.groupby('POP',sort=False)[pc].mean().reset_index()
        out['N_GWAS']=out.POP.map(d.groupby('POP').size())
    else:
        if not a.pca_weights or not a.projection_snps or not a.frequency_allele:
            raise ValueError('EAF mode requires --pca-weights, --projection-snps and --frequency-allele')
        if not {'POP','FILE','N_GWAS'}<=set(d) or d.POP.duplicated().any():raise ValueError('Invalid EAF manifest')
        w=pd.read_csv(a.pca_weights,sep=r'\s+')
        wanted=Path(a.projection_snps).read_text().split()
        if len(wanted)!=len(set(wanted)):raise ValueError('Duplicate projection SNPs')
        w=w.set_index(w.columns[1],drop=False).reindex(wanted)
        if w.iloc[:,5:].isna().any().any() or w.index.duplicated().any():raise ValueError('Missing/duplicate projection weights')
        weights=w.iloc[:,6:6+a.pcs].to_numpy(float)
        if weights.shape[1]!=a.pcs or not np.isfinite(weights).all():raise ValueError('Invalid PC weights')
        rows=[]
        for row in d.itertuples(index=False):
            path=Path(row.FILE)
            if not path.is_absolute():path=Path(a.input).resolve().parent/path
            f=pd.read_csv(path,sep='\t',dtype={'SNP':str})
            if f.SNP.duplicated().any():raise ValueError('Duplicate EAF SNPs')
            f=f.set_index('SNP').reindex(wanted)
            # Strict allele agreement avoids silent ambiguity and requires frequencies
            # already harmonized to the EXACT target projection effect allele.
            if not np.array_equal(f[a.frequency_allele].to_numpy(),w.iloc[:,5].to_numpy()):
                raise ValueError('EAF alleles must match target PCA weight alleles exactly; harmonize first')
            af=pd.to_numeric(f.EAF,errors='coerce').to_numpy()
            if not np.isfinite(af).all() or np.any((af<0)|(af>1)):raise ValueError('All projection SNPs require discovery EAF; no reference-frequency filling')
            mean=(2*af)@weights
            rows.append(dict(POP=row.POP,N_GWAS=row.N_GWAS,**dict(zip(pc,mean))))
        out=pd.DataFrame(rows)
    if out.POP.duplicated().any() or not np.isfinite(out.N_GWAS).all() or (out.N_GWAS<=0).any():raise ValueError('Invalid group sizes')
    out['pca_space']=a.pca_space;out['source']=a.source;out['kind']='discovery'
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(dest,sep='\t',index=False)
    Path(str(dest)+'.json').write_text(json.dumps(vars(a),indent=2)+'\n')
    print(f'DONE discovery centres: {dest}')

if __name__=='__main__':main()
