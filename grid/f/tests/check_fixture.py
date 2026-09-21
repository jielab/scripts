"""Independent numerical check, including deliberately shuffled sample IDs."""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
r=Path(sys.argv[1])
pca=pd.read_csv(r/'pca.tsv.gz',sep='\t').set_index('IID').sort_index()
med=pd.read_csv(r/'centers.tsv',sep='\t').iloc[:,1:].to_numpy()
D=np.linalg.norm(med[:,None,:]-med[None,:,:],axis=2)
a=np.linalg.solve(D,np.ones(4))
dis=np.linalg.norm(pca.to_numpy()[:,None,:]-med[None,:,:],axis=2)
q=a/dis; q/=q.sum(axis=1)[:,None]
X=np.column_stack([np.ones(len(pca)),pca.to_numpy()]); z=[]
for pop in ['AFR','EAS','EUR','SAS']:
    d=pd.read_csv(r/'scores/height/chr22/csx.pgs.gz',sep='\t').set_index('eid').reindex(pca.index)[f'csx.{pop}'].to_numpy()
    res=d-X@np.linalg.lstsq(X,d,rcond=None)[0]
    z.append((res-res.mean())/res.std(ddof=1))
    w=pd.read_csv(r/f'gwas/height.{pop}/gwas/height.{pop}.chr22.csx.gz',sep='\t')
    assert len(w)==120 and w.SNP.is_unique and np.isfinite(w.BETA).all()
    assert set(w.CHR)=={22} and np.isfinite(w.BP).all() and w.BP.min()==100001
expected=(np.array(z).T*q).sum(axis=1)
out=pd.read_csv(r/'scores/height/chr22/disco.pgs.gz',sep='\t').set_index('IID').reindex(pca.index).disco
assert np.allclose(expected,out,atol=1e-9),np.max(abs(expected-out))
print('PASS: official Disco scores match independent distance interpolation and PCA residualization')
