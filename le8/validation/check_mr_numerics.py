"""Check real MR regional optimizer starts and optional CUDA eigh from audit files."""
import hashlib,importlib.util,json,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
root=Path(__file__).resolve().parents[1]
def load(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
m=load('mr_numeric_check',root/'f/c2.mr_link2.py');r=load('resources_check',root/'f/0.resources.py')
w=Path(sys.argv[1]);rows=[]
for file in sorted((w/'serial-cache').rglob('*_numerical_audit.json')):
 audit=json.loads(file.read_text());trait=file.name.removeprefix('protein_').removesuffix('_numerical_audit.json')
 if not audit.get('optimizations'):continue
 matrices=np.load(w/'eigen'/audit['eigh']['cache_key']/'eigen.npz');lam=matrices['values'];u=matrices['vectors'];ld=(u*lam)@u.T
 def betas(kind):
  f=file.parent.parent/'prepared'/f'protein_{trait}.{kind}.ref.tsv';data=pd.read_csv(f,sep='\t').set_index('SNP');answer=[]
  for snp,a1,a2 in audit['ordered_alleles']:
   v=data.loc[snp];z=v.BETA/v.SE;b=z/np.sqrt(v.N+z*z);pair=(v.EA,v.NEA);comp=tuple(v.translate(str.maketrans('ACGT','TGCA')) for v in pair)
   if pair==(a1,a2) or comp==(a1,a2):answer.append(b)
   elif pair==(a2,a1) or comp==(a2,a1):answer.append(-b)
   else:raise ValueError('Unexpected allele mismatch')
  return np.array(answer)
 x,y=betas('exposure'),betas('outcome');nX,nY=audit['N_exposure'],audit['N_outcome'];guess=audit['sigma_initial']
 threshold=sorted(lam[lam>0],reverse=True)[np.argmin(np.cumsum(sorted(lam[lam>0],reverse=True))/sum(lam[lam>0])<=.99)]
 take=lam>=threshold
 for factor in [1.,.5,2.]:
  result=m.mr_link2(lam[take],u[:,take],x,y,nX,nY,guess[0]*factor,guess[1]*factor)
  rows.append(dict(trait=trait,backend='cpu',initial_factor=factor,**{k:(v.item() if hasattr(v,'item') else v) for k,v in result.items()}))
  if factor==1:np.testing.assert_allclose(result['alpha'],audit['optimizations'][0]['result']['alpha'],atol=1e-6,rtol=0)
 start=time.monotonic();gl,gu,diag=r.validated_eigh(ld,'cuda');seconds=time.monotonic()-start
 gthreshold=sorted(gl[gl>0],reverse=True)[np.argmin(np.cumsum(sorted(gl[gl>0],reverse=True))/sum(gl[gl>0])<=.99)];gt=gl>=gthreshold
 result=m.mr_link2(gl[gt],gu[:,gt],x,y,nX,nY,*guess)
 rows.append(dict(trait=trait,backend='cuda_eigh',eigh_seconds=seconds,initial_factor=1.,eigen_diagnostic=diag,
   **{k:(v.item() if hasattr(v,'item') else v) for k,v in result.items()}))
 np.testing.assert_allclose(result['alpha'],audit['optimizations'][0]['result']['alpha'],atol=1e-6,rtol=0)
(w/'numerical-starts-cuda.json').write_text(json.dumps(rows,indent=2))
print('PASS',len(rows),'actual regional fits; 3 CPU starts + CUDA float64 eigen backend per trait')
for trait in sorted({v['trait'] for v in rows}):
 z=[v for v in rows if v['trait']==trait];print(trait,'alpha_range',max(v['alpha'] for v in z)-min(v['alpha'] for v in z),'ha_loglik_range',max(v['ha_loglik'] for v in z)-min(v['ha_loglik'] for v in z),'all_converged',all(v[k] for v in z for k in ['optim_ha_success','optim_alpha_h0_success','optim_sigma_y_h0_success']))
