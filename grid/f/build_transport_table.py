#!/usr/bin/env python3
"""Build pairwise cross-ancestry effect-transportability observations."""
from __future__ import annotations
import argparse,gzip,itertools,math,re
from pathlib import Path
import numpy as np,pandas as pd
POPS0=["AFR","EAS","EUR","SAS"]
COMP=str.maketrans("ACGT","TGCA")
def comp(x): return str(x).translate(COMP)
def read(path): return pd.read_csv(path,sep='\t',compression='infer',low_memory=False,dtype={'SNP':str,'CHR':str})
def centers(path,pops,npc=10):
 d=pd.read_csv(path,sep=None,engine='python',compression='infer')
 pc=sorted([c for c in d.columns if re.fullmatch(r'PC\d+',str(c),re.I)],key=lambda x:int(re.findall(r'\d+',x)[0]))[:npc]
 pcoll=next((c for c in ['pop','POP','super_pop','ancestry','population','race'] if c in d.columns),None)
 if pcoll is None or len(pc)<2: raise SystemExit(f'Cannot read population PC centers from {path}: {list(d.columns)}')
 d[pcoll]=d[pcoll].astype(str).str.upper(); out={}
 for p in pops:
  x=d[d[pcoll]==p]
  if len(x): out[p]=x[pc].apply(pd.to_numeric,errors='coerce').median().to_numpy(float)
 if len(out)<2: raise SystemExit(f'Fewer than two population centers in {path}')
 return out
def align(d,ref):
 x=d.merge(ref,on='SNP',how='inner',suffixes=('','_REF'))
 a1=x.A1.astype(str).str.upper();a2=x.A2.astype(str).str.upper();r1=x.A1_REF.astype(str).str.upper();r2=x.A2_REF.astype(str).str.upper()
 same=(a1==r1)&(a2==r2); swap=(a1==r2)&(a2==r1); cs=(a1.map(comp)==r1)&(a2.map(comp)==r2); cw=(a1.map(comp)==r2)&(a2.map(comp)==r1)
 ok=same|swap|cs|cw; x=x[ok].copy(); flip=(swap|cw)[ok].to_numpy()
 x['BETA_ALIGNED']=pd.to_numeric(x.BETA,errors='coerce').to_numpy()*np.where(flip,-1,1)
 e=pd.to_numeric(x.get('EAF'),errors='coerce'); x['EAF_ALIGNED']=np.where(flip,1-e,e)
 return x

def external_age(path):
 if not path:return None
 d=pd.read_csv(path,sep=None,engine='python',compression='infer',low_memory=False)
 sc=next((c for c in ['SNP','rsid','RSID','id','variant_id'] if c in d.columns),None)
 ac=next((c for c in d.columns if any(z in str(c).lower() for z in ['age_gen','allele_age','geva_age','mean_age','age'])),None)
 if sc is None or ac is None: raise SystemExit(f'Cannot find SNP/age in {path}')
 age=pd.to_numeric(d[ac],errors='coerce'); unit='generations'
 if 'year' in ac.lower(): age=age/29.0; unit='years_to_generations'
 return pd.DataFrame({'SNP':d[sc].astype(str),'external_age_gen':age}).drop_duplicates('SNP'),unit

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--trait',required=True); ap.add_argument('--pops',default='AFR,EAS,EUR,SAS'); ap.add_argument('--chrs',required=True); ap.add_argument('--sumstats-dir',required=True); ap.add_argument('--arg-dir',required=True); ap.add_argument('--ldscore-dir',required=True); ap.add_argument('--centers',required=True); ap.add_argument('--external-age',default=''); ap.add_argument('--max-snps-per-chr',type=int,default=0); ap.add_argument('--out',required=True); a=ap.parse_args()
 pops=[x.upper() for x in re.split('[,; ]+',a.pops) if x]; chrs=[x for x in re.split('[,; ]+',a.chrs) if x]; cen=centers(a.centers,pops); ext=external_age(a.external_age); extdf=ext[0] if ext else None
 Path(a.out).parent.mkdir(parents=True,exist_ok=True); first=True; total=0
 for chrom in chrs:
  ds={}
  for p in pops:
   f=Path(a.sumstats_dir)/f'{a.trait}.{p}.chr{chrom}.tsv.gz'
   if not f.exists(): continue
   d=read(f); d['SNP']=d.SNP.astype(str); d['A1']=d.A1.astype(str).str.upper();d['A2']=d.A2.astype(str).str.upper(); d['CHR']=str(chrom); d['BP']=pd.to_numeric(d.BP,errors='coerce'); d=d.dropna(subset=['SNP','A1','A2','BETA','SE']).drop_duplicates('SNP')
   ds[p]=d
  if len(ds)<2: continue
  priority=[p for p in ['EUR','AFR','EAS','SAS'] if p in ds]
  long=pd.concat([ds[p][['SNP','A1','A2','BP']].assign(_prio=i) for i,p in enumerate(priority)],ignore_index=True).sort_values('_prio')
  ref=long.drop_duplicates('SNP')[['SNP','A1','A2','BP']]
  wide=ref.rename(columns={'A1':'REF_A1','A2':'REF_A2','BP':'REF_BP'}).copy()
  refalign=ref.rename(columns={'A1':'A1_REF','A2':'A2_REF','BP':'BP_REF'})
  for p,d in ds.items():
   z=align(d,refalign)
   z=z[['SNP','BETA_ALIGNED','SE','EAF_ALIGNED','P','N','BP']].rename(columns={c:f'{c}_{p}' for c in ['BETA_ALIGNED','SE','EAF_ALIGNED','P','N','BP']})
   wide=wide.merge(z,on='SNP',how='left')
  if a.max_snps_per_chr>0: wide=wide.sort_values('REF_BP').head(a.max_snps_per_chr)
  # ARG features: first rsID, then position fallback.
  afile=Path(a.arg_dir)/f'chr{chrom}.variants.tsv.gz'; arg=read(afile) if afile.exists() else pd.DataFrame()
  if len(arg):
   arg['SNP']=arg.SNP.astype(str); arg['bp']=pd.to_numeric(arg.bp,errors='coerce'); arg=arg.drop_duplicates('SNP')
   keep=[c for c in arg.columns if c not in {'chr','allele0','allele1'}]
   wide=wide.merge(arg[keep],on='SNP',how='left')
   miss=wide['arg_age_gen'].isna() if 'arg_age_gen' in wide else pd.Series(True,index=wide.index)
   if miss.any():
    apos=arg.drop_duplicates('bp').set_index('bp'); ix=pd.to_numeric(wide.loc[miss,'REF_BP'],errors='coerce');
    for c in [x for x in arg.columns if x not in {'chr','SNP','bp','allele0','allele1'}]:
      if c not in wide: wide[c]=np.nan
      wide.loc[miss,c]=ix.map(apos[c])
  if extdf is not None: wide=wide.merge(extdf,on='SNP',how='left')
  # LD score joins.
  for p in pops:
   lf=Path(a.ldscore_dir)/p/f'chr{chrom}.ldscore.tsv.gz'
   if lf.exists():
    ld=read(lf)[['SNP','ldscore']].drop_duplicates('SNP').rename(columns={'ldscore':f'ldscore_{p}'})
    wide=wide.merge(ld,on='SNP',how='left')
  rows=[]
  for pa,pb in itertools.combinations(pops,2):
   ba=f'BETA_ALIGNED_{pa}';bb=f'BETA_ALIGNED_{pb}';sa=f'SE_{pa}';sb=f'SE_{pb}'
   if not all(c in wide for c in [ba,bb,sa,sb]): continue
   x=wide[wide[ba].notna()&wide[bb].notna()&wide[sa].notna()&wide[sb].notna()].copy();
   if not len(x):continue
   va=pd.to_numeric(x[sa],errors='coerce')**2; vb=pd.to_numeric(x[sb],errors='coerce')**2; svar=va+vb; q=(pd.to_numeric(x[ba])-pd.to_numeric(x[bb]))**2/svar.replace(0,np.nan); y=np.log1p(np.maximum(q-1,0))
   ea=pd.to_numeric(x.get(f'EAF_ALIGNED_{pa}'),errors='coerce'); eb=pd.to_numeric(x.get(f'EAF_ALIGNED_{pb}'),errors='coerce')
   # ARG frequencies are orientation-free after conversion to MAF.
   if f'af_{pa}' in x: ea=ea.fillna(pd.to_numeric(x[f'af_{pa}'],errors='coerce'))
   if f'af_{pb}' in x: eb=eb.fillna(pd.to_numeric(x[f'af_{pb}'],errors='coerce'))
   ma=np.minimum(ea,1-ea); mb=np.minimum(eb,1-eb)
   la=pd.to_numeric(x.get(f'ldscore_{pa}'),errors='coerce'); lb=pd.to_numeric(x.get(f'ldscore_{pb}'),errors='coerce')
   dc=float(np.linalg.norm(cen[pa]-cen[pb])) if pa in cen and pb in cen else np.nan
   dcol=f'div_{pa}_{pb}' if f'div_{pa}_{pb}' in x else f'div_{pb}_{pa}'
   local=pd.to_numeric(x.get(dcol),errors='coerce')
   age=pd.to_numeric(x.get('arg_age_gen'),errors='coerce')
   if 'external_age_gen' in x: age=age.fillna(pd.to_numeric(x.external_age_gen,errors='coerce'))
   z=pd.DataFrame({'trait':a.trait,'chr':str(chrom),'bp':x.REF_BP,'SNP':x.SNP,'A1':x.REF_A1,'A2':x.REF_A2,'pop_a':pa,'pop_b':pb,'beta_a':x[ba],'beta_b':x[bb],'se_a':x[sa],'se_b':x[sb],'transport_heterogeneity':y,'sampling_var':svar,'log_sampling_var':np.log(svar),'maf_a':ma,'maf_b':mb,'mean_maf':(ma+mb)/2,'delta_maf':abs(ma-mb),'ldscore_a':la,'ldscore_b':lb,'mean_ldscore':(la+lb)/2,'delta_ldscore':abs(la-lb),'global_pca_distance':dc,'local_arg_divergence':local,'local_global_ratio':local/(dc+1e-8),'arg_age_gen':age,'log_arg_age':np.log1p(age)})
   for c in ['mutation_branch_gen','root_time_gen','carrier_frequency','lineage_breadth','lineage_entropy']:
    z[c]=pd.to_numeric(x.get(c),errors='coerce')
   rows.append(z)
  if rows:
   z=pd.concat(rows,ignore_index=True); z.to_csv(a.out,sep='\t',index=False,compression='gzip',mode='wt' if first else 'at',header=first,float_format='%.9g'); first=False; total+=len(z)
 if first: raise SystemExit('No pairwise transport observations were generated')
 print(f'rows={total} output={a.out}')
if __name__=='__main__':main()
