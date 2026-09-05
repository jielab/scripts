#!/usr/bin/env python3
"""Fit baseline vs evolutionary transportability models with blocked out-of-fold prediction."""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
BASE=['global_pca_distance','mean_maf','delta_maf','mean_ldscore','delta_ldscore','log_sampling_var']
EVOL=['log_arg_age','mutation_branch_gen','root_time_gen','carrier_frequency','lineage_breadth','lineage_entropy','local_arg_divergence','local_global_ratio']

def prepare(train,valid,features):
 med=train[features].median(numeric_only=True); tr=train[features].fillna(med).to_numpy(float); va=valid[features].fillna(med).to_numpy(float)
 mu=np.nanmean(tr,0); sd=np.nanstd(tr,0); sd[~np.isfinite(sd)|(sd<1e-12)]=1
 return (tr-mu)/sd,(va-mu)/sd,med.to_dict(),mu,sd

def ridge(X,y,w,alpha):
 X=np.column_stack([np.ones(len(X)),X]); sw=np.sqrt(w); Xw=X*sw[:,None]; yw=y*sw
 pen=np.eye(X.shape[1])*alpha; pen[0,0]=0
 try:return np.linalg.solve(Xw.T@Xw+pen,Xw.T@yw)
 except np.linalg.LinAlgError:return np.linalg.pinv(Xw.T@Xw+pen)@(Xw.T@yw)
def pred(X,b):return b[0]+X@b[1:]
def metrics(y,p):
 ok=np.isfinite(y)&np.isfinite(p); y=y[ok];p=p[ok]
 rmse=float(np.sqrt(np.mean((y-p)**2))); mae=float(np.mean(abs(y-p))); r=float(spearmanr(y,p).statistic) if len(y)>2 else np.nan
 den=np.sum((y-y.mean())**2); r2=float(1-np.sum((y-p)**2)/den) if den>0 else np.nan
 return {'n':int(len(y)),'rmse':rmse,'mae':mae,'spearman':r,'r2':r2}
def fit_oof(d,features,groups,alpha):
 p=np.full(len(d),np.nan); rec=[]
 for g in sorted(pd.unique(groups)):
  te=np.asarray(groups==g); tr=~te
  if tr.sum()<max(50,3*len(features)) or te.sum()<10: continue
  Xtr,Xte,_,_,_=prepare(d.loc[tr],d.loc[te],features); y=d['transport_heterogeneity'].to_numpy(float); w=d['model_weight'].to_numpy(float)
  b=ridge(Xtr,y[tr],w[tr],alpha); p[te]=pred(Xte,b); rec.append({'validation_group':str(g),**metrics(y[te],p[te])})
 return p,rec
def full_fit(d,features,alpha):
 X,_,med,mu,sd=prepare(d,d,features); y=d.transport_heterogeneity.to_numpy(float);w=d.model_weight.to_numpy(float); b=ridge(X,y,w,alpha)
 return b,med,mu,sd
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--input',required=True); ap.add_argument('--out-dir',required=True); ap.add_argument('--ridge-alpha',type=float,default=10); ap.add_argument('--conservation-min',type=float,default=.05); ap.add_argument('--conservation-max',type=float,default=.995); ap.add_argument('--seed',type=int,default=20260904); a=ap.parse_args()
 out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
 use=['trait','chr','bp','SNP','A1','A2','pop_a','pop_b','transport_heterogeneity','sampling_var']+BASE+EVOL
 d=pd.read_csv(a.input,sep='\t',compression='infer',usecols=lambda x:x in use,low_memory=False,dtype={'chr':str,'SNP':str})
 d['transport_heterogeneity']=pd.to_numeric(d.transport_heterogeneity,errors='coerce');d=d[np.isfinite(d.transport_heterogeneity)].reset_index(drop=True)
 if len(d)<200: raise SystemExit(f'Too few transport observations: {len(d)}')
 # Winsorize only the training outcome tail; preserve rank and zero mass.
 hi=d.transport_heterogeneity.quantile(.995); d.transport_heterogeneity=d.transport_heterogeneity.clip(upper=hi)
 sv=pd.to_numeric(d.sampling_var,errors='coerce'); w=1/np.sqrt(sv.clip(lower=sv[sv>0].quantile(.01) if (sv>0).any() else 1)); lo,hiw=w.quantile([.01,.99]); d['model_weight']=w.clip(lo,hiw).fillna(1); d.model_weight/=d.model_weight.mean()
 # Remove unusable columns; require at least one ARG-specific predictor.
 bfeat=[x for x in BASE if x in d and d[x].notna().mean()>.01]
 efeat=[x for x in EVOL if x in d and d[x].notna().mean()>.01]
 if len(bfeat)<3: raise SystemExit(f'Baseline MAF/LD/global-distance features incomplete: {bfeat}')
 if not any(x in efeat for x in ['log_arg_age','local_arg_divergence','root_time_gen']): raise SystemExit(f'ARG age/local genealogy unavailable: {efeat}')
 ffeat=bfeat+efeat
 chroms=d.chr.astype(str).nunique()
 if chroms>=3:
  groups=d.chr.astype(str).to_numpy(); validation='leave_one_chromosome_out'
 else:
  bp=pd.to_numeric(d.bp,errors='coerce').fillna(np.arange(len(d))); groups=((bp//5_000_000).astype(int)%5).astype(str).to_numpy(); validation='five_genomic_block_folds_pilot'
 pb,rb=fit_oof(d,bfeat,groups,a.ridge_alpha); pf,rf=fit_oof(d,ffeat,groups,a.ridge_alpha)
 y=d.transport_heterogeneity.to_numpy(float); mb=metrics(y,pb); mf=metrics(y,pf)
 selected='evolutionary_full' if mf['rmse']<mb['rmse'] else 'baseline'
 ps=pf if selected=='evolutionary_full' else pb
 d['pred_baseline_oof']=pb; d['pred_full_oof']=pf; d['pred_selected_oof']=ps
 # Variant prior is based exclusively on held-out predictions.
 v=(d.groupby(['trait','chr','bp','SNP','A1','A2'],dropna=False,as_index=False)
      .agg(predicted_heterogeneity=('pred_selected_oof','median'),observed_heterogeneity=('transport_heterogeneity','median'),n_pairs=('transport_heterogeneity','size')))
 ph=np.maximum(pd.to_numeric(v.predicted_heterogeneity,errors='coerce').fillna(0),0);v['conservation']=np.exp(-.5*ph).clip(a.conservation_min,a.conservation_max);v['selected_model']=selected
 v.to_csv(out/'variant_conservation.tsv.gz',sep='\t',index=False,compression='gzip',float_format='%.9g')
 d.to_csv(out/'pair_predictions.tsv.gz',sep='\t',index=False,compression='gzip',float_format='%.9g')
 cv=[]
 for name,rec in [('baseline',rb),('evolutionary_full',rf)]:
  for x in rec:cv.append({'model':name,'validation':validation,**x})
 cv += [{'model':'baseline','validation':'all_oof',**mb},{'model':'evolutionary_full','validation':'all_oof',**mf}]
 pd.DataFrame(cv).to_csv(out/'model_cv.tsv',sep='\t',index=False)
 models={}
 coeff=[]
 for name,features in [('baseline',bfeat),('evolutionary_full',ffeat)]:
  b,med,mu,sd=full_fit(d,features,a.ridge_alpha);models[name]={'features':features,'intercept':float(b[0]),'coefficients':[float(x) for x in b[1:]],'medians':{k:float(v) for k,v in med.items()},'means':[float(x) for x in mu],'scales':[float(x) for x in sd]}
  coeff.append(pd.DataFrame({'model':name,'feature':['intercept']+features,'standardized_coefficient':b}))
 pd.concat(coeff).to_csv(out/'coefficients.tsv',sep='\t',index=False)
 summary={'n_rows':len(d),'n_variants':len(v),'validation':validation,'baseline_features':bfeat,'evolutionary_features':efeat,'ridge_alpha':a.ridge_alpha,'baseline_oof':mb,'evolutionary_oof':mf,'selected_model':selected,'rmse_improvement':mb['rmse']-mf['rmse']}
 (out/'transport_model.json').write_text(json.dumps({'summary':summary,'models':models},indent=2)+'\n')
 print(json.dumps(summary))
if __name__=='__main__':main()
