#!/usr/bin/env python3
"""Apply a trusted frozen panome run to a new, unlabeled CSV cohort; no refitting."""
import argparse
import json
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
import torch
from representation import Autoencoder,encode


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-dir',type=Path,required=True)
    ap.add_argument('--phenotype',required=True);ap.add_argument('--omics',required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();root=args.run_dir
    cfg=json.loads((root/'manifest.json').read_text())['config'];idcol=cfg['id_col']
    p=pd.read_csv(args.phenotype,dtype={idcol:str});raw=pd.read_csv(args.omics,dtype={idcol:str})
    for table in [p,raw]:
        if idcol not in table or table[idcol].isna().any() or table[idcol].duplicated().any():raise ValueError('Missing / duplicate person IDs')
    if not raw[idcol].isin(p[idcol]).all():raise ValueError('Every omics person must have phenotype metadata')
    p=p.set_index(idcol).loc[raw[idcol]].reset_index()
    pre=joblib.load(root/'s2_preprocess/preprocessor.joblib')
    original=(root/'s1_prepare/features.txt').read_text().splitlines()
    needed=np.asarray(original)[pre.keep]
    missing=set(needed)-set(raw)
    if missing:raise ValueError(f'Missing required features: {sorted(missing)[:20]}; harmonize assay features explicitly')
    x=np.full((len(raw),len(original)),np.nan,dtype='float32');x[:,pre.keep]=raw[list(needed)].to_numpy(dtype='float32')
    if (np.mean(~np.isfinite(x[:,pre.keep]),axis=1)>=cfg['sample_missing']).any():raise ValueError('Some new people fail training-defined sample missingness threshold')
    x,_=pre.transform(x,p)
    checkpoint=torch.load(root/'s3_representation/autoencoder.pt',map_location='cpu',weights_only=True)
    torch.set_num_threads(cfg['cores']);model=Autoencoder(checkpoint['p'],checkpoint['latent'],checkpoint['hidden'])
    model.load_state_dict(checkpoint['state_dict']);model.eval();z=encode(model,x,'cpu')
    gm=joblib.load(root/'s4_graph/graph_model.joblib');zz=gm['scaler'].transform(z).astype('float32');index=gm['index'];k=gm['neighbors']
    if hasattr(index,'query'):ids,dist=index.query(zz,k=k)
    else:dist,ids=index.kneighbors(zz,n_neighbors=k)
    weights=np.exp(-dist/np.maximum(dist[:,-1:],1e-8));weights/=weights.sum(axis=1,keepdims=True)
    soft=np.stack([(weights*(gm['labels'][ids]==s)).sum(axis=1) for s in range(gm['labels'].max()+1)],axis=1)
    state=soft.argmax(axis=1);diffusion=(weights[:,:,None]*gm['psi'][ids]).sum(axis=1)
    cp=joblib.load(root/'s5_predict/clinical_preprocessor.joblib');c=cp['transformer'].transform(p)[:,cp['nonconstant']]
    pca=joblib.load(root/'s3_representation/pca.joblib').transform(x)
    pw=joblib.load(root/'s5_predict/pwas_score.joblib');score=x[:,pw['indices']]@pw['beta']
    hard=np.eye(soft.shape[1])[state][:,1:]
    blocks={'clinical':c,'clinical_pwas':np.c_[c,score],'clinical_pca':np.c_[c,pca],
            'clinical_ae':np.c_[c,z],'clinical_state':np.c_[c,hard],'clinical_soft_state':np.c_[c,soft[:,1:]],
            'clinical_diffusion':np.c_[c,diffusion],'clinical_panome':np.c_[c,z,soft[:,1:]],'clinical_elasticnet':np.c_[c,x]}
    result=p[[idcol]].copy();result['state']=state+1;result['novel']=dist[:,0]>gm['novelty_threshold']
    for j in range(soft.shape[1]):result[f'state_weight_{j+1}']=soft[:,j]
    for j in range(z.shape[1]):result[f'AE{j+1}']=z[:,j]
    for name,block in blocks.items():
        artifact=joblib.load(root/f's5_predict/{name}.joblib');b=artifact['scaler'].transform(block)[:,artifact['active']]
        kwargs={'alpha':artifact['alpha']} if name=='clinical_elasticnet' else {}
        fitted=artifact['model'];result[name+'_log_hazard']=fitted.predict(b,**kwargs)
        survival=fitted.predict_survival_function(b,**kwargs)
        for horizon in cfg['horizons']:
            if horizon<=survival[0].domain[1]:result[f'{name}_net_risk_{horizon:g}y']=[1-fn(horizon) for fn in survival]
    if args.output.exists():raise FileExistsError(f'Refusing to overwrite {args.output}')
    args.output.parent.mkdir(parents=True,exist_ok=True);result.to_csv(args.output,index=False)
    print(f'Projected {len(result)} people with frozen training models: {args.output}')

if __name__=='__main__':main()
