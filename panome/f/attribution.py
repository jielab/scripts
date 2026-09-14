"""Person-specific integrated gradients of the AE-Cox log hazard (not causal effects)."""
import joblib
import numpy as np
import pandas as pd
import torch
from representation import Autoencoder


def attribute(x,p,a,root,out):
    if a.attribution_samples==0:return
    torch.set_num_threads(a.cores)
    checkpoint=torch.load(root/'s3_representation/autoencoder.pt',map_location='cpu',weights_only=True)
    model=Autoencoder(checkpoint['p'],checkpoint['latent'],checkpoint['hidden'])
    model.load_state_dict(checkpoint['state_dict']);model.eval()
    fit=joblib.load(root/'s5_predict/clinical_ae.joblib')
    coefficients=np.zeros(len(fit['active']));coefficients[fit['active']]=fit['model'].coef_
    beta=torch.as_tensor((coefficients/fit['scaler'].scale_)[-a.latent:],dtype=torch.float32)
    rng=np.random.default_rng(a.seed);ix=np.flatnonzero(p.split.values=='test')
    ix=np.sort(rng.choice(ix,min(a.attribution_samples,len(ix)),replace=False))
    features=(root/'s2_preprocess/features.txt').read_text().splitlines();rows=[];completeness=[]
    for start in range(0,len(ix),64):
        subset=ix[start:start+64];target=torch.as_tensor(np.array(x[subset]),dtype=torch.float32)
        total=torch.zeros_like(target)
        for step in range(a.ig_steps):
            point=(target*((step+.5)/a.ig_steps)).requires_grad_(True)
            value=(model.encoder(point)*beta).sum()
            total+=torch.autograd.grad(value,point)[0].detach()/a.ig_steps
        contribution=(target*total).numpy()
        with torch.no_grad():delta=((model.encoder(target)-model.encoder(torch.zeros_like(target)))*beta).sum(axis=1).numpy()
        for local,index in enumerate(subset):
            selected=np.argsort(-np.abs(contribution[local]))[:20]
            for j in selected:rows.append(dict(person_id=p.iloc[index][a.id_col],feature=features[j],integrated_gradient=contribution[local,j]))
            completeness.append(dict(person_id=p.iloc[index][a.id_col],log_hazard_delta=float(delta[local]),
                                     attribution_sum=float(contribution[local].sum()),approximation_error=float(contribution[local].sum()-delta[local])))
    pd.DataFrame(rows).to_csv(out/'person_ae_cox_attributions.csv',index=False)
    pd.DataFrame(completeness).to_csv(out/'attribution_completeness.csv',index=False)
