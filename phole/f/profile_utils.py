"""Shared positive scaling and projection of frozen protein effect profiles."""
from __future__ import annotations
import numpy as np
from graph_model import training_scales

def aggregate(values,target):
    names=sorted(target.PROT.unique()); out=np.full((len(values),len(names)),np.nan)
    for j,name in enumerate(names):
        x=np.asarray(values[:,target.PROT.eq(name).to_numpy()],float)
        count=np.isfinite(x).sum(axis=1)
        out[:,j]=np.divide(np.nansum(x,axis=1),count,out=np.full(len(x),np.nan),where=count>0)
    return out,names

def normalise_cohort(beta,se,observed,cis,source,target,excluded_blocks):
    mask=observed&~cis&np.isfinite(beta)&np.isfinite(se)&(se>0)
    training=mask&~source.LD_BLOCK.astype(str).isin(set(map(str,excluded_blocks))).to_numpy()[:,None]
    scales=training_scales(beta,se,training); centres=np.zeros(beta.shape[1]); usable=training.sum(axis=0)>=10
    for j in range(beta.shape[1]):
        if usable[j]: centres[j]=float(np.median(np.asarray(beta)[training[:,j],j]/scales[j]))
    values=np.where(mask&usable[None,:],np.asarray(beta)/scales[None,:]-centres[None,:],np.nan)
    errors=np.where(mask&usable[None,:],np.asarray(se)/scales[None,:],np.nan)
    grouped,names=aggregate(values,target); error_grouped,_=aggregate(errors,target)
    # Mean SE is conservative for positively correlated repeated probes. It is
    # used for a screening threshold, never as an independent-assay meta-analysis.
    return dict(values=grouped,errors=error_grouped,proteins=names,scales=scales,centres=centres,usable=usable)

def project_profile(values,basis,mask=None):
    values=np.asarray(values,float); basis=np.asarray(basis,float)
    good=np.isfinite(values)&np.all(np.isfinite(basis),axis=1)
    if mask is not None: good &= np.asarray(mask,bool)
    out=np.full(values.shape,np.nan)
    if good.sum()<3: return out
    x=basis[good]
    # Orthogonal projection is invariant to rotation of the frozen basis.
    coef,_,rank,_=np.linalg.lstsq(x,values[good],rcond=1e-8)
    if int(good.sum())-int(rank)<3: return out
    out[good]=values[good]-x@coef
    return out

def grouped_basis(v,target,names):
    return np.column_stack([np.ones(len(names)),np.stack([np.mean(v[target.PROT.eq(n).to_numpy()],axis=0) for n in names])])

def grouped_geometry(features,target,names):
    """Freeze primary coordinate features; these carry no measured external effect."""
    result=[]
    for feature in features:
        grouped,order=aggregate(feature,target)
        if list(order)!=list(names): raise ValueError('Geometry and protein order differ')
        result.append(grouped.astype(np.float32))
    return np.stack(result,axis=2)
