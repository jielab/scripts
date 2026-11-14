#!/usr/bin/env python3
"""Masked, signed graph prediction with grouped nested holdouts.

Only observed, allele-harmonized beta/SE enter the objective. Genomic distance
is a competing predictor, never assumed to be a physical law. Matrix completion
is explicitly distinguished from prediction of a completely unseen protein.
"""
from __future__ import annotations
import argparse, hashlib, json, os, re, time
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

def atom_json(path,value):
    tmp=Path(str(path)+'.tmp'); tmp.write_text(json.dumps(value,indent=2,allow_nan=False,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))); tmp.replace(path)
def atom_npy(path,value):
    tmp=Path(str(path)+'.tmp.npy'); np.save(tmp,value); tmp.replace(path)
def frame(path,rows):
    pd.DataFrame(rows).to_csv(path,sep='\t',index=False,na_rep='.')
def log(s): print(time.strftime('%H:%M:%S')+' '+s,flush=True)
def checksum(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def token(value): return hashlib.sha256(json.dumps(value,sort_keys=True,default=str).encode()).hexdigest()
def load_matrix(folder):
    folder=Path(folder); meta=json.loads((folder/'manifest.json').read_text())
    if meta.get('schema',meta.get('format'))!='phole-signed-matrix-v2' or not meta.get('directional') or meta.get('input_scope')!='full':
        raise ValueError('Signed v2 full-GWAS matrices are required; MAGMA and selected .sig leads are not signed input')
    source=pd.read_csv(folder/'source.tsv',sep='\t',dtype={'CHR':str,'LD_BLOCK':str})
    target=pd.read_csv(folder/'target.tsv',sep='\t',dtype={'PROT':str,'ASSAY':str,'CHR':str})
    beta=np.load(folder/'beta.npy',mmap_mode='r'); se=np.load(folder/'se.npy',mmap_mode='r')
    mask=np.load(folder/'observed.npy',mmap_mode='r').astype(bool)
    cis=np.load(folder/'cis.npy',mmap_mode='r').astype(bool)
    if beta.shape!=(len(source),len(target)) or any(a.shape!=beta.shape for a in (se,mask,cis)):
        raise ValueError('Matrix and row/column metadata shapes differ')
    mask &= np.isfinite(beta)&np.isfinite(se)&(se>0)
    known_path=folder/'target_cis_known.npy'
    if known_path.is_file():
        known=np.load(known_path).astype(bool)
        if known.shape!=(len(target),): raise ValueError('target_cis_known shape mismatch')
        mask &= known[None,:]
    elif 'ANNOTATION_STATUS' in target:
        mask &= target.ANNOTATION_STATUS.eq('single_gene').to_numpy()[None,:]
    return beta,se,mask,cis,source,target,meta

def target_groups(target):
    """Keep known gene aliases/proteoforms together without merging columns.

    PROT_GROUP is annotation supplied by prepare, never inferred from effects or
    from nearby coordinates. Legacy inputs retain their exact PROT grouping.
    If only one repeated assay has a group label, propagate that same label to
    the other assays of its canonical protein; conflicting labels are an error.
    """
    proteins=target.PROT.astype(str)
    if proteins.str.strip().str.upper().isin(['','.','NAN','NONE','NA']).any():
        raise ValueError('Target PROT identifiers must be present for grouped evaluation')
    if 'PROT_GROUP' not in target: return proteins.rename('PROT_GROUP')
    supplied=target.PROT_GROUP.map(lambda value: None if pd.isna(value) or str(value).strip().upper() in ('','.','NAN','NONE','NA') else str(value).strip())
    lookup={}
    for protein in proteins.unique():
        labels=set(supplied[proteins.eq(protein)].dropna())
        if len(labels)>1: raise ValueError(f'Conflicting PROT_GROUP labels for {protein}')
        lookup[protein]=next(iter(labels)) if labels else protein
    return proteins.map(lookup).rename('PROT_GROUP')

def target_group_weights(target):
    groups=target_groups(target)
    return 1/groups.map(groups.value_counts()).to_numpy(float)

def grouped_folds(source,target,k,seed):
    def assignment(values,salt):
        groups=sorted(set(map(str,values))); rng=np.random.default_rng(seed+salt); rng.shuffle(groups)
        lookup={g:i%k for i,g in enumerate(groups)}
        return np.array([lookup[str(g)] for g in values],dtype=np.int8)
    groups=target_groups(target)
    if source.LD_BLOCK.nunique()<k or groups.nunique()<k:
        raise ValueError('Fewer LD/block or protein groups than folds; reduce --folds')
    rf=assignment(source.LD_BLOCK,17); cf=assignment(groups,43)
    return (rf[:,None]+cf[None,:])%k,rf,cf

def inner_group_holdout(source,target,seed,fraction=.2):
    rows,ru=pd.factorize(source.LD_BLOCK.astype(str),sort=True)
    cols,cu=pd.factorize(target_groups(target),sort=True)
    draw=np.random.default_rng(seed).random((len(ru),len(cu)))<fraction
    return draw[rows[:,None],cols[None,:]]

def target_intervals(row):
    """Read explicit interval triples or unambiguous parallel input columns."""
    explicit=row.get('TARGET_INTERVALS')
    if explicit is not None and not pd.isna(explicit) and str(explicit).strip() not in ('','.','nan'):
        values=[]
        for part in str(explicit).split(';'):
            match=re.fullmatch(r'([^:;]+):(\d+)-(\d+)',part.strip())
            if not match: return []
            values.append((match[1],float(match[2]),float(match[3])))
    else:
        chromosomes=str(row.get('CHR','.')).split(';')
        starts=str(row.get('START','.')).split(';'); ends=str(row.get('END','.')).split(';')
        if len(starts)!=len(ends): return []
        if len(chromosomes)==1: chromosomes*=len(starts)
        if len(chromosomes)!=len(starts): return []
        try: values=[(c,float(s),float(e)) for c,s,e in zip(chromosomes,starts,ends)]
        except (TypeError,ValueError): return []
    primary={*(str(i) for i in range(1,23)),'X','Y','MT'}
    if any(c not in primary or not np.isfinite(s) or not np.isfinite(e) or s<1 or e<s for c,s,e in values):
        return []
    return values

def geometry(source,target):
    ch=source.CHR.astype(str).to_numpy(); pos=source.POS.to_numpy(float)
    tc=target.CHR.astype(str).to_numpy()
    st=pd.to_numeric(target.START,errors='coerce').to_numpy(float); en=pd.to_numeric(target.END,errors='coerce').to_numpy(float)
    known=np.isfinite(st)&np.isfinite(en)&(st<=en)
    same=(ch[:,None]==tc[None,:])&known[None,:]
    dist=np.maximum(np.maximum(st[None,:]-pos[:,None],pos[:,None]-en[None,:]),0)
    dist=np.where(same,dist,np.inf)
    # Preserve the scalar path and numeric behavior, but use the nearest
    # same-chromosome interval for PAR genes and complete multi-interval targets.
    # CHR/START/END cannot reconstruct arbitrary deduplicated chromosome lists;
    # prepare therefore writes explicit triples in TARGET_INTERVALS.
    for j,(_,row) in enumerate(target.iterrows()):
        explicit=row.get('TARGET_INTERVALS')
        has_explicit=explicit is not None and not pd.isna(explicit) and str(explicit).strip() not in ('','.','nan')
        if not has_explicit and not any(';' in str(row.get(key,'')) for key in ('CHR','START','END')): continue
        same[:,j]=False; dist[:,j]=np.inf
        for chromosome,start,end in target_intervals(row):
            hit=ch==chromosome
            delta=np.maximum(np.maximum(start-pos,pos-end),0)
            same[:,j] |= hit
            dist[:,j]=np.minimum(dist[:,j],np.where(hit,delta,np.inf))
    return [same.astype(np.float32)]+[(same*np.exp(-dist/b)).astype(np.float32) for b in (1e6,1e7,1e8)]

def training_scales(beta,se,mask):
    # Do not inspect a held-out effect when estimating a column's scale.
    scale=np.ones(beta.shape[1])
    for j in range(beta.shape[1]):
        v=np.asarray(beta[mask[:,j],j],float); e=np.asarray(se[mask[:,j],j],float)
        if len(v)<3: continue
        scale[j]=max(float(np.median(np.abs(v)))*1.4826,float(np.median(e)),1e-12)
    return scale

def fit_model(beta,se,mask,features,rank=0,use_geometry=False,iterations=30,ridge=5.,seed=1,column_weights=None):
    n,m=beta.shape; scale=training_scales(beta,se,mask)
    y=np.where(mask,np.asarray(beta)/scale[None,:],0.).astype(np.float64)
    e=np.asarray(se)/scale[None,:]
    # Weight bound prevents a small number of extremely precise cells dominating.
    w=np.where(mask,np.clip(1/np.maximum(e,1e-12)**2,0.05,20),0.)
    if column_weights is not None: w*=np.asarray(column_weights)[None,:]
    w[~np.isfinite(w)]=0
    eye=np.eye(rank)
    weight_sum=max(float(w.sum()),1.)

    def one_start(start_seed):
        row=np.zeros(n); col=np.zeros(m); theta=np.zeros(len(features))
        rng=np.random.default_rng(start_seed)
        u=rng.normal(0,.02,(n,rank)); v=rng.normal(0,.02,(m,rank))
        g=np.zeros((n,m)); lr=u@v.T if rank else np.zeros((n,m)); old=None; converged=False
        for it in range(iterations):
            row=np.sum(w*(y-col[None,:]-g-lr),axis=1)/(w.sum(axis=1)+ridge)
            col=np.sum(w*(y-row[:,None]-g-lr),axis=0)/(w.sum(axis=0)+ridge)
            if use_geometry:
                residual=y-row[:,None]-col[None,:]-lr
                for k,x in enumerate(features):
                    residual_k=residual-g+theta[k]*x
                    new=np.sum(w*x*residual_k)/(np.sum(w*x*x)+ridge)
                    g+=(new-theta[k])*x; theta[k]=new
            if rank:
                z=y-row[:,None]-col[None,:]-g
                for i in range(n):
                    keep=w[i]>0
                    if keep.sum()<rank: u[i]=0; continue
                    vv=v[keep]; wi=w[i,keep]
                    u[i]=np.linalg.solve((vv.T*wi)@vv+ridge*eye,vv.T@(wi*z[i,keep]))
                for j in range(m):
                    keep=w[:,j]>0
                    if keep.sum()<rank: v[j]=0; continue
                    uu=u[keep]; wj=w[keep,j]
                    v[j]=np.linalg.solve((uu.T*wj)@uu+ridge*eye,uu.T@(wj*z[keep,j]))
                lr=u@v.T
            pred=row[:,None]+col[None,:]+g+lr
            data_loss=float(np.sum(w*(y-pred)**2)/weight_sum)
            penalty=float(ridge*(np.sum(row**2)+np.sum(col**2)+np.sum(theta**2)+
                                np.sum(u**2)+np.sum(v**2))/weight_sum)
            objective=data_loss+penalty
            if old is not None and abs(old-objective)<=1e-6*max(1.,abs(old)):
                converged=True; break
            old=objective
        return dict(row=row,col=col,theta=theta,u=u,v=v,scale=scale,loss=objective,
                    data_loss=data_loss,ridge_penalty=penalty,iterations=it+1,converged=converged,
                    rank=rank,use_geometry=use_geometry,selected_seed=start_seed)

    # Low-rank ALS is nonconvex. Choose among deterministic starts using only
    # the same penalized training objective used for convergence; never a test cell.
    seeds=[int(seed)+offset for offset in range(3 if rank else 1)]
    fits=[one_start(start_seed) for start_seed in seeds]
    best=min(fits,key=lambda model:model['loss'])
    best['start_seeds']=np.array(seeds,dtype=np.int64)
    best['start_objectives']=np.array([model['loss'] for model in fits])
    return best

def predict(model,features):
    p=model['row'][:,None]+model['col'][None,:]
    if model['use_geometry']:
        for t,x in zip(model['theta'],features): p=p+t*x
    if model['rank']: p=p+model['u']@model['v'].T
    return p*model['scale'][None,:]

def loss_on(beta,se,pred,mask,scale,column_weights=None):
    valid=mask&np.isfinite(pred)
    if not np.any(valid): return float('inf')
    weights=np.clip(1/np.maximum(np.asarray(se)/scale[None,:],1e-12)**2,.05,20)
    if column_weights is not None: weights*=np.asarray(column_weights)[None,:]
    z=(np.asarray(beta)-pred)/scale[None,:]
    return float(np.sum(weights[valid]*z[valid]**2)/np.sum(weights[valid]))

def metrics(y,p,se):
    keep=np.isfinite(y)&np.isfinite(p)&np.isfinite(se)&(se>0); y=y[keep]; p=p[keep]; se=se[keep]
    if len(y)<3: return dict(n=len(y),rmse=None,r2=None,pearson=None,sign_agreement=None,rmse_se=None)
    total=np.sum((y-y.mean())**2); good=np.abs(y/se)>=3
    return dict(n=len(y),rmse=float(np.sqrt(np.mean((y-p)**2))),r2=float(1-np.sum((y-p)**2)/total) if total else None,
       pearson=float(np.corrcoef(y,p)[0,1]) if np.std(y)>0 and np.std(p)>0 else None,
       sign_agreement=float(np.mean(np.sign(y[good])==np.sign(p[good]))) if good.any() else None,
       rmse_se=float(np.sqrt(np.mean(np.minimum(((y-p)/se)**2,1e6)))))

def block_bootstrap_gain(beta,se,pred,baseline,mask,source,target,seed,repeats=200):
    """Two-way source-block/target-group resampling; conditional on fixed panel."""
    rb,ru=pd.factorize(source.LD_BLOCK.astype(str),sort=True)
    cb,cu=pd.factorize(target_groups(target),sort=True)
    weights=target_group_weights(target)
    # Use the actual held-out loss: truncation can hide a small number of
    # disastrous predictions and falsely turn a worse model into an improvement.
    y=np.asarray(beta,dtype=float); uncertainty=np.asarray(se,dtype=float)
    base_error=np.zeros(y.shape); selected_error=np.zeros(y.shape)
    base_error[mask]=(y[mask]-np.asarray(baseline,dtype=float)[mask])**2
    selected_error[mask]=(y[mask]-np.asarray(pred,dtype=float)[mask])**2
    difference=np.zeros(y.shape)
    difference[mask]=base_error[mask]/uncertainty[mask]**2-selected_error[mask]/uncertainty[mask]**2
    if not np.isfinite(difference[mask]).all():
        raise ValueError('Nonfinite held-out loss; predictive-structure evidence cannot be computed')
    value=difference*weights[None,:]
    count=mask*weights[None,:]
    d=np.zeros((len(ru),beta.shape[1])); n=np.zeros_like(d)
    np.add.at(d,rb,value); np.add.at(n,rb,count)
    dg=np.zeros((len(ru),len(cu))); ng=np.zeros_like(dg)
    np.add.at(dg.T,cb,d.T); np.add.at(ng.T,cb,n.T)
    observed=float(dg.sum()/max(ng.sum(),1)); rng=np.random.default_rng(seed); samples=[]
    for _ in range(repeats):
        r=rng.multinomial(len(ru),np.full(len(ru),1/len(ru)))
        c=rng.multinomial(len(cu),np.full(len(cu),1/len(cu)))
        denominator=float(r@ng@c)
        if denominator: samples.append(float(r@dg@c)/denominator)
    low,high=np.quantile(samples,[.025,.975])
    denominator=max(float(count.sum()),1.)
    base_rmse=float(np.sqrt(np.sum(base_error*weights[None,:])/denominator))
    selected_rmse=float(np.sqrt(np.sum(selected_error*weights[None,:])/denominator))
    raw_base_rmse=float(np.sqrt(np.mean(base_error[mask])))
    raw_selected_rmse=float(np.sqrt(np.mean(selected_error[mask])))
    raw_improved=selected_rmse<base_rmse and raw_selected_rmse<raw_base_rmse
    return dict(comparison='selected minus hotspot improvement',loss_reduction=observed,ci_low=float(low),ci_high=float(high),
                predictive_structure_supported=bool(low>0 and observed>0 and raw_improved),bootstrap_repeats=len(samples),
                selected_raw_beta_rmse=raw_selected_rmse,hotspot_raw_beta_rmse=raw_base_rmse,
                selected_protein_weighted_beta_rmse=selected_rmse,hotspot_protein_weighted_beta_rmse=base_rmse,
                raw_beta_rmse_improved=bool(raw_improved),loss='untruncated squared beta/SE residual, inverse assay count per target group',
                n_source_blocks=len(ru),n_target_groups=len(cu),
                scope='two-way source-block/target-group bootstrap; internal completion conditional on discovery panel')

def fit_pipeline(matrix_dir,out,folds=5,ranks=(0,2,4,8),iterations=30,ridge=5.,seed=20261105,replace=False,fixed_rank=None,fixed_geometry=None):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    beta,se,observed,cis,source,target,meta=load_matrix(matrix_dir)
    groups=target_groups(target)
    if len(source)<10 or groups.nunique()<10: raise ValueError('At least 10 loci and 10 distinct protein/gene groups are required')
    valid=observed&~cis
    if valid.sum()<100: raise ValueError('Too few observed trans edges')
    input_files=['manifest.json','beta.npy','se.npy','observed.npy','cis.npy','source.tsv','target.tsv']
    known_path=Path(matrix_dir)/'target_cis_known.npy'
    if known_path.is_file(): input_files.append('target_cis_known.npy')
    input_hashes={name:checksum(Path(matrix_dir)/name) for name in input_files}
    signature=token(dict(code=checksum(__file__),input_hashes=input_hashes,folds=folds,ranks=ranks,
        iterations=iterations,ridge=ridge,seed=seed,fixed_rank=fixed_rank,fixed_geometry=fixed_geometry))
    mf=out/'fit.manifest.json'
    if mf.exists() and not replace:
        saved=json.loads(mf.read_text())
        required=['prediction.oof.npy','baseline.oof.npy','distance.oof.npy','residual.oof.npy','folds.npy','model.npz',
                  'metrics.tsv','selection.tsv','fold-metrics.tsv','prediction-gain.tsv','target-splits.tsv','source-splits.tsv']
        complete=saved.get('signature')==signature and all((out/f).is_file() for f in required)
        if complete and saved.get('model_hash')==checksum(out/'model.npz') and saved.get('residual_hash')==checksum(out/'residual.oof.npy'):
            log('Reuse matching completed fit '+str(out)); return saved
        raise ValueError('Existing fit differs or is incomplete; use --replace TRUE or a new output')
    cellfold,rf,cf=grouped_folds(source,target,folds,seed); feats=geometry(source,target)
    column_weights=target_group_weights(target)
    prediction=np.full(beta.shape,np.nan,dtype=np.float32); baseline=np.full(beta.shape,np.nan,dtype=np.float32)
    distance=np.full(beta.shape,np.nan,dtype=np.float32); selection=[]; fold_metrics=[]; selected_ranks=[]; selected_geometry=[]
    grid=[(0,False),(0,True)]+[(r,True) for r in sorted(set(ranks)) if r>0 and r<min(beta.shape)-2]
    checkpoint_dir=out/'checkpoints'; checkpoint_dir.mkdir(exist_ok=True)
    for f in range(folds):
        outer=valid&(cellfold==f); train=valid&~outer
        if not outer.any(): continue
        checkpoint=checkpoint_dir/f'fold-{f}.npz'
        if checkpoint.is_file() and not replace:
            saved=np.load(checkpoint)
            if str(saved['signature'])==signature and len(saved['prediction'])==int(outer.sum()):
                prediction[outer]=saved['prediction']; baseline[outer]=saved['baseline']; distance[outer]=saved['distance']
                selected_ranks.append(int(saved['rank'])); selected_geometry.append(bool(saved['geometry']))
                selection.extend(json.loads(str(saved['selection'])))
                for name,pp in [('zero',np.zeros(int(outer.sum()))),('hotspot',saved['baseline']),('distance',saved['distance']),('selected',saved['prediction'])]:
                    fold_metrics.append(dict(fold=f,model=name,**metrics(np.asarray(beta)[outer],pp,np.asarray(se)[outer])))
                log(f'Reused completed fold {f+1}/{folds}'); continue
        inner=train&inner_group_holdout(source,target,seed+211+f); inner_train=train&~inner
        if fixed_rank is None:
            choices=[]
            for r,g in grid:
                mod=fit_model(beta,se,inner_train,feats,r,g,max(8,iterations//2),ridge,seed+100*f,column_weights)
                val=loss_on(beta,se,predict(mod,feats),inner,mod['scale'],column_weights)
                choices.append((val,r,g)); selection.append(dict(fold=f,rank=r,geometry=g,inner_loss=val))
            _,rank,use_geom=min(choices,key=lambda x:(x[0],x[1],x[2]))
        else: rank,use_geom=int(fixed_rank),bool(fixed_geometry) if fixed_geometry is not None else True
        selected_ranks.append(rank); selected_geometry.append(use_geom); log(f'fold {f+1}/{folds}: rank={rank} geometry={use_geom}')
        # All benchmark models are trained on the identical outer training cells.
        hotspot=fit_model(beta,se,train,feats,0,False,iterations,ridge,seed+f,column_weights)
        base_pred=predict(hotspot,feats); baseline[outer]=base_pred[outer]
        geo=fit_model(beta,se,train,feats,0,True,iterations,ridge,seed+f,column_weights)
        geo_pred=predict(geo,feats); distance[outer]=geo_pred[outer]
        chosen=hotspot if rank==0 and not use_geom else geo if rank==0 else fit_model(beta,se,train,feats,rank,use_geom,iterations,ridge,seed+f,column_weights)
        pr=predict(chosen,feats); prediction[outer]=pr[outer]
        for name,pp in [('zero',np.zeros(beta.shape)),('hotspot',base_pred),('distance',geo_pred),('selected',pr)]:
            fold_metrics.append(dict(fold=f,model=name,**metrics(np.asarray(beta)[outer],pp[outer],np.asarray(se)[outer])))
        tmp=Path(str(checkpoint)+'.tmp.npz')
        np.savez_compressed(tmp,signature=signature,prediction=prediction[outer],baseline=baseline[outer],distance=distance[outer],
                            rank=rank,geometry=use_geom,selection=json.dumps([r for r in selection if r.get('fold')==f]))
        tmp.replace(checkpoint)
        # Atomic fold checkpoint; later rerun can reuse only a full matching fit.
        atom_npy(out/'prediction.oof.partial.npy',prediction)
    completed=valid&np.isfinite(prediction)
    pooled=[]
    for name,pp in [('zero',np.zeros(beta.shape)),('hotspot',baseline),('distance',distance),('selected',prediction)]:
        pooled.append(dict(scope='grouped_out_of_fold_trans',model=name,**metrics(np.asarray(beta)[completed],pp[completed],np.asarray(se)[completed])))
    family_counts=Counter(zip(selected_ranks,selected_geometry))
    final_rank,final_geometry=min(family_counts,key=lambda x:(-family_counts[x],x[0],x[1]))
    final=fit_model(beta,se,valid,feats,final_rank,final_geometry,iterations,ridge,seed+999,column_weights)
    atom_npy(out/'prediction.oof.npy',prediction); atom_npy(out/'baseline.oof.npy',baseline); atom_npy(out/'distance.oof.npy',distance)
    atom_npy(out/'residual.oof.npy',np.where(completed,np.asarray(beta)-prediction,np.nan).astype(np.float32))
    atom_npy(out/'folds.npy',cellfold)
    np.savez_compressed(out/'model.npz',**final,row_fold=rf,target_fold=cf)
    frame(out/'metrics.tsv',pooled); frame(out/'fold-metrics.tsv',fold_metrics); frame(out/'selection.tsv',selection or [dict(rank=final_rank,frozen=True)])
    gain=block_bootstrap_gain(beta,se,prediction,baseline,completed,source,target,seed+1701)
    frame(out/'prediction-gain.tsv',[gain])
    frame(out/'target-splits.tsv',target[['PROT','ASSAY']].assign(PROT_GROUP=groups,FOLD=cf).to_dict('records'))
    frame(out/'source-splits.tsv',source[['LOCUS','CHR','POS','LD_BLOCK']].assign(FOLD=rf).to_dict('records'))
    result=dict(schema='phole-signed-fit-v2',signature=signature,matrix_dir=str(Path(matrix_dir).resolve()),rank=final_rank,
       input_hashes=input_hashes,model_hash=checksum(out/'model.npz'),residual_hash=checksum(out/'residual.oof.npy'),
       rank_grid=[int(r) for r in ranks],
       n_target_proteins=int(target.PROT.nunique()),n_target_groups=int(groups.nunique()),
       target_group_definition='Provided PROT_GROUP (known gene aliases/proteoforms); PROT grouping when annotation group labels are absent; separate assay columns retained',
       selected_ranks=selected_ranks,folds=folds,iterations=iterations,ridge=ridge,seed=seed,n_observed_trans=int(completed.sum()),
       evaluation='matrix completion with entire source-block x protein-group cells held out; not cold-start protein prediction',
       panel_selection=meta.get('panel_selection','primary-selected loci; internal metrics conditioned on the panel'),
       scaling='training cells only, per fold',loss='masked inverse-variance weighted, weights clipped; cis excluded from trans loss',
       optimizer='three deterministic ALS starts for rank>0; minimize training loss plus ridge penalty on row, column, geometry and factors',
       causal_claim=False,backend='numpy CPU',external_rank_frozen=fixed_rank is not None,final_geometry=final_geometry,
       synthetic=bool(meta.get('synthetic',False)))
    result['predictive_structure']=gain
    atom_json(mf,result); (out/'fit.done').write_text(signature+'\n')
    (out/'prediction.oof.partial.npy').unlink(missing_ok=True)
    return result

def fit_main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--matrix-dir',required=True); p.add_argument('--out',required=True)
    p.add_argument('--folds',type=int,default=5); p.add_argument('--rank-grid',default='0,2,4,8'); p.add_argument('--iterations',type=int,default=30)
    p.add_argument('--ridge',type=float,default=5.); p.add_argument('--seed',type=int,default=20261105); p.add_argument('--replace',default='TRUE')
    p.add_argument('--fixed-rank',type=int); a=p.parse_args()
    if a.folds<3 or a.iterations<2 or a.ridge<=0: p.error('folds>=3, iterations>=2, ridge>0 required')
    fit_pipeline(a.matrix_dir,a.out,a.folds,tuple(map(int,a.rank_grid.split(','))),a.iterations,a.ridge,a.seed,a.replace.upper()=='TRUE',a.fixed_rank)
