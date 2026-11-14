#!/usr/bin/env python3
"""Freeze unexplained locus profiles as hypotheses for external replication."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from graph_model import load_matrix, atom_json, checksum, frame, token, fit_model, geometry, target_group_weights
from profile_utils import normalise_cohort, project_profile, grouped_basis, grouped_geometry

def aggregate_proteins(values,target):
    groups=sorted(target.PROT.unique()); out=np.full((len(values),len(groups)),np.nan)
    for j,p in enumerate(groups):
        a=values[:,target.PROT.eq(p).to_numpy()]
        count=np.isfinite(a).sum(axis=1)
        # Repeated probes are one biological target, not independent confirmations.
        out[:,j]=np.divide(np.nansum(a,axis=1),count,out=np.full(len(a),np.nan),where=count>0)
    return out,groups

def discover(matrix_dir,model_dir,out,min_targets=10,min_abs_z=3.,max_candidates=50):
    out=Path(out); out.mkdir(parents=True,exist_ok=True); model_dir=Path(model_dir)
    beta,se,observed,cis,source,target,meta=load_matrix(matrix_dir)
    fit_meta=json.loads((model_dir/'fit.manifest.json').read_text())
    if not fit_meta.get('input_hashes'):
        raise ValueError('Fit has no input provenance; rerun the current fit implementation')
    for name,expected in fit_meta['input_hashes'].items():
        if checksum(Path(matrix_dir)/name)!=expected: raise ValueError(f'Model/matrix provenance mismatch: {name}')
    for name,key in [('model.npz','model_hash'),('residual.oof.npy','residual_hash')]:
        if not fit_meta.get(key) or checksum(model_dir/name)!=fit_meta[key]:
            raise ValueError(f'Model output differs from fit receipt: {name}')
    residual=np.load(model_dir/'residual.oof.npy'); keep=observed&~cis&np.isfinite(residual)
    z=np.where(keep,residual/np.maximum(np.asarray(se),1e-30),np.nan)
    profiles,prots=aggregate_proteins(z,target)
    rows=[]
    for i,r in source.iterrows():
        v=profiles[i]; good=np.isfinite(v); n=int(good.sum()); large=int(np.sum(np.abs(v[good])>=min_abs_z))
        near=target[(target.CHR.astype(str)==str(r.CHR))&(pd.to_numeric(target.START,errors='coerce')-500000<=r.POS)&
                    (pd.to_numeric(target.END,errors='coerce')+500000>=r.POS)].GENE.dropna().unique() if 'GENE' in target else []
        rows.append(dict(LOCUS=r.LOCUS,SOURCE_INDEX=i,CHR=str(r.CHR),POS=int(r.POS),LD_BLOCK=str(r.LD_BLOCK),N_PROTEINS=n,
            N_LARGE_RESIDUAL=large,RESIDUAL_RMS_Z=float(np.sqrt(np.mean(np.minimum(v[good]**2,400)))) if n else None,
            RESIDUAL_MEDIAN_ABS_Z=float(np.median(np.abs(v[good]))) if n else None,
            NEARBY_MEASURED_GENES=';'.join(sorted(map(str,near))) or '.',PROGRAM_KIND='residual_locus',
            STATUS='descriptive_only',rank=0))
    all_rows=pd.DataFrame(rows).sort_values(['N_LARGE_RESIDUAL','RESIDUAL_RMS_Z'],ascending=False)
    eligible=all_rows[(all_rows.N_PROTEINS>=min_targets)&(all_rows.N_LARGE_RESIDUAL>=min_targets)]
    # One hypothesis per block reduces repeated nomination of LD proxies.
    eligible=eligible.drop_duplicates('LD_BLOCK').head(max_candidates).copy()
    eligible.insert(0,'CANDIDATE',[f'BH{i+1:04d}' for i in range(len(eligible))])
    eligible['rank']=np.arange(1,len(eligible)+1); eligible['STATUS']='discovery_candidate'
    # Re-estimate the ordinary background WITHOUT any nominated source block.
    # This prevents a candidate programme from defining the subspace against
    # which its own cross-cohort residual is evaluated.
    excluded=eligible.LD_BLOCK.astype(str).unique().tolist()
    normal=normalise_cohort(beta,se,observed,cis,source,target,excluded)
    training=observed&~cis&~source.LD_BLOCK.astype(str).isin(excluded).to_numpy()[:,None]
    training &= normal['usable'][None,:]
    features=geometry(source,target)
    if len(eligible):
        clean=fit_model(beta,se,training,features,rank=fit_meta['rank'],
             use_geometry=fit_meta['final_geometry'],iterations=fit_meta['iterations'],ridge=fit_meta['ridge'],
             seed=fit_meta['seed']+431, column_weights=target_group_weights(target))
    else:
        saved=np.load(model_dir/'model.npz'); clean={'v':saved['v']}
    names=normal['proteins']; basis=grouped_basis(clean['v'],target,names)
    geo=grouped_geometry(features,target,names)
    initial=eligible.copy(); projected=[]; children=[]; keep_ids=[]
    for ix,c in eligible.iterrows():
        i=int(c.SOURCE_INDEX); pr=project_profile(normal['values'][i],np.column_stack([basis,geo[i]]))
        zz=pr/np.maximum(normal['errors'][i],1e-30)
        child=np.isfinite(zz)&(np.abs(zz)>=min_abs_z)
        passed=int(child.sum())>=min_targets
        initial.loc[ix,'BACKGROUND_PROJECTION_TARGETS']=int(child.sum())
        initial.loc[ix,'ROBUSTNESS_PASS']=passed
        if passed:
            keep_ids.append(ix); projected.append(pr); children.append(child)
            eligible.loc[ix,'OOF_N_LARGE_RESIDUAL']=c.N_LARGE_RESIDUAL
            eligible.loc[ix,'N_LARGE_RESIDUAL']=int(child.sum())
            eligible.loc[ix,'PROGRAM_KIND']='frozen_background_residual'
    initial.to_csv(out/'initial-candidates.tsv',sep='\t',index=False,na_rep='.')
    eligible=eligible.loc[keep_ids].copy()
    eligible.to_csv(out/'candidates.tsv',sep='\t',index=False,na_rep='.')
    all_rows.to_csv(out/'all-locus-residuals.tsv',sep='\t',index=False,na_rep='.')
    indices=eligible.SOURCE_INDEX.to_numpy(int)
    np.savez_compressed(out/'candidate_programs.npz',profiles=np.array(projected).reshape(len(eligible),len(names)),
                        primary_values=normal['values'][indices],primary_se=normal['errors'][indices],
                        background_basis=basis,geometry_profiles=geo,children=np.array(children,dtype=bool).reshape(len(eligible),len(names)),
                        excluded_blocks=np.array(excluded,dtype=str),min_abs_z=float(min_abs_z),
                        proteins=np.array(names,dtype=str),source_index=indices,candidates=eligible.CANDIDATE.to_numpy(dtype=str))
    # Give the low-rank basis a separate, explicitly non-identifiable description.
    mod=np.load(model_dir/'model.npz'); u=mod['u']; v=mod['v']; modules=[]
    if u.shape[1]:
        qu,ru=np.linalg.qr(u,mode='reduced'); qv,rv=np.linalg.qr(v,mode='reduced')
        a,s,b=np.linalg.svd(ru@rv.T,full_matrices=False); left=qu@a; right=qv@b.T
        for h in range(len(s)):
            sw=left[:,h]**2; sw/=max(sw.sum(),1e-30); tw=right[:,h]**2; tw/=max(tw.sum(),1e-30)
            peak=int(np.argmax(sw)); r=source.iloc[peak]
            local=(source.CHR.astype(str)==str(r.CHR)).to_numpy()&(np.abs(source.POS.to_numpy()-r.POS)<=1000000)
            prot_weights=target.assign(WEIGHT=tw).groupby('PROT').WEIGHT.sum().to_numpy()
            modules.append(dict(MODULE=h+1,PEAK_LOCUS=r.LOCUS,CHR=str(r.CHR),POS=int(r.POS),
                EFFECTIVE_LOCI=float(1/np.sum(sw**2)),EFFECTIVE_PROTEINS=float(1/np.sum(prot_weights**2)),
                LOCAL_WEIGHT=float(sw[local].sum()),EXPLAINED_FRACTION=float(s[h]**2/np.sum(s**2)),
                INTERPRETATION='statistical latent basis; rotation/sign not an identified gene or protein'))
        np.savez_compressed(out/'latent_basis.npz',source_loading=left,target_loading=right,singular=s)
    pd.DataFrame(modules,columns=['MODULE','PEAK_LOCUS','CHR','POS','EFFECTIVE_LOCI','EFFECTIVE_PROTEINS','LOCAL_WEIGHT','EXPLAINED_FRACTION','INTERPRETATION']).to_csv(out/'latent_modules.tsv',sep='\t',index=False)
    # Cis instruments can be candidate maternal anchors, but overlapping genes remain ambiguous.
    anchors=[]
    pv=np.load(Path(matrix_dir)/'p.npy',mmap_mode='r') if (Path(matrix_dir)/'p.npy').is_file() else None
    for i,r in source.iterrows():
        js=np.flatnonzero(observed[i]&cis[i])
        js=[j for j in js if (float(pv[i,j])<=5e-8 if pv is not None and np.isfinite(pv[i,j]) else abs(float(beta[i,j]/se[i,j]))>=5.451)]
        names=set(target.iloc[js].PROT) if js else set()
        for j in js:
            anchors.append(dict(LOCUS=r.LOCUS,PROT=target.iloc[j].PROT,ASSAY=target.iloc[j].ASSAY,
                BETA=float(beta[i,j]),SE=float(se[i,j]),EA=r.EA,NEA=r.NEA,N_POSSIBLE_MOTHERS=len(names),
                ORIENTATION_TO_INCREASE_PROTEIN=1 if beta[i,j]>0 else -1,STATUS='cis anchor, not proven mediation'))
    pd.DataFrame(anchors,columns=['LOCUS','PROT','ASSAY','BETA','SE','EA','NEA','N_POSSIBLE_MOTHERS','ORIENTATION_TO_INCREASE_PROTEIN','STATUS']).to_csv(out/'maternal-anchors.tsv',sep='\t',index=False)
    result=dict(schema='phole-discovery-v2',matrix_dir=str(Path(matrix_dir).resolve()),model_dir=str(model_dir.resolve()),
        source_hash=meta.get('source_hash'),n_candidates=len(eligible),n_latent_modules=len(modules),
        min_targets=min_targets,min_abs_z=min_abs_z,max_candidates=max_candidates,
        selection='primary-only out-of-fold nomination, followed by leave-candidate-block-out background projection; one nomination per block',
        initial_candidates=len(initial),background_method='frozen primary subspace trained excluding all initially nominated source blocks; candidate-specific coordinate features projected out',
        source_signature=checksum(Path(matrix_dir)/'source.tsv'),target_signature=checksum(Path(matrix_dir)/'target.tsv'),
        candidate_status='hypotheses requiring independent profile replication, not discovered causal black holes',
        latent_identifiability='factor rotations and signs do not identify a unique unmeasured gene',
        program_signature=checksum(out/'candidate_programs.npz'),candidate_table_signature=checksum(out/'candidates.tsv'),
        synthetic=bool(meta.get('synthetic',False)))
    atom_json(out/'discovery.manifest.json',result); return result

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--matrix-dir',required=True); p.add_argument('--model-dir',required=True); p.add_argument('--out',required=True)
    p.add_argument('--min-targets',type=int,default=10); p.add_argument('--min-abs-z',type=float,default=3.); p.add_argument('--max-candidates',type=int,default=50)
    a=p.parse_args()
    if a.min_targets<3 or a.min_abs_z<=0 or a.max_candidates<1: p.error('Invalid discovery thresholds')
    print(json.dumps(discover(a.matrix_dir,a.model_dir,a.out,a.min_targets,a.min_abs_z,a.max_candidates),indent=2))
if __name__=='__main__': main()
