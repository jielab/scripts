#!/usr/bin/env python3
"""Explicit host data adapter. Reuses supplied survival outcomes; never invents dates."""
from pathlib import Path
from dataclasses import asdict,replace
import argparse,json
import numpy as np
import pandas as pd
from research_io import id_table,sha256,write_tables
from run_abm_twin import load_config
from abm_twin_core import ROLES,tie_hash,json_write


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--cohort-file',required=True);a.add_argument('--omics-file',required=True)
    a.add_argument('--config',required=True);a.add_argument('--output',required=True)
    a.add_argument('--time-col',default='time');a.add_argument('--event-col',default='event')
    a.add_argument('--prevalent-col',default='known_cad_at_index');a.add_argument('--clinical',default='')
    a.add_argument('--group-col');a.add_argument('--assume-unrelated',action='store_true');a.add_argument('--role-col')
    x=a.parse_args();cfg=load_config(x.config);d=id_table(x.cohort_file);m=id_table(x.omics_file)
    out=Path(x.output)
    if out.exists() and any(out.iterdir()):raise FileExistsError('Prepared output already exists')
    if not {x.time_col,x.event_col,x.prevalent_col}<=set(d):
        raise ValueError('Provide already-derived survival time/event/prevalence from the host outcome code; absent fields are not inferred')
    ids=d.index[d.index.isin(m.index)];d=d.loc[ids];m=m.loc[ids].drop(columns=['eid'])
    time=pd.to_numeric(d[x.time_col],errors='coerce').to_numpy();event=pd.to_numeric(d[x.event_col],errors='coerce').to_numpy()
    prev=pd.to_numeric(d[x.prevalent_col],errors='coerce').to_numpy()
    eligible=(prev==0)&np.isfinite(time)&(time>0)&np.isin(event,[0,1])
    cls=[v.strip() for v in x.clinical.split(',') if v.strip()]
    if set(cls)&{x.time_col,x.event_col,x.prevalent_col}:raise ValueError('Outcome leakage into clinical inputs')
    C=d[cls].apply(pd.to_numeric,errors='raise').to_numpy(float) if cls else np.empty((len(d),0))
    X=np.c_[m.apply(pd.to_numeric,errors='raise').to_numpy(float),C]
    names=list(m.columns)+['clinical::'+v for v in cls]
    if len(set(names))!=len(names):raise ValueError('Duplicate predictor names')
    if x.group_col:
        if x.group_col not in d or d[x.group_col].isna().any():raise ValueError('Missing family-group IDs')
        groups=d[x.group_col].astype(str).to_numpy()
    elif x.assume_unrelated:groups=np.asarray(ids,str)
    else:raise ValueError('Specify family IDs or explicitly --assume-unrelated')
    if x.role_col:
        role=d[x.role_col].astype(str).to_numpy()
    else:
        v=np.array([tie_hash(g,cfg.seed)/2**64 for g in groups]);role=np.full(len(v),'test',dtype='U25')
        edges=[0,.50,.60,.68,.76,.84]
        for j,r in enumerate(ROLES):role[(v>=edges[j])&(v<edges[j+1])]=r
    if not set(role)<=set(ROLES)|{'test'}:raise ValueError('Invalid roles; use exact explicit build/tune_model/tune_gate/calibration_fit/calibration_audit/test')
    rg=pd.DataFrame({'g':groups,'r':role})
    if rg.groupby('g').r.nunique().max()>1:raise ValueError('Family crosses roles')
    groups=np.asarray(groups,str);role=np.asarray(role,str)
    X,time,event,groups,role,ids=[v[eligible] for v in (X,time,event,groups,role,np.asarray(ids,str))]
    event=event.astype(int);dev=role!='test'
    if set(role)!=set(ROLES)|{'test'}:raise ValueError('Each development role and test must have eligible people')
    cfg=replace(cfg,always_include_indices=tuple(range(len(m.columns),X.shape[1])))
    out.mkdir(parents=True);local=out/'_LOCAL_ONLY';local.mkdir(mode=0o700)
    np.savez_compressed(local/'development.npz',X=X[dev],time=time[dev],event=event[dev],eid=ids[dev],group=groups[dev],role=role[dev],feature_names=np.asarray(names))
    np.savez_compressed(local/'test_X.npz',X=X[~dev],eid=ids[~dev],group=groups[~dev],feature_names=np.asarray(names))
    np.savez_compressed(local/'test_truth.npz',time=time[~dev],event=event[~dev],eid=ids[~dev])
    json_write(out/'config.json',asdict(cfg))
    json_write(out/'input_provenance.json',dict(cohort_sha256=sha256(x.cohort_file),omics_sha256=sha256(x.omics_file),
               time_column=x.time_col,event_column=x.event_col,prevalent_column=x.prevalent_col,
               group_assumption=x.group_col or 'declared_unrelated',
               outcome_definition='unchanged supplied host survival fields; death handled according to supplied event/censoring encoding'))
    counts=pd.DataFrame({'role':role}).groupby('role').size().reset_index(name='N')
    write_tables(out/'aggregate',{'cohort_roles':counts},'ABM twin',cfg.layer,cfg.trait,'prepared')
    print(json.dumps(dict(output=str(out),eligible_N=len(X),development_N=int(dev.sum()),test_N=int((~dev).sum()),input_features=X.shape[1],molecular_panel_budget=cfg.panel_size)))
if __name__=='__main__':main()
