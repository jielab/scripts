#!/usr/bin/env python3
"""GRID entry point, twin workflow, evaluation, and aggregate publication.

Model classes stay in twin_core so saved joblib models and their frozen hashes
remain compatible. This module also supplies the browser's numeric IO helpers.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import uuid
import zipfile
import xml.etree.ElementTree as ET
from dataclasses import asdict, replace
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
from threadpoolctl import threadpool_limits

from twin_core import (
    Config, TwinModel, assert_roles, assign_roles, hash_number, id_hash,
    state_hash, validate_ids,
)
from abm_routes import (
    VERSION as GRID_ABM_VERSION,
    apply_abm_routes, evaluate_abm_routes, fit_abm_routes,
    load_abm_routes_config, write_source_provenance,
    freeze_matched_controls,
)

# GRID_ABM_FINAL_20261017
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / 'config/twin.json'
DEFAULT_ABM_ROUTES_CONFIG = ROOT / 'config/abm_routes.json'


def abm_routes_config():
    return load_abm_routes_config(os.getenv('GRID_ABM_ROUTES_CONFIG', str(DEFAULT_ABM_ROUTES_CONFIG)))


# Numeric input and aggregate publication

IDENTIFIER_COLUMNS={'eid','iid','fid','participantid','individualid','subjectid','patientid','queryid','targetid','donorid','familyid','groupid','group','id','id1','id2','sampleid','ukbid','ukbeid'}


def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
    return h.hexdigest()


def json_write(path,value):
    def safe(v):
        if isinstance(v,pd.DataFrame):return {'type':'dataframe_metadata','rows':len(v),'columns':list(map(str,v.columns)),'contents_omitted':True}
        if isinstance(v,pd.Series):return {'type':'series_metadata','rows':len(v),'name':str(v.name),'contents_omitted':True}
        if isinstance(v,pd.Timestamp):return v.isoformat()
        if v is pd.NA or v is pd.NaT:return None
        if isinstance(v,dict):return {str(k):safe(x) for k,x in v.items()}
        if isinstance(v,(list,tuple)):return [safe(x) for x in v]
        if isinstance(v,np.ndarray):return safe(v.tolist())
        if isinstance(v,np.generic):return safe(v.item())
        if isinstance(v,Path):return str(v)
        if isinstance(v,float) and not np.isfinite(v):return None
        return v
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+'.tmp.'+uuid.uuid4().hex)
    tmp.write_text(json.dumps(safe(value),ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');os.replace(tmp,p)


def read_table(path):
    p=Path(path)
    if p.suffix=='.rds':
        try:import rdata
        except ImportError as e:raise RuntimeError('RDS input needs rdata>=1.1; use the existing grid environment or supply a TSV') from e
        d=rdata.read_rds(p)
        if isinstance(d,dict):
            candidates=[x for x in d.values() if isinstance(x,pd.DataFrame)]
            if len(candidates)!=1:raise ValueError('Expected exactly one data.frame in RDS')
            d=candidates[0]
        if not isinstance(d,pd.DataFrame):raise ValueError('Expected a data.frame in RDS')
        return d
    if p.suffix=='.parquet':return pd.read_parquet(p)
    sep=',' if '.csv' in p.name else '\t'
    # Explicit string IDs preserve values such as 00012; pandas otherwise infers integers.
    head=pd.read_csv(p,sep=sep,header=None,nrows=1).iloc[0].astype(str)
    if head.duplicated().any():raise ValueError('Duplicate source columns')
    dtype={c:str for c in head if re.sub(r'[^a-z0-9]','',c.lower()) in IDENTIFIER_COLUMNS or c=='ancestry'}
    return pd.read_csv(p,sep=sep,dtype=dtype,low_memory=False)


def identified(d):
    d=pd.DataFrame(d).copy()
    if 'eid' not in d:
        c=next((c for c in ('IID','#IID','ID') if c in d),None)
        if c is None:raise ValueError('Participant table needs eid or IID')
        d=d.rename(columns={c:'eid'})
    if d.eid.isna().any():raise ValueError('Missing IDs')
    if pd.api.types.is_numeric_dtype(d.eid):
        v=d.eid.to_numpy(float)
        if not np.isfinite(v).all() or np.any(v!=np.floor(v)):raise ValueError('Noninteger numeric IDs')
        d['eid']=d.eid.astype('int64').astype(str)
    else:d['eid']=d.eid.astype(str).str.strip()
    if d.eid.duplicated().any() or d.eid.eq('').any():raise ValueError('Duplicate/blank IDs')
    return d


def write_aggregate(directory,tables,metadata):
    out=Path(directory);out.mkdir(parents=True,exist_ok=True);manifest=[]
    for name,value in tables.items():
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*',name):raise ValueError('Unsafe table name')
        d=pd.DataFrame(value).copy()
        if not len(d.columns):d=pd.DataFrame(columns=['status'])
        for c in d.select_dtypes(include='object'):
            d[c]=d[c].map(lambda v:json.dumps(v,ensure_ascii=False,default=str) if isinstance(v,(dict,list,tuple)) else v)
        path=out/(name+'.tsv.gz');tmp=out/(name+'.tmp.gz')
        d.to_csv(tmp,sep='\t',index=False,float_format='%.17g',compression='gzip');os.replace(tmp,path)
        manifest.append(dict(name=name,path=path.name,sha256=sha256(path),rows=len(d),columns=list(map(str,d.columns)),**metadata))
    json_write(out/'table_manifest.json',manifest)
    return manifest


def xlsx_values(path,max_bytes=256*1024**2):
    """Read cached Excel *values* only. Never evaluate formula text or execute macros.
    Supports the current GRID ordinary XLSX aggregates. Empty/status sheets retained.
    This parser is for importing legacy evidence, not editing spreadsheets.
    """
    ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    def text(e):return ''.join(t.text or '' for t in e.findall('.//s:t',ns))
    with zipfile.ZipFile(path) as z:
        if sum(i.file_size for i in z.infolist())>max_bytes:raise ValueError('Expanded workbook exceeds configured limit')
        strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():strings=[text(si) for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)]
        relationships={e.attrib['Id']:e.attrib['Target'] for e in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        w=ET.fromstring(z.read('xl/workbook.xml'));out={}
        for sheet in w.findall('s:sheets/s:sheet',ns):
            rid=sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'];target=relationships[rid]
            target=target.lstrip('/') if target.startswith('/') else 'xl/'+target
            if '..' in Path(target).parts:raise ValueError('Unsafe XLSX relationship')
            rows=[]
            for row in ET.fromstring(z.read(target)).findall('s:sheetData/s:row',ns):
                vals=[]
                for cell in row.findall('s:c',ns):
                    letters=re.sub('[0-9]','',cell.get('r','A1'));col=0
                    for char in letters:col=26*col+ord(char)-64
                    while len(vals)<col:vals.append(None)
                    v=cell.find('s:v',ns);v=v.text if v is not None else None;typ=cell.get('t')
                    if typ=='s' and v is not None:v=strings[int(v)]
                    elif typ=='inlineStr':v=text(cell)
                    elif typ=='b' and v is not None:v=v=='1'
                    elif v is not None and typ not in ('str','e'):
                        try:v=float(v)
                        except ValueError:pass
                    vals[col-1]=v
                rows.append(vals)
            if not rows:out[sheet.attrib['name']]=pd.DataFrame();continue
            width=max(map(len,rows));header=rows[0]+[None]*(width-len(rows[0]));names=[str(v) if v is not None else f'column_{j+1}' for j,v in enumerate(header)]
            if len(set(names))!=len(names):raise ValueError('Duplicate worksheet columns')
            out[sheet.attrib['name']]=pd.DataFrame([r+[None]*(width-len(r)) for r in rows[1:]],columns=names)
        return out

# Prepared cohorts and benchmark adapters

class PreparedSourceChanged(ValueError):
    """A prepared feature cache must be rebuilt from current source files."""


def load_source(prepared_path,trait,data_file=None,schema_file=None,*,allow_unfinished_features=False):
    if data_file:
        if not schema_file:raise ValueError('Explicit external features require --schema JSON')
        spec=json.loads(Path(schema_file).read_text());d=identified(read_table(data_file))
        groups=spec['feature_groups'];metadata=spec.get('metadata',{})
        metadata.update(input_source=str(data_file),input_sha256=sha256(data_file),matrix_provenance='provider_declared')
    else:
        path=Path(prepared_path)
        if not path.is_file():raise FileNotFoundError(str(path)+'; reuse an existing prepared.joblib, or supply --data-file and --schema')
        # Only a trusted local file produced by the user's GRID pipeline. Never load
        # arbitrary downloaded pickle/joblib files through the public viewer.
        prepared=joblib.load(path)
        if not {'tables','feature_groups','cohort'}<=set(prepared):raise ValueError('Unrecognized GRID prepared cache')
        allowed_states=(None,'features_ready','prepared') if allow_unfinished_features else (None,'features_ready')
        if prepared.get('state') not in allowed_states:raise ValueError('Prepared cache has not completed genetic feature preparation')
        if trait not in prepared['tables']:raise ValueError('Trait absent from prepared GRID cache')
        d=identified(prepared['tables'][trait]);groups=prepared['feature_groups'][trait];cohort=identified(prepared['cohort'])
        if 'split' not in d or 'split' not in cohort:raise ValueError('Missing frozen outer split')
        c=cohort.set_index('eid');v=d.set_index('eid')
        if set(v.index)!=set(c.index):raise ValueError('Trait input must cover the entire frozen original cohort')
        canon=lambda s:s.astype(str).replace({'training':'train','testing':'test'})
        if not np.array_equal(canon(v['split']),canon(c.loc[v.index,'split'])):raise ValueError('Trait split differs from stored original cohort')
        if 'family_id' in c and ('family_id' not in v or not np.array_equal(v.family_id.astype(str),c.loc[v.index,'family_id'].astype(str))):
            raise ValueError('Trait reference groups differ from stored original cohort')
        metadata=dict(prepared.get('metadata',{}));metadata.update(input_source=str(path),input_sha256=sha256(path),matrix_provenance='existing_GRID_prepared_cache',outer_cohort_N=len(c),prepared_settings=prepared.get('settings',{}))
        for src in metadata.get('sources',[]):
            if not isinstance(src,dict) or 'path' not in src:continue
            source=Path(src['path'])
            if not source.is_file():raise FileNotFoundError('Prepared source missing: '+str(source))
            stat=source.stat()
            if ('size' in src and stat.st_size!=src['size']) or ('mtime_ns' in src and stat.st_mtime_ns!=src['mtime_ns']):
                raise PreparedSourceChanged('Prepared source changed; rebuild features on the frozen outer roster: '+str(source))
    return d,groups,metadata


def numeric_cohort(d,groups,config,roles_file=None):
    d=identified(d)
    # Genetic roles are explicit, not inferred from names or from test associations.
    allowed={'covariates','csx','ancestry','frequency','evolution','permuted_evolution','disco','genetic_profile'}
    if set(groups)-allowed:raise ValueError('Unknown feature groups: '+str(set(groups)-allowed))
    def union(keys):return list(dict.fromkeys(str(c) for k in keys for c in groups.get(k,[])))
    # Disco remains an explicit comparator; do not leak its fitted output into the
    # primary profile merely because it was available in the source matrix.
    genetic=union(['csx','ancestry','frequency','evolution','genetic_profile']);clinical=union(['covariates']);declared_clinical=clinical.copy()
    # Declared adjustment covariates never define the primary twin geometry,
    # including PCs declared in both source groups.
    genetic=[c for c in genetic if c not in declared_clinical]
    if set(genetic)&set(declared_clinical):raise AssertionError('Clinical covariate in primary geometry')
    comparator_extra=[c for c in union(['disco']) if c not in genetic and c not in clinical]
    if not genetic:raise ValueError('No genetic profile features supplied')
    features=genetic+clinical
    forbidden={'y','outcome','event','time','target','eid','family_id','split','role','ancestry'}
    if set(features)&forbidden:raise ValueError('Identity/outcome columns may not be predictors')
    if not set(features+comparator_extra)<=set(d):raise ValueError('Declared predictor missing from prepared table')
    # Extra Disco score(s) are outside primary model/QC; add as separate calibrated
    # comparators later rather than distort genetic missingness or twin geometry.
    X=d[features].apply(pd.to_numeric,errors='raise').to_numpy(dtype=float,na_value=np.nan)
    if 'split' not in d or 'y' not in d:raise ValueError('Source requires frozen split and observed phenotype column y')
    split=d.split.astype(str).replace({'training':'train','testing':'test'})
    if not set(split)<= {'train','test'} or set(split)!= {'train','test'}:raise ValueError('Keep the existing train/test outer roster')
    ids,group=validate_ids(d.eid,d.family_id if 'family_id' in d else None)
    if set(group[split=='train'])&set(group[split=='test']):raise ValueError('Supplied reference groups cross outer split')
    # Test outcomes are split off here and will not be passed to the fit API.
    dev=(split=='train').to_numpy();test=~dev;y=pd.to_numeric(d.y,errors='raise').to_numpy(dtype=float,na_value=np.nan)
    observed=np.isfinite(y);fit=dev&observed
    config=replace(config,grouping_status=config.grouping_status if 'family_id' in d else 'individual_unverified')
    if roles_file:
        r=identified(read_table(roles_file)).set_index('eid')
        if 'role' not in r or not set(ids[fit])<=set(r.index):raise ValueError('Roles file does not cover the development roster')
        roles=r.loc[ids[fit],'role'].astype(str).to_numpy()
        aliases={'calibration_audit':'audit'};roles=np.asarray([aliases.get(v,v) for v in roles])
    else:roles=assign_roles(group[fit],config)
    assert_roles(group[fit],roles)
    raw=dict(X=X[fit],y=y[fit],ids=ids[fit],groups=group[fit],roles=roles,features=features,
             genetic_indices=list(range(len(genetic))),covariate_indices=[features.index(c) for c in declared_clinical],
             csx_indices=[features.index(c) for c in groups.get('csx',[]) if c in features],disco_indices=[])
    query=dict(X=X[test],ids=ids[test],groups=group[test],features=features)
    ancestry=d['ancestry'].fillna('UNASSIGNED').astype(str).to_numpy() if 'ancestry' in d else np.full(len(d),'UNASSIGNED')
    truth=pd.DataFrame(dict(eid=ids[test],y=y[test],ancestry=ancestry[test]))
    flow=dict(original_N=len(d),outer_train_N=int(dev.sum()),outer_test_N=int(test.sum()),development_missing_y_N=int(np.sum(dev&~observed)),test_missing_y_N=int(np.sum(test&~observed)),test_roster_hash=id_hash(ids[test]),grouping_status=config.grouping_status,
              raw_SNP_QC='only for explicit variant_dosage input; summaries cannot reconstruct unseen SNP combinations',genetic_columns=genetic,clinical_columns=declared_clinical)
    audit_truth=pd.DataFrame(dict(eid=ids[fit][roles=='audit'],y=y[fit][roles=='audit'],ancestry=ancestry[fit][roles=='audit']))
    raw['y']=raw['y'].copy();raw['y'][roles=='audit']=np.nan
    extra=dict(d=d,fit=fit,test=test,role=roles,disco_columns=comparator_extra,audit_truth=audit_truth)
    return raw,query,truth,flow,config,extra


def calibrate_disco(extra,model):
    """Same-budget local calibration of supplied DiscoDivas scores, not a rerun of Disco."""
    from twin_core import Predictor,Calibration,loss
    from sklearn.impute import SimpleImputer
    from sklearn.preprocessing import StandardScaler
    cfg=model.config;columns=extra['disco_columns']
    if not columns:return None,None,dict(method='Disco_reference',status='source_scores_unavailable')
    d=extra['d'];fit=extra['fit'];test=extra['test'];r=extra['role']
    # Preserve covariates and original score values. Explicitly separate from raw
    # published DiscoDivas-untuned versus locally calibrated outcome predictions.
    cov=[model.features[i] for i in model.prep.clinical]
    names=cov+columns;x=d[names].apply(pd.to_numeric,errors='raise').to_numpy(dtype=float,na_value=np.nan)
    yy=pd.to_numeric(d.y,errors='raise').to_numpy(dtype=float,na_value=np.nan)[fit];xd=x[fit];xt=x[test]
    b=(r=='build');t=(r=='tune_model');c=(r=='calibration_fit')
    qc=model.prep.qc(d.loc[fit,model.features].to_numpy(dtype=float,na_value=np.nan))
    usable=qc&np.isfinite(xd[:,len(cov):]).all(1)
    if not (usable[b].all() and usable[t].all() and usable[c].all()):
        # A different training/calibration eligibility changes comparability;
        # report unavailability rather than quietly dropping rows for Disco.
        # QC failures in the primary do not count as missing benchmark scores.
        missing=qc&~np.isfinite(xd[:,len(cov):]).all(1)
        if np.any(missing&(b|t|c)):return None,None,dict(method='Disco_reference',status='incomplete_development_scores')
    if not np.isfinite(xt[:,len(cov):]).all():return None,None,dict(method='Disco_reference',status='incomplete_test_scores')
    imputer=SimpleImputer(strategy='median',keep_empty_features=True).fit(xd[b&qc]);scale=StandardScaler().fit(imputer.transform(xd[b&qc]));z=scale.transform(imputer.transform(xd));best=None
    for parameter in (cfg.logistic_c_grid if cfg.trait_type=='binary' else cfg.ridge_grid):
        m=Predictor().fit(z[b],yy[b],qc[b].astype(float),cfg,parameter);p=m.predict(z[t]);sc=float(np.mean(loss(yy[t][qc[t]],p[qc[t]],cfg.trait_type=='binary')))
        if best is None or sc<best[0]:best=(sc,m)
    m=best[1];cal=Calibration().fit(m.predict(z[c][qc[c]]),yy[c][qc[c]],cfg.trait_type=='binary')
    p=cal.predict(m.predict(scale.transform(imputer.transform(xt))))
    bundle=dict(columns=names,imputer=imputer,scale=scale,model=m,calibration=cal)
    return p,bundle,dict(method='Disco_reference',status='local_train_only_calibration',scope='same build/tune/cal roles; not a new run of upstream DiscoDivas')


def frozen_comparison_roster(root):
    """Read the immutable preparation roster; never require trained outputs."""
    root=Path(root)
    cache=root/'prsformer/cache'
    sources=[]
    data=cache/'data.tsv.gz'
    if data.is_file():
        metadata=json.loads((cache/'prepare.json').read_text())
        if hashlib.sha256(data.read_bytes()).hexdigest()!=metadata.get('sha256',{}).get('data.tsv.gz'):
            raise ValueError('Frozen roster data checksum mismatch')
        sources.append((str(data),identified(read_table(data))))
    published=root/'prsformer/split.rds'
    if published.is_file():sources.append((str(published),identified(read_table(published))))
    split=cache/'roster/prsformer.split.tsv.gz'
    families=cache/'roster/families.tsv.gz'
    if split.is_file() or families.is_file():
        if not split.is_file() or not families.is_file():raise ValueError('Incomplete frozen roster files')
        a=identified(read_table(split));b=identified(read_table(families))
        if set(a.eid)!=set(b.eid):raise ValueError('Frozen split/family IDs differ')
        sources.append((str(split),a.merge(b,on='eid',validate='one_to_one')))
    if not sources:raise FileNotFoundError('PRSformer frozen preparation roster missing under '+str(root)+'. Run first: ./grid.sh prsformer --traits height,ldl,t2dm')
    roster=None
    for source,current in sources:
        if not {'eid','split','group'}<=set(current):raise ValueError('Frozen roster lacks eid, split or group: '+source)
        current=current[['eid','split','group']].copy()
        if current.isna().any().any():raise ValueError('Missing frozen roster values: '+source)
        # TSV object strings and RDS StringDtype/factors must compare by value.
        current['split']=current.split.astype(str)
        current['group']=current.group.astype(str)
        if not set(current.split)<={'train','validation','test'} or not {'train','test'}<=set(current.split):
            raise ValueError('Invalid frozen roster partitions: '+source)
        if current.groupby('group').split.nunique().max()!=1:raise ValueError('Frozen family crosses outer roles')
        current=current.sort_values('eid').reset_index(drop=True)
        if roster is not None and not current.equals(roster):raise ValueError('Frozen comparison rosters disagree: '+source)
        roster=current
    return roster


def verify_existing_prsformer_roster(root,full_cohort,query_ids):
    r=frozen_comparison_roster(root);c=identified(full_cohort)
    if not {'eid','split'}<=set(r) or set(r.eid)!=set(c.eid):raise ValueError('Existing PRSformer outer cohort differs from current GRID cohort')
    rt=r.set_index('eid');test=set(rt.index[rt.split.astype(str).eq('test')])
    if test!=set(map(str,query_ids)):raise ValueError('Existing PRSformer test roster differs from ABM frozen test roster; do not compare or silently refit')
    if 'group' in rt and 'family_id' in c:
        if not np.array_equal(rt.loc[c.eid,'group'].astype(str),c.family_id.astype(str)):
            raise ValueError('Existing PRSformer reference groups differ from ABM groups')
    return dict(status='frozen_comparison_roster_matches_ABM',reused=True,test_N=len(test),test_roster_hash=id_hash(query_ids))


def verified_prsformer(root,trait,full_cohort,truth,config,covariates,endpoint_definition):
    """Use existing PRSformer outcome-scale predictions only with strict lineage.
    Called during evaluation, never to choose ABM roles/thresholds/parameters.
    """
    root=Path(root);roster=root/'prsformer/split.rds';path=root/'prsformer/scores'/trait/'predictions.rds'
    if not roster.is_file() or not path.is_file():return None,dict(method='PRSformer',status='saved_prediction_or_roster_unavailable')
    r=identified(read_table(roster));c=identified(full_cohort);v=identified(read_table(path))
    if not {'eid','split'}<=set(r) or set(r.eid)!=set(c.eid):raise ValueError('PRSformer outer cohort differs')
    a=r.set_index('eid').loc[c.eid,'split'].eq('test').to_numpy();b=c.split.astype(str).replace({'testing':'test'}).eq('test').to_numpy()
    if not np.array_equal(a,b):raise ValueError('PRSformer test split differs')
    if not {'prediction','outcome','split','covariate_names','endpoint_definition'}<=set(v):raise ValueError('PRSformer lacks required prediction/endpoint/covariate provenance')
    if not v.split.eq('test').all() or set(v.eid)!=set(c.loc[b,'eid']):raise ValueError('PRSformer scores do not cover the full original test roster')
    if v.endpoint_definition.isna().any() or set(v.endpoint_definition.astype(str))!={endpoint_definition}:raise ValueError('PRSformer endpoint definition differs')
    if v.covariate_names.isna().any() or set(v.covariate_names.astype(str))!={','.join(covariates)}:raise ValueError('PRSformer covariate names differ')
    z=v.set_index('eid').loc[truth.eid];cd=c.set_index('eid').loc[truth.eid]
    if not np.allclose(pd.to_numeric(z.outcome).to_numpy(),truth.y.to_numpy(),equal_nan=True,atol=1e-7,rtol=1e-7):raise ValueError('PRSformer outcome/units differ')
    for name in covariates:
        col='covariate.'+name
        if col not in z or name not in cd or not np.allclose(pd.to_numeric(z[col]),pd.to_numeric(cd[name]),equal_nan=True,atol=1e-8,rtol=1e-7):
            raise ValueError('PRSformer covariate values differ: '+name)
    p=pd.to_numeric(z.prediction,errors='raise').to_numpy()
    return p,dict(method='PRSformer',status='same_roster_endpoint_and_covariates_verified',source_sha256=sha256(path),scope='upstream model may have a larger training-label budget; audit comparison is not claimed')

# Frozen-model evaluation

METRIC_DIRECTIONS={"MSE":"lower","RMSE":"lower","MAE":"lower","prediction_R2":"higher","MSE_over_train_variance":"lower","AUC":"higher","AUPRC":"higher","Brier":"lower","LogLoss":"lower"}


def metrics(y,p,binary,training_variance,weights=None):
    y=np.asarray(y,float);p=np.asarray(p,float);w=np.ones(len(y)) if weights is None else np.asarray(weights,float)
    good=w>0;y,p,w=y[good],p[good],w[good];n=float(w.sum())
    if not n:return {k:np.nan for k in (('AUC','AUPRC','Brier','LogLoss') if binary else ('MSE','RMSE','MAE','prediction_R2','MSE_over_train_variance'))}
    if not np.isfinite(y).all() or not np.isfinite(p).all():raise ValueError('Nonfinite prediction inside a locked common evaluation mask')
    error=y-p;mse=float(np.average(error**2,weights=w))
    if binary:
        if np.any((p<0)|(p>1)):raise ValueError('Binary benchmarks require probabilities, not PRS values/logits')
        pp=np.clip(p,1e-7,1-1e-7);two=len(np.unique(y))==2
        return dict(AUC=float(roc_auc_score(y,p,sample_weight=w)) if two else np.nan,
                    AUPRC=float(average_precision_score(y,p,sample_weight=w)) if two else np.nan,
                    Brier=mse,LogLoss=float(np.average(-(y*np.log(pp)+(1-y)*np.log1p(-pp)),weights=w)))
    var=float(np.average((y-np.average(y,weights=w))**2,weights=w))
    return dict(MSE=mse,RMSE=math.sqrt(mse),MAE=float(np.average(np.abs(error),weights=w)),
                prediction_R2=1-mse/var if var>0 else np.nan,
                MSE_over_train_variance=mse/training_variance if training_variance>0 else np.nan)


class FixedMetricCache:
    """Pre-sort frozen prediction ties once; bootstrap changes weights, not scores.
    Weighted AUC/AP then use O(n) cumulative sums instead of O(n log n) sorting
    in every bootstrap draw. Losses and outcome moments are also precomputed.
    """
    def __init__(self,y,p,binary,training_variance):
        self.y=np.asarray(y,float);self.p=np.asarray(p,float);self.binary=binary;self.tv=training_variance
        if not np.isfinite(self.y).all() or not np.isfinite(self.p).all():raise ValueError('Invalid common-mask metric input')
        self.err2=(self.y-self.p)**2;self.abserr=np.abs(self.y-self.p)
        if binary:
            if not np.isin(self.y,[0,1]).all() or np.any((self.p<0)|(self.p>1)):raise ValueError('Invalid binary probability metric')
            values,self.bin=np.unique(self.p,return_inverse=True);self.nb=len(values)
            pp=np.clip(self.p,1e-7,1-1e-7);self.logloss=-(self.y*np.log(pp)+(1-self.y)*np.log1p(-pp))
    def __call__(self,w):
        w=np.asarray(w,float);total=w.sum()
        keys=('AUC','AUPRC','Brier','LogLoss') if self.binary else ('MSE','RMSE','MAE','prediction_R2','MSE_over_train_variance')
        if not total:return {k:np.nan for k in keys}
        mse=float(w@self.err2/total)
        if self.binary:
            pos=np.bincount(self.bin,weights=w*self.y,minlength=self.nb)
            neg=np.bincount(self.bin,weights=w*(1-self.y),minlength=self.nb)
            np_,nn=pos.sum(),neg.sum();auc=ap=np.nan
            if np_>0 and nn>0:
                auc=float(pos@(np.cumsum(neg)-.5*neg)/(np_*nn))
                tp=np.cumsum(pos[::-1]);fp=np.cumsum(neg[::-1]);den=tp+fp
                precision=np.divide(tp,den,out=np.zeros_like(tp),where=den>0)
                ap=float(pos[::-1]@precision/np_)
            return dict(AUC=auc,AUPRC=ap,Brier=mse,LogLoss=float(w@self.logloss/total))
        mean=float(w@self.y/total);var=float(w@((self.y-mean)**2)/total)
        return dict(MSE=mse,RMSE=math.sqrt(mse),MAE=float(w@self.abserr/total),prediction_R2=1-mse/var if var>0 else np.nan,MSE_over_train_variance=mse/self.tv if self.tv>0 else np.nan)


def composition_control(selected,eligible,risk,ancestry,ids,edges,seed):
    """Outcome-blind descriptive control; batch-composition dependent, NOT deployment."""
    bins=np.searchsorted(edges,risk);mask=np.zeros(len(ids),bool)
    for pop in np.unique(ancestry):
        for b in np.unique(bins[ancestry==pop]):
            scope=(ancestry==pop)&(bins==b);ix=np.flatnonzero(eligible&scope)
            order=sorted(ix,key=lambda i:hash_number(ids[i],seed));count=int(np.sum(selected&scope))
            mask[np.array(order[:count],int)]=True
    return mask


def evaluate(model,prediction,truth,*,external=None,bootstrap=None):
    """truth: eid,y and optionally ancestry. All test subjects remain in the flow."""
    cfg=model.config
    if prediction['frozen_hash']!=state_hash(model):raise ValueError('Frozen model identity changed')
    truth=pd.DataFrame(truth).copy();truth['eid']=truth.eid.astype(str)
    ids=prediction['eid'];groups=prediction['group'];n=len(ids)
    if truth.eid.duplicated().any() or set(truth.eid)!=set(ids):raise ValueError('Truth must cover exactly the frozen query roster')
    truth=truth.set_index('eid').loc[ids]
    y=pd.to_numeric(truth.y,errors='raise').to_numpy(float);known=np.isfinite(y)
    if cfg.trait_type=='binary' and not np.isin(y[known],[0,1]).all():raise ValueError('Binary endpoint must be coded 0/1')
    pop=truth.ancestry.fillna('UNASSIGNED').astype(str).to_numpy() if 'ancestry' in truth else np.full(n,'UNASSIGNED')
    qc=prediction['diagnostics']['technical_qc'];base_scope=known&qc
    if base_scope.sum()<2:raise ValueError('Insufficient evaluable technical-QC passing outcomes')
    p={k:np.asarray(v,float) for k,v in prediction['predictions'].items()};availability=[]
    for name,values in (external or {}).items():
        if name in p:raise ValueError('External method collides with an internal method name')
        values=np.asarray(values,float)
        if len(values)!=n or not np.isfinite(values[base_scope]).all():
            availability.append(dict(method=name,status='unavailable_incomplete_frozen_roster'));continue
        if cfg.trait_type=='binary' and np.any((values[base_scope]<0)|(values[base_scope]>1)):
            raise ValueError(f'{name}: expected probabilities, not uncalibrated scores')
        p[name]=values;availability.append(dict(method=name,status='same_roster_prediction_available'))
    for name,v in p.items():
        if not np.isfinite(v[base_scope]).all():raise ValueError(f'{name}: incomplete common-mask predictions; do not silently use pairwise intersections')
    selected=prediction['selected'];masks={'all':base_scope,'selected':base_scope&selected,'rejected':base_scope&~selected}
    for name,v in prediction['masks'].items():masks[name]=base_scope&v
    # Count matching uses predictor values only, before looking at y. Missing-y
    # proportions are then reported rather than balancing using test outcomes.
    mask_random=composition_control(selected,qc,np.nan_to_num(p[model.fallback]),pop,ids,model.risk_edges,cfg.seed+471)
    masks['risk_ancestry_matched_random']=base_scope&mask_random
    pool=np.flatnonzero(qc);risk=np.nan_to_num(p[model.fallback]);center=.5 if cfg.trait_type=='binary' else model.prediction_center
    order=sorted(pool,key=lambda i:(-abs(risk[i]-center),hash_number(ids[i],cfg.seed+472)))
    cc=np.zeros(n,bool);cc[np.asarray(order[:int(selected.sum())],int)]=True
    masks['confidence_count_matched']=base_scope&cc
    scopes=['ALL',*sorted(set(pop))];rows=[];lookup={};binary=cfg.trait_type=='binary'
    for ancestry in scopes:
        am=np.ones(n,bool) if ancestry=='ALL' else pop==ancestry
        total=int(am.sum());cases=int(np.sum(am&known&(y==1))) if binary else None
        for subset,mask in masks.items():
            use=mask&am;count=int(use.sum())
            for name,values in p.items():
                met=metrics(y[use],values[use],binary,model.training_variance)
                row=dict(trait=cfg.trait,ancestry=ancestry,subset=subset,method=name,N=count,original_N=total,
                         actual_coverage=count/total if total else np.nan,selected_policy_coverage=float(np.sum(am&selected)/total) if total else np.nan,
                         cases=int(np.sum(y[use]==1)) if binary else np.nan,controls=int(np.sum(y[use]==0)) if binary else np.nan,
                         case_coverage=float(np.sum(y[use]==1)/cases) if binary and cases else np.nan,
                         outcome_variance=float(np.var(y[use])) if count else np.nan,
                         status='estimated' if count else 'empty',evaluation_scope=prediction['evaluation_scope'],**met)
                rows.append(row);lookup[(ancestry,subset,name)]=met
    primary=cfg.primary;full='ABM_full';comparisons=[
        ('selection_effect',full,'selected',full,'all'),
        ('training_effect_same_subset',primary,'selected',full,'selected'),
        ('combined_effect',primary,'selected',full,'all'),
        ('training_effect_all',primary,'all',full,'all'),
        ('pruning_specificity',primary,'selected','ABM_random','selected'),
        ('beyond_baseline_confidence',full,'selected',full,'confidence_count_matched'),
        ('beyond_risk_ancestry_composition',full,'selected',full,'risk_ancestry_matched_random'),
        ('fallback_policy_all','ABM_policy','all',model.fallback,'all')]
    benchmark_names=[name for name in p if name not in {primary,full,'ABM_soft','ABM_random','Clinical','ABM_policy'}]
    for name in benchmark_names:
        comparisons.extend([(f'benchmark_selection:{name}',name,'selected',name,'all'),
                            (f'benchmark_same_subset:{name}',primary,'selected',name,'selected'),
                            (f'benchmark_combined:{name}',primary,'selected',name,'all')])
    B=int(cfg.bootstrap if bootstrap is None else bootstrap)
    if B<0:raise ValueError('Bootstrap count must be nonnegative')
    metric_names=list(next(iter(lookup.values())))
    contrasts=[]
    # Cluster resampling keeps each selected group nested inside each original draw.
    for ancestry in scopes:
        am=np.ones(n,bool) if ancestry=='ALL' else pop==ancestry
        scope_indices=np.flatnonzero(am);ug,gi=np.unique(groups[am],return_inverse=True)
        if not len(ug):continue
        needed=set((method,subset) for _,a,b,c,d in comparisons for method,subset in [(a,b),(c,d)])
        scoped_masks={subset:masks[subset][am] for _,subset in needed}
        yy=y[am];vv={k:v[am] for k,v in p.items()}
        cached={(method,subset):FixedMetricCache(yy[scoped_masks[subset]],vv[method][scoped_masks[subset]],binary,model.training_variance)
                for method,subset in needed}
        collected={(label,metric):[] for label,*_ in comparisons for metric in metric_names}
        rng=np.random.default_rng(cfg.seed+702)
        for _ in range(B):
            gw=np.bincount(rng.integers(0,len(ug),len(ug)),minlength=len(ug));w=gw[gi]
            bm={}
            for method,subset in needed:
                use=scoped_masks[subset];bm[(method,subset)]=cached[(method,subset)](w[use])
            for label,a,s,b,t in comparisons:
                for metric in metric_names:
                    collected[(label,metric)].append(bm[(a,s)][metric]-bm[(b,t)][metric])
        for label,a,s,b,t in comparisons:
            for metric in metric_names:
                samples=np.asarray(collected[(label,metric)]);samples=samples[np.isfinite(samples)]
                enough=len(samples)>=max(20,math.ceil(.8*B));bounds=np.quantile(samples,[.025,.975]) if enough else [np.nan,np.nan]
                contrasts.append(dict(trait=cfg.trait,ancestry=ancestry,contrast=label,method=a,subset=s,reference=b,reference_subset=t,metric=metric,
                                      delta=lookup[(ancestry,s,a)][metric]-lookup[(ancestry,t,b)][metric],lower=bounds[0],upper=bounds[1],valid_boot=len(samples),requested_boot=B,
                                      direction=METRIC_DIRECTIONS[metric],grouping_status=cfg.grouping_status,
                                      inference='conditional_on_frozen_models_and_thresholds',multiplicity='primary contrast prespecified; others exploratory',
                                      status='interval_estimated' if enough else 'point_only_or_insufficient_bootstrap'))
    confidence_rows=[];confidence=prediction['diagnostics']['pre_prediction_confidence']
    for ancestry in scopes:
        am=np.ones(n,bool) if ancestry=='ALL' else pop==ancestry
        for lo in range(0,100,10):
            mask=base_scope&am&(confidence>=lo)&((confidence<lo+10) if lo<90 else confidence<=100)
            if mask.sum()<10:continue
            for name in (primary,full,model.fallback):
                row=dict(trait=cfg.trait,ancestry=ancestry,method=name,confidence_lo=lo,confidence_hi=lo+10,N=int(mask.sum()),
                         mean_confidence=float(confidence[mask].mean()),matching_supported_fraction=float(prediction['diagnostics']['matching_supported'][mask].mean()),
                         **metrics(y[mask],p[name][mask],binary,model.training_variance))
                confidence_rows.append(row)
    state_metric='Brier' if binary else 'MSE'
    tab=pd.DataFrame(contrasts);check=tab[(tab.ancestry=='ALL')&(tab.metric==state_metric)]
    selected_known=masks['selected'];enough_info=selected_known.sum()>=30
    if binary:enough_info &= min(np.sum(y[selected_known]==v) for v in (0,1))>=10
    def supported(label):
        d=check[check.contrast==label]
        return bool(enough_info and len(d)==1 and d.iloc[0].valid_boot>=200 and d.iloc[0].upper<0)
    status=dict(trait=cfg.trait,task=cfg.trait_type,evaluation_scope=prediction['evaluation_scope'],
                original_N=n,known_outcome_N=int(known.sum()),technical_QC_excluded=int((~qc).sum()),
                selected_N=int(selected.sum()),selected_with_observed_outcome_N=int(selected_known.sum()),actual_policy_coverage=float(selected.mean()),
                selection_loss_supported=supported('selection_effect'),training_loss_supported=supported('training_effect_same_subset'),
                beyond_confidence_loss_supported=supported('beyond_baseline_confidence'),primary_loss_metric=state_metric,
                grouping_status=cfg.grouping_status,clinical_release=False,external_validation=False,
                dataset_has_been_inspected=True,roster_hash=id_hash(ids),frozen_hash=prediction['frozen_hash'])
    flow=pd.DataFrame(dict(reason=prediction['reason'],outcome_observed=known)).groupby(['reason','outcome_observed'],dropna=False).size().reset_index(name='N')
    tables=dict(metrics=pd.DataFrame(rows),contrasts=tab,confidence_validation=pd.DataFrame(confidence_rows),
                cohort_flow=flow,status=pd.DataFrame([status]),comparator_availability=pd.DataFrame(availability))
    for d in tables.values():
        d['run_id']=prediction['frozen_hash'];d['trait']=cfg.trait;d['grouping_status']=cfg.grouping_status
    return tables,status

# Aggregate figures

def draw_report(tables,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    m=tables['metrics'];c=tables['confidence_validation'];d=tables['contrasts']
    score='Brier' if 'Brier' in m else 'MSE'
    # One chart per file, with no fixed palette.
    allrows=m[(m.ancestry=='ALL')&m['subset'].isin(['all','selected','rejected'])]
    fig,ax=plt.subplots(figsize=(10,6))
    for name,z in allrows.groupby('method',sort=False):
        z=z.set_index('subset').reindex(['all','selected','rejected']);ax.plot(z.index,z[score],marker='o',label=name)
    ax.set(xlabel='Identical evaluation subsets',ylabel=score+' (lower is better)',title='Full / selected / rejected — same-person model comparison')
    ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'same_subset.png',dpi=170);plt.close(fig)
    if len(c):
        fig,ax=plt.subplots(figsize=(9,5))
        for name,z in c[c.ancestry=='ALL'].groupby('method',sort=False):
            z=z.sort_values('mean_confidence');ax.plot(z.mean_confidence,z[score],marker='o',label=name)
        ax.set(xlabel='Pre-outcome confidence (reference support, 0–100)',ylabel=score,title='Confidence versus held-out error')
        ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'confidence_error.png',dpi=170);plt.close(fig)
    selected=d[(d.ancestry=='ALL')&(d.metric==score)&d.contrast.isin(['selection_effect','training_effect_same_subset','combined_effect','beyond_baseline_confidence','beyond_risk_ancestry_composition'])]
    fig,ax=plt.subplots(figsize=(10,5))
    for i,row in enumerate(selected.itertuples()):
        if np.isfinite(row.lower) and np.isfinite(row.upper) and row.lower<=row.delta<=row.upper:
            ax.errorbar(row.delta,i,xerr=[[row.delta-row.lower],[row.upper-row.delta]],fmt='o',capsize=3)
        else:
            ax.plot(row.delta,i,marker='o')
            if np.isfinite(row.lower) and np.isfinite(row.upper):ax.plot([row.lower,row.upper],[i,i])
    ax.axvline(0,linestyle='--');ax.set_yticks(np.arange(len(selected)),selected.contrast)
    ax.set(xlabel='Paired Δ'+score+' (negative favours the candidate)',title='Selection and training are different effects')
    fig.tight_layout();fig.savefig(out/'effect_decomposition.png',dpi=170);plt.close(fig)
    q=m[(m.ancestry=='ALL')&m['subset'].str.startswith('support_quality:')]
    fig,ax=plt.subplots(figsize=(9,5))
    for name,z in q.groupby('method',sort=False):
        z=z.sort_values('actual_coverage');ax.plot(z.actual_coverage,z[score],marker='o',label=name)
    ax.set(xlabel='Actual coverage of the original test roster',ylabel=score,title='Frozen twin-confidence coverage')
    ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'coverage.png',dpi=170);plt.close(fig)

# Analysis commands

def config_file(path,trait=None):
    d=json.loads(Path(path).read_text()) if path else {}
    for name in ('role_fractions','coverages','ridge_grid','logistic_c_grid'):
        if name in d:d[name]=tuple(d[name])
    cfg=Config(**d)
    if trait:cfg=replace(cfg,trait=trait,trait_type='binary' if trait=='t2dm' else 'continuous')
    cfg.validate();return cfg


def npz(path):
    with np.load(path,allow_pickle=False) as z:return {k:z[k] for k in z.files}


def export_projection(p,path):
    d=pd.DataFrame(dict(eid=p['eid'],group=p['group'],selected=p['selected'],reason=p['reason'],**p['diagnostics']))
    for name,v in p['predictions'].items():d['prediction::'+name]=v
    for key in ('selected_support','selected_policy','selected_predictability'):
        if key in p:d[key]=p[key]
    for name,v in p.get('selector_scores_final',{}).items():d['selector::'+name]=v
    d['expected_squared_error']=p['expected_squared_error'];d['model_hash']=p['frozen_hash']
    d.to_csv(path,sep='\t',index=False,float_format='%.17g',compression='gzip')


def save_fit(model,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);private=out/'_LOCAL_ONLY';private.mkdir(mode=0o700,exist_ok=True)
    joblib.dump(model,private/'model.joblib',compress=3)
    digest=state_hash(model)
    abm=getattr(model,'grid_abm',None)
    json_write(out/'MODEL_FROZEN.json',dict(model_hash=digest,config=asdict(model.config),test_outcomes_used=False,clinical_release=False,
        core_support_selector=model.config.selector,final_predictability_selector=abm.config.predictability_selector if abm is not None else None,
        abm_version=GRID_ABM_VERSION,abm_routes_config=asdict(abm.config) if abm is not None else None,
        abm_routes_audit=abm.audit if abm is not None else None))
    pd.DataFrame(model.donors).to_csv(private/'donors.tsv.gz',sep='\t',index=False,float_format='%.17g')
    feature_rows=[]
    for pos,original in enumerate(model.prep.keep):
        is_genetic=original in set(model.prep.profile_keep)
        gp=list(model.prep.profile_keep).index(original) if is_genetic else None
        feature_rows.append(dict(feature=model.features[original],role='genetic_profile' if is_genetic else 'clinical_covariate',training_mean=model.prep.mean[pos],training_sd=model.prep.scale[pos],imputation_median=model.prep.median[pos],distance_weight=model.geometry.weights[gp] if is_genetic else 0.,full_model_coefficient=model.predictors['ABM_full'].coefficients[pos]))
    tables=dict(training_summary=pd.DataFrame([model.training_summary]),roles=pd.DataFrame(model.role_summary),tuning=pd.DataFrame(model.tuning),OOF_audit=pd.DataFrame(model.oof_audit),profile_features=pd.DataFrame(feature_rows))
    if abm is not None:
        tables['route_tuning']=pd.DataFrame(abm.tuning)
        tables['abm_training_audit']=pd.DataFrame([{k:v for k,v in abm.audit.items() if not isinstance(v,(list,dict,tuple))}])
    write_aggregate(out/'aggregate',tables,dict(trait=model.config.trait,module='Twin training',run_id=digest))


def fit_cmd(args):
    d=npz(args.data);cfg=config_file(args.config)
    needed={'X','y','ids','groups','roles','features','genetic_indices','covariate_indices','csx_indices'}
    if not needed<=set(d) or 'test' in set(d['roles']):raise ValueError('Fit accepts development-only arrays')
    out=Path(args.output)
    if out.exists():raise FileExistsError('Choose a fresh fit directory; old models are never overwritten')
    if abm_routes_config().publication_mode:raise ValueError('publication_mode requires the full abm entry and comparator preflight')
    model=TwinModel(cfg).fit(**{k:d[k] for k in needed})
    model=fit_abm_routes(model,d,abm_routes_config())
    save_fit(model,out);write_source_provenance(out,ROOT,abm_routes_config())
    print(json.dumps(dict(model=str(out/'_LOCAL_ONLY/model.joblib'),hash=state_hash(model),abm=GRID_ABM_VERSION)))


def project_cmd(args):
    model=joblib.load(args.model);d=npz(args.data)
    if set(d)!= {'X','ids','groups','features'}:raise ValueError('Project accepts only predictors/identity, never query outcomes')
    p=model.project(**d);p=apply_abm_routes(model,p,d['X'],d['ids'],d['groups']);out=Path(args.output)
    if out.exists():raise FileExistsError(out)
    out.parent.mkdir(parents=True,exist_ok=True);joblib.dump(p,out,compress=3)
    export_projection(p,out.with_suffix('.confidence.tsv.gz'));print(json.dumps(dict(query_N=len(p['eid']),selected_N=int(p['selected'].sum()))))


def benchmark_manifest(path,ids,truth,trait,endpoint):
    if not path:return {},[]
    path=Path(path);spec=json.loads(path.read_text());items=spec['methods'];external={};audit=[]
    for item in items:
        name=item['name'];source=path.parent/item['path']
        if item.get('trait')!=trait or item.get('endpoint_definition')!=endpoint or item.get('value_type')!='outcome_scale_prediction':
            raise ValueError('Comparator endpoint/units manifest mismatch: '+name)
        if item.get('test_roster_sha256')!=id_hash(ids):raise ValueError('Comparator does not declare this exact frozen test roster: '+name)
        if not item.get('predictions_frozen_without_test_labels',False):raise ValueError('Comparator must declare out-of-training prediction provenance')
        if not source.exists():audit.append(dict(method=name,status='file_unavailable'));continue
        if 'sha256' not in item or sha256(source)!=item['sha256']:raise ValueError('Comparator file checksum missing/mismatch: '+name)
        d=identified(read_table(source))
        if set(d.eid)!=set(ids):audit.append(dict(method=name,status='different_or_incomplete_test_roster'));continue
        d=d.set_index('eid').loc[ids];p=pd.to_numeric(d[item['prediction_column']],errors='raise').to_numpy(dtype=float,na_value=np.nan)
        if name in external:raise ValueError('Duplicate benchmark method')
        if item.get('outcome_column'):
            yy=pd.to_numeric(d[item['outcome_column']],errors='raise').to_numpy(dtype=float,na_value=np.nan)
            if not np.allclose(yy,truth.set_index('eid').loc[ids,'y'],equal_nan=True):raise ValueError('Comparator outcome definitions/units differ')
        external[name]=p;audit.append(dict(method=name,status='provider_declared_frozen_prediction_verified_file_and_roster',source_sha256=item['sha256'],training_budget=item.get('training_budget','not_recorded')))
    return external,audit


def evaluate_cmd(args):
    model=joblib.load(args.model);p=joblib.load(args.predictions);truth=identified(read_table(args.truth))
    external,availability=benchmark_manifest(args.comparators,p['eid'],truth,model.config.trait,args.endpoint)
    tables,status=evaluate(model,p,truth,external=external,bootstrap=args.bootstrap)
    if 'selected_predictability' in p:tables.update(evaluate_abm_routes(model,p,truth,bootstrap=args.bootstrap))
    if availability:tables['comparator_provenance']=pd.DataFrame(availability)
    out=Path(args.output)
    if out.exists():raise FileExistsError('Choose a new report directory')
    write_aggregate(out,tables,dict(trait=model.config.trait,module='Twin evaluation',run_id=p['frozen_hash']))
    json_write(out/'status.json',status)
    draw_report(tables,out/'figures');print(json.dumps(status))


def publication_preflight(cfg,route_cfg,raw,query,prsformer,comparator,out_root):
    """Fail before fitting when publication prerequisites are not verifiable."""
    if not route_cfg.publication_mode:return
    from twin_core import NumericProfile,FrozenRule
    import tempfile
    failures=[]
    if cfg.grouping_status not in {'family_verified','kinship_verified','verified_family','verified_groups'}:
        failures.append('relatedness/grouping not verified')
    if set(raw['groups'])&set(query['groups']):failures.append('group crosses outer split')
    if cfg.upstream_gwas_overlap not in {'verified_nonoverlap','verified_no_overlap'}:
        failures.append('GWAS/target non-overlap not explicitly verified')
    if not comparator or comparator.get('blockers'):failures.append('comparator provenance not verified')
    if not prsformer.get('reused') or prsformer.get('test_roster_hash')!=id_hash(query['ids']):
        failures.append('PRSformer/ABM test hash mismatch')
    if cfg.validation_scope not in {'exploratory_previously_inspected_test','external_untouched','new_untouched'}:
        failures.append('validation scope not explicitly declared')
    try:
        prep=NumericProfile().fit(raw['X'][raw['roles']=='build'],raw['genetic_indices'],raw['covariate_indices'],cfg)
        cal=raw['roles']=='calibration_fit';ok=prep.qc(raw['X'][cal]);ids=raw['ids'][cal]
        if ok.sum()<cfg.min_role_n:raise ValueError('Insufficient calibration QC rows')
        # Pre-fit feasibility check uses identities, never audit/test outcomes.
        score=np.array([hash_number(x,cfg.seed)/2**64 for x in ids])
        for q in route_cfg.coverages:FrozenRule().fit(score,ok,ids,q,cfg.seed).apply(score,ok,ids)
    except (ValueError,IndexError) as exc:failures.append('six coverage masks unavailable: '+str(exc))
    try:
        Path(out_root).mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='.provenance-check-',dir=out_root) as dest:
            payload=write_source_provenance(dest,ROOT,route_cfg)
            if len(payload['files_sha256'])<8:raise ValueError('missing core source files')
    except (OSError,ValueError) as exc:failures.append('source provenance unavailable: '+str(exc))
    if failures:raise ValueError('publication_mode preflight failed: '+'; '.join(failures))


def run(args):
    traits=args.traits.split(',')
    if not traits or len(set(traits))!=len(traits) or not set(traits)<= {'height','ldl','t2dm'}:raise ValueError('Use unique height,ldl,t2dm traits')
    route_cfg=abm_routes_config()
    prepared_path=args.prepared or str(Path(args.cache_dir)/'prepared.joblib')
    if getattr(args,'prepare',False):
        if args.data_file:raise ValueError('--prepare uses frozen comparison inputs; do not combine with --data-file')
        from abm_prepare import prepare_abm
        prepare_abm(args)
        return
    if not args.data_file and not args.check:
        from abm_prepare import prepare_abm
        prepare_abm(args)
    require_comparators=getattr(args,'require_comparators',False) or route_cfg.publication_mode
    if args.analysis_root and require_comparators:
        command=[sys.executable,str(ROOT/'f/check_comparators.py'),'--analysis-root',args.analysis_root,
                 '--prsformer-upstream',os.getenv('GRID_PRSFORMER_UPSTREAM','/mnt/f/software/PRSformer'),'--strict']
        receipt=subprocess.run(command,capture_output=True,text=True)
        if receipt.returncode:
            raise ValueError('Published comparator preflight failed: '+(receipt.stdout or receipt.stderr).strip())
        comparator_check=json.loads(receipt.stdout)
        print('[GRID preflight] Published comparator provenance passed',flush=True)
    else:comparator_check=dict(status='deferred_until_comparison',blockers=[])
    # Validate every requested trait before fitting the first one. These input
    # contexts keep evaluation truth separate from the development-only fit API.
    contexts=[]
    for trait in traits:
        cfg=config_file(args.config,trait)
        d,groups,meta=load_source(args.prepared or str(Path(args.cache_dir)/'prepared.joblib'),trait,
                                   args.data_file.replace('{trait}',trait) if args.data_file else None,
                                   args.schema.replace('{trait}',trait) if args.schema else None)
        raw,query,truth,flow,cfg,extra=numeric_cohort(d,groups,cfg,args.roles_file)
        endpoint=meta.get('endpoint_definition') or meta.get('prepared_settings',{}).get(trait+'_col',trait)
        if trait=='t2dm' and endpoint=='auto':endpoint='baseline_t2dm_from_Yr2e_Yt2e'
        endpoint=args.endpoint or endpoint
        prsformer_reuse=verify_existing_prsformer_roster(args.analysis_root,d,query['ids']) if args.analysis_root else dict(status='analysis_root_unavailable',reused=False)
        if args.analysis_root and require_comparators:
            cov=meta.get('prepared_settings',{}).get('covariates','')
            cov=cov.split(',') if isinstance(cov,str) and cov else list(groups.get('covariates',[]))
            _,comparison=verified_prsformer(args.analysis_root,trait,d,truth,cfg,cov,endpoint)
            if comparison['status']!='same_roster_endpoint_and_covariates_verified':
                raise ValueError('PRSformer comparison preflight failed: '+comparison['status'])
            prsformer_reuse['provenance']=comparison
        publication_preflight(cfg,route_cfg,raw,query,prsformer_reuse,comparator_check,args.out_root)
        print(json.dumps(dict(trait=trait,input='validated',endpoint=endpoint,prsformer_reuse=prsformer_reuse,**flow)),flush=True)
        if not args.check:contexts.append((trait,cfg,d,groups,meta,raw,query,truth,flow,extra,endpoint,prsformer_reuse))
    if args.check:return
    from abm_resume import run_trait
    for context in contexts:
        run_trait(sys.modules[__name__], args, route_cfg, context, comparator_check)
        if args.analysis_root:
            from workflow_status import readiness
            readiness(args.analysis_root,args.traits,out_root=args.out_root,prepared=prepared_path,show_abm=False)
    if args.analysis_root and Path(args.out_root).resolve().is_relative_to(Path(args.analysis_root).resolve()):
        print('[GRID] Updating the final evidence index for completed ABM results',flush=True)
        subprocess.run([sys.executable,str(ROOT/'f/browser_index.py'),'--analysis-root',args.analysis_root],check=True)


def demo(args):
    root=Path(args.output)
    if root.exists():raise FileExistsError(root)
    root.mkdir(parents=True);rng=np.random.default_rng(219)
    n=2400;ids=np.asarray([f'SYNTHETIC_{i}' for i in range(n)]);pop=np.where(np.arange(n)%4,'EUR','OTHER')
    g=rng.normal(size=(n,6));age=rng.normal(size=n);sex=rng.integers(0,2,n)
    # Heteroscedasticity exercises confidence evaluation; no tests require an ABM win.
    y=170+4*g[:,0]+.5*g[:,1]+2*sex+rng.normal(size=n)*(3+np.abs(g[:,2]))
    if args.binary:y=rng.binomial(1,1/(1+np.exp(-(-1+.6*g[:,0]+.2*age))))
    d=pd.DataFrame(g,columns=[f'csx_{i}' for i in range(6)]);d['age']=age;d['sex']=sex
    d['eid']=ids;d['family_id']=np.asarray([f'F_{i//2}' for i in range(n)]);d['ancestry']=pop;d['y']=y
    d['split']=np.where(np.arange(n)<n//2,'train','test');d.loc[n//2:n//2+4,[f'csx_{i}' for i in range(6)]]=np.nan
    data=root/'synthetic.tsv.gz';d.to_csv(data,sep='\t',index=False)
    schema=root/'schema.json';json_write(schema,dict(feature_groups=dict(csx=[f'csx_{i}' for i in range(6)],covariates=['age','sex']),metadata=dict(synthetic_only=True)))
    cfg=replace(Config(),folds=3,repeats=1,min_role_n=12,min_class_n=3,k=12,min_twins=3,min_group_ess=2,gate_iterations=8,bootstrap=25,ridge_grid=(10.,100.),logistic_c_grid=(.1,))
    config=root/'config.json';json_write(config,asdict(cfg))
    trait='t2dm' if args.binary else 'height'
    run(argparse.Namespace(traits=trait,prepared=None,cache_dir='unused',data_file=str(data),schema=str(schema),roles_file=None,config=str(config),check=False,out_root=str(root/'abm'),analysis_root=None,run_id='synthetic',endpoint='SIMULATED_ONLY',comparators=None,bootstrap=25))
    (root/'SYNTHETIC_ONLY.txt').write_text('Synthetic execution check, not UK Biobank evidence or performance validation.\n')


def analysis_main(argv=None):
    a=argparse.ArgumentParser(description='GRID twin: run/fit/project/evaluate/demo with separated outcome access.');sp=a.add_subparsers(dest='command',required=True)
    r=sp.add_parser('abm');r.add_argument('--traits',default='height,ldl,t2dm');r.add_argument('--cache-dir',default='/mnt/d/analysis/grid/cache/abm');r.add_argument('--prepared');r.add_argument('--data-file');r.add_argument('--schema');r.add_argument('--roles-file')
    r.add_argument('--require-comparators',action='store_true',help='Require completed, verified PRSformer comparisons before fitting (also enabled by publication_mode)')
    r.add_argument('--prepare',action='store_true',help='Prepare ABM features from the frozen preparation roster, without waiting for PRSformer results')
    r.add_argument('--prepare-config',help='JSON input path/covariate overrides for ABM preparation')
    r.add_argument('--config',default=str(DEFAULT_CONFIG),help='Configuration JSON (default: config/twin.json beside grid.sh)');r.add_argument('--out-root',default='/mnt/d/analysis/grid/abm');r.add_argument('--analysis-root',default='/mnt/d/analysis/grid');r.add_argument('--run-id');r.add_argument('--new-run',action='store_true',help='Explicitly fit a new ABM run instead of resuming verified work');r.add_argument('--check',action='store_true');r.add_argument('--endpoint');r.add_argument('--comparators');r.add_argument('--bootstrap',type=int);r.add_argument('--threads',type=int,default=4,help='Numeric worker threads (default: 4)')
    f=sp.add_parser('fit');f.add_argument('--data',required=True);f.add_argument('--config',required=True);f.add_argument('--output',required=True);f.add_argument('--threads',type=int,default=4)
    p=sp.add_parser('project');p.add_argument('--model',required=True);p.add_argument('--data',required=True);p.add_argument('--output',required=True);p.add_argument('--threads',type=int,default=4)
    e=sp.add_parser('evaluate');e.add_argument('--model',required=True);e.add_argument('--predictions',required=True);e.add_argument('--truth',required=True);e.add_argument('--output',required=True);e.add_argument('--comparators');e.add_argument('--endpoint',required=True);e.add_argument('--bootstrap',type=int);e.add_argument('--threads',type=int,default=4)
    d=sp.add_parser('demo');d.add_argument('--output',required=True);d.add_argument('--binary',action='store_true');d.add_argument('--threads',type=int,default=2)
    args=a.parse_args(argv)
    if args.threads<1:a.error('--threads must be positive')
    action=dict(abm=run,fit=fit_cmd,project=project_cmd,evaluate=evaluate_cmd,demo=demo)[args.command]
    os.umask(0o077)
    with threadpool_limits(args.threads):action(args)

# Public command dispatch

HELP='''Modules:
    prsformer   PRSformer benchmark
    prsformer2  PRSformer with separate prediction calibration
    abm         Twin-based ABM analysis (Clinical, CSx-linear, MoE, BlockNet)
    final       Final comparison and result index
    cleanup     Preview unused large caches; add --apply to release
    shiny       Result viewer

Usage:
    ./grid.sh [abm] --traits height,ldl,t2dm
    ./grid.sh module1,module2 --traits height,ldl
    Default module: abm; workflow (or all) runs prsformer -> prsformer2 -> abm -> final
    Comma-separated modules run in order with shared options; stop on failure.

Options:
    --traits CSV          [height,ldl,t2dm]
    --module-only         Accepted; abm is already standalone
    --predictor NAME      [blocknet] or moe
    --representation NAME [csx_ldblocks_fixed] or csx_scores
    --device NAME         [auto], cpu or cuda
    --config FILE_OR_DIR  [config/abm.json]
    --analysis-root DIR   [/mnt/d/analysis/grid]
    --threads N           [4]
    --dry-run             Preview without running
    <module> --help       Module options
    -h, --help            Help

Example:
    # Run both benchmarks sequentially:
    ./grid.sh prsformer,prsformer2 --traits height,ldl
    # Run modules one at a time, in this order:
    ./grid.sh prsformer --traits height,ldl,t2dm
    ./grid.sh prsformer2 --traits height,ldl,t2dm
    ./grid.sh abm --module-only --traits height,ldl,t2dm
    ./grid.sh final --traits height,ldl,t2dm
    ./grid.sh shiny
'''


def r_runtime(explicit=None):
    """Find an R with the viewer dependencies, keeping system/Conda libs separate."""
    override=explicit or os.getenv('GRID_RSCRIPT') or os.getenv('R_BIN')
    candidates=[override] if override else [shutil.which('Rscript'),'/usr/bin/Rscript']
    probe="p<-c('shiny','DT','DBI','RSQLite','ggplot2','jsonlite');quit(status=if(all(vapply(p,requireNamespace,logical(1),quietly=TRUE)))0L else 1L)"
    for candidate in dict.fromkeys(x for x in candidates if x):
        command=[candidate]
        resolved=shutil.which(candidate)
        if resolved and Path(resolved).resolve()==Path('/usr/bin/Rscript').resolve():
            command=['env','-u','R_ENVIRON_USER','-u','R_LIBS_USER',candidate]
        try:
            if subprocess.run([*command,'--vanilla','-e',probe],capture_output=True).returncode==0:return command
        except OSError:
            pass
    raise RuntimeError('No usable Shiny R runtime. Install R_PACKAGES.txt dependencies or set GRID_RSCRIPT.')


def main(argv=None):
    args=list(sys.argv[1:] if argv is None else argv)
    if args and args[0] in ('-h','--help','help'):print(HELP);return
    command = args.pop(0) if args and not args[0].startswith('-') else 'abm'
    allowed = {'prsformer', 'prsformer2', 'abm', 'workflow', 'all', 'abm-assess',
               'final', 'shiny', 'cleanup', 'self-check'}
    modules = [name.strip() for name in command.split(',')]
    # Validate the entire list before starting any potentially expensive module.
    for module in modules:
        if module not in allowed:
            raise ValueError('Unknown module: ' + (module or '<empty>') + '; see ./grid.sh --help')
    if len(modules) > 1:
        for index, module in enumerate(modules, 1):
            print(f'[GRID] MODULE {index}/{len(modules)} {module}', flush=True)
            # Keep module runtimes isolated, just like consecutive shell commands.
            # Forward --dry-run and --help so each child handles its own preview.
            subprocess.run([sys.executable, str(ROOT / 'f/twin.py'), module, *args], check=True)
        return
    command = modules[0]
    if command=='self-check':
        from self_check import main as self_check_main
        self_check_main(args);return
    if command == 'cleanup':
        from cache_cleanup import main as cleanup_main
        cleanup_main(args)
        return
    dry='--dry-run' in args;args=[x for x in args if x!='--dry-run'];env=None
    quiet=any(x in ('-h','--help') for x in args)
    if command in ('prsformer','prsformer2','abm','workflow','all') and not dry and not quiet:
        # Parse first: invalid options must not create caches or run self-checks.
        if command in ('abm', 'workflow', 'all'):
            from abm_core.pipeline import parser
            preflight = parser()
        else:
            from abm_core.cold_start import benchmark_parser
            preflight = benchmark_parser(command, command)
            if command == 'prsformer2':
                preflight.add_argument('--check', action='store_true')
                preflight.add_argument('--upstream-dir')
        options = preflight.parse_args(args)
        from abm_core.cold_start import benchmark_config
        benchmark_config(options)
        os.environ.setdefault('GRID_CACHE_ROOT', str(Path(options.analysis_root) / 'cache'))
        from automatic_check import ensure_self_check
        ensure_self_check()
    if command=='prsformer':
        from prsformer_resume import dispatch
        dispatch(ROOT, args, dry=dry)
        return
    elif command=='abm':
        from abm_core.pipeline import main as abm_main
        abm_main([*args,*(['--dry-run'] if dry else [])])
        return
    elif command in ('workflow', 'all'):
        from abm_core.workflow import main as workflow_main
        workflow_main([*args,*(['--dry-run'] if dry else [])])
        return
    elif command=='prsformer2':
        from abm_core.prsformer_adapter import main as prsformer_main
        prsformer_main([*args,*(['--dry-run'] if dry else [])])
        return
    elif command=='abm-assess':
        from abm_core.deploy import main as abm_assess_main
        if dry:raise ValueError("Use abm-assess --help to inspect its input contract")
        abm_assess_main(args)
        return
    elif command=='final':cmd=[sys.executable,str(ROOT/'f/final_dispatch.py'),*args]
    elif command=='shiny':
        import argparse
        p=argparse.ArgumentParser();p.add_argument('--analysis-root',default='/mnt/d/analysis/grid');p.add_argument('--db');p.add_argument('--traits',default='height,ldl,t2dm');p.add_argument('--port',type=int,default=3839);p.add_argument('--host',default='127.0.0.1');p.add_argument('--r-bin');p.add_argument('--review',action='store_true');p.add_argument('--reindex',action='store_true')
        a=p.parse_args(args);db=Path(a.db or str(Path(a.analysis_root)/'final/evidence.sqlite'))
        if not 1024<=a.port<=65535:raise ValueError('Invalid port')
        if a.host not in ('127.0.0.1','::1','localhost') and os.getenv('GRID_SHINY_ALLOW_NETWORK')!='TRUE':raise ValueError('Network binding requires GRID_SHINY_ALLOW_NETWORK=TRUE and your own authentication')
        if a.reindex or not db.is_file():
            ix=[sys.executable,str(ROOT/'f/browser_index.py'),'--analysis-root',a.analysis_root,'--out',str(db.parent),'--traits',a.traits]
            print(shlex.join(ix),flush=True)
            if not dry:subprocess.run(ix,check=True)
        runtime=[a.r_bin or os.getenv('GRID_RSCRIPT') or os.getenv('R_BIN','Rscript')] if dry else r_runtime(a.r_bin)
        cmd=[*runtime,str(ROOT/'shiny/app.R'),'--db',str(db),'--port',str(a.port),'--host',a.host,*(['--review'] if a.review else [])]
    else:raise ValueError('Unknown command: '+command)
    print('[GRID] '+shlex.join(cmd),flush=True)
    if not dry:subprocess.run(cmd,check=True,env=env)

if __name__=='__main__':
    try:main()
    except subprocess.CalledProcessError as exc:
        print('ERROR: '+str(exc),file=sys.stderr);sys.exit(exc.returncode if exc.returncode>0 else 128-exc.returncode)
    except (ValueError,OSError,RuntimeError) as exc:
        print('ERROR: '+str(exc),file=sys.stderr);sys.exit(2)

