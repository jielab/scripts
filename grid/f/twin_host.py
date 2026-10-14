"""Reuse existing GRID prepared feature matrices; never rerun GWAS or rewrite splits."""
from __future__ import annotations
from dataclasses import replace
from pathlib import Path
import json,sys,importlib.util
import joblib,numpy as np,pandas as pd
from twin_core import validate_ids,assign_roles,assert_roles,id_hash
from twin_io import identified,read_table,sha256


def load_source(prepared_path,trait,data_file=None,schema_file=None):
    if data_file:
        if not schema_file:raise ValueError('Explicit external features require --schema JSON')
        spec=json.loads(Path(schema_file).read_text());d=identified(read_table(data_file))
        groups=spec['feature_groups'];metadata=spec.get('metadata',{})
        metadata.update(input_source=str(data_file),input_sha256=sha256(data_file),matrix_provenance='provider_declared')
    else:
        path=Path(prepared_path)
        if not path.is_file():raise FileNotFoundError(str(path)+'; reuse the existing prepared.joblib, or explicitly run the preserved upstream prepare/evolution stages first')
        # Only a trusted local file produced by the user's GRID pipeline. Never load
        # arbitrary downloaded pickle/joblib files through the public viewer.
        prepared=joblib.load(path)
        if not {'tables','feature_groups','cohort'}<=set(prepared):raise ValueError('Unrecognized GRID prepared cache')
        if prepared.get('state') not in (None,'features_ready'):raise ValueError('Prepared cache has not completed genetic feature preparation')
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
                raise ValueError('Prepared source changed; explicitly verify/rebuild features, preserving the outer roster: '+str(source))
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
    # Distances exclude covariates; a PC used as both adjustment and genetic input
    # is one numeric column with its explicit genetic identity, not a duplicate.
    clinical=[c for c in clinical if c not in genetic]
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
    extra=dict(d=d,fit=fit,test=test,role=roles,disco_columns=comparator_extra)
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
