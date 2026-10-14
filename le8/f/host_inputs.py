"""Prepare current ABM inputs with fixed outcome definitions and fresh family roles."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]


def host_config(args):
    import shlex
    import cohort_inputs as host
    values=['--Y',args.Y,'--biom',args.biom,'--seed',str(args.seed),'--r-bin',args.r_bin]
    for key,flag in [('ukb_phe','ukb-phe'),('group_file','group-file'),('group_col','group-col'),
                     ('outer_roster','outer-roster'),('end_date','end-date'),
                     ('diagnosis_col','diagnosis-col'),('shared_covariates','covariates')]:
        value=getattr(args,key,None)
        if value is not None:
            values += ['--'+flag,str(value)]
    # Scope and random seed are owned by the public invocation.
    extra=shlex.split(args.host_args)
    if any(t.split('=')[0] in {'--Y','--biom','--seed'} for t in extra):
        raise ValueError('Set outcome, layer and seed using the public options')
    a=host.configure(host.parser().parse_args(values+extra))
    a.cores=args.cores
    return host,a


def check_group_roles(groups,roles):
    groups=np.asarray(groups,str);roles=np.asarray(roles,str)
    order=np.argsort(groups,kind='stable')
    if np.any((groups[order][1:]==groups[order][:-1])&(roles[order][1:]!=roles[order][:-1])):
        raise ValueError('Roles divide a family; check the supplied grouping and roster')


def clinical_matrix(p, columns, categorical, training):
    """Build-only vocabulary, no label-derived preprocessing or rare-level pooling.

    Numeric imputation/scale and constant-column removal remain fold-local inside
    twin. Dummy vocabulary is schema metadata; held-out levels stay explicitly missing.
    Technical factors are fixed adjustment inputs, outside the molecular budget.
    """
    arrays=[];names=[];audit=[]
    for col in columns:
        s=p[col]
        if col in categorical or not pd.api.types.is_numeric_dtype(s):
            values=s.astype('string');levels=sorted(values[training].dropna().unique())
            if not levels:
                raise ValueError('No training levels for '+col)
            valid=values.isin(levels).to_numpy()
            for level in levels:
                x=values.eq(level).fillna(False).to_numpy(dtype=float)
                x[~valid]=np.nan;arrays.append(x);names.append('clinical::'+col+'='+str(level))
            audit.append(dict(column=col,encoding='training-vocabulary one-hot',levels=list(map(str,levels)),
                              unknown_or_missing_N=int((~valid).sum())))
        else:
            arrays.append(pd.to_numeric(s,errors='raise').to_numpy(float));names.append('clinical::'+col)
            audit.append(dict(column=col,encoding='numeric'))
    if not arrays:
        raise ValueError('Clinical covariates required')
    return np.column_stack(arrays),names,audit


def prepare_host(args, work):
    host,a=host_config(args)
    print(f'[LE8] START ABM inputs: biom={a.biom}; omics={a.omics_file}; '
          f'met_input={a.met_input}; transform={a.transform}; '
          f'covariates={a.covariates}; extra_adjustment={a.residualize or "none"}',flush=True)
    prepared=work/'host_inputs_LOCAL_ONLY';prepared.mkdir(parents=True,mode=0o700)
    p=host.prepare(a,prepared)
    raw=np.load(prepared/'raw.npy');features=(prepared/'features.txt').read_text().splitlines()
    p,outcome_audit=host.outcomes(p,a)
    eligible=p.eligible.to_numpy(bool);eligible_p=p.loc[eligible].reset_index(drop=True)
    roles=host.assign_roles(eligible_p,a)
    groups=host.family_groups(eligible_p,a)
    check_group_roles(groups,roles)
    ids=eligible_p[a.id_col].astype(str).to_numpy()
    columns=list(dict.fromkeys(host.words(a.covariates)+host.words(a.residualize)))
    C,cn,encoding=clinical_matrix(eligible_p,columns,host.words(a.categorical),roles=='build')
    x=raw[eligible]
    if a.transform=='log1p':
        x=np.log1p(x)
    # Deviation preserves prevalent participants and does not exclude future cases.
    from run_deviation import fixed_roles
    all_groups=host.family_groups(p,a)
    dr=fixed_roles(all_groups,a.seed)
    mapping={'build':'train','tune_model':'tune','tune_gate':'tune',
             'calibration_fit':'calibrate','calibration_audit':'calibrate','test':'test'}
    inherited=dict(zip(groups,[mapping[r] for r in roles]))
    dr=np.asarray([inherited.get(g,r) for g,r in zip(all_groups,dr)])
    check_group_roles(all_groups,dr)
    dc,dcn,de=clinical_matrix(p,columns,host.words(a.categorical),dr=='train')
    baseline=host.parsed_date(p[a.baseline_col],a.baseline_col)
    known=p.prevalent.astype(float).to_numpy()
    known[baseline.isna().to_numpy()]=np.nan
    if a.disease_evidence_col:
        unknown=(pd.to_numeric(p[a.disease_evidence_col],errors='raise')>0)&p[a.diagnosis_col].isna()
        known[unknown.to_numpy()]=np.nan
    endpoint=pd.DataFrame({'eid':p[a.id_col].astype(str),'known_cad_at_index':known,'index_date':baseline})
    if args.objective_file:
        from research_io import id_table
        extra=id_table(args.objective_file).reset_index(drop=True)
        # A supplied baseline status takes precedence; future-event columns never do.
        fields=[c for c in extra if c in {'eid','known_cad_at_index','index_date','objective_date',
                 'objective_cad','objective_modality','intervening_cad','intervening_treatment'}]
        for col in set(fields)&set(endpoint)-{'eid'}:
            endpoint=endpoint.drop(columns=col)
        endpoint=endpoint.merge(extra[fields],on='eid',how='left',validate='one_to_one')
    audit=dict(outcome=outcome_audit,split_source='explicit_outer_roster' if a.outer_roster else 'fresh_outcome_blind_family_assignment',
               roster_sha256=hashlib.sha256(json.dumps(list(zip(ids.tolist(),roles.tolist())),separators=(',',':')).encode()).hexdigest(),
               clinical_encoding=encoding,
               deviation_encoding=de,clinical_columns=columns,technical_adjustment='fixed covariate inputs; no full-cohort residualization',
               family_status='verified_supplied_groups' if a.group_col else 'unavailable_individual_ID_groups_only',
               independent_family_support_verified=bool(a.group_col),
               validation_status='exploratory',input_config=vars(a))
    return dict(X=np.c_[x,C],time=eligible_p.time.to_numpy(),event=eligible_p.event.to_numpy(),
                eid=ids,family_group=groups,role=roles,feature_names=features+cn,clinical_feature_names=cn), \
           dict(M=raw,C=dc,covariate_names=dcn,features=features,endpoint=endpoint,
                eid=np.asarray(p[a.id_col],str),group=np.asarray(all_groups,str),role=np.asarray(dr,str)),audit
