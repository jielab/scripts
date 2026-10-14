"""Boundary adapter: prepared numeric host arrays -> frozen TwinABM experiment.

This adapter does NOT reinterpret UKB endpoints, load RDS, build PRSs, change the
original split, or replace the production mixed-type/tail preprocessing. It is
called after the host has constructed an explicit numeric baseline-only schema.
See docs/HOST_INTEGRATION.md for the production integration contract.
"""
from dataclasses import asdict,replace
from pathlib import Path
import numpy as np
from abm_twin_core import Config,ROLES,json_write,fingerprint
from run_abm_twin import fit_command,write_result
from abm_twin_evaluate import evaluate


def run_prepared_numeric_cohort(*, X, time, event, eid, family_group, role,
                               feature_names, clinical_feature_names,
                               output, config:Config=Config(),
                               confirmatory_test=False):
    """Keep the existing roster; accept exactly five development roles plus test.

    clinical_feature_names must be explicit so a protein-only experiment cannot
    silently be described as a clinical-plus-protein baseline comparison.
    confirmatory_test is metadata ONLY; this function cannot verify that the
    user has never inspected that test's labels/results in earlier experiments.
    """
    out=Path(output)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError('Use a fresh versioned directory; never overwrite previous ABM results')
    X=np.asarray(X,float);ids=np.asarray(eid).astype(str);groups=np.asarray(family_group).astype(str)
    roles=np.asarray(role).astype(str);names=list(map(str,feature_names))
    if set(roles)!=set(ROLES)|{'test'}:raise ValueError('Missing/unknown roles; the adapter never creates a new test split')
    if not clinical_feature_names or not set(clinical_feature_names)<=set(names):
        raise ValueError('Declare included numeric baseline clinical features; no silently weaker comparator')
    if len(set(names))!=len(names):raise ValueError('Duplicate feature names')
    forbidden={'time','event','outcome','target','prevalent_status','role','split','diagnosis_date','t2e','Yt2e'}
    if forbidden & set(names) or any(name.lower().endswith(('.t2e','.yt2e','.bi2e')) for name in names):
        raise ValueError('Outcome or split field in X')
    config=replace(config,always_include_indices=tuple(names.index(str(n)) for n in clinical_feature_names))
    dev=roles!='test';test=~dev
    if set(groups[dev]) & set(groups[test]):raise ValueError('Family overlap with test')
    out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out/'development_LOCAL_ONLY.npz',X=X[dev],time=np.asarray(time)[dev],event=np.asarray(event)[dev],
        eid=ids[dev],group=groups[dev],role=roles[dev],feature_names=np.asarray(names))
    # Test outcomes are in a separate file and never passed to fit/project.
    np.savez_compressed(out/'test_truth_LOCAL_ONLY.npz',time=np.asarray(time)[test],event=np.asarray(event)[test],eid=ids[test])
    np.savez_compressed(out/'test_X_LOCAL_ONLY.npz',X=X[test],eid=ids[test],group=groups[test],feature_names=np.asarray(names))
    model=fit_command(out/'development_LOCAL_ONLY.npz',out/'fit',config)
    pred=model.predict(X[test],ids[test],groups[test],names)
    import joblib
    joblib.dump(pred,out/'test_predictions_LOCAL_ONLY.joblib',compress=3)
    import pandas as pd
    pd.DataFrame({'eid':pred['eid'],'selected':pred['selected'],'reason':pred['reason'],
                  'candidate_risk':pred['probability'][model.config.primary],
                  'fallback_risk':pred['probability']['elasticnet_full'],
                  'policy_risk':pred['probability']['fallback_policy'],**pred['diagnostics']}).to_csv(
                  out/'test_predictions_LOCAL_ONLY.confidence.tsv.gz',sep='\t',index=False)
    result=evaluate(model,pred,np.asarray(time)[test],np.asarray(event)[test],ids[test])
    result['status']['test_provenance']='user_declared_untouched' if confirmatory_test else 'previously_inspected_or_unverified_exploratory'
    result['status']['clinical_release']=False
    write_result(out/'test_report',result)
    json_write(out/'host_contract.json',dict(
        original_roster_preserved=True,clinical_feature_names=list(clinical_feature_names),
        numeric_schema=names,all_ids_hash=fingerprint(ids),test_ids_hash=fingerprint(ids[test]),
        UKB_phenotype_definition='supplied by host; not reconstructed by this adapter',
        preprocessing='numeric reference core; production categorical/technical/tail adapter requires validation',
        config=asdict(config),confirmed_independent_external_validation=False))
    return result
