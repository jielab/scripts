"""Public preparation entry points must not invent missing family/endpoint data."""
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'f'))
from abm_twin_core import ROLES


def prepare_command(tmp_path, mode, grouping=None, invalid_group=False):
    rng = np.random.default_rng(317)
    if mode == 'abm':
        roles = np.repeat([*ROLES, 'test'], 40)
    else:
        roles = np.repeat(['train', 'tune', 'calibrate', 'test'], [120, 40, 40, 40])
    n = len(roles)
    ids = np.array([f'person_{i:04d}' for i in range(n)])
    groups = np.array([f'family_{i//2:04d}' for i in range(n)])
    cohort = pd.DataFrame(dict(eid=ids, age=rng.normal(55, 4, n), sex=rng.integers(0, 2, n),
        time=np.where(np.arange(n)%3 == 0, 5., 15.), event=(np.arange(n)%3 == 0).astype(int),
        known_cad_at_index=np.zeros(n, int), role=roles))
    if grouping == 'family':
        cohort['family'] = groups
        if invalid_group:
            cohort.loc[0, 'family'] = np.nan
    G = rng.normal(size=(n, 3)); M = .5*G + rng.normal(size=(n, 3))
    omics = pd.DataFrame(M, columns=['A','B','C']); omics.insert(0,'eid',ids)
    cohort_file = tmp_path/'cohort.tsv'; cohort.to_csv(cohort_file,sep='\t',index=False)
    omics_file = tmp_path/'omics.tsv'; omics.to_csv(omics_file,sep='\t',index=False)
    config = tmp_path/'config.json'
    config.write_text(json.dumps({} if mode=='abm' else dict(reference_min_n=20,folds=3,bootstrap=1)))
    output = tmp_path/'prepared'
    command = ['bash',str(ROOT/'le8.sh'),mode,'prepare','--cohort-file',str(cohort_file),
        '--omics-file',str(omics_file),'--clinical','age,sex','--role-col','role',
        '--config',str(config),'--output',str(output)]
    if grouping == 'family': command += ['--group-col','family']
    if grouping == 'declared': command += ['--assume-unrelated']
    if mode == 'deviation':
        pgs = pd.DataFrame(G,columns=['A.SCORE_SUM','B.SCORE_SUM','C.SCORE_SUM']);pgs.insert(0,'eid',ids)
        pgs_file = tmp_path/'pgs.tsv'; pgs.to_csv(pgs_file,sep='\t',index=False)
        command += ['--pgs-root',str(pgs_file)]
    return command, output, ids, groups, roles


def run(command):
    return subprocess.run(command,cwd=ROOT,text=True,capture_output=True,timeout=60,
        env={**os.environ,'PYTHON':sys.executable,'PYTHONDONTWRITEBYTECODE':'1',
             'SCRIPT_VERBOSE':'1','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})


@pytest.mark.parametrize('mode',['abm','deviation'])
@pytest.mark.parametrize('grouping',[None,'declared','family'])
def test_prepare_optional_family_preserves_roster_and_unverified_status(tmp_path,mode,grouping):
    command,out,ids,groups,roles = prepare_command(tmp_path,mode,grouping)
    result = run(command)
    assert result.returncode == 0, result.stdout + result.stderr
    config = json.loads((out/'config.json').read_text())
    expected_unit = 'family' if grouping == 'family' else 'individual_unverified'
    assert config['grouping_unit'] == expected_unit
    provenance_file = 'input_provenance.json' if mode == 'abm' else 'preparation.json'
    audit = json.loads((out/provenance_file).read_text())
    assert audit['independent_family_support_verified'] == (grouping == 'family')
    assert audit['family_status'] == ('verified_supplied_groups' if grouping == 'family'
                                     else 'unavailable_individual_ID_groups_only')
    assert audit['assume_unrelated'] == (grouping == 'declared')
    if grouping is None: assert audit['group_assumption'] == 'individual_ID_only_unverified'
    dev = roles != 'test'
    with np.load(out/'_LOCAL_ONLY/development.npz',allow_pickle=False) as d:
        np.testing.assert_array_equal(d['eid'],ids[dev])
        np.testing.assert_array_equal(d['role'],roles[dev])
        np.testing.assert_array_equal(d['group'],(groups if grouping == 'family' else ids)[dev])
        if mode == 'deviation':
            assert not d['endpoint_eligible'].any() and np.isnan(d['y']).all()
            from molecular_deviation import DeviationConfig,DeviationStudy
            tr = d['role']=='train'
            model = DeviationStudy(DeviationConfig(**config)).fit_reference(
                d['M'][tr],d['G'][tr],d['C'][tr],d['eid'][tr],d['group'][tr],d['reference'][tr],
                d['feature_names'],d['covariate_names'])
            model.fit_disease_models({k:d[k] for k in d.files})
            assert model.classifiers == {} and len(model.reference_table)==3
            assert model.classification_status == 'objective_endpoint_unavailable'
    with np.load(out/'_LOCAL_ONLY/test_X.npz',allow_pickle=False) as q:
        np.testing.assert_array_equal(q['eid'],ids[~dev])
        assert not {'y','time','event'} & set(q.files)
        if mode == 'deviation':
            projected=model.project(q['M'],q['G'],q['C'],q['eid'],q['group'],q['feature_names'],q['covariate_names'])
            assert np.isfinite(projected['Z']).all() and projected['probability']=={}


@pytest.mark.parametrize('mode',['abm','deviation'])
def test_supplied_but_incomplete_family_is_not_silently_ignored(tmp_path,mode):
    command,*_ = prepare_command(tmp_path,mode,'family',invalid_group=True)
    result = run(command)
    assert result.returncode != 0
    assert 'Missing family' in result.stdout + result.stderr
