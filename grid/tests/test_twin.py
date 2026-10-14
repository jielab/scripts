from pathlib import Path
import sys,inspect
from dataclasses import replace
import numpy as np,pandas as pd,pytest,joblib
from threadpoolctl import threadpool_limits
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'f'))
from twin_core import *
from twin_evaluate import metrics,evaluate

@pytest.fixture(scope='module')
def fitted():
    rng=np.random.default_rng(482);n=1100;X=rng.normal(size=(n,14));y=3*X[:,0]-X[:,1]+rng.normal(size=n)*2
    roles=np.repeat(ROLE_NAMES,[650,150,100,100,100]);ids=np.asarray([f'd{i}' for i in range(n)])
    X[:20,:10]=np.nan
    cfg=Config(folds=3,repeats=2,k=12,min_twins=3,min_group_ess=2,gate_iterations=6,min_role_n=10,min_class_n=2,bootstrap=3,ridge_grid=(10.,100.))
    with threadpool_limits(1):m=TwinModel(cfg).fit(X,y,ids,ids,roles,[f'x{i}' for i in range(14)],list(range(10)),list(range(10,14)),[0,1])
    q=rng.normal(size=(80,14));q[:10,:10]=np.nan;qi=np.asarray([f'q{i}' for i in range(80)])
    p=m.project(q,qi,qi,m.features)
    return m,X,y,ids,roles,q,qi,p


def test_all_missing_genetics_abstain_keep_denominator(fitted):
    m,*_,q,ids,p=fitted
    assert len(p['eid'])==80
    assert not p['diagnostics']['technical_qc'][:10].any()
    assert not p['diagnostics']['pre_prediction_confidence'][:10].any()
    assert not p['selected'][:10].any()
    assert np.isnan(p['predictions']['ABM_policy'][:10]).all()


def test_missing_donors_never_enter_banks(fitted):
    m,*_=fitted
    for bank in m.banks.values():assert not set(bank.ids)&{f'd{i}' for i in range(20)}
    assert not m.donors['retained'][:20].any()
    assert not m.donors['random_retained'][:20].any()
    assert not m.donors['soft_weight'][:20].any()
    assert m.donors['technical_excluded'][:20].all()
    assert not m.donors['profile_discordance_excluded'][:20].any()


def test_dummy_count_cannot_dilute_genetic_qc():
    rng=np.random.default_rng(2);g=rng.normal(size=(100,10));q=g[:2].copy();q[0]=np.nan;q[1,:3]=np.nan
    for n in (1,50,150):
        p=NumericProfile().fit(np.c_[g,rng.normal(size=(100,n))],range(10),range(10,10+n),Config())
        query=np.c_[q,np.zeros((2,n))]
        np.testing.assert_array_equal(p.missing(query)[0],[1,.3]);assert not p.qc(query).any()


def test_projection_is_batch_independent(fitted):
    m,*_,q,ids,p=fitted
    for i in (0,10,50):
        one=m.project(q[i:i+1],ids[i:i+1],ids[i:i+1],m.features)
        for name,v in p['predictions'].items():np.testing.assert_allclose(one['predictions'][name],v[i:i+1],atol=1e-11)
        for name,v in p['diagnostics'].items():np.testing.assert_allclose(one['diagnostics'][name],v[i:i+1],atol=1e-11)
        assert one['selected'][0]==p['selected'][i]


def test_query_outcome_rejected(fitted):
    m,*_,q,ids,p=fitted
    assert not {'y','event','time','outcome'}&set(inspect.signature(m.project).parameters)
    with pytest.raises(TypeError):m.project(q,ids,ids,m.features,y=np.zeros(len(ids)))


def test_original_schema_and_overlap_guards(fitted):
    m,*_,q,ids,p=fitted
    with pytest.raises(ValueError,match='Feature names/order'):m.project(q,ids,ids,m.features[::-1])
    with pytest.raises(ValueError,match='overlaps'):m.project(q[:1],m.ids[:1],m.groups[:1],m.features)
    with pytest.raises(ValueError,match='Audit'):m.project(q[:1],ids[:1],ids[:1],m.features,audit=True)


def test_saved_model_roundtrip(fitted,tmp_path):
    m,*_,q,ids,p=fitted;f=tmp_path/'m.joblib';joblib.dump(m,f);loaded=joblib.load(f)
    pred=loaded.project(q,ids,ids,loaded.features)
    np.testing.assert_allclose(pred['diagnostics']['pre_prediction_confidence'],p['diagnostics']['pre_prediction_confidence'])
    assert state_hash(m)==state_hash(loaded)


def test_twin_profile_not_scalar_risk():
    z=np.array([[1.,-1.],[-1.,1.]]);geo=Geometry().fit(z,np.ones(2),Config())
    assert np.all(z@np.ones(2)==0)
    assert np.linalg.norm(geo.transform(z)[0]-geo.transform(z)[1])>0


def test_distant_query_has_no_twin_but_fallback(fitted):
    m,*_,q,ids,p=fitted;xx=q[20:21].copy();xx[:,:10]=10000
    out=m.project(xx,ids[:1],ids[:1],m.features)
    assert out['diagnostics']['technical_qc'][0]
    assert not out['diagnostics']['matching_supported'][0]
    assert out['diagnostics']['pre_prediction_confidence'][0]==0
    assert out['predictions']['ABM_policy'][0]==out['predictions'][m.fallback][0]


def test_same_group_not_a_twin():
    cfg=Config(k=3,min_twins=1,min_group_ess=1,min_role_n=2)
    b=ReferenceBank().fit(np.arange(8.).reshape(4,2),np.array(['a','b','c','d']),np.array(['F','F','F','G']),np.zeros(4),cfg);b.radius=100
    d=b.query(np.array([[0.,1.]]),np.array(['q']),np.array(['F']))
    assert d['matched_count'][0]==1


def test_pruned_and_random_composition(fitted):
    m,*_=fitted;d=m.donors
    assert d['retained'].sum()==d['random_retained'].sum()
    assert np.all(d['retained']<=d['donor_eligible'])
    assert np.all(d['random_retained']<=d['donor_eligible'])


def test_frozen_rule_never_forces_unsupported():
    ids=np.array([str(i) for i in range(10)]);s=np.arange(10.)
    r=FrozenRule().fit(s,np.zeros(10,bool),ids,1,8)
    assert not r.apply(s,np.zeros(10,bool),ids).any()


def test_dosage_qc_not_faked_from_summary():
    x=np.random.default_rng(2).uniform(0,2,(100,5));p=NumericProfile().fit(x,range(5),[],replace(Config(),profile_kind='variant_dosage'))
    q=x[:1].copy();q[0,0]=3
    assert not p.qc(q)[0]


def test_r2_is_not_squared_correlation():
    y=np.arange(10.);p=-y
    met=metrics(y,p,False,np.var(y))
    assert met['prediction_R2']<0
    assert np.corrcoef(y,p)[0,1]**2>.999


def test_eval_separates_effects_and_counts(fitted):
    m,*_,q,ids,p=fitted;y=np.arange(len(ids))*.02
    with threadpool_limits(1):tables,status=evaluate(m,p,pd.DataFrame({'eid':ids,'y':y}),bootstrap=3)
    assert status['original_N']==80 and status['technical_QC_excluded']==10
    assert {'selection_effect','training_effect_same_subset','combined_effect'}<=set(tables['contrasts'].contrast)
    assert not tables['contrasts'].lower.notna().any() # 3 draws cannot produce a credible CI
    assert not status['clinical_release'] and not status['external_validation']
    assert not status['training_loss_supported']


def test_incomplete_external_not_silently_intersected(fitted):
    m,*_,q,ids,p=fitted;v=np.zeros(80);v[35]=np.nan
    tables,status=evaluate(m,p,pd.DataFrame({'eid':ids,'y':np.arange(80)}),external={'external':v},bootstrap=0)
    assert 'external' not in set(tables['metrics'].method)
    assert tables['comparator_availability'].iloc[0].status=='unavailable_incomplete_frozen_roster'


def test_truth_roster_mismatch_rejected(fitted):
    m,*_,q,ids,p=fitted
    with pytest.raises(ValueError,match='roster'):evaluate(m,p,pd.DataFrame({'eid':ids[:-1],'y':np.zeros(79)}),bootstrap=0)


def test_roles_do_not_use_y():
    ids=np.asarray([f'family{i//2}' for i in range(1000)]);a=assign_roles(ids,Config());b=assign_roles(ids,Config())
    np.testing.assert_array_equal(a,b)
    assert all(a[i]==a[i+1] for i in range(0,1000,2))


def test_audit_labels_do_not_fit_model():
    rng=np.random.default_rng(55);n=700;X=rng.normal(size=(n,6));y=X[:,0]+rng.normal(size=n)
    ids=np.array([f'a{i}' for i in range(n)]);role=np.repeat(ROLE_NAMES,[400,100,80,60,60]);cfg=Config(folds=2,repeats=1,k=8,min_twins=2,min_group_ess=1,min_role_n=10,gate_iterations=3,ridge_grid=(10.,))
    with threadpool_limits(1):
        a=TwinModel(cfg).fit(X,y,ids,ids,role,[str(i) for i in range(6)],list(range(4)),[4,5],[0])
        y[role=='audit']+=1e6
        b=TwinModel(cfg).fit(X,y,ids,ids,role,[str(i) for i in range(6)],list(range(4)),[4,5],[0])
    assert state_hash(a)==state_hash(b)


def test_binary_path_preserves_probability_scale():
    rng=np.random.default_rng(77);n=900;X=rng.normal(size=(n,6));y=rng.binomial(1,expit(X[:,0]))
    ids=np.array([f'b{i}' for i in range(n)]);role=np.repeat(ROLE_NAMES,[500,120,100,90,90]);cfg=Config(trait='t2dm',trait_type='binary',folds=2,repeats=1,min_class_n=2,min_role_n=8,k=8,min_twins=2,min_group_ess=1,gate_iterations=3,logistic_c_grid=(.1,))
    with threadpool_limits(1):m=TwinModel(cfg).fit(X,y,ids,ids,role,[str(i) for i in range(6)],range(4),[4,5],[0])
    qi=np.array(['bq1','bq2']);p=m.project(X[:2],qi,qi,m.features)
    for v in p['predictions'].values():assert np.all((v>=0)&(v<=1))


@pytest.mark.parametrize('binary',[False,True])
def test_cached_bootstrap_metrics_equal_reference(binary):
    from twin_evaluate import FixedMetricCache
    rng=np.random.default_rng(891);n=120
    y=rng.binomial(1,.3,n).astype(float) if binary else rng.normal(size=n)
    p=rng.choice([.1,.3,.7,.9],n) if binary else rng.normal(size=n)
    cache=FixedMetricCache(y,p,binary,2.)
    for _ in range(5):
        w=rng.integers(0,4,n).astype(float)
        a=cache(w);b=metrics(y,p,binary,2.,w)
        for key in a:np.testing.assert_allclose(a[key],b[key],rtol=1e-12,atol=1e-12)
