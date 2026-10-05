"""Supplied S7 numerical contracts executed directly against production kernels."""
import copy
import importlib.util
from pathlib import Path
import pickle
import os
import sys
import numpy as np
import pandas as pd
import pytest

path = Path(__file__).resolve().parents[1] / 'f' / 'c1.abm.py'
if 'c1_abm' in sys.modules:
    m = sys.modules['c1_abm']
else:
    spec = importlib.util.spec_from_file_location('c1_abm', path)
    m = importlib.util.module_from_spec(spec); sys.modules[spec.name] = m; spec.loader.exec_module(m)


def test_rare_continuous_signal_survives():
    x = np.zeros((10000, 1)); x[:30, 0] = 20.
    p = m.S7TailTransform().fit(x)
    z = p.transform(x)
    assert np.var(z['bulk']) == 0
    assert np.var(z['tail_up']) > 0
    assert np.all(z['tail_up'][:30] > 0) and np.all(z['tail_up'][30:] == 0)
    assert p.audit_.rare_signal_preserved_in_tail.iloc[0]


def test_rare_binary_not_winsorized():
    x = np.zeros((10000, 1)); x[:30] = 1
    p = m.S7TailTransform().fit(x); z = p.transform(x)
    np.testing.assert_array_equal(x, z['bulk'])
    assert p.feature_types_ == ('binary',)


def test_explicit_dosage_preserved():
    x = np.array([[0], [.1], [1.], [1.8], [2.]])
    p = m.S7TailTransform().fit(x, ['dosage'])
    np.testing.assert_allclose(p.transform(x)['bulk'][:, 0], x[:, 0]/2)


def test_training_constant_and_missing_dropped():
    x = np.c_[np.arange(10.), np.ones(10), np.full(10, np.nan)]
    p = m.S7TailTransform().fit(x)
    assert p.keep_.tolist() == [0]
    assert p.audit_.status.tolist() == ['retained','constant_raw','all_missing']


@pytest.mark.parametrize('x', [np.zeros((4, 1)), np.full((4, 1), np.nan)])
def test_no_usable_features_fails(x):
    with pytest.raises(ValueError): m.S7TailTransform().fit(x)


def test_projection_cannot_refit():
    rng = np.random.default_rng(5); x = rng.normal(size=(500, 3))
    p = m.S7TailTransform().fit(x); before = pickle.dumps(p)
    p.transform(np.full((5, 3), 1e10))
    assert before == pickle.dumps(p)


def test_transform_row_batch_invariance():
    x = np.random.default_rng(4).normal(size=(100, 4))
    p = m.S7TailTransform().fit(x[:60]); a = p.transform(x[60:])
    b = p.transform(x[60:][::-1])
    for name in a: np.testing.assert_array_equal(a[name], b[name][::-1])
    for name in a: np.testing.assert_array_equal(a[name][:1], p.transform(x[60:61])[name])


def test_missing_is_not_extreme():
    x = np.arange(100.)[:, None]; p = m.S7TailTransform().fit(x)
    z = p.transform(np.array([[np.nan], [np.inf], [-np.inf]]))
    assert not z['observed'].any()
    assert not z['tail_up'].any() and not z['tail_down'].any()


def test_lower_and_upper_tails_have_separate_identity():
    p = m.S7TailTransform().fit(np.arange(100.)[:, None])
    z = p.transform(np.array([[-100.], [200.]]))
    assert z['tail_down'][0, 0] > 0 and z['tail_up'][0, 0] == 0
    assert z['tail_up'][1, 0] > 0 and z['tail_down'][1, 0] == 0


def test_numerical_tail_cap_is_flagged():
    p = m.S7TailTransform().fit(np.arange(100.)[:, None])
    z = p.transform(np.array([[-1e300], [1e300]]))
    assert z['tail_cap_hit'].all()
    assert np.isfinite(z['tail_up']).all() and np.isfinite(z['tail_down']).all()


def test_linear_comparators_get_same_information():
    x = np.arange(80.).reshape(20, 4)
    p = m.S7TailTransform().fit(x)
    assert p.linear_design(x).shape == (20, 16)


def test_wrong_source_feature_count_rejected():
    p = m.S7TailTransform().fit(np.arange(10.)[:,None])
    with pytest.raises(ValueError): p.transform(np.ones((4, 2)))


def test_binary_schema_violation_fails():
    p = m.S7TailTransform().fit(np.array([[0.],[1.]]))
    with pytest.raises(ValueError): p.transform(np.array([[2.]]))


def borrow(d, **kw):
    d = np.array(d, dtype=float)
    return m.fixed_scale_borrow(d, np.ones_like(d), np.ones_like(d), np.full(len(d), .2),
                               radius=kw.pop('radius', 1.), **kw)


def test_fixed_kernel_can_focus_without_copy1():
    r = borrow([[0.] + [1.]*99], temperature=.25)
    assert r.weights[0,0] > .95
    assert r.borrow_fraction[0] < .05
    assert .2 < r.probability[0] < .3


def test_absolute_distance_reduces_support():
    near = borrow([[.1, .2, .3]])
    far = borrow([[10., 20., 30.]])
    assert far.support[0] < near.support[0]
    assert far.borrow_fraction[0] < near.borrow_fraction[0]
    assert abs(far.probability[0]-.2) < 1e-10


@pytest.mark.parametrize('mode', ['empty','masked','zero_weight','infinite'])
def test_no_donors_falls_back(mode):
    d = np.ones((2, 0 if mode == 'empty' else 3))
    w = np.zeros_like(d) if mode == 'zero_weight' else np.ones_like(d)
    if mode == 'infinite': d[:] = np.inf
    allowed = np.zeros_like(d, dtype=bool) if mode == 'masked' else None
    r = m.fixed_scale_borrow(d, np.ones_like(d), w, np.array([.1,.7]), radius=1, allowed=allowed)
    np.testing.assert_array_equal(r.probability, [.1,.7])
    assert not r.borrow_fraction.any() and not r.weights.any()


def test_related_donors_do_not_inflate_group_ess():
    r = borrow([[0.,0.,0.,0.]], donor_groups=np.array([['A','A','B','B']]))
    np.testing.assert_allclose(r.effective_donors, [4.])
    np.testing.assert_allclose(r.effective_groups, [2.])
    assert r.borrow_fraction[0] < borrow([[0.,0.,0.,0.]]).borrow_fraction[0]


def test_donor_order_invariance():
    d = np.array([[.1,.8,.3]]); y = np.array([[1.,0.,1.]]); w = np.array([[2.,1.,3.]])
    a = m.fixed_scale_borrow(d,y,w,np.array([.3]),radius=.7)
    b = m.fixed_scale_borrow(d[:,::-1],y[:,::-1],w[:,::-1],np.array([.3]),radius=.7)
    np.testing.assert_allclose(a.probability,b.probability,rtol=0,atol=1e-15)
    np.testing.assert_allclose(a.weights,b.weights[:,::-1],rtol=0,atol=1e-15)


def test_masked_label_change_has_no_effect():
    d = np.array([[.1,.2]]); y = np.array([[1.,0.]])
    kw = dict(distances=d, donor_ipcw=np.ones_like(d), prior=np.array([.2]), radius=1., allowed=np.array([[False,True]]))
    a=m.fixed_scale_borrow(labels=y,**kw); y[0,0]=0; b=m.fixed_scale_borrow(labels=y,**kw)
    np.testing.assert_array_equal(a.probability,b.probability)


def test_invalid_distance_rejected():
    with pytest.raises(ValueError): borrow([[-1.,1.]])
    with pytest.raises(ValueError): borrow([[np.nan,1.]])


def metadata():
    return pd.DataFrame({'eid':['00001','00002'], 'age':[50,60]})


def test_projection_family_file_without_internal_column():
    p=metadata(); mapping=pd.DataFrame({'eid':['00002','00001','extra'], 'group':['002','001','003']})
    z=m.attach_projection_groups(p,mapping=mapping)
    assert z.eid.tolist() == ['00001','00002']
    assert z['.le8_family'].tolist() == ['001','002']


@pytest.mark.parametrize('mode', ['missing','duplicate','conflict'])
def test_bad_family_files_fail(mode):
    p=metadata(); mp=pd.DataFrame({'eid':['00001','00002'], 'group':['A','B']})
    kw={}
    if mode=='missing': mp=mp.iloc[:1]
    elif mode=='duplicate': mp=pd.concat([mp,mp.iloc[:1]])
    else: p['family']=['A','C']; kw['source_group_col']='family'
    with pytest.raises(ValueError): m.attach_projection_groups(p,mapping=mp,**kw)


def test_projection_requires_explicit_family_information():
    with pytest.raises(ValueError): m.attach_projection_groups(metadata())
    z=m.attach_projection_groups(metadata(),allow_individual_fallback=True)
    assert z['.le8_family'].tolist() == metadata().eid.tolist()


def test_noncontext_development_person_is_blocked():
    with pytest.raises(ValueError): m.assert_external_roster(['train_not_context'],['F9'],{'train_not_context'},{'F1'})


def test_development_family_is_blocked():
    with pytest.raises(ValueError): m.assert_external_roster(['new'],['F1'],{'old'},{'F1'})


def test_external_person_and_family_pass():
    m.assert_external_roster(['new'],['new_family'],{'old'},{'old_family'})


def test_nonstring_id_rejected_instead_of_losing_zeros():
    p=metadata();p.eid=[1,2]
    with pytest.raises(ValueError): m.attach_projection_groups(p,allow_individual_fallback=True)
