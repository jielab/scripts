import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd

from abm_core.model import ABMAssessor
from abm_core import reporting, pipeline
from abm_core.losses import row_loss
from abm_core.heads import CountHeads
from abm_core.support import count_twins
from abm_core.contracts import ordered_ids_sha256


class LossAndContrastTests(unittest.TestCase):
    def test_logloss_fit_unavailable_roundtrip_observed_units(self):
        for name,expected in [('logloss',.916290731874155),('brier',.34)]:
            with self.subTest(loss=name),tempfile.TemporaryDirectory() as td:
                frames=[]
                for role in ['gate','head','policy','test']:
                    frames.append(pd.DataFrame(dict(id=[role+str(i) for i in range(40)],
                        group_id=[role+str(i) for i in range(40)],y=np.tile([0.,1.],20),
                        p_base=.2,p_candidate=.4,n_twins=40,supported=True,ancestry='EUR')))
                gate,head,policy,test=frames
                manifest=dict(baseline_id='base',candidate_id='candidate',trait='t2dm',**{
                    k:'a'*64 for k in ('bank_sha256','geometry_sha256','predictor_bundle_sha256','role_manifest_sha256')})
                model=ABMAssessor.fit(gate,head,policy,manifest=manifest,task='binary',loss_name=name,
                    head_options=dict(n_bins=1,n_boot=30,min_calibration_groups=2,prior_strength=0))
                p=Path(td)/'assessor.json';model.save(p);model=ABMAssessor.load(p)
                assessed=model.assess(test.drop(columns='y'),manifest=manifest).merge(test[['id','y']],on='id')
                cfg=dict(head_options=dict(seed=4,min_calibration_groups=2),evaluation=dict(bootstrap_replicates=30))
                table=reporting.evaluation_calibration(assessed,model,cfg,'test')
                row=table[(table.stratum=='all')&(table.objective=='expected_loss')].iloc[0]
                self.assertAlmostEqual(row.observed_estimate,expected,places=14)
                self.assertAlmostEqual(row.expected_mean,expected,places=14)
                self.assertAlmostEqual(row.bias_estimate,0,places=14)
                diag=reporting.fit_diagnostics(gate,head,[], 'binary',policy.drop(columns='y'),loss_name=name)
                self.assertEqual(reporting.diagnostic_loss(diag),name)
                unavailable=ABMAssessor.fit(gate.assign(supported=False),head.assign(supported=False),policy,
                    manifest=manifest,task='binary',loss_name=name)
                unavailable.save(p);self.assertEqual(ABMAssessor.load(p).heads.loss_name,name)
                reporting.freeze_matches(assessed.drop(columns='y'),Path(td),model.heads.seed)
                self.assertEqual(CountHeads.unavailable('binary',loss_name=name).loss_name,name)
        for task,loss in [('continuous','brier'),('binary','mse'),('other','mse')]:
            with self.assertRaises(ValueError):CountHeads.unavailable(task,loss_name=loss)
        with self.assertRaisesRegex(ValueError,'unknown schema'):
            reporting.diagnostic_loss(dict(task='binary',models={}))

    def test_paired_count_contrasts_equal_predictions_and_unequal_groups(self):
        # Target=(y-p_base)^2=1. A predicts 0, C=D predict 1: D-A=1, D-C=0.
        d=pd.DataFrame(dict(id=list('abcdef'),group_id=['g0','g1','g1','g2','g2','g2'],
                            y=1.,p_base=0.,p_candidate=1.,assessment_available=True))
        for key,pred in [('A',0.),('B',.5),('C',1.),('D',1.)]:
            for obj in ('expected_loss','expected_gain'):
                d[f'diagnostic_{key}_{obj}']=pred;d[f'diagnostic_{key}_{obj}_selected']=True
        bundle=dict(schema_version='abm.diagnostics.1',task='continuous',loss_name='mse',models={'frozen':True})
        cfg=dict(evaluation=dict(bootstrap_replicates=400),head_options=dict(seed=19))
        result=reporting.count_incremental_contrasts(bundle,d,cfg)
        self.assertEqual(len(result),4)
        zero=result[result.contrast=='D_vs_C']
        np.testing.assert_array_equal(zero[['mse_reduction','ci_low','ci_high']].to_numpy(),np.zeros((2,3)))
        self.assertEqual(((d.y-d.p_base)**2-d.diagnostic_A_expected_loss).pow(2).mean(),1.)
        # Per-person delta = [1,2,2,4,4,4], original group sizes 1,2,3.
        delta=np.array([1.,2.,2.,4.,4.,4.]);d['diagnostic_A_expected_loss']=1+np.sqrt(delta)
        d['diagnostic_B_expected_loss']=1.
        row=reporting.count_incremental_contrasts(bundle,d,cfg).iloc[0]
        self.assertAlmostEqual(row.mse_reduction,17/6)
        rng=np.random.default_rng(19);draws=[]
        for _ in range(400):
            counts=np.bincount(rng.integers(0,3,3),minlength=3)
            draws.append(counts@np.array([1.,4.,12.])/(counts@np.array([1.,2.,3.])))
        np.testing.assert_allclose([row.ci_low,row.ci_high],np.quantile(draws,[.025,.975]))
        d['diagnostic_D_expected_loss_selected']=False
        r=reporting.count_incremental_contrasts(bundle,d,cfg)
        self.assertTrue(pd.isna(r.iloc[1].selection_difference))
        r=reporting.count_incremental_contrasts(bundle,d.assign(group_id='one'),cfg)
        self.assertTrue(r.mse_reduction.isna().all())

    def test_support_cache_refreshes_evidence_without_recomputing_distances(self):
        rng=np.random.default_rng(9);n=80
        frame=pd.DataFrame(dict(id=[f'i{i}' for i in range(n)],group_id=[f'g{i}' for i in range(n)],
            role=['build']*40+['tune_model']*20+['test']*20,y=0.))
        cfg= json.loads((pipeline.ROOT/'config/abm.json').read_text())
        cfg['support'].update(radius_k_ref=5,min_count=1,min_group_ess=1)
        x=rng.normal(size=(n,4,1));rows=np.arange(n)
        with tempfile.TemporaryDirectory() as td:
            a,saved=pipeline.prepare_support(frame,x,rows,['scores']*4,{},Path(td),cfg,'cpu',False)
            with patch.object(pipeline,'count_twins',side_effect=AssertionError('must reuse')):
                b,reused=pipeline.prepare_support(frame.assign(y=100),x,rows,['scores']*4,{},Path(td),cfg,'cuda',True)
            np.testing.assert_array_equal(a.n_twins,b.n_twins)
            self.assertTrue(b.groups_verified.all());self.assertEqual(reused['executed_device'],'cpu')
            self.assertTrue(reused['cache_reused']);self.assertTrue(reused['bank']['groups_verified'])
            self.assertTrue(b.loc[b.supported,'support_reason'].eq('supported_verified_groups').all())

    def test_cuda_boundary_ties_groups_and_frozen_policy(self):
        import torch
        if not torch.cuda.is_available():self.skipTest('No CUDA: execution parity not verified')
        donors=np.array([[1.,0.],[-1.,0.],[0.,1.],[0.,-1.],[1.+1e-8,0.],[1.-1e-8,0.],[0.,0.]])
        queries=np.array([[0.,0.],[0.,0.]])
        opts=dict(donor_ids=['z','b','a','m','x','y','self'],query_ids=['q','self'],
            donor_groups=['g0','g1','g2','g3','g4','g5','g1'],query_groups=['g1','same'],
            radius=1.,k_used=3,min_count=2,min_group_ess=2,query_batch=1,donor_batch=2)
        cpu=count_twins(donors,queries,device='cpu',**opts)
        gpu=count_twins(donors,queries,device='cuda',**opts)
        pd.testing.assert_frame_equal(cpu.table.drop(columns='device_used'),gpu.table.drop(columns='device_used'),rtol=1e-10,atol=1e-10)
        np.testing.assert_array_equal(cpu.neighbor_indices,gpu.neighbor_indices)
        # The identical support enters the same frozen count-only heads/rules.
        from test_model import make_role,manifest
        gate,head,policy=[make_role(p,20) for p in ("gate","head","policy")]
        for f in (gate,head,policy):f["n_twins"]=f.n_twins.map({50:4,100:5})
        model=ABMAssessor.fit(gate,head,policy,manifest=manifest(),head_options={'n_boot':30,'n_bins':2,'min_calibration_groups':2})
        for table in (cpu.table,gpu.table):
            table['p_base']=.2;table['p_candidate']=.1;table['group_id']=['qg','selfg']
        a=model.assess(cpu.table,manifest=manifest());b=model.assess(gpu.table,manifest=manifest())
        np.testing.assert_array_equal(a.selected,b.selected);np.testing.assert_array_equal(a.p_policy,b.p_policy)
        self.assertTrue(a.selected.any());self.assertTrue((~a.selected).any())
