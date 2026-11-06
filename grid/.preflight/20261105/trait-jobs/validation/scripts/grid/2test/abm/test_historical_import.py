import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import joblib
import numpy as np
import pandas as pd
from import_legacy_comparators import import_comparators
from abm_core.model import file_sha256,write_json
from twin_core import id_hash,state_hash


def fixture(root,order=None):
    run=root/'source';local=run/'_LOCAL_ONLY';local.mkdir(parents=True)
    roster=np.array(['001','002','003','004']);order=np.arange(4) if order is None else np.array(order)
    joblib.dump(dict(eid=roster,frozen_hash='frozen'),local/'predictions.joblib')
    joblib.dump(dict(eid=roster[order],predictions={'Disco_reference':np.array([.2,.3,.8,.7])[order]}),local/'comparator_predictions.joblib')
    pd.DataFrame(dict(eid=roster,y=[0,0,1,1],ancestry='EUR')).iloc[order].to_csv(local/'test_truth.tsv.gz',sep='\t',index=False)
    features=np.array(['csx.AFR','csx.EAS','csx.EUR','csx.SAS','age','sex','PC1','PC2'])
    np.savez(local/'test_X.npz',ids=roster[order],groups=roster[order],X=np.ones((4,8)),features=features)
    ids=np.array(['d1','d2','d3']);roles=np.array(['build','tune_model','calibration_fit'])
    np.savez(local/'development.npz',ids=ids,groups=ids,roles=roles,features=features)
    pd.DataFrame(dict(eid=ids,group=ids,role=roles)).to_csv(local/'development_roles.tsv.gz',sep='\t',index=False)
    joblib.dump(dict(columns=['age','sex','PC1','PC2','disco'],imputer=None,scale=None,
                     model=SimpleNamespace(),calibration=SimpleNamespace()),local/'disco_calibration.joblib')
    write_json(run/'protocol.json',dict(abm_version='GRID-ABM-FINAL-2026-10-17',config=dict(trait='t2dm',trait_type='binary'),
        endpoint_definition='baseline_t2dm_from_Yr2e_Yt2e'))
    write_json(run/'SOURCE_FROZEN.json',dict(format='GRID-source-provenance-final'))
    write_json(run/'COMPLETE.json',dict(status='completed',test_roster_hash=id_hash(roster),model_hash='frozen'))
    seal(run)
    return run


def seal(run):
    write_json(run/'RESUME.json',dict(stages={'exports':{str(p.relative_to(run)):file_sha256(p)
        for p in run.rglob('*') if p.is_file() and p.name!='RESUME.json'}}))


class HistoricalImportTests(unittest.TestCase):
    def test_eight_feature_baseline_stays_historical(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);run=fixture(root)
            with np.load(run/'_LOCAL_ONLY/test_X.npz') as q:features=q['features'].tolist()
            model=SimpleNamespace(features=features,prep=SimpleNamespace(keep=np.arange(8)),
                columns={'CSx_reference':np.arange(8)},calibration={'CSx_reference':SimpleNamespace()})
            path=run/'fit/_LOCAL_ONLY/model.joblib';path.parent.mkdir(parents=True);joblib.dump(model,path)
            p=joblib.load(run/'_LOCAL_ONLY/predictions.joblib');p['frozen_hash']=state_hash(model)
            p['predictions']={'CSx_reference':np.array([.1,.2,.9,.8])}
            joblib.dump(p,run/'_LOCAL_ONLY/predictions.joblib')
            complete=json.loads((run/'COMPLETE.json').read_text());complete['model_hash']=state_hash(model)
            write_json(run/'COMPLETE.json',complete);seal(run)
            rows,methods=import_comparators(root,run,root/'historical/import')
            row=methods[methods.old_name.eq('CSx_reference')].iloc[0]
            self.assertEqual(row.n_features,8);self.assertEqual(row.status,'available')
            self.assertFalse((root/'abm').exists())

    def test_read_only_alignment_and_actual_disco_identity(self):
        with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
            left,right=Path(a),Path(b);ra=fixture(left);rb=fixture(right,[3,1,0,2])
            before={str(p):file_sha256(p) for p in ra.rglob('*') if p.is_file()}
            with patch('sklearn.linear_model.LogisticRegression.fit',side_effect=AssertionError('never fit')):
                da,ma=import_comparators(left,ra,left/'historical/import')
                db,mb=import_comparators(right,rb,right/'historical/import')
            self.assertAlmostEqual(da.brier.item(),.065)
            self.assertEqual(da.brier.item(),db.brier.item())
            self.assertEqual(ma.method.item(),'csx-disco-legacy/raw-interpolation')
            self.assertEqual(ma.n_features.item(),5)
            self.assertEqual(before,{str(p):file_sha256(p) for p in ra.rglob('*') if p.is_file()})

    def test_duplicate_length_and_roster_fail(self):
        for change in ('duplicate','length','roster'):
            with self.subTest(change=change),tempfile.TemporaryDirectory() as td:
                root=Path(td);run=fixture(root);p=run/'_LOCAL_ONLY/comparator_predictions.joblib';d=joblib.load(p)
                if change=='duplicate':d['eid'][0]=d['eid'][1]
                elif change=='length':d['predictions']['Disco_reference']=np.array([.2])
                else:d['eid'][0]='other'
                joblib.dump(d,p);seal(run)
                with self.assertRaises(ValueError):import_comparators(root,run,root/'out')

    def test_missing_artifacts_report_paths_without_training(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);_,methods=import_comparators(root,'missing',root/'out')
            self.assertEqual(methods.status.item(),'unavailable_missing_artifacts')
            self.assertIn(str(root/'missing/_LOCAL_ONLY/comparator_predictions.joblib'),methods.missing.item())
