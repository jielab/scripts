import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT/'le8/f'), str(ROOT/'0f')]
import le8_scheduler as scheduler
import workflow_jobs as common


class Followups(unittest.TestCase):
    def test_labels_and_policy(self):
        for number in (1, 25):
            self.assertEqual(scheduler.job_label({'key':'dx.compare'},number,[]),f'le8.dates.compare.{number}')
            self.assertEqual(scheduler.job_label({'key':'timing.cvd_cad.prot'},number,[]),
                             f'le8.cvd_cad.prot.any1.timing.{number}')
        for key in ('timing.cvd_cad.prot','dx.compare'):
            self.assertEqual(scheduler.dependency_condition({'key':key},'any1.cvd_cad.prot.analysis'),'ended')
            self.assertEqual(scheduler.dependency_condition({'key':key},'aod.cvd_cad'),'done')
        self.assertEqual(scheduler.dependency_condition({'key':'any1.cvd_cad.prot.analysis'},'aod.cvd_cad'),'done')

    def test_failed_parent_worker_and_submission(self):
        with tempfile.TemporaryDirectory() as temp, patch.dict(os.environ, {'submit':'local','SUBMIT':'local'}):
            root=Path(temp); output=root/'compared'
            jobs=[common.Job('any1.cvd_cad.prot.primary',[sys.executable,'-c','raise SystemExit(7)'],[],
                             dict(cores=1,queue='33',gpus=0),{},temp),
                  common.Job('timing.cvd_cad.prot',[sys.executable,'-c','pass'],['any1.cvd_cad.prot.primary'],
                             dict(cores=1,queue='33',gpus=0),{},temp),
                  common.Job('dx.compare',[sys.executable,'-c',f'from pathlib import Path; Path({str(output)!r}).write_text("ok")'],
                             ['any1.cvd_cad.prot.primary','timing.cvd_cad.prot'],dict(cores=1,queue='33',gpus=0),{},temp)]
            plan=scheduler.write_plan(jobs,root/'jobs',mode='local',name='le8',analysis_root=root,modules=[])
            self.assertEqual(scheduler.worker(plan,jobs[0].key),7)
            self.assertEqual(scheduler.worker(plan,jobs[1].key),0)
            self.assertEqual(scheduler.worker(plan,jobs[2].key),0)
            self.assertEqual(output.read_text(),'ok')
            path,data=scheduler.load_plan(plan)
            failed=json.loads(common._state_path(path,jobs[0].key).read_text())
            self.assertEqual(failed['status'],'EXIT')
            with self.assertRaises(RuntimeError):
                scheduler.check_parent(path,data,{'key':'science'},jobs[0].key,failed)
            captured=[]
            with patch.object(common,'submit_graph',side_effect=lambda tasks,*a,**kw:captured.extend(tasks)):
                common._run_hpc(path,data,3,False,False,dependency_condition=scheduler.dependency_condition)
            self.assertEqual(captured[1]['depends'],[(jobs[0].key,'ended')])
            self.assertEqual(captured[2]['depends'],[(jobs[0].key,'ended'),(jobs[1].key,'ended')])
            captured.clear()
            with patch.object(common,'submit_graph',side_effect=lambda tasks,*a,**kw:captured.extend(tasks)):
                common._run_hpc(path,data,3,False,False)
            self.assertTrue(all(condition=='done' for task in captured for _,condition in task['depends']))

    def test_running_missing_and_killed_parent(self):
        with tempfile.TemporaryDirectory() as temp:
            plan=Path(temp)/'plan.json'; data={'mode':'hpc'}; job={'key':'dx.compare'}
            for status in ('RUN','PEND',None):
                with patch.object(common,'_probe',return_value={'42':status}), self.assertRaises(RuntimeError):
                    scheduler.check_parent(plan,data,job,'science',{'status':'RUN','job_id':'42'})
            with patch.object(common,'_probe',return_value={'42':'EXIT'}):
                scheduler.check_parent(plan,data,job,'science',{'status':'RUN','job_id':'42'})
            with self.assertRaises(RuntimeError):
                scheduler.check_parent(plan,data,job,'science',{})

    def test_report_excludes_stale_failed_results(self):
        import pandas as pd
        import dx_compare
        import timing_study
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); aod=root/'trait/aod/main/results'; aod.mkdir(parents=True)
            pd.DataFrame([{'dx':d,'endpoint_sha256':d} for d in ('good','bad')]).to_csv(aod/'dx_counts.tsv',sep='\t',index=False)
            for date,state in [('good','DONE'),('bad','FAILED')]:
                base=root/date; scope=base/'c1_correlate'
                for p in [base/'qc',base/'logs',scope/'qc/evaluation']: p.mkdir(parents=True,exist_ok=True)
                (base/'qc/dx.json').write_text(json.dumps({'dx':date,'trait':'trait','endpoint_sha256':date}))
                (base/'logs/analysis.status.json').write_text(json.dumps({'status':state,'started':'2026','steps':[]}))
                (scope/'qc/C1_FINAL_INDEX.json').write_text(json.dumps({'status':'completed','evaluation':True}))
                (scope/'qc/evaluation/COMPLETE.json').write_text('{}')
                pd.DataFrame([{'model':'glmnet','AUC':.8}]).to_csv(scope/'qc/evaluation/performance.tsv',sep='\t',index=False)
            with patch.object(dx_compare,'scope_for_read',side_effect=lambda layer,date:root/date), patch.object(dx_compare,'reliability_tables',return_value=({},{})):
                result=dx_compare.report(root,['trait'],['prot'],['good','bad'])
            self.assertEqual(result['performance'].dx.tolist(),['good'])
            self.assertEqual(len(result['availability']),2)
            with patch.object(timing_study,'scope_for_read',return_value=root/'bad'):
                for status in ('FAILED','EXIT','RUNNING','BLOCKED'):
                    (root/'bad/logs/analysis.status.json').write_text(json.dumps({'status':status}))
                    self.assertIsNone(timing_study.native_predictions(root,'bad','glmnet',{'bad':'bad'}))


if __name__=='__main__':
    unittest.main(verbosity=2)
