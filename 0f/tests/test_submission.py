"""No real LSF calls: test fake clients, receipts, graph semantics and local jobs."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import lsf_submit as lsf
import submit_jobs
import command_jobs
from submission_io import atomic_json
from workflow_jobs import Job, write_plan, run_plan


class SubmissionTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.cmd=self.root/'work.cmd';self.cmd.write_text('#!/bin/bash\ntrue\n')
        self.rows=[]
    def tearDown(self):self.tmp.cleanup()
    def result(self,code=0,out='Job <101>',err=''):
        return subprocess.CompletedProcess(['bsub'],code,out,err)
    def attempt(self,**kw):
        return lsf.submit_attempt(['bsub'],path=self.cmd,persist=self.rows.append,**kw)
    def test_intent_precedes_submission(self):
        def fake(*a,**k):
            self.assertEqual(self.rows[-1]['status'],'submitting');return self.result()
        with patch.object(lsf,'run_bsub',side_effect=fake) as client:r=self.attempt()
        self.assertEqual(r['job_id'],'101');self.assertEqual(r['status'],'submitted');self.assertEqual(client.call_count,1)
    def test_rejection(self):
        with patch.object(lsf,'run_bsub',return_value=self.result(1,'','rejected')):r=self.attempt()
        self.assertEqual(r['status'],'failed')
    def test_nonzero_with_id_is_unknown(self):
        with patch.object(lsf,'run_bsub',return_value=self.result(1)):r=self.attempt()
        self.assertEqual((r['status'],r['job_id']),('unknown','101'))
    def test_missing_and_multiple_ids_are_unknown(self):
        for output in ('accepted','Job <1> Job <2>'):
            with patch.object(lsf,'run_bsub',return_value=self.result(out=output)):r=self.attempt()
            self.assertEqual(r['status'],'unknown')
    def test_timeout_keeps_partial_id(self):
        with patch.object(lsf,'run_bsub',side_effect=subprocess.TimeoutExpired(['bsub'],1,output=b'Job <99>')) as client:r=self.attempt()
        self.assertEqual((r['status'],r['job_id']),('unknown','99'));self.assertEqual(client.call_count,1)
    def test_hash_mismatch_never_submits(self):
        with patch.object(lsf,'run_bsub') as client:r=self.attempt(expected_sha256='wrong')
        self.assertFalse(client.called);self.assertEqual(r['status'],'failed')
    def test_failed_intent_write_never_submits(self):
        def broken(row):raise OSError('disk full')
        with patch.object(lsf,'run_bsub') as client:
            with self.assertRaises(OSError):lsf.submit_attempt(['bsub'],path=self.cmd,persist=broken)
        self.assertFalse(client.called)
    def test_failed_ack_write_leaves_intent(self):
        receipt=self.root/'receipt.json'
        def persist(row):
            if row['status']=='submitting':atomic_json(receipt,row)
            else:raise OSError('disk full')
        with patch.object(lsf,'run_bsub',return_value=self.result()) as client:
            with self.assertRaises(OSError):lsf.submit_attempt(['bsub'],path=self.cmd,persist=persist)
        self.assertEqual(client.call_count,1);self.assertEqual(json.loads(receipt.read_text())['status'],'submitting')
    def test_graph_edges_and_resume(self):
        tasks=[dict(key='a',path=self.cmd),dict(key='b',path=self.cmd,depends=[('a','done')]),dict(key='c',path=self.cmd,depends=[('b','ended')])]
        calls=[];saved={}
        def fake(argv,**kwargs):calls.append(argv);return self.result(out=f'Job <{len(calls)}>')
        with patch.object(lsf,'run_bsub',side_effect=fake):
            lsf.submit_graph(tasks,saved,persist=lambda k,r:saved.update({k:r}),receipt='test')
            lsf.submit_graph(tasks,saved,persist=lambda k,r:saved.update({k:r}),receipt='test')
        self.assertEqual(len(calls),3);self.assertEqual(calls[1][-1],'done(1)');self.assertEqual(calls[2][-1],'ended(2)')
    def test_graph_unknown_blocks_all(self):
        with patch.object(lsf,'run_bsub') as client:
            with self.assertRaises(RuntimeError):lsf.submit_graph([dict(key='a',path=self.cmd)],{'a':{'status':'unknown'}},persist=lambda *a:None,receipt='test')
        self.assertFalse(client.called)
    def test_batch_merge_conditions(self):
        for condition in ('done','ended'):
            calls=[]
            def fake(argv,**kw):calls.append(argv);return self.result(out=f'Job <{len(calls)}>')
            with patch.object(lsf,'run_bsub',side_effect=fake):
                ids=submit_jobs.run([self.cmd],'hpc',merge=self.cmd,merge_condition=condition)
            self.assertEqual(ids,['1','2']);self.assertEqual(calls[-1][-1],f'{condition}(1)')
    def test_legacy_grid_uncertain_blocks(self):
        state=self.root/'state.json';atomic_json(state,dict(schema='ld.submission.1',workers={},uncertain=True))
        with patch.object(lsf,'run_bsub') as client:
            with self.assertRaises(RuntimeError):submit_jobs.run([self.cmd],'hpc',state=state,state_schema='ld.submission.1')
        self.assertFalse(client.called)
    def test_capacity_wait_uses_one_client(self):
        fake=self.root/'client.py';counter=self.root/'count'
        fake.write_text(f"from pathlib import Path\nimport time,sys\np=Path({str(counter)!r});p.write_text(p.read_text()+'x' if p.exists() else 'x')\nprint('Pending job threshold reached. Retrying in 1 seconds...',file=sys.stderr,flush=True)\ntime.sleep(1.5)\nprint('Job <42>')\n")
        result=lsf.run_bsub([sys.executable,str(fake)],capture_output=True,text=True,timeout=.5)
        self.assertIn('42',result.stdout);self.assertEqual(counter.read_text(),'x')
    def test_real_timeout_kills_client_once(self):
        with self.assertRaises(subprocess.TimeoutExpired) as exc:
            lsf.run_bsub([sys.executable,'-c',"import time;print('Job <43>',flush=True);time.sleep(2)"],capture_output=True,text=True,timeout=.5)
        self.assertIn('43',lsf.decoded(exc.exception.stdout))
    def test_local_failure_stops_queued_work(self):
        self.cmd.write_text('#!/bin/bash\nexit 7\n');marker=self.root/'unexpected';later=self.root/'later.cmd'
        later.write_text(f'#!/bin/bash\ntouch {marker}\n')
        with self.assertRaises(RuntimeError):submit_jobs.run([self.cmd,later],'local',jobs=1)
        self.assertFalse(marker.exists())
    def test_local_array_uses_shared_pool(self):
        marker=self.root/'indices';self.cmd.write_text(f'#!/bin/bash\n#BSUB -J test[1-3]%1\n#BSUB -n 1\n#BSUB -R "rusage[mem=1MB/host]"\necho "$LSB_JOBINDEX" >> {marker}\n')
        result=command_jobs.local_submit({},self.cmd)
        self.assertEqual(result['status'],'submitted');self.assertEqual(marker.read_text().split(),['1','2','3'])
    def test_metag_receipt_blocks_unknown_retry(self):
        with patch.dict(os.environ,{'submit':'hpc'}),patch.object(lsf,'run_bsub',side_effect=subprocess.TimeoutExpired(['bsub'],1)) as client:
            with self.assertRaises(command_jobs.SubmissionError):command_jobs.submit_cmd({},self.cmd)
            with self.assertRaises(command_jobs.SubmissionError):command_jobs.submit_cmd({},self.cmd)
        self.assertEqual(client.call_count,1)
    def test_accepted_script_cannot_be_submitted_twice(self):
        with patch.dict(os.environ,{'submit':'hpc'}),patch.object(lsf,'run_bsub',return_value=self.result()) as client:
            command_jobs.submit_cmd({},self.cmd)
            with self.assertRaises(command_jobs.SubmissionError):command_jobs.submit_cmd({},self.cmd)
        self.assertEqual(client.call_count,1)
    def test_dag_engine_uses_shared_graph_and_resumes(self):
        jobs=[Job('a',['true'],resources={'cores':1,'memory_gib':1,'hours':1}),Job('b',['true'],depends=['a'],resources={'cores':1,'memory_gib':1,'hours':1})]
        plan=write_plan(jobs,self.root/'dag',mode='hpc',name='test');calls=[]
        def fake(argv,**kwargs):calls.append(argv);return self.result(out=f'Job <{len(calls)}>')
        with patch.object(lsf,'run_bsub',side_effect=fake):
            run_plan(plan,submit='hpc');run_plan(plan,submit='hpc')
        self.assertEqual(len(calls),2);self.assertEqual(calls[-1][-1],'done(1)')
    def test_invalid_timeout_rejected(self):
        for value in ('nan','inf','0','-1','bad'):
            with self.assertRaises(ValueError):lsf.submission_timeout(value)

    def test_local_jobs_obey_central_ceiling(self):
        active = peak = 0
        lock = threading.Lock()
        def fake(path, index=None):
            nonlocal active, peak
            with lock:
                active += 1
                peak = max(peak, active)
            time.sleep(.05)
            with lock:
                active -= 1
            return dict(cmd=str(path), index=index, returncode=0, status='DONE')
        with patch.object(submit_jobs, 'local_one', side_effect=fake):
            results = submit_jobs.run_local([self.cmd]*24, jobs=100)
        self.assertEqual(len(results), 24)
        self.assertEqual(peak, 8)

    def test_hpc_has_no_artificial_job_count_or_walltime_limit(self):
        jobs = [Job(f'job{i}', ['true'], resources={'cores':1, 'memory_gib':1, 'hours':1})
                for i in range(20)]
        plan = write_plan(jobs, self.root/'unlimited', mode='hpc', name='test')
        calls = []
        def fake(argv, **kwargs):
            calls.append(argv)
            return self.result(out=f'Job <{len(calls)}>')
        with patch.object(lsf, 'run_bsub', side_effect=fake):
            run_plan(plan, submit='hpc', jobs=1)
        self.assertEqual(len(calls), 20)
        self.assertTrue(all('-w' not in call for call in calls))
        for cmd in Path(plan).parent.joinpath('commands').glob('*.cmd'):
            self.assertNotIn('#BSUB -W', cmd.read_text())


if __name__=='__main__':unittest.main()
