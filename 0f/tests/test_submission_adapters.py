"""Legacy receipts remain usable while project adapters share one dispatcher."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'0f'))
import lsf_submit
import command_jobs
from submission_io import receipt_lock


def module(name,project,file):
    sys.path.insert(0,str(ROOT/project/'f'))
    spec=importlib.util.spec_from_file_location(name,ROOT/project/'f'/file)
    mod=importlib.util.module_from_spec(spec);sys.modules[name]=mod;spec.loader.exec_module(mod)
    return mod


class AdapterTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.cmd=self.root/'child.cmd';self.cmd.write_text('#!/bin/bash\ntrue\n')
    def tearDown(self):self.tmp.cleanup()
    def test_allele_retry_keeps_unselected_parent(self):
        batch=module('test_allele_batch','allele3','batch_runner.py')
        state={'jobs':{'parent':{'status':'submitted','job_id':7}}}
        plan={'jobs':[{'id':'child','deps':['parent'],'throttle_after':None}]}
        commands={'child':{'path':str(self.cmd),'sha256':hashlib.sha256(self.cmd.read_bytes()).hexdigest()}}
        with patch.object(lsf_submit,'run_bsub',return_value=subprocess.CompletedProcess(['bsub'],0,'Job <8>','')) as client:
            batch.submit_hpc(plan,self.root,commands,state)
            batch.submit_hpc(plan,self.root,commands,state)
        self.assertEqual(client.call_count,1)
        self.assertEqual(client.call_args.args[0][0], 'bsub')
        self.assertEqual(client.call_args.args[0][-2:], ['-w','done(7)'])
        self.assertEqual(state['jobs']['child']['job_id'],8)
    def test_gu_keeps_journal_contract(self):
        gu=module('test_gu_hpc','gu','hpc.py');(self.root/'jobs').mkdir()
        (self.root/'jobs/prepare.cmd').write_text('#!/bin/bash\ntrue\n')
        journal=self.root/'events.jsonl'
        with journal.open('a') as out,patch.object(gu,'bsub_args',return_value=['bsub']),patch.object(lsf_submit,'run_bsub',return_value=subprocess.CompletedProcess(['bsub'],0,'Job <9>','')):
            job=gu.submit_one({'run_dir':str(self.root)},'prepare','test',None,None,0,out,'attempt')
        rows=[json.loads(line) for line in journal.read_text().splitlines()]
        self.assertEqual(job,'9');self.assertEqual([r['state'] for r in rows],['submitting','submitted'])
    def test_alphagenome_snapshots_shared_submission_engine(self):
        spec=importlib.util.spec_from_file_location('test_alpha',ROOT/'allele3/0.alphag_2vcf.py')
        alpha=importlib.util.module_from_spec(spec);spec.loader.exec_module(alpha)
        self.assertNotIn('allele3/f/batch_runner.py',alpha.EMBEDDED_RESOURCES)
        alpha.prepare_runtime(self.root)
        base=alpha.NATIVE_SOURCE.parents[2]
        for relative in ('allele3/f/batch_runner.py','0f/lsf_submit.py','0f/submission_io.py','0f/workflow_lock.py'):
            self.assertEqual((base/relative).read_bytes(),(ROOT/relative).read_bytes())
        result=subprocess.run([sys.executable,'-c',"import sys;sys.path.insert(0,'allele3/f');import batch_runner;assert batch_runner.submit_graph.__module__=='lsf_submit'"],cwd=base,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
    def test_receipt_lock_stops_concurrent_dispatch(self):
        with receipt_lock(self.cmd.with_suffix('.submission.json')),patch.object(lsf_submit,'run_bsub') as client:
            with self.assertRaises(RuntimeError):command_jobs.submit_cmd({},self.cmd)
        self.assertFalse(client.called)


if __name__=='__main__':unittest.main()
