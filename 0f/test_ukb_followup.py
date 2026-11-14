"""Canonical date propagation and frozen-config contracts; no cohort data."""
import copy
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]


class FollowupTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / 'scripts'
        for project in ('0f', 'le8', 'allele3', 'flux', 'gu', 'ukb'):
            shutil.copytree(ROOT / project, self.root / project,
                            ignore=shutil.ignore_patterns('__pycache__', '.pytest_cache'))
        self.source = self.root / '0f/phenotype.R'
        self.original = self.source.read_text()
        self.expression('as.Date("2022-12-31")')
        spec = importlib.util.spec_from_file_location('fixture_followup', self.root / '0f/ukb_followup.py')
        self.adapter = importlib.util.module_from_spec(spec); spec.loader.exec_module(self.adapter)

    def expression(self, expression):
        lines = self.original.splitlines(True)
        lines = ['date_follow_end <- '+expression+'\n' if x.startswith('date_follow_end <-') else x for x in lines]
        self.source.write_text(''.join(lines))

    def run_code(self, code):
        result = subprocess.run([sys.executable, '-c', code], cwd=self.root,
                                capture_output=True, text=True, timeout=90)
        self.assertEqual(result.returncode, 0, result.stdout+'\n'+result.stderr)

    def test_expression_cache_and_failures(self):
        real = subprocess.run
        with patch.object(self.adapter.subprocess, 'run', wraps=real) as spy:
            first = self.adapter.resolve_followup(environ={})
            self.assertEqual(first['date'], '2022-12-31')
            self.assertEqual(first['source'], str(self.source))
            first['date'] = 'mutated'
            self.assertEqual(self.adapter.resolve_followup(environ={})['date'], '2022-12-31')
            self.assertEqual(spy.call_count, 1)
            self.expression('as.Date("2022-12-30") + 1L')
            self.assertEqual(self.adapter.resolve_followup(environ={})['date'], '2022-12-31')
            self.assertEqual(spy.call_count, 2)
            self.expression('as.Date("2022-12-29")')
            self.assertEqual(self.adapter.resolve_followup(environ={})['date'], '2022-12-29')
        for expression in ('as.Date(NA)', '"2022-12-31"', 'as.Date(character())'):
            self.expression(expression)
            with self.assertRaises(RuntimeError): self.adapter.resolve_followup(environ={})
        self.source.write_text(self.original.replace('date_follow_end <-', 'unused_date <-', 1))
        with self.assertRaises(RuntimeError): self.adapter.resolve_followup(environ={})
        with self.assertRaises(RuntimeError): self.adapter.resolve_followup(rscript='/nonexistent/Rscript', environ={})
        with self.assertRaisesRegex(ValueError, 'no longer supported'):
            self.adapter.resolve_followup(environ={'DATE_FOLLOW_END': '2020-01-01'})

    def test_python_consumers_and_saved_configs(self):
        self.run_code(r'''
import sys, importlib.util, json, tempfile, copy
from pathlib import Path
from types import SimpleNamespace
sys.path[:0]=[str(Path('le8/f').resolve()),str(Path('0f').resolve())]
def load(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
le8=load('le8_common','le8/f/0.common.py')
host=load('le8_host','le8/f/0.host.py')
try:host.cohort_parser().parse_args(['--end-date','2020-01-01'])
except SystemExit as exc:assert exc.code==2
else:raise AssertionError('Host parser accepted --end-date')

values,sources,env=le8.shared_analysis_settings(SimpleNamespace(),environ={})
assert values['end_date']=='2022-12-31' and 'DATE_FOLLOW_END' not in env
assert sources['end_date']['definition']['variable']=='date_follow_end'
sys.path.insert(0,str(Path('allele3/f').resolve()))
planner=load('planner','allele3/f/burden_plan.py')
args=planner.parser().parse_args(['--submit','local','--memory-limit-gb','8'])
cfg=planner.resolve(args)[0]
assert cfg['phenotypes']['end_date']=='2022-12-31'
worker=load('burden','allele3/f/burden.py')
import pandas as pd
people=pd.DataFrame({'date_attend':['2020-01-01']*4,'fod_icd10_cvd_cad':['2023-01-01']*4,
 'date_lost':[None,'2021-01-01',None,None],'date_death':[None,None,'2021-06-01',None]})
derived,audit=worker.phenotype_derive_outcome(people,'cvd_cad',worker.load_config())
assert derived['censor_date'].dt.strftime('%Y-%m-%d').tolist()==['2022-12-31','2021-01-01','2021-06-01','2022-12-31']
assert derived['event'].sum()==0
assert audit['followup_provenance']['date']=='2022-12-31'

p=Path('saved.json');p.write_text(json.dumps(cfg))
assert worker.load_config(p)['phenotypes']['end_date']=='2022-12-31'
args.config=str(p.resolve());assert planner.resolve(args)[0]['phenotypes']['end_date']=='2022-12-31'
for bad in [{'phenotypes':{'end_date':'2020-01-01'}},{'phenotypes':{'outcomes':{'not_selected':{'end_date':'2020-01-01'}}}}]:
 p.write_text(json.dumps(bad))
 for call in [lambda:worker.load_config(p),lambda:planner.resolve(args)]:
  try:call()
  except ValueError:pass
  else:raise AssertionError('Accepted override')
p.write_text(json.dumps(cfg))
sys.path.insert(0,str(Path('flux/f').resolve()))
common=load('common','flux/f/common.py');data=load('data','flux/f/data.py')
f=common.load_config();assert f['data']['admin_end']=='2022-12-31'
assert data.Importer(None,f,Path('.')).admin_end==data.day('2022-12-31')
y=Path('saved.yml');common.save_config(f,y);assert common.load_config(y)['data']['admin_end']=='2022-12-31'
for call in [lambda:common.load_config(overrides={'data':{'admin_end':'2020-01-01'}}),lambda:data.Importer(None,{'data':{'admin_end':'2020-01-01'}},Path('.'))]:
 try:call()
 except ValueError:pass
 else:raise AssertionError('Accepted FLUX override')
y.write_text('data:\n  admin_end: 2020-01-01\n')
try:common.load_config(y)
except ValueError:pass
else:raise AssertionError('Accepted YAML override')
common.save_config(f,y)
canonical=Path('0f/phenotype.R');canonical.write_text(canonical.read_text().replace('as.Date("2022-12-31")','as.Date("2022-12-30") + 1L',1))
assert common.load_config(y)['data']['admin_end']=='2022-12-31'
assert worker.load_config(p)['phenotypes']['end_date']=='2022-12-31'
canonical.write_text(canonical.read_text().replace('as.Date("2022-12-30") + 1L','as.Date("2021-12-31")',1))
for call in [lambda:worker.load_config(p),lambda:planner.resolve(args),lambda:common.load_config(y)]:
 try:call()
 except ValueError as exc:assert 'Stale' in str(exc),str(exc)
 else:raise AssertionError('Accepted stale plan')
''')

    def test_native_r_and_rng(self):
        # Execute the actual canonical-reader expression and getter without loading
        # unrelated model packages or accessing phenotype input paths.
        code = r'''
args <- commandArgs(TRUE)
e <- new.env(parent=globalenv())
e$.le8_canonical_phenotype_file <- normalizePath("0f/phenotype.R")
expressions <- parse("le8/f/0.common.R")
for (x in expressions) {
 if (is.call(x) && identical(x[[1]],as.name("<-")) &&
     is.symbol(x[[2]]) && as.character(x[[2]]) %in% c("le8_follow_end",".le8_canonical_phe")) {
   set.seed(918); before <- .Random.seed
   eval(x,envir=e)
   stopifnot(identical(before,.Random.seed))
 }
}
stopifnot(as.character(e$le8_follow_end())=="2022-12-31")
Sys.setenv(DATE_FOLLOW_END="2020-01-01")
stopifnot(inherits(try(e$le8_follow_end(),silent=TRUE),"try-error"))
Sys.unsetenv("DATE_FOLLOW_END")
# The GU extra source preserves RNG; its original subsequent set.seed is unchanged.
g <- new.env(parent=globalenv());g$.gu_canonical_phenotype_file <- e$.le8_canonical_phenotype_file
for(x in parse("gu/f/phenotype.R")) if(is.call(x) && identical(x[[1]],as.name("<-")) && identical(x[[2]],as.name("date_follow_end"))) {
 set.seed(77);before <- .Random.seed;eval(x,envir=g);stopifnot(identical(before,.Random.seed),as.character(g$date_follow_end)=="2022-12-31")
}
'''
        result=subprocess.run(['Rscript','--vanilla','-e',code],cwd=self.root,capture_output=True,text=True,timeout=60)
        self.assertEqual(result.returncode,0,result.stderr)

    def test_help_without_r_and_removed_cli(self):
        env=dict(os.environ,R_BIN='/nonexistent/Rscript')
        for entry in ['le8/f/0.common.py','le8/f/c1_correlate.py','flux/f/flux.py']:
            result=subprocess.run([sys.executable,entry,'--help'],cwd=self.root,env=env,capture_output=True,text=True,timeout=180)
            self.assertEqual(result.returncode,0,(entry,result.stderr))
            result=subprocess.run([sys.executable,entry,'--end-date','2020-01-01'],cwd=self.root,env=env,capture_output=True,text=True,timeout=180)
            self.assertNotEqual(result.returncode,0,entry)
            self.assertIn('--end-date',result.stderr)


if __name__ == '__main__': unittest.main()
