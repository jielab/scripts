import sys,types,json,tempfile,subprocess,unittest,os
from pathlib import Path
from unittest.mock import patch
ROOT=Path(os.environ.get('GRID_TEST_ROOT', Path(__file__).resolve().parents[1]))
sys.path[:0]=[str(ROOT/'f'),str(ROOT.parent/'0f')]
import trait_workflow
from abm import count
source=Path(sys.argv.pop(1)) if len(sys.argv)>1 else ROOT/'f/hpc_workflow.py'
planner=types.ModuleType('bundle_planner');planner.__file__=str(ROOT/'f/hpc_workflow.py');sys.modules[planner.__name__]=planner
exec(compile(source.read_text(),planner.__file__,'exec'),planner.__dict__)
class Bundles(unittest.TestCase):
 def plan(self,mode,module='all'):
  args=planner.parser().parse_args(['--entry','grid','--workflow-module',module,'--traits','height,ldl,t2dm','--queue','auto','--job-cores','40','--replace','FALSE'])
  planner.normalize_runs(args)
  with patch.object(planner,'resolve_mode',return_value=mode):return planner.build_jobs(args,[])
 def test_both_modes_same_bundles(self):
  for mode in ('local','hpc'):
   jobs=self.plan(mode);by={j.key:j for j in jobs};self.assertEqual(len(jobs),11)
   for trait in ('height','ldl','t2dm'):
    j=by[trait+'.all'];steps=json.loads(j.env['GRID_TRAIT_STEPS'])
    self.assertEqual([s['argv'][2] for s in steps],['radius_twins','evaluate'])
    self.assertEqual(j.depends,[trait+'.abm.cpu_geometry'])
    self.assertEqual(j.env['GRID_REPLACE'],'FALSE')
    self.assertEqual(j.resources['gpus'],0)
    self.assertGreaterEqual(j.resources['memory_gib'],256 if mode == 'hpc' else 1)
    self.assertEqual(planner.scheduler_name('grid',j.key),'grid.'+trait)
   j=by['abm.cpu_final'];self.assertEqual(j.depends,['height.all','ldl.all','t2dm.all'])
   self.assertEqual([s['argv'][2] for s in json.loads(j.env['GRID_TRAIT_STEPS'])],['architecture','final'])
   self.assertEqual(json.loads(j.env['GRID_TRAIT_STEPS'])[0]['argv'][4],'height,ldl,t2dm')
   if os.environ.get('GRID_TEST_CHECK_INPUTS') == '1':
    planner.check_upstream_files(jobs)
 def test_standalone_modules_unchanged(self):
  for mode in ('local','hpc'):
   for module in ('prepare','twins','evaluate','architecture','final'):
    self.assertTrue(all('GRID_TRAIT_STEPS' not in j.env for j in self.plan(mode,module)))
 def test_native_checks_and_reuse(self):
  with tempfile.TemporaryDirectory() as td:
   d=Path(td);a=d/'result';a.write_text('valid')
   count.commit_stage(d,'twins',[a]);self.assertEqual(count.verify_stage(d,'twins')['status'],'complete')
   self.assertTrue(count.reusable_stage(types.SimpleNamespace(replace=False),d,'twins'))
   self.assertFalse(count.reusable_stage(types.SimpleNamespace(replace=True),d,'twins'))
   a.unlink();self.assertFalse(count.reusable_stage(types.SimpleNamespace(replace=False),d,'twins'))
   with self.assertRaises(ValueError):count.verify_stage(d,'twins')
 def test_local_execution_order_and_failure(self):
  with tempfile.TemporaryDirectory() as td:
   receipt=Path(td)/'ready';output=Path(td)/'out'
   def step(key,code):return dict(key=key,argv=[sys.executable,'-c',code],env={})
   steps=[step('first',f'from pathlib import Path;Path({str(receipt)!r}).write_text("done")'),step('second',f'from pathlib import Path;assert Path({str(receipt)!r}).read_text()=="done";Path({str(output)!r}).touch()')]
   self.assertEqual(trait_workflow.run_steps(steps,dict(os.environ)),0);self.assertTrue(output.exists())
   output.unlink();steps[0]=step('first','raise SystemExit(7)')
   self.assertEqual(trait_workflow.run_steps(steps,dict(os.environ)),7);self.assertFalse(output.exists())
if __name__=='__main__':unittest.main()
