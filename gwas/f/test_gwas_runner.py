"""Small synthetic jobs verify chromosome scheduling without running GWAS."""
import json
from pathlib import Path
import tempfile
import unittest
import gwas_runner as runner


class ChromosomeExecutionTest(unittest.TestCase):
    def make_job(self, root, name, inputs=(), fail=False):
        output=root/name
        code=("import pathlib,time,sys; "
              "start=time.time(); time.sleep(0.15); "
              "pathlib.Path(sys.argv[1]).write_text(str(start)+' '+str(time.time())); "
              f"sys.exit({1 if fail else 0})")
        return dict(name=name,cmd=['python3','-c',code,str(output)],
                    inputs=[str(root/x) for x in inputs],outputs=[str(output)],covars=None)

    def test_parallel_dependencies_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            jobs=[self.make_job(root,'null-model')]
            for chrom in ('chr1','chrX','chr2'):
                jobs += [self.make_job(root,chrom+'.view',['null-model']),
                         self.make_job(root,chrom+'.association',[chrom+'.view'])]
            jobs += [self.make_job(root,'merge-results',[c+'.association' for c in ('chr1','chrX','chr2')])]
            plan=root/'plan.json'
            plan.write_text(json.dumps(dict(config={'jobs':2},replace=False,jobs=jobs)))
            runner.execute(plan)
            times={j['name']:list(map(float,(root/j['name']).read_text().split())) for j in jobs}
            self.assertLess(times['chr1.view'][0],times['chrX.view'][1])
            self.assertLess(times['chrX.view'][0],times['chr1.view'][1])
            for chrom in ('chr1','chrX','chr2'):
                self.assertGreaterEqual(times[chrom+'.view'][0],times['null-model'][1])
                self.assertGreaterEqual(times[chrom+'.association'][0],times[chrom+'.view'][1])
                self.assertGreaterEqual(times['merge-results'][0],times[chrom+'.association'][1])
            events=sorted([(start,1) for start,end in times.values()]+[(end,-1) for start,end in times.values()])
            active=0
            for _,delta in events:
                active+=delta
                self.assertLessEqual(active,2)
            before={j['name']:(root/j['name']).stat().st_mtime_ns for j in jobs}
            runner.execute(plan)
            self.assertEqual(before,{j['name']:(root/j['name']).stat().st_mtime_ns for j in jobs})

    def test_failure_blocks_merge(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            jobs=[self.make_job(root,'chr1.association',fail=True),
                  self.make_job(root,'chrX.association'),self.make_job(root,'merge-results')]
            plan=root/'plan.json'
            plan.write_text(json.dumps(dict(config={'jobs':2},replace=False,jobs=jobs)))
            with self.assertRaises(RuntimeError):
                runner.execute(plan)
            self.assertFalse((root/'merge-results').exists())
            self.assertFalse((root/'chr1.association.done.json').exists())


if __name__=='__main__':
    unittest.main()
