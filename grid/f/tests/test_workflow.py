#!/usr/bin/env python3
"""Integration checks with synthetic PCs; no UKB data or real PLINK work."""
import csv
import gzip
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


class WorkflowTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='grid-workflow-')
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.bin = self.work / 'bin'
        self.bin.mkdir()
        self.env = dict(os.environ, GRID_CONDA_ENV='/nonexistent-test-env',
                        PATH=f'{self.bin}:{os.environ["PATH"]}')
        self.pca = self.work / 'pca' / 'ukb.discodivas.pca.tsv.gz'
        self.pca.parent.mkdir()
        self.phe = self.work / 'phe.rds'
        self.out = self.work / 'out'
        self.args = ['--pca-file', str(self.pca), '--phe-file', str(self.phe),
                     '--output-root', str(self.out), '--jobs', '1']

    def run_cmd(self, *args, ok=True):
        result = subprocess.run(list(map(str, args)), env=self.env,
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if ok:
            self.assertEqual(result.returncode, 0, result.stdout)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout)
        return result.stdout

    def fixture(self):
        with open('/mnt/d/files/DiscoDivas/med.g1000.4pop.tsv') as f:
            centers = list(csv.DictReader(f, delimiter='\t'))
        rows = []
        for i in range(80):
            center = centers[i % 4]
            values = [float(center[f'PC{j}']) + (i % 7 - 3) * (j + 1)
                      if j <= 10 else (i % 11 - 5) * j for j in range(1, 21)]
            rows.append([str(i + 1), *values])
        with gzip.open(self.pca, 'wt') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['IID', *[f'PC{j}' for j in range(1, 21)]])
            writer.writerows(rows)
        self.run_cmd('Rscript', '-e',
                     'saveRDS(data.frame(eid=as.character(1:80), '
                     'ethnicity=rep(c(4,5,1,3),20)), commandArgs(TRUE)[1])', self.phe)
        self.stub_plink(fail=True)

    def stub_plink(self, fail):
        script = self.bin / 'plink2'
        if fail:
            script.write_text('#!/bin/sh\necho "UNEXPECTED PLINK" >&2\nexit 99\n')
        else:
            script.write_text('''#!/usr/bin/env python3
import sys
from pathlib import Path
a=sys.argv
out=Path(a[a.index('--out')+1])
n=int(a[a.index('--score-col-nums')+1].split('-')[1])-6
Path(str(out)+'.sscore').write_text('#IID\\t'+'\\t'.join(f'PC{i}_SUM' for i in range(1,n+1))+'\\n1\\t'+'\\t'.join(['1']*n)+'\\n')
Path(str(out)+'.sscore.vars').write_text('rs1\\n')
''')
        script.chmod(0o755)

    def test_existing_projection_resume_and_disco(self):
        self.fixture()
        original = self.pca.read_bytes()
        first = self.run_cmd('bash', ROOT / 'grid.sh', 'pca', *self.args)
        self.assertIn('SKIP projection', first)
        self.assertEqual(self.pca.read_bytes(), original)
        distance = self.pca.parent / 'ukb_reference_distances.tsv.gz'
        anc = self.pca.parent / 'ukb.ancestry.auto.tsv.gz'
        self.assertTrue(distance.stat().st_size > 0)
        with gzip.open(anc, 'rt') as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        self.assertEqual(len(rows), 80)
        self.assertIn('ancestry_probability', rows[0])
        stamp = distance.stat().st_mtime_ns
        second = self.run_cmd('bash', ROOT / 'pca.sh', *self.args)
        self.assertIn('SKIP reference distances/QC', second)
        self.assertIn('SKIP ancestry', second)
        anc.unlink()
        third = self.run_cmd('bash', ROOT / 'pca.sh', *self.args)
        self.assertIn('SKIP reference distances/QC', third)
        self.assertEqual(distance.stat().st_mtime_ns, stamp)
        self.assertTrue(anc.exists())
        scores = self.out / 'height' / 'scores'
        scores.mkdir(parents=True)
        with gzip.open(scores / 'csx.tsv.gz', 'wt') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['eid', 'CSX_AFR', 'CSX_EAS', 'CSX_EUR', 'CSX_SAS'])
            writer.writerows([[str(i), i, i * 2, i * 3, i * 4] for i in range(1, 81)])
        self.run_cmd('bash', ROOT / 'disco.sh', *self.args)
        self.assertTrue((scores / 'disco_zero.tsv.gz').stat().st_size > 0)
        log = (self.out / 'height' / 'log' / 'disco.log').read_text()
        self.assertNotIn('Error', log)

    def test_projection_persists_when_distance_fails(self):
        # Stub only the expensive genotype scoring; use real R combination/QC.
        self.fixture()
        self.pca.unlink()
        self.stub_plink(fail=False)
        imp = self.work / 'imp'
        imp.mkdir()
        for chrom in range(1, 23):
            for ext in ('pgen', 'psam', 'pvar.zst' if chrom % 2 else 'pvar'):
                (imp / f'chr{chrom}.{ext}').write_text('fixture\n')
        failure = self.run_cmd('bash', ROOT / 'pca.sh', *self.args,
                              '--dir-imp', imp, '--med-file', self.work / 'missing.tsv', ok=False)
        self.assertIn('Missing/empty file', failure)
        self.assertTrue(self.pca.stat().st_size > 0)
        self.run_cmd('python3', ROOT / 'f/pca_cache.py', self.pca, '20')
        self.stub_plink(fail=True)
        # A real completed projection must be reusable even with no genotypes.
        resumed = self.run_cmd('bash', ROOT / 'pca.sh', *self.args,
                              '--dir-imp', self.work / 'absent')
        self.assertIn('SKIP projection', resumed)

    def test_entrypoints_and_grid_steps(self):
        for script in ('pca.sh', 'csx.sh', 'disco.sh', 'grid.sh'):
            self.run_cmd('bash', ROOT / script, '--help')
        for module in ('ancestry', 'csx', 'disco', 'all', 'arg', 'ld', 'transport', 'fit', 'weights', 'score'):
            self.run_cmd('bash', ROOT / 'grid.sh', module, ok=False)
        for module in ('grid', 'eval'):
            self.run_cmd('bash', ROOT / 'grid.sh', module, '--help')
        # GRID validates prerequisites internally before any expensive computation.
        result = self.run_cmd('bash', ROOT / 'grid.sh', 'grid', '--output-root', self.out,
                              '--arg-dir', self.work / 'absent', '--chrs', '22', ok=False)
        self.assertIn('missing GRID ARG-Needle data', result)

    def test_grid_runs_internal_steps_in_order_and_stops_on_failure(self):
        app = self.work / 'app'
        (app / 'f').mkdir(parents=True)
        for name in ('grid.sh', 'f/grid.sh', 'f/common.sh', 'f/grid_common.sh', 'f/environment.sh'):
            shutil.copy2(ROOT / name, app / name)
        order = self.work / 'steps.txt'
        self.env['TEST_STEP_LOG'] = str(order)
        steps = ['arg', 'ld', 'transport', 'fit', 'weights', 'score']
        for step in steps:
            (app / 'f' / f'{step}.sh').write_text(
                f'#!/bin/bash\necho {step} >> "$TEST_STEP_LOG"\n')
        self.run_cmd('bash', app / 'grid.sh', 'grid', '--output-root', self.out)
        self.assertEqual(order.read_text().splitlines(), steps)
        order.unlink()
        (app / 'f' / 'fit.sh').write_text('#!/bin/bash\necho fit >> "$TEST_STEP_LOG"\nexit 7\n')
        self.run_cmd('bash', app / 'grid.sh', 'grid', '--output-root', self.out, ok=False)
        self.assertEqual(order.read_text().splitlines(), steps[:4])


if __name__ == '__main__':
    unittest.main(verbosity=2)
