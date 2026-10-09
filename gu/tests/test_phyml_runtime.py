"""Regression coverage for build caching, durable resume and timeout deferral."""
import csv
import importlib.util
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import unittest
from concurrent.futures import ProcessPoolExecutor
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('phyml_runtime_test', ROOT / 'f/phyml.py')
phyml = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = phyml
spec.loader.exec_module(phyml)
common = sys.modules['gu_0_common']


def table(path, records, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(records[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(records)


def add_timeout(args):
    out, index = args
    out = Path(out)
    (out / 'loci').mkdir(parents=True)
    Path(str(out / 'loci/haplotypes.phy') + '.phyml.timeout.json').write_text(json.dumps(dict(
        rc=124, timeout_seconds=7200, elapsed_seconds=7200, bootstrap_completed=index,
        n_sequences=100, n_sites=200, recorded_at='test', reason='time_limit')))
    phyml.update_timeout_list(out, 'test', f'1:{index}:A:G', True)


class RuntimeTests(unittest.TestCase):
    def setUp(self):
        # The result layer intentionally rejects /tmp as a permanent work root.
        self.temp = tempfile.TemporaryDirectory(prefix='.runtime-test-', dir=ROOT)
        self.root = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)
        self.addCleanup(lambda: shutil.rmtree(common.gu_native_view_root(self.root), ignore_errors=True))

    def fixture(self):
        out = self.root / 'phyml/test/chr1_test'
        (out / 'loci').mkdir(parents=True)
        refs = self.root / 'refs'
        refs.mkdir()
        (refs / 'chr1.vcf.gz').write_text('reference')
        target = self.root / 'target'
        target.mkdir()
        (target / 'chr1.pvar').write_text('genotypes')
        (target / 'samples.txt').write_text('sample panel')
        lead = dict(locus_id='1:2:A:G', chr='1', lead_pos='2')
        leads = self.root / 'leads.tsv'
        table(leads, [lead])
        bed = out / (out.name + '.bed')
        bed.write_text('1\t1\t3\t1:2:A:G\n')
        cmd = out / (out.name + '.cmd')
        argv = [str(ROOT / 'gu.sh'), 'phyml', '--loci', str(bed), '--target-dir', str(target / 'chr'),
                '--target', 'test', '--replace-phyml', 'FALSE', '--memory-cap', '32G']
        cmd.write_text(f'export GU_PHYML_LEAD_TABLE={shlex.quote(str(leads))}\nexec {shlex.join(argv)}\n')
        final = out / 'final'
        table(final / 'gwas_lead.tsv', [lead])
        summary = [dict(status='insufficient_high_LD_markers', reason='too few', lineage=x) for x in ('Neanderthal', 'Denisovan')]
        table(final / 'gwas_loci.tsv', summary)
        table(final / 'trees.tsv', [dict(tree_status='not_run')])
        for name in ['gwas_haplotypes', 'gwas_copies', 'loci', 'evidence_trees', 'haplotypes', 'haplotype_samples', 'skipped_loci']:
            table(final / (name + '.tsv'), [], ['locus_id'])
        (final / 'gwas_parameters.json').write_text(json.dumps(dict(workflow=phyml.WORKFLOW)))
        cmd.with_suffix('.log').write_text(f'archaic reference root={refs}\nsample metadata={target / "samples.txt"}\n'
                                         f'analysis_unit={out.name} status=complete\n')
        return out, cmd, refs

    def test_build_cache_concurrency_and_invalidation(self):
        source = self.root / 'chr1.pvar'
        gold = self.root / 'gold.tsv'
        source.write_text('37\n')
        gold.write_text('sentinels\n')
        env = dict(os.environ, CHECK_GRCH_SNP_LIST=str(gold), GU_BUILD_CHECK_CACHE_DIR=str(self.root / 'cache'),
                   TEST_CALLS=str(self.root / 'calls'), TEST_SOURCE=str(source))
        script = f'''set -euo pipefail
source {shlex.quote(str(ROOT / 'f/0.build-check.sh'))}
check_GRCH() {{ echo checked >> "$TEST_CALLS"; sleep 0.15; [[ $(cat "$1") == "$2" ]]; }}
phe_zcat() {{ cat "$@"; }}
phe_check_grch_tabix_rows() {{ return 1; }}
gu_check_grch_pvar_positions() {{ return 1; }}
gu_check_build_cached "$TEST_SOURCE" "${{TEST_BUILD:-37}}" pfile
'''
        workers = [subprocess.Popen(['bash', '-c', script], env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE) for _ in range(6)]
        outputs = [p.communicate() for p in workers]
        self.assertTrue(all(p.returncode == 0 for p in workers), outputs)
        self.assertEqual((self.root / 'calls').read_text().count('checked'), 1)
        source.write_text('37\n')  # Same size/content, new mtime: invalidate.
        subprocess.run(['bash', '-c', script], env=env, check=True, capture_output=True)
        gold.write_text('changed sentinels\n')
        subprocess.run(['bash', '-c', script], env=env, check=True, capture_output=True)
        for _ in range(2):
            result = subprocess.run(['bash', '-c', script], env=dict(env, TEST_BUILD='38'), capture_output=True)
            self.assertNotEqual(result.returncode, 0)
        self.assertEqual((self.root / 'calls').read_text().count('checked'), 5)

    def test_old_missing_lead_is_rechecked_once_after_allele_fix(self):
        out, cmd, refs = self.fixture()
        summary = [dict(status='lead_absent_or_ambiguous', reason='old ordered allele match', lineage=x)
                   for x in ('Neanderthal', 'Denisovan')]
        table(out / 'final/gwas_loci.tsv', summary)
        self.assertTrue(phyml.process('seal', cmd, refs))
        self.assertFalse(phyml.process('check', cmd, refs))
        params = out / 'final/gwas_parameters.json'
        params.write_text(json.dumps(dict(workflow=phyml.WORKFLOW, lead_matching_rule=phyml.LEAD_MATCHING_RULE)))
        self.assertTrue(phyml.process('seal', cmd, refs))
        self.assertTrue(phyml.process('check', cmd, refs))

    def test_receipt_survives_compact_and_restore(self):
        out, cmd, refs = self.fixture()
        self.assertTrue(phyml.process('seal', cmd, refs))
        common.gu_compact_results(self.root, ['phyml'], out)
        self.assertFalse((out / 'final').exists())
        self.assertTrue(cmd.exists())
        self.assertTrue(cmd.with_suffix('.log').exists())
        self.assertTrue((out / '.phyml.locus.complete.json').exists())
        with tarfile.open(out / 'phyml.raw.tar.gz') as bundle:
            self.assertIn(cmd.name, bundle.getnames())
            self.assertIn(cmd.with_suffix('.log').name, bundle.getnames())
            self.assertIn('.phyml.locus.complete.json', bundle.getnames())
        common.gu_restore_results(self.root, 'phyml', out)
        self.assertTrue(phyml.process('check', cmd, refs))
        (self.root / 'target/chr1.pvar').write_text('different genotypes')
        self.assertFalse(phyml.process('check', cmd, refs))

    def test_actual_timeout_kills_children_and_defers_resume(self):
        out, cmd, refs = self.fixture()
        phy = out / 'loci/haplotypes.phy'
        phy.write_text('4 4\nA AAAA\nB AAAT\nC AATT\nD TTTT\n')
        binary = self.root / 'bin/phyml'
        binary.parent.mkdir()
        binary.write_text(f'#!{sys.executable}\nimport sys,time,subprocess\nfrom pathlib import Path\n'
                          'p=sys.argv[sys.argv.index("-i")+1]\n'
                          'Path(p+"_phyml_tree.txt").write_text("(A,B,C,D);\\n")\n'
                          'child=subprocess.Popen([sys.executable,"-c","import time;time.sleep(30)"])\n'
                          'Path(p+".child").write_text(str(child.pid))\ntime.sleep(30)\n')
        binary.chmod(0o755)
        result = subprocess.run([sys.executable, str(ROOT / 'f/phyml.py'), 'run', '--phy', str(phy),
                                 '--scope', 'test', '--cpus', '1', '--timeout', '0.5'],
                                env=dict(os.environ, PATH=str(binary.parent) + ':' + os.environ['PATH']), capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 124, result.stderr.decode())
        self.assertTrue(phy.exists())
        self.assertFalse(Path(str(phy) + '_phyml_tree.txt').exists())
        status = phyml.rows(Path(str(phy) + '.phyml.run.status.tsv'))[0]
        self.assertEqual(status['status'], 'TIMEOUT')
        recorded = Path(str(phy) + '.phyml.timeout.json').read_bytes()
        again = subprocess.run([sys.executable, str(ROOT / 'f/phyml.py'), 'run', '--phy', str(phy),
                                '--scope', 'test', '--cpus', '1', '--timeout', '7200'],
                               env=dict(os.environ, PATH=str(binary.parent) + ':' + os.environ['PATH'], PHYML_RETRY_TIMEOUTS='FALSE'),
                               capture_output=True, timeout=5)
        self.assertEqual(again.returncode, 124)
        self.assertEqual(Path(str(phy) + '.phyml.timeout.json').read_bytes(), recorded)
        child = Path('/proc') / Path(str(phy) + '.child').read_text() / 'stat'
        self.assertTrue(not child.exists() or child.read_text().split()[2] == 'Z')
        summary = [dict(status='tree_failed', reason='tree_timeout', lineage=x) for x in ('Neanderthal', 'Denisovan')]
        table(out / 'final/gwas_loci.tsv', summary)
        table(out / 'final/trees.tsv', [dict(tree_status='failed')])
        phyml.update_timeout_list(out, 'test', '1:2:A:G', True)
        self.assertTrue(phyml.process('seal', cmd, refs))
        with patch.dict(os.environ, {'PHYML_RETRY_TIMEOUTS': 'FALSE'}):
            self.assertTrue(phyml.process('check', cmd, refs))
        with patch.dict(os.environ, {'PHYML_RETRY_TIMEOUTS': 'TRUE'}):
            self.assertFalse(phyml.process('check', cmd, refs))
        common.gu_compact_results(self.root, ['phyml'], out)
        common.gu_restore_results(self.root, 'phyml', out)
        self.assertTrue(phyml.process('check', cmd, refs))
        summary[0]['reason'] = 'numerical_failure'
        table(out / 'final/gwas_loci.tsv', summary)
        with self.assertRaises(ValueError):
            phyml.process('seal', cmd, refs)

    def test_timeout_list_parallel_writers_and_recovery(self):
        dataset = self.root / 'phyml/test'
        with ProcessPoolExecutor(max_workers=4) as workers:
            list(workers.map(add_timeout, [(str(dataset / f'unit{i}'), i) for i in range(8)]))
        self.assertEqual(len(phyml.rows(dataset / 'phyml.timeouts.tsv')), 8)
        phyml.update_timeout_list(dataset / 'unit2', 'test', '1:2:A:G')
        self.assertEqual(len(phyml.rows(dataset / 'phyml.timeouts.tsv')), 7)
        self.assertFalse((dataset / 'unit2/phyml.timeout.json').exists())

    def test_numerical_failure_reports_cause_without_accepting_partial_tree(self):
        phy = self.root / 'haplotypes.phy'
        phy.write_text('4 3\nA GGG\nB GGG\nC GGG\nD GGG\n')
        binary = self.root / 'bin/phyml'
        binary.parent.mkdir()
        binary.write_text(f'#!{sys.executable}\nimport sys\nfrom pathlib import Path\n'
                          'p=sys.argv[sys.argv.index("-i")+1]\n'
                          'Path(p+"_phyml_tree.txt").write_text("(A,B,C,D);\\n")\n'
                          'print(". Failed to invert the matrix.")\n'
                          'print(". Cannot work out eigen vectors.")\nsys.exit(1)\n')
        binary.chmod(0o755)
        result = subprocess.run([sys.executable, str(ROOT / 'f/phyml.py'), 'run', '--phy', str(phy),
                                 '--scope', 'test', '--cpus', '1', '--timeout', '10'],
                                env=dict(os.environ, PATH=str(binary.parent) + ':' + os.environ['PATH']),
                                capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 1)
        self.assertIn('numerical_model_fit_failed', result.stderr.decode())
        self.assertNotIn('interrupted run', result.stderr.decode())
        self.assertEqual(phyml.rows(Path(str(phy) + '.phyml.run.status.tsv'))[0]['status'], 'FAILED')
        self.assertFalse(Path(str(phy) + '.phyml.complete.json').exists())
        self.assertFalse(Path(str(phy) + '_phyml_tree.txt').exists())
        # The explicit failure remains diagnostic after partial outputs are removed.
        self.assertIn('numerical_model_fit_failed', phyml.completion_error(phy, 100))

    def test_interrupted_tree_uses_recorded_runtime_not_idle_time(self):
        out, cmd, refs = self.fixture()
        phy = out / 'loci/haplotypes.phy'
        phy.write_text('4 4\nA AAAA\nB AAAT\nC AATT\nD TTTT\n')
        binary = self.root / 'bin/phyml'
        binary.parent.mkdir()
        launched = self.root / 'launched'
        binary.write_text(f'#!/bin/sh\ntouch {shlex.quote(str(launched))}\nexit 2\n')
        binary.chmod(0o755)
        status = Path(str(phy) + '.phyml.run.status.tsv')
        log = Path(str(phy) + '.phyml.log')
        env = dict(os.environ, PATH=str(binary.parent) + ':' + os.environ['PATH'], PHYML_RETRY_TIMEOUTS='FALSE')
        args = [sys.executable, str(ROOT / 'f/phyml.py'), 'run', '--phy', str(phy), '--scope', 'test', '--cpus', '1', '--timeout', '7200']
        for observed_seconds, expected_rc in [(7300, 124), (10, 2)]:
            phyml.state(phy, 'test', 'RUNNING', '')
            log.write_text('incomplete tree\n')
            start = time.time() - 10000
            os.utime(phy, (start - 1, start - 1))
            os.utime(status, (start, start))
            os.utime(log, (start + observed_seconds, start + observed_seconds))
            result = subprocess.run(args, env=env, capture_output=True, timeout=10)
            self.assertEqual(result.returncode, expected_rc, result.stderr.decode())
            self.assertEqual(launched.exists(), observed_seconds < 7200)
        # An explicit retry bypasses the historical-budget deferral.
        launched.unlink()
        phyml.state(phy, 'test', 'RUNNING', '')
        os.utime(phy, (start - 1, start - 1))
        os.utime(status, (start, start))
        os.utime(log, (start + 7300, start + 7300))
        result = subprocess.run(args, env=dict(env, PHYML_RETRY_TIMEOUTS='TRUE'), capture_output=True, timeout=10)
        self.assertEqual(result.returncode, 2)
        self.assertTrue(launched.exists())


if __name__ == '__main__':
    unittest.main()
