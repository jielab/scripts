#!/usr/bin/env python3
"""Exercise the real copy CLI with isolated sources/destinations; no bsub calls."""
import contextlib
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


ENTRY = Path(__file__).resolve().parents[2] / 'github_copy.sh'


class GithubCopyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='github-copy-test-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / 'analysis'
        self.source.mkdir()
        self.env = {
            'PATH': os.environ['PATH'], 'HOME': os.environ['HOME'],
            'submit': 'hpc', 'DIR_DAT': str(self.root),
            'UTILITY_PYTHON': sys.executable, 'PYTHONDONTWRITEBYTECODE': '1',
        }
        self.write('allele3/runs/current/primary/primary_manifest.json', '{"status":"complete"}')
        self.write('allele3/runs/current/primary/primary_enrichment.tsv', 'family\toe\nERV1\t1.2\n')

    def write(self, relative, data):
        p = self.source / relative
        p.parent.mkdir(parents=True, exist_ok=True)
        if isinstance(data, bytes):
            p.write_bytes(data)
        else:
            p.write_text(data)
        return p

    def run_cli(self, *args, env=None, ok=True):
        result = subprocess.run(
            ['bash', str(ENTRY), *args], cwd=self.root,
            env={**self.env, **(env or {})}, capture_output=True, text=True, timeout=30,
        )
        if ok:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def copy(self, *args):
        return self.run_cli('allele3', '--check', 'size,type', '--submit', 'local',
                            '--dir-dest', str(self.root / 'github'), *args)

    def files(self):
        return {p.relative_to(self.root / 'github').as_posix()
                for p in (self.root / 'github/allele3').rglob('*') if p.is_file()}

    @contextlib.contextmanager
    def loaded_copy(self):
        embedded = ENTRY.read_text().split("<<'PY'\n", 1)[1].rsplit('\nPY', 1)[0]
        definitions = embedded.split('\nsignal.signal(signal.SIGTERM', 1)[0]
        namespace = {'__name__': 'github_copy_test'}
        with patch.dict(os.environ, {**self.env, 'ANALYSIS_ROOT': str(self.source), 'DIR_RUNTIME_TMP': str(self.root)}):
            exec(compile(definitions, str(ENTRY), 'exec'), namespace)
            namespace['SOURCE'] = self.source
            namespace['DEST'] = self.root / 'github'
            namespace['CHECKS'] = frozenset({'size', 'type'})
            with contextlib.redirect_stdout(io.StringIO()):
                yield namespace

    def test_hpc_defaults_and_read_only_preview(self):
        r = self.run_cli('allele3', '--check', 'size,type', '--dry-run')
        self.assertIn('Source: ' + str(self.source), r.stdout)
        self.assertIn('Destination: ' + str(self.root / 'github'), r.stdout)
        self.assertIn('Selected: 1 files', r.stdout)
        self.assertFalse((self.root / 'github').exists())
        self.assertFalse((self.source / 'github_copy').exists())

    def test_hpc_plan_preserves_exact_selection_and_worker_roots(self):
        self.run_cli('allele3', '--check', 'size,type', '--write-cmd')
        paths = list((self.source / 'github_copy/jobs').glob('*/plan.json'))
        self.assertEqual(len(paths), 1)
        plan = json.loads(paths[0].read_text())
        self.assertEqual(plan['mode'], 'hpc')
        job, = plan['jobs']
        self.assertEqual(job['env']['GITHUB_COPY_SOURCE_ROOT'], str(self.source))
        self.assertEqual(job['env']['GITHUB_COPY_DEST_ROOT'], str(self.root / 'github'))
        self.assertEqual(job['argv'][2:5], ['allele3', '--check', 'size,type'])
        self.assertEqual(job['argv'][-3:], ['--submit', 'hpc', '--execute'])
        self.assertEqual(job['resources']['cores'], 10)
        self.assertFalse((self.root / 'github').exists())

    def test_local_default_is_preserved(self):
        r = self.run_cli('allele3', '--submit', 'local', '--check', 'size,type', '--dry-run')
        self.assertIn('Destination: ' + str(self.root / 'github/analysis'), r.stdout)

    def test_source_and_destination_overrides(self):
        custom = self.root / 'custom-source'
        (custom / 'allele3').mkdir(parents=True)
        (custom / 'allele3/result.tsv').write_text('value\n3\n')
        target = self.root / 'explicit-destination'
        r = self.run_cli('allele3', '--check', 'size,type', '--submit', 'local',
                         '--dir', str(custom), '--dir-dest', str(target),
                         env={'GITHUB_COPY_DEST_ROOT': str(self.root / 'ignored-destination')})
        self.assertIn('COPY OK: 1 files', r.stdout)
        self.assertEqual((target / 'allele3/result.tsv').read_text(), 'value\n3\n')
        self.assertFalse((self.root / 'ignored-destination').exists())
        r = self.run_cli('allele3', '--dry-run', '--check', 'size,type', env={
            'GITHUB_COPY_SOURCE_ROOT': str(custom), 'GITHUB_COPY_DEST_ROOT': str(target)})
        self.assertIn('Source: ' + str(custom), r.stdout)
        self.assertIn('Destination: ' + str(target), r.stdout)

    def test_selection_copies_completed_results_without_changing_sources(self):
        wanted = {
            'allele3/runs/current/primary/primary_enrichment.tsv': 'family\toe\nERV1\t1.2\n',
            'allele3/runs/current/summary/results.tsv.gz': gzip.compress(b'family\tp\nERVK\t0.03\n'),
            'allele3/runs/current/genotypes/genotype_concordance.png': b'figure-preview',
            'allele3/runs/current/genotypes/cohort_retention.tsv': 'cohort\tn\nEUR\t500\n',
            'allele3/runs/current/genotypes/genotype_validation.tsv': 'variant\tstatus\nrs1\tpass\n',
            'allele3/runs/current/candidates/targets.tsv': 'chrom\tpos\n1\t123\n',
            'allele3/stores/gnomad_scan/summary/enrichment_summary.tsv': 'family\toe\nERVK\t1.4\n',
        }
        for name, value in wanted.items():
            self.write(name, value)
        for stage, manifest in [('summary', 'summary_manifest.json'),
                                ('genotypes', 'genotypes_manifest.json'),
                                ('candidates', 'candidate_manifest.json')]:
            self.write(f'allele3/runs/current/{stage}/{manifest}', '{"status":"complete"}')
        for name in ['allele3/jobs/job_ids.tsv', 'allele3/verification/old_result.tsv',
                     'allele3/stores/archaic_overlap/1kg/counts/chr1.regions.bed',
                     'allele3/stores/archaic_overlap/1kg/tracks/liftover.unmapped.bed',
                     'allele3/stores/archaic_overlap/1kg/enrichment_block_counts.tsv.gz',
                     'allele3/stores/gnomad_scan/summary/architecture.global_strata.tsv.gz',
                     'allele3/stores/gnomad_scan/shards/chr1/result.tsv',
                     'allele3/runs/current/summary/cache/result.tsv',
                     'allele3/runs/current/genotypes/sample_burden.tsv',
                     'allele3/runs/current/genotypes/balanced_samples.tsv',
                     'allele3/runs/current/candidates/cache/result.tsv',
                     'allele3/runs/current/primary/data.pgen',
                     'allele3/runs/current/primary/data.pgen.gz',
                     'allele3/.spatial.history/old/result.png']:
            self.write(name, 'not a selected result\n')
        self.write('allele3/runs/current/spatial/spatial_manifest.json', '{"status":"running"}')
        self.write('allele3/runs/current/spatial/unfinished.tsv', 'partial\n1\n')
        self.write('allele3/runs/current/primary/empty.tsv', '')
        large = self.write('allele3/runs/current/primary/large.tsv', '')
        with large.open('wb') as h:
            h.truncate(10_000_001)
        (self.source / 'allele3/linked.tsv').symlink_to(
            self.source / 'allele3/runs/current/primary/primary_enrichment.tsv')
        before = {str(p): (p.stat().st_size, p.stat().st_mtime_ns, hashlib.sha256(p.read_bytes()).hexdigest())
                  for p in self.source.rglob('*') if p.is_file()}
        other = self.root / 'github/other/keep.txt'
        other.parent.mkdir(parents=True)
        other.write_text('untouched')
        r = self.copy()
        self.assertIn('COPY OK:', r.stdout)
        self.assertEqual(self.files(), set(wanted))
        for name in wanted:
            self.assertEqual((self.source / name).read_bytes(), (self.root / 'github' / name).read_bytes())
        for name, expected in before.items():
            p = Path(name)
            self.assertEqual((p.stat().st_size, p.stat().st_mtime_ns,
                              hashlib.sha256(p.read_bytes()).hexdigest()), expected)
        self.assertEqual(other.read_text(), 'untouched')
        manifest = (self.root / 'github/allele3.files.lst').read_text()
        self.assertIn('# Checks: size,type', manifest)
        self.assertEqual({line for line in manifest.splitlines() if line and not line.startswith('#')}, set(wanted))

    def test_repeat_copy_refreshes_only_selected_project(self):
        self.copy()
        stale = self.root / 'github/allele3/stale.tsv'
        stale.write_text('old')
        p = self.write('allele3/runs/current/primary/primary_enrichment.tsv', 'family\toe\nERV1\t2\n')
        self.copy()
        self.assertFalse(stale.exists())
        self.assertEqual((self.root / 'github' / p.relative_to(self.source)).read_bytes(), p.read_bytes())

    def test_folder_scope_preserves_other_destination_results(self):
        self.copy()
        self.write('allele3/runs/current/summary/summary_manifest.json', '{"status":"complete"}')
        self.write('allele3/runs/current/summary/new.tsv', 'summary\n4\n')
        self.write('allele3/logs/pipeline.log', 'pipeline still running\n')
        self.write('allele3/runs/current/pipeline_status.json', '{"status":"running","stages":[]}')
        self.copy('--folder', 'runs/current/summary')
        self.assertEqual(len(self.files()), 2)
        r = self.run_cli('allele3', '--check', 'size,type', '--folder', 'runs/current/summary', '--dry-run')
        self.assertIn('Selected: 1 files', r.stdout)

    def test_scoped_unfinished_stage_is_skipped(self):
        self.write('allele3/runs/current/spatial/spatial_manifest.json', '{"status":"running"}')
        self.write('allele3/runs/current/spatial/result.tsv', 'partial\n1\n')
        r = self.run_cli('allele3', '--check', 'size,type', '--folder', 'runs/current/spatial', '--dry-run')
        self.assertIn('Selected: 0 files', r.stdout)
        self.assertIn('Excluded unfinished_stage: 1', r.stdout)

    def test_explicit_auxiliary_folder_scopes_are_still_excluded(self):
        for folder in ('jobs', 'verification', 'runs/current/summary/verification',
                       'stores/archaic_overlap/1kg/counts', 'stores/archaic_overlap/1kg/tracks'):
            with self.subTest(folder=folder):
                self.write('allele3/' + folder + '/helper.tsv', 'helper\n1\n')
                r = self.run_cli('allele3', '--check', 'size,type', '--folder', folder, '--dry-run')
                self.assertIn('Selected: 0 files', r.stdout)
                self.assertIn('Excluded intermediate_directory: 1', r.stdout)

    def test_compressed_tables_honor_optional_id_check(self):
        self.write('allele3/with_ids.tsv.gz', gzip.compress(b'sample_id\tvalue\n1234567\t1\n'))
        r = self.run_cli('allele3', '--check', 'size,type', '--dry-run')
        self.assertIn('Selected: 2 files', r.stdout)
        r = self.run_cli('allele3', '--check', 'size,type,id', '--dry-run')
        self.assertIn('Selected: 1 files', r.stdout)
        self.assertIn('Excluded ukb_id_column: 1', r.stdout)

    def test_invalid_checks_and_overlapping_roots_fail(self):
        self.run_cli('allele3', '--check', 'size,bogus', '--dry-run', ok=False)
        self.run_cli('allele3', '--check', 'size,type', '--dir-dest', str(self.source / 'copy'),
                     '--dry-run', ok=False)
        self.assertFalse((self.root / 'github').exists())

    def test_existing_git_checkout_is_not_replaced(self):
        target = self.root / 'github/allele3'
        (target / '.git').mkdir(parents=True)
        (target / 'keep').write_text('existing checkout')
        r = self.run_cli('allele3', '--check', 'size,type', '--submit', 'local',
                         '--dir-dest', str(self.root / 'github'), ok=False)
        self.assertIn('contains .git', r.stderr)
        self.assertEqual((target / 'keep').read_text(), 'existing checkout')

    def test_important_logs_prioritize_current_job_and_keep_only_recent_executions(self):
        self.write('allele3/runs/current/pipeline_status.json', json.dumps({
            'status': 'running', 'updated': '2026-11-06T01:00:00Z', 'batch_id': 'batch-a',
            'configuration': {'unused': 'do not export detailed verification data'},
            'stages': [{'stage': 'primary', 'status': 'complete', 'worker_exit_code': 0},
                       {'stage': 'spatial', 'status': 'running', 'lsf_job_id': '101'}],
        }))
        for i in range(1, 5):
            path = self.write(f'allele3/jobs/job{i}/logs/run.{100+i}.out', f'execution {i}\n')
            os.utime(path, (i * 100, i * 100))
        error = self.write('allele3/jobs/job4/logs/run.104.err', 'important warning\n')
        os.utime(error, (400, 400))
        for path in ('allele3/jobs/job4/plan.json', 'allele3/jobs/job4/helper.tsv',
                     'allele3/jobs/verification/logs/run.999.out',
                     'allele3/jobs/debug/logs/run.998.out',
                     'allele3/jobs/job4/logs/trace.log',
                     'allele3/verification/logs/run.log',
                     'allele3/runs/.hidden/pipeline_status.json'):
            self.write(path, 'auxiliary data')
        self.write('allele3/jobs/job4/logs/run.104.log', '')
        linked = self.source / 'allele3/jobs/job4/logs/linked.log'
        linked.symlink_to(self.source / 'allele3/jobs/job4/logs/run.104.out')
        self.copy()
        logs = {name for name in self.files() if '/run_logs/' in name}
        expected = {
            'allele3/run_logs/runs-current.status.log',
            'allele3/run_logs/jobs-job1/run.101.out.log',
            'allele3/run_logs/jobs-job3/run.103.out.log',
            'allele3/run_logs/jobs-job4/run.104.out.log',
            'allele3/run_logs/jobs-job4/run.104.err.log',
        }
        self.assertEqual(logs, expected)
        self.assertEqual(len(self.files()), 6)
        for name in self.files():
            self.assertFalse({'jobs', 'verification'} & set(Path(name).parts))
        summary = (self.root / 'github/allele3/run_logs/runs-current.status.log').read_text()
        self.assertIn('primary: complete', summary)
        self.assertIn('worker_exit_code=0', summary)
        self.assertIn('spatial: running | job=101', summary)
        self.assertNotIn('do not export detailed verification data', summary)
        warning = (self.root / 'github/allele3/run_logs/jobs-job4/run.104.err.log').read_text()
        self.assertIn('# Source: allele3/jobs/job4/logs/run.104.err', warning)
        self.assertIn('# Snapshot UTC:', warning)
        self.assertTrue(warning.endswith('important warning\n'))

    def test_long_logs_retain_head_tail_and_leave_source_unchanged(self):
        body = 'FIRST configuration\n' + '运行进度 checkpoint\n' * 90_000 + 'LAST fatal diagnostic\n'
        source = self.write('allele3/logs/pipeline.log', body)
        before = (source.stat().st_mtime_ns, hashlib.sha256(source.read_bytes()).hexdigest())
        self.copy()
        copied = self.root / 'github/allele3/run_logs/root-logs/pipeline.log.log'
        text = copied.read_text()
        self.assertLessEqual(copied.stat().st_size, 1_000_000)
        self.assertIn('TRUNCATED', text)
        self.assertIn('FIRST configuration', text)
        self.assertIn('LAST fatal diagnostic', text)
        self.assertEqual((source.stat().st_mtime_ns, hashlib.sha256(source.read_bytes()).hexdigest()), before)

    def test_appending_log_does_not_invalidate_frozen_copy(self):
        source = self.write('allele3/logs/pipeline.out', 'started\n')
        with self.loaded_copy() as namespace:
            selected = namespace['scan'](['allele3'])
            with source.open('a') as handle:
                handle.write('later progress\n')
            namespace['copy_files'](['allele3'], selected)
        snapshot = (self.root / 'github/allele3/run_logs/root-logs/pipeline.out.log').read_text()
        self.assertTrue(snapshot.endswith('started\n'))
        self.assertNotIn('later progress', snapshot)
        self.assertEqual(source.read_text(), 'started\nlater progress\n')

    def test_log_id_checks_and_malformed_status_do_not_block_results(self):
        self.write('allele3/runs/current/pipeline_status.json', '{partially written')
        self.write('allele3/logs/participants.log', 'sample_id\tvalue\n1234567\t2\n')
        preview = self.run_cli('allele3', '--check', 'size,type', '--dry-run')
        self.assertIn('Selected: 2 files', preview.stdout)
        checked = self.run_cli('allele3', '--check', 'size,type,id', '--dry-run')
        self.assertIn('Selected: 1 files', checked.stdout)
        self.assertIn('Skipped status log:', checked.stdout)
        self.assertIn('Skipped execution log:', checked.stdout)

    def test_log_snapshots_share_project_size_budget(self):
        self.write('allele3/logs/pipeline.log', 'started\n')
        with self.loaded_copy() as namespace:
            namespace['MAX_PROJECT_BYTES'] = 100
            selected = namespace['scan'](['allele3'])
        self.assertEqual(set(selected), {'allele3/runs/current/primary/primary_enrichment.tsv'})
        self.assertEqual(namespace['LOG_SNAPSHOTS'], {})

    def test_stage_restart_aborts_before_replacing_destination(self):
        self.copy()
        original = (self.root / 'github/allele3.files.lst').read_bytes()
        with self.loaded_copy() as namespace:
            selected = namespace['scan'](['allele3'])
            self.write('allele3/runs/current/primary/primary_manifest.json', '{"status":"running"}')
            with self.assertRaises(namespace['SourceChanged']):
                namespace['copy_files'](['allele3'], selected)
        self.assertEqual((self.root / 'github/allele3.files.lst').read_bytes(), original)
        self.assertEqual(len(self.files()), 1)


if __name__ == '__main__':
    unittest.main(verbosity=2)
