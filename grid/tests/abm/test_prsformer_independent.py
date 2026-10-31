"""Independent trait scheduling, publication, resume and shared-cache consumers."""
import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import numpy as np
import pandas as pd

import prsformer
import prsformer_resume as original
from abm_core import prsformer_adapter as calibrated
from abm_core.model import write_json, file_sha256
from benchmark_paths import prsformer2_current, prsformer_publication
from browser_index import Index
from dosage_storage import _benchmark_evaluations
from final_compare import current_tables

ROOT = Path(__file__).resolve().parents[2]
TRAITS = ['height', 'ldl', 't2dm']


def cache_fixture(root):
    cache = root / 'prsformer/cache'
    cache.mkdir(parents=True)
    ids = [str(i) for i in range(101, 107)]
    d = pd.DataFrame(dict(eid=ids, group=ids, target='EUR',
        split=['train', 'validation', 'test', 'test', 'test', 'test'],
        height=[170., 171., 172., 173., 174., 175.],
        ldl=[2., 3., 4., 5., 6., 7.], t2dm=[0., 1., 0., 1., 0., 1.]))
    d.to_csv(cache / 'data.tsv.gz', sep='\t', index=False)
    pd.DataFrame(dict(SNP=['rs1'])).to_csv(cache / 'variants.tsv.gz', sep='\t', index=False)
    np.save(cache / 'genotypes.npy', np.zeros((len(d), 1), dtype=np.float16))
    write_json(cache / 'prepare.json', dict(signature='fixture', shape=[len(d), 1],
        options=dict(covariates='none', height_col='height', ldl_col='ldl', t2dm_col='t2dm'),
        sha256={name: file_sha256(cache / name) for name in ('data.tsv.gz', 'variants.tsv.gz')}))
    return cache


def evaluated_run(run, cache, trait):
    run.mkdir(parents=True, exist_ok=True)
    (run / 'model.pt').write_bytes(('model for ' + trait).encode())
    (run / 'training.pt').write_bytes(('checkpoint for ' + trait).encode())
    d = pd.read_csv(cache / 'data.tsv.gz', sep='\t', dtype={'eid': str})
    d = d.loc[d.split.eq('test')].copy()
    prediction = [.2, .7, .3, .8] if trait == 't2dm' else d[trait] + .1
    pd.DataFrame(dict(eid=d.eid, target=d.target, split='test', trait=trait,
        y=d[trait], baseline=0., prediction=prediction, genetic_score=prediction)).to_csv(
            run / 'test_predictions.tsv.gz', sep='\t', index=False)
    pd.DataFrame([dict(trait=trait, target='EUR', split='test',
        metric='AUC' if trait == 't2dm' else 'predictive_R2', value=.8, n=len(d))]).to_csv(
            run / 'metrics.tsv', sep='\t', index=False)
    pd.DataFrame(dict(trait=trait, epoch=[1, 1], split=['train', 'validation'], loss=[.5, .6])).to_csv(
        run / 'history.tsv', sep='\t', index=False)
    meta = json.loads((cache / 'prepare.json').read_text())
    model, geno = (run / 'model.pt').stat(), (cache / 'genotypes.npy').stat()
    write_json(run / 'evaluation.json', dict(format_version=1, preparation_signature=meta['signature'],
        inputs=dict(data_sha256=meta['sha256']['data.tsv.gz'], variants_sha256=meta['sha256']['variants.tsv.gz'],
                    genotypes=dict(size=geno.st_size, mtime_ns=geno.st_mtime_ns, shape=meta['shape'])),
        outputs_sha256={name: file_sha256(run / name) for name in
                       ('metrics.tsv', 'history.tsv', 'test_predictions.tsv.gz')},
        model=dict(size=model.st_size, mtime_ns=model.st_mtime_ns)))


def snapshot(folder):
    return {p.relative_to(folder): (p.read_bytes(), p.stat().st_mtime_ns)
            for p in folder.rglob('*') if p.is_file() and not p.name.endswith('.lock')}


def value(command, flag):
    return command[command.index(flag) + 1]


class IndependentTraitsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.cache = cache_fixture(self.root)
        self.frozen = snapshot(self.cache)
        self.output = io.StringIO()

    def test_original_failure_resume_skip_and_report_repair_are_trait_local(self):
        calls, fail = [], ['ldl']
        legacy = self.root / 'prsformer/run'
        legacy.mkdir()
        (legacy / 'model.pt').write_bytes(b'old joint model must remain untouched')
        old = snapshot(legacy)

        def execute(command, **kwargs):
            trait = value(command, '--traits')
            run = Path(value(command, '--prsformer-run-dir'))
            base = Path(value(command, '--prsformer-publication-dir'))
            self.assertEqual(run, self.root / 'prsformer' / trait / 'run')
            self.assertEqual(base, run.parent)
            self.assertEqual(value(command, '--prsformer-root'), str(self.root))
            self.assertEqual(kwargs['env']['PRSFORMER_KEEP_GENOTYPES'], '1')
            mode = value(command, '--prsformer-mode') if '--prsformer-mode' in command else 'all'
            calls.append((trait, mode, '--resume' in command))
            if mode != 'report':
                run.mkdir(parents=True, exist_ok=True)
                (run / 'model.pt').write_bytes(trait.encode())
                (run / 'training.pt').write_bytes(b'unfinished')
                if fail == [trait]:
                    fail.clear()
                    raise subprocess.CalledProcessError(7, command)
                evaluated_run(run, self.cache, trait)
            if mode in ('all', 'report'):
                prsformer.main(['publish', '--run-dir', str(run), '--cache-dir', str(self.cache),
                    '--output-root', str(self.root), '--publication-dir', str(base),
                    '--score-dir', value(command, '--prsformer-score-dir'), '--traits', trait,
                    *(['--replace'] if mode == 'report' else [])])

        with patch('abm_core.cold_start.benchmark_device', return_value='cuda'), \
             patch('abm_core.cold_start.ensure_inputs') as prepare, \
             patch.object(original, 'subprocess', SimpleNamespace(run=execute)), \
             contextlib.redirect_stdout(self.output):
            argv = ['--analysis-root', str(self.root), '--traits', ','.join(TRAITS)]
            with self.assertRaises(subprocess.CalledProcessError):
                original.dispatch(ROOT, argv)
            self.assertEqual(calls, [('height', 'all', False), ('ldl', 'all', False)])
            self.assertFalse((self.root / 'prsformer/t2dm').exists())
            saved = snapshot(self.root / 'prsformer/height')
            calls.clear()
            original.dispatch(ROOT, argv)
            self.assertEqual(calls, [('ldl', 'train', True), ('ldl', 'report', False), ('t2dm', 'all', False)])
            self.assertEqual(saved, snapshot(self.root / 'prsformer/height'))
            calls.clear()
            # Separate invocations reuse exactly the same per-trait results.
            for trait in reversed(TRAITS):
                original.dispatch(ROOT, ['--analysis-root', str(self.root), '--traits', trait])
            self.assertEqual(calls, [])
            for trait in TRAITS:
                self.assertEqual(original.plan(self.root, trait), 'skip')
                roster, scores = prsformer_publication(self.root, trait)
                self.assertTrue(roster.is_file() and scores.is_file())
            (self.root / 'prsformer/height/performance.png').unlink()
            original.dispatch(ROOT, argv)
            self.assertEqual(calls, [('height', 'report', False)])
            self.assertEqual(prepare.call_count, 6)
        self.assertEqual(snapshot(legacy), old)
        self.assertEqual(snapshot(self.cache), self.frozen)
        tables, _ = current_tables(self.root, TRAITS)
        self.assertEqual(set(pd.concat(tables).training_mode), {'independent_trait'})

    def test_calibrated_traits_share_cache_and_resume_without_touching_other_models(self):
        self._assert_calibrated_traits_resume(keep_shared=False)

    def test_calibrated_traits_retain_shared_cache_after_resume_and_skip(self):
        self._assert_calibrated_traits_resume(keep_shared=True)

    def _assert_calibrated_traits_resume(self, keep_shared):
        calls, fail = [], ['ldl']
        def execute(command, **kwargs):
            trait, run = value(command, '--traits'), Path(value(command, '--out-dir'))
            calls.append((trait, '--resume' in command))
            self.assertEqual(run.parent.parent, self.root / 'prsformer2' / trait)
            self.assertEqual(value(command, '--genotypes'), str(self.cache / 'genotypes.npy'))
            (run / 'model.pt').write_bytes(trait.encode())
            (run / 'training.pt').write_bytes(b'unfinished')
            if fail == [trait]:
                fail.clear()
                raise subprocess.CalledProcessError(7, command)
            evaluated_run(run, self.cache, trait)
        # Set the policy explicitly: the public shell entry exports retention,
        # while direct adapter calls may release verified shared inputs.
        with patch.dict(os.environ, GRID_KEEP_SHARED_INPUTS='1' if keep_shared else '0'), \
             patch.object(calibrated, 'benchmark_device', return_value='cuda'), \
             patch.object(calibrated, 'ensure_inputs') as prepare, \
             patch.object(calibrated, 'partition_calibration', side_effect=lambda roles, seed: roles), \
             patch.object(calibrated, 'subprocess', SimpleNamespace(run=execute)), \
             patch('dosage_storage.release_dosages') as cleanup, \
             contextlib.redirect_stdout(self.output):
            argv = ['--analysis-root', str(self.root), '--traits', ','.join(TRAITS)]
            with self.assertRaises(subprocess.CalledProcessError):
                calibrated.main(argv)
            self.assertEqual(calls, [('height', False), ('ldl', False)])
            cleanup.assert_not_called()
            saved_run = prsformer2_current(self.root, 'height')
            saved = snapshot(saved_run)
            calls.clear()
            calibrated.main(argv)
            self.assertEqual(calls, [('ldl', True), ('t2dm', False)])
            self.assertEqual(snapshot(saved_run), saved)
            self.assertEqual(cleanup.call_count, 0 if keep_shared else 1)
            calls.clear()
            pointers = {trait: prsformer2_current(self.root, trait) for trait in TRAITS}
            for index, trait in enumerate(reversed(TRAITS), 2):
                calibrated.main(['--analysis-root', str(self.root), '--traits', trait])
                self.assertEqual(prsformer2_current(self.root, trait), pointers[trait])
                self.assertEqual(cleanup.call_count, 0 if keep_shared else index)
            self.assertEqual(calls, [])
            self.assertEqual(prepare.call_count, 5)
        self.assertEqual(snapshot(self.cache), self.frozen)
        tables, status = current_tables(self.root, TRAITS)
        d = pd.concat(tables)
        self.assertEqual(len(d), 3)
        self.assertEqual(set(d.trait), set(TRAITS))
        self.assertTrue(d.training_mode.eq('independent_trait').all())
        Index(self.root, self.root / 'final').run(TRAITS)
        catalogue = pd.read_csv(self.root / 'final/catalogue.tsv', sep='\t')
        self.assertEqual(set(catalogue.loc[catalogue.module.eq('prsformer2'), 'trait']), set(TRAITS))
        self.assertFalse(catalogue.source.str.contains('_LOCAL_ONLY').any())
        # Cleanup must account for every evaluation, including nested trait runs.
        self.assertEqual(len(_benchmark_evaluations(self.root, self.cache)), 3)
        unfinished = self.root / 'prsformer/ldl/run'
        unfinished.mkdir(parents=True)
        (unfinished / 'training.pt').write_bytes(b'unfinished')
        self.assertIsNone(_benchmark_evaluations(self.root, self.cache))
        (unfinished / 'training.pt').unlink()
        unfinished = self.root / 'prsformer2/ldl/interrupted/_LOCAL_ONLY'
        unfinished.mkdir(parents=True)
        (unfinished / 'training.pt').write_bytes(b'unfinished')
        self.assertIsNone(_benchmark_evaluations(self.root, self.cache))

    def test_partial_independent_layout_does_not_fall_back_to_joint_results(self):
        evaluated_run(self.root / 'prsformer/run', self.cache, 'ldl')
        evaluated_run(self.root / 'prsformer/height/run', self.cache, 'height')
        legacy = self.root / 'prsformer2/old_joint'
        evaluated_run(legacy / '_LOCAL_ONLY', self.cache, 'ldl')
        write_json(legacy / 'RUN_MANIFEST.json', {})
        (legacy / '_LOCAL_ONLY/roles.tsv.gz').write_bytes(b'roles')
        calibrated.publish(legacy, self.cache, ['ldl'])
        (self.root / 'prsformer2/height').mkdir()
        tables, status = current_tables(self.root, TRAITS)
        self.assertEqual(pd.concat(tables).trait.tolist(), ['height'])
        self.assertTrue(any(s['trait'] == 'ldl' and s['method'] == 'PRSformer' and s['status'] == 'not_completed' for s in status))
        self.assertIsNone(prsformer2_current(self.root, 'ldl'))
        evaluated_run(self.root / 'prsformer/height/run', self.cache, 'ldl')
        with self.assertRaisesRegex(ValueError, 'missing requested traits'):
            current_tables(self.root, ['height'])

    def test_joint_only_outputs_are_ignored_before_any_independent_run_exists(self):
        from twin import frozen_comparison_roster, verified_prsformer
        from main import add_prsformer
        evaluated_run(self.root / 'prsformer/run', self.cache, 'height')
        base = self.root / 'prsformer'
        (base / 'split.rds').write_bytes(b'old roster must never be read')
        old_scores = base / 'scores/height'
        old_scores.mkdir(parents=True)
        (old_scores / 'predictions.rds').write_bytes(b'old prediction must never be read')
        (base / 'performance.tsv').write_text('method\tvalue\nold_joint\t1\n')
        nested = self.root / 'benchmark/prsformer'
        nested.mkdir(parents=True)
        (nested / 'performance.tsv').write_text('method\tvalue\nold_joint\t1\n')
        legacy = self.root / 'prsformer2/old_joint'
        evaluated_run(legacy / '_LOCAL_ONLY', self.cache, 'height')
        write_json(legacy / 'RUN_MANIFEST.json', {})
        (legacy / '_LOCAL_ONLY/roles.tsv.gz').write_bytes(b'roles')
        calibrated.publish(legacy, self.cache, ['height'])
        tables, status = current_tables(self.root, TRAITS)
        self.assertEqual(tables, [])
        self.assertTrue(all(row['status'] == 'not_completed' for row in status if row['method'] != 'DiscoDivas'))
        self.assertIsNone(prsformer2_current(self.root, 'height'))
        self.assertEqual(len(frozen_comparison_roster(self.root)), 6)
        prediction, _ = verified_prsformer(self.root, 'height', None, None, None, [], None)
        self.assertIsNone(prediction)
        _, report = add_prsformer(SimpleNamespace(prsformer_root=self.root), 'height', None,
                                 None, pd.DataFrame())
        self.assertEqual(report['status'], 'not_run_on_common_split')
        Index(self.root, self.root / 'final').run(TRAITS)
        catalogue = pd.read_csv(self.root / 'final/catalogue.tsv', sep='\t')
        self.assertTrue(catalogue.empty)

    def test_current_pointer_cannot_redirect_to_joint_or_other_trait_directory(self):
        base = self.root / 'prsformer2/height'
        base.mkdir(parents=True)
        for name in ('../old_joint', '../ldl/model', '/tmp/old_joint'):
            with self.subTest(name=name):
                write_json(base / 'CURRENT.json', {'run': name})
                with self.assertRaisesRegex(ValueError, 'inside its trait directory'):
                    prsformer2_current(self.root, 'height')


if __name__ == '__main__':
    unittest.main()
