"""Frozen plans must detect same-size changes while jobs are queued."""
import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import audit_integrity
import hpc_cli as cli
import hpc_common
import hpc_query
import hpc_compare


class PlanRuntimeIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        with patch.dict(os.environ, {'MEMORY_DEFAULT_GIB': '32'}):
            self.cfg = hpc_common.load_config()
        self.cfg.update({key: str(self.root / key) for key in hpc_common.PATH_KEYS})
        self.cfg.update(root=str(self.root), classifier_engine='kraken', classify_workers=1)
        self.cfg['db_dir'] = str(Path(self.cfg['library_dir']) / 'db')
        self.cfg['kunpeng_db_dir'] = str(Path(self.cfg['library_dir']) / 'kunpeng')
        self.cfg['prepared_dir'] = str(Path(self.cfg['library_dir']) / 'prepared')
        self.cfg['reference_plan'] = str(Path(self.cfg['library_dir']) / 'reference_plan.json')
        self.cfg['taxonomy_dir'] = str(Path(self.cfg['metadata_dir']) / 'taxonomy')
        self.db = Path(self.cfg['db_dir'])
        self.tax = self.db / 'taxonomy'
        self.tax.mkdir(parents=True)
        for name in ('hash.k2d', 'taxo.k2d', 'opts.k2d'):
            (self.db / name).write_bytes(b'fixture!')
        for name in ('nodes.dmp', 'names.dmp'):
            (self.tax / name).write_bytes(b'fixture!')
        self.snapshot = {'identity': 'fixture', 'db_dir': str(self.db),
                         'files': {n: {} for n in ('hash.k2d', 'taxo.k2d', 'opts.k2d')},
                         'taxonomy_dir': str(self.tax),
                         'taxonomy_files': {n: {} for n in ('nodes.dmp', 'names.dmp')}}
        self.fastq = self.root / 'sample.fastq'
        self.fastq.write_text('@read\nACGT\n+\nIIII\n')
        self.sample = {'sample_id': 'data/sample', 'dataset': 'data', 'sample': 'sample',
                       'r1': str(self.fastq), 'r2': None}

    def _generated(self):
        with patch.object(cli, 'engine_snapshot', return_value=self.snapshot), \
             patch.object(hpc_query, 'discover_samples', return_value=[self.sample]), \
             contextlib.redirect_stdout(io.StringIO()):
            path = cli.generate_split(self.cfg)
        return path, hpc_common.load_config(path.parent / 'config.json')

    def _assert_worker_refuses(self, path, cfg, message):
        args = argparse.Namespace(command='split-worker', allow_local=True,
                                  job_plan=str(path), batch_id=1, offset=0)
        with patch.object(cli, 'engine_snapshot', return_value=self.snapshot), \
             patch.object(hpc_query, 'run_batch') as run, \
             contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(RuntimeError, message):
                cli.run_worker(args, cfg)
        run.assert_not_called()

    def test_same_size_fastq_replaced_after_generation_is_refused(self):
        path, cfg = self._generated()
        self.fastq.write_text('@read\nTGCA\n+\nIIII\n')
        self._assert_worker_refuses(path, cfg, 'FASTQ changed after job generation')

    def test_same_size_database_replaced_after_generation_is_refused(self):
        path, cfg = self._generated()
        (self.db / 'hash.k2d').write_bytes(b'changed!')
        self._assert_worker_refuses(path, cfg, 'runtime state changed after job generation')

    def test_summary_refuses_changed_taxonomy_for_kraken_and_kunpeng_only_cases(self):
        for engine in ('kraken', 'kunpeng_full'):
            with self.subTest(engine=engine):
                (self.db / 'hash_1.k2d').write_bytes(b'fixture!')
                snapshot = dict(self.snapshot)
                if engine != 'kraken':
                    snapshot.update(identity='converted', source_database_identity='fixture',
                                    files={n: {} for n in ('hash_1.k2d', 'taxo.k2d', 'opts.k2d')})
                with patch.object(cli, 'engine_snapshot', return_value=snapshot), \
                     patch.object(hpc_query, 'discover_samples', return_value=[self.sample]), \
                     contextlib.redirect_stdout(io.StringIO()):
                    path = cli.generate_compare(self.cfg, methods=[engine], repeats=1, limit=1)
                cfg = hpc_common.load_config(path.parent / 'config.json')
                target = self.tax / 'names.dmp'
                target.write_bytes(b'changed!')
                args = argparse.Namespace(command='compare-summary-worker', allow_local=True, job_plan=str(path))
                with patch.object(cli, 'engine_snapshot', return_value=snapshot) as get_snapshot, \
                     patch.object(hpc_compare, 'summarize_cases') as summarize:
                    with self.assertRaisesRegex(RuntimeError, 'runtime state changed before comparison summary'):
                        cli.run_worker(args, cfg)
                summarize.assert_not_called()
                self.assertEqual(get_snapshot.call_args.args[0]['classifier_engine'], engine)

    def test_optional_seal_is_bound_to_library_and_database_identity(self):
        seal = self.root / 'seal.json'
        seal.write_text('{}')
        cfg = dict(self.cfg, integrity_seal=str(seal))
        (self.db / 'database.json').write_text(json.dumps({'identity': 'fixture'}))
        verified = {'paths': {key: cfg[key] for key in (
            'reference_raw', 'metadata_dir', 'library_dir', 'db_dir', 'kunpeng_db_dir')},
            'indexes': {'kraken': {'status': 'complete', 'directory': str(self.db), 'identity': 'fixture'}}}
        with patch.object(audit_integrity, 'verify_seal', return_value=verified):
            observed = cli._integrity_seal(cfg, required_indexes=['kraken'])
            self.assertEqual(observed['path'], str(seal))
            verified['paths']['db_dir'] = str(self.root / 'different-db')
            with self.assertRaisesRegex(RuntimeError, 'different db_dir'):
                cli._integrity_seal(cfg, required_indexes=['kraken'])
            verified['paths']['db_dir'] = cfg['db_dir']
            verified['indexes']['kraken']['identity'] = 'other'
            with self.assertRaisesRegex(RuntimeError, 'database identity'):
                cli._integrity_seal(cfg, required_indexes=['kraken'])

    def test_references_only_seal_cannot_authorize_unsealed_query_index(self):
        seal = self.root / 'seal.json'
        seal.write_text('{}')
        cfg = dict(self.cfg, integrity_seal=str(seal))
        verified = {'paths': {key: cfg[key] for key in (
            'reference_raw', 'metadata_dir', 'library_dir', 'db_dir', 'kunpeng_db_dir')},
            'indexes': {'kraken': {'status': 'not_complete'}}}
        with patch.object(audit_integrity, 'verify_seal', return_value=verified):
            with self.assertRaisesRegex(RuntimeError, 'does not cover'):
                cli._integrity_seal(cfg, required_indexes=['kraken'])

    def test_seal_changes_are_rejected_before_use(self):
        seal = self.root / 'seal.json'
        seal.write_text('{}')
        cfg = dict(self.cfg, integrity_seal=str(seal))
        with patch.object(audit_integrity, 'verify_seal') as verify:
            with self.assertRaisesRegex(RuntimeError, 'changed after job generation'):
                cli._integrity_seal(cfg, {'path': str(seal), 'sha256': '0' * 64})
        verify.assert_not_called()


if __name__ == '__main__':
    unittest.main()
