"""Use account quotas, including Lustre scratch, instead of filesystem totals."""
import json
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'f'))
import ld_pipeline as ld


class StorageTests(unittest.TestCase):
    def setUp(self):
        self.pipeline = object.__new__(ld.Pipeline)
        self.pipeline.min_free = 0

    def check(self, path, outputs):
        with patch.object(ld.Path, 'stat', return_value=SimpleNamespace(st_gid=os.getgid(), st_mode=0)), \
             patch.object(ld.shutil, 'disk_usage', return_value=SimpleNamespace(free=170 * 1024**4)), \
             patch.object(ld.subprocess, 'run', side_effect=[SimpleNamespace(stdout=x) for x in outputs]) as run:
            value = self.pipeline.check_free_space(path)
            return value, run.call_args_list

    @staticmethod
    def lustre():
        return json.dumps({'filesystems': [{'source': 'server:/lustre1', 'fstype': 'lustre', 'target': '/scratch'}]})

    def test_scratch_user_quota_limits_unlimited_group_and_large_filesystem(self):
        value, calls = self.check('/scratch/test', [self.lustre(),
            '/scratch 553656764 5368709120 5410652160 - 629 100000 110000 -\n',
            '/scratch 553656764 0 0 - 629 0 0 -\n'])
        self.assertEqual(value, (5368709120 - 553656764) * 1024)
        self.assertEqual(calls[1].args[0][:3], ['lfs', 'quota', '-u'])
        self.assertEqual(calls[2].args[0][:3], ['lfs', 'quota', '-g'])

    def test_stricter_group_quota_wins(self):
        value, _ = self.check('/scratch/test', [self.lustre(),
            '/scratch 100 1000 1100 - 0 0 0 -', '/scratch 100* 150 160 - 0 0 0 -'])
        self.assertEqual(value, 50 * 1024)

    def test_quota_failure_does_not_fall_back_to_filesystem_free(self):
        with self.assertRaisesRegex(ld.Stop, 'Cannot verify scratch user/group quota'):
            self.check('/scratch/test', [self.lustre(), 'unavailable quota information'])

    def test_reserve_uses_quota_remaining(self):
        self.pipeline.min_free = 100 * 10**9
        with self.assertRaisesRegex(ld.Stop, 'Free disk below MIN_FREE_GB'):
            self.check('/scratch/test', [self.lustre(),
                '/scratch 900000 1000000 1100000 - 0 0 0 -', '/scratch 900000 0 0 - 0 0 0 -'])

    def test_work_gpfs_queries_real_filesystem_name(self):
        mount = json.dumps({'filesystems': [{'source': 'fs01', 'fstype': 'gpfs', 'target': '/work'}]})
        quota = ('mmlsquota::HEADER:filesystemName:blockUsage:blockQuota:blockLimit:blockInDoubt\n'
                 'mmlsquota:::fs01:10:100:120:5\n')
        value, calls = self.check('/work/test', [mount, quota])
        self.assertEqual(value, 85 * 1024)
        self.assertEqual(calls[1].args[0][-1], 'fs01')

    def test_moved_cache_resolves_old_output_symlink(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            old = root / 'old/ld'
            old.mkdir(parents=True)
            target = root / 'scratch/analysis/ld'
            target.mkdir(parents=True)
            (target / 'test.pgen').write_bytes(b'data')
            (old / '38').symlink_to(target, target_is_directory=True)
            record = dict(path=str(old / '38/test.pgen'), size=4, mtime_ns=1)
            with patch.object(ld, 'configured_path', return_value=str(old)), \
                 patch.dict(os.environ, SCRATCH_ANALYSIS_ROOT=str(target.parent)):
                result = ld.relocated_output_stamp(record)
            self.assertEqual(result, dict(record, path=str(target / 'test.pgen')))

    def test_verified_move_records_precision_loss_but_detects_later_changes(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / 'scratch'
            target.mkdir()
            old = root / 'old'
            old.symlink_to(target, target_is_directory=True)
            path = target / 'test.pgen'
            path.write_bytes(b'data')
            after = ld.stamp(path)
            before = dict(after, path=str(old / 'test.pgen'), mtime_ns=after['mtime_ns'] + 123)
            (target / 'LD_OUTPUT_MIGRATION.json').write_text(json.dumps(dict(status='complete',
                checksum_verified=True, files={before['path']: dict(before=before, after=after)})))
            self.assertEqual(ld.relocated_output_stamp(before), after)
            os.utime(path, ns=(after['mtime_ns'], after['mtime_ns'] + 2000000000))
            self.assertNotEqual(ld.relocated_output_stamp(before), ld.stamp(path))


if __name__ == '__main__':
    unittest.main()
