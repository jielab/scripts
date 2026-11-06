import fcntl
import gzip
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from cache_cleanup import clean_csx
from csx.ld_cache import cached_ld_file, fingerprint


class CacheCleanupTests(unittest.TestCase):
    def test_ld_staging_released_on_success_and_exception(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            source = root / 'reference.hdf5'
            source.write_bytes(b'permanent reference')
            with patch.dict(os.environ, {'GRID_CSX_LD_CACHE': str(root / 'cache'),
                                         'GRID_CSX_KEEP_LD_CACHE': '0'}):
                with cached_ld_file(source) as staged:
                    self.assertEqual(staged.read_bytes(), source.read_bytes())
                self.assertFalse(staged.exists())
                with self.assertRaises(RuntimeError):
                    with cached_ld_file(source) as staged:
                        raise RuntimeError('injected error')
                self.assertFalse(staged.exists())
                self.assertEqual(source.read_bytes(), b'permanent reference')

    def test_existing_ld_and_merged_scores_preview_release_and_active_locks(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            source = root / 'reference.hdf5'
            source.write_bytes(b'permanent reference')
            ld = root / 'cache/csx-ld'
            ld.mkdir(parents=True)
            cached = ld / 'ld-test.hdf5'
            cached.write_bytes(source.read_bytes())
            cached.with_suffix('.json').write_text(json.dumps(fingerprint(source)))
            run = root / 'cache/tables/key/scores/run'
            run.mkdir(parents=True)
            (run / 'config.json').write_text('{}')
            with gzip.open(run / 'AFR.tsv.gz', 'wb') as stream:
                stream.write(b'eid\tCSX_AFR\n001\t0.5\n')
            score = run / 'AFR.chr1.sscore'
            score.write_bytes(b'chromosome score')
            unfinished = run / 'EUR.chr1.sscore'
            unfinished.write_bytes(b'unmerged score')
            size = cached.stat().st_size + score.stat().st_size
            self.assertEqual(clean_csx(root), size)
            with cached.with_suffix('.lock').open('a') as lock:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                self.assertEqual(clean_csx(root, apply=True), len(b'chromosome score'))
                self.assertTrue(cached.exists())
            self.assertEqual(clean_csx(root, apply=True), len(b'permanent reference'))
            self.assertTrue(source.exists())
            self.assertTrue(unfinished.exists())
            self.assertFalse(score.exists())
            self.assertEqual(clean_csx(root, apply=True), 0)

    def test_corrupt_merge_and_unbound_posterior_are_retained(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            run = root / 'cache/tables/key/meta/run'
            run.mkdir(parents=True)
            (run / 'config.json').write_text('{}')
            (run / 'combined.gz').write_bytes(b'interrupted merge')
            score = run / 'chr1.sscore'
            score.write_bytes(b'original scores')
            posterior = root / 'cache/tables/other/posterior_scores/run'
            (posterior / 'chr1').mkdir(parents=True)
            (posterior / 'config.json').write_text('{}')
            moments = posterior / 'chr1/moments.h5'
            moments.write_bytes(b'unpublished moments')
            self.assertEqual(clean_csx(root, apply=True), 0)
            self.assertTrue(score.exists())
            self.assertTrue(moments.exists())


if __name__ == '__main__':
    unittest.main()
