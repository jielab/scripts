import sys
from pathlib import Path
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hpc_query import discover_samples


class DiscoveryTests(unittest.TestCase):
    def test_incoming_downloads_are_skipped_only_when_requested(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data = root / 'dataset'
            data.mkdir()
            (root / 'empty').mkdir()
            for name in ['done_1.fastq.gz', 'done_2.fastq.gz', 'partial_1.fastq.gz',
                         'partial_2.fastq.gz', 'missing_1.fastq.gz']:
                (data / name).write_bytes(b'fixture')
            (data / 'partial_1.fastq.gz.aria2').touch()
            cfg = {'reads_raw': str(root)}
            result = discover_samples(cfg, ready_only=True)
            self.assertEqual([v['sample_id'] for v in result], ['dataset/done'])
            with self.assertRaises(RuntimeError):
                discover_samples(cfg)
            with self.assertRaises(RuntimeError):
                discover_samples(cfg, only=['dataset/partial'], ready_only=True)
            with self.assertRaises(RuntimeError):
                discover_samples(cfg, only=['dataset/missing'], ready_only=True)


if __name__ == '__main__':
    unittest.main()
