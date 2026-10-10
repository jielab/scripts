"""Dataset destinations and recovery after moving an existing result root."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('gu_path_test', ROOT / 'f/0.common.py')
common = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = common
spec.loader.exec_module(common)


class OutputPaths(unittest.TestCase):
    def resolve(self, *args):
        return subprocess.run(['bash', '-c', 'source "$1"; shift; gu_result_root "$@"',
                               'test', str(ROOT / 'f/0.common.sh'), *args],
                              capture_output=True, text=True)

    def test_dataset_defaults_and_literal_postfix(self):
        for args, expected in [
            (('1kg',), '/mnt/d/analysis/gu/1kg'),
            (('ukb',), '/mnt/d/analysis/gu/ukb'),
            (('ukb', '', '003'), '/mnt/d/analysis/gu/ukb003'),
            (('ukb', '/mnt/d/analysis/gu/ukb///', '003'), '/mnt/d/analysis/gu/ukb003'),
            (('ukb-cohort-test', '/mnt/d/analysis/gu/ukb003', ''), '/mnt/d/analysis/gu/ukb003'),
        ]:
            with self.subTest(args=args):
                result = self.resolve(*args)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout.strip(), expected)

    def test_invalid_paths_and_postfix_do_not_create_directories(self):
        for directory, postfix in [('/tmp/result', ''), ('/', '003'),
                                   ('/mnt/d/analysis/gu/ukb', '/003'),
                                   ('/mnt/d/analysis/gu/ukb', '../other')]:
            with self.subTest(directory=directory, postfix=postfix):
                self.assertNotEqual(self.resolve('ukb', directory, postfix).returncode, 0)

    def test_relocation_preserves_archives_and_rebases_recovered_commands(self):
        with tempfile.TemporaryDirectory() as directory:
            old = Path(directory) / 'results'
            run = old / 'phyml/1kg/chr1'
            run.mkdir(parents=True)
            (run / 'final').mkdir()
            command = run / 'chr1.cmd'
            command.write_text(f'export GU_ANALYSIS_ROOT={old}\nexec gu.sh --loci {run}/final/loci.tsv\n')
            (run / 'final/loci.tsv').write_text('chr\tstart\tend\n1\t0\t100\n')
            common.gu_archive_run(old, run, 'phyml')
            new = old / '1kg'
            new.mkdir()
            (old / 'phyml').rename(new / 'phyml')
            (new / '.gu-path-relocations.json').write_text(json.dumps([str(old)]))
            run = new / 'phyml/1kg/chr1'
            archive = run / 'phyml.raw.tar.gz'
            digest = common.gu_result_digest(archive)
            restored = Path(directory) / 'recovered'
            restored.mkdir()
            common.gu_extract_native(archive, restored, run, new)
            text = (restored / 'chr1.cmd').read_text()
            self.assertIn(f'GU_ANALYSIS_ROOT={new}', text)
            self.assertIn(str(run / 'final/loci.tsv'), text)
            common.gu_rebase_paths(restored, old, new)
            self.assertEqual(text, (restored / 'chr1.cmd').read_text())
            self.assertEqual(digest, common.gu_result_digest(archive))


if __name__ == '__main__':
    unittest.main()
