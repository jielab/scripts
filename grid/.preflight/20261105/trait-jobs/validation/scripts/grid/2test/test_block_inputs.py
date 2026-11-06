import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'f'))
import block_inputs as blocks


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


class CachePathsTest(unittest.TestCase):
    def test_relocation_reuses_identity_but_content_change_invalidates(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            for directory in ('old', 'new'):
                (root/directory).mkdir()
                for pop in blocks.POPS:
                    (root/directory/f'{pop}.37.bed').write_text('1\t0\t100\n')
            names = ','.join(f'{pop}.37.bed' for pop in blocks.POPS)
            old = blocks.resolve_block_files(root/'old', names)
            new = blocks.resolve_block_files(root/'new', names)
            self.assertEqual(old['signature'], new['signature'])
            self.assertEqual(old['layout_id'], new['layout_id'])
            self.assertNotEqual(old['files'], new['files'])
            (root/'new'/'EUR.37.bed').write_text('1\t0\t101\n')
            self.assertNotEqual(old['signature'], blocks.resolve_block_files(root/'new', names)['signature'])

    def test_noncanonical_boundary_names_are_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            for pop in blocks.POPS:
                (root / f'{pop}.37.bed').write_text('1\t0\t100\n')
            with self.assertRaisesRegex(ValueError, 'supported GRCh37'):
                blocks.resolve_block_files(root, 'AFR.37.bed,EAS.b37.bed,EUR.37.bed,SAS.37.bed')


if __name__ == '__main__':
    unittest.main()
