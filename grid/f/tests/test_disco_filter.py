"""Missing PCA samples must be removed from every aligned Disco input."""
import contextlib
import io
from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline_io import disco_inputs


class DiscoFilterTest(unittest.TestCase):
    def test_filter_and_alignment(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pcs = {f'PC{i}': [1., 2., 3.] for i in range(1, 6)}
            pd.DataFrame({'IID': ['003', '001', 'pca_only'], **pcs}).to_csv(root/'pca.tsv', sep='\t', index=False)
            pops = ['AFR', 'EAS', 'EUR', 'SAS']
            pd.DataFrame({'POP': pops, **{f'PC{i}': [10., 20., 30., 40.] for i in range(1, 6)}}).to_csv(root/'centers.tsv', sep='\t', index=False)
            for k,pop in enumerate(pops):
                z = pd.DataFrame({'eid': ['003', 'missing', '001'], f'CSX_{pop}': [3.+k, 8.+k, 1.+k]})
                z.sample(frac=1, random_state=k).to_csv(root/f'csx.{pop}.tsv.gz', sep='\t', index=False)
            def run():
                disco_inputs(root/'pca.tsv', root/'centers.tsv', root, root/'out', 5)
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf): run()
            self.assertIn('scored=3; missing PCA=1; retained=2', buf.getvalue())
            for name in ['pca']+pops:
                z = pd.read_csv(root/'out'/f'{name}.tsv', sep='\t', dtype={'IID': str})
                self.assertEqual(z.IID.tolist(), ['001', '003'])
                if name in pops:
                    self.assertEqual(z.PRS.tolist(), [1.+pops.index(name), 3.+pops.index(name)])
            # The permanent merged file alone must suffice, including withdrawals.
            merged = pd.DataFrame({'eid':['001','003','-1'], **{f'csx.{p}':[1.,3.,8.] for p in pops}})
            merged.to_csv(root/'csx.pgs.gz', sep='\t', index=False)
            for pop in pops: (root/f'csx.{pop}.tsv.gz').unlink()
            with contextlib.redirect_stdout(io.StringIO()): run()
            self.assertEqual(pd.read_csv(root/'out'/'AFR.tsv',sep='\t',dtype={'IID':str}).IID.tolist(),['001','003'])
            # No-PCA overlap is still a useful failure, not an empty score file.
            pca = pd.read_csv(root/'pca.tsv', sep='\t', dtype={'IID': str})
            pca.IID = ['other1', 'other2', 'other3']
            pca.to_csv(root/'pca.tsv', sep='\t', index=False)
            with self.assertRaisesRegex(ValueError, 'Fewer than two'):
                with contextlib.redirect_stdout(io.StringIO()): run()


if __name__ == '__main__':
    unittest.main()
