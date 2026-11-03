"""Discover finalized GU publications in the genotype-adjacent HPC layout."""
import gzip
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'f'))
from herv_archaic import discover

class GenLayoutTests(unittest.TestCase):
    def publish(self, root):
        d = root / 'ibdmix/chrX'
        (d / 'samples').mkdir(parents=True)
        (d / 'final').mkdir()
        (d / 'samples/ALL.txt').write_text('0001\n')
        (d / '.complete').touch()
        (d / 'run.meta.tsv').write_text('genome_build\tb37\ncoordinate_system\t0-based-half-open\nrefs\tAltai\nlod_cut\t4\nlen_cut\t50000\npipeline_version\ttest\n')
        with gzip.open(d / 'final/all_archaic_refs.segments.tsv.gz', 'wt') as f:
            f.write('ID\tchrom\tstart\tend\tanc\tslod\n')

    def test_gen_build_and_chunk_keep_stable_dataset_names(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            self.publish(root / '1kg/37')
            self.publish(root / 'ukb/37/cohort013')
            records = discover(root, '1kg,ukb013', chromosomes='X', references=('Altai',))
            self.assertEqual([r['dataset'] for r in records], ['1kg.chrX', 'ukb013.chrX'])
            self.assertEqual([r['n_samples'] for r in records], [1, 1])
            self.assertEqual(records[1]['source'], str(root / 'ukb/37/cohort013/ibdmix'))

    def test_legacy_explicit_root_still_works(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            self.publish(root / '1kg')
            record = discover(root, '1kg', chromosomes='X', references=('Altai',))[0]
            self.assertEqual(record['source'], str(root / '1kg/ibdmix'))

    def test_missing_result_names_expected_gen_build(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            (root / '1kg/37').mkdir(parents=True)
            with self.assertRaisesRegex(ValueError, str(root / '1kg/37/ibdmix')):
                discover(root, '1kg', chromosomes='X', references=('Altai',))
