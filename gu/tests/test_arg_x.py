import importlib.util
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'f'))
import arg_x
spec = importlib.util.spec_from_file_location('needle_x', Path(__file__).resolve().parents[1] / 'f/arg.needle_x.py')
needle = importlib.util.module_from_spec(spec)
spec.loader.exec_module(needle)


class XInputs(unittest.TestCase):
    def test_par_boundaries(self):
        for build, start, end in [('37', 2699521, 154931043), ('38', 2781480, 155701382)]:
            self.assertTrue(arg_x.nonpar(start, build))
            self.assertTrue(arg_x.nonpar(end, build))
            self.assertFalse(arg_x.nonpar(start-1, build))
            self.assertFalse(arg_x.nonpar(end+1, build))

    def test_gt(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            (p/'sex').write_text('#IID SEX\nm 1\nf 2\n')
            header = ('##fileformat=VCFv4.2\n##contig=<ID=X,length=155270560>\n'
                      '##INFO=<ID=AA,Number=1,Type=String,Description="ancestral">\n'
                      '##FORMAT=<ID=GT,Number=1,Type=String,Description="genotype">\n'
                      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tm\tf\n')
            for male, female, ok in [('0','0|1',True), ('0|0','0|1',False), ('0','0/1',False), ('0','1',False)]:
                (p/'x.vcf').write_text(header + f'X\t3000000\ta\tA\tC\t.\tPASS\tAA=A\tGT\t{male}\t{female}\n')
                result = subprocess.run([sys.executable, str(Path(__file__).resolve().parents[1]/'f/arg.validate_vcf.py'), str(p/'x.vcf'), '--chr','X','--x-psam',str(p/'sex')], capture_output=True, text=True)
                self.assertEqual(result.returncode == 0, ok, result.stderr)

    def test_pack_preserves_real_copies(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'chrX'
            psam = Path(d) / 'source.psam'
            psam.write_text('#IID SEX\nm1 1\nf1 2\nm2 1\n')
            keep = Path(d) / 'keep'
            keep.write_text('#IID\nm1\nf1\nm2\n')
            Path(str(p)+'.sample').write_text('ID_1 ID_2 missing\n0 0 0\nm1 m1 0\nf1 f1 0\nm2 m2 0\n')
            Path(str(p)+'.haps').write_text('X a 3000000 A C 0 - 1 0 1 -\nX b 3000100 A C 1 - 0 1 0 -\n')
            needle.pack(p, psam, '37', keep)
            self.assertEqual(Path(str(p)+'.haps').read_text().splitlines()[0], 'X a 3000000 A C 0 1 0 1')
            self.assertEqual(Path(str(p)+'.haplotypes.tsv').read_text().splitlines(), ['eid\thap','m1\t0','f1\t0','f1\t1','m2\t0'])
            self.assertEqual(len(Path(str(p)+'.sample').read_text().splitlines()), 4)

    def test_odd_panel_requires_explicit_policy(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            (p/'s').write_text('#IID SEX\nm 1\nf 2\n')
            (p/'k').write_text('#IID\nm\nf\n')
            with self.assertRaises(ValueError):
                needle.panel(p/'k', p/'s', 'error', p/'q')
            self.assertIn('m', (p/'k').read_text())
            needle.panel(p/'k', p/'s', 'drop-last', p/'q')
            self.assertEqual((p/'k').read_text(), '#IID\nf\n')

    def test_zarr_removes_padding_and_rejects_changing_ploidy(self):
        import numpy as np
        import zarr
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            src = zarr.open_group(str(p/'src'), mode='w')
            src.create_dataset('sample_id', data=np.array(['male','female']))
            src.create_dataset('call_genotype', data=np.array([[[0,-2],[1,0]], [[-1,-2],[0,1]]], dtype='i1'), chunks=(1,2,2))
            identities = arg_x.haploid_zarr(p/'src', p/'dest')
            self.assertEqual(identities, [('male',1),('female',1),('female',2)])
            np.testing.assert_array_equal(zarr.open(str(p/'dest'))['call_genotype'][:,:,0], [[0,1,0],[-1,0,1]])
            src['call_genotype'][1,0,1] = 0
            with self.assertRaises(ValueError):
                arg_x.haploid_zarr(p/'src', p/'bad')

    def test_tree_identity(self):
        import tskit
        tables = tskit.TableCollection(100)
        for _ in range(3):
            tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
        root = tables.nodes.add_row(time=1)
        for node in range(3):
            tables.edges.add_row(0,100,root,node)
        tables.sort()
        ts = arg_x.restore_individuals(tables.tree_sequence(), [('m',1),('f',1),('f',2)])
        self.assertEqual(ts.num_individuals, 2)
        self.assertEqual(list(ts.individual(0).nodes), [0])
        self.assertEqual(list(ts.individual(1).nodes), [1,2])
        self.assertEqual(ts.individual(0).metadata['sample'], 'm')


if __name__ == '__main__':
    unittest.main()
