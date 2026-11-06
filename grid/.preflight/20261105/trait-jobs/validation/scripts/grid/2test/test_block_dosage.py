"""Use real tiny PGENs to check supplementation across different target panels."""
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
import pgenlib

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'f'))
import block_dosage as dosage
import block_inputs as blocks


def write_pgen(root, ids, variants, values):
    root.mkdir()
    (root / 'chr1.psam').write_text('#IID\n' + '\n'.join(ids) + '\n')
    with (root / 'chr1.pvar').open('w') as stream:
        stream.write('#CHROM\tPOS\tID\tREF\tALT\n')
        for snp, bp, ref, alt in variants:
            stream.write(f'1\t{bp}\t{snp}\t{ref}\t{alt}\n')
    with pgenlib.PgenWriter(os.fsencode(root / 'chr1.pgen'), len(ids),
                           variant_ct=len(variants), nonref_flags=False, dosage_present=True) as writer:
        for row in values:
            writer.append_dosages(np.array(row, dtype=np.float32))


class SupplementTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        write_pgen(self.root / 'hap', ['01', '02', '03'], [('rs2', 20, 'C', 'T')],
                   [[.25, 1.5, -9]])
        # Both variant order and sample order differ from the frozen cache.
        write_pgen(self.root / 'typ', ['02', 'unused', '03', '01'],
                   [('rs3', 30, 'A', 'G'), ('unused', 10, 'A', 'C'), ('rs2', 20, 'C', 'T')],
                   [[0, 2, 2, 1], [0, 1, 2, 0], [1.5, 0, -9, .25]])
        self.config = dict(cache_root=str(self.root / 'cache'), dosage_supplement_dir=str(self.root / 'typ'))
        self.prep = dict(options=dict(dir_gen=str(self.root / 'hap')))
        self.cohort = pd.DataFrame(dict(eid=['03', '01', '02']))
        self.variants = pd.DataFrame(dict(CHR=[1], BP=[5], SNP=['rs1'], REF=['A'], ALT=['G']))
        self.original = np.array([[.25], [1], [2]], dtype=np.float16)
        self.missing = pd.DataFrame(dict(CHR=[1, 1], BP=[20, 30], SNP=['rs2', 'rs3'],
                                         A1=['C', 'G'], A2=['T', 'A'], BETA=[.5, 2.]))

    def run_supplement(self):
        return dosage.supplement(self.config, self.prep, self.cohort, self.variants,
                                 self.original, self.missing, blocks.sha256)

    def test_configured_panel_aligns_real_dosages_and_reuses_cache(self):
        variants, combined, receipt = self.run_supplement()
        expected = np.array([[.25, -1, 2], [1, .25, 1], [2, 1.5, 0]], dtype=np.float32)
        self.assertEqual(variants.SNP.tolist(), ['rs1', 'rs2', 'rs3'])
        np.testing.assert_array_equal(combined[:], expected)
        np.testing.assert_array_equal(combined[np.ix_([2, 0], [2, 0, 1])], expected[np.ix_([2, 0], [2, 0, 1])])
        np.testing.assert_array_equal(self.original[:, 0], [.25, 1, 2])
        self.assertEqual(receipt['identity']['source_directory'], str(self.root / 'typ'))
        self.assertTrue(all(Path(s['path']).parent == self.root / 'typ' for s in receipt['identity']['sources']))
        bed = pd.DataFrame(dict(CHR=[1], start=[0], end=[100], block_local=[0]))
        np.testing.assert_allclose(blocks.direct_assigned_score(combined, [0, 1, 2], self.missing, bed, variants),
                                   [4, 2.875, .25])  # REF effect, ALT effect, and missing=-1.
        with patch.object(pgenlib, 'PgenReader', side_effect=AssertionError('unexpected rebuild')):
            _, reused, saved = self.run_supplement()
        self.assertEqual(saved, receipt)
        np.testing.assert_array_equal(reused[:], expected)

    def test_legacy_source_is_used_when_override_is_absent(self):
        self.config.pop('dosage_supplement_dir')
        self.missing = self.missing.iloc[:1]
        _, combined, receipt = self.run_supplement()
        np.testing.assert_array_equal(combined[:][:, [1]], [[-1], [.25], [1.5]])
        self.assertEqual(receipt['identity']['source_directory'], str(self.root / 'hap'))

    def test_released_supplement_rebuilds_from_pgen(self):
        _, combined, receipt = self.run_supplement()
        expected = combined[:].copy()
        combined.supplement._mmap.close()
        folder = self.root / 'cache/dosage_supplements' / receipt['signature']
        (folder / 'genotypes.released.json').write_text('{}')
        (folder / 'genotypes.npy').unlink()
        _, rebuilt, saved = self.run_supplement()
        np.testing.assert_array_equal(rebuilt[:], expected)
        self.assertEqual(saved, receipt)
        self.assertFalse((folder / 'genotypes.released.json').exists())

    def test_wrong_panel_fails_with_actionable_source_before_cache_write(self):
        self.config.pop('dosage_supplement_dir')
        with self.assertRaisesRegex(ValueError, r'1 weighted SNPs unavailable.*hap.*rs3.*dosage_supplement_dir'):
            self.run_supplement()
        self.assertFalse((self.root / 'cache').exists())

    def test_source_directory_changes_cache_identity(self):
        _, _, first = self.run_supplement()
        shutil.copytree(self.root / 'typ', self.root / 'other')
        self.config['dosage_supplement_dir'] = str(self.root / 'other')
        _, _, second = self.run_supplement()
        self.assertNotEqual(first['signature'], second['signature'])

    def test_coordinate_and_allele_errors_are_still_rejected(self):
        for column, value, error in [('BP', 21, 'CHR:BP mismatches'), ('A1', 'A', 'allele mismatch'),
                                      ('A2', 'G', 'allele mismatch')]:
            with self.subTest(column=column):
                previous = self.missing.loc[0, column]
                self.missing.loc[0, column] = value
                with self.assertRaisesRegex(ValueError, error):
                    self.run_supplement()
                self.missing.loc[0, column] = previous
        self.assertFalse((self.root / 'cache').exists())

    def test_missing_cohort_id_does_not_publish_a_supplement(self):
        self.cohort.loc[0, 'eid'] = 'absent'
        with self.assertRaisesRegex(ValueError, 'cohort IDs missing'):
            self.run_supplement()
        self.assertEqual(list((self.root / 'cache').rglob('COMPLETE.json')), [])


if __name__ == '__main__':
    unittest.main()
