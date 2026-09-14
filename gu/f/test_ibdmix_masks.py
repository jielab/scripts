"""Local mask preflight must fail before any network access or computation."""
import gzip
import hashlib
import urllib.request
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import ibdmix_workflow as ibd


class LocalMaskTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)/'mask-override'
        self.raw = self.root
        self.raw.mkdir(parents=True)
        self.archaic = Path(self.tmp.name)/'archaic'/'vcf'
        self.mask = self.archaic.parent/'mask'
        network = patch.object(urllib.request, 'urlopen', side_effect=AssertionError('Network forbidden'))
        self.network = network.start()
        self.addCleanup(network.stop)

    def write(self, path, data=b'1\t0\t100\n'):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(gzip.compress(data) if path.suffix in {'.gz', '.bgz'} else data)
        return path

    def shared(self):
        self.write(self.raw/'common/1kg.strict_mask.autosomes.bed')
        self.write(self.raw/'common/genomicSuperDups.txt.gz', b'0\tchr1\t10\t20\n')

    def test_all_missing_reported_without_download_or_directory_creation(self):
        root = Path(self.tmp.name)/'absent'
        with self.assertRaises(FileNotFoundError) as error:
            ibd.check_mask_inputs(root, ['1', '2'], ['Chagyr', 'Denisova25'], self.archaic)
        message = str(error.exception)
        for name in ['1kg.strict_mask.autosomes.bed', 'genomicSuperDups.txt.gz',
                     'Chagyr/chr1_mask.bed.gz', 'Chagyr/chr2_mask.bed.gz',
                     'Denisova25.chr1.hg19', 'Denisova25.chr2.hg19']:
            self.assertIn(name, message)
        self.assertFalse(root.exists())
        self.network.assert_not_called()

    def test_local_and_completed_cache_are_both_reused(self):
        self.shared()
        local = self.write(self.mask/'Chagyr/chr1_mask.bed.gz')
        cached = self.write(self.raw/'Chagyr/chr2_mask.bed.gz')
        self.assertEqual(ibd.check_mask_inputs(self.root, ['1', '2'], ['Chagyrskaya'], self.archaic), 4)
        self.assertEqual(ibd.published_reference_mask('Chagyr', '1', self.raw, self.archaic), local)
        self.assertEqual(ibd.published_reference_mask('Chagyr', '2', self.raw, self.archaic), cached)
        self.network.assert_not_called()

    def test_aria2_preallocation_is_not_a_complete_file(self):
        self.shared()
        local = self.write(self.mask/'Chagyr/chr1_mask.bed.gz')
        Path(str(local)+'.aria2').touch()
        with self.assertRaises(FileNotFoundError):
            ibd.check_mask_inputs(self.root, ['1'], ['Chagyr'], self.archaic)
        with self.assertRaises(FileNotFoundError):
            ibd.published_reference_mask('Chagyr', '1', self.raw, self.archaic)
        # A separate, completed cache remains usable while aria2 is downloading.
        cached = self.write(self.raw/'Chagyr/chr1_mask.bed.gz')
        self.assertEqual(ibd.published_reference_mask('Chagyr', '1', self.raw, self.archaic), cached)

    def test_empty_and_part_files_do_not_satisfy_mask(self):
        path = self.raw/'Chagyr/chr1_mask.bed.gz'
        self.write(Path(str(path)+'.part.123'))
        path.touch()
        with self.assertRaises(FileNotFoundError):
            ibd.published_reference_mask('Chagyr', '1', self.raw, self.archaic)

    def test_selected_refs_only_and_x_uses_correct_mask(self):
        self.write(self.raw/'common/genomicSuperDups.txt.gz')
        x = self.write(self.mask/'Denisova25/Denisova25.chrX.hg19.L35MQ25.map35_100.GCcov.noSimpleRepeat.noIndel.bed.gz')
        self.assertEqual(ibd.check_mask_inputs(self.root, ['23'], ['Den25'], self.archaic), 2)
        self.assertEqual(ibd.published_reference_mask('Denisova25', 'X', self.raw, self.archaic), x)
        with self.assertRaises(FileNotFoundError):
            ibd.check_mask_inputs(self.root, ['1'], ['Denisova25'], self.archaic)

    def test_custom_empty_exclusion_bed_is_valid_but_absent_is_not(self):
        custom = Path(self.tmp.name)/'custom'
        self.write(custom/'Altai/chr1.bed', b'')
        self.assertEqual(ibd.check_mask_inputs(self.root, ['1'], ['Altai'], custom_masks=custom), 1)
        with self.assertRaises(FileNotFoundError):
            ibd.check_mask_inputs(self.root, ['1', '2'], ['Altai'], custom_masks=custom)

    def test_prepare_fails_before_axt_or_modern_vcf_work(self):
        with patch.object(ibd, 'local_axt', side_effect=AssertionError('Too late')):
            with self.assertRaises(FileNotFoundError):
                ibd.prepare_masks(self.root, '1', 'absent.vcf', 'absent.fa', 'absent-upstream',
                                  Path(self.tmp.name)/'out', refs=['Chagyr'], archaic_root=self.archaic)
        self.network.assert_not_called()

    def test_unified_root_needs_no_axt_cache_and_preserves_minimal_masks(self):
        self.shared()
        minimal = self.write(self.raw/'Altai/AltaiNea.map35_50.MQ30.Cov.indels.TRF.bed.bgz')
        self.write(self.raw/'Chagyr/chr1_mask.bed.gz')
        self.assertEqual(ibd.check_mask_inputs(self.root, ['1'], ['Altai', 'Chagyr']), 4)
        self.assertEqual(ibd.reference_mask_resource(ibd.MINIMAL_MASK_NAMES['Altai'], self.root), minimal)
        self.assertFalse((self.root/'raw').exists())
        self.assertFalse((self.root/'ibdmix-cache').exists())

    def test_mask_root_override_has_precedence(self):
        self.write(self.mask/'Chagyr/chr1_mask.bed.gz')
        override = self.write(self.root/'Chagyr/chr1_mask.bed.gz')
        self.assertEqual(ibd.published_reference_mask('Chagyr', '1', self.root, self.archaic), override)

    def test_missing_axt_or_checksum_never_downloads(self):
        axt_root = Path(self.tmp.name)/'axt'
        with self.assertRaises(FileNotFoundError):
            ibd.local_axt('1', 'panTro2', axt_root)
        self.assertFalse(axt_root.exists())
        self.write(axt_root/'chr1.hg19.panTro2.synNet.axt.gz')
        with self.assertRaises(FileNotFoundError):
            ibd.local_axt('1', 'panTro2', axt_root)
        self.network.assert_not_called()

    def test_local_axt_checks_md5_and_rejects_aria2(self):
        axt_root = Path(self.tmp.name)/'axt'
        axt = self.write(axt_root/'chr1.hg19.panTro2.synNet.axt.gz')
        digest = hashlib.md5(axt.read_bytes()).hexdigest()
        md5 = axt_root/'md5sum.txt'
        md5.write_text(f'{digest}  {axt.name}\n')
        self.assertEqual(ibd.local_axt('1', 'panTro2', axt_root), axt)
        partial = Path(str(axt)+'.aria2')
        partial.touch()
        with self.assertRaises(FileNotFoundError):
            ibd.local_axt('1', 'panTro2', axt_root)
        partial.unlink()
        axt.write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, 'Publisher checksum mismatch'):
            ibd.local_axt('1', 'panTro2', axt_root)
        self.network.assert_not_called()

    def test_local_axt_accepts_existing_alternate_checksum_directory(self):
        axt_root = Path(self.tmp.name)/'axt'
        axt = self.write(axt_root/'chr1.hg19.panTro2.synNet.axt.gz')
        md5 = axt_root/'checksums/panTro2.md5sum.txt'
        md5.parent.mkdir()
        md5.write_text(f'{hashlib.md5(axt.read_bytes()).hexdigest()}  {axt.name}\n')
        self.assertEqual(ibd.local_axt('1', 'panTro2', axt_root), axt)


if __name__ == '__main__':
    unittest.main()
