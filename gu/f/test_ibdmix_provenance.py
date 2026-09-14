"""Resume survives resource moves but still rejects incompatible results."""
import unittest
import os
from pathlib import Path
import subprocess
import tempfile
from ibdmix_provenance import differences

RECORD = '''pipeline_version\t2026-09-14.1
profile\tmulti_reference
refs\tAltai Denisova
lod_cut\t4
mask_root\t/mnt/i/refGen/archaic/37/mask
software\\t/mnt/d/scripts/gu/f/ibdmix.sh:23367:100
software\\t/mnt/d/scripts/gu/f/ibdmix_workflow.py:30024:100
software\\t/mnt/d/software/gu/ibdmix/build/src/ibdmix:419088:100
sample_file\\t/mnt/i/refGen/1kg/37/samples.txt:57720:100
modern_source\tchr8\tformat\tvcf
modern_source\tchr8\tvcf\\t/mnt/i/refGen/1kg/37/vcf/chr8.vcf.gz:861233457:100
mask_manifest_sha256\tC8\told-metadata
excluded_mask_sha256\tC8\tAltai\t'''+ 'a'*64 + '\nexcluded_mask_sha256\tC8\tDenisova\t' + 'b'*64 + '\n'


class ProvenanceTests(unittest.TestCase):
    def test_recorded_move_and_nonsemantic_script_edits_are_compatible(self):
        moved = RECORD.replace('/mnt/i/refGen/', '/mnt/e/refGen/')
        moved = moved.replace('23367:100', '24000:200').replace('30024:100', '31000:200')
        moved = moved.replace('old-metadata', 'new-metadata')
        self.assertEqual(differences(RECORD, moved), '')

    def test_effective_changes_remain_incompatible(self):
        changes = [
            ('2026-09-14.1', '2026-09-15.1'),
            ('lod_cut\t4', 'lod_cut\t5'),
            ('861233457:100', '861233458:100'),
            ('861233457:100', '861233457:101'),
            ('419088:100', '419088:101'),
            ('57720:100', '57721:100'),
            ('format\tvcf', 'format\tpfile'),
            ('a'*64, 'c'*64),
            ('/mnt/i/refGen/1kg/', '/mnt/e/refGen/another-cohort/'),
        ]
        for old, new in changes:
            with self.subTest(change=(old, new)):
                self.assertTrue(differences(RECORD, RECORD.replace(old, new)))

    def test_missing_checksums_cannot_hide_manifest_change(self):
        incomplete = '\n'.join(line for line in RECORD.splitlines()
                               if not line.startswith('excluded_mask_sha256\tC8\tAltai'))
        with self.assertRaises(ValueError):
            differences(incomplete, incomplete.replace('old-metadata', 'new-metadata'))

    def test_literal_and_real_stat_tabs_compare_equally(self):
        self.assertEqual(differences(RECORD, RECORD.replace('\\t', '\t')), '')

    def test_shell_resume_keeps_results_and_rejects_changed_inputs(self):
        source = Path(__file__).with_name('ibdmix.sh').read_text()
        function = source.split('check_run_provenance(){', 1)[1].split('completed_output_ok(){', 1)[0]
        function = 'check_run_provenance(){' + function
        for compatible in (True, False):
            with self.subTest(compatible=compatible), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                for name in ('tmp', 'genotype', 'raw', 'segments', 'final'):
                    (root/name).mkdir()
                (root/'run.meta.tsv').write_text(RECORD)
                candidate = RECORD.replace('/mnt/i/refGen/', '/mnt/e/refGen/')
                if not compatible:
                    candidate = candidate.replace('lod_cut\t4', 'lod_cut\t5')
                (root/'candidate.tsv').write_text(candidate)
                (root/'raw/keep.txt').write_text('existing result')
                env = dict(os.environ, dirout=tmp, IBDMIX_REPLACE='0',
                           F=str(Path(__file__).parent.resolve()))
                result = subprocess.run(['bash', '-c', 'set -e\n' + function +
                    '\nrun_record(){ cat "$dirout/candidate.tsv"; }\ncheck_run_provenance'],
                    env=env, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0 if compatible else 1, result.stderr)
                self.assertEqual((root/'raw/keep.txt').read_text(), 'existing result')
                self.assertEqual((root/'run.meta.tsv').read_text(), candidate if compatible else RECORD)


if __name__ == '__main__':
    unittest.main()
