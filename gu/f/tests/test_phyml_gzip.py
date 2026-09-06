"""Gzip inventories: producers, legacy readers and normalize packaging."""
import argparse
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from comm import read_tsv_rows, resolve_tsv_path, write_tsv_rows
import dual_lead
import normalize_results
import phyml_evidence as evidence
import phyml_region_scan as region


class GzipInventories(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)

    def test_roundtrip_and_legacy(self):
        plain = self.root / 'example.unfiltered.tsv'
        plain.write_text('id\tseq\nlegacy\tACGT\n')
        compressed = plain.with_name(plain.name + '.gz')
        self.assertEqual(evidence.read_tsv(compressed)[0]['id'], 'legacy')
        rows = [{'id': '样本\t"quoted"\nline', 'seq': 'ACGT' * 50000}]
        region.write_tsv(plain, ['id', 'seq'], iter(rows))
        self.assertFalse(plain.exists())
        self.assertTrue(compressed.exists())
        for reader in (region.read_tsv, evidence.read_tsv, dual_lead.rows):
            self.assertEqual(reader(plain), rows)
            self.assertEqual(reader(compressed), rows)
        plain.write_text('id\tseq\nstale\tAAAA\n')
        self.assertEqual(read_tsv_rows(plain), rows)
        self.assertEqual(resolve_tsv_path(plain), compressed)
        write_tsv_rows(self.root / 'empty.unfiltered.tsv', ['id'], [])
        self.assertEqual(read_tsv_rows(self.root / 'empty.unfiltered.tsv.gz'), [])
        write_tsv_rows(self.root / 'ordinary.tsv', ['id'], [{'id': 'plain'}])
        self.assertTrue((self.root / 'ordinary.tsv').is_file())

    def test_failed_write_preserves_previous_outputs(self):
        plain = self.root / 'example.unfiltered.tsv'
        write_tsv_rows(plain, ['id'], [{'id': 'old'}])
        compressed = plain.with_name(plain.name + '.gz')
        previous = compressed.read_bytes()
        plain.write_text('id\nlegacy\n')

        def broken_rows():
            yield {'id': 'new'}
            raise RuntimeError('interrupted producer')

        with self.assertRaisesRegex(RuntimeError, 'interrupted producer'):
            write_tsv_rows(plain, ['id'], broken_rows())
        self.assertEqual(compressed.read_bytes(), previous)
        self.assertEqual(plain.read_text(), 'id\nlegacy\n')
        self.assertEqual(list(self.root.glob('.*.tmp')), [])

    def test_region_production_and_downstream(self):
        out = self.root / 'phyml/1kg/test'
        final, locus = out / 'final', out / 'loci/L1'
        write_tsv_rows(final / 'loci.tsv', ['locus_id', 'chr', 'core_start', 'core_end'],
                       [{'locus_id': 'L1', 'chr': '6', 'core_start': 100, 'core_end': 116}])
        write_tsv_rows(locus / 'sites.tsv', ['chr', 'pos', 'id', 'ref', 'alt'],
                       [{'chr': '6', 'pos': 101 + i, 'id': f'rs{i}', 'ref': 'A', 'alt': 'G'}
                        for i in range(16)])
        write_tsv_rows(locus / 'haplotypes.tsv', ['hap_id', 'n', 'copies', 'seq'],
                       [{'hap_id': 'H1', 'n': 2, 'copies': 'S1:1;S1:2', 'seq': 'ACGT' * 4}])
        write_tsv_rows(locus / 'archaic.tsv', ['archaic', 'seq'],
                       [{'archaic': 'Vindija', 'seq': 'ACGT' * 4}])
        parser = argparse.ArgumentParser()
        region.add_args(parser)
        args = parser.parse_args(['--out', str(out)])
        region.prepare(args)
        region.summarize(args)
        self.assertEqual(len(list(final.glob('*.unfiltered.tsv.gz'))), 8)
        self.assertFalse(list(final.glob('*.unfiltered.tsv')))
        families = evidence.read_tsv(final / 'region_candidate_families.unfiltered.tsv.gz')
        self.assertEqual(len(families), 1)
        self.assertIn('region_common_tree_linked', families[0])
        self.assertIn('L1', evidence.final_rows_by_locus(out, 'region_scan_loci.unfiltered.tsv.gz'))
        auth = {'_evidence_file': str(final / 'evidence_haplotypes.tsv'), 'locus_id': 'L1'}
        _, _, hp, regional = dual_lead.resolve_sequence_files(auth)
        self.assertTrue(regional)
        self.assertEqual(dual_lead.rows(hp)[0]['hap_id'], 'H1')
        old = hp.with_suffix('')
        with gzip.open(hp, 'rb') as source:
            old.write_bytes(source.read())
        hp.unlink()
        self.assertEqual(dual_lead.resolve_sequence_files(auth)[2], old)
        self.assertEqual(evidence.read_tsv(hp)[0]['hap_id'], 'H1')
        outputs = json.loads((final / 'region_scan_parameters.json').read_text())['outputs']
        self.assertTrue(all(name.endswith('.gz') for name in outputs.values() if 'unfiltered' in name))

    def test_normalize_packages_gzip_and_legacy_once(self):
        src = self.root / 'phyml/1kg/test/final'
        src.mkdir(parents=True)
        legacy = src / 'legacy.unfiltered.tsv'
        legacy.write_text('id\nlegacy\n')
        modern = src / 'modern.unfiltered.tsv'
        write_tsv_rows(modern, ['id'], [{'id': 'current'}])
        modern.write_text('id\nstale\n')
        package = self.root / 'normalize'
        package.mkdir()
        for _ in range(2):
            counts = normalize_results.package_core_results(self.root, package)
            self.assertEqual(counts['phyml'], 2)
            dest = package / 'phyml/1kg/test/final'
            self.assertFalse(list(dest.glob('*unfiltered.tsv')))
            self.assertEqual(read_tsv_rows(dest / legacy.name), [{'id': 'legacy'}])
            self.assertEqual(read_tsv_rows(dest / modern.name), [{'id': 'current'}])
        self.assertTrue(legacy.exists())


if __name__ == '__main__':
    unittest.main()
