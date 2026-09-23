"""Regression checks for independent reference burdens and missing-run handling."""
import csv
import json
import sqlite3
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import introgression_density as density


class ReferenceSummaryTests(unittest.TestCase):
    def setUp(self):
        self.con = sqlite3.connect(':memory:')
        self.con.execute('CREATE TABLE segments (dataset_id, genome_build, method, source, source_class, sample_id, chr, start, end)')
        self.con.executemany('INSERT INTO segments VALUES (?,?,?,?,?,?,?,?,?)', [
            ('test', 'GRCh37', 'ibdmix', ref, lineage, sample, chrom, left, right)
            for ref, lineage, sample, chrom, left, right in [
                ('Altai', 'Neanderthal', 'a', '1', 0, 20),
                ('Altai', 'Neanderthal', 'a', '1', 10, 30),
                ('Altai', 'Neanderthal', 'a', '2', 0, 50),
                ('Altai', 'Neanderthal', 'a', 'X', 0, 50),
                ('Chagyr', 'Neanderthal', 'a', '1', 20, 45),
                ('Chagyrskaya', 'Neanderthal', 'a', '1', 40, 60),
                ('Vindija', 'Neanderthal', 'a', '1', 0, 90),
                ('Denisova', 'Denisovan', 'a', '1', 0, 10),
                ('Denisova25', 'Denisovan', 'a', '1', 0, 15),
            ]])
        self.con.commit()
        self.lengths = {'1': 100, '2': 100, 'X': 100}
        self.targets = {'1': {'a', 'b'}}

    def tearDown(self):
        self.con.close()

    def summary(self, ref):
        return density.reference_summary(self.con, 'test', 'GRCh37', ['a', 'b', 'c'], self.lengths, self.targets, ref)

    def test_independent_unions_autosomes_and_alias(self):
        for ref, expected in zip(density.SUMMARY_REFERENCES, [30, 40, 90, 10, 15]):
            with self.subTest(reference=ref):
                a, b, c = self.summary(ref)
                self.assertEqual(a['archaic_bp'], expected)
                self.assertEqual(a['tested_bp'], 200)
                self.assertEqual(a['n_chromosomes'], 1)
                self.assertEqual(a['coverage_pct'], expected / 2)
                self.assertEqual(b['archaic_bp'], 0)  # Tested, no positive calls.
                self.assertEqual(c['archaic_bp'], '')  # No certified run.

    def test_cache_certifies_references_and_excludes_background_only(self):
        for background_only in [False, True]:
            with self.subTest(background_only=background_only), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                meta = root / 'run.meta.tsv'
                meta.write_text('modern_vcf\tfixture.vcf\nrefs\tAltai Denisova\nbackground_filter\t1\nexport_denisovan\t' + ('0' if background_only else '1') + '\n')
                (root / 'samples').mkdir()
                (root / 'samples' / 'ALL.txt').write_text('a\nb\n')
                database = root / 'test.sqlite'
                dest = sqlite3.connect(database)
                self.con.backup(dest)
                dest.execute('CREATE TABLE method_runs (dataset_id, genome_build, method, status, chr, evidence_eligible, raw_file, availability_note)')
                dest.execute('INSERT INTO method_runs VALUES (?,?,?,?,?,?,?,?)', ('test', 'GRCh37', 'ibdmix', 'complete', '1', 1, str(meta), 'audited IBDmix profile=test'))
                dest.execute('CREATE TABLE sample_populations (dataset_id, sample_id, population, super_population)')
                dest.executemany('INSERT INTO sample_populations VALUES (?,?,?,?)', [('test', s, 'YRI', 'AFR') for s in ['a', 'b']])
                dest.commit()
                dest.close()
                output = root / 'cache'
                with patch.dict(density.CHROM_LENGTHS, {'37': self.lengths}):
                    density.prepare(database, output, bin_bp=100)
                with (output / 'test.GRCh37' / 'current.tsv').open() as stream:
                    pointer = next(csv.DictReader(stream, delimiter='\t'))['directory']
                cache = output / 'test.GRCh37' / pointer
                with (cache / 'archaic_summary.tsv').open() as stream:
                    rows = {(r['reference'], r['sample_id']): r for r in csv.DictReader(stream, delimiter='\t')}
                self.assertEqual(rows['Altai', 'a']['archaic_bp'], '30')
                self.assertEqual(rows['Altai', 'b']['archaic_bp'], '0')
                self.assertEqual(rows['Denisova', 'a']['archaic_bp'], '' if background_only else '10')
                self.assertEqual(rows['Vindija', 'a']['archaic_bp'], '')
                self.assertEqual(rows['Denisova25', 'a']['archaic_bp'], '')
                with (cache / 'neanderthal_summary.tsv').open() as stream:
                    self.assertEqual(next(csv.DictReader(stream, delimiter='\t'))['neanderthal_bp'], '30')
                self.assertEqual(json.loads((cache / 'manifest.json').read_text())['schema'], 10)


if __name__ == '__main__':
    unittest.main()
