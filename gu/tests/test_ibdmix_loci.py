"""Regression tests for cohort window support, not raw-call multiplicity."""
import json
from pathlib import Path
import sqlite3
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'f'))
import ibdmix_loci as loci


class WindowSupportTests(unittest.TestCase):
    def profile(self, raw, chrom='1', lineage='Neanderthal', n=1000, width=50000):
        raw = np.array(raw, dtype=float)
        starts = np.arange(len(raw)) * width
        return dict(source_class=lineage, chr=chrom, n_tested=n, starts=starts,
                    ends=starts+width, raw=raw, covered_bp=raw*2*n*width)

    def test_ten_carriers_in_1000_people_do_not_qualify(self):
        members = {str(i) for i in range(1000)}
        targets = {('Neanderthal', '1'): {'Altai': members}}
        rows = [('Neanderthal', '1', str(i), 0, 150000, 'Altai', 10) for i in range(10)]
        p = loci.window_profiles(sorted(rows), targets, {'1': 1000000}, 50000)
        np.testing.assert_allclose(p[0]['raw'][:3], .005)
        called, _ = loci.call_loci(p, 'test', 'GRCh37')
        self.assertEqual(called, [])
        self.assertEqual(p[0]['n_tested'], 1000)  # Includes all 990 zero-call people.

    def test_union_deduplicates_reference_hits_and_partial_windows(self):
        targets = {('Neanderthal', '1'): {'Altai': {'A','B'}, 'Vindija': {'A','B'}}}
        rows = [('Neanderthal','1','A',0,75000,'Altai',10),
                ('Neanderthal','1','A',0,75000,'Altai',10),
                ('Neanderthal','1','A',25000,100000,'Vindija',10),
                ('Neanderthal','1','B',75000,100000,'Altai',10),
                ('Neanderthal','1','OUTSIDE_ROSTER',0,100000,'Altai',10)]
        p = loci.window_profiles(rows, targets, {'1': 100000}, 50000)[0]
        np.testing.assert_equal(p['covered_bp'], [50000, 75000])
        np.testing.assert_allclose(p['raw'], [.25, .375])

    def test_integral_matches_per_base_brute_force_including_terminal_bin(self):
        rng = np.random.default_rng(42)
        length, width = 137, 17
        members = {str(i) for i in range(8)}
        targets = {('Neanderthal', '1'): {'Altai': members, 'Vindija': members}}
        expected = np.zeros((8, length), dtype=bool)
        rows = []
        for i in range(8):
            for _ in range(13):
                lo, hi = sorted(rng.integers(0, length+1, size=2))
                rows.append(('Neanderthal','1',str(i),int(lo),int(hi),'Altai',5))
                rows.append(('Neanderthal','1',str(i),int(lo),int(hi),'Vindija',5))
                expected[i,lo:hi] = True
        p = loci.window_profiles(sorted(rows), targets, {'1': length}, width)[0]
        exact = [expected[:,x:x+width].sum() for x in range(0,length,width)]
        np.testing.assert_equal(p['covered_bp'], exact)

    def test_isolated_peak_removed_and_zero_gap_splits_regions(self):
        p = self.profile([0,.06,0,0,.05,.05,.05,0,.05,.05,.05]+[0]*89)
        called, _ = loci.call_loci([p], 'test', 'GRCh37')
        self.assertEqual([(r['start'],r['end']) for r in called], [(200000,350000),(400000,550000)])
        self.assertTrue(all(r['core_windows'] == 3 for r in called))

    def test_raw_gate_prevents_smoothing_from_filling_a_valley(self):
        p = self.profile([0,.05,.05,0,.05,.05,0]+[0]*93)
        called, _ = loci.call_loci([p], 'test', 'GRCh37')
        self.assertEqual(len(called), 2)
        self.assertEqual(called[0]['end'], 150000)
        self.assertEqual(called[1]['start'], 200000)

    def test_genome_background_adapts_and_flat_background_is_not_a_peak(self):
        p = self.profile([.02]*100)
        called, bg = loci.call_loci([p], 'test', 'GRCh37')
        self.assertAlmostEqual(bg['Neanderthal'], .02)
        self.assertEqual(called, [])
        p['raw'][20:23] = .07
        p['covered_bp'] = p['raw'] * 2 * p['n_tested'] * 50000
        called, bg = loci.call_loci([p], 'test', 'GRCh37')
        self.assertEqual(len(called), 1)
        self.assertAlmostEqual(called[0]['core_threshold'], 2*bg['Neanderthal'])

    def test_chromosomes_and_lineages_never_merge_and_no_count_cap(self):
        raw = ([0]*10 + [.05]*2 + [0]*10)*130
        profiles = [self.profile(raw), self.profile([0]*10+[.05]*2+[0]*10, '2'),
                    self.profile([0]*10+[.05]*2+[0]*10, '1', 'Denisovan')]
        called, _ = loci.call_loci(profiles, 'test', 'GRCh37')
        self.assertEqual(len(called), 132)
        self.assertEqual(len({r['segment_code'] for r in called}), 132)
        self.assertEqual([r['rank'] for r in called], list(range(1,133)))

    def test_tested_zero_call_chromosome_contributes_background(self):
        targets = {('Neanderthal', c): {'Altai': {'A'}} for c in ('1','2')}
        rows = [('Neanderthal','1','A',0,100000,'Altai',10)]
        p = loci.window_profiles(rows, targets, {'1':100000,'2':100000}, 50000)
        _, bg = loci.call_loci(p, 'test', 'GRCh37')
        self.assertAlmostEqual(bg['Neanderthal'], .25)
        np.testing.assert_equal(p[1]['raw'], [0,0])

    def test_annotation_boundaries_are_half_open_and_reference_hits_deduplicated(self):
        p = self.profile([0]*10 + [.05]*2 + [0]*10)
        called, _ = loci.call_loci([p], 'test', 'GRCh37')
        targets = {('Neanderthal','1'): {'Altai': {'A'}, 'Vindija': {'A'}}}
        rows = [('Neanderthal','1','A',450000,500000,'Altai',100),
                ('Neanderthal','1','A',500000,600000,'Altai',5),
                ('Neanderthal','1','A',500000,600000,'Vindija',6),
                ('Neanderthal','1','A',600000,650000,'Altai',100)]
        loci.annotate_loci(rows, called, targets)
        self.assertEqual(called[0]['n_carriers'],1)
        self.assertEqual(called[0]['n_references'],2)
        self.assertEqual(called[0]['n_calls'],2)
        self.assertEqual(called[0]['max_lod'],6)

    def test_roster_provenance_cache_and_control_exclusion(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            database = root/'gu.sqlite'
            run = root/'ibdmix'/'chr1'
            roster = run/'samples'/'C1'/'ALL.txt'
            roster.parent.mkdir(parents=True)
            roster.write_text('\n'.join(str(i) for i in range(1000)))
            metadata = run/'run.meta.tsv'
            metadata.write_text('modern_vcf\tinput.vcf\nrefs\tAltai Altai.2013\n')
            with sqlite3.connect(database) as con:
                con.executescript('''CREATE TABLE method_runs(dataset_id,genome_build,chr,method,status,evidence_eligible,raw_file);
                CREATE TABLE segments_evidence(dataset_id,genome_build,method,reference_role,source_class,chr,sample_id,start,end,source,score);''')
                con.execute('INSERT INTO method_runs VALUES(?,?,?,?,?,?,?)',
                            ('test','GRCh37','1','ibdmix','complete',1,str(metadata)))
                for i in range(60):
                    con.execute('INSERT INTO segments_evidence VALUES(?,?,?,?,?,?,?,?,?,?,?)',
                                ('test','GRCh37','ibdmix','scientific_reference','Neanderthal','1',str(i),100000,250000,'Altai',10))
                con.execute('INSERT INTO segments_evidence VALUES(?,?,?,?,?,?,?,?,?,?,?)',
                            ('test','GRCh37','ibdmix','ibdmix_reproduction_only','Neanderthal','1','999',0,249250621,'Altai.2013',100))
            path = loci.prepare(database, root/'cache', 'test', 'GRCh37')
            result = json.loads(path.read_text())
            self.assertEqual(len(result['loci']), 1)
            self.assertEqual(result['loci'][0]['n_tested'], 1000)
            self.assertEqual(result['loci'][0]['reference_names'], 'Altai')
            self.assertAlmostEqual(result['loci'][0]['peak_fraction'], .03)
            # A cache hit must not rescan calls; changed zero-call roster must invalidate it.
            with patch.object(loci, 'window_profiles', side_effect=AssertionError('cache missed')):
                self.assertEqual(loci.prepare(database,root/'cache','test','GRCh37'), path)
            roster.write_text(roster.read_text()+'\n1000')
            new = loci.prepare(database,root/'cache','test','GRCh37')
            self.assertNotEqual(new,path)
            self.assertEqual(json.loads(new.read_text())['loci'][0]['n_tested'],1001)
            roster.unlink()
            missing = json.loads(loci.prepare(database,root/'cache','test','GRCh37').read_text())
            self.assertEqual(missing['loci'],[])
            self.assertTrue(missing['metadata']['warnings'])


if __name__ == '__main__':
    unittest.main()
