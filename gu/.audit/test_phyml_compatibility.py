import argparse
import math
import sys
import unittest
import tempfile
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'f'))
from phyml_thresholds import apply_profile, ils_probability
from phyml_evidence import add_common_arguments, validate_args, best_contiguous_cluster
from phyml_anchor import find_anchor_row, read_tsv
from phyml_evidence import write_tsv
from phyml_risk import summarize_risk


class CompatibilityTests(unittest.TestCase):
    def parse(self, *args):
        p = argparse.ArgumentParser()
        add_common_arguments(p)
        a = p.parse_args(['--out', '/tmp/test', *args])
        validate_args(a)
        return a

    def test_old_excel_ils_with_stable_tail(self):
        self.assertAlmostEqual(ils_probability(111969) / 1.05e-8, 1, delta=.01)
        self.assertGreater(ils_probability(15735), .1)  # observed ABO motif
        self.assertLess(ils_probability(44293), .1)    # COVID positive control
        self.assertGreater(ils_probability(644347), 0) # old 1-CDF printed 0

    def test_profiles_and_explicit_override(self):
        a = self.parse()
        self.assertEqual(a.min_candidate_copies, 11)
        self.assertEqual(a.max_ils_probability, .1)
        self.assertEqual(self.parse('--min-candidate-copies', '17').min_candidate_copies, 17)
        old = self.parse('--threshold-profile', 'strict')
        self.assertEqual((old.min_candidate_copies, old.min_diagnostic_match_prop,
                          old.min_candidate_purity, old.max_diagnostic_gap_bp), (2, .8, .8, 50000))

    def test_legacy_no_added_gap_gate(self):
        self.assertEqual(best_contiguous_cluster([100, 110, 90000], 0), [100, 110, 90000])
        self.assertEqual(best_contiguous_cluster([100, 110, 90000], 50000), [100, 110])

    def test_alias_requires_build_and_both_alleles(self):
        aliases = read_tsv(Path(__file__).resolve().parents[1] / 'f/phyml_markers.tsv')
        locus = dict(start=66186240, end=66423003, name='rs6525167', locus_id='rs6525167', chrom='X')
        row = ['X', '66423003', 'X:66423003:G:A', 'G', 'A', '1']
        with patch('phyml_anchor.query_rows', return_value=(['s'], [row])):
            self.assertEqual(find_anchor_row(Path('x'), 'X', locus, aliases, 'GRCh37')[1], row)
            self.assertIsNone(find_anchor_row(Path('x'), 'X', locus, aliases, 'GRCh38')[1])
            row[4] = 'T'
            self.assertIsNone(find_anchor_row(Path('x'), 'X', locus, aliases, 'GRCh37')[1])

    def test_duplicate_ids_remain_ambiguous(self):
        locus = dict(start=1, end=4, name='rs1', locus_id='rs1', chrom='1')
        row = ['1', '2', 'rs1;rs2', 'G', 'A', '1']
        with patch('phyml_anchor.query_rows', return_value=(['s'], [row, row])):
            with self.assertRaises(SystemExit):
                find_anchor_row(Path('x'), '1', locus)

    def test_risk_uses_its_own_lineage_tree_and_not_archaic_anchor(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            def save(name, records):
                write_tsv(root / 'final' / name, list(records[0]), records)
            save('evidence_loci.tsv', [dict(locus_id='rs9349320', name='rs9349320', genome_build='GRCh37',
                 evidence_lineage='Neanderthal', introgression_call='supported_candidate')])
            save('anchor_copies.tsv', [dict(locus_id='rs9349320', sample='s', haplotype='1', hap_id='H1', allele='G')])
            save('evidence_haplotypes.tsv', [dict(locus_id='rs9349320', hap_id='H1', diagnostic_lineage='Denisovan', diagnostic_candidate_pass='1')])
            save('evidence_trees.tsv', [dict(locus_id='rs9349320', expected_lineage='Denisovan', candidate_clade_pass='1', candidate_clade_bootstrap='97')])
            row = summarize_risk(root, Path(__file__).resolve().parents[1] / 'f/phyml_markers.tsv')[0]
            self.assertEqual((row['status'], row['lineage_call'], row['lineage_tree_bootstrap']),
                             ('risk_concordant', 'supported_candidate', '97'))


if __name__ == '__main__':
    unittest.main()
