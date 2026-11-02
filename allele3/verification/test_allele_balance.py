import itertools
import math
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"f"))
from allele_balance import rarefaction, summarize_row


class TestBalance(unittest.TestCase):
    def test_exact_joint_rarefaction_against_all_subsamples(self):
        for counts in ([3,2,1], [2,2,1,1], [4,1], [1,1,1,1]):
            copies=[i for i,c in enumerate(counts) for _ in range(c)]
            for n in range(1,len(copies)+1):
                supports=[len({copies[i] for i in selected}) for selected in itertools.combinations(range(len(copies)),n)]
                p,e=rarefaction(counts,n)
                self.assertAlmostEqual(p,sum(k>=3 for k in supports)/len(supports))
                self.assertAlmostEqual(e,sum(supports)/len(supports))

    def row(self, counts=(90,5,5), ref="A"):
        return dict(site_id="1:1:A",chrom="1",pos="1",ref=ref,family="ERV1",group="case",population="EUR",
                    evaluable="1",represented_status="all_expected_alts_encoded",record_status="ok",encoded_alts="C,G",
                    A_AC=str(counts[0]),C_AC=str(counts[1]),G_AC=str(counts[2]),T_AC="NA",AN=str(sum(counts)),call_rate="1")

    def test_ref_absent_is_only_two_states(self):
        r=summarize_row(self.row((0,60,40)),2,100,.9)
        self.assertEqual(r["at_least_three_supported"],0)
        self.assertEqual(r["rank1_base"],"C")
        self.assertEqual(r["p_rarefied_ge3_each_ac_ge1"],0)

    def test_three_actual_bases_and_ties(self):
        r=summarize_row(self.row(),5,100,.9)
        self.assertEqual(r["at_least_three_supported"],1)
        self.assertEqual(r["p_rarefied_ge3_each_ac_ge1"],1)
        self.assertAlmostEqual(r["entropy_nats"],-.9*math.log(.9)-.1*math.log(.05))

    def test_unencoded_target_not_absent(self):
        row=self.row();row["represented_status"]="expected_alt_not_encoded";row["evaluable"]="0"
        self.assertEqual(summarize_row(row,2,100,.9)["status"],"unavailable")

    def test_invalid_denominator_rejected(self):
        row=self.row();row["AN"]="101"
        with self.assertRaises(ValueError): summarize_row(row,2,100,.9)

    def test_too_few_copies_not_extrapolated(self):
        self.assertEqual(rarefaction([2,1,1],5),(None,None))

    def test_odd_an_retained(self):
        r=summarize_row(self.row((93,5,5)),2,100,.9)
        self.assertEqual(r["AN"],103)


if __name__=="__main__": unittest.main()
