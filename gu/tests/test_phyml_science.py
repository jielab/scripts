#!/usr/bin/env python3
"""Scientific regression cases for allele identity and saved-tree interpretation."""
import importlib.util
import io
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
import hashlib

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("phyml_science_test", ROOT / "f/phyml.py")
p = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = p
spec.loader.exec_module(p)
audit = p.load_module("phyml.audit.py")


class AlleleIdentity(unittest.TestCase):
	def test_build_reference_switch_keeps_biological_effect(self):
		row = dict(lead_pos=123, ref="A", alt="G", effect_allele="G", beta_j="0.5", index_snp="8:123:A:G")
		f = ["8", "123", "rs123", "G", "A", "0|1", "1|1"]
		lead = p.exact_lead(row, [f], 2)
		self.assertEqual(lead["risk_allele"], "G")
		self.assertEqual(lead["copies"], [1, 0, 0, 0])
		self.assertEqual(lead["allele_order"], "swapped")
		lead = p.exact_lead(dict(row, beta_j="-0.5"), [f], 2)
		self.assertEqual(lead["copies"], [0, 1, 1, 1])

	def test_third_allele_and_ambiguous_records_still_fail(self):
		row = dict(lead_pos=123, ref="A", alt="G", effect_allele="G", beta_j="0.5", index_snp="rs123")
		for records in ([['8','123','rs123','A','C','0|1']],
			[['8','123','rs123','A','G,C','0|1']],
			[['8','123','rs123','A','G','0|1']] * 2):
			with self.assertRaises(p.SkipLocus): p.exact_lead(row, records, 1)


class TreeInterpretation(unittest.TestCase):
	def setUp(self):
		self.haps = [dict(hap_id="H1", role="risk", seq="AA"), dict(hap_id="H2", role="nonrisk", seq="GG")]
		self.tree = "((H2,Altai,Chagyr,Vindija)81,(Denisova,Denisova25)90,H1,Ancestral);"

	def test_nonrisk_signal_is_separate_from_original_risk_call(self):
		d = p.allele_topology_audit(self.tree, self.haps, "Neanderthal")
		self.assertEqual(d["supported_allele_role"], "nonrisk")
		self.assertEqual(d["nonrisk_tree_bootstrap"], 81)
		self.assertIsNone(p.risk_clade(self.tree, {"H1"}, {"H1", "H2"}))
		self.assertIsNone(p.allele_topology_audit(self.tree, self.haps, "Neanderthal", False)["nonrisk_tree_pass"])

	def test_unrooted_complement_and_extra_control(self):
		other_root = "(H2,Altai,Chagyr,Vindija,((Denisova,Denisova25)90,H1,Ancestral)81);"
		self.assertEqual(p.allele_topology_audit(other_root, self.haps, "Neanderthal")["nonrisk_tree_bootstrap"], 81)
		with_extra = "((H2,H1,Altai,Chagyr,Vindija)99,Denisova,Denisova25,Ancestral);"
		self.assertEqual(p.allele_topology_audit(with_extra, self.haps, "Neanderthal")["either_allele_tree_pass"], 0)

	def test_lineage_stage_is_recomputed(self):
		tree = dict(tree_status="complete", tree_newick="((H1,Denisova,Denisova25)90,H2,Altai,Chagyr,Vindija,Ancestral);")
		ident = dict(status="tree_not_supported", tree_pass=0, tree_bootstrap=None, n_sites=20, topology_status="not_supported")
		summaries, _, _ = p.lineage_results(ident, tree, [], self.haps, {})
		self.assertEqual(summaries[0]["topology_status"], "not_supported")
		self.assertEqual(summaries[1]["topology_status"], "supported")
		self.assertEqual(summaries[1]["tree_pass"], 1)

	def test_partial_reference_affinity_is_not_strict_support(self):
		tree = "((H1,Vindija)95,H2,Altai,Chagyr,Denisova,Denisova25,Ancestral);"
		self.assertEqual(p.allele_topology_audit(tree, self.haps, "Neanderthal")["either_allele_tree_pass"], 0)
		edges = audit.affinity_edges(p, tree, self.haps, "Neanderthal")
		self.assertEqual(len(edges), 1)
		self.assertIn("not_introgression_call", edges[0]["evidence"])


class ArchiveIntegrity(unittest.TestCase):
	def test_corrupt_member_is_rejected_without_extraction(self):
		with tempfile.TemporaryDirectory() as directory:
			path = Path(directory) / "test.tar.gz"
			name = "final/gwas_loci.tsv"
			manifest = dict(files=[dict(path=name, sha256=hashlib.sha256(b"original").hexdigest())])
			with tarfile.open(path, "w:gz") as archive:
				for key, value in ((name, b"changed"), ("GU-MANIFEST.json", json.dumps(manifest).encode())):
					info = tarfile.TarInfo(key); info.size = len(value)
					archive.addfile(info, io.BytesIO(value))
			with self.assertRaisesRegex(ValueError, "checksum mismatch"): audit.archive_tables(path)
			self.assertFalse((Path(directory) / "final").exists())


if __name__ == "__main__": unittest.main()
