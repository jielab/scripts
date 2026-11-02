"""Export audited nucleotide identities without reopening chains or genotypes."""
from contextlib import redirect_stdout
import csv
import gzip
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"f"))
import mapped_targets as mt
from analysis_integrity import inventory


def mapped_row(pos=10, pos37=20, ref="A", alts="C,G", ref37="A", strand="+", site_id=None):
    bases37 = {mt.oriented(b, strand) for b in {ref, *alts.split(",")}}
    return dict(site_id=site_id or f"1:{pos}:{ref}", chrom="1", pos=str(pos), ref=ref, alts=alts,
                mapping_status="mapped", chrom38="1", pos38=str(pos), ref38=ref, chrom37="1", pos37=str(pos37),
                ref37=ref37, alts37=",".join(sorted(bases37-{ref37})), strand37to38=strand,
                source_ref38=ref, source_to38="+", reference_changed=str(int(mt.oriented(ref37, strand) != ref)))


class MappedTargetsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.time = self.root/"time"
        self.time.mkdir()
        self.targets = self.root/"targets.tsv"
        self.out = self.root/"export"

    def tearDown(self):
        self.temp.cleanup()

    def fixture(self, rows, targets=None):
        fields = sorted(set().union(*(row.keys() for row in rows)))
        with gzip.open(self.time/"mapped_catalog.tsv.gz", "wt") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t")
            writer.writeheader(); writer.writerows(rows)
        (self.time/"other_time_output.tsv").write_text("status\ncomplete\n")
        receipt = dict(status="complete_descriptive_all_allele_panel", version="herv-all-allele-time-2.0-nucleotide-counts",
                       parameters=dict(pgen_prefix="/absent/original-panel", fasta37="/absent/ref37", chain38to37="/absent/chain"),
                       output_inventory=inventory(self.time, "manifest.json"))
        (self.time/"manifest.json").write_text(json.dumps(receipt))
        targets = targets if targets is not None else [dict(site_id=r["site_id"], chrom=r["chrom"], pos=r["pos"], ref=r["ref"], alts=r["alts"],
                                                           family="ERV1", group="case", match_id="pair1") for r in rows]
        fields = sorted(set().union(*(row.keys() for row in targets)))
        with self.targets.open("w") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t")
            writer.writeheader(); writer.writerows(targets)

    def run_export(self, **changes):
        args = mt.parser().parse_args(["--time-dir", str(self.time), "--targets", str(self.targets), "--out", str(self.out)])
        for key, value in changes.items():
            setattr(args, key, value)
        with redirect_stdout(io.StringIO()) as screen:
            result = mt.run(args)
        self.assertIn("[dataset]", screen.getvalue())
        self.assertIn("[build-check]", screen.getvalue())
        return result

    def table(self, name):
        with (self.out/name).open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_plus_strand_keeps_every_requested_base_and_sorts_alt_labels(self):
        self.fixture([mapped_row(alts="G,C")])
        result = self.run_export()
        row = self.table("targets_grch37.tsv")[0]
        self.assertEqual((row["chrom"], row["pos"], row["ref"], row["alts"], row["build"]), ("1", "20", "A", "C,G", "GRCh37"))
        self.assertEqual(row["site_id"], "1:10:A")
        self.assertEqual({r["base38"]:r["base37"] for r in self.table("base_mapping.tsv")}, {"A":"A", "C":"C", "G":"G"})
        self.assertEqual(result["status"], "complete")
        self.assertEqual(result["n_base_mappings"], 3)

    def test_negative_strand_complements_ref_and_all_requested_alts(self):
        self.fixture([mapped_row(ref37="T", strand="-")])
        self.run_export()
        row = self.table("targets_grch37.tsv")[0]
        self.assertEqual((row["ref"], row["alts"], row["reference_changed"]), ("T", "C,G", "0"))
        self.assertEqual(row["requested_bases"], "C,G,T")
        mapping = self.table("base_mapping.tsv")
        self.assertEqual({r["base38"]:r["base37"] for r in mapping}, {"A":"T", "C":"G", "G":"C"})
        self.assertTrue(all(r["complemented"] == "1" for r in mapping))

    def test_true_reference_switch_retains_old_ref_as_new_alt(self):
        self.fixture([mapped_row(ref37="C")])
        self.run_export()
        row = self.table("targets_grch37.tsv")[0]
        self.assertEqual((row["ref"], row["alts"], row["reference_changed"]), ("C", "A,G", "1"))
        mapping = {r["base38"]:r for r in self.table("base_mapping.tsv")}
        self.assertEqual((mapping["A"]["role38"], mapping["A"]["role37"]), ("REF_label", "ALT"))
        self.assertEqual((mapping["C"]["role38"], mapping["C"]["role37"]), ("ALT", "REF_label"))

    def test_reference_outside_requested_bases_is_only_an_assembly_label(self):
        self.fixture([mapped_row(ref37="T")])
        self.run_export()
        row = self.table("targets_grch37.tsv")[0]
        self.assertEqual((row["ref"], row["alts"], row["grch37_ref_among_requested_bases"]), ("T", "A,C,G", "0"))
        self.assertEqual(row["requested_bases"], "A,C,G")
        self.assertNotIn("T", row["requested_bases"].split(","))
        self.assertEqual({r["base37"] for r in self.table("base_mapping.tsv")}, {"A", "C", "G"})

    def test_actual_producer_map_targets_semantics_match_negative_ref_switch(self):
        import herv_time
        class FakeMapper:
            def map(self, *args):
                return dict(chrom38="1", pos38=10, ref38="A", chrom37="1", pos37=20, ref37="G",
                            strand37to38="-", source_ref38="A", source_to38="+", reference_changed=1)
        native = herv_time.map_targets([dict(site_id="1:10:A", chrom="1", pos=10, ref="A", alts="C,G")], FakeMapper())
        self.fixture(native)
        self.run_export()
        row = self.table("targets_grch37.tsv")[0]
        self.assertEqual((row["ref"], row["alts"]), ("G", "C,T"))
        mapping = {r["base38"]:r for r in self.table("base_mapping.tsv")}
        self.assertEqual((mapping["A"]["base37"], mapping["A"]["role37"]), ("T", "ALT"))

    def test_candidate_subset_uses_requested_labels_only_including_ref(self):
        target = dict(site_id="1:10:A", chrom="1", pos=10, ref="A", alts="C")
        self.fixture([mapped_row()], [target])
        self.run_export()
        self.assertEqual(self.table("targets_grch37.tsv")[0]["alts"], "C")
        self.assertEqual(self.table("targets_grch37.tsv")[0]["requested_bases"], "A,C")
        self.assertEqual({r["base38"] for r in self.table("base_mapping.tsv")}, {"A", "C"})

    def test_requested_subset_excludes_structural_ref_not_originally_requested(self):
        target = dict(site_id="1:10:A", chrom="1", pos=10, ref="A", alts="C")
        self.fixture([mapped_row(ref37="T")], [target])
        self.run_export()
        row = self.table("targets_grch37.tsv")[0]
        self.assertEqual((row["ref"], row["alts"], row["requested_bases"]), ("T", "A,C", "A,C"))
        self.assertEqual({r["base37"] for r in self.table("base_mapping.tsv")}, set(row["requested_bases"].split(",")))

    def test_missing_ambiguous_and_reference_mismatch_are_exclusions_not_zeros(self):
        rows = [mapped_row(), dict(site_id="1:11:A", chrom="1", pos="11", ref="A", alts="C,G", mapping_status="ambiguous_mapping")]
        targets = [dict(site_id="1:10:A",chrom="1",pos=10,ref="A",alts="C,G"),
                   dict(site_id="1:11:A",chrom="1",pos=11,ref="A",alts="C,G"),
                   dict(site_id="1:12:A",chrom="1",pos=12,ref="A",alts="C,G")]
        self.fixture(rows, targets)
        result = self.run_export()
        self.assertEqual((result["status"], result["n_exported"], result["n_excluded"]), ("complete_with_exclusions",1,2))
        reasons = {r["site_id"]:r["reason"] for r in self.table("exclusions.tsv")}
        self.assertEqual(reasons, {"1:11:A":"source_mapping_unavailable", "1:12:A":"not_in_time_catalog"})
        self.assertEqual(len(self.table("targets_grch37.tsv")), 1)

    def test_target_reference_change_is_not_confused_with_assembly_reference_switch(self):
        target = dict(site_id="1:10:A", chrom="1", pos=10, ref="C", alts="A,G")
        self.fixture([mapped_row(ref37="C")], [target])
        result = self.run_export()
        self.assertEqual(result["status"], "not_available")
        self.assertEqual(self.table("exclusions.tsv")[0]["reason"], "source38_reference_mismatch")
        self.assertEqual(self.table("targets_grch37.tsv"), [])

    def test_missing_requested_nucleotide_cannot_be_dropped_silently(self):
        target = dict(site_id="1:10:A", chrom="1", pos=10, ref="A", alts="C,G,T")
        self.fixture([mapped_row()], [target])
        result = self.run_export()
        self.assertEqual(result["status"], "not_available")
        self.assertEqual(self.table("exclusions.tsv")[0]["reason"], "requested_bases_missing_from_source")

    def test_identical_source_and_candidate_duplicates_deduplicate(self):
        row = mapped_row()
        self.fixture([row, dict(row)])
        result = self.run_export()
        self.assertEqual(result["n_exported"], 1)
        self.assertEqual(result["identical_candidate_duplicates"], 1)
        self.assertEqual(result["stream_counts"]["identical_source_duplicates"], 1)

    def test_conflicting_source_site_id_is_an_error(self):
        self.fixture([mapped_row(), mapped_row(pos37=30)])
        with self.assertRaisesRegex(ValueError, "Conflicting duplicate source site_id"):
            self.run_export()

    def test_conflicting_candidate_site_id_is_an_error(self):
        candidates = [dict(site_id="one",chrom="1",pos=10,ref="A",alts="C,G"),
                      dict(site_id="one",chrom="1",pos=11,ref="A",alts="C,G")]
        self.fixture([mapped_row()], candidates)
        with self.assertRaisesRegex(ValueError, "Conflicting duplicate candidate site_id"):
            self.run_export()

    def test_distinct_original_sites_colliding_in_37_are_both_excluded(self):
        self.fixture([mapped_row(pos=10), mapped_row(pos=11)])
        result = self.run_export()
        self.assertEqual(result["status"], "not_available")
        self.assertEqual(result["n_excluded"], 2)
        self.assertTrue(all(r["reason"] == "mapped_coordinate_collision" for r in self.table("exclusions.tsv")))

    def test_inconsistent_source_allele_labels_fail(self):
        row = mapped_row(strand="-",ref37="T"); row["alts37"]="A,C"
        self.fixture([row])
        with self.assertRaisesRegex(ValueError, "Source alts37 disagrees"):
            self.run_export()

    def test_inventory_tampering_and_unlisted_mapping_fail(self):
        self.fixture([mapped_row()])
        (self.time/"other_time_output.tsv").write_text("changed\n")
        with self.assertRaisesRegex(ValueError, "Output identity changed"):
            self.run_export()
        self.fixture([mapped_row()])
        path = self.time/"manifest.json"; receipt=json.loads(path.read_text())
        receipt["output_inventory"].pop("mapped_catalog.tsv.gz"); path.write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError, "absent from the verified time output inventory"):
            self.run_export()

    def test_missing_unfinished_or_wrong_producer_manifest_fails(self):
        self.fixture([mapped_row()])
        path = self.time/"manifest.json"; original=json.loads(path.read_text())
        path.unlink()
        with self.assertRaises(FileNotFoundError): self.run_export()
        path.write_text(json.dumps(dict(original,status="running")))
        with self.assertRaisesRegex(ValueError, "Nonterminal/unusable"): self.run_export()
        path.write_text(json.dumps(dict(original,version="unrelated-producer")))
        with self.assertRaisesRegex(ValueError, "producer manifest"): self.run_export()

    def test_cap_and_build_declaration(self):
        self.fixture([mapped_row(), mapped_row(pos=11,pos37=21)])
        with self.assertRaisesRegex(ValueError, "More than 1"): self.run_export(max_targets=1)
        candidate = dict(site_id="1:10:A",chrom="1",pos=10,ref="A",alts="C,G",build="GRCh37")
        self.fixture([mapped_row()], [candidate])
        with self.assertRaisesRegex(ValueError, "must be explicitly GRCh38"): self.run_export()
        candidate.pop("build"); candidate["target_build"] = "GRCh37"
        self.fixture([mapped_row()], [candidate])
        with self.assertRaisesRegex(ValueError, "must be explicitly GRCh38"): self.run_export()

    def test_output_protects_time_directory_and_candidate_inputs(self):
        self.fixture([mapped_row()])
        before = (self.time/"manifest.json").read_bytes()
        for destination in (self.time, self.time/"nested", self.root):
            with self.subTest(destination=destination):
                with self.assertRaisesRegex(ValueError, "overlaps protected"):
                    self.run_export(out=str(destination))
        self.assertEqual((self.time/"manifest.json").read_bytes(), before)

    def test_cli_no_usable_targets_has_audited_nonzero_exit(self):
        self.fixture([dict(site_id="1:10:A",chrom="1",pos=10,ref="A",alts="C,G",mapping_status="unmapped")])
        cp = subprocess.run([sys.executable, str(ROOT/"f/mapped_targets.py"), "--time-dir", str(self.time),
                             "--targets", str(self.targets), "--out", str(self.out)], capture_output=True, text=True)
        self.assertEqual(cp.returncode, 2, cp.stdout+cp.stderr)
        self.assertEqual(json.loads((self.out/"manifest.json").read_text())["status"], "not_available")
        self.assertEqual(len(self.table("exclusions.tsv")), 1)


if __name__ == "__main__":
    unittest.main()
