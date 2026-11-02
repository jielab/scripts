"""Bounded data provenance checks without access to any real cohort files."""
from contextlib import redirect_stdout
import csv
import gzip
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"f"))
import data_sources as ds

try:
    import pysam
except ImportError:
    pysam = None


def source(path, build="38", kind="vcf", identifier="test", evidence="aggregate_counts"):
    return dict(id=identifier, type=kind, path=str(path), build=build, evidence=evidence, role="Synthetic audit fixture")


def write_vcf(path, *, reference=None, contig_length=None, sample=False, base="A", n=24):
    text = "##fileformat=VCFv4.2\n"
    if reference:
        text += f"##reference={reference}\n"
    if contig_length:
        text += f"##contig=<ID=1,length={contig_length}>\n"
    text += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" + ("\tFORMAT\ts1\ts2" if sample else "") + "\n"
    for pos in range(1, n+1):
        text += f"1\t{pos}\t.\t{base}\t{'G' if base != 'G' else 'T'}\t.\tPASS\t." + ("\tGT\t0/1\t0/0" if sample else "") + "\n"
    path.write_text(text)


class DataSources(unittest.TestCase):
    def audit(self, entries, **kwargs):
        stream = io.StringIO()
        with redirect_stdout(stream):
            result = ds.audit_entries(entries, **kwargs)
        self.assertIn("[dataset]", stream.getvalue())
        self.assertIn("[build-check]", stream.getvalue())
        return result

    def test_registry_is_explicit_and_ancestry_models_stay_separate(self):
        entries = {row["id"]: row for row in ds.load_config(ROOT/"config/data_sources.tsv")}
        self.assertEqual(entries["gu_1kg"]["build"], "metadata")
        self.assertEqual(entries["ancient_dates"]["build"], "not_applicable")
        for model in ("ancient_ancestry_3way", "ancient_ancestry_4way"):
            self.assertEqual(entries[model]["build"], "not_applicable")
            self.assertIn("never concatenate", entries[model]["notes"])
        self.assertEqual(entries["ukb_imputed"]["evidence"], "imputed_genotypes")
        self.assertEqual(entries["1kg_sites"]["evidence"], "aggregate_counts")

    def registry(self, root, paths):
        config = root/"registry.tsv"
        with config.open("w") as handle:
            writer = csv.DictWriter(handle, ds.REQUIRED, delimiter="\t")
            writer.writeheader()
            for i, path in enumerate(paths):
                record = source(path, identifier=f"fixture{i}")
                record["role"] = "Retain {dir_data} in non-path descriptions"
                writer.writerow(record)
        return config

    def test_registry_hpc_roots_expand_only_declared_tokens_and_keep_chrom_for_units(self):
        with tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ, {
                "DIR_DATA":"/hpc/data cohort", "DIR_WORK":"/hpc/work", "ANALYSIS_ROOT":"/scratch/date/analysis"}):
            config = self.registry(Path(tmp), ["{dir_data}/gen/chr{chrom}.vcf.gz",
                "{dir_work}/metadata.tsv", "{analysis_root}/gu/1kg"])
            entries = ds.load_config(config)
            self.assertEqual([r["path"] for r in entries], ["/hpc/data cohort/gen/chr{chrom}.vcf.gz",
                "/hpc/work/metadata.tsv", "/scratch/date/analysis/gu/1kg"])
            self.assertEqual(str(ds.units(entries[0], ["X"])[0][1]), "/hpc/data cohort/gen/chrX.vcf.gz")
            self.assertEqual(entries[0]["role"], "Retain {dir_data} in non-path descriptions")

    def test_registry_local_defaults_and_correct_actual_gnomad_location(self):
        with patch.dict(os.environ, {}, clear=True):
            entries = {r["id"]: r for r in ds.load_config(ROOT/"config/data_sources.tsv")}
        self.assertEqual(entries["gnomad_scan"]["path"], "/mnt/f/annot/gnomad/allele3")
        self.assertEqual(entries["gnomad_sites"]["path"], "/mnt/f/annot/gnomad/sites")
        self.assertEqual(entries["allofus"]["path"], "/mnt/f/annot/allele3/allofus/hg38_allofus.txt.gz")
        self.assertEqual(entries["gu_1kg"]["path"], "/mnt/f/analysis/gu/1kg")
        self.assertEqual(entries["allofus"]["evidence"], "aggregate_frequencies")

    def test_old_concrete_registry_paths_and_actual_manifest_inputs_stay_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ, {
                "DIR_DATA":"/hpc/data", "DIR_WORK":"/hpc/work", "ANALYSIS_ROOT":"/hpc/analysis"}):
            root = Path(tmp)
            paths = ["/mnt/f/annot/allele3/legacy.vcf", "/mnt/d/explicit/chr{chrom}.vcf", "$DIR_DATA/literal.vcf"]
            entries = ds.load_config(self.registry(root, paths))
            self.assertEqual([r["path"] for r in entries], paths)
            scan = root/"scan"; scan.mkdir()
            (scan/"run.json").write_text(json.dumps(dict(build="GRCh38", inputs=paths[:1])))
            row = self.audit([source(scan, kind="scan")])["datasets"][0]
            self.assertEqual(row["details"]["metadata"]["recorded_input_paths"], paths[:1])

    def test_registry_unknown_placeholders_and_format_expressions_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            for expression in ("{unknown}", "{DIR_DATA}", "{dir_data.real}", "{dir_data[0]}",
                               "{dir_data!r}", "{dir_work:20}", "{chrom!r}", "{dir_data"):
                with self.subTest(expression=expression), self.assertRaises(ValueError):
                    ds.load_config(self.registry(Path(tmp), [expression+"/chr1.vcf"]))
            with self.assertRaisesRegex(ValueError, "Only a literal"):
                ds.check_entries([source("{dir_data}/not_loaded.vcf")])

    def test_fasta_defaults_follow_shared_environment_and_preserve_explicit_paths(self):
        with patch.dict(os.environ, {"DIR_DATA":"/hpc/data"}, clear=True):
            args = ds.parser().parse_args(["audit", "--out", "/tmp/not_written"])
            self.assertEqual(args.fasta37, "/hpc/data/gen/1hgp/37/GRCH37.fasta.gz")
            self.assertEqual(args.fasta38, "/hpc/data/gen/1hgp/38/GRCH38.fasta.gz")
            args = ds.parser().parse_args(["audit", "--out", "/tmp/not_written", "--fasta38", "/explicit/local/reference.fa"])
            self.assertEqual(args.fasta38, "/explicit/local/reference.fa")
        with patch.dict(os.environ, {"DIR_DATA":"/hpc/data", "ALLELE_FASTA37":"/custom/b37.fa", "ALLELE_FASTA":"/custom/b38.fa"}, clear=True):
            args = ds.parser().parse_args(["audit", "--out", "/tmp/not_written"])
            self.assertEqual((args.fasta37, args.fasta38), ("/custom/b37.fa", "/custom/b38.fa"))

    def test_actual_header_distinguishes_sites_and_individuals(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); sites, people = root/"1kg.vcf", root/"gnomad.vcf"
            write_vcf(sites); write_vcf(people, sample=True)
            result = self.audit([source(sites, identifier="sites"), source(people, identifier="people")])
            a, b = result["datasets"]
            self.assertEqual((a["sample_count"], a["genotype_scope"]), (0, "sites_only"))
            self.assertEqual((b["sample_count"], b["genotype_scope"]), (2, "individuals"))
            self.assertIn("no individual columns", result["source_plan"][0]["constraint"])

    def test_path_name_does_not_validate_a_build_and_header_conflict_blocks(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/"hg38_GRCh38.vcf"
            write_vcf(path)
            result = self.audit([source(path)])
            self.assertEqual(result["datasets"][0]["build_status"], "declared_only")
            write_vcf(path, reference="GRCh37")
            result = self.audit([source(path)])
            self.assertTrue(result["has_conflicts"])
            self.assertEqual(result["datasets"][0]["build_status"], "conflict")

    def test_canonical_contig_length_conflict_is_not_ignored(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/"input.vcf"
            write_vcf(path, contig_length=249250621)
            result = self.audit([source(path, "38")])
            self.assertTrue(result["has_conflicts"])
            self.assertEqual(result["datasets"][0]["contig_evidence"]["matches"]["37"], ["1"])

    def test_conflicting_duplicate_contig_headers_are_blocking(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/"input.vcf"
            path.write_text("##contig=<ID=1,length=249250621>\n##contig=<ID=chr1,length=248956422>\n")
            row = self.audit([source(path)])["datasets"][0]
            self.assertEqual(row["status"], "conflict")
            self.assertTrue(row["conflict"])

    def test_missing_sources_are_deferred_without_global_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); vcf = root/"existing.vcf"
            write_vcf(vcf)
            result = self.audit([source(vcf, identifier="ready"), source(root/"missing.vcf", identifier="absent")])
            self.assertFalse(result["has_conflicts"])
            self.assertEqual([r["status"] for r in result["datasets"]], ["audited", "deferred_missing"])

    def test_allofus_is_sniffed_without_assuming_idx_or_schema(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/"data.txt.gz"
            with gzip.open(path, "wt") as handle:
                handle.write("chrom\tpos\tvariant\n1\t100\tA_G\n")
            Path(str(path)+".idx").write_bytes(b"some index")
            row = self.audit([source(path, kind="allofus_text")])["datasets"][0]
            self.assertEqual(row["status"], "requires_schema")
            self.assertEqual(row["details"]["metadata"]["index_status"], "unclassified_idx_not_tabix")
            self.assertIsNone(row["sample_count"])
            self.assertEqual(row["build_status"], "declared_only")

    def test_supplied_allofus_header_reports_af_only_and_does_not_validate_build_from_name(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/"hg38_allofus.txt.gz"
            with gzip.open(path, "wt") as handle:
                handle.write("#"+"\t".join(ds.ALLOFUS_COLUMNS)+"\n")
                for alt in ("A", "C", "G"):
                    handle.write("chr1\t10001\t10001\tT\t"+alt+"\t0.01\t0.02\t0\t0.01\t0\t0.01\t\t0.02\t0\n")
                handle.write("chr1\t10110\t10110\tACCCT\tA\t0.01\t0.02\t0\t0.01\t0\t0.01\t\t0.02\t0\n")
            Path(str(path)+".idx").write_text("unclassified original index")
            result = self.audit([source(path, kind="allofus_text", evidence="aggregate_frequencies")])
            row = result["datasets"][0]; meta = row["details"]["metadata"]
            self.assertEqual(row["status"], "audited_af_only")
            self.assertEqual(row["build_status"], "declared_only")
            self.assertEqual(meta["schema_status"], "annovar-allofus-hg38")
            self.assertEqual(len(meta["populations"]), 8)
            self.assertNotIn("MAX", meta["populations"])
            self.assertFalse(meta["counts_available"])
            self.assertEqual(meta["additional_metrics"], ["gvs_max_af"])
            self.assertEqual(meta["sniff"]["column_counts"], [14]*5)
            self.assertEqual(meta["index_status"], "unclassified_idx_not_tabix")
            self.assertIn("AC>=5", result["source_plan"][0]["constraint"])

    def test_acount_does_not_create_a_genotype_or_wgs_source(self):
        with tempfile.TemporaryDirectory() as temp:
            prefix = Path(temp)/"chr1"
            Path(str(prefix)+".acount").write_text("#CHROM\tID\tREF\tALT\tALT_CTS\tOBS_CT\n")
            result = self.audit([source(prefix, "37", kind="pgen", evidence="imputed_genotypes")])
            row = result["datasets"][0]
            self.assertFalse(row["genotype_files_present"])
            self.assertEqual(row["status"], "deferred_missing")
            self.assertTrue(row["details"]["fileset"]["acount_present"])
            self.assertIn("Neither .acount nor imputed PGEN is WGS", result["source_plan"][0]["constraint"])

    def test_pgen_fileset_presence_is_not_genotype_decoding(self):
        with tempfile.TemporaryDirectory() as temp:
            prefix = Path(temp)/"chr1"
            Path(str(prefix)+".pgen").write_bytes(b"synthetic presence only")
            Path(str(prefix)+".psam").write_text("#IID\ns1\n")
            Path(str(prefix)+".pvar").write_text("##reference=GRCh37\n#CHROM\tPOS\tID\tREF\tALT\n1\t1\t.\tA\tC\n")
            row = self.audit([source(prefix, "37", kind="pgen")])["datasets"][0]
            self.assertTrue(row["genotype_files_present"])
            self.assertEqual(row["genotype_scope"], "individual_fileset_present_not_decoded")
            self.assertFalse(row["details"]["fileset"]["genotype_content_verified"])

    def test_provisional_pvar_ref_is_not_used_to_decide_build(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/"p.pvar"
            path.write_text("#CHROM\tPOS\tID\tREF\tALT\tINFO\n1\t1\t.\tA\tC\tPR\n")
            row = self.audit([source(path, "37", kind="pvar")])["datasets"][0]
            self.assertEqual(row["details"]["metadata"]["sampled_snvs"], 0)
            self.assertEqual(row["details"]["metadata"]["provisional_ref_sites_skipped"], 1)

    def test_gu_build_comes_from_metadata_and_dates_are_not_genome_coordinates(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); gu = root/"ibdmix/chr1"
            gu.mkdir(parents=True)
            (gu/"run.meta.tsv").write_text("genome_build\\tb37\ncoordinate_system\t0-based-half-open\npipeline_version\tunit-test\n")
            dates = root/"qpAdm.tsv"; dates.write_text("IID\tAge\tscore\ns1\t1000\t0.1\n")
            result = self.audit([source(root, "metadata", "gu", "gu"), source(dates, "na", "metadata", "dates")], chromosomes=["1"])
            self.assertEqual(result["datasets"][0]["build_status"], "metadata_supports_37")
            self.assertEqual(result["datasets"][1]["build_status"], "not_applicable")

    def test_header_read_is_bounded_even_before_chrom_line(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/"input.vcf"
            path.write_text("##comment=keep\n"*10 + "##reference=GRCh37\nthis must not be read\n")
            row = self.audit([source(path)], max_records=5)["datasets"][0]
            self.assertEqual(row["status"], "deferred_header_limit")
            self.assertEqual(row["details"]["metadata"]["lines_read"], 5)
            self.assertFalse(row["conflict"])

    def test_scan_records_actual_inputs_without_guessing_vcf_names(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)/"scan"; root.mkdir()
            paths = ["/actual/path/unusually_named_source.vcf.bgz"]
            (root/"run.json").write_text(json.dumps(dict(status="complete", options=dict(build="GRCh38"), inputs=paths)))
            row = self.audit([source(root, "38", "scan")])["datasets"][0]
            self.assertEqual(row["details"]["metadata"]["recorded_input_paths"], paths)
            self.assertFalse(row["details"]["metadata"]["scan_identity_checked"])

    def test_output_cannot_replace_source_directory(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); path = root/"input.vcf"
            write_vcf(path); before = path.read_bytes()
            with self.assertRaisesRegex(ValueError, "overlaps protected"):
                self.audit([source(path)], out=root)
            self.assertEqual(path.read_bytes(), before)

    def test_cli_publishes_conflicts_and_returns_nonzero_but_missing_returns_zero(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); path = root/"wrong.vcf"; write_vcf(path, reference="GRCh37")
            config = root/"sources.tsv"
            with config.open("w") as handle:
                writer = csv.DictWriter(handle, ds.REQUIRED, delimiter="\t")
                writer.writeheader()
                writer.writerow(source(path))
                writer.writerow(source(root/"absent", identifier="missing"))
            command = [sys.executable, str(ROOT/"f/data_sources.py"), "audit", "--config", str(config),
                       "--fasta37", str(root/"absent37.fa"), "--fasta38", str(root/"absent38.fa")]
            out = root/"audit"
            cp = subprocess.run(command+["--out", str(out)], capture_output=True, text=True)
            self.assertEqual(cp.returncode, 2, cp.stdout+cp.stderr)
            report = json.loads((out/"audit.json").read_text())
            self.assertEqual(report["n_conflicts"], 1)
            self.assertEqual(len(report["datasets"]), 2)
            self.assertTrue((out/"audit.tsv").is_file())
            self.assertTrue((out/"source_plan.tsv").is_file())
            missing = subprocess.run(command+["--datasets", "missing", "--out", str(root/"missing_audit")], capture_output=True, text=True)
            self.assertEqual(missing.returncode, 0, missing.stdout+missing.stderr)


@unittest.skipIf(pysam is None, "pysam is required for real indexed FASTA tests")
class IndexedFastaTests(unittest.TestCase):
    audit = DataSources.audit
    def make_reference(self, path, sequence):
        path.write_text(">1\n"+sequence+"\n")
        pysam.faidx(str(path))

    def test_dual_fasta_matches_are_ambiguous_and_not_validated(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); f37, f38, path = root/"b37.fa", root/"b38.fa", root/"v.vcf"
            self.make_reference(f37, "A"*100); self.make_reference(f38, "A"*100); write_vcf(path)
            row = self.audit([source(path)], fasta37=f37, fasta38=f38)["datasets"][0]
            self.assertEqual(row["build_status"], "ambiguous_two_builds")
            self.assertEqual(row["ref_check"]["joint"]["both_match"], 24)
            self.assertNotIn("validated", row["build_status"])

    def test_discriminating_ref_sample_support_and_conflict(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); f37, f38, path = root/"b37.fa", root/"b38.fa", root/"v.vcf"
            self.make_reference(f37, "A"*100); self.make_reference(f38, "C"*100); write_vcf(path)
            yes = self.audit([source(path, "37")], fasta37=f37, fasta38=f38)
            self.assertFalse(yes["has_conflicts"])
            self.assertEqual(yes["datasets"][0]["build_status"], "sample_supports_supplied_reference_37")
            no = self.audit([source(path, "38")], fasta37=f37, fasta38=f38)
            self.assertTrue(no["has_conflicts"])
            self.assertEqual(no["references"]["37"]["assembly_identity"], "unverified_custom_or_partial_reference")

    def test_reference_indexes_are_never_created_implicitly(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); ref, path = root/"noindex.fa", root/"v.vcf"
            ref.write_text(">1\nAAAA\n"); write_vcf(path)
            result = self.audit([source(path)], fasta38=ref)
            self.assertEqual(result["references"]["38"]["status"], "deferred_missing_fasta_or_index")
            self.assertFalse(Path(str(ref)+".fai").exists())

    def test_supplied_fasta_label_conflict_is_reported_before_querying(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); ref, path = root/"wrong37.fa", root/"v.vcf"
            ref.write_text(">1\nAAAA\n")
            Path(str(ref)+".fai").write_text("1\t248956422\t3\t4\t5\n")
            write_vcf(path)
            result = self.audit([source(path, "37")], fasta37=ref)
            self.assertEqual(result["references"]["37"]["status"], "reference_label_conflict")
            self.assertTrue(result["has_conflicts"])


if __name__ == "__main__":
    unittest.main()
