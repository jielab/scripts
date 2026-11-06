"""Relocated scan masks preserve filtering and immutable source provenance."""
import copy
import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
import herv_pipeline
from analysis_integrity import validate_manifest
from scan_input_paths import resolve_mask_options
from scan_store import fingerprint


class ScanMaskRelocationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.mask = self.root / "hpc-mask.bed"
        self.mask.write_text("chr1\t599\t600\n")
        self.old = "/unavailable-local-drive/original-mask.bed"
        record = dict(fingerprint(self.mask), path=self.old, mtime_ns=1)
        self.run = dict(status="complete", signature="fixture", shards=["chr1"], options={
            "build": "GRCh38", "architecture_low_mappability_bed": self.old,
            "architecture_annotations": {"files": {"architecture_low_mappability_bed": record}}})
        self.overrides = {"low_mappability": str(self.mask)}

    def test_copy_with_different_mtime_is_accepted_without_mutating_provenance(self):
        original = copy.deepcopy(self.run)
        resolved = resolve_mask_options(self.run, self.overrides, verify=True)
        self.assertEqual(resolved["architecture_low_mappability_bed"], str(self.mask))
        self.assertEqual(self.run, original)
        self.assertEqual(resolve_mask_options(self.run)["architecture_low_mappability_bed"], self.old)

    def test_changed_mask_is_rejected_even_with_same_size(self):
        self.mask.write_text("chr1\t699\t700\n")
        with self.assertRaisesRegex(ValueError, "differs from original"):
            resolve_mask_options(self.run, self.overrides, verify=True)

    def test_missing_original_filter_or_fingerprint_cannot_be_overridden(self):
        for field in ("architecture_low_mappability_bed", "architecture_annotations"):
            run = copy.deepcopy(self.run)
            del run["options"][field]
            with self.subTest(field=field), self.assertRaises(ValueError):
                resolve_mask_options(run, self.overrides, verify=True)

    def test_preflight_does_not_read_mask_data_but_worker_rejects_missing_copy(self):
        self.mask.unlink()
        with patch("scan_store.fingerprint", side_effect=AssertionError("preflight opened mask")):
            resolved = resolve_mask_options(self.run, self.overrides)
        self.assertEqual(resolved["architecture_low_mappability_bed"], str(self.mask))
        with self.assertRaisesRegex(ValueError, "missing"):
            resolve_mask_options(self.run, self.overrides, verify=True)

    def fixture(self):
        modern = self.root / "scan"
        shard = modern / "chr1"
        shard.mkdir(parents=True)
        (modern / "run.json").write_text(json.dumps(self.run))
        (shard / "scan.json").write_text(json.dumps({"complete": True, "scope": ["1"], "fields_read": ["lcr", "segdup"]}))
        sites = shard / "multiallelic_sites.tsv.gz"
        with gzip.open(sites, "wt") as handle:
            handle.write("chrom\tpos\tref\talts_by_af\talt_ac\talt_an\tcontext\n")
            for pos in (500, 600, 5000):
                handle.write(f"1\t{pos}\tA\tC,G\t10,5\t100,100\tAAA\n")
        (shard / "COMPLETE.json").write_text(json.dumps({"signature": "fixture", "outputs": {sites.name: {"size": sites.stat().st_size}}}))
        annotation = self.root / "joined.tsv"
        annotation.write_text("0\tchr1\t100\t1000\tfixture#LTR/ERV1\t0\t+\t100\t1000\t0\t1\t900\t0\tfixture\n")
        return modern, annotation

    def test_pipeline_preflight_and_worker_use_same_override_for_both_catalog_modes(self):
        modern, annotation = self.fixture()
        for stage in ("candidates", "alleles"):
            args = herv_pipeline.parser().parse_args([stage, "--scan-dir", str(modern),
                "--out", str(self.root / "result"), "--erv-joined", str(annotation),
                "--low-mappability-bed", str(self.mask)])
            step = herv_pipeline.plan(args, herv_pipeline.resolve(args), [stage])[0]
            self.assertIn(self.mask, step.inputs)
            self.assertNotIn(Path(self.old), step.inputs)
            self.assertEqual(herv_pipeline.availability(step)["missing_inputs"], [])
            self.assertEqual(step.command[step.command.index("--low-mappability-bed") + 1], str(self.mask))
            self.assertIn("scan_input_paths.py", step.modules)

    def test_real_candidate_cli_excludes_relocated_mask_and_preserves_source_files(self):
        modern, annotation = self.fixture()
        before = {p: p.read_bytes() for p in modern.rglob("*.json")}
        out = self.root / "candidates"
        command = [sys.executable, str(ROOT / "f/herv_candidates.py"), "--modern-dir", str(modern),
                   "--erv-joined", str(annotation), "--out", str(out), "--controls-per-case", "1",
                   "--low-mappability-bed", str(self.mask)]
        result = subprocess.run(command, capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        manifest = validate_manifest(out, "candidate_manifest.json", require_evidence=True)
        self.assertEqual(manifest["n_cases"], 1)
        self.assertEqual(manifest["n_controls"], 1)
        self.assertEqual(manifest["pass1"]["excluded_difficult_regions"], 1)
        with (out / "targets.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({int(row["pos"]) for row in rows}, {500, 5000})
        self.assertTrue(all(row["strict_filter_complete"] == "1" for row in rows))
        self.assertEqual(before, {p: p.read_bytes() for p in before})


if __name__ == "__main__":
    unittest.main()
