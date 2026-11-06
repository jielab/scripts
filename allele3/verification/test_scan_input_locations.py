"""Moved raw inputs retain their original scan identity and remain verifiable."""
import contextlib
import copy
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
import herv_pipeline as pipeline
import scan
import scan_store as store


class ExternalInputLocations(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.modern = self.root / "scan"
        self.modern.mkdir()
        self.current = self.root / "current"
        self.current.mkdir()
        self.old = Path("/unavailable-host/inputs")
        self.files = {}
        for name in ("chr1.vcf.bgz", "reference.fa", "joined.tsv"):
            target = self.current / name
            target.write_text("same original bytes " + name)
            self.files[name] = dict(store.fingerprint(target), path=str(self.old / name), mtime_ns=1)
        identity = dict(inputs=[self.files["chr1.vcf.bgz"]], fasta=self.files["reference.fa"])
        self.raw = dict(status="complete", signature=store.stable_hash(identity), identity=identity,
                        shards=[str(self.modern / "chr1")], options=dict(
                            fasta=str(self.old / "reference.fa"), erv_joined=str(self.old / "joined.tsv"),
                            architecture_annotations=dict(files={"erv_joined": self.files["joined.tsv"]})))
        self.manifest = self.modern / "run.json"
        self.manifest.write_text(json.dumps(self.raw))
        self.before = self.manifest.read_bytes()

    def configure(self):
        return store.configure_input_locations(self.modern, {str(self.old): str(self.current)})

    def test_register_read_and_verify_without_mutating_original_scan(self):
        self.assertEqual(self.configure()["relocated_inputs"], 3)
        with patch.object(store, "fingerprint", side_effect=AssertionError("metadata reader opened raw input")):
            run = store.load_run_metadata(self.modern)
        self.assertEqual(run["identity"], self.raw["identity"])
        self.assertEqual(run["options"]["fasta"], str(self.current / "reference.fa"))
        self.assertEqual(store.resolve_run_path(run, self.old / "chr1.vcf.bgz.tbi"), self.current / "chr1.vcf.bgz.tbi")
        store.verify_recorded_input(run, self.files["reference.fa"])
        store.verify_recorded_input(run, run["options"]["architecture_annotations"]["files"]["erv_joined"])
        self.assertEqual(self.manifest.read_bytes(), self.before)

    def test_different_content_cannot_be_registered(self):
        path = self.current / "reference.fa"
        path.write_text("X" * path.stat().st_size)
        with self.assertRaisesRegex(ValueError, "differs from the completed scan"):
            self.configure()
        self.assertFalse((self.modern / "input_locations.json").exists())

    def test_changes_after_registration_are_rejected(self):
        self.configure()
        run = store.load_run_metadata(self.modern)
        target = self.current / "reference.fa"
        before = target.stat()
        os.utime(target, ns=(before.st_atime_ns, before.st_mtime_ns + 1000000000))
        with self.assertRaisesRegex(ValueError, "Recorded scan input changed"):
            store.verify_recorded_input(run, self.files["reference.fa"])

    def test_wrong_scan_incomplete_mapping_and_conflicting_receipt_are_rejected(self):
        self.configure()
        path = self.modern / "input_locations.json"
        original = json.loads(path.read_text())
        for kind in ("scan", "missing", "changed"):
            value = copy.deepcopy(original)
            if kind == "scan": value["run_signature"] = "wrong"
            elif kind == "missing": value["files"].pop()
            else: value["files"][0]["destination"]["edge_sha256"] = "different"
            path.write_text(json.dumps(value))
            with self.subTest(kind=kind), self.assertRaises(ValueError):
                store.load_run_metadata(self.modern)

    def test_annotation_cache_can_validate_its_external_source_after_a_move(self):
        self.configure()
        cache = self.modern / "annotations"
        cache.mkdir()
        output = cache / "chr1.tsv.gz"
        output.write_bytes(b"existing cached annotation")
        # The cache supplies the store context for a source outside the store.
        (cache / "manifest.json").write_text(json.dumps(dict(source=self.files["joined.tsv"],
            files={output.name: store.fingerprint(output)})))
        scan.validate_existing_joined_cache(self.current / "joined.tsv", cache)
        self.assertEqual(self.manifest.read_bytes(), self.before)


class DefaultAllExecution(unittest.TestCase):
    def test_native_default_stops_after_failure_and_records_remaining_modules(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out"
            args = pipeline.parser().parse_args(["--out", str(out)])
            steps = [pipeline.Stage(name, [], [], out / name, (), ()) for name in pipeline.DEFAULT_STAGES]
            def plan(args, paths, names):
                if names == ["summary"]: raise ValueError("Summary inputs not finalized")
                return [step for step in steps if step.name in names]
            with patch.object(pipeline, "plan", side_effect=plan), \
                    patch.object(pipeline, "run_stage", return_value=dict(stage="primary", status="complete")) as run, \
                    contextlib.redirect_stdout(io.StringIO()):
                code = pipeline.execute(args)
            self.assertEqual(code, 2)
            self.assertEqual(run.call_count, 1)
            status = json.loads((out / "pipeline_status.json").read_text())
            self.assertEqual([row["stage"] for row in status["stages"]], list(pipeline.DEFAULT_STAGES))
            self.assertEqual(status["stages"][1]["status"], "failed")
            self.assertTrue(all(row["status"] == "not_available" for row in status["stages"][2:]))


if __name__ == "__main__":
    unittest.main()
