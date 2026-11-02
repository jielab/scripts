import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import audit_integrity as audit
from hpc_common import fingerprint, sha256


class IntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.cfg = {key: str(self.root / name) for key, name in (
            ("reference_raw", "raw"), ("reads_raw", "reads"), ("metadata_dir", "metadata"),
            ("library_dir", "lib"), ("work_dir", "work"))}
        self.db = self.root / "lib/db"
        (self.db / "taxonomy").mkdir(parents=True)
        files = {}
        for name in ("hash.k2d", "taxo.k2d", "opts.k2d"):
            path = self.db / name
            path.write_bytes((name * 4).encode())
            files[name] = {"state": {"size": path.stat().st_size}, "sha256": sha256(path)}
        taxonomy = {}
        for name in ("nodes.dmp", "names.dmp", "merged.dmp"):
            path = self.db / "taxonomy" / name
            path.write_bytes(name.encode())
            taxonomy[name] = {"state": {"size": path.stat().st_size}, "sha256": sha256(path)}
        context = {"files": files, "taxonomy_files": taxonomy}
        receipt = {"schema": "metag-hpc-database-v1", "status": "complete", "context": context,
                   "identity": fingerprint(context), "files": files, "taxonomy_files": taxonomy,
                   "db_dir": str(self.db), "taxonomy_dir": str(self.db / "taxonomy")}
        (self.db / "database.json").write_text(json.dumps(receipt))
        self.seal_path = self.root / "work/seal.json"

    def test_report_is_explicitly_not_hash_validation(self):
        with patch.object(audit, "sha256", side_effect=AssertionError("unexpected hash")):
            result = audit.report(self.cfg)
        self.assertFalse(result["content_hashes_checked"])
        self.assertEqual(result["historical_integrity"], "not_established")
        self.assertEqual(result["files_count"], 7)

    def test_known_hash_mismatch_prevents_publication(self):
        target = self.db / "hash.k2d"
        old = target.read_bytes()
        target.write_bytes(b"X" + old[1:])
        with self.assertRaisesRegex(RuntimeError, "Historical recorded SHA256 mismatch"):
            audit.seal(self.cfg, self.seal_path, allow_local=True)
        self.assertFalse(self.seal_path.exists())

    def test_unknown_legacy_files_are_baseline_only(self):
        result = audit.seal(self.cfg, self.seal_path, allow_local=True)
        self.assertEqual(result["baseline_only_files"], 1)
        saved = json.loads(self.seal_path.read_text())
        receipt = next(row for row in saved["context"]["files"] if row["path"].endswith("database.json"))
        self.assertEqual(receipt["historical_evidence"], audit.BASELINE_ONLY)
        self.assertEqual(audit.verify_seal(self.seal_path)["status"], "pass")
        self.assertTrue(audit.verify_seal(self.seal_path, full=True)["content_hashes_checked"])

    def test_equal_size_modification_is_caught_by_quick_verification(self):
        audit.seal(self.cfg, self.seal_path, allow_local=True)
        target = self.db / "hash.k2d"
        previous = target.stat()
        old = target.read_bytes()
        target.write_bytes(b"Y" + old[1:])
        os.utime(target, ns=(previous.st_atime_ns, previous.st_mtime_ns + 1_000_000_000))
        with self.assertRaisesRegex(RuntimeError, "size/mtime/ctime changed"):
            audit.verify_seal(self.seal_path)

    def test_manifest_mutation_is_rejected(self):
        audit.seal(self.cfg, self.seal_path, allow_local=True)
        saved = json.loads(self.seal_path.read_text())
        saved["context"]["files"][0]["sha256"] = "0" * 64
        self.seal_path.write_text(json.dumps(saved))
        with self.assertRaisesRegex(RuntimeError, "fingerprint"):
            audit.verify_seal(self.seal_path)

    def test_new_index_shard_is_detected(self):
        audit.seal(self.cfg, self.seal_path, allow_local=True)
        (self.db / "hash_extra.k2d").write_bytes(b"extra")
        with self.assertRaisesRegex(RuntimeError, "directory membership changed"):
            audit.verify_seal(self.seal_path)

    def test_seal_never_overwrites_original_or_existing_seal(self):
        with self.assertRaisesRegex(RuntimeError, "outside protected"):
            audit.seal(self.cfg, self.db / "new-seal.json", allow_local=True)
        audit.seal(self.cfg, self.seal_path, allow_local=True)
        before = self.seal_path.read_bytes()
        with self.assertRaisesRegex(RuntimeError, "overwrite"):
            audit.seal(self.cfg, self.seal_path, allow_local=True)
        self.assertEqual(self.seal_path.read_bytes(), before)

    def test_hashes_each_selected_file_once(self):
        paths = [row["path"] for row in audit.collect_files(self.cfg)["files"]]
        with patch.object(audit, "sha256", wraps=sha256) as hashing:
            audit.seal(self.cfg, self.seal_path, allow_local=True)
        self.assertEqual(sorted(str(call.args[0]) for call in hashing.call_args_list), sorted(paths))

    def test_full_rehash_detects_corruption_when_file_states_appear_unchanged(self):
        audit.seal(self.cfg, self.seal_path, allow_local=True)
        target = self.db / "hash.k2d"
        saved_state = audit._state(target)
        original_state = audit._state
        target.write_bytes(b"Z" + target.read_bytes()[1:])
        # Model bit rot below the filesystem metadata layer: stat appears unchanged.
        def unchanged_state(path):
            return saved_state if Path(path) == target else original_state(path)
        with patch.object(audit, "_state", side_effect=unchanged_state):
            self.assertEqual(audit.verify_seal(self.seal_path)["status"], "pass")
            with self.assertRaisesRegex(RuntimeError, "Sealed SHA256 mismatch"):
                audit.verify_seal(self.seal_path, full=True)

    def _reference_fixture(self):
        raw = Path(self.cfg["reference_raw"])
        raw.mkdir()
        fasta = raw / "GCF_000000001.1_fixture.fna"
        fasta.write_text(">one\nACGTACGT\n")
        raw_rows = [{"path": fasta.name, "accession": "GCF_000000001.1",
                     "state": {"size": fasta.stat().st_size}}]
        meta = Path(self.cfg["metadata_dir"])
        (meta / "taxonomy").mkdir(parents=True)
        metadata_files = {}
        for name in ("taxonomy/nodes.dmp", "taxonomy/names.dmp", "taxonomy/merged.dmp", "reference_taxids.tsv"):
            path = meta / name
            path.write_text(name + "\n")
            metadata_files[name] = {"state": {"size": path.stat().st_size}}
        metadata_context = {"files": metadata_files, "raw_inventory": fingerprint(raw_rows), "custom_map": None}
        metadata = {"schema": "metag-hpc-reference-v1", "status": "complete", "context": metadata_context,
                    "identity": fingerprint(metadata_context)}
        (meta / "metadata.json").write_text(json.dumps(metadata))
        lib = Path(self.cfg["library_dir"])
        manifests = lib / "reference_manifests"
        manifests.mkdir()
        manifest_path = manifests / "task000001.json"
        manifest_path.write_text(json.dumps({"task_id": 1, "files": [{**raw_rows[0], "index": 0, "taxid": 10}]}))
        task = {"task_id": 1, "files": [0], "manifest": {
            "path": str(manifest_path.relative_to(lib)), "state": {"size": manifest_path.stat().st_size},
            "sha256": sha256(manifest_path)}}
        plan_context = {"tasks": [task], "reference_raw": str(raw),
                        "metadata_identity": metadata["identity"], "raw_inventory": fingerprint(raw_rows)}
        plan = {"schema": "metag-hpc-reference-v1", "status": "complete", "context": plan_context,
                "tasks": [task], "identity": fingerprint(plan_context)}
        (lib / "reference_plan.json").write_text(json.dumps(plan))
        prepared = lib / "prepared/task000001"
        prepared.mkdir(parents=True)
        prepared_files = {}
        for name in ("sequences.fna", "seqid2taxid.map"):
            path = prepared / name
            path.write_text(name + "\n")
            prepared_files[name] = {"state": {"size": path.stat().st_size}}
        context = {"task": task, "files": prepared_files}
        receipt = {"schema": "metag-hpc-reference-v1", "status": "complete", "context": context,
                   "identity": fingerprint(context), "plan_identity": plan["identity"], "task_id": 1}
        (prepared / "receipt.json").write_text(json.dumps(receipt))
        return fasta

    def test_references_scope_covers_legacy_raw_prepared_and_metadata(self):
        fasta = self._reference_fixture()
        audit.seal(self.cfg, self.seal_path, "references", allow_local=True)
        saved = json.loads(self.seal_path.read_text())
        rows = saved["context"]["files"]
        raw_row = next(row for row in rows if row["path"] == str(fasta))
        self.assertEqual(raw_row["historical_evidence"], audit.BASELINE_ONLY)
        roles = {role for row in rows for role in row["roles"]}
        self.assertTrue({"prepared_reference", "reference_metadata", "raw_reference", "reference_task_manifest"} <= roles)
        self.assertEqual(audit.verify_seal(self.seal_path, full=True)["status"], "pass")
        (fasta.parent / "extra.fna").write_text(">new\nTGCA\n")
        with self.assertRaisesRegex(RuntimeError, "directory membership changed"):
            audit.verify_seal(self.seal_path)


if __name__ == "__main__":
    unittest.main()
