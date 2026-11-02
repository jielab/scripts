"""Reject short native index writes before accepting a complete database."""
import json
from pathlib import Path
import struct
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import hpc_reference as reference
from hpc_common import fingerprint, sha256


class NativeIndexStructureTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.db = self.root / "lib/db"
        (self.db / "taxonomy").mkdir(parents=True)
        self.cfg = {"library_dir": str(self.root / "lib"),
                    "metadata_dir": str(self.root / "metadata"),
                    "reference_raw": str(self.root / "raw")}
        self.hash = self.db / "hash.k2d"
        self.taxonomy = self.db / "taxo.k2d"
        self.hash.write_bytes(struct.pack("<QQQQ", 8, 4, 29, 3) + bytes(8 * 4))
        self.taxonomy.write_bytes(b"K2TAXDAT" + struct.pack("<QQQ", 6, 4, 5)
                                 + bytes(6 * 56 + 4 + 5))
        opts = bytearray(64)
        struct.pack_into("<QQ", opts, 0, 35, 31)
        opts[32] = 1
        (self.db / "opts.k2d").write_bytes(opts)
        for name in reference.TAX_FILES:
            (self.db / "taxonomy" / name).write_text(name + "\n")

    def _marker(self, capacity=8):
        def record(path):
            return {"state": {"size": path.stat().st_size}, "sha256": sha256(path)}
        files = {name: record(self.db / name) for name in ("hash.k2d", "taxo.k2d", "opts.k2d")}
        taxonomy = {name: record(self.db / "taxonomy" / name) for name in reference.TAX_FILES}
        context = {"parameters": dict(reference.PARAMETERS), "files": files,
                   "taxonomy_files": taxonomy}
        if capacity is not None:
            context["resource_estimate"] = {"capacity": capacity}
        marker = {"schema": reference.DB_SCHEMA, "status": "complete", "context": context,
                  "identity": fingerprint(context), "files": files, "taxonomy_files": taxonomy,
                  "db_dir": str(self.db), "taxonomy_dir": str(self.db / "taxonomy"),
                  "parameters": dict(reference.PARAMETERS)}
        (self.db / "database.json").write_text(json.dumps(marker))
        return marker

    def test_valid_native_layout(self):
        result = reference._native_index_structure(self.db, expected_capacity=8)
        self.assertEqual(result, {"capacity": 8, "used": 4, "key_bits": 29,
                                 "value_bits": 3, "taxonomy_nodes": 6,
                                 "hash_bytes": 64, "taxonomy_bytes": 377})

    def test_validation_reads_only_two_headers(self):
        original = Path.open
        reads = []

        class HeaderOnly:
            def __init__(self, handle):
                self.handle = handle

            def __enter__(self):
                return self

            def __exit__(self, *args):
                self.handle.close()

            def read(self, size=-1):
                if size != 32:
                    raise AssertionError(f"Unexpected payload read: {size}")
                reads.append(size)
                return self.handle.read(size)

            def fileno(self):
                return self.handle.fileno()

        def open_header(path, *args, **kwargs):
            return HeaderOnly(original(path, *args, **kwargs))

        with patch.object(Path, "open", open_header):
            reference._native_index_structure(self.db)
        self.assertEqual(reads, [32, 32])

    def test_truncated_headers_are_rejected(self):
        for path in (self.hash, self.taxonomy):
            original = path.read_bytes()
            with self.subTest(path=path.name):
                path.write_bytes(original[:31])
                with self.assertRaisesRegex(RuntimeError, "Truncated native index header"):
                    reference._native_index_structure(self.db)
                path.write_bytes(original)

    def test_truncated_or_extra_hash_payload_is_rejected(self):
        original = self.hash.read_bytes()
        for content in (original[:-1], original + b"x"):
            with self.subTest(length=len(content)):
                self.hash.write_bytes(content)
                with self.assertRaisesRegex(RuntimeError, "Native hash length mismatch"):
                    reference._native_index_structure(self.db)

    def test_invalid_hash_headers_are_rejected(self):
        for fields in ((0, 0, 29, 3), (8, 9, 29, 3), (8, 4, 30, 3),
                       (8, 4, 32, 0), (8, 4, 0, 32)):
            with self.subTest(fields=fields):
                self.hash.write_bytes(struct.pack("<QQQQ", *fields) + bytes(32))
                with self.assertRaisesRegex(RuntimeError, "Invalid native hash header"):
                    reference._native_index_structure(self.db)

    def test_expected_capacity_must_match(self):
        for expected in (9, 0, True, "8"):
            with self.subTest(expected=expected):
                with self.assertRaisesRegex(RuntimeError, "differs from the build estimate"):
                    reference._native_index_structure(self.db, expected_capacity=expected)

    def test_invalid_taxonomy_magic_is_rejected(self):
        self.taxonomy.write_bytes(b"badmagic" + self.taxonomy.read_bytes()[8:])
        with self.assertRaisesRegex(RuntimeError, "Invalid native taxonomy magic"):
            reference._native_index_structure(self.db)

    def test_invalid_taxonomy_counts_are_rejected(self):
        for fields in ((1, 4, 5), (6, 0, 5), (6, 4, 0)):
            with self.subTest(fields=fields):
                self.taxonomy.write_bytes(b"K2TAXDAT" + struct.pack("<QQQ", *fields)
                                          + bytes(345))
                with self.assertRaisesRegex(RuntimeError, "Invalid native taxonomy header"):
                    reference._native_index_structure(self.db)

    def test_truncated_or_extra_taxonomy_payload_is_rejected(self):
        original = self.taxonomy.read_bytes()
        for content in (original[:-1], original + b"x"):
            with self.subTest(length=len(content)):
                self.taxonomy.write_bytes(content)
                with self.assertRaisesRegex(RuntimeError, "Native taxonomy length mismatch"):
                    reference._native_index_structure(self.db)

    def test_hash_value_bits_must_fit_taxonomy(self):
        self.hash.write_bytes(struct.pack("<QQQQ", 8, 4, 30, 2) + bytes(32))
        with self.assertRaisesRegex(RuntimeError, "taxid bits cannot encode"):
            reference._native_index_structure(self.db)

    def test_symlinks_are_rejected(self):
        for path in (self.hash, self.taxonomy):
            target = path.with_suffix(".target")
            path.rename(target)
            path.symlink_to(target)
            with self.subTest(path=path.name):
                with self.assertRaisesRegex(RuntimeError, "regular non-symlink"):
                    reference._native_index_structure(self.db)
            path.unlink()
            target.rename(path)

    def test_snapshot_accepts_legacy_valid_receipt_without_rewriting_identity(self):
        for capacity in (8, None):
            expected = self._marker(capacity)
            before = (self.db / "database.json").read_bytes()
            with self.subTest(capacity=capacity):
                with patch.object(reference, "sha256", side_effect=AssertionError("full scan forbidden")):
                    self.assertEqual(reference.database_snapshot(self.cfg), expected)
                self.assertEqual((self.db / "database.json").read_bytes(), before)

    def test_snapshot_rejects_short_hash_even_with_matching_receipt_size_and_sha(self):
        self.hash.write_bytes(self.hash.read_bytes()[:-4])
        self._marker()
        with self.assertRaisesRegex(RuntimeError, "Native hash length mismatch"):
            reference.database_snapshot(self.cfg)

    def test_snapshot_rejects_truncated_taxonomy_with_matching_receipt(self):
        self.taxonomy.write_bytes(self.taxonomy.read_bytes()[:-4])
        self._marker()
        with self.assertRaisesRegex(RuntimeError, "Native taxonomy length mismatch"):
            reference.database_snapshot(self.cfg)

    def test_snapshot_rejects_capacity_mismatch_against_build_resource(self):
        self._marker(capacity=9)
        with self.assertRaisesRegex(RuntimeError, "differs from the build estimate"):
            reference.database_snapshot(self.cfg)


if __name__ == "__main__":
    unittest.main()
