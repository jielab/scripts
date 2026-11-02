"""Query scratch must coexist with, but never overlap, protected data."""
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hpc_query import _roots


class QueryRootTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        base = self.root / "scratch" / "analysis" / "metag"
        self.cfg = {
            "work_dir": str(base),
            "output_dir": str(base / "split_hpc"),
            "metadata_dir": str(base / "metadata_hpc"),
            "library_dir": str(base / "lib_hpc"),
            "reads_raw": str(self.root / "scratch" / "metag" / "raw"),
            "reference_raw": str(self.root / "data" / "annot" / "metag" / "raw"),
            "tools_dir": str(self.root / "data" / "tools" / "metag"),
        }

    def test_default_shared_base_allows_distinct_query_subtree(self):
        out, work = _roots(self.cfg)
        self.assertEqual(out, Path(self.cfg["output_dir"]))
        self.assertEqual(work, Path(self.cfg["work_dir"]))
        self.assertTrue(out.is_dir())
        self.assertTrue(work.is_dir())

    def test_query_scratch_cannot_overlap_any_protected_root(self):
        query_root = Path(self.cfg["work_dir"]) / "queries"
        for key in ("reads_raw", "reference_raw", "metadata_dir", "library_dir", "tools_dir"):
            for protected in (query_root, query_root / "asset", query_root.parent):
                with self.subTest(protected_key=key, protected=str(protected)):
                    cfg = dict(self.cfg, **{key: str(protected)})
                    with self.assertRaisesRegex(RuntimeError, "overlaps protected"):
                        _roots(cfg)

    def test_output_cannot_overlap_any_protected_root(self):
        output = Path(self.cfg["output_dir"])
        for key in ("reads_raw", "reference_raw", "metadata_dir", "library_dir", "tools_dir"):
            for protected in (output, output / "asset", output.parent):
                with self.subTest(protected_key=key, protected=str(protected)):
                    cfg = dict(self.cfg, **{key: str(protected)})
                    with self.assertRaisesRegex(RuntimeError, "overlaps protected"):
                        _roots(cfg)

    def test_output_and_query_scratch_cannot_overlap(self):
        query_root = Path(self.cfg["work_dir"]) / "queries"
        for output in (query_root, query_root / "outputs", query_root.parent):
            with self.subTest(output=str(output)):
                with self.assertRaisesRegex(RuntimeError, "output and scratch roots overlap"):
                    _roots(dict(self.cfg, output_dir=str(output)))

    def test_query_scratch_symlink_is_rejected(self):
        work = Path(self.cfg["work_dir"])
        work.mkdir(parents=True)
        target = self.root / "external"
        target.mkdir()
        (work / "queries").symlink_to(target, target_is_directory=True)
        with self.assertRaisesRegex(RuntimeError, "Symlink/alias in query work root"):
            _roots(self.cfg)
        self.assertFalse(list(target.iterdir()))


if __name__ == "__main__":
    unittest.main()
