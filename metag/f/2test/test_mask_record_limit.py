"""Reject overlong single records without confusing record and file lengths."""
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import hpc_reference as ref


class MaskRecordLimitTests(unittest.TestCase):
    def emit(self, text, masking=True):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "ref.fna"
            source.write_text(text)
            cfg = {"library_dir": str(root / "lib"), "metadata_dir": str(root / "meta"),
                   "reference_raw": str(root), "masking": masking}
            rows = [{"path": "ref.fna", "index": 1, "taxid": 10, "state": {"size": source.stat().st_size}}]
            stats = {}
            with patch.object(ref, "MAX_MASK_RECORD", 30):
                ref._emit_fasta(cfg, rows, io.BytesIO(), io.BytesIO(), stats)
            return stats

    def test_limit_is_per_record_and_total_file_can_be_larger(self):
        stats = self.emit(">a\n" + "A" * 20 + "\n>b\n" + "C" * 20 + "\n")
        self.assertEqual((stats["bases"], stats["longest_record"]), (40, 20))

    def test_masked_record_limit_is_checked_across_wrapped_lines(self):
        with self.assertRaisesRegex(RuntimeError, "signed-position limit"):
            self.emit(">a\n" + "A" * 20 + "\n" + "C" * 20 + "\n")

    def test_unmasked_record_does_not_inherit_masker_integer_limit(self):
        self.assertEqual(self.emit(">a\n" + "A" * 40 + "\n", masking=False)["bases"], 40)
