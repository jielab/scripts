"""Small field-projection regressions; no participant data or network required."""
from pathlib import Path
import gzip
import subprocess
import sys
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "ukb/f/date_check_io.py"


class FieldProjectionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="phe-fields-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.raw = self.root / "raw.tsv"
        self.raw.write_text(
            "eid\tp23460_i0\tp23461_i0\tp23460_i1\tp53_i2\tp53_i3\n"
            "1001\t8\t2\t10\t2018-01-01\t2022-01-01\n"
        )
        self.mapping = self.root / "met.lst"
        self.out = self.root / "out.tsv.gz"

    def run_mapping(self, text, expected=0):
        self.mapping.write_text(text)
        p = subprocess.run(
            [sys.executable, str(SCRIPT), "fields", "--input", str(self.raw),
             "--mapping", str(self.mapping), "--output", str(self.out)],
            text=True, capture_output=True,
        )
        self.assertEqual(p.returncode, expected, p.stdout + p.stderr)
        if not expected:
            with gzip.open(self.out, "rt") as f:
                return f.read().splitlines()
        return p

    def test_met_header_ratio_and_baseline(self):
        lines = self.run_mapping(
            "data_field\tmet_name\tfull_name\n"
            "p23460\tAla\tAlanine\n"
            "p23461\tGln\tGlutamine\n"
            "p23460/p23461\tAla2Gln\tAlanine to glutamine\n"
        )
        self.assertEqual(lines, ["eid\tp23460_i0\tp23461_i0", "1001\t8\t2"])

    def test_ratio_only_extracts_both_operands_once(self):
        lines = self.run_mapping("p23460 / p23461\tratio\n")
        self.assertEqual(lines[0], "eid\tp23460_i0\tp23461_i0")

    def test_all_instances_unchanged(self):
        lines = self.run_mapping("p53\tdate_attend\tall\n")
        self.assertEqual(lines[0], "eid\tp53_i2\tp53_i3")

    def test_baseline_and_single_instance_mapping_compatibility(self):
        a = self.run_mapping("p23460\tAla\ti0\n")
        b = self.run_mapping("p23460_i0\tAla\ti0\n")
        self.assertEqual(a, b)
        self.assertEqual(a[0], "eid\tp23460_i0")

    def test_invalid_mapping_keeps_existing_output(self):
        self.run_mapping("p23460\tAla\ti0\n")
        before = self.out.read_bytes()
        for expr in ("p23460+1", "p23460/p23461/p53", "p23460_i0/p23461", "typo"):
            with self.subTest(expression=expr):
                self.run_mapping(expr + "\tbad\n", expected=2)
                self.assertEqual(self.out.read_bytes(), before)
                self.assertFalse(list(self.root.glob(".*.tmp")))

    def test_malformed_input_keeps_existing_output(self):
        self.run_mapping("p23460\tAla\n")
        before = self.out.read_bytes()
        self.raw.write_text(self.raw.read_text() + "1002\t8\n")
        self.run_mapping("p23460\tAla\n", expected=2)
        self.assertEqual(self.out.read_bytes(), before)
        self.assertFalse(list(self.root.glob(".*.tmp")))


if __name__ == "__main__":
    unittest.main()
