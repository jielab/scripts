"""Exercise the real console subprocess path used by research workflows."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
CONSOLE = ROOT / "runtime" / "pipeline_console.py"
RESEARCH = ROOT / "f" / "research.py"
LABELS = (
    "research-balance", "research-space", "research-domains",
    "research-cohort", "research-archaic-alleles",
    "research-data", "research-allofus", "research-mfi", "research-panel", "research-mapped-targets", "pca-trajectory",
)


class TestResearchLogging(unittest.TestCase):
    def invoke_console(self, log_dir, label, command):
        result = subprocess.run(
            [sys.executable, str(CONSOLE), "--log-dir", str(log_dir),
             "--label", label, "--interval", "300", "--", *command],
            capture_output=True, text=True, timeout=30,
            env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"),
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def test_all_research_labels_preserve_raw_stdout_and_stderr(self):
        with tempfile.TemporaryDirectory(prefix="allele3-research-logs-") as tmp:
            log_dir = Path(tmp) / "logs"
            payload = b"raw stdout\rprogress\nraw stderr\n"
            command = [sys.executable, "-c",
                       "import sys; "
                       "sys.stdout.buffer.write(b'raw stdout\\rprogress\\n'); "
                       "sys.stderr.buffer.write(b'raw stderr\\n')"]
            for label in LABELS:
                with self.subTest(label=label):
                    result = self.invoke_console(log_dir, label, command)
                    logs = list(log_dir.glob(label + "_*.log"))
                    self.assertEqual(len(logs), 1)
                    self.assertEqual(logs[0].read_bytes(), payload)
                    self.assertIn("full log=" + str(logs[0]), result.stdout)

    def test_balance_transaction_keeps_sibling_logs_at_reported_paths(self):
        with tempfile.TemporaryDirectory(prefix="allele3-balance-log-") as tmp:
            root = Path(tmp)
            source, out = root / "genotype_validation.tsv", root / "balance"
            # Logs must remain outside the transactionally replaced result.
            log_dir = root / "balance.logs"
            rows = []
            for pos, counts in ((1, (90, 5, 5)), (2, (0, 60, 40))):
                rows.append(dict(
                    site_id=f"1:{pos}:A", chrom="1", pos=pos, ref="A",
                    family="ERV1", group="case", population="EUR",
                    n_selected=50, n_called=50, AN=100, call_rate=1,
                    A_AC=counts[0], C_AC=counts[1], G_AC=counts[2], T_AC="NA",
                    encoded_alts="C,G", evaluable=1,
                    represented_status="all_expected_alts_encoded", record_status="ok",
                ))
            with source.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            command = [sys.executable, str(RESEARCH), "balance",
                       "--genotypes", str(source), "--out", str(out)]
            saved_logs = set()
            # A refresh moves the previous result to history. Both live log
            # paths must survive these real successful balance transactions.
            for attempt in range(2):
                result = self.invoke_console(log_dir, "research-balance", command)
                manifest = json.loads((out / "balance_manifest.json").read_text())
                self.assertEqual(manifest["status"], "complete")
                current_logs = set(log_dir.glob("research-balance_*.log"))
                self.assertEqual(len(current_logs), attempt + 1)
                self.assertTrue(saved_logs <= current_logs)
                for log in current_logs - saved_logs:
                    self.assertIn("full log=" + str(log), result.stdout)
                    self.assertTrue(log.is_file())
                saved_logs = current_logs
            with (out / "allele_balance.tsv").open() as handle:
                values = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([r["n_observed_bases"] for r in values], ["3", "2"])
            self.assertEqual([r["at_least_three_supported"] for r in values], ["1", "0"])
            self.assertFalse(list((root / ".balance.history").rglob("*.log")))


if __name__ == "__main__":
    unittest.main()
