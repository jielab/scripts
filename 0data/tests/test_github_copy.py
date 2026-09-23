"""Synthetic end-to-end publishing regressions; never reads real UKB data."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile


SCRIPT = Path(__file__).resolve().parents[1] / "github_copy.sh"
PREFIX = "panome/cvd_cad/prot/v4_reference/"


class PublishingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="github-copy-test-")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / "source"
        self.dest = self.root / "dest"
        (self.source / "panome").mkdir(parents=True)
        self.dest.mkdir()
        self.env = dict(os.environ, GITHUB_COPY_SOURCE_ROOT=str(self.source),
                        GITHUB_COPY_DEST_ROOT=str(self.dest))

    def write(self, name, content):
        path = self.source / PREFIX / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def run_copy(self, mode, ok=True):
        args = ["bash", str(SCRIPT), mode]
        if mode == "scan":
            args.append("panome")
        result = subprocess.run(args, env=self.env, capture_output=True, text=True)
        if ok:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        # Scanner diagnostics must never echo individual values.
        self.assertNotIn("9000001", result.stdout + result.stderr)
        return result

    def manifest_entries(self):
        return {line for line in (self.dest / "github_files.lst").read_text().splitlines()
                if line and not line.startswith("#")}

    def sentinel(self):
        target = self.dest / "panome" / "existing.txt"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("keep until all checks pass")
        return target

    def python_wrapper(self, body):
        directory = self.root / "bin"
        directory.mkdir()
        wrapper = directory / "python3"
        wrapper.write_text("#!/usr/bin/env bash\n" + body)
        wrapper.chmod(0o755)
        self.env["PATH"] = str(directory) + os.pathsep + self.env["PATH"]
        self.env["REAL_PYTHON"] = sys.executable

    def test_compound_identifiers_and_single_person_are_excluded(self):
        headers = ["eid_A", "eid_B", "reference_eid", "case_eid", "control_eid",
                   "EID.A", "referenceEid", "participantID", "sample_id",
                   "person_id", "IID", "FID", "ID", "id_1", "ID_A", "donor_id"]
        for i, header in enumerate(headers):
            self.write(f"result_{i}.csv", f"{header},risk\n9000001,0.2\n")
        self.write("same_risk_pairs.csv", "eid_A,eid_B,predicted_risk_A,predicted_risk_B,profile_distance,largest_contrasts\n9000001,9000002,0.2,0.2,1.5,synthetic\n")
        self.write("reference_utilization.csv", "reference_eid,test_matches,local_net_risk\n9000001,2,0.1\n")
        self.write("summary.csv", "model,n,AUC\nclinical,100,0.7\n")
        self.run_copy("scan")
        self.assertEqual(self.manifest_entries(), {PREFIX + "summary.csv"})

    def test_identifier_free_outputs_and_configuration_are_allowed(self):
        files = {
            "test_predictions.csv": "time,event,AE1,AE2,state_weight_1,risk\n2,1,0.2,0.5,0.8,0.1\n",
            "reference_bank.csv": "feature_id,gene_id,mean\np1,12345,0.3\n",
            "manifest.json": json.dumps({"id_col": "eid", "features": ["age", "sex"]}),
            "empty.csv": "eid_A,risk\n,0.2\nNA,0.3\n",
            "features.txt": "id1\nid2\nid3\neid1\neid2\n",
            "metric_features.csv": "feature,mean_abs_beta\nid1,0.1\nid2,0.2\n",
            "molecular_blocks.csv": "module,feature,loading\nM1,id1,0.1\nM2,id2,0.2\n",
        }
        for name, content in files.items():
            self.write(name, content)
        self.run_copy("scan")
        self.assertEqual(self.manifest_entries(), {PREFIX + name for name in files})
        self.run_copy("sync")
        for name, content in files.items():
            self.assertEqual((self.dest / PREFIX / name).read_text(), content)

    def test_nested_and_column_oriented_json(self):
        fixtures = {
            "nested.json": {"records": [{"reference_eid": 9000001, "risk": 0.2}]},
            "columns.json": {"eid_A": [9000001], "risk": [0.2]},
            "mapping.json": {"9000001": {"risk": 0.2}},
            "rows.json": [["eid_A", "risk"], [9000001, 0.2]],
        }
        for name, value in fixtures.items():
            self.write(name, json.dumps(value))
        self.write("records.jsonl", '{"id_col":"eid"}\n{"reference_eid":9000001}\n')
        self.run_copy("scan")
        self.assertEqual(self.manifest_entries(), set())

    def test_all_rows_and_text_formats(self):
        self.write("late.csv", "metric,value\n" + "auc,0.7\n" * 12000 + "eid_A,risk\n9000001,0.2\n")
        self.write("tab.tsv", "reference_eid\trisk\n9000001\t0.2\n")
        self.write("report.html", "<table><tr><th>eid_A</th><th>risk</th></tr><tr><td>9000001</td><td>0.2</td></tr></table>")
        self.write("data.yaml", "reference_eid: 9000001\nrisk: 0.2\n")
        self.write("report.md", "| eid_A | risk |\n| --- | --- |\n| 9000001 | 0.2 |\n")
        self.write("figure.svg", '<svg><text>reference_eid: 9000001</text></svg>')
        self.run_copy("scan")
        self.assertEqual(self.manifest_entries(), set())

    def test_hidden_xlsx_sheet_and_far_column(self):
        path = self.write("report.xlsx", "")
        ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
        with zipfile.ZipFile(path, "w") as z:
            z.writestr("xl/workbook.xml", f'<workbook xmlns="{ns}"><sheets><sheet name="Data" state="hidden"/></sheets></workbook>')
            z.writestr("xl/sharedStrings.xml", f'<sst xmlns="{ns}"><si><t>reference_eid</t></si></sst>')
            z.writestr("xl/worksheets/sheet2.xml", f'<worksheet xmlns="{ns}"><sheetData><row r="2"><c r="AA2" t="s"><v>0</v></c></row><row r="3"><c r="AA3"><v>9000001</v></c></row></sheetData></worksheet>')
        self.run_copy("scan")
        self.assertEqual(self.manifest_entries(), set())

    def test_uninspectable_formats_fail_closed(self):
        self.write("bad.xlsx", "not a zip")
        path = self.write("broken_string.xlsx", "")
        with zipfile.ZipFile(path, "w") as z:
            z.writestr("xl/worksheets/sheet1.xml", '<worksheet><sheetData><row><c r="A1" t="s"><v>0</v></c></row></sheetData></worksheet>')
        self.write("bad.json", "{")
        path = self.write("bad.csv", "")
        path.write_bytes(b"\xff\x80invalid")
        self.run_copy("scan")
        self.assertEqual(self.manifest_entries(), set())

    def test_sync_rejects_content_changed_after_scan_without_touching_destination(self):
        path = self.write("summary.csv", "model,AUC\nclinical,0.7\n")
        self.run_copy("scan")
        keep = self.sentinel()
        path.write_text("reference_eid,risk\n9000001,0.2\n")
        result = self.run_copy("sync", ok=False)
        self.assertIn("UKB SYNC BLOCKED", result.stderr)
        self.assertTrue(keep.exists())
        self.assertFalse((self.dest / PREFIX / "summary.csv").exists())

    def test_manually_added_sensitive_manifest_entry_cannot_bypass_sync(self):
        self.write("summary.csv", "model,AUC\nclinical,0.7\n")
        self.run_copy("scan")
        self.write("same_risk_pairs.csv", "eid_A,risk\n9000001,0.2\n")
        with (self.dest / "github_files.lst").open("a") as manifest:
            manifest.write(PREFIX + "same_risk_pairs.csv\n")
        keep = self.sentinel()
        self.run_copy("sync", ok=False)
        self.assertTrue(keep.exists())
        self.assertFalse((self.dest / PREFIX / "same_risk_pairs.csv").exists())

    def test_sync_publishes_only_the_inspected_snapshot(self):
        content = "model,AUC\nclinical,0.7\n"
        path = self.write("summary.csv", content)
        self.run_copy("scan")
        self.env["MUTATE_SOURCE"] = str(path)
        self.env["STAGE_RECORD"] = str(self.root / "stage-path")
        self.python_wrapper('''"$REAL_PYTHON" "$@"
result=$?
if [[ ${2-} == /tmp/github_copy.stage.* ]]; then
    printf '%s' "$2" > "$STAGE_RECORD"
    printf 'reference_eid,risk\\n9000001,0.2\\n' > "$MUTATE_SOURCE"
fi
exit "$result"
''')
        self.run_copy("sync")
        self.assertEqual((self.dest / PREFIX / "summary.csv").read_text(), content)
        self.assertIn("reference_eid", path.read_text())
        self.assertFalse(Path((self.root / "stage-path").read_text()).exists())

    def test_scanner_failure_preserves_manifest_and_destination(self):
        self.write("summary.csv", "model,AUC\nclinical,0.7\n")
        self.run_copy("scan")
        before = (self.dest / "github_files.lst").read_bytes()
        keep = self.sentinel()
        self.python_wrapper("exit 23\n")
        self.run_copy("scan", ok=False)
        self.assertEqual((self.dest / "github_files.lst").read_bytes(), before)
        self.run_copy("sync", ok=False)
        self.assertTrue(keep.exists())

    def test_sync_removes_stale_files_only_in_selected_project(self):
        self.write("summary.csv", "model,AUC\nclinical,0.7\n")
        keep = self.sentinel()
        unrelated = self.dest / "grid" / "keep.txt"
        unrelated.parent.mkdir()
        unrelated.write_text("untouched")
        self.run_copy("scan")
        self.run_copy("sync")
        self.assertFalse(keep.exists())
        self.assertEqual(unrelated.read_text(), "untouched")


if __name__ == "__main__":
    unittest.main()
