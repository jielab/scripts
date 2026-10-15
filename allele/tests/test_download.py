"""Downloader command/manifest regressions; no network requests or large files."""

import json
import os
from pathlib import Path
import subprocess
import shutil
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class DownloadTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.script = ROOT / "gnomad_download.sh"
        self.capture = self.base / "captured.json"
        stub = self.base / "aria2c"
        stub.write_text(
            "#!/usr/bin/env python3\n"
            "import json, os, sys\n"
            "from pathlib import Path\n"
            "args=sys.argv[1:]\n"
            "source=next(x.split('=',1)[1] for x in args if x.startswith('--input-file='))\n"
            "Path(os.environ['ALLELE_DOWNLOAD_CAPTURE']).write_text(json.dumps("
            "{'args':args, 'source':source, 'manifest':Path(source).read_text()}))\n"
            "sys.exit(int(os.environ.get('ALLELE_DOWNLOAD_TEST_EXIT','0')))\n"
        )
        stub.chmod(0o755)
        self.env = {
            **os.environ,
            "PATH": str(self.base) + os.pathsep + os.environ["PATH"],
            "ALLELE_CONFIG": "/dev/null",
            "GNOMAD_OUT": str(self.base / "downloads"),
            "GNOMAD_PROXY": "http://127.0.0.1:7897",
            "ALLELE_DOWNLOAD_CAPTURE": str(self.capture),
        }

    def run_download(self, *args):
        return subprocess.run(
            [str(self.script), *args],
            env=self.env,
            cwd=self.base,
            text=True,
            capture_output=True,
        )

    def read_capture(self):
        data = json.loads(self.capture.read_text())
        self.assertFalse(Path(data["source"]).exists(), "temporary manifest leaked")
        return data

    def test_all_selects_48_files_with_checksums_and_resume(self):
        result = self.run_download("--chr", "all")
        self.assertEqual(result.returncode, 0, result.stderr)
        data = self.read_capture()
        self.assertEqual(
            data["manifest"], (ROOT / "gnomad.genomes.v4.1.aria2.txt").read_text()
        )
        lines = data["manifest"].splitlines()
        self.assertEqual(len(lines), 48 * 3)
        self.assertEqual(len({line for line in lines if line.startswith("  out=")}), 48)
        self.assertIn("chr22.vcf.bgz", lines[1])
        for i in range(0, len(lines), 3):
            self.assertEqual(len(lines[i].split("\t")), 2)
            self.assertRegex(lines[i + 2], r"^  checksum=md5=[0-9a-f]{32}$")
        for flag in [
            "-c",
            "--always-resume=true",
            "--allow-overwrite=false",
            "--check-integrity=true",
        ]:
            self.assertIn(flag, data["args"])
        self.assertIn("--https-proxy=http://127.0.0.1:7897", data["args"])
        self.assertIn("--no-proxy=", data["args"])

    def test_chromosome_selection_includes_only_requested_pairs(self):
        result = self.run_download(
            "--chr", "1,22,1", "--proxy", "http://localhost:8888"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        data = self.read_capture()
        names = [
            line.strip().removeprefix("out=")
            for line in data["manifest"].splitlines()
            if line.startswith("  out=")
        ]
        expected = {
            f"gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz{suffix}"
            for chrom in (1, 22)
            for suffix in ("", ".tbi")
        }
        self.assertEqual(set(names), expected)
        self.assertEqual(len(names), 4)
        self.assertIn("--https-proxy=http://localhost:8888", data["args"])

    def test_direct_mode_overrides_inherited_proxies(self):
        self.env["https_proxy"] = "http://localhost:9999"
        result = self.run_download("--chr", "22", "--proxy", "none")
        self.assertEqual(result.returncode, 0, result.stderr)
        args = self.read_capture()["args"]
        for flag in ["--https-proxy=", "--http-proxy=", "--all-proxy="]:
            self.assertIn(flag, args)

    def test_invalid_input_never_starts_downloader(self):
        for args in [
            ("--chr", "23"),
            ("--chr", "22,"),
            ("--chr", "1,,2"),
            ("--chr", "22", "--proxy"),
            ("--chr", "22", "--proxy", "invalid"),
            ("--chr", "22", "--proxy", "http://localhost:65536"),
            ("--chr", "22", "--jobs", "0"),
            ("--chr", "22", "--out"),
            ("--proxy", "none"),
        ]:
            with self.subTest(args=args):
                result = self.run_download(*args)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(self.capture.exists())

    def test_downloader_failure_is_preserved_and_manifest_removed(self):
        self.env["ALLELE_DOWNLOAD_TEST_EXIT"] = "7"
        result = self.run_download("--chr", "22")
        self.assertEqual(result.returncode, 7)
        self.read_capture()

    def test_copied_script_runs_without_project_or_manifest(self):
        remote = self.base / "another computer"
        remote.mkdir()
        self.script = remote / "gnomad_download.sh"
        shutil.copy2(ROOT / "gnomad_download.sh", self.script)
        output = remote / "download files"
        result = self.run_download(
            "--proxy",
            "none",
            "--chr",
            "1,2,3,4,5",
            "--out",
            str(output),
            "--jobs",
            "1",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        data = self.read_capture()
        self.assertEqual(len(data["manifest"].splitlines()), 10 * 3)
        self.assertIn("--dir=" + str(output), data["args"])
        self.assertIn("--max-concurrent-downloads=1", data["args"])
        self.assertTrue(output.is_dir())

    def test_installer_no_longer_accepts_download_command(self):
        result = subprocess.run(
            [str(ROOT / "install.sh"), "--download-gnomad", "22"],
            env=self.env,
            text=True,
            capture_output=True,
        )
        self.assertEqual(result.returncode, 2)
        self.assertFalse(self.capture.exists())


if __name__ == "__main__":
    unittest.main()
