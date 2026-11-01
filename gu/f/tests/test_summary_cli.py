"""Regression tests for denominators, zero calls, and public command planning.

CLI planning tests mock only the expensive IBDmix mask preflight; they never
execute inference and are not a substitute for data/hardware validation.
"""
from __future__ import annotations
import importlib.util
import os
from pathlib import Path
import sqlite3
import shlex
import subprocess
import sys
import tempfile
import unittest

GU = Path(__file__).resolve().parents[2]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


norm = load("summary_test_normalize", GU / "f/results_normalize.py")


class SummaryTests(unittest.TestCase):
    def test_union_excludes_untested_and_retains_tested_zero(self):
        con = sqlite3.connect(":memory:")
        con.execute("CREATE TABLE segments(dataset_id,genome_build,method,source,sample_id,chr,start,end)")
        con.executemany("INSERT INTO segments VALUES('d','GRCh37','ibdmix','Altai',?,?,?,?)", [
            ("positive", "1", 10, 30), ("positive", "1", 20, 50),
            ("positive", "2", 0, 100), ("untested", "1", 0, 100),
        ])
        rows = norm.reference_summary(con, "d", "GRCh37", ["positive", "zero", "untested"],
                                      {"1": 100, "2": 100}, {"1": {"positive", "zero"}}, "Altai")
        self.assertEqual([r["archaic_bp"] for r in rows], [40, 0, ""])
        self.assertEqual(rows[0]["coverage_pct"], 20)
        self.assertEqual(rows[0]["coordinate_coverage_pct"], 40)
        tagged = [dict(r, reference="Altai") for r in rows]
        aggregate = next(r for r in norm.population_burdens(tagged, {}, "d", "GRCh37", "autosomes",
                                                         {"positive", "zero", "untested"}, []) if r["group"] == "ALL")
        self.assertEqual((aggregate["n_tested"], aggregate["n_zero"], aggregate["n_untested"]), (2, 1, 1))
        self.assertEqual(aggregate["complete_chromosomes"], 0)
        self.assertEqual(aggregate["mean_mb"], 20 / 1e6)
        con.close()

    def test_x_is_clipped_to_nonpar_physical_span(self):
        con = sqlite3.connect(":memory:")
        con.execute("CREATE TABLE segments(dataset_id,genome_build,method,source,sample_id,chr,start,end)")
        length = 155270560
        con.execute("INSERT INTO segments VALUES('d','GRCh37','ibdmix','Altai','male','X',-100,?)", (length + 100,))
        row = norm.reference_summary(con, "d", "GRCh37", ["male"], {"X": length}, {"X": {"male"}}, "Altai", True)[0]
        span = sum(b-a for a,b in norm.x_nonpar_intervals("GRCh37", length))
        self.assertEqual(row["archaic_bp"], span)
        self.assertEqual(row["tested_bp"], span)
        self.assertEqual(row["coverage_pct"], 100)
        self.assertEqual(row["coordinate_coverage_pct"], 100)
        con.close()

    def test_same_males_require_complete_autosomes_and_matching_reference(self):
        def row(sample, reference, n, bp, span, copies):
            return dict(sample_id=sample, reference=reference, n_chromosomes=n,
                        archaic_bp=bp, haploid_bp=span, tested_bp=span*copies,
                        coverage_pct=100*bp/(span*copies))
        autos = [row("m1", "Altai", 22, 200, 10000, 2), row("m2", "Altai", 21, 100, 9000, 2),
                 row("m3", "Denisova", 22, 100, 10000, 2)]
        xs = [row(s, "Altai", 1, 2, 1000, 1) for s in ("m1", "m2", "m3")]
        result = next(r for r in norm.paired_male_burdens(xs, autos, {}, "d", "GRCh37") if r["group"] == "ALL")
        self.assertEqual(result["n_same_males"], 1)
        self.assertAlmostEqual(result["x_to_autosome_reported_coverage_ratio"], .2)
        self.assertAlmostEqual(result["x_to_autosome_coordinate_coverage_ratio"], .1)

class PublicPlanningTests(unittest.TestCase):
    def test_chromosome_ranges_and_invalid_input(self):
        shell = 'source "$1"; gu_normalize_chrs "$2"'
        p = subprocess.run(["bash", "-c", shell, "test", str(GU / "f/0.common.sh"), "1-22,X,chr1,23"], capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(p.stdout.strip().split(), [*map(str, range(1, 23)), "X"])
        for bad in ("22-1", "0-3", "1-24", "X-Y", "foo"):
            p = subprocess.run(["bash", "-c", shell, "test", str(GU / "f/0.common.sh"), bad], capture_output=True, text=True)
            self.assertNotEqual(p.returncode, 0, bad)

    def plan(self, tmp, method, gpu=None, refs=None):
        tmp = Path(tmp)
        bindir = tmp / "bin"
        bindir.mkdir()
        wrapper = bindir / "python3"
        wrapper.write_text(f"#!{sys.executable}\nimport os,sys\nif len(sys.argv)>2 and sys.argv[1].endswith('ibdmix.py') and sys.argv[2]=='check-masks': sys.exit(0)\nos.execv({sys.executable!r}, [{sys.executable!r},*sys.argv[1:]])\n")
        wrapper.chmod(0o755)
        argdir = tmp / "arg"
        argdir.mkdir()
        for chrom in ("1", "2"):
            (argdir / f"chr{chrom}.trees").write_text("planning fixture only\n")
        env = dict(os.environ, GU_ENV_REEXEC="1", SCRIPT_CONSOLE_ACTIVE="1", GU_KEEP_SCRATCH="1",
                   PATH=str(bindir)+os.pathsep+os.environ["PATH"], GU_SOFT=str(tmp / "soft"))
        for key in ("GU_CMD_WORKER", "GU_RESULT_MODULE", "GU_ANALYSIS_ROOT", "GU_PUBLISHED_ROOT", "GU_FILTERED_TARGET", "GU_UKB_DIRECT", "GU_UKB_PGEN_PREPARED", "GU_CHRS", "GU_KEEP", "GU_MALE_ONLY"):
            env.pop(key, None)
        args = ["bash", str(GU / "gu.sh"), method, "--grch", "37", "--target", "1kg",
                "--chr", "1,2", "--dir-target", str(tmp / "gen/chr"), "--dir-out", str(tmp / "out"),
                "--run-cmd", "FALSE", "--auto-final", "FALSE", "--foreground", "TRUE"]
        if method.startswith("ibdmix"):
            args += ["--archaic-path", str(tmp / "archaic")]
            if refs is not None:
                args += ["--archaic-gen", refs]
        else:
            args += ["--dir-arg", str(argdir)]
        if gpu is not None:
            args += ["--GPU", gpu]
        p = subprocess.run(args, env=env, capture_output=True, text=True, timeout=30)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        return sorted((tmp / "out").rglob("*.cmd"))

    def test_ibdmix_backend_is_explicit_and_layout_unified(self):
        for method, gpu in [("ibdmix", "FALSE"), ("ibdmix", "TRUE")]:
            with self.subTest(method=method, gpu=gpu), tempfile.TemporaryDirectory(dir=os.environ.get("GU_TEST_OUTPUT_ROOT", GU.parent)) as tmp:
                cmds = self.plan(tmp, method, gpu)
                self.assertEqual(len(cmds), 2)
                for cmd in cmds:
                    self.assertEqual(cmd.relative_to(tmp).parts[1], "ibdmix")
                    content = cmd.read_text()
                    self.assertRegex(content, r"--GPU\s+" + gpu)
                    self.assertNotIn(".2013", content)
                    argv = shlex.split(next(x for x in content.splitlines() if x.startswith("exec ")))
                    self.assertEqual(argv[argv.index("--archaic-gen") + 1], "Altai,Chagyr,Denisova,Denisova25,Vindija")
                    self.assertIn("--dir-target", argv)

    def test_optional_reference_requires_explicit_selection(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get("GU_TEST_OUTPUT_ROOT", GU.parent)) as tmp:
            cmds = self.plan(tmp, "ibdmix", "FALSE", refs="Altai.2013,Denisova.2013")
            for cmd in cmds:
                argv = shlex.split(next(x for x in cmd.read_text().splitlines() if x.startswith("exec ")))
                self.assertEqual(argv[argv.index("--archaic-gen") + 1], "Altai.2013,Denisova.2013")

    def test_trace_plans_one_joint_request(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get("GU_TEST_OUTPUT_ROOT", GU.parent)) as tmp:
            cmds = self.plan(tmp, "trace")
            self.assertEqual(len(cmds), 1)
            self.assertIn("--chr 1\\ 2", cmds[0].read_text())

    def test_single_ibdmix_module_rejects_other_method_names(self):
        p = subprocess.run(["bash", str(GU / "gu.sh"), "ibdmix" + "-GPU", "--GPU", "TRUE"],
                           capture_output=True, text=True, timeout=10)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("unknown method", p.stderr)


if __name__ == "__main__":
    unittest.main()
