"""Portable CLI/LSF integration regressions; no UKB data or real scheduler needed."""
from __future__ import annotations
import gzip
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
import dispatch


class PortableRuntime(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="gwas-portable-")
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.env = {k: v for k, v in os.environ.items() if not k.startswith(("GWAS_", "LSB_", "WORKFLOW_", "CONDA_"))
                    and k not in {"submit", "SUBMIT", "DIR_WORK", "DIR_DATA", "DIR_SCRATCH", "DIR_SOFTWARE", "SCRIPT_CONSOLE_ACTIVE", "SCRIPT_LOG_DIR", "CHECK_GRCH_SNP_LIST", "PGS_PLINK_MEMORY_MB", "jobs_N"}}
        self.env.update(DIR_WORK=str(self.work / "work"), DIR_DATA=str(self.work / "data"),
                        DIR_SCRATCH=str(self.work / "scratch"), SHARED_PYTHON=sys.executable,
                        WORKFLOW_PYTHON=sys.executable, GWAS_PLAN_PYTHON=sys.executable,
                        GWAS_PLAN_DIR=str(self.work / "plans"), SCRIPT_CONSOLE_PASSTHROUGH="1",
                        PYTHONDONTWRITEBYTECODE="1")

    def call(self, entry, *args, env=None, success=True):
        result = subprocess.run(["bash", str(ROOT / entry), *map(str, args)], cwd=self.work,
                                env=env or self.env, capture_output=True, text=True, timeout=30)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def plan(self, entry, module, *args, env=None):
        result = self.call(entry, module, "--submit", "hpc", "--plan-only", *args, env=env)
        match = re.search(r"^\[PLAN\] (.+) \(1 jobs; hpc\)$", result.stdout, re.M)
        self.assertIsNotNone(match, result.stdout)
        path = Path(match.group(1))
        return path, json.loads(path.read_text())

    def tool(self, name, body):
        path = self.work / "bin" / name
        path.parent.mkdir(exist_ok=True)
        path.write_text("#!" + sys.executable + "\n" + body)
        path.chmod(0o755)
        self.env["PATH"] = str(path.parent) + os.pathsep + self.env["PATH"]
        return path

    def genotype_fixture(self):
        geno = self.work / "genotypes with spaces"
        geno.mkdir()
        for chrom in ["1", "22"]:
            for ext in ["pgen", "pvar", "psam"]:
                (geno / f"chr{chrom}.{ext}").write_text("fixture\n")
        plink = self.tool("plink-fixture", '''import pathlib,sys
a=sys.argv[1:]
out=pathlib.Path(a[a.index('--out')+1])
with pathlib.Path(str(out)+'.calls').open('a') as h: h.write('call\\n')
pathlib.Path(str(out)+'.acount').write_text('#CHROM ID REF ALT ALT_CTS OBS_CT\\n1 rsKeep A G 20 100\\n1 rsDrop A C 0 100\\n')
''')
        return geno, plink

    def test_help_relocatable(self):
        for entry in ["format.sh", "gwas.sh", "compare.sh"]:
            with self.subTest(entry=entry):
                result = self.call(entry, "--help")
                self.assertIn("--submit auto|local|hpc", result.stdout)

    def test_hpc_planning_never_reads_raw_fifo(self):
        raw = self.work / "must-not-open.tsv.gz"
        os.mkfifo(raw)
        path, plan = self.plan("format.sh", "format", "--raw-file", raw, "--jobs", "3", "--tool-threads", "4")
        job = plan["jobs"][0]
        self.assertEqual(job["resources"]["cores"], 12)
        self.assertEqual(job["resources"]["queue"], "33")
        self.assertEqual(job["env"]["GWAS_TOOL_THREADS"], "4")
        self.assertIn(str(raw), job["argv"])
        self.assertEqual(job["argv"][-2:], ["--foreground", "TRUE"])
        cmd = path.parent / "commands/format.cmd"
        self.assertEqual(subprocess.run(["bash", "-n", str(cmd)]).returncode, 0)
        blocked = subprocess.run(["bash", str(cmd)], env=self.env, capture_output=True, text=True)
        self.assertNotEqual(blocked.returncode, 0)
        self.assertIn("login node", blocked.stderr)

    def test_hpc_paths_and_legacy_root(self):
        env = dict(self.env)
        for k in ["DIR_WORK", "DIR_DATA", "DIR_SCRATCH"]:
            env.pop(k)
        _, plan = self.plan("format.sh", "magma", env=env)
        settings = plan["jobs"][0]["env"]
        self.assertEqual(settings["GWAS_DATA_ROOT"], "/work/sph-huangj/data")
        self.assertEqual(settings["GWAS_FILES_ROOT"], "/work/sph-huangj/data/files")
        _, plan = self.plan("format.sh", "format", "--dir-data-root", str(self.work / "legacy"))
        self.assertEqual(plan["jobs"][0]["env"]["GWAS_FILES_ROOT"], str(self.work / "legacy/files"))
        self.assertEqual(plan["jobs"][0]["env"]["GWAS_SOFTWARE_ROOT"], str(self.work / "legacy/software"))

    def test_genotype_resource_budget(self):
        _, plan = self.plan("gwas.sh", "regenie", "--jobs", "3", "--threads", "2", "--step2-threads", "4", "--queue", "ser")
        r = plan["jobs"][0]["resources"]
        self.assertEqual(r["cores"], 12)
        self.assertGreaterEqual(r["memory_gib"], 3 * 24000 / 1024)
        self.call("gwas.sh", "regenie", "--submit", "hpc", "--plan-only", "--cores", "2", success=False)
        self.call("gwas.sh", "regenie", "--submit", "hpc", "--plan-only", "--memory-gb", "16", success=False)
        self.call("format.sh", "pgs", "--submit", "hpc", "--plan-only", "--memory-gb", "8", success=False)

    def test_invalid_combination_rejected_before_submission(self):
        self.call("gwas.sh", "regenie,acount", "--submit", "hpc", "--plan-only", success=False)

    def test_default_eight_workers_use_larger_queue(self):
        _, plan = self.plan("gwas.sh", "regenie")
        r = plan['jobs'][0]['resources']
        self.assertEqual(r['cores'], 32)
        self.assertEqual(r['memory_gib'], 192)
        self.assertEqual(r['queue'], '38')
        result = self.call("gwas.sh", "regenie", "--submit", "hpc", "--plan-only", "--queue", "33", success=False)
        self.assertIn('single-node budget', result.stderr)

    def test_auto_queue_covers_large_worker_core_count(self):
        _, plan = self.plan("gwas.sh", "regenie", "--jobs", "12", "--threads", "8")
        r = plan['jobs'][0]['resources']
        self.assertEqual(r['cores'], 96)
        self.assertEqual(r['queue'], '1t81c')

    def test_magma_combination_preserves_execution_default(self):
        from unittest.mock import patch
        with patch.dict(os.environ, {**self.env, "submit": "hpc"}, clear=True):
            self.assertFalse(dispatch.build("format.sh", ["magma,lead"])[1])
            self.assertTrue(dispatch.build("format.sh", ["format"])[1])
            self.assertFalse(dispatch.build("format.sh", ["--grch", "38", "magma"])[1])
            self.assertTrue(dispatch.build("format.sh", ["--label", "magma", "format"])[1])

    def test_existing_allocation_is_checked(self):
        env = dict(self.env, LSB_JOBID="123", LSB_DJOB_NUMPROC="2", LSB_QUEUE="ser")
        result = self.call("gwas.sh", "regenie", "--submit", "hpc", env=env, success=False)
        self.assertIn("Insufficient allocated CPU", result.stderr)

    def test_local_acount_plan_execution_and_reuse(self):
        geno, plink = self.genotype_fixture()
        args = ["acount", "--submit", "local", "--dir-imputed", geno, "--chr", "1,22", "--plink2", plink, "--jobs", "2"]
        self.call("gwas.sh", *args, "--run-cmd", "FALSE")
        self.assertFalse(list(geno.glob("*.calls")))
        self.call("gwas.sh", *args, "--run-cmd", "TRUE")
        for chrom in ["1", "22"]:
            self.assertEqual((geno / f"chr{chrom}.extract").read_text(), "rsKeep\n")
        self.call("gwas.sh", *args, "--run-cmd", "TRUE")
        self.assertTrue(all(p.read_text() == "call\n" for p in geno.glob("*.calls")))

    def test_hpc_generated_worker_executes_inside_mock_allocation(self):
        geno, plink = self.genotype_fixture()
        path, plan = self.plan("gwas.sh", "acount", "--dir-imputed", geno, "--chr", "1,22", "--plink2", plink, "--jobs", "2")
        cmd = path.parent / "commands/acount.cmd"
        result = subprocess.run(["bash", str(cmd)], env=dict(self.env, LSB_JOBID="456", LSB_DJOB_NUMPROC="40"),
                                text=True, capture_output=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        state = json.loads((path.parent / "status/acount.json").read_text())
        self.assertEqual(state["status"], "DONE")
        self.assertEqual((geno / "chr22.extract").read_text(), "rsKeep\n")

    def test_failed_submission_receipt_preserved(self):
        self.tool("bsub", "import sys\nsys.stdin.read()\nprint('fixture submission rejected',file=sys.stderr)\nsys.exit(7)\n")
        result = self.call("format.sh", "format", "--submit", "hpc", "--run-cmd", "TRUE", success=False)
        self.assertIn("submission", result.stderr)
        receipts = list((self.work / "plans").glob("*/submission.json"))
        if not receipts:
            receipts = list((self.work / "plans").glob("*/submissions.json"))
        self.assertTrue(receipts, result.stdout + result.stderr)
        self.assertFalse(list((self.work / "plans").glob("*/status/*.json")))

    def test_runtime_environment_applies_locally(self):
        marker = self.work / "activated"
        runtime = self.work / "runtime.sh"
        runtime.write_text("printf 'yes' > " + repr(str(marker)) + "\n")
        geno, plink = self.genotype_fixture()
        self.call("gwas.sh", "acount", "--submit", "local", "--runtime-env", runtime,
                  "--dir-imputed", geno, "--chr", "1", "--plink2", plink, "--run-cmd", "FALSE")
        self.assertEqual(marker.read_text(), "yes")

    @unittest.skipUnless(importlib.util.find_spec("pysam"), "pysam required for real BGZF/TBI fixture")
    def test_real_format_bgzf_index_and_resume(self):
        # These adapters call real HTSlib via pysam, not placeholder gzip/index files.
        self.tool("bgzip", '''import sys,gzip,pysam,tempfile,pathlib
a=sys.argv[1:]
if '-t' in a:
    with gzip.open(a[-1],'rb') as f: f.read()
else:
    with tempfile.TemporaryDirectory() as d:
        p=pathlib.Path(d)/'in'; p.write_bytes(sys.stdin.buffer.read())
        pysam.tabix_compress(str(p),str(p)+'.gz',force=True)
        sys.stdout.buffer.write(pathlib.Path(str(p)+'.gz').read_bytes())
''')
        self.tool("tabix", '''import sys,pysam
a=sys.argv[1:]
if '-l' in a:
    with pysam.TabixFile(a[-1]) as f: print('\\n'.join(f.contigs))
else:
    def opt(flag,default): return int(a[a.index(flag)+1]) if flag in a else default
    pysam.tabix_index(a[-1],seq_col=opt('-s',1)-1,start_col=opt('-b',2)-1,end_col=opt('-e',2)-1,
                     line_skip=opt('-S',0),zerobased=False,force='-f' in a,csi='-C' in a)
''')
        raw = self.work / "tiny.tsv.gz"
        with gzip.open(raw, "wt") as f:
            f.write("SNP\tCHR\tPOS\tEA\tNEA\tEAF\tBETA\tSE\tP\tN\n"
                    "rs2\t2\t200\tG\tA\t0.2\t0.1\t0.05\t0.04\t1000\n"
                    "rs1\t1\t100\tT\tC\t0.4\t-0.2\t0.1\t0.04\t1000\n")
        output = self.work / "formatted"
        args = ["format", "--submit", "local", "--raw-file", raw, "--gwas", "tiny",
                "--dir-out", output, "--grch", "38", "--hm3", "FALSE", "--thin", "FALSE",
                "--rsid", "FALSE", "--run-cmd", "TRUE", "--jobs", "1"]
        self.call("format.sh", *args)
        final = output / "common/tiny/gwas/tiny.gz"
        with gzip.open(final, "rt") as f:
            lines = f.read().splitlines()
        self.assertEqual(len(lines), 3)
        self.assertEqual([x.split("\t")[0] for x in lines[1:]], ["rs1", "rs2"])
        import pysam
        with pysam.TabixFile(str(final)) as tb:
            rows = list(tb.fetch("1", 99, 100))
        self.assertEqual(len(rows), 1)
        self.assertIn("rs1", rows[0])
        old = final.stat().st_mtime_ns
        result = self.call("format.sh", *args)
        self.assertEqual(final.stat().st_mtime_ns, old)
        self.assertIn("SKIP", result.stdout + result.stderr)
        hpc_output = self.work / "hpc-formatted"
        path, _ = self.plan("format.sh", "format", "--raw-file", raw, "--gwas", "tiny",
                            "--dir-out", hpc_output, "--grch", "38", "--hm3", "FALSE",
                            "--thin", "FALSE", "--rsid", "FALSE", "--jobs", "1")
        result = subprocess.run(["bash", str(path.parent / "commands/format.cmd")],
                                env=dict(self.env, LSB_JOBID="789", LSB_DJOB_NUMPROC="40"),
                                text=True, capture_output=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        with gzip.open(hpc_output / "common/tiny/gwas/tiny.gz", "rt") as f:
            self.assertEqual(f.read().splitlines(), lines)


if __name__ == "__main__":
    unittest.main()
