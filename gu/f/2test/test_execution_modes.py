"""Exercise public wrappers and saved workers without launching inference."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


GU = Path(__file__).resolve().parents[2]


class ExecutionModes(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="gu-mode-qa-", dir=GU.parent)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        excluded = ("GU_", "IBDMIX_", "ARG_", "TRACE_", "PHYML_", "DIR_", "OS_", "LSB_", "WORKFLOW_", "GRID_")
        self.env = {k: v for k, v in os.environ.items()
                    if not k.startswith(excluded) and k not in {"submit", "SUBMIT", "ANALYSIS_ROOT", "SCRIPTS_ROOT", "DATE_SCRATCH", "jobs_N", "SHARED_PYTHON"}}
        self.env.update(GU_SUBMIT_PYTHON=sys.executable, GU_HPC_PYTHON=sys.executable,
                        SCRIPT_CONSOLE_ACTIVE="1", PYTHONDONTWRITEBYTECODE="1")
        bindir = self.root / "bin"
        bindir.mkdir()
        bsub = bindir / "bsub"
        bsub.write_text("#!/bin/bash\necho 'test must not submit jobs' >&2\nexit 99\n")
        bsub.chmod(0o755)
        self.env["PATH"] = str(bindir) + os.pathsep + self.env["PATH"]
        rscript = bindir / "Rscript"
        rscript.write_text("#!/bin/bash\nexit 0\n")
        rscript.chmod(0o755)
        self.env["GU_RESULTS_RSCRIPT"] = str(rscript)
        self.env["GU_REF_ROOT"] = str(self.root / "genotypes")


    def run_shell(self, entry, args, env=None):
        return subprocess.run(["bash", str(GU / entry), *map(str, args)],
                              env=env or self.env, text=True, capture_output=True, timeout=30)

    def ibdmix_command(self, mode, env=None):
        selected = env or self.env
        reference = Path(selected["GU_REF_ROOT"]) / "1kg/37/pfile"
        reference.mkdir(parents=True, exist_ok=True)
        for suffix in (".pgen", ".pvar", ".psam"):
            (reference / ("chrX" + suffix)).write_text("fixture\n")
        out = self.root / mode / "cohort"
        result = self.run_shell("gu.sh", ["ibdmix", "--submit", mode, "--chr", "X",
                               "--target", "1kg", "--dir-out", out, "--run-cmd", "FALSE"], env)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        commands = list(out.parent.glob("cohort.jobs/*/*.cmd"))
        self.assertEqual(len(commands), 1)
        return commands[0]

    def test_explicit_local_wins_when_bsub_is_on_path(self):
        command = self.ibdmix_command("local")
        text = command.read_text()
        self.assertIn("export submit=local SUBMIT=local", text)
        self.assertNotIn('${LSB_JOBID:', text)
        self.assertIn("--memory-cap 32G", text)

    def test_hpc_worker_requires_allocation(self):
        command = self.ibdmix_command("hpc")
        text = command.read_text()
        self.assertIn("export submit=hpc SUBMIT=hpc", text)
        self.assertIn('${LSB_JOBID:', text)
        result = subprocess.run(["bash", str(command)], env=self.env, text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Submit this .cmd with bsub", result.stderr)

    def test_saved_worker_keeps_explicit_paths_in_a_fresh_shell(self):
        conda = self.root / "custom-conda"
        conda_source = conda / "etc/profile.d/conda.sh"
        conda_source.parent.mkdir(parents=True)
        conda_source.write_text('conda() { if [[ $2 == /* ]]; then export CONDA_PREFIX="$2"; else export CONDA_PREFIX="$GU_CONDA_BASE/envs/$2"; fi; }\n')
        env = {**self.env, "submit": "hpc", "SUBMIT": "hpc",
               "DIR_WORK": str(self.root / "work"), "DIR_DATA": str(self.root / "data"),
               "DIR_SCRATCH": str(self.root / "scratch"), "DIR_SOFTWARE": str(self.root / "software"),
               "GU_CONDA_BASE": str(conda), "GU_REF_ROOT": str(self.root / "custom-genotypes")}
        command = self.ibdmix_command("hpc", env)
        # Execute only the actual generated environment prefix, never the final
        # inference command. The fake conda initializer does no installations.
        prefix = "\n".join(line for line in command.read_text().splitlines() if not line.startswith("exec "))
        prefix += '\nprintf "%s\\n" "$DIR_DATA" "$DIR_SCRATCH" "$GU_REF_ROOT" "$GU_CONDA_BASE"\n'
        clean = {"PATH": self.env["PATH"], "HOME": self.env["HOME"], "LSB_JOBID": "123"}
        result = subprocess.run(["bash", "-c", prefix], env=clean, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(result.stdout.splitlines()[-4:],
                         [env[k] for k in ("DIR_DATA", "DIR_SCRATCH", "GU_REF_ROOT", "GU_CONDA_BASE")])

    def test_direct_hpc_inference_stops_before_runtime_or_payload(self):
        for entry, args in (("gu.sh", ["trace"]), ("gu.sh", ["phyml"]),
                            ("arg.sh", ["build"]), ("arg.sh", ["prep_gen"])):
            with self.subTest(entry=entry, args=args):
                result = self.run_shell(entry, [*args, "--submit", "hpc"])
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("allocated LSF job", result.stderr)
                self.assertIn("hpc plan", result.stderr)
                self.assertNotIn("per-process virtual-memory", result.stderr)

    def test_hpc_help_is_available_without_allocation(self):
        for entry, args in (("gu.sh", []), ("arg.sh", ["build"])):
            with self.subTest(entry=entry):
                result = self.run_shell(entry, [*args, "--submit", "hpc", "--help"])
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_existing_advanced_hpc_entry_plans_all_non_ibdmix_methods(self):
        gen = self.root / "genotypes"
        gen.mkdir()
        for extension, content in (("pgen", "fixture"), ("pvar", "#CHROM\tPOS\tID\tREF\tALT\n"),
                                   ("psam", "#IID SEX\na 1\nb 2\n")):
            (gen / ("chr1." + extension)).write_text(content)
        runtime = self.root / "runtime.env"
        runtime.write_text("# Plan-only fixture\n")
        loci = self.root / "loci.jma.cojo"
        loci.write_text("Chr SNP bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r\n")
        for method in ("arg", "arg-trace", "trace", "phyml"):
            with self.subTest(method=method):
                out = self.root / ("plan-" + method)
                args = ["hpc", "plan", method, "--submit", "hpc", "--queue", "ser",
                        "--run-dir", out, "--env-file", runtime, "--dir-target", gen / "chr", "--chr", "1"]
                if method in {"arg", "arg-trace"}: args += ["--dir-gen", gen]
                if method == "trace": args += ["--dir-arg", self.root / "trees"]
                if method == "phyml": args += ["--loci", loci, "--archaic-path", self.root / "archaic"]
                result = self.run_shell("gu.sh", args)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                plan = json.loads((out / "plan.json").read_text())
                self.assertEqual(plan["method"], method)
                self.assertEqual(len(plan["tasks"]), 1)
                self.assertFalse((out / "state/submission.jsonl").exists())


if __name__ == "__main__":
    unittest.main()
