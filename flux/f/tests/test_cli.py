"""Public CLI and scheduler tests: no UKB reads, GPUs, or real submissions."""
from __future__ import annotations

import builtins
from contextlib import redirect_stderr, redirect_stdout
import copy
import io
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import flux as cli
import submit
from common import ROOT, load_config


class CLITests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="flux cli ")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.output = self.directory / "analysis result"

    def run_main(self, arguments):
        with patch.object(sys, "argv", ["flux.py", *arguments]), redirect_stdout(io.StringIO()):
            return cli.main()

    def plan(self, stage="all", machine="local", device="cpu", extra=()):
        return self.run_main([stage, "--machine", machine, "--device", device,
                              "--plan-only", "--dir-out", str(self.output),
                              "--dir-ukb", str(self.directory / "UKB does not exist"), *extra])

    def snapshot(self):
        snapshots = list((self.output / "logs").glob("plan_*.yml"))
        self.assertEqual(len(snapshots), 1)
        return yaml.safe_load(snapshots[0].read_text())

    def test_no_arguments_submits_all_hpc_gpu_without_local_analysis(self):
        with patch.dict(os.environ, {}, clear=True), \
             patch.object(submit, "plan") as planner, \
             patch.object(cli, "execute", side_effect=AssertionError("Local analysis")):
            self.run_main([])
        cfg, args, machine = planner.call_args.args
        self.assertEqual((args.stage, machine, args.device), ("all", "hpc", "gpu"))
        self.assertFalse(args.plan_only)

    def test_plan_only_never_imports_analysis_or_invokes_subprocess(self):
        original_import = builtins.__import__
        prohibited = {"torch", "numpy", "pandas", "data", "train", "model", "evaluate", "rollout", "demo"}

        def import_guard(name, *args, **kwargs):
            if name.split(".")[0] in prohibited:
                raise AssertionError(f"Plan-only attempted heavy import: {name}")
            return original_import(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=import_guard), \
             patch.object(submit.subprocess, "run", side_effect=AssertionError("Plan submitted or ran a worker")):
            self.plan(machine="hpc", device="gpu", extra=("--gpus", "4"))
        self.assertEqual({p.name for p in self.output.iterdir()}, {"logs"})
        self.assertEqual({p.stem for p in (self.output / "logs").glob("*.cmd")}, set(cli.STAGES))
        self.assertEqual(self.snapshot()["resources"]["gpus"], 4)

    def test_public_shell_plan_handles_spaces_without_data_or_environment_install(self):
        environment = os.environ.copy()
        environment.update(FLUX_PYTHON=sys.executable, FLUX_ENV=str(self.directory / "missing environment"))
        result = subprocess.run(["bash", str(ROOT / "flux.sh"), "all", "--machine", "local",
                                 "--device", "cpu", "--plan-only", "--dir-out", str(self.output),
                                 "--dir-ukb", str(self.directory / "missing UKB")],
                                env=environment, text=True, capture_output=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse((self.directory / "missing environment").exists())
        self.assertEqual({p.name for p in self.output.iterdir()}, {"logs"})
        for command in (self.output / "logs").glob("*.cmd"):
            script = command.read_text()
            self.assertNotIn("#BSUB", script)
            tokens = shlex.split(script.splitlines()[-1])
            self.assertEqual(tokens[tokens.index("--machine") + 1], "local")
            self.assertEqual(tokens[tokens.index("--device") + 1], "cpu")
            snapshot = Path(tokens[tokens.index("--config") + 1])
            self.assertEqual(yaml.safe_load(snapshot.read_text())["paths"]["out"], str(self.output))

    def test_hpc_gpu_and_cpu_resources_match_stage(self):
        self.plan(machine="hpc", device="gpu", extra=("--gpus", "3", "--cores", "6"))
        train = (self.output / "logs/train.cmd").read_text()
        readout = (self.output / "logs/readout.cmd").read_text()
        prepare = (self.output / "logs/prepare.cmd").read_text()
        self.assertIn('#BSUB -gpu "num=3:mode=exclusive_process"', train)
        self.assertIn('#BSUB -gpu "num=1:mode=exclusive_process"', readout)
        self.assertNotIn("#BSUB -gpu", prepare)
        self.assertIn("OMP_NUM_THREADS=6 MKL_NUM_THREADS=6", train)
        self.assertIn('span[hosts=1] rusage[mem=128GB/host]', train)
        self.assertIn("--device cpu", prepare)
        self.assertIn("--device gpu", train)

    def test_demo_preserves_small_defaults_user_overrides_and_worker_output(self):
        custom = self.directory / "custom.yml"
        custom.write_text("train:\n  max_steps: 2\npaths:\n  out: ${FLUX_TEST_OUTPUT}\n")
        with patch.dict(os.environ, {"FLUX_TEST_OUTPUT": str(self.output)}):
            self.run_main(["demo", "--machine", "hpc", "--device", "gpu", "--gpus", "2",
                           "--plan-only", "--config", str(custom)])
        cfg = self.snapshot()
        self.assertEqual(cfg["model"]["deterministic_dim"], 24)
        self.assertEqual(cfg["train"]["epochs"], 3)
        self.assertEqual(cfg["train"]["max_steps"], 2)
        self.assertEqual(cfg["paths"]["out"], str(self.output))
        command = (self.output / "logs/demo.cmd").read_text()
        self.assertIn('#BSUB -gpu "num=2:mode=exclusive_process"', command)
        self.assertIn("--device gpu", command)
        planned = next((self.output / "logs").glob("plan_*.yml"))
        # Re-planning the worker's command must retain the same small demo and
        # output path; previously it silently fell back to ROOT/demo_out.
        self.run_main(["demo", "--worker", "--machine", "local", "--device", "gpu",
                       "--plan-only", "--config", str(planned)])
        snapshots = list((self.output / "logs").glob("plan_*.yml"))
        self.assertEqual(len(snapshots), 2)
        for snapshot in snapshots:
            resolved = yaml.safe_load(snapshot.read_text())
            self.assertEqual(resolved["paths"]["out"], str(self.output))
            self.assertEqual(resolved["model"]["deterministic_dim"], 24)

    def test_local_plan_does_not_apply_hpc_queue_policy(self):
        fake_root = self.directory / "repo/flux"
        helper = fake_root.parent / "0f/qiming_resources.py"
        helper.parent.mkdir(parents=True)
        helper.write_text('def cpu_slots(*args):\n    raise AssertionError("HPC policy called for local plan")\n')
        with patch.object(submit, "ROOT", fake_root):
            self.plan(stage="prepare", machine="local", extra=("--cores", "200"))
        self.assertNotIn("#BSUB", (self.output / "logs/prepare.cmd").read_text())

    def test_hpc_plan_uses_repository_cpu_slot_policy(self):
        fake_root = self.directory / "repo/flux"
        helper = fake_root.parent / "0f/qiming_resources.py"
        helper.parent.mkdir(parents=True)
        helper.write_text('def cpu_slots(queue, cores, memory, hours):\n    return "testqueue", 30\n')
        with patch.object(submit, "ROOT", fake_root):
            self.plan(stage="prepare", machine="hpc", extra=("--cores", "8"))
        command = (self.output / "logs/prepare.cmd").read_text()
        self.assertIn("#BSUB -q testqueue", command)
        self.assertIn("#BSUB -n 30", command)
        self.assertIn("OMP_NUM_THREADS=8", command)

    def test_invalid_resources_and_cpu_multigpu_fail_before_plan_creation(self):
        for options in (("--cores", "0"), ("--gpus", "0"), ("--gpus", "2")):
            with self.subTest(options=options), self.assertRaises(ValueError):
                self.plan(stage="train", machine="local", device="cpu", extra=options)
        for key, value in (("hours", "-.inf"), ("cpu_memory_gib", "-1"),
                           ("gpu_memory_gib", ".nan"), ("cores", "1.5")):
            custom = self.directory / "invalid.yml"
            custom.write_text(f"resources:\n  {key}: {value}\n")
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.plan(stage="prepare", extra=("--config", str(custom)))
        self.assertFalse(self.output.exists())

    def test_repeated_plans_keep_immutable_configuration_and_valid_minutes(self):
        custom = self.directory / "resource.yml"
        custom.write_text("resources:\n  hours: 1.999\n  gpus: 1.0\n")
        self.plan(machine="hpc", extra=("--config", str(custom), "--cores", "2"))
        first_path = next((self.output / "logs").glob("plan_*.yml"))
        first_content = first_path.read_bytes()
        self.plan(machine="hpc", extra=("--config", str(custom), "--cores", "3"))
        self.assertEqual(first_path.read_bytes(), first_content)
        self.assertIsInstance(yaml.safe_load(first_content)["resources"]["gpus"], int)
        self.assertEqual(len(list((self.output / "logs").glob("plan_*.yml"))), 2)
        command = (self.output / "logs/train.cmd").read_text()
        self.assertIn("#BSUB -W 2:00", command)
        self.assertNotIn("1:60", command)

    def test_submissions_form_stage_dependency_chain_without_executing_workers(self):
        cfg = load_config()
        cfg["paths"]["out"] = str(self.output)
        args = SimpleNamespace(stage="all", device="gpu", plan_only=False, resume=False, replace=False)
        submitted = []

        def bsub(command, **kwargs):
            submitted.append((command, kwargs["stdin"].read()))
            return SimpleNamespace(stdout=f"Job <{100 + len(submitted)}> is submitted to queue <test>.", stderr="")

        with patch.object(submit.shutil, "which", return_value="/mock/bsub"), \
             patch.object(submit.subprocess, "run", side_effect=bsub), redirect_stdout(io.StringIO()):
            submit.plan(cfg, args, "hpc")
        self.assertEqual(len(submitted), len(cli.STAGES))
        self.assertEqual(submitted[0][0], ["bsub"])
        for index, (command, _) in enumerate(submitted[1:], start=1):
            self.assertEqual(command, ["bsub", "-w", f"done({100 + index})"])
        self.assertFalse((self.output / "data").exists())
        self.assertEqual(len((self.output / "logs/jobs.tsv").read_text().splitlines()), 7)

    def test_unknown_options_and_invalid_enums_are_rejected(self):
        for arguments in (["train", "--unknown-option"], ["train", "--machine", "cloud"],
                          ["train", "--device", "cuda"], ["nonexistent-stage"]):
            with self.subTest(arguments=arguments), redirect_stderr(io.StringIO()), \
                 self.assertRaises(SystemExit) as caught:
                cli.parser().parse_args(arguments)
            self.assertEqual(caught.exception.code, 2)


class ShellRuntimeTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="flux runtime ")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.bin = self.directory / "bin"
        self.bin.mkdir()

    def environment(self):
        environment = os.environ.copy()
        for key in ("LSB_JOBID", "LSF_ENVDIR", "LSF_SERVERDIR", "DIR_WORK", "DIR_UKB",
                    "ANALYSIS_ROOT", "FLUX_ENV", "FLUX_MACHINE"):
            environment.pop(key, None)
        environment.update(FLUX_PYTHON=sys.executable, FLUX_ROOT=str(self.directory / "standalone/flux"),
                           FLUX_RUNTIME=str(ROOT / "f/runtime.sh"), PATH=str(self.bin))
        return environment

    def detect(self, machine="auto", markers=None):
        environment = self.environment()
        environment.update(markers or {})
        environment["FLUX_TEST_MACHINE"] = machine
        result = subprocess.run(["/bin/bash", "-c",
                                 'set -euo pipefail; source "$FLUX_RUNTIME"; '
                                 'flux_runtime --machine "$FLUX_TEST_MACHINE"; '
                                 'printf "%s\\n" "$FLUX_MACHINE" "$DIR_UKB" "$ANALYSIS_ROOT"'],
                                env=environment, text=True, capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout.splitlines()

    def test_auto_defaults_local_and_detects_lsf_without_submission(self):
        self.assertEqual(self.detect(), ["local", "/mnt/d/phe/ukb", "/mnt/d/analysis"])
        self.assertEqual(self.detect(markers={"LSB_JOBID": "123"})[0], "hpc")
        self.assertEqual(self.detect(markers={"LSF_ENVDIR": "/mock/lsf"})[0], "hpc")
        self.assertEqual(self.detect(markers={"LSF_SERVERDIR": "/mock/lsf/server"})[0], "hpc")
        bsub = self.bin / "bsub"
        bsub.write_text("#!/bin/sh\nexit 99\n")
        bsub.chmod(0o755)
        self.assertEqual(self.detect()[0], "hpc")

    def test_explicit_machine_overrides_auto_detection(self):
        self.assertEqual(self.detect("local", {"LSB_JOBID": "123"})[0], "local")
        self.assertEqual(self.detect("hpc"), ["hpc", "/work/sph-huangj/phe/ukb", "/work/sph-huangj/analysis"])

    def test_installer_check_never_creates_environment_or_installs(self):
        capture = self.directory / "check arguments.txt"
        python = self.bin / "python mock"
        python.write_text('#!/bin/bash\nprintf "%s\\n" "$@" > "$FLUX_TEST_CAPTURE"\n')
        python.chmod(0o755)
        environment = os.environ.copy()
        environment.update(FLUX_PYTHON=str(python), FLUX_ENV=str(self.directory / "environment"),
                           FLUX_TEST_CAPTURE=str(capture), DIR_WORK=str(self.directory / "work"))
        result = subprocess.run(["bash", str(ROOT / "install.sh"), "--machine", "local", "--check"],
                                env=environment, text=True, capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(capture.read_text().strip(), str(ROOT / "f/check_install.py"))
        self.assertFalse((self.directory / "environment").exists())
        self.assertFalse((self.directory / "work").exists())

    def test_installer_rejects_unknown_options_without_installation(self):
        result = subprocess.run(["bash", str(ROOT / "install.sh"), "--unknown-option"],
                                text=True, capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 2)
        self.assertIn("Unknown option", result.stderr)


if __name__ == "__main__":
    unittest.main()

