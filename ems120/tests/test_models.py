"""Runtime-policy tests without importing Torch, patient data, or pretrained models.

Extract the production helper definitions with AST so these checks run on a
login node and in CI even when the scientific Conda environment is absent.
"""
import ast
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
import types
import unittest
from unittest import mock
import uuid


SOURCE = Path(__file__).resolve().parents[1] / "f" / "ems120.py"
TREE = ast.parse(SOURCE.read_text(encoding="utf-8"))


def helpers(*names, **globals_):
    selected = [node for node in TREE.body if isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name in names]
    assert len(selected) == len(names), names
    env = dict(os=os, json=json, Path=Path, re=re, shutil=shutil, tempfile=tempfile, uuid=uuid,
               ems_log=lambda message: None)
    env.update(globals_)
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(SOURCE), "exec"), env)
    return env


def fake_xgboost(device="cuda:0", version="2.1.0", compiled=True, fail=False):
    class Regressor:
        def __init__(self, **params):
            self.params = params

        def fit(self, x, y):
            if fail:
                raise RuntimeError("synthetic CUDA driver failure")
            return self

        def get_booster(self):
            return types.SimpleNamespace(save_config=lambda: json.dumps({"learner": {"generic_param": {"device": device}}}))

    return types.SimpleNamespace(__version__=version, build_info=lambda: {"USE_CUDA": compiled}, XGBRegressor=Regressor)


class GPUPolicyTests(unittest.TestCase):
    def env(self, policy="required", cuda=True):
        return helpers("_xgboost_is_cuda", "check_phone_gpu", PHONE_GPU=policy, CPU_CORES=4,
                       _XGBOOST_GPU_CHECKED=False,
                       torch=types.SimpleNamespace(cuda=types.SimpleNamespace(is_available=lambda: cuda)),
                       np=types.SimpleNamespace(float32=float, asarray=lambda value, dtype: value))

    def test_required_checks_actual_booster_and_does_not_accept_fallback(self):
        env = self.env()
        with mock.patch.dict(sys.modules, {"xgboost": fake_xgboost(device="cpu")}):
            with self.assertRaisesRegex(RuntimeError, "fell back to CPU"):
                env["check_phone_gpu"]()
        self.assertFalse(env["_XGBOOST_GPU_CHECKED"])

    def test_required_rejects_old_cpu_only_or_broken_cuda(self):
        for module, message in [(fake_xgboost(version="1.7.6"), ">= 2.0"),
                                (fake_xgboost(compiled=False), "no CUDA support"),
                                (fake_xgboost(fail=True), "driver failure")]:
            with self.subTest(message=message), mock.patch.dict(sys.modules, {"xgboost": module}):
                with self.assertRaisesRegex(RuntimeError, message):
                    self.env()["check_phone_gpu"]()

    def test_cuda_success_is_reused(self):
        env = self.env()
        module = fake_xgboost()
        with mock.patch.dict(sys.modules, {"xgboost": module}):
            self.assertTrue(env["check_phone_gpu"]())
            module.XGBRegressor = lambda **kwargs: self.fail("GPU probe repeated")
            self.assertTrue(env["check_phone_gpu"]())

    def test_explicit_auto_and_off_policies(self):
        with mock.patch.dict(sys.modules, {"xgboost": fake_xgboost(device="cpu")}):
            self.assertFalse(self.env("auto")["check_phone_gpu"]())
        with mock.patch.dict(sys.modules, {"xgboost": None}):
            self.assertFalse(self.env("off")["check_phone_gpu"]())

    def test_required_never_uses_gpu_when_torch_cannot_see_allocation(self):
        with mock.patch.dict(sys.modules, {"xgboost": fake_xgboost()}):
            with self.assertRaisesRegex(RuntimeError, "CUDA is unavailable"):
                self.env(cuda=False)["check_phone_gpu"]()


class CacheTests(unittest.TestCase):
    def test_training_failure_reaches_caller_instead_of_nan_results(self):
        def fail_training(*args):
            raise RuntimeError("training failed")
        env = helpers("ensure_phone_sco3_model", "predict_phone_sco3",
                      _batch_size=lambda kind, value: value, _is_digit8=lambda value: True,
                      PHONE_SCO3_MODEL_FILE="/nonexistent/ems120/model.pkl",
                      PHONE_SCO3_FEATURE_FILE="/nonexistent/ems120/features.pkl",
                      PHONE_SCO3_REPORT_FILE="/nonexistent/ems120/report.csv", TRAIN_FILE="training.xlsx",
                      train_phone_sco3_ml=fail_training)
        with self.assertRaisesRegex(RuntimeError, "training failed"):
            env["predict_phone_sco3"](["12345678"])

    def test_interrupted_phone_pair_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            model = Path(directory) / "model.pkl"
            Path(str(model) + ".incomplete").write_text("{}")
            env = helpers("_load_phone_model", _PHONE_CACHE=None, _PHONE_CACHE_KEY=None)
            with self.assertRaisesRegex(RuntimeError, "interrupted"):
                env["_load_phone_model"](str(model), str(Path(directory) / "features.pkl"))

    def test_loaded_cpu_winner_is_bounded_and_cached_without_changing_method(self):
        class CPUWinner:
            def __init__(self):
                self.params = {"n_jobs": -1, "model__n_jobs": -1, "n_estimators": 500}
            def get_params(self, deep=True):
                return dict(self.params)
            def set_params(self, **kwargs):
                self.params.update(kwargs)
                return self
        winner = CPUWinner()
        loads = []
        def load(path):
            loads.append(path)
            return winner if str(path).endswith("model.pkl") else ["feature"]
        env = helpers("_bound_phone_model", "_load_phone_model", _PHONE_CACHE=None, _PHONE_CACHE_KEY=None,
                      CPU_CORES=6, PHONE_GPU="required", joblib=types.SimpleNamespace(load=load))
        with tempfile.TemporaryDirectory() as directory:
            model, feature = Path(directory) / "model.pkl", Path(directory) / "features.pkl"
            model.write_text("model")
            feature.write_text("features")
            self.assertIs(env["_load_phone_model"](model, feature)[0], winner)
            self.assertIs(env["_load_phone_model"](model, feature)[0], winner)
        self.assertEqual(len(loads), 2)
        self.assertEqual(winner.params, {"n_jobs": 6, "model__n_jobs": 6, "n_estimators": 500})

    def test_cached_xgboost_is_retargeted_without_retraining(self):
        class CachedXGB:
            __module__ = "xgboost.sklearn"
            def __init__(self):
                self.params = {"device": "cuda:4", "n_jobs": -1, "n_estimators": 600}
            def get_params(self, deep=True):
                return dict(self.params)
            def set_params(self, **kwargs):
                self.params.update(kwargs)
                return self
        env = helpers("_bound_phone_model", CPU_CORES=8, check_phone_gpu=lambda: True)
        model = CachedXGB()
        self.assertIs(env["_bound_phone_model"](model), model)
        self.assertEqual(model.params, {"device": "cuda", "n_jobs": 8, "n_estimators": 600})
        env["check_phone_gpu"] = lambda: False
        with mock.patch.dict(sys.modules, {"xgboost": fake_xgboost()}):
            env["_bound_phone_model"](model)
        self.assertEqual(model.params["device"], "cpu")

    def test_atomic_file_failure_keeps_previous_cache(self):
        env = helpers("_atomic_write")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cache.json"
            path.write_text("complete")
            def broken_writer(temporary):
                Path(temporary).write_text("partial")
                raise RuntimeError("interrupted write")
            with self.assertRaisesRegex(RuntimeError, "interrupted write"):
                env["_atomic_write"](path, broken_writer)
            self.assertEqual(path.read_text(), "complete")
            self.assertEqual(list(Path(directory).iterdir()), [path])

    def test_transformer_failed_publication_restores_previous_model(self):
        env = helpers("_atomic_write", "_atomic_json", "_publish_transformer")
        class Bundle:
            def save_pretrained(self, directory):
                Path(directory, "weights").write_text("new")
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "hfl"
            target.mkdir()
            (target / "weights").write_text("old")
            original_replace = os.replace
            def fail_publication(source, destination):
                if Path(destination) == target and ".staging-" in Path(source).name:
                    raise OSError("failed model publication")
                return original_replace(source, destination)
            with mock.patch("os.replace", side_effect=fail_publication):
                with self.assertRaisesRegex(OSError, "failed model publication"):
                    env["_publish_transformer"](Bundle(), Bundle(), target, {"labels.json": {"a": 1}})
            self.assertEqual((target / "weights").read_text(), "old")
            self.assertEqual(list(Path(directory).iterdir()), [target])

    def test_transformer_publication_contains_all_metadata(self):
        env = helpers("_atomic_write", "_atomic_json", "_publish_transformer", "_transformer_files_ready")
        class Bundle:
            def save_pretrained(self, directory):
                for name in ("model.safetensors", "config.json", "tokenizer_config.json"):
                    Path(directory, name).write_text("complete")
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "geo"
            env["_publish_transformer"](Bundle(), Bundle(), target, {"labels.json": {"a": 1}})
            self.assertTrue(env["_transformer_files_ready"](target))
            self.assertEqual(json.loads((target / "labels.json").read_text()), {"a": 1})

    def test_batch_override_is_validated_without_changing_default(self):
        env = helpers("_positive_int", "_batch_size")
        with mock.patch.dict(os.environ, {}, clear=True):
            self.assertEqual(env["_batch_size"]("DX", 1024), 1024)
            os.environ["EMS120_DX_BATCH_SIZE"] = "256"
            self.assertEqual(env["_batch_size"]("DX", 1024), 256)
            for invalid in ("0", "-1", "3.5", "nan"):
                os.environ["EMS120_DX_BATCH_SIZE"] = invalid
                with self.assertRaises(ValueError):
                    env["_batch_size"]("DX", 1024)


if __name__ == "__main__":
    unittest.main()
