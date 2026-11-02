"""A moved raw cache retains its explicit analysis/kinship identity."""
import importlib.util
import json
from pathlib import Path
import shutil

import pandas as pd
import pytest

import kinship
import prsformer
from prsformer_benchmark import publish
from shared_cohort import sha256


def test_custom_cache_publication_verifies_truth_hashes_and_analysis_namespace(tmp_path, monkeypatch):
    # Reuse the existing evaluated-file fixture, not a trained GPU model. This
    # exercises the actual publication verifier and its strict identity gates.
    source = Path(__file__).parent / "abm/test_prsformer_independent.py"
    spec = importlib.util.spec_from_file_location("prsformer_receipt_fixture", source)
    fixture = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fixture)
    root = tmp_path / "analysis"
    cache = fixture.cache_fixture(root)
    external = tmp_path / "external_cache"
    shutil.move(cache, external)
    run = root / "prsformer/height/run"
    fixture.evaluated_run(run, external, "height")
    observed_roots = []
    def active(path):
        observed_roots.append(Path(path))
        return None
    monkeypatch.setattr(kinship, "active", active)
    receipt = publish(root, "height", cache=external)
    assert observed_roots and set(observed_roots) == {root}
    assert receipt["schema"] == "grid.target_benchmark.1"
    assert receipt["method_label"] == "PRSformer independent-trait adaptation"
    assert all(sha256(path) == digest for path, digest in receipt["files"].items())
    truth = pd.read_csv(receipt["test_truth_file"], sep="\t")
    pred = pd.read_csv(receipt["predictions_file"], sep="\t")
    assert set(pred.id) == set(truth.id) and "y" not in pred
    assert pred.role.eq("test").all()
    assert receipt["test_labels_used_only_after_prediction_freeze"]
    monkeypatch.setattr(kinship, "active", lambda path: {"signature": "newly-activated-family-correction"})
    with pytest.raises(ValueError, match="predates the active kinship"):
        prsformer.validate_evaluation(run, external, analysis_root=root)
