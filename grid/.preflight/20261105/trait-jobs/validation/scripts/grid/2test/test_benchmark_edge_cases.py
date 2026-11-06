"""Regression tests for unused-role scaling and producer-cohort resolution."""
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
sys.path.insert(0, str(Path(__file__).parent))

from target_benchmark import cohort_source, fit_score_model, main as target_main, standardize_scores
from csx_blocks import main as block_main
from test_block_producer import fixture as block_fixture


def test_singleton_unused_roles_do_not_change_valid_fit_or_test_scaling():
    rng = np.random.default_rng(411)
    scores = rng.normal(size=(62, 4))
    side = rng.normal(size=(62, 2))
    roles = np.array(["target_validation"] * 30 + ["test"] * 30 +
                     ["excluded", "development_unused"])
    yfit = scores[:30] @ np.array([2., -1., .5, .3]) + side[:30].sum(1)
    all_rows = fit_score_model(scores, side, yfit, roles, False, "per_split")
    relevant_rows = fit_score_model(scores[:60], side[:60], yfit, roles[:60], False, "per_split")
    assert np.isfinite(all_rows[1]).all()
    for full, relevant in zip(all_rows[:3], relevant_rows[:3]):
        np.testing.assert_allclose(full[:60], relevant, rtol=0, atol=0)
    assert all_rows[3] == relevant_rows[3]
    scaled, record = standardize_scores(scores, roles, "per_split")
    assert set(record["statistics"]) == {"target_validation", "test"}
    expected_unused = (scores[60:] - scores[:30].mean(0)) / scores[:30].std(0)
    np.testing.assert_allclose(scaled[60:], expected_unused)
    np.testing.assert_allclose(scaled[30:60].mean(0), 0, atol=1e-14)
    np.testing.assert_allclose(scaled[30:60].std(0), 1, atol=1e-14)


def test_targetfit_and_blockfit_use_bound_custom_producer_cohort(tmp_path, monkeypatch):
    config_path, config = block_fixture(tmp_path, "ALL.37.bed", n=240)
    source = Path(config["prsformer_cache"]) / "data.tsv.gz"
    frame = pd.read_csv(source, sep="\t", dtype={"eid": str})
    genotype = np.load(Path(config["prsformer_cache"]) / "genotypes.npy")
    rng = np.random.default_rng(422)
    frame["split"] = np.where(np.arange(len(frame)) < 180, "train", "test")
    frame["age"] = rng.normal(size=len(frame))
    frame["sex"] = rng.integers(0, 2, len(frame))
    frame["PC1"] = rng.normal(size=len(frame))
    frame["PC2"] = rng.normal(size=len(frame))
    frame["height"] = 100 + 4 * genotype[:, 0] - 2 * genotype[:, 1] + .2 * frame.age
    frame.to_csv(source, sep="\t", index=False)
    scores = pd.DataFrame(np.column_stack([genotype, rng.normal(size=(len(frame), 2))]),
                          columns=["csx.AFR", "csx.EAS", "csx.EUR", "csx.SAS"])
    scores.insert(0, "id", frame.eid)
    score_path = tmp_path / "independent_population_scores.tsv"
    scores.to_csv(score_path, sep="\t", index=False)
    output = tmp_path / "analysis"

    def no_raw_repreparation(*args, **kwargs):
        raise AssertionError("A valid configured cohort must not trigger default-cache preparation")

    import shared_cohort
    monkeypatch.setattr(shared_cohort, "ensure_inputs", no_raw_repreparation)
    common = ["--analysis-root", str(output), "--trait", "height", "--block-config", str(config_path),
              "--validation-n", "120"]
    target_main([*common, "--scores", str(score_path), "--models", "populations",
                 "--prs-standardization", "validation"])
    target_receipt = json.loads((output / "csx/height/targetfit/COMPLETE.json").read_text())
    assert Path(target_receipt["metadata"]["source"]) == source.resolve()
    assert target_receipt["files"][str(source.resolve())] == shared_cohort.sha256(source)
    assert target_receipt["roster"]["actual_validation_N"] == 120
    assert target_receipt["roster"]["test_N"] == 60

    block_main(["blockridge", *common, "--device", "cpu", "--threads", "1"])
    block_receipts = list((output / "csx_blocks").glob("*/height/blockridge/COMPLETE.json"))
    assert len(block_receipts) == 1
    receipt = json.loads(block_receipts[0].read_text())
    assert Path(receipt["cohort_metadata"]["source"]) == source.resolve()
    assert receipt["contract"]["cohort_source_sha256"] == shared_cohort.sha256(source)
    assert receipt["contract"]["roster"] == target_receipt["roster"]
    performance = pd.read_csv(block_receipts[0].parent / "performance.tsv", sep="\t").iloc[0]
    assert performance.N_total == 60 and performance.build_N + performance.tune_N == 120
    assert not (output / "prsformer/cache").exists()

    # Explicit sources remain authoritative even when a different cache config
    # is supplied. Prepared inputs continue through load_cohort's bound reader.
    args = SimpleNamespace(cohort="declared.tsv", prepared="declared.joblib", block_config=str(config_path),
                           analysis_root=str(output))
    assert cohort_source(args) == "declared.tsv"
    args.cohort = None
    assert cohort_source(args) is None
