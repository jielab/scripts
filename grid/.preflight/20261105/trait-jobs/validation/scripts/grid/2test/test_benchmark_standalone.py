"""Scientific integration tests: genuine fits, fixed families, and sealed test Y."""
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
from shared_cohort import small_validation_roster, frozen_roster, load_cohort
from target_benchmark import evaluate, fit_score_model


def fixture(folder, binary=False, n=600):
    rng = np.random.default_rng(81)
    score = rng.normal(size=(n, 4))
    side = np.column_stack([rng.normal(55, 4, n), rng.integers(0, 2, n), rng.normal(size=(n,2))])
    linear = score @ np.array([2.5, -1.3, .6, .25]) + side @ np.array([.2, .3, .1, -.2])
    y = rng.binomial(1, 1/(1+np.exp(-linear+11))) if binary else linear + 7
    d = pd.DataFrame(dict(id=[str(i) for i in range(n)], group_id=[f"g{i//2}" for i in range(n)],
        split=np.where(np.arange(n)<400, "train", "test"), y=y,
        age=side[:,0], sex=side[:,1], PC1=side[:,2], PC2=side[:,3], ancestry="EUR"))
    s = pd.DataFrame(score, columns=["csx.AFR", "csx.EAS", "csx.EUR", "csx.SAS"])
    s.insert(0, "id", d.id)
    s["disco"] = score @ np.array([.1,.2,.4,.3])
    cohort, scores = folder/"cohort.tsv", folder/"scores.tsv"
    d.to_csv(cohort, sep="\t", index=False)
    s.to_csv(scores, sep="\t", index=False)
    args = SimpleNamespace(analysis_root=str(folder/"analysis"), scores=str(scores), score_dir=str(folder),
        cohort=str(cohort), prepared=None, group_col=None, validation_n=200, seed=20261102,
        covariates="age,sex,PC1,PC2", models="populations", prs_standardization="validation")
    return d, s, args


def test_family_budget_assignment_is_label_and_row_order_blind(tmp_path):
    frame, _, _ = fixture(tmp_path)
    a = small_validation_roster(frame, 121)
    frame["y"] = frame.y[::-1].to_numpy()
    b = small_validation_roster(frame.sample(frac=1, random_state=2), 121)
    pd.testing.assert_frame_equal(a.sort_values("id").reset_index(drop=True), b.sort_values("id").reset_index(drop=True))
    assert a.role.eq("target_validation").sum() == 120
    assert a.groupby("group_id").role.nunique().max() == 1
    assert a.groupby("group_id").block_role.nunique().max() == 1
    assert set(a.loc[a.role.eq("test"), "id"]) == set(frame.loc[frame.split.eq("test"), "id"])
    frame.loc[0, "group_id"] = frame.loc[500, "group_id"]
    with pytest.raises(ValueError, match="crosses"):
        small_validation_roster(frame, 121)


def test_actual_four_beta_fit_then_test_label_permutation_cannot_change_model(tmp_path):
    frame, score, args = fixture(tmp_path)
    perf, receipt = evaluate(args, "height")
    assert perf.iloc[0].prediction_r2 == pytest.approx(1., abs=1e-12)
    out = Path(args.analysis_root)/"csx/height/targetfit"
    fitted = json.loads((out/"_LOCAL_ONLY/frozen_models.json").read_text())
    sd = np.array(fitted["populations"]["score_scaling"]["statistics"]["target_validation"]["sd"])
    beta = np.array(fitted["populations"]["model"]["beta"][:4]) / sd
    np.testing.assert_allclose(beta, [2.5,-1.3,.6,.25], atol=1e-10)
    assert len(pd.read_csv(out/"population_betas.tsv",sep="\t")) == 4
    predictions = pd.read_csv(receipt["predictions_file"],sep="\t")
    original_model = (out/"_LOCAL_ONLY/frozen_models.json").read_bytes()
    test = frame.split.eq("test")
    frame.loc[test,"y"] = frame.loc[test,"y"].to_numpy()[::-1] + 50
    frame.to_csv(args.cohort,sep="\t",index=False)
    changed, new = evaluate(args,"height")
    assert (out/"_LOCAL_ONLY/frozen_models.json").read_bytes() == original_model
    pd.testing.assert_frame_equal(predictions,pd.read_csv(new["predictions_file"],sep="\t"))
    assert changed.iloc[0].prediction_r2 < 0
    assert len(pd.read_csv(new["test_truth_file"],sep="\t")) == 200
    assert (out/"COMPLETE.json").is_file()


def test_binary_and_disco_own_standalone_metrics_and_bound_receipt(tmp_path):
    _, _, args = fixture(tmp_path, binary=True)
    perf, receipt = evaluate(args,"t2dm")
    assert perf.iloc[0].N == 200 and 0 <= perf.iloc[0].auc <= 1
    assert 0 <= perf.iloc[0].brier <= 1
    disco, receipt = evaluate(args,"t2dm",method="disco")
    from disco_publication import read_performance
    saved, status = read_performance(args.analysis_root,"t2dm")
    assert status == "completed" and len(saved) == 1
    assert disco.iloc[0].validation_N == 200
    with open(receipt["predictions_file"],"ab") as stream:
        stream.write(b"tampered")
    assert read_performance(args.analysis_root,"t2dm")[0] is None


def test_public_score_benchmark_entry_runs_without_grid_model_cache(tmp_path):
    _, _, args = fixture(tmp_path)
    for script, module in (("2.csx.sh","targetfit"),("3.disco.sh","evaluate")):
        result = subprocess.run(["bash",str(ROOT/script),module,"--analysis-root",args.analysis_root,
            "--traits","height","--cohort",args.cohort,"--scores",args.scores,
            "--validation-n","200","--covariates",args.covariates,
            "--prs-standardization","validation","--models","populations"],
            check=True,capture_output=True,text=True,env=dict(os.environ,WORKFLOW_WORKER="1",submit="local"))
        assert "untouched test N=200" in result.stdout
    assert not (Path(args.analysis_root)/"prsformer").exists()
    assert not (Path(args.analysis_root)/"cache/abm/prepared.joblib").exists()


def test_whole_group_roster_cannot_be_reselected_by_later_method(tmp_path):
    frame, _, _ = fixture(tmp_path)
    frozen_roster(frame,tmp_path/"a",200,31)
    with pytest.raises(ValueError,match="roster/settings changed"):
        frozen_roster(frame,tmp_path/"a",220,31)


@pytest.mark.parametrize("method,filename", [("csx", "1csx.scores.rds"), ("disco", "2disco.scores.rds")])
@pytest.mark.parametrize("explicit", [False, True])
def test_subset_reads_matching_score_and_preserves_genomewide_performance(tmp_path, monkeypatch, method, filename, explicit):
    _, scores, args = fixture(tmp_path)
    args.scores = None
    args.score_subset = "chr22"
    expected = tmp_path / "height" / "chr22" / filename
    expected.parent.mkdir(parents=True)
    scores.to_csv(expected, sep="\t", index=False)
    if explicit:
        args.scores = str(expected)
        args.score_subset = ""
    import target_benchmark
    def read_score(path):
        assert Path(path) == expected
        return pd.read_csv(path, sep="\t", dtype={"id": str})
    monkeypatch.setattr(target_benchmark, "table", read_score)
    base = (Path(args.analysis_root) / "csx/height/targetfit" if method == "csx" else
            Path(args.analysis_root) / "disco/height")
    base.mkdir(parents=True)
    previous = base / "performance.tsv"
    previous.write_text("previous genome-wide result\n")
    performance, receipt = evaluate(args, "height", method=method)
    assert previous.read_text() == "previous genome-wide result\n"
    assert (base / "chr22/performance.tsv").is_file()
    assert performance.iloc[0].chromosome_subset == "chr22"
    assert receipt["scores"] == str(expected.resolve())
