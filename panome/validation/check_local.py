"""Local migration regressions. Synthetic data only unless --real-input is set."""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import tempfile
import numpy as np
import pandas as pd
import pyreadr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "f"))
from panome import configure, parser
from io_data import read_table, prepare, generate_demo
from preprocess import outcomes, MolecularPreprocessor, MetadataDesign
from common import words
from survival import CensoringKM, metrics


def reader_check():
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "metadata.rds"
        subprocess.run(["Rscript", "-e", 'args <- commandArgs(TRUE); '
                        'saveRDS(data.frame(eid=c("001","002"), '
                        'plate=c(5e-320,1e-319)),args[1])', str(path)], check=True)
        expected = pyreadr.read_r(str(path))[None].reset_index(drop=True)
        pd.testing.assert_frame_equal(expected, read_table(path), check_exact=True)
    print("PASS: binary RDS preserves exact subnormal values and leading-zero IDs")
    people = pd.DataFrame({"time": [1.] * 10 + [12.] * 20, "event": [1] * 10 + [0] * 20})
    y, w = CensoringKM().fit(people, 10).labels_weights(people)
    np.testing.assert_array_equal(w, np.ones(30))
    assert np.isclose(metrics(y, w, np.full(30, .2))["Brier_IPCW"], np.mean((y-.2)**2))
    print("PASS: IPCW matches ordinary Brier without early censoring")


def real_check():
    a = configure(parser().parse_args(["--max-samples", "1600"]))
    with tempfile.TemporaryDirectory(prefix="panome_real_input_") as tmp:
        out = Path(tmp)
        p = prepare(a, out)
        raw = np.load(out / "raw.npy")
        p, audit = outcomes(p, a)
        eligible = p.eligible.to_numpy()
        raw, p = raw[eligible], p.loc[eligible].reset_index(drop=True)
        from pipeline import split_people
        split = split_people(p, a)
        build = split == "build"
        keep = (np.mean(~np.isfinite(raw[build]), axis=0) < a.feature_missing) & (np.nanstd(raw[build], axis=0) > 1e-8)
        qc = np.mean(~np.isfinite(raw[:, keep]), axis=1) <= a.sample_missing
        raw, p, build = raw[qc], p.loc[qc].reset_index(drop=True), build[qc]
        prep = MolecularPreprocessor(a.feature_missing, words(a.residualize), words(a.categorical), a.transform)
        prep.fit(raw[build], p.loc[build], allowed=keep)
        x, _ = prep.transform(raw, p)
        clinical = MetadataDesign(words(a.covariates), words(a.categorical)).fit(p.loc[build])
        assert np.isfinite(x).all() and np.isfinite(clinical.transform(p)).all()
        # The public audit CSV must never be reloaded as training metadata.
        assert p["prot.plate"].dtype.kind in "fi"
        print(json.dumps(dict(input_audit=json.loads((out / "input_audit.json").read_text()),
                              outcomes=audit, after_qc=len(p), features=x.shape[1],
                              plate_levels=p["prot.plate"].nunique()), indent=2))
    print("PASS: real RDS + raw PPP proteins + sample QC + residualization; no full training")


def projection_check(run):
    import joblib
    from types import SimpleNamespace
    from pipeline import project
    bundle = joblib.load(run / "model_bundle.joblib")
    a = SimpleNamespace(**bundle["config"])
    p, omics = generate_demo(a)
    expected = pd.read_csv(run / "test_individuals.csv", dtype={a.id_col: str}).head(12)
    ids = expected[a.id_col].tolist()
    cols = list(dict.fromkeys([a.id_col] + words(a.covariates) + words(a.residualize) +
                             ([a.group_col] if a.group_col else [])))
    p = p.set_index(a.id_col).loc[ids].reset_index()[cols]
    omics = omics.set_index(a.id_col).loc[ids].reset_index()
    names = json.loads((run / "prediction_columns.json").read_text())
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        phe, assays, output = tmp / "phe.csv", tmp / "omics.csv", tmp / "pred.csv"
        p.to_csv(phe, index=False)
        omics.to_csv(assays, index=False)
        project(run, phe, assays, output)
        pred = pd.read_csv(output, dtype={a.id_col: str})
        np.testing.assert_allclose(pred[names], expected[names], atol=1e-6, rtol=1e-6)
        assert pred.nearest_reference.tolist() == expected.nearest_reference.tolist()
        p.iloc[::-1].to_csv(phe, index=False)
        omics.iloc[::-1].to_csv(assays, index=False)
        project(run, phe, assays, output)
        reversed_pred = pd.read_csv(output).iloc[::-1].reset_index(drop=True)
        np.testing.assert_allclose(reversed_pred[names], pred[names], atol=1e-6, rtol=1e-6)
        p.iloc[:1].to_csv(phe, index=False)
        omics.iloc[:1].to_csv(assays, index=False)
        project(run, phe, assays, output)
        np.testing.assert_allclose(pd.read_csv(output)[names], pred[names].iloc[:1], atol=1e-6, rtol=1e-6)
    matches = pd.read_csv(run / "test_reference_matches.csv")
    np.testing.assert_allclose(matches.groupby(a.id_col).copy_weight.sum(), 1)
    split = pd.read_csv(run / "split.csv").merge(generate_demo(a)[0][[a.id_col, a.group_col]], on=a.id_col)
    assert split.groupby(a.group_col).split.nunique().max() == 1
    print(f"PASS: {len(names)} models reproduce frozen test predictions without outcomes; batch/order invariance, normalized matches, and family separation")


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--real-input", action="store_true")
    cli.add_argument("--projection-run", type=Path)
    args = cli.parse_args()
    reader_check()
    if args.real_input:
        real_check()
    if args.projection_run:
        projection_check(args.projection_run)
