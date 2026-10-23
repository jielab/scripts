"""Validate complete, frozen-test PRSformer outputs before aggregation."""
import numpy as np
import pandas as pd

from .evaluation import prediction_metrics


def performance(run, cache, traits, method):
    from prsformer import validate_evaluation
    validate_evaluation(run, cache)
    result = pd.read_csv(run / "test_predictions.tsv.gz", sep="\t", dtype={"eid": str})
    cohort = pd.read_csv(cache / "data.tsv.gz", sep="\t", dtype={"eid": str})
    cohort = cohort.loc[cohort.split.eq("test")].set_index("eid")
    if cohort.empty or not cohort.index.is_unique:
        raise ValueError("PRSformer prepared test roster must be nonempty and unique")
    if not result.split.eq("test").all() or result.duplicated(["eid", "trait"]).any():
        raise ValueError("PRSformer predictions must be unique held-out test rows")
    if not set(traits) <= set(result.trait):
        raise ValueError("PRSformer evaluation is missing requested traits")
    rows = []
    for trait in traits:
        group = result.loc[result.trait.eq(trait)].set_index("eid")
        if len(group) != len(cohort) or set(group.index) != set(cohort.index):
            raise ValueError("PRSformer predictions have an incomplete test roster: " + trait)
        expected = pd.to_numeric(cohort.loc[group.index, trait]).to_numpy(float)
        if not np.allclose(group.y, expected, rtol=0, atol=1e-10, equal_nan=True):
            raise ValueError("PRSformer test outcomes differ from prepared inputs: " + trait)
        rows.append(dict(method=method, trait=trait, stage="test",
                         **prediction_metrics(group.y.to_numpy(), group.prediction.to_numpy(),
                                              "binary" if trait == "t2dm" else "continuous")))
    return pd.DataFrame(rows)
