"""Offline test double only; never imported by the production pipeline."""

import hashlib
from types import SimpleNamespace

import anndata
import numpy as np
import pandas as pd

SCORERS = ["ATAC", "DNASE", "CHIP_TF", "CHIP_HISTONE", "CAGE", "PROCAP"]


def vector(variant, scorer):
    seed = int.from_bytes(
        hashlib.sha256(f"{variant}|{scorer}".encode()).digest()[:8], "little"
    )
    return np.random.default_rng(seed).normal(size=6).astype("f4")


class FakeAtlas:
    calls = 0
    fail_at = None

    def scorer_metadata(self):
        return {
            s: SimpleNamespace(
                is_signed=True, track_metadata=pd.DataFrame(index=range(6))
            )
            for s in SCORERS
        }

    def query_variants(self, variants, *, requested_scorers, **kwargs):
        self.calls += 1
        if self.fail_at == self.calls:
            raise ValueError("simulated interrupted query")
        # Different scorers return different orderings, as the real API may do.
        answer = {}
        for i, scorer in enumerate(requested_scorers):
            order = list(reversed(variants)) if i % 2 else variants
            answer[scorer] = anndata.AnnData(
                X=np.stack([vector(v, scorer) for v in order]),
                obs=pd.DataFrame(
                    {"variant": order}, index=[str(j) for j in range(len(order))]
                ),
                var=pd.DataFrame(index=[f"track{j}" for j in range(6)]),
            )
        return answer
