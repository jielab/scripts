"""Bounded CPU-only resampling; preserve the original bootstrap RNG stream."""
import numpy as np
from joblib import Parallel, delayed, parallel_config


def _metrics(true, predicted, indices, classes):
    true, predicted = true[indices], predicted[indices]
    support = np.bincount(true, minlength=classes)
    denominator = support + np.bincount(predicted, minlength=classes)
    correct = true == predicted
    hits = np.bincount(true[correct], minlength=classes)
    f1 = np.divide(2 * hits, denominator,
                   out=np.zeros(classes, dtype=float), where=denominator > 0)
    return (float(correct.mean()), float(f1[denominator > 0].mean()),
            float(np.dot(f1, support) / len(indices)))


def bootstrap_metrics(y_true, y_pred, n_boot, seed, cores):
    labels, encoded = np.unique(np.concatenate([y_true, y_pred]), return_inverse=True)
    n = len(y_true)
    if not n or n_boot < 1: raise ValueError('Bootstrap requires observations and positive replicates')
    rng = np.random.default_rng(int(seed))
    workers = max(1, min(int(cores), int(n_boot)))
    # Threads share label arrays; each replicate only allocates one index vector.
    with parallel_config(backend='threading'):
        rows = Parallel(n_jobs=workers, pre_dispatch=workers)(
            delayed(_metrics)(encoded[:n], encoded[n:], rng.integers(0, n, size=n), len(labels))
            for _ in range(int(n_boot)))
    return dict(zip(('accuracy', 'macro_f1', 'weighted_f1'), np.asarray(rows).T))
