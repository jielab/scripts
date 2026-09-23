"""Classical PCA geometry and quality-weighted real-person coverage selection."""
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

class MolecularGeometry:
    def fit(self, x, importance, dimensions, seed):
        d = min(dimensions, len(x)-2, x.shape[1]-1)
        self.unsupervised = PCA(d, svd_solver="randomized", random_state=seed).fit(x)
        importance = np.maximum(importance, 0)
        self.has_supervised_signal = bool(importance.max() > 0)
        self.weight = (np.sqrt(.02+.98*importance/max(importance.max(), 1e-12))
                       if self.has_supervised_signal else np.ones_like(importance)).astype("float32")
        self.supervised = PCA(d, svd_solver="randomized", random_state=seed+1).fit(x*self.weight)
        self.scale_u = np.sqrt(np.maximum(self.unsupervised.explained_variance_, 1e-6)).astype("float32")
        self.scale_s = np.sqrt(np.maximum(self.supervised.explained_variance_, 1e-6)).astype("float32")
        self.normalizer = np.sqrt(2*d)
        return self

    def transform(self, x):
        return np.c_[self.unsupervised.transform(x)/self.scale_u,
                     self.supervised.transform(x*self.weight)/self.scale_s].astype("float32")/self.normalizer

def diverse(z, pool, size, quality, seed):
    """Quality-weighted farthest-first real exemplars, not synthetic centroids."""
    pool = np.asarray(pool)
    size = min(int(size), len(pool))
    if size <= 0:
        return np.array([], int)
    rng = np.random.default_rng(seed)
    q = np.asarray(quality)[pool]
    q = .25+.75*pd.Series(q).rank(pct=True).to_numpy()
    # Start near a representative centre; cap leverage of very distant points.
    cloud = z[pool]
    central = np.sum((cloud-np.median(cloud, axis=0))**2, axis=1)
    density_cap = np.quantile(central, .95)+1e-8
    q *= np.minimum(1, density_cap/np.maximum(central, 1e-8))
    first = int(np.argmax(q/(1+central)))
    chosen, distance = [], np.full(len(pool), np.inf)
    available = np.ones(len(pool), bool)
    current = first
    for _ in range(size):
        chosen.append(pool[current]); available[current] = False
        delta = np.sum((cloud-cloud[current])**2, axis=1)
        distance = np.minimum(distance, delta)
        score = distance*q
        score[~available] = -1
        current = int(np.argmax(score+rng.uniform(0, 1e-12, len(score))))
    return np.asarray(chosen, int)

