"""Reference selection, local outcome borrowing and query-only diagnostics."""
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import MiniBatchKMeans


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


def nearest(bank, query, k, ids=None, bank_ids=None, groups=None, bank_groups=None):
    """Exact deterministic neighbors, excluding self/relatives; batch invariant."""
    extra = 1
    if bank_groups is not None:
        extra += int(pd.Series(bank_groups).value_counts().max())
    request = min(len(bank), k+extra)
    model = NearestNeighbors(n_neighbors=request, algorithm="brute", metric="euclidean", n_jobs=1).fit(bank)
    all_j, all_d = [], []
    for start in range(0, len(query), 512):
        distances, indices = model.kneighbors(query[start:start+512])
        for pos, (dist, idx) in enumerate(zip(distances, indices), start):
            keep = np.ones(len(idx), bool)
            if ids is not None:
                keep &= bank_ids[idx] != ids[pos]
            if groups is not None and bank_groups is not None:
                keep &= bank_groups[idx] != groups[pos]
            ix = np.flatnonzero(keep)
            ix = ix[np.lexsort((idx[ix], dist[ix]))][:k]
            if len(ix) < k:
                raise ValueError("Insufficient independent references after self/family exclusion")
            all_j.append(idx[ix]); all_d.append(dist[ix])
    return np.asarray(all_j), np.asarray(all_d)


def kernel(dist):
    bandwidth = np.maximum(dist[:, -1:], 1e-6)
    w = np.exp(-.5*(dist/bandwidth)**2)
    return w/w.sum(1, keepdims=True)


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


def choose_panel(z, quality, size, mode, seed):
    y = quality.horizon_label.to_numpy()
    known = quality.known_label.to_numpy(bool)
    gain = quality.fit_gain.fillna(-1e6).to_numpy()
    pool = np.flatnonzero(known)
    if mode == "topfit":
        return pool[np.argsort(quality.OOF_logloss.to_numpy()[pool], kind="stable")[:size]]
    if mode == "all":
        return pool
    if mode == "reliable":
        pool = np.flatnonzero(quality.reliable_candidate.to_numpy(bool))
    size = min(size, len(pool))
    # Match observed class prevalence; do not construct a 50/50 case-control panel.
    event_fraction = y[known].mean()
    ncase = min(int(round(size*event_fraction)), int(np.sum(y[pool] == 1)))
    if size >= 2 and np.sum(y[pool] == 1) and ncase == 0:
        ncase = 1
    ncontrol = min(size-ncase, int(np.sum(y[pool] == 0)))
    ncase = min(size-ncontrol, int(np.sum(y[pool] == 1)))
    picked = []
    for label, count in [(0, ncontrol), (1, ncase)]:
        candidates = pool[y[pool] == label]
        if mode == "random":
            selected = np.random.default_rng(seed+label).choice(candidates, count, replace=False)
        else:
            score = gain if mode == "reliable" else np.ones(len(y))
            selected = diverse(z, candidates, count, score, seed+label)
        picked.extend(selected)
    return np.sort(picked)


class ReferencePanel:
    def fit(self, z, x, ids, groups, y, w, anchors, donor_k, prior_strength):
        if len(anchors) < 2:
            raise ValueError("Fewer than two usable prototypes")
        self.z, self.x = z[anchors].copy(), x[anchors].copy()
        self.ids, self.anchor_indices = ids[anchors], np.asarray(anchors)
        self.groups = groups[anchors] if groups is not None else None
        self.labels = y[anchors]
        self.outcome_weight = w[anchors]
        self.prior = float(np.average(y, weights=w))
        donor = np.flatnonzero(w > 0)
        bg = groups[donor] if groups is not None else None
        largest_group = int(pd.Series(bg).value_counts().max()) if bg is not None else 1
        self.donor_k = min(donor_k, len(donor)-largest_group-1)
        if self.donor_k < 10:
            raise ValueError("Insufficient independent donors for local risk")
        jj, dd = nearest(z[donor], self.z, self.donor_k, self.ids, ids[donor], self.groups, bg)
        donors = donor[jj]
        base_weight = kernel(dd)
        weights = base_weight*w[donors]
        weights /= weights.sum(1, keepdims=True)
        ess = 1/np.sum(weights**2, axis=1)
        local = np.sum(weights*y[donors], axis=1)
        self.local_risk = (local*ess+prior_strength*self.prior)/(ess+prior_strength)
        self.donor_ess = ess
        self.donor_events = np.sum(y[donors], axis=1)
        self.donor_indices, self.donor_weights = donors, weights
        return self

    def subset(self, indices):
        import copy
        result = copy.copy(self)
        lookup = {int(v): j for j, v in enumerate(self.anchor_indices)}
        take = np.asarray([lookup[int(v)] for v in indices])
        for name in ["z", "x", "ids", "anchor_indices", "labels", "outcome_weight",
                     "local_risk", "donor_ess", "donor_events", "donor_indices", "donor_weights"]:
            setattr(result, name, getattr(self, name)[take].copy())
        if self.groups is not None:
            result.groups = self.groups[take].copy()
        return result

    def predict(self, z, x, observed, ids, groups, k, mode="local", diagnostics=True):
        max_group = int(pd.Series(self.groups).value_counts().max()) if self.groups is not None else 1
        k = min(k, len(self.z)-max_group)
        if k < 1:
            raise ValueError("Panel cannot support independent matching")
        jj, dd = nearest(self.z, z, k, ids, self.ids, groups, self.groups)
        weights = kernel(dd)
        if mode == "label":
            weights *= self.outcome_weight[jj]
            weights /= weights.sum(1, keepdims=True)
        values = self.labels if mode == "label" else self.local_risk
        p = np.sum(weights*values[jj], axis=1)
        if not diagnostics:
            return p
        # Reconstruct the measured molecular profile, never the missing values
        # when calculating this diagnostic. It is not genotype imputation r2.
        reconstruction = np.zeros_like(x)
        for i in range(k):
            reconstruction += weights[:, i:i+1]*self.x[jj[:, i]]
        mse = np.sum(np.where(observed, (x-reconstruction)**2, 0), axis=1)/np.maximum(observed.sum(1), 1)
        dispersion = np.sqrt(np.sum(weights*(values[jj]-p[:, None])**2, axis=1))
        info = pd.DataFrame(dict(nearest_reference=self.ids[jj[:, 0]],
                                 nearest_distance=dd[:, 0], reconstruction_RMSE=np.sqrt(mse),
                                 reference_ESS=1/np.sum(weights**2, axis=1),
                                 donor_ESS=np.sum(weights*self.donor_ess[jj], axis=1),
                                 donor_event_support=np.sum(weights*self.donor_events[jj], axis=1),
                                 reference_disagreement=dispersion))
        return p, info, jj, weights


class AcceptanceRule:
    def fit(self, info, quantile, min_ess):
        self.quantile, self.min_ess = quantile, min_ess
        self.thresholds = {c: float(np.quantile(info[c], quantile)) for c in
                           ["nearest_distance", "reconstruction_RMSE", "reference_disagreement"]}
        return self

    def apply(self, info, missing, unknown):
        reasons = [[] for _ in range(len(info))]
        for col, threshold in self.thresholds.items():
            for j in np.flatnonzero(info[col].to_numpy() > threshold):
                reasons[j].append(col)
        for key, bad in [("low_reference_ESS", info.reference_ESS.to_numpy() < self.min_ess),
                         ("excess_missingness", missing), ("unknown_category", unknown)]:
            for j in np.flatnonzero(bad):
                reasons[j].append(key)
        return np.asarray([not r for r in reasons]), np.asarray([";".join(r) for r in reasons])


def make_profiles(x_build, geometry, features, seed, module_file=None):
    """Frozen descriptive axes; module signs do not assert causality/pathway activity."""
    if module_file:
        table = pd.read_csv(module_file, sep="\t", dtype=str)
        if not {"feature", "module"} <= set(table):
            raise ValueError("--module-file needs feature and module TSV columns")
        table = table[table.feature.isin(features)].copy()
        if table.empty:
            raise ValueError("No module-file features match measured feature names exactly")
        membership = [(str(m), [features.index(f) for f in t.feature.unique()]) for m, t in table.groupby("module")]
    else:
        loading = geometry.unsupervised.components_.T.copy()
        loading /= np.maximum(np.linalg.norm(loading, axis=1, keepdims=True), 1e-8)
        count = min(12, max(2, len(features)//15))
        labels = MiniBatchKMeans(count, random_state=seed, n_init=3, batch_size=512).fit_predict(loading)
        membership = [(f"data_block_{i+1:02d}", np.flatnonzero(labels == i).tolist()) for i in range(count)]
    weights = np.zeros((len(features), len(membership)), dtype="float32")
    rows = []
    for j, (name, idx) in enumerate(membership):
        if len(idx) == 1:
            v = np.ones(1)
        else:
            v = PCA(1, svd_solver="randomized", random_state=seed).fit(x_build[:, idx]).components_[0]
            if v[np.argmax(np.abs(v))] < 0:
                v = -v
        weights[idx, j] = v
        for f, b in zip(idx, v):
            rows.append(dict(module=name, feature=features[f], loading=float(b)))
    score = x_build @ weights
    scale = np.maximum(score.std(0), 1e-6)
    weights /= scale
    return weights, [m[0] for m in membership], pd.DataFrame(rows)
