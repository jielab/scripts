"""Auditable COPY and risk borrowing. Y enters values, never matching scores."""
import hashlib
import copy
import numpy as np
import pandas as pd
from scipy.special import softmax
from reference import diverse
from survival import loss
from attention import numpy_head_probabilities, numpy_mixture_scores


def stable_seed(eid, salt=0):
    return int.from_bytes(hashlib.blake2b((str(eid)+":"+str(salt)).encode(), digest_size=8).digest(), "little")


def select_references(z, quality, size, mode, seed, utility=None, balanced=False):
    """Actual participants; fit scores are never interpreted as data validity."""
    known = quality.known_label.to_numpy(bool)
    y = quality.horizon_label.to_numpy(int)
    pool = np.flatnonzero(known)
    if mode == "all":
        return pool
    if mode == "topfit":
        return pool[np.argsort(quality.OOF_logloss.to_numpy()[pool], kind="stable")[:size]]
    if mode in ["reliable", "utility"]:
        pool = np.flatnonzero(known & quality.reliable_candidate.to_numpy(bool))
    size = min(size, len(pool))
    if size < 2:
        return np.array([], int)
    fraction = .5 if balanced else float(y[known].mean())
    counts = [size-int(round(size*fraction)), int(round(size*fraction))]
    if len(np.unique(y[pool])) == 2:
        counts[1] = max(1, min(size-1, counts[1])); counts[0] = size-counts[1]
    counts = [min(counts[v], int(np.sum(y[pool] == v))) for v in [0, 1]]
    for v in [0, 1]:
        counts[v] += min(size-sum(counts), int(np.sum(y[pool] == v))-counts[v])
    gain = quality.fit_gain.fillna(-1e6).to_numpy()
    result = []
    for label, count in enumerate(counts):
        candidates = pool[y[pool] == label]
        if not count:
            continue
        if mode == "random":
            take = np.random.default_rng(seed+label).choice(candidates, count, replace=False)
        elif mode == "stratified_fit":
            take = candidates[np.argsort(quality.OOF_logloss.to_numpy()[candidates], kind="stable")[:count]]
        else:
            score = np.ones(len(y)) if mode == "diversity" else gain.copy()
            if mode == "utility" and utility is not None:
                # Rank-based blend keeps scale/large-exposure donors from dominating.
                score = (pd.Series(score).rank(pct=True).to_numpy() +
                         2*pd.Series(utility).rank(pct=True).to_numpy())
            take = diverse(z, candidates, count, score, seed+label)
        result.extend(take)
    return np.sort(result)


def risk_from_weights(weights, labels, prior, strength):
    sum_sq = np.sum(weights**2, axis=1)
    ess = np.divide(1., sum_sq, out=np.zeros_like(sum_sq), where=sum_sq>0)
    prior_weight = np.divide(strength, ess+strength, out=np.ones_like(ess), where=(ess+strength)>0)
    coefficients = weights*(1-prior_weight[:, None])
    p = np.sum(coefficients*labels, axis=1)+prior_weight*prior
    return p, coefficients, prior_weight, ess


class ReferenceBank:
    def __init__(self, z, x, observed, ids, groups, y, w, indices, encoder=None,
                 inclusion_correction=False, seed=2026):
        self.indices = np.asarray(indices, int)
        self.z, self.x, self.observed = z[indices], x[indices], observed[indices]
        self.ids, self.groups = ids[indices], None if groups is None else groups[indices]
        self.y, self.w = y[indices].astype(float), w[indices].astype(float)
        self.prior = float(np.average(y, weights=w))
        self.seed = seed
        self.gate_weight = self.gate_bias = self.temperature = None
        self.cross_mode = "metric"
        self.query_weight = self.key_weight = None
        if encoder is not None:
            self.cross_mode = encoder.cross_mode
            self.query_weight = encoder.query.weight.detach().cpu().numpy().copy()
            self.key_weight = encoder.key.weight.detach().cpu().numpy().copy()
            self.gate_weight = encoder.gate.weight.detach().cpu().numpy().copy()
            self.gate_bias = encoder.gate.bias.detach().cpu().numpy().copy()
            t = encoder.temperature.detach().cpu().numpy()
            self.temperature = .02+1.98/(1+np.exp(-t))
        self.correction = np.ones(len(indices))
        if inclusion_correction:
            # Only a class sampling correction; does not correct within-class selection.
            for label in [0, 1]:
                ix = self.y == label
                if ix.any():
                    self.correction[ix] = np.mean(y[w > 0] == label)/np.mean(ix)
        permutation = np.random.default_rng(seed+841).permutation(len(indices))
        self.permuted_y, self.permuted_w = self.y[permutation], self.w[permutation]

    def head_scores(self, q):
        h = len(self.temperature); d = q.shape[1]//h
        scale = q.shape[1]**.5
        qq = ((q*scale)@self.query_weight.T).reshape(-1,h,d)
        kk = ((self.z*scale)@self.key_weight.T).reshape(-1,h,d)
        return np.einsum("bhd,rhd->brh",qq,kk,optimize=True)/(d**.5*self.temperature)

    def head_gate(self, q):
        return softmax((q@self.gate_weight.T+self.gate_bias).astype("float64"),axis=1)

    def scores(self, q, kernel="attention"):
        if self.gate_weight is None or kernel == "euclidean":
            return -np.maximum(np.sum(q*q, 1)[:, None]+np.sum(self.z*self.z, 1)[None, :]-2*q@self.z.T, 0)
        if self.cross_mode == "qkv":
            return numpy_mixture_scores(self.head_scores(q),self.head_gate(q))
        h = len(self.temperature); d = q.shape[1]//h
        qq, zz = q.reshape(-1, h, d), self.z.reshape(-1, h, d)
        distances = np.maximum(np.sum(qq*qq, 2)[:, None, :]+np.sum(zz*zz, 2)[None, :, :]
                               -2*np.einsum("bhd,rhd->brh", qq, zz, optimize=True), 0)
        gate = softmax(q@self.gate_weight.T+self.gate_bias, axis=1)
        return -np.sum(distances*gate[:, None, :]/self.temperature, axis=2)

    def match(self, q, ids, groups, k=20, temperature=1., strength=2., kernel="attention",
              random_candidates=False, permuted=False, literal=False, batch_size=256):
        k = min(int(k), len(self.ids))
        if k < 1:
            raise ValueError("Empty reference panel")
        all_j, all_w, all_distance = [], [], []
        all_heads, all_gates = [], []
        multihead = self.cross_mode=="qkv" and kernel in ["attention","equal_heads"]
        for begin in range(0, len(q), batch_size):
            query = q[begin:begin+batch_size]
            forbidden = ids[begin:begin+batch_size, None] == self.ids[None, :]
            if groups is not None and self.groups is not None:
                forbidden |= groups[begin:begin+batch_size, None] == self.groups[None, :]
            if multihead:
                logits = self.head_scores(query)/temperature
                logits[forbidden] = -np.inf
                gate = self.head_gate(query)
                if kernel=="equal_heads": gate[:] = 1/gate.shape[1]
                score = numpy_mixture_scores(logits,gate)
            else:
                score = self.scores(query,kernel)/temperature
                score[forbidden] = -np.inf
            if random_candidates:
                jj = []
                for i, eid in enumerate(ids[begin:begin+batch_size]):
                    valid = np.flatnonzero(~forbidden[i])
                    selected = np.random.default_rng(stable_seed(eid, self.seed)).choice(valid, min(k,len(valid)), replace=False)
                    padding = np.flatnonzero(forbidden[i])[:k-len(selected)]
                    jj.append(np.r_[selected,padding])
                jj = np.asarray(jj)
            else:
                jj = np.argsort(-score, axis=1, kind="stable")[:, :k]
            values = np.take_along_axis(score,jj,axis=1)
            ww = self.permuted_w if permuted else self.w
            measure = ww[jj]*self.correction[jj]
            if multihead:
                chosen = np.take_along_axis(logits,jj[:,:,None],axis=1)
                with np.errstate(divide='ignore'):
                    chosen = chosen+np.log(measure)[:,:,None]
                head_weights = numpy_head_probabilities(chosen)
                weights = np.sum(head_weights*gate[:,None,:],axis=2)
                all_heads.append(head_weights.astype('float32')); all_gates.append(gate)
                # A geometric support diagnostic independent of logit offsets.
                distance = np.sum((query-self.z[jj[:,0]])**2,axis=1)
                distance[weights.sum(1)==0] = np.inf
            else:
                finite = np.isfinite(values)
                maximum = values.max(1, keepdims=True)
                maximum = np.where(np.isfinite(maximum),maximum,0)
                weights = np.where(finite,np.exp(values.astype("float64")-maximum),0)*measure
                total = weights.sum(1, keepdims=True)
                np.divide(weights,total,out=weights,where=total>0)
                distance = -values.max(1)*temperature
            all_j.append(jj); all_w.append(weights); all_distance.append(distance)
        jj, weights = np.concatenate(all_j), np.concatenate(all_w)
        labels = (self.permuted_y if permuted else self.y)[jj]
        if literal:
            if k != 1:
                raise ValueError("Literal COPY needs k=1")
            strength = 0
            weights = (weights > 0).astype(float)  # literal binary copying is exact
        p, coefficient, pw, ess = risk_from_weights(weights, labels, self.prior, strength)
        if literal:
            p[ess==0] = np.nan  # no donor's Y exists to COPY; do not invent one
        diagnostics = dict(jj=jj, weights=weights, coefficient=coefficient, prior_weight=pw,
            ESS=ess, nearest_distance=np.concatenate(all_distance), labels=labels,
            entropy=-np.sum(weights*np.log(np.maximum(weights, 1e-30)), axis=1),
            case_weight=np.sum(weights*labels, axis=1))
        if multihead:
            diagnostics["head_weights"] = np.concatenate(all_heads)
            diagnostics["head_gate"] = np.concatenate(all_gates)
            diagnostics["head_risk"] = np.sum(diagnostics["head_weights"]*labels[:,:,None],axis=1)
        return p, diagnostics

    def reconstruct(self, detail):
        jj, ww = detail["jj"], detail["weights"]
        reconstructed = np.zeros((len(jj), self.x.shape[1]), dtype="float32")
        denominator = np.zeros_like(reconstructed)
        # A missing donor assay is never silently treated as a measured value.
        for rank in range(jj.shape[1]):
            weight = ww[:, rank:rank+1]*self.observed[jj[:, rank]]
            reconstructed += weight*self.x[jj[:, rank]]; denominator += weight
        np.divide(reconstructed, denominator, out=reconstructed, where=denominator > 0)
        return reconstructed, denominator > 0

    def describe(self, x, observed, p, detail):
        reconstruction, covered = self.reconstruct(detail)
        valid = observed & covered
        rmse = np.sqrt(np.sum(np.where(valid, (x-reconstruction)**2, 0), axis=1)/np.maximum(valid.sum(1), 1))
        ww, labels, jj = detail["weights"], detail["labels"], detail["jj"]
        variance = np.sum(ww*(labels-detail["case_weight"][:, None])**2, axis=1)
        most_weighted = np.argmax(ww, axis=1)
        return pd.DataFrame(dict(nearest_reference=np.where(detail["ESS"]>0,self.ids[jj[:,0]],""),
            highest_weight_reference=np.where(detail["ESS"]>0,self.ids[jj[np.arange(len(jj)), most_weighted]],""),
            nearest_distance=detail["nearest_distance"], reconstruction_RMSE=rmse,
            reference_ESS=detail["ESS"], attention_entropy=detail["entropy"],
            independent_references=np.sum(ww>0,axis=1),
            normalized_entropy=detail["entropy"]/max(np.log(jj.shape[1]), 1),
            reference_disagreement=np.sqrt(variance), case_weight=detail["case_weight"],
            case_references=np.sum(labels*(ww>0),axis=1), prior_weight=detail["prior_weight"],
            reconstructed_assay_fraction=covered.mean(1)))


def tune_bank(bank, z, ids, groups, y, w, ks, temperatures, priors, name,
              kernel="attention", literal=False):
    rows, best = [], None
    for k in ks:
        if k > len(bank.ids):
            continue
        for temp in temperatures:
            _, cached = bank.match(z, ids, groups, k=k, temperature=temp, strength=0.,
                                   kernel=kernel, literal=literal)
            for strength in priors:
                spec = dict(k=k, temperature=temp, strength=strength, kernel=kernel, literal=literal)
                p = risk_from_weights(cached["weights"], cached["labels"], bank.prior, strength)[0]
                ll = float(np.average(loss(y, p), weights=w))
                rows.append(dict(model=name, **spec, tune_logloss=ll))
                if best is None or ll < best[0]:
                    best = ll, spec
    if best is None:
        raise ValueError(f"No usable neighbor setting for {name}")
    return best[1], rows


def reference_utility(bank, z, ids, groups, y, w, k=32):
    """Tuning-set deletion utility, NOT an OOF subject-fit score or causality."""
    p, d = bank.match(z, ids, groups, k=min(k, len(bank.ids)), strength=2)
    summed, exposure = np.zeros(len(bank.ids)), np.zeros(len(bank.ids))
    ww, jj, labels = d["weights"], d["jj"], d["labels"]
    for r in range(jj.shape[1]):
        reduced = ww.copy(); reduced[:, r] = 0
        reduced /= reduced.sum(1, keepdims=True)
        removed = risk_from_weights(reduced, labels, bank.prior, 2)[0]
        delta = w*(loss(y, removed)-loss(y, p))
        np.add.at(summed, jj[:, r], delta)
        np.add.at(exposure, jj[:, r], w)
    value = summed/(exposure+20)  # shrink rarely used donor estimates toward zero
    return value, pd.DataFrame(dict(reference_eid=bank.ids, tuning_exposure=exposure,
                                   shrunk_deletion_utility=value, total_loss_difference=summed))


class SupportRule:
    """X-only support plus a tuning-trained expected-error diagnostic.

    Expected error is evaluated on a separate audit set. It is not a posterior
    uncertainty interval, guarantee for a person, or evidence of cohort invalidity.
    """
    columns = ["nearest_distance", "reconstruction_RMSE", "reference_ESS",
               "normalized_entropy", "prior_weight", "reconstructed_assay_fraction"]

    def fit(self, info, risk, y, w, quantile=.95, min_ess=3):
        from sklearn.ensemble import HistGradientBoostingRegressor
        self.quantile, self.min_ess = quantile, min_ess
        self.thresholds = {k: float(info[k].quantile(quantile)) for k in
                           ["nearest_distance", "reconstruction_RMSE"]}
        known = w > 0
        self.error_model = HistGradientBoostingRegressor(max_iter=80, max_leaf_nodes=5,
            min_samples_leaf=max(10, min(100, int(known.sum()/15))), l2_regularization=10,
            early_stopping=False, random_state=583)
        self.error_model.fit(info[self.columns].to_numpy()[known], (y[known]-risk[known])**2,
                             sample_weight=w[known])
        self.error_threshold = float(np.quantile(self.error_score(info), quantile))
        return self

    def error_score(self, info):
        return np.maximum(self.error_model.predict(info[self.columns].to_numpy()), 0)

    def apply(self, info, missing, unknown, max_missing=.2, use_error=True):
        reasons = [[] for _ in range(len(info))]
        tests = {key: info[key].to_numpy() > threshold for key, threshold in self.thresholds.items()}
        tests.update(low_reference_ESS=info.reference_ESS.to_numpy() < self.min_ess,
                     no_independent_reference=info.reference_ESS.to_numpy() == 0,
                     excess_missingness=np.asarray(missing)>max_missing, unknown_category=unknown)
        if use_error:
            tests["high_expected_error"] = self.error_score(info)>self.error_threshold
        for key, bad in tests.items():
            for j in np.flatnonzero(bad):
                reasons[j].append(key)
        return np.array([not v for v in reasons]), np.array([";".join(v) for v in reasons])
