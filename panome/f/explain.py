"""Individual reference evidence, blockwise mosaic matching and copy stability."""
import gzip
import json
import numpy as np
import pandas as pd
from reference import ReferencePanel, kernel
from survival import metrics


class MosaicReference:
    """Different real reference people may explain different molecular blocks.

    PC1 block coordinates are descriptive and learned without Y. Local outcome
    borrowing is supervised. There is no chromosome order or recombination HMM.
    """
    def fit(self, profiles, x, ids, groups, y, w, anchors, names, a):
        self.names, self.panels = names, []
        self.k = min(a.match_k[0], len(anchors)-1)
        for j in range(profiles.shape[1]):
            panel = ReferencePanel().fit(profiles[:, j:j+1], x, ids, groups, y, w,
                anchors, a.donor_neighbors, a.prior_strength)
            self.panels.append(panel)
        self.weights = np.ones(len(self.panels))/len(self.panels)
        return self

    def block_predictions(self, profiles, ids, groups):
        predictions, nearest_ids, nearest_distance = [], [], []
        for j, panel in enumerate(self.panels):
            from reference import nearest
            jj, dd = nearest(panel.z, profiles[:, j:j+1], self.k,
                             ids, panel.ids, groups, panel.groups)
            ww = kernel(dd)
            predictions.append(np.sum(ww*panel.local_risk[jj], axis=1))
            nearest_ids.append(panel.ids[jj[:, 0]])
            nearest_distance.append(dd[:, 0])
        return np.array(predictions).T, np.array(nearest_ids).T, np.array(nearest_distance).T

    def tune(self, values, y, w):
        null = np.average(y, weights=w)
        baseline = np.average((y-null)**2, weights=w)
        scores = np.array([np.average((y-values[:, j])**2, weights=w) for j in range(values.shape[1])])
        gain = np.maximum(baseline-scores, 0)
        # Equal weights remain a comparator; small shrinkage avoids all weight
        # collapsing onto one noisy block in a finite tuning sample.
        self.weights = (.1/len(gain)+.9*gain/gain.sum()) if gain.sum() > 0 else np.ones(len(gain))/len(gain)
        return pd.DataFrame(dict(module=self.names, tuning_Brier=scores,
                                 positive_Brier_gain=gain, copy_weight=self.weights))


def individual_evidence(out, id_col, ids, x, panel, matches, weights, features,
                        profiles, module_names, mosaic_predictions, mosaic_ids,
                        mosaic_distances, mosaic_weights, max_features=5):
    borrowed = np.zeros_like(x)
    for k in range(matches.shape[1]):
        borrowed += weights[:, k:k+1]*panel.x[matches[:, k]]
    with gzip.open(out/"test_feature_explanations.csv.gz", "wt", encoding="utf-8") as handle:
        import csv
        writer = csv.writer(handle)
        writer.writerow([id_col,"evidence_type","feature","target_z","reference_mixture_z","difference_z"])
        for i, eid in enumerate(ids):
            difference = x[i]-borrowed[i]
            agreement = np.minimum(np.abs(x[i]), np.abs(borrowed[i]))/(1+np.abs(difference))
            agreement[np.sign(x[i]) != np.sign(borrowed[i])] = -1
            for label, order in [("shared_extreme", np.argsort(-agreement)[:max_features]),
                                 ("mismatch", np.argsort(-np.abs(difference))[:max_features])]:
                for f in order:
                    writer.writerow([eid,label,features[f],float(x[i,f]),float(borrowed[i,f]),float(difference[f])])
    with gzip.open(out/"test_person_cards.jsonl.gz", "wt", encoding="utf-8") as handle:
        for i, eid in enumerate(ids):
            top = np.argsort(-np.abs(x[i]-borrowed[i]))[:max_features]
            card = dict(eid=str(eid), reference_people=[dict(eid=str(panel.ids[j]), weight=float(w),
                            local_net_risk=float(panel.local_risk[j])) for j,w in zip(matches[i],weights[i])],
                        molecular_mismatches=[dict(feature=features[j],target_z=float(x[i,j]),
                            reference_z=float(borrowed[i,j])) for j in top],
                        modules=[dict(module=name, coordinate=float(profiles[i,j]),
                            reference_eid=str(mosaic_ids[i,j]), distance=float(mosaic_distances[i,j]),
                            local_net_risk=float(mosaic_predictions[i,j]), mixture_weight=float(mosaic_weights[j]))
                            for j,name in enumerate(module_names)])
            handle.write(json.dumps(card, allow_nan=False)+"\n")
    nmodule = len(module_names)
    pd.DataFrame({id_col: np.repeat(ids, nmodule), "module": np.tile(module_names, len(ids)),
                  "target_coordinate": profiles.ravel(), "reference_eid": mosaic_ids.ravel(),
                  "module_distance": mosaic_distances.ravel(), "module_local_risk": mosaic_predictions.ravel(),
                  "mixture_weight": np.tile(mosaic_weights, len(ids))}).to_csv(
                      out/"test_mosaic_matches.csv.gz", index=False, compression="gzip")
    utilization = pd.Series(panel.ids[matches.ravel()]).value_counts()
    pd.DataFrame({"reference_eid":panel.ids,
                  "test_matches":utilization.reindex(panel.ids, fill_value=0).to_numpy(),
                  "local_net_risk":panel.local_risk}).to_csv(out/"reference_utilization.csv", index=False)


def masked_copy_stability(raw, p, prep, geometry, panel, spec, original_matches,
                          ids, groups, a, out):
    """Hold out query assays, rerun matching, assess the masked values only.

    These are empirical stress tests, not individual accuracy guarantees.
    Use deterministic sample/feature draws; do not consult test outcomes.
    """
    rng = np.random.default_rng(a.seed+400)
    take = np.sort(rng.choice(len(raw), min(a.explanation_samples, len(raw)), replace=False))
    if not len(take):
        return pd.DataFrame()
    original_x, original_observed = prep.transform(raw[take], p.iloc[take])
    rows = []
    count = max(1, int(len(prep.keep)*a.mask_fraction))
    def coordinates(values):
        if spec.get("space") == "full_proteome":
            return values
        zz = geometry.transform(values)
        return zz[:, :geometry.unsupervised.n_components_] if spec.get("space") == "unsupervised" else zz
    original_p, _, original_jj, _ = panel.predict(coordinates(original_x), original_x, original_observed,
        ids[take], groups[take] if groups is not None else None, spec["k"], spec["mode"])
    for repeat in range(a.mask_repeats):
        mask = np.zeros_like(original_observed)
        for j in range(len(take)):
            available = np.flatnonzero(original_observed[j])
            mask[j, rng.choice(available, min(count,len(available)), replace=False)] = True
        hidden = raw[take].copy()
        local_i, local_j = np.where(mask)
        hidden[local_i, prep.keep[local_j]] = np.nan
        hx, observed = prep.transform(hidden, p.iloc[take])
        hz = coordinates(hx)
        pred, _, jj, ww = panel.predict(hz, hx, observed, ids[take], groups[take] if groups is not None else None,
                                       spec["k"], spec["mode"])
        reconstructed = np.zeros_like(hx)
        for k in range(jj.shape[1]):
            reconstructed += ww[:,k:k+1]*panel.x[jj[:,k]]
        for j, idx in enumerate(take):
            mse = float(np.mean((original_x[j,mask[j]]-reconstructed[j,mask[j]])**2))
            null_mse = float(np.mean(original_x[j,mask[j]]**2))
            a_set, b_set = set(original_jj[j]), set(jj[j])
            rows.append(dict(eid=ids[idx], repeat=repeat+1, masked_assays=int(mask[j].sum()),
                             masked_RMSE=np.sqrt(mse), mean_reference_RMSE=np.sqrt(null_mse),
                             reconstruction_skill_vs_zero=1-mse/null_mse if null_mse>0 else np.nan,
                             match_Jaccard=len(a_set & b_set)/len(a_set | b_set),
                             absolute_raw_risk_change=float(abs(pred[j]-original_p[j]))))
    return pd.DataFrame(rows)
