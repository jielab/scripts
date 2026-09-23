"""Individual evidence and falsifiable explanation checks; no causal claims."""
from pathlib import Path
import gzip
import json
import numpy as np
import pandas as pd
from borrowing import risk_from_weights, stable_seed
from neural import encode, device_for
from common import log, dump


def write_match_table(path, idcol, ids, bank, detail, extra=None):
    jj = detail["jj"]
    data = {idcol:np.repeat(ids, jj.shape[1]), "rank":np.tile(np.arange(jj.shape[1])+1, len(ids)),
        "reference_eid":bank.ids[jj.ravel()], "reference_observed_Y":detail["labels"].ravel(),
        "normalized_copy_weight":detail["weights"].ravel(),
        "risk_coefficient":detail["coefficient"].ravel(),
        "raw_risk_contribution":(detail["coefficient"]*detail["labels"]).ravel(),
        "reference_IPCW":bank.w[jj.ravel()]}
    if "head_weights" in detail:
        for h in range(detail["head_gate"].shape[1]):
            head = detail["head_weights"][:,:,h]
            gate = np.repeat(detail["head_gate"][:,h,None],jj.shape[1],axis=1)
            coefficient = (1-detail["prior_weight"][:,None])*gate*head
            data[f"head_{h+1}_within_head_weight"] = head.ravel()
            data[f"head_{h+1}_gate"] = gate.ravel()
            data[f"head_{h+1}_raw_contribution"] = (coefficient*detail["labels"]).ravel()
    if extra:
        data.update({k:v.ravel() for k,v in extra.items()})
    table = pd.DataFrame(data)
    table.loc[table.normalized_copy_weight>0].to_csv(path, index=False, compression="infer")


def write_evidence(out, a, ids, x, observed, bank, d, calibrated, features, membership,
                   token_names, calibrator, strength):
    raw = np.sum(d["coefficient"]*d["labels"], axis=1)+d["prior_weight"]*bank.prior
    ww, labels = d["weights"], d["labels"]
    deletion = np.zeros_like(ww)
    for r in range(ww.shape[1]):
        reduced = ww.copy(); reduced[:, r] = 0
        sums = reduced.sum(1, keepdims=True)
        np.divide(reduced, sums, out=reduced, where=sums>0)
        risk = risk_from_weights(reduced, labels, bank.prior, strength)[0]
        risk = np.where(sums[:, 0]>0, risk, bank.prior)
        deletion[:, r] = calibrated-calibrator.predict(risk)
    flip_raw = d["coefficient"]*(1-2*labels)
    write_match_table(out/"test_reference_matches.csv.gz", a.id_col, ids, bank, d,
        {"calibrated_delta_delete_fixed_neighbor":deletion, "raw_delta_flip_donor_Y":flip_raw})
    pd.DataFrame({a.id_col:ids, "reference_contribution_sum":np.sum(d["coefficient"]*labels, axis=1),
        "prior_weight":d["prior_weight"], "build_prior":bank.prior,
        "prior_contribution":d["prior_weight"]*bank.prior, "reconstructed_raw_risk":raw,
        "reconstructed_calibrated_risk":calibrator.predict(raw),
        "calibrated_risk":calibrated, "max_absolute_reference_deletion_delta":np.abs(deletion).max(1)
        }).to_csv(out/"individual_risk_decomposition.csv", index=False)
    if not np.allclose(calibrator.predict(raw), calibrated, atol=1e-8):
        raise AssertionError("Reference contributions fail to reconstruct deployed probability")
    select = sorted(range(len(ids)), key=lambda i:stable_seed(ids[i], a.seed))[:a.explanation_samples]
    reconstruction, _ = bank.reconstruct({k:v[select] for k,v in d.items()}) if select else (None, None)
    with gzip.open(out/"individual_explanations.jsonl.gz", "wt", encoding="utf-8") as handle:
        for pos, i in enumerate(select):
            overlap = observed[i]
            mismatch = np.where(overlap, np.abs(x[i]-reconstruction[pos]), -np.inf)
            common = np.where(overlap & (x[i]*reconstruction[pos]>0),
                              np.minimum(np.abs(x[i]), np.abs(reconstruction[pos])), -np.inf)
            def describe(score):
                ix = np.argsort(-score, kind="stable")[:8]
                return [dict(feature=features[j], token=token_names[membership[j]],
                    query_standardized=float(x[i,j]), reference_reconstruction=float(reconstruction[pos,j]))
                    for j in ix if np.isfinite(score[j])]
            row = dict(eid=str(ids[i]), raw_borrowed_risk=float(raw[i]), calibrated_net_risk=float(calibrated[i]),
                prior_contribution=float(d["prior_weight"][i]*bank.prior), reference_ESS=float(d["ESS"][i]),
                shared_abnormal_assays=describe(common), largest_mismatches=describe(mismatch),
                interpretation="descriptive matching evidence; attention/contributions are not causal protein effects")
            handle.write(json.dumps(row, ensure_ascii=False)+"\n")
    count = pd.Series(bank.ids[d["jj"][d["weights"]>0]]).value_counts()
    utilization = pd.DataFrame({"reference_eid":bank.ids, "observed_Y":bank.y})
    utilization["test_matches"] = utilization.reference_eid.map(count).fillna(0).astype(int)
    utilization.to_csv(out/"reference_utilization.csv", index=False)


def main_prediction(bundle, x, observed, ids, groups, random=False):
    name = bundle["primary"]
    spec = bundle["specs"][name]
    model = bundle["encoders"][spec["encoder"]]
    cfg = bundle["config"]
    z, _, reconstructed = encode(model, x, observed, device_for(cfg["device"]), cfg["batch_size"], True)
    model.cpu()
    kw = {k:v for k,v in spec.items() if k not in ["encoder", "space"]}
    p, detail = bundle["banks"][name].match(z, ids, groups, **kw, random_candidates=random)
    return p, detail, reconstructed


def subset_indices(ids, count, seed):
    return np.array(sorted(range(len(ids)), key=lambda i:stable_seed(ids[i], seed))[:count], int)


def jaccard(a, b):
    return np.array([len(set(x)&set(y))/len(set(x)|set(y)) for x,y in zip(a,b)])


def masked_validation(bundle, x, observed, ids, groups, out, a):
    selected = subset_indices(ids, a.explanation_samples, a.seed)
    x, observed, ids = x[selected], observed[selected], ids[selected]
    groups = None if groups is None else groups[selected]
    original, orig_detail, _ = main_prediction(bundle, x, observed, ids, groups)
    bank = bundle["banks"][bundle["primary"]]
    rows, feature_stats = [], {}
    for repeat in range(a.mask_repeats):
        hidden = np.zeros_like(observed)
        for i, eid in enumerate(ids):
            valid = np.flatnonzero(observed[i])
            n = max(1, int(round(len(valid)*a.mask_fraction)))
            if len(valid):
                hidden[i, np.random.default_rng(stable_seed(eid, a.seed+repeat+981)).choice(valid, min(n,len(valid)), replace=False)] = True
        query = x.copy(); query[hidden] = 0
        mask = observed & ~hidden
        pred, detail, decoder = main_prediction(bundle, query, mask, ids, groups)
        _, random_detail, _ = main_prediction(bundle, query, mask, ids, groups, random=True)
        copied, cover = bank.reconstruct(detail)
        random_copy, _ = bank.reconstruct(random_detail)
        pca = bundle["geometry"].unsupervised
        pca_copy = pca.inverse_transform(pca.transform(query))
        methods = {"reference_attention":copied, "decoder":decoder, "random_references":random_copy,
                   "pca_reconstruction":pca_copy, "zero_build_mean":np.zeros_like(x)}
        nj = jaccard(orig_detail["jj"], detail["jj"])
        for name, estimate in methods.items():
            square = np.where(hidden, (estimate-x)**2, 0)
            n = hidden.sum(1)
            for i, eid in enumerate(ids):
                rows.append(dict(eid=eid, repeat=repeat+1, method=name, masked_count=int(n[i]),
                    masked_RMSE=float(np.sqrt(square[i].sum()/max(n[i],1))),
                    match_Jaccard=float(nj[i]) if name=="reference_attention" else np.nan,
                    absolute_raw_risk_change=float(abs(pred[i]-original[i])) if name=="reference_attention" else np.nan,
                    copied_assay_coverage=float(cover[i,hidden[i]].mean()) if name=="reference_attention" else np.nan))
            stats = feature_stats.setdefault(name, np.zeros((4, x.shape[1])))
            stats[0] += hidden.sum(0); stats[1] += np.where(hidden,x,0).sum(0)
            stats[2] += np.where(hidden,x*x,0).sum(0); stats[3] += square.sum(0)
        log("DONE", "masked_validation", f"repeat={repeat+1}/{a.mask_repeats}, n={len(ids)}")
    pd.DataFrame(rows).to_csv(out/"masked_reconstruction.csv", index=False)
    table = []
    for name, stats in feature_stats.items():
        n, total, total_sq, sse = stats
        ss = total_sq-total**2/np.maximum(n,1)
        r2 = np.full(len(n), np.nan)
        ok = (n>=3) & (ss>1e-8); r2[ok] = 1-sse[ok]/ss[ok]
        for j, feature in enumerate(bundle["retained_features"]):
            table.append(dict(method=name, feature=feature, masked_count=int(n[j]),
                masked_R2=r2[j], masked_RMSE=np.sqrt(sse[j]/n[j]) if n[j] else np.nan))
    pd.DataFrame(table).to_csv(out/"masked_feature_metrics.csv", index=False)


def module_perturbations(bundle, x, observed, ids, groups, out, a):
    ix = subset_indices(ids, min(a.explanation_samples, a.perturbation_samples), a.seed+3)
    if not len(ix):
        return
    x, observed, ids = x[ix], observed[ix], ids[ix]
    groups = None if groups is None else groups[ix]
    baseline, detail, _ = main_prediction(bundle, x, observed, ids, groups)
    cal = bundle["calibrators"][bundle["primary"]]
    original = cal.predict(baseline)
    rows = []
    for module, name in enumerate(bundle["token_names"]):
        query, mask = x.copy(), observed.copy()
        mask[:, bundle["membership"]==module] = False
        query[:, bundle["membership"]==module] = 0
        changed, other, _ = main_prediction(bundle, query, mask, ids, groups)
        calibrated_changed = cal.predict(changed)
        overlap = jaccard(detail["jj"], other["jj"])
        for i, eid in enumerate(ids):
            rows.append(dict(eid=eid, token=name, original_risk=original[i],
                risk_after_module_mask=calibrated_changed[i],
                delta_original_minus_masked=original[i]-calibrated_changed[i],
                match_Jaccard=overlap[i], interpretation="missing-module sensitivity; not a treatment effect"))
        log("DONE", "module_perturbation", name)
    pd.DataFrame(rows).to_csv(out/"module_perturbations.csv.gz", index=False, compression="gzip")


def same_risk_pairs(ids, risk, profile, names, caliper=.01, max_pairs=100):
    order = np.argsort(risk, kind="stable")
    used, rows = set(), []
    # Deterministic, no outcome inspection; cap candidates to avoid an N^2 matrix.
    for pos in np.linspace(0, max(0,len(order)-1), min(len(order), max_pairs*20), dtype=int):
        i = order[pos]
        if i in used: continue
        lo = np.searchsorted(risk[order], risk[i]-caliper)
        hi = np.searchsorted(risk[order], risk[i]+caliper, side="right")
        candidate = order[lo:hi]
        candidate = np.array([j for j in candidate if j!=i and j not in used], int)
        if not len(candidate): continue
        if len(candidate)>500: candidate = candidate[np.linspace(0,len(candidate)-1,500,dtype=int)]
        distance = np.linalg.norm(profile[candidate]-profile[i], axis=1)
        j = candidate[np.argmax(distance)]; used.update([i,j])
        contrast = np.argsort(-np.abs(profile[i]-profile[j]))[:3]
        rows.append(dict(eid_A=ids[i], eid_B=ids[j], predicted_risk_A=risk[i], predicted_risk_B=risk[j],
            profile_distance=float(np.max(distance)), largest_contrasts=";".join(names[k] for k in contrast)))
        if len(rows)>=max_pairs: break
    return pd.DataFrame(rows, columns=["eid_A","eid_B","predicted_risk_A","predicted_risk_B","profile_distance","largest_contrasts"])
