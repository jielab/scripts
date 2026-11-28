#!/usr/bin/env python3
"""Train-frozen block neighbourhoods and conservative OOF residual borrowing.

The full-space consensus geometry uses every variable LD-block score. Exact
tiled CPU search counts all eligible radius members; k only limits prediction.
Explicit PCA/global-score sensitivity routes retain their declared geometry.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import joblib
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.ensemble import HistGradientBoostingRegressor

PRIMARY = "csx.population"
RAW_COLUMNS = [f"csx.{p}.raw_score" for p in ("EUR", "AFR", "EAS", "SAS")]
COVERAGES = (.2, .4, .5, .6, .8, 1.)


def table(path):
    return pd.read_csv(path, sep="\t", dtype={"id": str, "group_id": str, "ancestry": str})


def align(frame, ids, name):
    if "id" not in frame or frame.id.isna().any() or frame.id.duplicated().any():
        raise ValueError(f"{name}: unique nonmissing id required")
    if not set(ids) <= set(frame.id):
        raise ValueError(f"{name}: does not cover requested roster")
    return frame.set_index("id").loc[list(ids)].reset_index()


def regional_feature_schema(manifest,width):
    """Ordered scientific identities, excluding paths and mutable QC counts."""
    frame=table(manifest) if isinstance(manifest,(str,Path)) else pd.DataFrame(manifest).copy()
    if len(frame)!=width or width<1 or not frame.columns.is_unique:
        raise ValueError("Regional feature manifest does not match the input width/schema")
    columns=[c for c in ("block_id","name","map","population","chrom38","start38","end38","chrom37","start37","end37") if c in frame]
    if not columns:
        raise ValueError("Regional feature manifest needs block identities or coordinates")
    schema=frame[columns].fillna("").astype(str).to_dict("records")
    identities=[json.dumps(row,sort_keys=True) for row in schema]
    if len(set(identities))!=width:
        raise ValueError("Duplicate regional feature identities")
    return schema


def development_outcomes(outcomes, cohort):
    """Filter identities/roles BEFORE reading any outcome values into fitting."""
    if outcomes.id.isna().any() or outcomes.id.duplicated().any():
        raise ValueError("Outcomes require unique nonmissing IDs")
    development_ids = cohort.loc[cohort.role.isin(["train", "tune"]), "id"]
    selected = align(outcomes.loc[outcomes.id.isin(development_ids)], development_ids, "development outcomes")
    values = pd.to_numeric(selected.y, errors="raise")
    if np.isinf(values).any():
        raise ValueError("Infinite development outcomes")
    return pd.Series(values.to_numpy(float), index=selected.id).reindex(cohort.id).to_numpy(float)


def check_roster(cohort):
    required = {"id", "group_id", "ancestry", "role"}
    if not required <= set(cohort) or cohort[list(required)].isna().any().any() or cohort.id.duplicated().any():
        raise ValueError("Roster needs unique id and complete group_id/ancestry/role")
    if not set(cohort.role) <= {"train", "tune", "calibration", "test"}:
        raise ValueError("Unknown roster role")
    if not {"train", "tune", "calibration", "test"} <= set(cohort.role):
        raise ValueError("All four disjoint roles are required")
    if cohort.groupby("group_id").role.nunique().max() != 1:
        raise ValueError("A family/group crosses roles")


def check_oof(cohort, predictions, observed_train):
    column = PRIMARY + ".oof_prediction"
    if column not in predictions or "oof_fold" not in predictions:
        raise ValueError("Primary group-held-out OOF predictions and oof_fold are required")
    rows = cohort.loc[observed_train, ["group_id"]].copy()
    rows["fold"] = predictions.loc[observed_train, "oof_fold"].to_numpy()
    folds = pd.to_numeric(rows.fold, errors="raise").to_numpy(float)
    if not np.isfinite(folds).all() or (folds < 0).any() or not np.equal(folds, np.floor(folds)).all():
        raise ValueError("Observed train donors need valid held-out fold numbers")
    if len(np.unique(folds)) < 2 or rows.groupby("group_id").fold.nunique().max() != 1:
        raise ValueError("OOF folds must be multiple and family-disjoint")
    oof = pd.to_numeric(predictions[column], errors="raise").to_numpy(float)
    if not np.isfinite(oof[observed_train]).all():
        raise ValueError("Missing primary donor OOF predictions; no in-sample replacement allowed")
    return oof


def batches(indices, size):
    for start in range(0, len(indices), size):
        yield indices[start:start+size]


def fit_projection(x, train_rows, maps, dimensions=16, seed=20261115,
                   normalization="signal", batch_size=512):
    """Bounded-memory randomized PCA; parameters depend on training X only.

One global scale per map preserves within-map relative regional amplitudes.
The explicit unit_block ablation equalizes individual block SDs.
"""
    n, width = len(train_rows), x.shape[1]
    if n < 3 or len(maps) != width:
        raise ValueError("Projection requires training rows and one map per feature")
    mean, m2, seen = np.zeros(width), np.zeros(width), 0
    for rows in batches(np.asarray(train_rows), batch_size):
        a = np.asarray(x[rows], float)
        if not np.isfinite(a).all():
            raise ValueError("Nonfinite genetic features")
        bn, bm = len(a), a.mean(0)
        delta = bm-mean
        m2 += ((a-bm)**2).sum(0) + delta**2*seen*bn/(seen+bn)
        mean += delta*bn/(seen+bn); seen += bn
    variance = m2/seen
    keep = np.flatnonzero(variance > 1e-14)
    if not len(keep):
        raise ValueError("No variable genetic features in training")
    labels = np.asarray(maps, str)[keep]
    unique = np.unique(labels)
    scale = np.empty(len(keep))
    for label in unique:
        where = labels == label
        if normalization == "unit_block":
            sd = np.sqrt(variance[keep[where]])
            scale[where] = np.maximum(sd, max(np.quantile(sd, .1), 1e-12))*np.sqrt(where.sum()*len(unique))
        else:
            # Each map has unit aggregate training variance, without equalizing
            # weak and strong blocks inside that map.
            scale[where] = np.sqrt(variance[keep[where]].sum()*len(unique))
    center = mean[keep]
    if dimensions is None:
        return dict(keep=keep, mean=center, scale=scale, components=None,
                    normalization=normalization, n_train=n, input_width=width,
                    dimensions=len(keep), captured_training_variance=1.,
                    representation="all_variable_blocks_without_projection",storage_dtype="float32")
    rank = min(int(dimensions), len(keep), n-1)
    if rank < 1:
        raise ValueError("Projection dimension must be positive")
    if len(keep) <= rank:
        components = np.eye(len(keep))
        explained = 1.
    else:
        sketch_rank = min(rank+8, len(keep), n-1)
        omega = np.random.default_rng(seed).standard_normal((len(keep), sketch_rank))
        sketch = np.empty((n, sketch_rank))
        offset = 0
        for rows in batches(np.asarray(train_rows), batch_size):
            a = (np.asarray(x[rows], float)[:, keep]-center)/scale
            sketch[offset:offset+len(rows)] = a@omega; offset += len(rows)
        q, _ = np.linalg.qr(sketch, mode="reduced")
        del sketch
        small = np.zeros((sketch_rank, len(keep)))
        offset = 0
        for rows in batches(np.asarray(train_rows), batch_size):
            a = (np.asarray(x[rows], float)[:, keep]-center)/scale
            small += q[offset:offset+len(rows)].T@a; offset += len(rows)
        _, singular, vt = np.linalg.svd(small, full_matrices=False)
        components = vt[:rank].T
        total = float(np.sum(variance[keep]/scale**2)*n)
        explained = float(np.square(singular[:rank]).sum()/total) if total else np.nan
    return dict(keep=keep, mean=center, scale=scale, components=components,
                normalization=normalization, n_train=n, input_width=width,
                dimensions=components.shape[1], captured_training_variance=explained)


def transform_projection(x, state, batch_size=512, output_path=None):
    if output_path is None:
        out = np.empty((len(x), state["dimensions"]), dtype=state.get("storage_dtype","float64"))
    else:
        out = np.lib.format.open_memmap(output_path, mode="w+", dtype="float32",
                                       shape=(len(x), state["dimensions"]))
    for rows in batches(np.arange(len(x)), batch_size):
        a = np.asarray(x[rows], float)[:, state["keep"]]
        if not np.isfinite(a).all():
            raise ValueError("Nonfinite query genetic features")
        normalized = (a-state["mean"])/state["scale"]
        out[rows] = normalized if state["components"] is None else normalized@state["components"]
    if isinstance(out,np.memmap):out.flush()
    return out


class IndexedFeatures:
    def __init__(self, array, rows):
        self.array, self.rows = array, np.asarray(rows, int)
        self.shape = (len(rows), array.shape[1])
    def __len__(self):
        return len(self.rows)
    def __getitem__(self, rows):
        return self.array[self.rows[rows]]


def eligible_nearest(tree, donor_groups, query_z, query_groups, k, workers=1):
    """Exact k eligible donors, adaptively expanding past same-family members."""
    n = len(donor_groups)
    out_i = np.full((len(query_z), k), -1, dtype=np.int64)
    out_d = np.full((len(query_z), k), np.inf)
    pending = np.arange(len(query_z))
    requested = min(n, k+16)
    while len(pending):
        distance, index = tree.query(query_z[pending], k=requested, workers=workers)
        distance = np.asarray(distance).reshape(len(pending), requested)
        index = np.asarray(index).reshape(len(pending), requested)
        again = []
        for local, row in enumerate(pending):
            good = donor_groups[index[local]] != query_groups[row]
            ii, dd = index[local][good][:k], distance[local][good][:k]
            out_i[row, :len(ii)], out_d[row, :len(ii)] = ii, dd
            if len(ii) < k and requested < n:
                again.append(row)
        pending = np.asarray(again, int)
        requested = min(n, max(requested+1, requested*2))
    return out_i, out_d


def radius_from_calibration(tree, donor_groups, z, groups, k=64, quantile=.9,
                            batch_size=128, workers=1):
    distances = []
    for rows in batches(np.arange(len(z)), batch_size):
        _, d = eligible_nearest(tree, donor_groups, z[rows], groups[rows], k, workers)
        distances.extend(d[:, -1])
    d = np.asarray(distances)
    if not len(d) or not np.isfinite(d).all():
        raise ValueError("Calibration queries need at least radius-k eligible training donors")
    radius = float(np.quantile(d, quantile))
    # Identical coordinates still have well-defined full-radius membership.
    return max(radius, 1e-10)


def neighbourhoods(donor_z, query_z, donor_groups, query_groups, donor_y,
                   donor_residual, radius, k_values=(32,64,128), shrinkage=20.,
                   min_count=20, min_group_ess=15., batch_size=128, workers=1,search_k=None):
    tree = cKDTree(donor_z)
    _, donor_group_code = np.unique(donor_groups, return_inverse=True)
    group_sizes = pd.Series(donor_groups).value_counts().to_dict()
    n = len(query_z)
    correction = np.zeros((n, len(k_values)))
    support = []
    raw = np.empty(n)
    for rows in batches(np.arange(n), batch_size):
        ii, dd = eligible_nearest(tree, donor_groups, query_z[rows], query_groups[rows],
                                  max(max(k_values),search_k or 0), workers)
        neighborhoods = tree.query_ball_point(query_z[rows], radius, workers=workers)
        for local, row in enumerate(rows):
            members = np.asarray(neighborhoods[local], int)
            members = members[donor_groups[members] != query_groups[row]]
            _, sizes = np.unique(donor_group_code[members], return_counts=True)
            count = len(members)
            ess = float(count**2/np.square(sizes.astype(float)).sum()) if count else 0.
            enough = count >= min_count and ess >= min_group_ess
            raw[row] = float(donor_y[members].mean()) if count else float(donor_y.mean())
            support.append(dict(n_support=count, n_support_groups=len(sizes), support_group_ess=ess,
                n_eligible=len(donor_groups)-group_sizes.get(query_groups[row],0),
                nearest_distance=float(dd[local,0]), radius=radius, supported=enough))
            for j, k in enumerate(k_values):
                take = (ii[local,:k] >= 0) & (dd[local,:k] <= radius)
                indices, distances = ii[local,:k][take], dd[local,:k][take]
                if not enough or not len(indices):
                    continue
                bandwidth = max(float(distances[-1]), 1e-10)
                w = np.exp(-.5*(distances/bandwidth)**2)
                _, code = np.unique(donor_group_code[indices], return_inverse=True)
                group_w = np.bincount(code, weights=w)
                local_ess = float(w.sum()**2/np.square(group_w).sum())
                correction[row,j] = float(w@donor_residual[indices]/w.sum())*local_ess/(local_ess+shrinkage)
    return pd.DataFrame(support), correction, raw



def projected_neighbourhoods(donor_z, query_z, donor_groups, query_groups, donor_y,
                             donor_residual, radius, k_values, shrinkage, min_count, min_group_ess,
                             batch_size, workers, donor_ids, evidence_rows=(), evidence_neighbors=0):
    """Count all eligible radius members; retain closest donors with stable ID tie breaks."""
    tree = cKDTree(donor_z)
    requested = set(map(int, evidence_rows))
    donor_ids = np.asarray(donor_ids, str)
    levels, group_codes, group_sizes = np.unique(donor_groups, return_inverse=True, return_counts=True)
    family_sizes = dict(zip(levels, group_sizes))
    ranks = np.argsort(np.argsort(donor_ids, kind="stable"), kind="stable")
    width = max(max(k_values), evidence_neighbors)
    correction = np.zeros((len(query_z), len(k_values)))
    rows_out, evidence = [], {}
    raw = np.full(len(query_z), float(np.mean(donor_y)))
    for rows in batches(np.arange(len(query_z)), batch_size):
        balls = tree.query_ball_point(query_z[rows], radius, workers=workers, return_sorted=False)
        empty = []
        for row, ball in zip(rows, balls):
            members = np.asarray(ball, int)
            members = members[donor_groups[members] != query_groups[row]]
            _, sizes = np.unique(group_codes[members], return_counts=True)
            n = len(members)
            ess = n*n / np.square(sizes.astype(float)).sum() if n else 0.
            supported = n >= min_count and ess >= min_group_ess
            eligible = len(donor_z) - family_sizes.get(query_groups[row], 0)
            nearest = np.inf
            chosen, distances = np.empty(0, int), np.empty(0)
            if n:
                raw[row] = np.mean(donor_y[members])
                delta = donor_z[members] - query_z[row]
                squared = np.einsum("ij,ij->i", delta, delta)
                # Include all ties at the top-k boundary before ordering by donor ID.
                pool = np.arange(n)
                if n > width:
                    pool = np.flatnonzero(squared <= np.partition(squared,width-1)[width-1])
                order = pool[np.lexsort((ranks[members[pool]], squared[pool]))[:width]]
                chosen, distances = members[order], np.sqrt(squared[order])
                nearest = distances[0]
            elif eligible:
                empty.append(row)
            if row in requested:
                evidence[int(row)] = dict(indices=chosen.copy(), distances=distances.copy(),
                                          supported=supported, n_support=n)
            rows_out.append(dict(n_support=n, n_support_groups=len(sizes), support_group_ess=ess,
                                 n_eligible=eligible, nearest_distance=nearest, radius=radius, supported=supported))
            if supported:
                for j, k in enumerate(k_values):
                    ix, ds = chosen[:k], distances[:k]
                    weights = np.exp(-.5*(ds/max(float(ds[-1]),1e-10))**2)
                    _, codes = np.unique(group_codes[ix], return_inverse=True)
                    mass = np.bincount(codes, weights=weights)
                    local_ess = weights.sum()**2 / np.square(mass).sum()
                    correction[row,j] = (weights @ donor_residual[ix] / weights.sum()
                                         * local_ess / (local_ess+shrinkage))
        if empty:
            _, nearest = eligible_nearest(tree, donor_groups, query_z[empty], query_groups[empty], 1, workers)
            for row, d in zip(empty, nearest[:,0]):
                rows_out[row]["nearest_distance"] = d
    return pd.DataFrame(rows_out), correction, raw, evidence


def add_counts(support, k, alpha):
    """Keep all-radius support and actual residual borrowing distinguishable per person."""
    count = support.n_support.to_numpy(int)
    support["support_fraction"] = np.divide(count, support.n_eligible.to_numpy(float),
                                             out=np.zeros(len(count)), where=support.n_eligible.to_numpy()>0)
    support["n_donors_used"] = np.where(support.supported.to_numpy(bool) & (alpha != 0),
                                        np.minimum(count, int(k)), 0)


def full_neighbourhoods(donor_z, query_z, donor_groups, query_groups, donor_y,
                       donor_residual, radius, k_values=(32,64,128), shrinkage=20.,
                       min_count=20, min_group_ess=15., batch_size=64,
                       donor_batch=4096, donor_ids=None, query_ids=None,
                       evidence_rows=(), search_k=None, workers=1):
    """Full-space all-radius sufficient statistics and local residual corrections."""
    from block_neighbors import ExactBlockIndex
    from parallel_runtime import worker_count, ordered_map, row_slices
    bank=ExactBlockIndex(donor_z,donor_groups,donor_ids,donor_batch)
    bank.validate_queries(query_z,query_groups,query_ids)
    estimate=bank.estimate(batch_size,batch_size,max(k_values))
    count=worker_count(workers,max(1,len(query_z)//batch_size),
        memory_per_worker_gib=.5+estimate["bounded_scan_working_bytes"]/2**30,
        reserve_gib=(donor_z.nbytes+query_z.nbytes)/2**30)
    if count>1:
        slices=row_slices(len(query_z),count)
        requested=np.asarray(list(evidence_rows),int)
        arguments=((donor_z,query_z[a:b],donor_groups,np.asarray(query_groups)[a:b],
            donor_y,donor_residual,radius,k_values,shrinkage,min_count,min_group_ess,
            batch_size,donor_batch,donor_ids,None if query_ids is None else np.asarray(query_ids)[a:b],
            requested[(requested>=a)&(requested<b)]-a,search_k,1) for a,b in slices)
        print(f"[abm] full-block query processes={count}; BLAS threads/process=1",flush=True)
        pieces=list(ordered_map(full_neighbourhoods,arguments,count))
        evidence={a+row:value for (a,b),piece in zip(slices,pieces) for row,value in piece[3].items()}
        return (pd.concat([p[0] for p in pieces],ignore_index=True),
                np.concatenate([p[1] for p in pieces]),np.concatenate([p[2] for p in pieces]),evidence)
    n=len(query_z);correction=np.zeros((n,len(k_values)));raw=np.empty(n)
    support=[];evidence={};requested=set(map(int,evidence_rows))
    y=np.asarray(donor_y,float);residual=np.asarray(donor_residual,float)
    if not np.isfinite(residual).all():raise ValueError("Nonfinite OOF donor residuals")
    for rows,ii,dd,stats in bank.batches(query_z,query_groups,max(max(k_values),search_k or 0),
            radius=radius,query_ids=query_ids,batch_size=batch_size,donor_y=y):
        counts=stats["n_support"]
        raw[rows]=np.divide(stats.pop("outcome_sum"),counts,out=np.full(len(rows),y.mean()),where=counts>0)
        second=stats.pop("outcome_square")
        variance=np.divide(second,counts,out=np.full(len(rows),np.nan),where=counts>0)-raw[rows]**2
        stats["raw_neighbor_variance"]=np.maximum(variance,0.)
        stats["radius"]=np.full(len(rows),radius)
        enough=(counts>=min_count)&(stats["support_group_ess"]>=min_group_ess)
        stats["supported"]=enough
        support.append(pd.DataFrame(stats))
        for local,row in enumerate(rows):
            if row in requested:
                evidence[int(row)]={"indices":ii[local].copy(),"distances":dd[local].copy(),
                                    "supported":bool(enough[local]),"n_support":int(counts[local])}
            if not enough[local]:continue
            for j,k in enumerate(k_values):
                take=(ii[local,:k]>=0)&(dd[local,:k]<=radius)
                indices=ii[local,:k][take];distances=dd[local,:k][take]
                if not len(indices):continue
                weights=np.exp(-.5*(distances/max(float(distances[-1]),1e-10))**2)
                _,code=np.unique(np.asarray(donor_groups)[indices],return_inverse=True)
                mass=np.bincount(code,weights=weights)
                ess=float(weights.sum()**2/np.square(mass).sum())
                correction[row,j]=float(weights@residual[indices]/weights.sum())*ess/(ess+shrinkage)
    return pd.concat(support,ignore_index=True),correction,raw,evidence


def evidence_selection(cohort,maximum=200):
    """Outcome-blind stable examples; prefer TEST identities for inspection."""
    rows=np.flatnonzero(cohort.role.eq("test")) if "role" in cohort and cohort.role.eq("test").any() else np.arange(len(cohort))
    order=sorted(rows,key=lambda i:hashlib.sha256(str(cohort.id.iloc[i]).encode()).digest())
    return np.asarray(order[:max(0,int(maximum))],int)


def write_twin_evidence(path,evidence,cohort,donor_ids,donor_groups,donor_y,donor_residual,
                        baseline,model,neighbors=16):
    """Bounded local-only donor evidence. Never attaches query outcomes."""
    columns=["query_id","query_role","donor_rank","donor_id","donor_group","distance","radius",
             "within_radius","query_supported","n_support","prediction_k","query_baseline",
             "donor_outcome","donor_oof_prediction","donor_oof_residual","used_for_correction",
             "normalized_kernel_weight","local_family_ess","shrink_factor","alpha",
             "residual_contribution","evidence_is_truncated"]
    records=[];k=int(model["k"]);radius=float(model["radius"])
    for row,item in evidence.items():
        ii,dd=item["indices"],item["distances"]
        valid=(ii[:k]>=0)&(dd[:k]<=radius)
        chosen=ii[:k][valid];distance=dd[:k][valid]
        weights=np.exp(-.5*(distance/max(float(distance[-1]),1e-10))**2) if len(chosen) else np.array([])
        if len(weights):
            weights/=weights.sum();_,codes=np.unique(np.asarray(donor_groups)[chosen],return_inverse=True)
            ess=float(1/np.square(np.bincount(codes,weights=weights)).sum())
        else:ess=0.
        shrink=ess/(ess+float(model["shrinkage"])) if ess else 0.
        lookup=dict(zip(chosen,weights))
        for rank,(index,d) in enumerate(zip(ii[:max(0,int(neighbors))],dd[:max(0,int(neighbors))]),1):
            if index<0:continue
            used=bool(model["alpha"] != 0 and item["supported"] and index in lookup)
            weight=float(lookup.get(index,0.))
            records.append(dict(query_id=str(cohort.id.iloc[row]),query_role=str(cohort.role.iloc[row]),
                donor_rank=rank,donor_id=str(donor_ids[index]),donor_group=str(donor_groups[index]),
                distance=float(d),radius=radius,within_radius=bool(d<=radius),query_supported=item["supported"],
                n_support=item["n_support"],prediction_k=k,query_baseline=float(baseline[row]),
                donor_outcome=float(donor_y[index]),donor_oof_prediction=float(donor_y[index]-donor_residual[index]),
                donor_oof_residual=float(donor_residual[index]),used_for_correction=used,
                normalized_kernel_weight=weight,local_family_ess=ess,shrink_factor=shrink,
                alpha=float(model["alpha"]),residual_contribution=float(model["alpha"])*shrink*weight*float(donor_residual[index]) if used else 0.,
                evidence_is_truncated=int(item["n_support"])>int(neighbors)))
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(records,columns=columns).to_csv(path,sep="\t",index=False)
    return len(records)


def clipped(prediction, binary):
    return np.clip(prediction, 0., 1.) if binary else np.asarray(prediction, float)


def group_mean_se(values, groups):
    values = np.asarray(values, float)
    _, code = np.unique(np.asarray(groups,str), return_inverse=True)
    n = np.bincount(code); sums = np.bincount(code, weights=values)
    if len(n) < 2:
        return np.inf
    mean = values.mean()
    return float(np.sqrt(len(n)/(len(n)-1)*np.square(sums-n*mean).sum()/len(values)**2))


def select_correction(y, baseline, corrections, groups, k_values,
                      lambdas=(0.,.25,.5,1.), binary=False):
    observed = np.isfinite(y)
    if observed.sum() < 20:
        raise ValueError("At least 20 observed tune outcomes required")
    y, base = y[observed], baseline[observed]
    errors0 = np.square(y-clipped(base,binary))
    candidates = []
    for j,k in enumerate(k_values):
        for weight in lambdas:
            loss = np.square(y-clipped(base+weight*corrections[observed,j],binary))
            delta = loss-errors0
            candidates.append(dict(k=int(k), weight=float(weight), mse=float(loss.mean()),
                delta=float(delta.mean()), delta_se=group_mean_se(delta,np.asarray(groups)[observed])))
    # Conservative one-SE rule on paired tune loss: when evidence is weak,
    # retain lambda=0. No candidate is chosen from test or calibration outcomes.
    admissible = [x for x in candidates if x["weight"]>0 and x["delta"]+x["delta_se"]<0]
    chosen = min(admissible,key=lambda x:(x["mse"],x["weight"],x["k"])) if admissible else dict(
        k=int(k_values[len(k_values)//2]),weight=0.,mse=float(errors0.mean()),delta=0.,delta_se=0.)
    return chosen, candidates


def error_features(cohort, baseline, support):
    columns = [c for c in ["age","sex"]+[f"PC{i}" for i in range(1,21)] if c in cohort]
    base = [np.asarray(baseline,float), np.square(np.asarray(baseline,float))]
    base += [pd.to_numeric(cohort[c],errors="coerce").to_numpy(float) for c in columns]
    # HistGradientBoosting handles missing covariates. Ancestry coordinates are
    # already represented by PCs; no test labels define categorical encodings.
    x = np.column_stack(base)
    distance = support.nearest_distance.to_numpy(float)
    distance[~np.isfinite(distance)] = np.nan
    extra = np.column_stack([np.log1p(support.n_support),np.log1p(support.support_group_ess),distance])
    return x, np.column_stack([x,extra]), ["baseline","baseline_squared"]+columns


def freeze_masks(scores, calibration, coverages=COVERAGES, larger=False, eligible=None):
    scores = np.asarray(scores,float)
    valid = np.ones(len(scores),bool) if eligible is None else np.asarray(eligible,bool)
    ranked = scores if larger else -scores
    ranked = np.where(valid & np.isfinite(ranked),ranked,-np.inf)
    cal = ranked[np.asarray(calibration,bool)]
    thresholds, masks = {}, {}
    for coverage in coverages:
        ix = max(0,len(cal)-int(np.ceil(coverage*len(cal))))
        threshold = float(np.sort(cal)[ix])
        label = f"{coverage:.2f}".replace(".","p")
        thresholds[label] = threshold if np.isfinite(threshold) else None
        masks[label] = valid & np.isfinite(scores) & (ranked>=threshold)
    return thresholds,masks


def predict_frozen(bundle,cohort_query,base_predictions,features=None,*,workers=1,batch_size=128,
                   evidence_path=None,evidence_queries=200,evidence_neighbors=16,feature_manifest=None):
    """Predict new rows with a trusted local fitted bundle; no outcomes argument.

New regional bundles require the training producer's frozen ordered manifest.
Query batching cannot change the projection, radius, heads or threshold values.
"""
    state=joblib.load(bundle) if isinstance(bundle,(str,Path)) else bundle
    if state.get("schema")!="grid.abm2.bundle.1":
        raise ValueError("Unknown ABM bundle")
    model=state["model"]
    query=cohort_query.copy().reset_index(drop=True)
    if "role" not in query:query["role"]="query"
    required=["id","group_id","ancestry"]
    if not set(required)<=set(query) or query[required].isna().any().any() or query.id.duplicated().any():
        raise ValueError("New queries require unique IDs and explicit groups/ancestry")
    donor_lookup=dict(zip(state["donor_ids"],state["donor_groups"]))
    for identity,group in zip(query.id.astype(str),query.group_id.astype(str)):
        if identity in donor_lookup and donor_lookup[identity]!=group:
            raise ValueError("A donor/query ID has conflicting family identity")
    base=align(base_predictions,query.id,"new primary predictions")
    prediction=pd.to_numeric(base[PRIMARY+".prediction"],errors="raise").to_numpy(float)
    binary=model["task"]=="binary"
    if not np.isfinite(prediction).all() or (binary and ((prediction<0)|(prediction>1)).any()):
        raise ValueError("Invalid new primary predictions")
    if model["geometry"]=="global_prs":
        x=base[RAW_COLUMNS].to_numpy(float)
    else:
        if features is None or features.shape!=(len(query),state["projection"]["input_width"]):
            raise ValueError("Regional query features must match the frozen input shape/order")
        if state.get("feature_schema") is not None:
            if feature_manifest is None:
                raise ValueError("This frozen regional model requires a feature_manifest")
            if regional_feature_schema(feature_manifest,features.shape[1])!=state["feature_schema"]:
                raise ValueError("Regional feature identities/order differ from the frozen model")
        x=features
    z=transform_projection(x,state["projection"])
    if model["geometry"]=="consensus":
        support,correction,raw,evidence=full_neighbourhoods(state["donor_z"],z,state["donor_groups"],query.group_id.to_numpy(str),
            state["donor_y"],state["donor_residual"],model["radius"],(model["k"],),model["shrinkage"],
            model["min_count"],model["min_group_ess"],batch_size,model.get("donor_batch",4096),
            state["donor_ids"],query.id.to_numpy(str),evidence_selection(query,evidence_queries) if evidence_path else (),
            model["search_k"],workers=workers)
    elif model.get("projected_search") == "exact_radius_id_ties":
        support,correction,raw,evidence=projected_neighbourhoods(state["donor_z"],z,state["donor_groups"],
            query.group_id.to_numpy(str),state["donor_y"],state["donor_residual"],model["radius"],(model["k"],),
            model["shrinkage"],model["min_count"],model["min_group_ess"],batch_size,workers,state["donor_ids"],
            evidence_selection(query,evidence_queries) if evidence_path else (), evidence_neighbors)
    else:
        evidence={}
        support,correction,raw=neighbourhoods(state["donor_z"],z,state["donor_groups"],query.group_id.to_numpy(str),
            state["donor_y"],state["donor_residual"],model["radius"],(model["k"],),model["shrinkage"],
            model["min_count"],model["min_group_ess"],batch_size,workers,model["search_k"])
    if evidence_path and (model["geometry"]=="consensus" or model.get("projected_search")):
        write_twin_evidence(evidence_path,evidence,query,state["donor_ids"],state["donor_groups"],
            state["donor_y"],state["donor_residual"],prediction,model,evidence_neighbors)
    add_counts(support,model["k"],model["alpha"])
    result=query[["id","group_id","ancestry","role"]].copy()
    result["abm.prediction"]=clipped(prediction+model["alpha"]*correction[:,0],binary)
    result["abm.global_residual.prediction"]=clipped(prediction+model["global_residual_alpha"]*model["global_residual_mean"],binary)
    result["abm.raw_radius_mean.prediction"]=clipped(raw,binary)
    for col in support:result[col]=support[col].to_numpy()
    for col in model["error_features"][2:]:
        if col not in query:query[col]=np.nan
    basic,extended,_=error_features(query,prediction,support)
    scores={"support":support.support_group_ess.to_numpy(float),"twin_count":support.n_support.to_numpy(float)}
    for label,matrix in (("risk_error",basic),("support_error",extended)):
        scores[label]=np.maximum(state["error_heads"][label].predict(matrix),0.)
        result[label+".expected_loss"]=scores[label]
    for label,thresholds in model["policies"].items():
        ranked=scores[label] if label in {"support","twin_count"} else -scores[label]
        eligible=support.supported.to_numpy(bool) if label=="support" else np.ones(len(query),bool)
        if label=="twin_count":eligible=support.n_support.to_numpy()>0
        for level,threshold in thresholds.items():
            result[f"selected_{label}_{level}"]=eligible&np.isfinite(ranked)&(True if threshold is None else ranked>=threshold)
    return result


def fit_abm(cohort, outcomes, predictions, output_dir, *, trait, features=None,
            feature_maps=None, geometry="global_prs", dimensions=16, k_values=(32,64,128),
            radius_k=64, radius_quantile=.9, seed=20261115, workers=1,
            batch_size=128, min_count=20, min_group_ess=15., shrinkage=20.,
            donor_batch=4096,evidence_queries=200,evidence_neighbors=16,feature_manifest=None):
    cohort = cohort.copy().reset_index(drop=True)
    check_roster(cohort)
    predictions = align(predictions,cohort.id,"benchmark predictions")
    baseline = pd.to_numeric(predictions[PRIMARY+".prediction"],errors="raise").to_numpy(float)
    if not np.isfinite(baseline).all():
        raise ValueError("Primary predictions must cover every cohort role")
    y = development_outcomes(outcomes,cohort)
    train = cohort.role.eq("train").to_numpy()
    donor = train & np.isfinite(y)
    tune = cohort.role.eq("tune").to_numpy()
    calibration = cohort.role.eq("calibration").to_numpy()
    binary = trait == "t2dm"
    if donor.sum()<max(radius_k,min_count):
        raise ValueError("Too few observed training donors for the declared radius/support rule")
    if binary and not np.isin(y[np.isfinite(y)],[0.,1.]).all():
        raise ValueError("T2DM is a binary baseline endpoint")
    oof = check_oof(cohort,predictions,donor)
    if binary and (np.any((baseline<0)|(baseline>1)) or np.any((oof[donor]<0)|(oof[donor]>1))):
        raise ValueError("Binary baseline and OOF predictions must be probabilities")
    if geometry == "global_prs":
        if not set(RAW_COLUMNS)<=set(predictions):
            raise ValueError("Four raw ancestry PRS columns required for global_prs geometry")
        x = predictions[RAW_COLUMNS].to_numpy(float)
        maps = np.asarray(["EUR","AFR","EAS","SAS"])
        dimensions = 4
    else:
        if features is None or feature_maps is None or len(features)!=len(cohort):
            raise ValueError("Regional geometry requires aligned features and map manifest")
        x,maps = features,np.asarray(feature_maps,str)
    feature_schema=(regional_feature_schema(feature_manifest,x.shape[1])
                    if geometry!="global_prs" and feature_manifest is not None else None)
    normalization = "unit_block" if geometry == "unit_block" else "signal"
    if geometry not in {"consensus","regional","global_prs","unit_block"}:
        raise ValueError("Unknown declared ABM geometry")
    out = Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    projection = fit_projection(x,np.flatnonzero(train),maps,None if geometry=="consensus" else dimensions,seed,normalization)
    z = transform_projection(x,projection,output_path=out/"geometry.npy" if geometry=="consensus" else None)
    dg = cohort.loc[donor,"group_id"].to_numpy(str)
    groups = cohort.group_id.to_numpy(str)
    dz = z[donor]
    evidence={};work_estimate={}
    if geometry=="consensus":
        from block_neighbors import ExactBlockIndex
        bank=ExactBlockIndex(dz,dg,cohort.loc[donor,"id"].to_numpy(str),donor_batch)
        work_estimate=bank.estimate(len(z),batch_size,max(k_values))
        print(f"[abm] full-block exact CPU search: {work_estimate['distance_pairs']:,} pairs; "
              f"~{work_estimate['approximate_distance_flops']/1e12:.3f} TFLOP arithmetic; "
              f"scan buffers ~{work_estimate['bounded_scan_working_bytes']/2**20:.1f} MiB (excludes features/BLAS)",flush=True)
        _,reference=bank.nearest(z[calibration],groups[calibration],radius_k,
            query_ids=cohort.loc[calibration,"id"].to_numpy(str),batch_size=batch_size,workers=workers)
        if not len(reference) or not np.isfinite(reference[:,-1]).all():
            raise ValueError("Calibration queries need at least radius-k eligible training donors")
        radius=max(float(np.quantile(reference[:,-1],radius_quantile)),1e-10)
    else:
        tree = cKDTree(dz)
        radius = radius_from_calibration(tree,dg,z[calibration],groups[calibration],radius_k,
                                        radius_quantile,batch_size,workers)
    print(f"[abm] donor={donor.sum():,}; all queries={len(z):,}; frozen dimensions={z.shape[1]}; "
          f"radius={radius:.6g}; worst-case search pairs={int(donor.sum())*len(z):,}; "
          f"embedding storage={z.nbytes/2**20:.1f} MiB; no dense pair matrix",flush=True)
    if geometry=="consensus":
        support,correction,raw,evidence = full_neighbourhoods(dz,z,dg,groups,y[donor],y[donor]-oof[donor],radius,
            k_values,shrinkage,min_count,min_group_ess,batch_size,donor_batch,
            cohort.loc[donor,"id"].to_numpy(str),cohort.id.to_numpy(str),evidence_selection(cohort,evidence_queries),workers=workers)
    else:
        support,correction,raw,evidence = projected_neighbourhoods(dz,z,dg,groups,y[donor],y[donor]-oof[donor],
            radius,k_values,shrinkage,min_count,min_group_ess,batch_size,workers,cohort.loc[donor,"id"].to_numpy(str),
            evidence_selection(cohort,evidence_queries),evidence_neighbors)
    chosen,candidates = select_correction(y[tune],baseline[tune],correction[tune],groups[tune],k_values,binary=binary)
    add_counts(support,chosen["k"],chosen["weight"])
    j = list(k_values).index(chosen["k"])
    prediction = clipped(baseline+chosen["weight"]*correction[:,j],binary)
    global_residual=float(np.mean(y[donor]-oof[donor]))
    global_chosen,_=select_correction(y[tune],baseline[tune],np.full((int(tune.sum()),1),global_residual),
                                     groups[tune],(0,),binary=binary)
    query = ~train
    result = cohort.loc[query,["id","group_id","ancestry","role"]].reset_index(drop=True)
    result["abm.prediction"] = prediction[query]
    result["abm.global_residual.prediction"] = clipped(baseline[query]+global_chosen["weight"]*global_residual,binary)
    result["abm.raw_radius_mean.prediction"] = clipped(raw[query],binary)
    for col in support:
        result[col] = support.loc[query,col].to_numpy()
    # Train error targets are genuinely OOF baseline losses. The primary
    # baseline predictor is swapped to its OOF values on training donors.
    error_base = baseline.copy();error_base[donor]=oof[donor]
    basic,extended,error_names = error_features(cohort,error_base,support)
    heads,policy = {},{}
    target = np.square(y[donor]-oof[donor])
    for label,matrix in (("risk_error",basic),("support_error",extended)):
        head = HistGradientBoostingRegressor(max_iter=80,max_leaf_nodes=7,l2_regularization=10.,
            min_samples_leaf=max(5,min(200,int(donor.sum()/20))),random_state=seed,early_stopping=False)
        head.fit(matrix[donor],target)
        expected = np.maximum(head.predict(matrix),0.)
        heads[label]=head
        result[label+".expected_loss"] = expected[query]
        thresholds,masks = freeze_masks(expected,calibration)
        policy[label]=thresholds
        for level,mask in masks.items():
            result[f"selected_{label}_{level}"] = mask[query]
    thresholds,masks = freeze_masks(support.support_group_ess,calibration,larger=True,
                                    eligible=support.supported.to_numpy(bool))
    policy["support"]=thresholds
    for level,mask in masks.items():
        result["selected_support_"+level] = mask[query]
    thresholds,masks = freeze_masks(support.n_support,calibration,coverages=(.2,),larger=True,
                                    eligible=support.n_support.to_numpy()>0)
    policy["twin_count"] = thresholds
    result["selected_twin_count_0p20"] = masks["0p20"][query]
    result.to_csv(out/"abm.predictions.tsv.gz",sep="\t",index=False)
    pd.DataFrame(candidates).to_csv(out/"tuning.tsv",sep="\t",index=False)
    model = dict(version="grid.abm2.2",trait=trait,primary=PRIMARY,task="binary" if binary else "continuous",
        geometry=geometry,projection_dimensions=z.shape[1],captured_training_variance=projection["captured_training_variance"],
        projection_fit_role="train_X_only",normalization=normalization,
        projection_quality=("All variable normalized LD-block coordinates; no PCA; float32 frozen coordinates and float64 exact search"
            if geometry=="consensus" else "Captured normalized training feature variance; randomized approximation, not full original-space distance preservation"),
        residual_contrast="OOF primary-residual borrowing compared with primary alone, global OOF-residual mean and all-radius raw outcome mean",
        neighbour_definition=("all eligible donors inside radius in full frozen regional score space"
            if geometry=="consensus" else "all eligible donors inside radius in frozen projected genetic-score space"),
        family_ess_definition="N_squared_over_sum_family_counts_squared; concentration diagnostic, not independence guarantee",
        donor_count=int(donor.sum()),radius=radius,radius_k=radius_k,radius_quantile=radius_quantile,
        radius_fit_role="calibration_X_only",alpha=float(chosen["weight"]),k=int(chosen["k"]),
        search_k=max(k_values),donor_batch=int(donor_batch),search_work=work_estimate,
        global_residual_mean=global_residual,global_residual_alpha=float(global_chosen["weight"]),
        raw_ablation="all-radius unadjusted raw donor mean; global raw donor mean if zero neighbours",
        correction_selection="tune paired loss improvement beyond one cluster standard error; zero is fallback",
        residual_source=PRIMARY+".oof_prediction",min_count=min_count,min_group_ess=min_group_ess,shrinkage=shrinkage,
        support_threshold_role="calibration_X_only",policies=policy,error_target="train group-OOF primary squared loss",
        error_features=error_names,support_error_extra_features=["log1p(n_support)","log1p(support_group_ess)","nearest_distance"],
        test_outcomes_used=False,calibration_outcomes_used=False,seed=seed,
        statement="Support and expected-error scores are empirical diagnostics, not individual error guarantees")
    model["feature_identity_validation"]=("fixed named global ancestry score columns" if geometry=="global_prs"
        else "ordered regional manifest required at prediction" if feature_schema is not None
        else "input-width contract; no regional manifest supplied")
    if geometry!="consensus":
        model["projected_search"]="exact_radius_id_ties"
    evidence_count=write_twin_evidence(out/"twins.tsv.gz",evidence,cohort,cohort.loc[donor,"id"].to_numpy(str),
            dg,y[donor],y[donor]-oof[donor],baseline,model,evidence_neighbors)
    model["twin_evidence"]={"file":"twins.tsv.gz","queries":len(evidence),"rows":evidence_count,
            "max_neighbors_per_query":int(evidence_neighbors),"selection":"stable identity hash, TEST preferred, no query Y",
            "scope":"bounded donor examples; all-radius n_support is computed independently of this export"}
    # model.json stays local with the modelling snapshots; public evaluation
    # exports aggregate metadata only and never embeds participant identities.
    (out/"model.json").write_text(json.dumps(model,indent=2,allow_nan=False)+"\n")
    bundle=dict(schema="grid.abm2.bundle.1",projection=projection,error_heads=heads,model=model,
        donor_z=dz,donor_y=y[donor],donor_residual=y[donor]-oof[donor],
        donor_groups=dg,donor_ids=cohort.loc[donor,"id"].to_numpy(str),feature_schema=feature_schema)
    joblib.dump(bundle,out/"abm.model.joblib",compress=3)
    print(f"[abm] frozen lambda={chosen['weight']}, k={chosen['k']}; {out/'abm.predictions.tsv.gz'}",flush=True)
    return result,model


def main(argv=None):
    argv=list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0]=="predict":
        p=argparse.ArgumentParser(description="Apply a trusted frozen ABM bundle without query outcomes")
        for key in ("bundle","cohort","predictions","output"):p.add_argument("--"+key,required=True)
        p.add_argument("--features-npy");p.add_argument("--features-ids");p.add_argument("--feature-manifest")
        p.add_argument("--workers",type=int,default=1);p.add_argument("--batch-size",type=int,default=128)
        p.add_argument("--evidence-output");p.add_argument("--evidence-queries",type=int,default=200)
        p.add_argument("--evidence-neighbors",type=int,default=16)
        a=p.parse_args(argv[1:]);query=table(a.cohort);features=None
        if a.features_npy:
            if not a.features_ids:p.error("--features-ids is required with regional features")
            array=np.load(a.features_npy,mmap_mode="r",allow_pickle=False);ids=table(a.features_ids)
            if len(ids)!=len(array) or ids.id.duplicated().any():raise ValueError("Invalid regional sample manifest")
            rows=pd.Series(np.arange(len(ids)),index=ids.id).reindex(query.id)
            if rows.isna().any():raise ValueError("Regional features do not cover new queries")
            features=IndexedFeatures(array,rows.to_numpy(int))
        result=predict_frozen(a.bundle,query,table(a.predictions),features,workers=a.workers,batch_size=a.batch_size,
            evidence_path=a.evidence_output,evidence_queries=a.evidence_queries,evidence_neighbors=a.evidence_neighbors,
            feature_manifest=a.feature_manifest)
        output=Path(a.output);output.parent.mkdir(parents=True,exist_ok=True)
        result.to_csv(output,sep="\t",index=False)
        return
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort",required=True);parser.add_argument("--outcomes",required=True)
    parser.add_argument("--predictions",required=True);parser.add_argument("--trait",required=True)
    parser.add_argument("--output-dir",required=True)
    parser.add_argument("--features-npy");parser.add_argument("--features-ids");parser.add_argument("--feature-manifest")
    parser.add_argument("--geometry",choices=["consensus","regional","global_prs","unit_block"],default="consensus")
    parser.add_argument("--dimensions",type=int,default=16);parser.add_argument("--k-values",default="32,64,128")
    parser.add_argument("--radius-k",type=int,default=64);parser.add_argument("--radius-quantile",type=float,default=.9)
    parser.add_argument("--seed",type=int,default=20261115);parser.add_argument("--workers",type=int,default=1)
    parser.add_argument("--batch-size",type=int,default=128)
    parser.add_argument("--donor-batch",type=int,default=4096)
    parser.add_argument("--evidence-queries",type=int,default=200);parser.add_argument("--evidence-neighbors",type=int,default=16)
    args=parser.parse_args(argv)
    cohort=table(args.cohort)
    features=maps=manifest=None
    if args.geometry!="global_prs":
        if not all([args.features_npy,args.features_ids,args.feature_manifest]):
            parser.error("Regional geometry needs all three feature inputs; explicitly choose global_prs otherwise")
        array=np.load(args.features_npy,mmap_mode="r",allow_pickle=False)
        ids=table(args.features_ids)
        if len(ids)!=len(array) or ids.id.duplicated().any():
            raise ValueError("Regional sample manifest does not match array")
        rows=pd.Series(np.arange(len(ids)),index=ids.id).reindex(cohort.id)
        if rows.isna().any():
            raise ValueError("Regional features do not cover cohort")
        features=IndexedFeatures(array,rows.to_numpy(int))
        manifest=table(args.feature_manifest)
        column="map" if "map" in manifest else "population"
        maps=manifest[column].astype(str).to_numpy()
    k_values=tuple(sorted(set(int(x) for x in args.k_values.split(","))))
    if not k_values or min(k_values)<1 or args.radius_k<1 or not 0<=args.radius_quantile<=1:
        parser.error("Invalid neighbour configuration")
    fit_abm(cohort,table(args.outcomes),table(args.predictions),args.output_dir,trait=args.trait,
        features=features,feature_maps=maps,geometry=args.geometry,dimensions=args.dimensions,
        k_values=k_values,radius_k=args.radius_k,radius_quantile=args.radius_quantile,
        seed=args.seed,workers=args.workers,batch_size=args.batch_size,donor_batch=args.donor_batch,
        evidence_queries=args.evidence_queries,evidence_neighbors=args.evidence_neighbors,feature_manifest=manifest)


if __name__=="__main__":
    main()
