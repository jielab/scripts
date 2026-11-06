"""Outcome-free, fixed-radius genetic-profile support for ABM.

``n_twins`` counts *every* eligible donor in a frozen radius.  It is not a
nearest-neighbour hyperparameter.  ``k_used`` only caps the donors retained for
local summaries; their individual and family/group ESS are reported separately.
These are neighbours in the caller's declared genetic representation, not
biological twins.  Fit every representation/scaler and radius on development
data, then freeze them before calling this module on an unseen person.

The NumPy and optional CUDA implementations run the same blockwise exhaustive
search in float64.  CUDA accelerates distances/search; it does not change the
scientific definition.  This reference implementation has no approximate
nearest-neighbour index and makes no unmeasured speed claim.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence
import time

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class TwinResult:
    table: pd.DataFrame
    # Indices refer to the caller's original donor array. Padding is -1 / inf.
    neighbor_indices: np.ndarray
    neighbor_distances: np.ndarray


def _strings(values: Sequence, n: int, name: str, *, unique: bool) -> np.ndarray:
    raw = np.asarray(values, dtype=object)
    if raw.ndim != 1 or len(raw) != n:
        raise ValueError(f"{name} must be a one-dimensional array of length {n}")
    if any(pd.isna(v) for v in raw):
        raise ValueError(f"{name} cannot contain missing identifiers")
    out = np.asarray([str(v) for v in raw], dtype=str)
    if any(not v.strip() for v in out):
        raise ValueError(f"{name} cannot contain empty identifiers")
    if unique and len(set(out.tolist())) != n:
        raise ValueError(f"{name} must be unique after conversion to strings")
    return out


def _validate_inputs(donor_z, query_z, donor_ids, query_ids, donor_groups,
                     query_groups, k_used, query_batch, donor_batch):
    # Do not make a float64 copy of an entire genome-wide matrix. Conversion
    # happens per block; memmaps and float32 storage remain usable.
    dz, qz = np.asarray(donor_z), np.asarray(query_z)
    if dz.ndim != 2 or qz.ndim != 2 or dz.shape[1] != qz.shape[1] or dz.shape[1] == 0:
        raise ValueError("donor_z/query_z must be 2D with the same positive feature count")
    if not np.issubdtype(dz.dtype, np.number) or not np.issubdtype(qz.dtype, np.number):
        raise ValueError("geometry arrays must be numeric")
    if np.iscomplexobj(dz) or np.iscomplexobj(qz):
        raise ValueError("geometry arrays must be real-valued")
    # Validate in bounded blocks rather than allocate a whole-matrix mask.
    for label, values in (("donor_z", dz), ("query_z", qz)):
        for start in range(0, len(values), 4096):
            if not np.isfinite(values[start:start + 4096]).all():
                raise ValueError(f"{label} contains non-finite values; impute using build-only parameters")
    for name, val in (("k_used", k_used), ("query_batch", query_batch), ("donor_batch", donor_batch)):
        if isinstance(val, bool) or not isinstance(val, (int, np.integer)) or val < 1:
            raise ValueError(f"{name} must be a positive integer")
    di = _strings(donor_ids, len(dz), "donor_ids", unique=True)
    qi = _strings(query_ids, len(qz), "query_ids", unique=True)
    dg = _strings(donor_groups, len(dz), "donor_groups", unique=False)
    qg = _strings(query_groups, len(qz), "query_groups", unique=False)
    return dz, qz, di, qi, dg, qg


def _backend(device: str):
    if device not in {"auto", "cpu", "cuda"}:
        raise ValueError("device must be 'auto', 'cpu', or 'cuda'")
    if device == "cpu":
        return "cpu", None
    try:
        import torch
    except ImportError:
        if device == "cuda":
            raise RuntimeError("device='cuda' requested but PyTorch is not installed") from None
        return "cpu", None
    if torch.cuda.is_available():
        return "cuda", torch
    if device == "cuda":
        raise RuntimeError("device='cuda' requested but no CUDA device is available")
    return "cpu", None


def _scan(dz, qz, di, qi, dg, qg, *, radius, k, device, query_batch, donor_batch):
    """Exhaustive bounded-memory scan; an infinite radius is internal only."""
    used_device, torch = _backend(device)
    nd, nq = len(dz), len(qz)
    # Canonical lexical ID order defines exact-distance tie breaking, regardless
    # of how the donor file was ordered or how it is divided into blocks.
    canonical = np.argsort(di, kind="stable")
    all_ids = {v: i for i, v in enumerate(sorted(set(di.tolist() + qi.tolist())))}
    all_groups = {v: i for i, v in enumerate(sorted(set(dg.tolist() + qg.tolist())))}
    did = np.asarray([all_ids[v] for v in di], dtype=np.int64)
    qid = np.asarray([all_ids[v] for v in qi], dtype=np.int64)
    dgroup = np.asarray([all_groups[v] for v in dg], dtype=np.int64)
    qgroup = np.asarray([all_groups[v] for v in qg], dtype=np.int64)
    result_idx = np.full((nq, k), -1, dtype=np.int64)
    result_d2 = np.full((nq, k), np.inf, dtype=np.float64)
    eligible_total = np.zeros(nq, dtype=np.int64)
    in_radius_total = np.zeros(nq, dtype=np.int64)
    nearest_d2 = np.full(nq, np.inf, dtype=np.float64)
    radius2 = radius * radius
    resident = resident_norm = None
    if used_device == "cuda" and nd:
        free, _ = torch.cuda.mem_get_info()
        needed = nd * dz.shape[1] * 8
        # Leave room for queries, distances, stable-sort workspaces and the UI.
        if needed < free * .65:
            resident = torch.empty((nd, dz.shape[1]), dtype=torch.float64, device="cuda")
            resident_norm = torch.empty(nd, dtype=torch.float64, device="cuda")
            for ds in range(0, nd, donor_batch):
                ix = canonical[ds:ds + donor_batch]
                resident[ds:ds + len(ix)] = torch.as_tensor(np.asarray(dz[ix], dtype=np.float64), device="cuda")
                resident_norm[ds:ds + len(ix)] = (resident[ds:ds + len(ix)] ** 2).sum(dim=1)
        tdid = torch.as_tensor(did[canonical], device="cuda")
        tdgroup = torch.as_tensor(dgroup[canonical], device="cuda")
    last_progress = time.monotonic()

    for qs in range(0, nq, query_batch):
        qe = min(qs + query_batch, nq)
        # A bounded writable copy also avoids PyTorch's read-only memmap warning.
        q = np.array(qz[qs:qe], dtype=np.float64, copy=True)
        m = len(q)
        best_i = np.full((m, k), -1, dtype=np.int64)
        best_d = np.full((m, k), np.inf, dtype=np.float64)
        counts = np.zeros(m, dtype=np.int64)
        neligible = np.zeros(m, dtype=np.int64)
        near = np.full(m, np.inf, dtype=np.float64)
        if used_device == "cuda":
            tq = torch.as_tensor(q, dtype=torch.float64, device="cuda")
            tqnorm = (tq * tq).sum(dim=1)
            tqid = torch.as_tensor(qid[qs:qe], device="cuda")
            tqgroup = torch.as_tensor(qgroup[qs:qe], device="cuda")
            tbest_i = torch.full((m, k), -1, dtype=torch.int64, device="cuda")
            tbest_d = torch.full((m, k), float("inf"), dtype=torch.float64, device="cuda")
            tcounts = torch.zeros(m, dtype=torch.int64, device="cuda")
            teligible = torch.zeros(m, dtype=torch.int64, device="cuda")
            tnear = torch.full((m,), float("inf"), dtype=torch.float64, device="cuda")
        else:
            qnorm = np.einsum("ij,ij->i", q, q)

        for ds in range(0, nd, donor_batch):
            ix = canonical[ds:ds + donor_batch]
            if used_device == "cuda":
                td = (resident[ds:ds + len(ix)] if resident is not None else
                      torch.as_tensor(np.asarray(dz[ix], dtype=np.float64), dtype=torch.float64, device="cuda"))
                tdnorm = resident_norm[ds:ds + len(ix)] if resident is not None else (td * td).sum(dim=1)
                distance2 = (tqnorm[:, None] + tdnorm[None, :] - 2 * (tq @ td.T)).clamp_min_(0)
                if not bool(torch.isfinite(distance2).all().item()):
                    raise ValueError("distance overflow: check the frozen representation/scaling")
                # Dot products can cancel near a radius boundary. Refine only
                # those few pairs with the same direct float64 CPU expression
                # used by the CPU path, then compare to the unexpanded radius.
                if np.isfinite(radius2):
                    tolerance = 64 * np.finfo(np.float64).eps * (tqnorm[:, None] + tdnorm[None, :] + radius2 + 1)
                    edge = torch.nonzero(torch.abs(distance2 - radius2) <= tolerance, as_tuple=False)
                    if edge.numel():
                        edge_np = edge.detach().cpu().numpy()
                        delta = q[edge_np[:, 0]] - np.asarray(dz[ix[edge_np[:, 1]]], dtype=np.float64)
                        exact = np.einsum("ij,ij->i", delta, delta)
                        distance2[edge[:, 0], edge[:, 1]] = torch.as_tensor(exact, device="cuda")
                ok = ((tqid[:, None] != tdid[ds:ds + len(ix)][None, :]) &
                      (tqgroup[:, None] != tdgroup[ds:ds + len(ix)][None, :]))
                teligible += ok.sum(dim=1)
                tnear = torch.minimum(tnear, torch.where(ok, distance2, float("inf")).amin(dim=1))
                inside = ok & (distance2 <= radius2)
                tcounts += inside.sum(dim=1)
                masked = torch.where(inside, distance2, float("inf"))
                order = torch.argsort(masked, dim=1, stable=True)[:, :min(k, len(ix))]
                chunk_d = torch.gather(masked, 1, order)
                chunk_i = torch.as_tensor(ix, dtype=torch.int64, device="cuda")[order]
                merged_d = torch.cat((tbest_d, chunk_d), dim=1)
                merged_i = torch.cat((tbest_i, chunk_i), dim=1)
                pick = torch.argsort(merged_d, dim=1, stable=True)[:, :k]
                tbest_d = torch.gather(merged_d, 1, pick)
                tbest_i = torch.gather(merged_i, 1, pick)
            else:
                d = np.asarray(dz[ix], dtype=np.float64)
                dnorm = np.einsum("ij,ij->i", d, d)
                distance2 = np.maximum(qnorm[:, None] + dnorm[None, :] - 2 * (q @ d.T), 0)
                if not np.isfinite(distance2).all():
                    raise ValueError("distance overflow: check the frozen representation/scaling")
                if np.isfinite(radius2):
                    tolerance = 64 * np.finfo(np.float64).eps * (qnorm[:, None] + dnorm[None, :] + radius2 + 1)
                    er, ec = np.nonzero(np.abs(distance2 - radius2) <= tolerance)
                    if len(er):
                        delta = q[er] - d[ec]
                        distance2[er, ec] = np.einsum("ij,ij->i", delta, delta)
                ok = ((qid[qs:qe, None] != did[ix][None, :]) &
                      (qgroup[qs:qe, None] != dgroup[ix][None, :]))
                neligible += ok.sum(axis=1)
                near = np.minimum(near, np.where(ok, distance2, np.inf).min(axis=1))
                inside = ok & (distance2 <= radius2)
                counts += inside.sum(axis=1)
                masked = np.where(inside, distance2, np.inf)
                order = np.argsort(masked, axis=1, kind="stable")[:, :min(k, len(ix))]
                chunk_d = np.take_along_axis(masked, order, axis=1)
                chunk_i = ix[order]
                merged_d = np.concatenate((best_d, chunk_d), axis=1)
                merged_i = np.concatenate((best_i, chunk_i), axis=1)
                pick = np.argsort(merged_d, axis=1, kind="stable")[:, :k]
                best_d = np.take_along_axis(merged_d, pick, axis=1)
                best_i = np.take_along_axis(merged_i, pick, axis=1)
        if used_device == "cuda":
            best_d, best_i = tbest_d.cpu().numpy(), tbest_i.cpu().numpy()
            counts, neligible, near = tcounts.cpu().numpy(), teligible.cpu().numpy(), tnear.cpu().numpy()
        best_i[~np.isfinite(best_d)] = -1
        result_idx[qs:qe], result_d2[qs:qe] = best_i, best_d
        eligible_total[qs:qe], in_radius_total[qs:qe] = neligible, counts
        nearest_d2[qs:qe] = near
        if time.monotonic() - last_progress >= 45:
            print(f"[GRID ABM FINAL] support {used_device}: {qe:,}/{nq:,} queries; {nd:,} donors; "
                  f"resident_bank={resident is not None}", flush=True)
            last_progress = time.monotonic()
    return (result_idx, np.sqrt(result_d2), eligible_total, in_radius_total,
            np.sqrt(nearest_d2), dgroup, used_device)


def count_twins(donor_z, query_z, *, donor_ids, query_ids, donor_groups,
                query_groups, radius: float, k_used: int = 128,
                min_count: int = 20, min_group_ess: float = 15.0,
                device: str = "auto", query_batch: int = 128,
                donor_batch: int = 4096, groups_verified: bool = False) -> TwinResult:
    """Count eligible non-self donors outside the query's declared group.

    Arrays contain the *same frozen representation*, in the same feature order;
    this function does not standardize, fit or consult any phenotype. All IDs
    must be unique within their respective input. Family/group identifiers must
    be supplied explicitly. Using IDs as provisional singleton groups is
    permitted only with ``groups_verified=False`` and is labelled accordingly.

    ``d_nearest`` is the nearest eligible donor, even if outside the radius.
    ``d_k_used`` is the farthest retained in-radius neighbour; ``n_used`` tells
    whether the cap was reached. Empty neighbourhood distances are infinity.
    ``ess_group`` sums kernel weights within a family before calculating ESS.
    A supported result with unverified groups is *provisional*, not a claim of
    kinship-verified independent support.
    """
    if not np.isfinite(radius) or radius <= 0 or not np.isfinite(radius * radius):
        raise ValueError("radius must be finite, positive, and safe to square")
    if isinstance(min_count, bool) or not isinstance(min_count, (int, np.integer)) or min_count < 1:
        raise ValueError("min_count must be a positive integer")
    if not np.isfinite(min_group_ess) or min_group_ess <= 0:
        raise ValueError("min_group_ess must be finite and positive")
    dz, qz, di, qi, dg, qg = _validate_inputs(
        donor_z, query_z, donor_ids, query_ids, donor_groups, query_groups,
        k_used, query_batch, donor_batch)
    idx, distances, ne, nt, nearest, group_code, used_device = _scan(
        dz, qz, di, qi, dg, qg, radius=radius, k=k_used, device=device,
        query_batch=query_batch, donor_batch=donor_batch)
    rows = []
    for i in range(len(qz)):
        keep = idx[i] >= 0
        ii, dd = idx[i, keep], distances[i, keep]
        weights = np.exp(-0.5 * (dd / radius) ** 2)
        sumw = float(weights.sum())
        individual_ess = sumw * sumw / float(np.square(weights).sum()) if len(weights) else 0.0
        if len(ii):
            _, inverse = np.unique(group_code[ii], return_inverse=True)
            group_weights = np.bincount(inverse, weights=weights)
            group_ess = sumw * sumw / float(np.square(group_weights).sum())
            n_groups = len(group_weights)
        else:
            group_ess, n_groups = 0.0, 0
        # ESS is a floating-point quantity; avoid failing an exact integer
        # boundary solely because equal weights round 15 to 14.999999999999998.
        enough_ess = group_ess >= min_group_ess - 1e-10 * max(1.0, min_group_ess)
        supported = bool(nt[i] >= min_count and enough_ess)
        if ne[i] == 0:
            reason = "no_eligible_donors"
        elif nt[i] == 0:
            reason = "no_twins_in_radius"
        elif nt[i] < min_count:
            reason = "too_few_twins"
        elif not enough_ess:
            reason = "low_group_ess"
        else:
            reason = "supported_verified_groups" if groups_verified else "supported_provisional_groups"
        rows.append({
            "id": qi[i], "n_eligible": int(ne[i]), "n_twins": int(nt[i]),
            "twins_per_10000": float(10000 * nt[i] / ne[i]) if ne[i] else np.nan,
            "n_used": int(len(ii)), "n_groups_used": int(n_groups),
            "ess_individual": individual_ess, "ess_group": group_ess,
            "d_nearest": float(nearest[i]), "d_k_used": float(dd[-1]) if len(dd) else np.inf,
            "supported": supported, "support_reason": reason,
            "groups_verified": bool(groups_verified), "radius": float(radius),
            "device_used": used_device,
        })
    columns = ["id", "n_eligible", "n_twins", "twins_per_10000", "n_used", "n_groups_used",
               "ess_individual", "ess_group", "d_nearest", "d_k_used", "supported",
               "support_reason", "groups_verified", "radius", "device_used"]
    return TwinResult(pd.DataFrame(rows, columns=columns), idx, distances)


def estimate_radius(build_z, calibration_z, *, donor_ids, query_ids,
                    donor_groups, query_groups, k_ref: int = 64,
                    quantile: float = 0.90, device: str = "auto",
                    query_batch: int = 128, donor_batch: int = 4096) -> float:
    """Freeze one global radius from development-only k-th-neighbour distances.

    The caller must provide a development calibration role with no test people.
    Each calibration query must have at least ``k_ref`` eligible build donors.
    This is one global threshold; it does not force every future query to have
    ``k_ref`` twins. The default 90th percentile is an inherited starting rule,
    not an empirically established optimal radius or an ancestry-neutral rule.
    """
    if not np.isfinite(quantile) or not 0 <= quantile <= 1:
        raise ValueError("quantile must lie between zero and one")
    dz, qz, di, qi, dg, qg = _validate_inputs(
        build_z, calibration_z, donor_ids, query_ids, donor_groups, query_groups,
        k_ref, query_batch, donor_batch)
    if not len(qz):
        raise ValueError("cannot estimate a radius from an empty calibration set")
    _, distances, ne, _, _, _, _ = _scan(
        dz, qz, di, qi, dg, qg, radius=np.inf, k=k_ref, device=device,
        query_batch=query_batch, donor_batch=donor_batch)
    if np.any(ne < k_ref):
        raise ValueError("every radius-calibration query must have at least k_ref eligible donors")
    radius = float(np.quantile(distances[:, k_ref - 1], quantile))
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("calibration produced a zero/non-finite radius; inspect representation or specify a positive frozen radius")
    return radius
