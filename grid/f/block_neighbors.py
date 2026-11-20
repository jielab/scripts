"""Exact, tiled CPU neighbours in the full frozen consensus-block space.

No N-query by N-donor matrix is allocated. Runtime remains quadratic in sample
count; tiling bounds memory, it does not turn exhaustive search into ANN. Family
counts are accumulated across donor-tile boundaries, independently of top-k.
"""
from __future__ import annotations

import time
import numpy as np
import pandas as pd


class ExactBlockIndex:
    def __init__(self, donor_z, donor_groups, donor_ids=None, donor_batch=4096):
        self.z = donor_z
        if len(donor_z.shape) != 2 or not len(donor_z) or donor_z.shape[1] < 1:
            raise ValueError("Nonempty two-dimensional donor features required")
        self.groups = np.asarray(donor_groups, str)
        self.ids = (np.asarray(donor_ids, str) if donor_ids is not None
                    else np.asarray([f"{i:020d}" for i in range(len(donor_z))]))
        if self.groups.shape != (len(donor_z),) or self.ids.shape != self.groups.shape:
            raise ValueError("Donor identity/group shape mismatch")
        if len(set(self.ids)) != len(self.ids) or any(not s.strip() for s in self.groups):
            raise ValueError("Unique donor IDs and nonempty groups required")
        if int(donor_batch) != donor_batch or donor_batch < 1:
            raise ValueError("Positive donor_batch required")
        self.donor_batch = int(donor_batch)
        levels, self.group_code = np.unique(self.groups, return_inverse=True)
        self.group_lookup = dict(zip(levels, range(len(levels))))
        self.group_sizes = np.bincount(self.group_code)
        # A family can span tiles; stream carry below makes its all-radius ESS exact.
        self.order = np.lexsort((self.ids, self.group_code))
        self.tie_rank = np.argsort(np.argsort(self.ids, kind="stable"), kind="stable")
        self.id_lookup = dict(zip(self.ids, self.groups))

    def estimate(self, n_queries, query_batch=64, k=128):
        d = int(self.z.shape[1]); qb = min(int(query_batch), int(n_queries))
        db = min(self.donor_batch, len(self.z))
        return dict(search="exact_exhaustive_float64_tiled_CPU",
                    donor_n=int(len(self.z)), query_n=int(n_queries), dimensions=d,
                    distance_pairs=int(len(self.z))*int(n_queries),
                    approximate_distance_flops=2*int(len(self.z))*int(n_queries)*d,
                    # Includes several distance/mask/top-k temporaries, not BLAS workspace.
                    bounded_scan_working_bytes=int(8*(db*d+qb*d+6*qb*db+8*qb*k)),
                    full_pair_matrix_allocated=False,
                    complexity="O(Nquery*Ndonor*B) compute; query/donor tiles bound scan memory")

    @staticmethod
    def _squared(q, d, qnorm, dnorm, radius2):
        distance = qnorm[:, None]+dnorm[None, :]-2*(q@d.T)
        if not np.isfinite(distance).all():
            raise ValueError("Nonfinite full-space distances")
        tolerance = 64*np.finfo(float).eps*(qnorm[:, None]+dnorm[None, :]+1)
        near = np.abs(distance) <= tolerance
        if radius2 is not None:
            near |= np.abs(distance-radius2) <= tolerance+64*np.finfo(float).eps*radius2
        rr, cc = np.nonzero(near)
        # Avoid an unbounded pair-by-feature temporary when many profiles coincide.
        for start in range(0, len(rr), 128):
            a, b = rr[start:start+128], cc[start:start+128]
            delta = q[a]-d[b]
            distance[a, b] = np.einsum("ij,ij->i", delta, delta)
        return np.maximum(distance, 0.)

    def _local_top(self, distance, indices, k):
        width = min(k, distance.shape[1])
        selected = np.argpartition(distance, width-1, axis=1)[:, :width]
        # Stable exact ties use the donor ID, never row/tile order.
        for r in range(len(distance)):
            cut = distance[r, selected[r]].max()
            strict = np.flatnonzero(distance[r] < cut)
            tied = np.flatnonzero(distance[r] == cut)
            need = width-len(strict)
            if need:
                tied = tied[np.argsort(self.tie_rank[indices[tied]], kind="stable")[:need]]
            selected[r] = np.r_[strict, tied[:need]]
        return np.take_along_axis(distance, selected, axis=1), indices[selected]

    def validate_queries(self, query_z, query_groups, query_ids=None):
        if len(query_z.shape) != 2 or query_z.shape[1] != self.z.shape[1]:
            raise ValueError("Query feature schema differs from donor features")
        groups = np.asarray(query_groups, str)
        if groups.shape != (len(query_z),) or any(not s.strip() for s in groups):
            raise ValueError("Query groups must be explicit")
        ids = None if query_ids is None else np.asarray(query_ids, str)
        if ids is not None:
            if ids.shape != groups.shape or len(set(ids)) != len(ids):
                raise ValueError("Query IDs must be unique and aligned")
            for identity, group in zip(ids, groups):
                if identity in self.id_lookup and self.id_lookup[identity] != group:
                    raise ValueError("Same ID has a conflicting family identity")
        return groups,ids

    def batches(self, query_z, query_groups, k, *, radius=None, query_ids=None,
                batch_size=64, donor_y=None):
        groups,ids=self.validate_queries(query_z,query_groups,query_ids)
        if int(k) != k or k < 1 or int(batch_size) != batch_size or batch_size < 1:
            raise ValueError("Positive integer k and batch_size required")
        if radius is not None and (not np.isfinite(radius) or radius <= 0):
            raise ValueError("Finite positive radius required")
        values = None if donor_y is None else np.asarray(donor_y, float)
        if values is not None and (values.shape != (len(self.z),) or not np.isfinite(values).all()):
            raise ValueError("Finite donor outcome vector required")
        radius2 = None if radius is None else float(radius)**2
        qcodes = np.asarray([self.group_lookup.get(s, -1) for s in groups], int)
        last = time.monotonic()
        for begin in range(0, len(query_z), int(batch_size)):
            end = min(begin+int(batch_size), len(query_z)); n = end-begin
            q = np.asarray(query_z[begin:end], dtype=np.float64)
            if not np.isfinite(q).all(): raise ValueError("Nonfinite query features")
            qnorm = np.einsum("ij,ij->i", q, q)
            best_d = np.full((n, int(k)), np.inf); best_i = np.full((n, int(k)), -1, int)
            count = np.zeros(n, np.int64); group_count = np.zeros(n, np.int64)
            group_square = np.zeros(n); outcome_sum = np.zeros(n); outcome_square = np.zeros(n)
            distance_sum = np.zeros(n)
            previous_code = -2; previous_count = np.zeros(n, np.int64)
            for start in range(0, len(self.z), self.donor_batch):
                ix = self.order[start:start+self.donor_batch]
                d = np.asarray(self.z[ix], dtype=np.float64)
                if not np.isfinite(d).all(): raise ValueError("Nonfinite donor features")
                dn = np.einsum("ij,ij->i", d, d)
                dd = self._squared(q, d, qnorm, dn, radius2)
                eligible = qcodes[begin:end, None] != self.group_code[ix][None, :]
                distance_sum += np.where(eligible, dd, 0.).sum(axis=1)
                if radius2 is not None:
                    inside = eligible & (dd <= radius2)
                    count += inside.sum(axis=1)
                    gc = self.group_code[ix]
                    offsets = np.r_[0, np.flatnonzero(gc[1:] != gc[:-1])+1]
                    counts = np.add.reduceat(inside.astype(np.int32), offsets, axis=1)
                    group_count += (counts > 0).sum(axis=1)
                    group_square += np.square(counts.astype(float)).sum(axis=1)
                    if gc[0] == previous_code:
                        group_square += 2*previous_count*counts[:, 0]
                        group_count -= (previous_count > 0) & (counts[:, 0] > 0)
                        if len(offsets) == 1:
                            carry = previous_count+counts[:, -1]
                        else:
                            carry = counts[:, -1].copy()
                    else:
                        carry = counts[:, -1].copy()
                    previous_count, previous_code = carry, gc[-1]
                    if values is not None:
                        outcome_sum += inside@values[ix]
                        outcome_square += inside@np.square(values[ix])
                dd[~eligible] = np.inf
                local_d, local_i = self._local_top(dd, ix, int(k))
                candidate_d = np.concatenate([best_d, local_d], axis=1)
                candidate_i = np.concatenate([best_i, local_i], axis=1)
                ranks = np.where(candidate_i >= 0, self.tie_rank[np.maximum(candidate_i, 0)], len(self.z))
                order = np.lexsort((ranks, candidate_d), axis=1)[:, :int(k)]
                best_d = np.take_along_axis(candidate_d, order, axis=1)
                best_i = np.take_along_axis(candidate_i, order, axis=1)
            best_i[~np.isfinite(best_d)] = -1
            nearest = np.sqrt(best_d)
            neligible = np.full(n, len(self.z), np.int64)
            known = qcodes[begin:end] >= 0
            neligible[known] -= self.group_sizes[qcodes[begin:end][known]]
            rms = np.sqrt(np.divide(distance_sum, neligible, out=np.full(n, np.nan), where=neligible > 0))
            ess = np.divide(count.astype(float)**2, group_square, out=np.zeros(n), where=group_square > 0)
            summary = dict(n_support=count, n_support_groups=group_count, support_group_ess=ess,
                           n_eligible=neligible, nearest_distance=nearest[:, 0],
                           rms_eligible_distance=rms,
                           nearest_rms_ratio=np.divide(nearest[:, 0], rms, out=np.zeros(n), where=rms > 0),
                           kth_rms_ratio=np.divide(nearest[:, -1], rms, out=np.zeros(n), where=rms > 0),
                           outcome_sum=outcome_sum, outcome_square=outcome_square)
            yield np.arange(begin, end), best_i, nearest, summary
            if time.monotonic()-last >= 40:
                print(f"[abm] exact full-block search {end:,}/{len(query_z):,} queries; {len(self.z):,} donors", flush=True)
                last = time.monotonic()

    def nearest(self, query_z, query_groups, k, *, query_ids=None, batch_size=64, workers=1):
        from parallel_runtime import worker_count, ordered_map, row_slices
        self.validate_queries(query_z,query_groups,query_ids)
        count = worker_count(workers, max(1, len(query_z)//batch_size),
                             memory_per_worker_gib=.5+self.estimate(batch_size,batch_size,k)["bounded_scan_working_bytes"]/2**30,
                             reserve_gib=self.z.nbytes/2**30)
        if count > 1:
            slices = row_slices(len(query_z), count)
            arguments = ((self, query_z[a:b], np.asarray(query_groups)[a:b], k,
                          None if query_ids is None else np.asarray(query_ids)[a:b], batch_size)
                         for a,b in slices)
            parts = list(ordered_map(_nearest_part, arguments, count))
            return np.concatenate([p[0] for p in parts]), np.concatenate([p[1] for p in parts])
        indices = np.full((len(query_z), int(k)), -1, int)
        distances = np.full((len(query_z), int(k)), np.inf)
        for rows, ix, dd, _ in self.batches(query_z, query_groups, k, query_ids=query_ids, batch_size=batch_size):
            indices[rows], distances[rows] = ix, dd
        return indices, distances


def _nearest_part(bank, queries, groups, k, ids, batch_size):
    return bank.nearest(queries, groups, k, query_ids=ids, batch_size=batch_size, workers=1)
