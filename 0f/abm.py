"""Shared ABM kernels for LE8 and GRID (no cohort or outcome policy).

Exact bounded-memory Euclidean retrieval, explicit ID/family/fold exclusion,
stable hash ties, and neural attention primitives. Callers choose the tie
tolerance: LE8 uses 1e-10; GRID uses exact distances. No approximate search,
mixed precision, or implicit CPU fallback is used.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path
import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

SOURCE_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


class BlockedKNN:
    def __init__(self, coordinates, ids, groups, tie_keys, device='cpu',
                 dtype='float64', query_block=256, donor_block=8192,
                 tie_tolerance=1e-10):
        x = np.asarray(coordinates)
        self.ids, self.groups = np.asarray(ids, str), np.asarray(groups, str)
        if x.ndim != 2 or not len(x) or not x.shape[1] or not np.isfinite(x).all():
            raise ValueError('Finite nonempty donor matrix required')
        if len(self.ids) != len(x) or len(set(self.ids)) != len(x) or len(self.groups) != len(x) or len(tie_keys) != len(x):
            raise ValueError('Invalid donor identities')
        if query_block < 1 or donor_block < 1 or not np.isfinite(tie_tolerance) or tie_tolerance < 0:
            raise ValueError('Positive blocks and nonnegative tie tolerance required')
        self.device = torch.device(device)
        if self.device.type not in {'cpu', 'cuda'}:
            raise ValueError('ABM device must be cpu or cuda[:index]')
        if self.device.type == 'cuda' and not torch.cuda.is_available():
            raise RuntimeError('CUDA explicitly requested, no CPU fallback')
        self.dtype = {'float64': torch.float64, 'float32': torch.float32}[dtype]
        self.order = np.argsort(np.asarray(tie_keys), kind='stable')
        self.x = torch.as_tensor(x[self.order].copy(), dtype=self.dtype, device=self.device)
        self.query_block, self.donor_block = int(query_block), int(donor_block)
        self.tie_tolerance = float(tie_tolerance)
        # Encode once, instead of allocating/comparing Python strings for every
        # query-by-donor block and copying those masks across PCIe.
        self.id_lookup = {s: i for i, s in enumerate(self.ids)}
        levels, codes = np.unique(self.groups, return_inverse=True)
        self.group_lookup = {s: i for i, s in enumerate(levels)}
        self.bid = torch.as_tensor(self.order, device=self.device)
        self.bg = torch.as_tensor(codes[self.order], device=self.device)

    def _rank(self, distance):
        return torch.round(distance / self.tie_tolerance) if self.tie_tolerance else distance

    @torch.no_grad()
    def query(self, coordinates, ids, groups, k, query_folds=None, bank_folds=None,
              *, radius=None, return_counts=False):
        q = np.ascontiguousarray(coordinates)
        ids, groups = np.asarray(ids, str), np.asarray(groups, str)
        if q.ndim != 2 or q.shape[1] != self.x.shape[1] or len(ids) != len(q) or len(groups) != len(q) or not np.isfinite(q).all():
            raise ValueError('Invalid query')
        if int(k) != k or k < 0 or (k == 0 and not return_counts):
            raise ValueError('Positive k required unless counting radius matches')
        if (query_folds is None) != (bank_folds is None):
            raise ValueError('Supply both query and bank folds')
        if radius is not None and (np.isnan(radius) or radius < 0):
            raise ValueError('Nonnegative radius required')
        bf = qf = None
        if bank_folds is not None:
            if len(bank_folds) != len(self.x) or len(query_folds) != len(q):
                raise ValueError('Fold shape mismatch')
            _, codes = np.unique(np.r_[bank_folds, query_folds], return_inverse=True)
            bf = torch.as_tensor(codes[:len(self.x)][self.order], device=self.device)
            qf = torch.as_tensor(codes[len(self.x):], device=self.device)
        qi = torch.as_tensor([self.id_lookup.get(s, -1) for s in ids], dtype=torch.long, device=self.device)
        qg = torch.as_tensor([self.group_lookup.get(s, -1) for s in groups], dtype=torch.long, device=self.device)
        keep = min(int(k), len(self.x))
        indices = np.full((len(q), keep), -1, dtype=np.int64)
        distances = np.full((len(q), keep), np.inf)
        counts = np.zeros(len(q), dtype=np.int64)
        for begin in range(0, len(q), self.query_block):
            end = min(begin + self.query_block, len(q))
            xx = torch.as_tensor(q[begin:end], dtype=self.dtype, device=self.device)
            best_d = xx.new_empty((len(xx), 0))
            best_i = torch.empty((len(xx), 0), dtype=torch.long, device=self.device)
            count = torch.zeros(len(xx), dtype=torch.long, device=self.device)
            for start in range(0, len(self.x), self.donor_block):
                stop = min(start + self.donor_block, len(self.x))
                d = torch.cdist(xx, self.x[start:stop], compute_mode='donot_use_mm_for_euclid_dist')
                allowed = (qi[begin:end, None] != self.bid[None, start:stop]) & (qg[begin:end, None] != self.bg[None, start:stop])
                if bf is not None:
                    allowed &= qf[begin:end, None] != bf[None, start:stop]
                d.masked_fill_(~allowed, torch.inf)
                if radius is not None:
                    count += (allowed & (d <= radius)).sum(1)
                if not keep:
                    continue
                rank = self._rank(d)
                width = min(keep, stop - start)
                # Select local k first. Resolve boundary ties by donor hash rank;
                # sort only <=2k candidates, rather than the entire donor block.
                cutoff = torch.topk(rank, width, dim=1, largest=False).values[:, -1:]
                strict = (rank < cutoff).sum(1, keepdim=True)
                tied = rank == cutoff
                eligible = (rank < cutoff) | (tied & (tied.cumsum(1) <= width - strict))
                take = torch.topk(rank.masked_fill(~eligible, torch.inf), width,
                                  dim=1, largest=False, sorted=False).indices
                all_d = torch.cat([best_d, d.gather(1, take)], dim=1)
                all_i = torch.cat([best_i, take + start], dim=1)
                by_id = torch.argsort(all_i, dim=1, stable=True)
                all_d, all_i = all_d.gather(1, by_id), all_i.gather(1, by_id)
                chosen = torch.argsort(self._rank(all_d), dim=1, stable=True)[:, :keep]
                best_d, best_i = all_d.gather(1, chosen), all_i.gather(1, chosen)
            if keep:
                dd = best_d.cpu().numpy()
                jj = self.order[best_i.cpu().numpy()]
                distances[begin:end] = dd
                indices[begin:end] = np.where(np.isfinite(dd), jj, -1)
            if radius is not None:
                counts[begin:end] = count.cpu().numpy()
        return (indices, distances, counts) if return_counts else (indices, distances)


def neighbors(bank, query, bank_ids, query_ids, bank_groups, query_groups, tie, k,
              device, query_block=256, donor_block=8192, bank_folds=None,
              query_folds=None, radius=None):
    """GRID-compatible exact-distance API; LE8 uses BlockedKNN's tolerance."""
    index = BlockedKNN(bank, bank_ids, bank_groups, tie, device=device,
                       query_block=query_block, donor_block=donor_block, tie_tolerance=0.)
    return index.query(query, query_ids, query_groups, k, query_folds, bank_folds,
                       radius=radius, return_counts=True)


def safe_softmax(logits, allowed, dim):
    masked = logits.masked_fill(~allowed, -torch.inf)
    maximum = masked.amax(dim=dim, keepdim=True)
    maximum = torch.where(torch.isfinite(maximum), maximum, torch.zeros_like(maximum))
    exp = torch.exp(masked - maximum)
    return exp / exp.sum(dim=dim, keepdim=True).clamp_min(1e-30)


def family_effective_size(weights, codes):
    """Exact family ESS with O(batch*k) space, including empty neighborhoods."""
    order = torch.argsort(codes, dim=1, stable=True)
    sorted_codes, sorted_w = codes.gather(1, order), weights.gather(1, order)
    first = torch.ones_like(sorted_codes, dtype=torch.bool)
    first[:, 1:] = sorted_codes[:, 1:] != sorted_codes[:, :-1]
    compact = first.long().cumsum(1) - 1
    mass = torch.zeros_like(weights).scatter_add(1, compact, sorted_w)
    ss = mass.square().sum(1)
    return torch.where(ss > 0, 1 / ss.clamp_min(1e-30), torch.zeros_like(ss))


def masked_mse(decoded, target, mask):
    # Avoid CUDA boolean indexing/nonzero and the mask.any() host synchronization.
    error = torch.where(mask, decoded - target, torch.zeros_like(decoded))
    return error.square().sum() / mask.sum().clamp_min(1)


class GradientTracker:
    """Track ever-nonzero gradients, retiring proven parameters once per epoch.

    Avoid launching a reduction for every already-proven parameter at every
    optimizer step. Flush only at epoch boundaries, not inside the hot loop.
    """
    def __init__(self, named_parameters, device):
        self.pending = list(named_parameters)
        self.flags = torch.zeros(len(self.pending), dtype=torch.bool, device=device)
        self.seen = set()

    def update(self):
        if self.pending:
            self.flags.logical_or_(torch.stack([
                torch.any(p.grad != 0) if p.grad is not None else self.flags.new_tensor(False)
                for _, p in self.pending]))

    def flush(self):
        flags = self.flags.cpu().tolist()
        self.seen.update(name for (name, _), flag in zip(self.pending, flags) if flag)
        self.pending = [pair for pair, flag in zip(self.pending, flags) if not flag]
        self.flags = self.flags.new_zeros(len(self.pending))
        return self.seen


class ContextBlock(nn.Module):
    def __init__(self, width, heads, dropout, mode="learned"):
        super().__init__()
        self.heads, self.width, self.mode = heads, width, mode
        self.norm1, self.norm2 = nn.LayerNorm(width), nn.LayerNorm(width)
        self.qkv = nn.Linear(width, 3 * width)
        self.out = nn.Linear(width, width)
        self.ffn = nn.Sequential(
            nn.Linear(width, 4 * width),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(4 * width, width),
        )
        self.dropout = nn.Dropout(dropout)
        self.attention_dropout = dropout
        self.intervention = (
            None  # evaluation-only intervention; never a fitted parameter
        )

    def forward(self, x, capture=False):
        b, n, _ = x.shape
        q, k, v = (
            self.qkv(self.norm1(x))
            .view(b, n, 3, self.heads, self.width // self.heads)
            .permute(2, 0, 3, 1, 4)
        )
        mode = self.intervention or self.mode
        weights = None
        if mode == "learned" and not capture:
            context = F.scaled_dot_product_attention(
                q, k, v, dropout_p=self.attention_dropout if self.training else 0.0
            )
        else:
            if mode == "uniform":
                weights = x.new_full((b, self.heads, n, n), 1 / n)
            elif mode == "identity":
                weights = torch.eye(n, device=x.device, dtype=x.dtype)[
                    None, None
                ].expand(b, self.heads, -1, -1)
            else:
                weights = torch.softmax(
                    q @ k.transpose(-2, -1) / math.sqrt(q.shape[-1]), dim=-1
                )
            context = (
                F.dropout(weights, p=self.attention_dropout, training=self.training) @ v
            )
        x = x + self.dropout(
            self.out(context.transpose(1, 2).reshape(b, n, self.width))
        )
        x = x + self.dropout(self.ffn(self.norm2(x)))
        return x, weights

