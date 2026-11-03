"""PHOLE CUDA kernels. CPU owns random draws; explicit CUDA never falls back."""
from __future__ import annotations
import os
import time
import numpy as np

TIMINGS = []


def resolve_device(request="auto"):
    if request not in {"auto", "cpu", "cuda"}:
        raise ValueError("device must be auto|cpu|cuda")
    # User policy: prefer available CUDA; retain explicit CPU for measured slow workloads.
    if request == "cpu":
        return "cpu"
    try:
        import torch
        available = torch.cuda.is_available()
    except ImportError:
        available = False
    if request == "cuda" and not available:
        raise RuntimeError("PHOLE CUDA requested but unavailable in this Python environment")
    return "cuda" if available else "cpu"


def preflight(request):
    device = resolve_device(request)
    info = dict(backend=device, numpy=np.__version__, gpu_name=None, gpu_peak_memory_mb=0)
    if device == "cuda":
        os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
        import torch
        torch.use_deterministic_algorithms(True)
        torch.set_float32_matmul_precision("highest")
        # Exercise the actual linalg runtime before reading/writing analysis data.
        x = torch.eye(4, device="cuda", dtype=torch.float32)
        torch.linalg.svd(x)
        torch.cuda.synchronize()
        torch.cuda.reset_peak_memory_stats()
        info.update(gpu_name=torch.cuda.get_device_name(), torch=torch.__version__, cuda=torch.version.cuda)
    return info


def memory_mb(device):
    if device != "cuda":
        return 0
    import torch
    torch.cuda.synchronize()
    return torch.cuda.max_memory_allocated() / 2**20


def randomized_svd_gpu(matrix, rank, seed=1, n_iter=2, oversample=8, device="cuda", omega=None):
    import torch
    if resolve_device(device) != "cuda":
        raise RuntimeError("randomized_svd_gpu requires CUDA")
    started = time.perf_counter()
    a = torch.as_tensor(np.asarray(matrix, dtype=np.float32, order="C"), device="cuda")
    m, n = a.shape
    k = max(1, min(int(rank), m, n))
    ell = min(k + int(oversample), m, n)
    # Same projection as the NumPy reference, independent of global RNG state.
    if omega is None:
        omega = np.random.default_rng(seed).normal(size=(n, ell))
    if np.shape(omega) != (n, ell):
        raise ValueError("Invalid SVD projection shape")
    with torch.no_grad():
        projection = torch.as_tensor(omega, device="cuda", dtype=a.dtype)
        q = torch.linalg.qr(a @ projection, mode="reduced").Q
        for _ in range(max(0, int(n_iter))):
            z = torch.linalg.qr(a.T @ q, mode="reduced").Q
            q = torch.linalg.qr(a @ z, mode="reduced").Q
        ub, s, vh = torch.linalg.svd(q.T @ a, full_matrices=False)
        u = q @ ub[:, :k]
        result = u.cpu().numpy(), s[:k].cpu().numpy(), vh[:k].cpu().numpy()
    TIMINGS.append(dict(operation="randomized_svd", backend="cuda", n=m, p=n, rank=k, seed=seed,
                        wall_seconds=time.perf_counter()-started, gpu_peak_memory_mb=memory_mb("cuda")))
    return result


class PairPredictor:
    """Resident float64 model; transfer only each position kernel and predictions.

    Accumulate all factor SSEs in one pass through the sampled pairs, avoiding
    rank separate gathers/copies. Table construction and statistics stay on CPU.
    """
    def __init__(self, row_base, alpha, beta, col_base, batch=250_000):
        import torch
        self.torch = torch
        self.batch = int(batch)
        if self.batch < 1:
            raise ValueError("prediction batch must be positive")
        self.row = torch.as_tensor(row_base, device="cuda", dtype=torch.float64)
        self.alpha = torch.as_tensor(alpha, device="cuda", dtype=torch.float64)
        self.beta = torch.as_tensor(beta, device="cuda", dtype=torch.float64)
        self.col = torch.as_tensor(col_base, device="cuda", dtype=torch.float64)
        self.seconds = 0.0
        self.calls = 0

    def evaluate(self, kernel, rows, cols, observed):
        torch = self.torch
        started = time.perf_counter()
        rank = self.alpha.shape[1]
        out = np.empty(len(rows), dtype=np.float64)
        with torch.no_grad():
            kg = torch.as_tensor(kernel, device="cuda", dtype=torch.float64)
            beta_query = self.beta @ kg.T
            col_query = kg @ self.col
            ss = torch.zeros(rank + 1, device="cuda", dtype=torch.float64)
            for start in range(0, len(rows), self.batch):
                stop = min(len(rows), start + self.batch)
                r = torch.as_tensor(rows[start:stop], device="cuda", dtype=torch.long)
                c = torch.as_tensor(cols[start:stop], device="cuda", dtype=torch.long)
                y = torch.as_tensor(observed[start:stop], device="cuda", dtype=torch.float64)
                baseline = self.row[r] + col_query[c]
                terms = self.alpha[r] * beta_query[:, c].T
                cumulative = baseline[:, None] + torch.cumsum(terms, dim=1)
                ss[0] += ((y-baseline)**2).sum()
                ss[1:] += ((y[:, None]-cumulative)**2).sum(dim=0)
                # Match the reference's sum-then-add for final predictions.
                out[start:stop] = (baseline + terms.sum(dim=1)).cpu().numpy()
            sse = ss.cpu().numpy()
        self.seconds += time.perf_counter() - started
        self.calls += 1
        return out, sse
