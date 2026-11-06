"""PRS-CSx block algebra; float64 CUDA via GRID's existing PyTorch runtime.

NumPy retains all RNG ownership and the GIG sampler remains unchanged on CPU.
Small blocks can use SciPy to avoid per-block CUDA transfers/synchronization.
LD blocks are cached within a bounded device budget. No CUDA failure is hidden.
"""
from __future__ import annotations
import argparse
from collections import OrderedDict
import json
import os
import time
import numpy as np


def gpu_min_block():
    """Zero preserves existing callers; new launchers explicitly select hybrid."""
    value = int(os.getenv("GRID_CSX_GPU_MIN_BLOCK", "0"))
    if value < 0:
        raise ValueError("GRID_CSX_GPU_MIN_BLOCK must be a nonnegative integer")
    return value


def linalg_identity(backend):
    # Routing changes execution, not the float64 block conditional or RNG order.
    # Keep the fitted-input identity so completed chromosomes remain reusable.
    return "grid-block-f64-v1"


class BlockLinalg:
    def __init__(self, backend="cpu", device=0):
        if backend not in {"auto", "cpu", "cuda"}:
            raise ValueError("backend must be auto|cpu|cuda")
        self.backend = "cpu"
        self.device = int(device)
        self.cache = OrderedDict()
        self.cache_bytes = 0
        self.calls = 0
        self.cpu_calls = 0
        self.cuda_calls = 0
        self.min_gpu_block = gpu_min_block()
        self.seconds = 0.0
        if backend == "cpu":
            return
        try:
            import torch
        except ImportError:
            if backend == "cuda":
                raise RuntimeError("CUDA requested: install GRID's torch environment")
            return
        if not torch.cuda.is_available():
            if backend == "cuda":
                raise RuntimeError("CUDA requested but unavailable")
            return
        if not 0 <= self.device < torch.cuda.device_count():
            raise ValueError("Invalid GRID_CSX_GPU device index")
        os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
        self.torch = torch
        self.dev = torch.device("cuda", self.device)
        torch.cuda.set_device(self.dev)
        torch.use_deterministic_algorithms(True)
        torch.linalg.cholesky(torch.eye(2, device=self.dev, dtype=torch.float64))
        torch.cuda.synchronize(self.dev)
        torch.cuda.reset_peak_memory_stats(self.dev)
        free, _ = torch.cuda.mem_get_info(self.dev)
        self.cache_limit = int(min(float(os.getenv("GRID_CSX_GPU_CACHE_GIB", "4"))*2**30, free*.4))
        if self.cache_limit < 0:
            raise ValueError("GRID_CSX_GPU_CACHE_GIB must be nonnegative")
        self.backend = "cuda"

    @property
    def is_cuda(self):
        return self.backend == "cuda"

    def preload_blocks(self, blocks):
        # Upload lazily once per resident block; do not exhaust VRAM by eagerly
        # copying every population's chromosome into device memory.
        return blocks

    def _resident(self, ld):
        key = id(ld)
        if key in self.cache:
            self.cache.move_to_end(key)
            return self.cache[key][1]
        size = ld.size * 8
        while self.cache and self.cache_bytes + size > self.cache_limit:
            _, (host, _) = self.cache.popitem(last=False)
            self.cache_bytes -= host.size * 8
        value = self.torch.as_tensor(ld, device=self.dev, dtype=self.torch.float64)
        if size <= self.cache_limit:
            self.cache[key] = (ld, value)  # Retain host identity against id reuse.
            self.cache_bytes += size
        return value

    def beta_block(self, ld, psi, beta_mrg, noise, sigma_over_n):
        started = time.perf_counter()
        if not self.is_cuda or ld.shape[0] < self.min_gpu_block:
            from scipy import linalg
            dinvt = ld + np.diag(1.0 / np.asarray(psi).ravel())
            chol = linalg.cholesky(dinvt)
            tmp = linalg.solve_triangular(chol, beta_mrg, trans="T") + np.sqrt(sigma_over_n)*noise
            beta = linalg.solve_triangular(chol, tmp, trans="N")
            quad = float((beta.T @ dinvt @ beta).item())
            self.cpu_calls += 1
        else:
            torch = self.torch
            with torch.no_grad():
                a = self._resident(ld).clone()
                pg = torch.as_tensor(psi, device=self.dev, dtype=torch.float64).flatten()
                a.diagonal().add_(1.0 / pg)
                b = torch.as_tensor(beta_mrg, device=self.dev, dtype=torch.float64)
                e = torch.as_tensor(noise, device=self.dev, dtype=torch.float64)
                # Upper R, A = R.T @ R, exactly the SciPy reference convention.
                chol = torch.linalg.cholesky(a, upper=True)
                tmp = torch.linalg.solve_triangular(chol.T, b, upper=False)
                tmp.add_(e * float(np.sqrt(sigma_over_n).item()))
                value = torch.linalg.solve_triangular(chol, tmp, upper=True)
                quad = float((value.T @ a @ value).item())
                beta = value.cpu().numpy()
            self.cuda_calls += 1
        if not np.isfinite(quad) or not np.isfinite(beta).all():
            raise FloatingPointError("Nonfinite PRS-CSx block update; refusing posterior output")
        self.calls += 1
        self.seconds += time.perf_counter() - started
        return beta, quad

    def synchronize(self):
        if self.is_cuda:
            self.torch.cuda.synchronize(self.dev)

    def info(self):
        result = dict(backend=self.backend, numpy=np.__version__, gpu_name=None, device=self.device if self.is_cuda else None,
                      block_calls=self.calls, cpu_block_calls=self.cpu_calls, cuda_block_calls=self.cuda_calls,
                      gpu_min_block=self.min_gpu_block, block_seconds=self.seconds, gpu_peak_memory_mb=0)
        if self.is_cuda:
            result.update(gpu_name=self.torch.cuda.get_device_name(self.dev), torch=self.torch.__version__,
                          cuda=self.torch.version.cuda, ld_cache_bytes=self.cache_bytes,
                          gpu_peak_memory_mb=self.torch.cuda.max_memory_allocated(self.dev)/2**20)
        return result


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--backend", choices=["auto", "cpu", "cuda"], default="auto")
    p.add_argument("--gpu-device", type=int, default=0)
    a = p.parse_args()
    print(json.dumps(BlockLinalg(a.backend, a.gpu_device).info()))
