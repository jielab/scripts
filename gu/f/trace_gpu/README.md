# Experimental TRACE CUDA backend

`./gu.sh trace-GPU` retains the installed TRACE CLI, initialization, forward/
backward recursion, ECM fitting and summarization. Two kernels run on CUDA:
Gamma emission log densities over trees/windows, and per-position transition
posteriors used in fitting. They use double precision and disable fused
multiply-add. ARG extraction, SciPy fitting/reductions, NPZ IO and segment
summarization remain CPU. This is not a GPU ARG builder or an upstream release.

Upstream: https://github.com/YulinZhang9806/trace (MIT; see LICENSE).
The installed Python and Cython source/binary hashes are pinned in the build
manifest. Changed source or missing CUDA causes an explicit failure. There is
no silent CPU fallback. A single inference process owns the device.

```bash
bash f/trace_gpu/build.sh
```

`TRACE_GPU_HOME` defaults to `/mnt/d/software/gu/trace-gpu`; `TRACE_GPU_DEVICE`
defaults to 0. The GU environment already has the required NumPy/SciPy/TRACE;
no PyTorch or CuPy dependency is added. `TRACE_SEED` defaults to 42. GPU outputs
use `<result-root>/trace-GPU/<unit>/`, defaulting to
`/mnt/d/analysis/gu/1kg/trace-GPU/`. Command lists and logs also use `trace-GPU`.
Use `GU_RESULT_MODULE=trace-GPU ./gu.sh final --dir-out <result-root>` to write
the corresponding summary under `<result-root>/final/trace-GPU/`.
The upstream package installation is unchanged.

## Shared CPU/GPU numerical correction

The installed upstream initialization calls
`gamma_sp.pdf(subncoal, self.emi2_a1, 1 / self.emi2_b1)`. SciPy interprets the
third positional argument as `loc`, not `scale`. On the real 1KG test this made
all cached null emissions infinite and all fitted posteriors NaN.

Both `gu.sh trace` and `gu.sh trace-GPU` now run through `../trace.py infer`,
which applies exactly the same one-line correction:
`gamma_sp.pdf(subncoal, a=self.emi2_a1, scale=1 / self.emi2_b1)`.
It is recorded as `gu-trace-gamma-scale-v1`. Old run metadata cannot silently
reuse results produced without this correction; use a new output directory or
explicit replacement. CPU/GPU comparisons should use separate output roots.
Inference outputs are checked for finite parameters, normalized finite
posteriors and valid ZIP data before being published.

Other upstream behavior, including forward/backward precision and masked-column
transition indexing, is preserved. The purpose of this backend is comparison
with the same TRACE model; it is not an independent validation of that model.

Historical validation compared kernels, complete fits, decoded posteriors and
classifications at 0.9. The comparisons also checked CLI NPZ outputs
and pipeline segment tables. Speed depends on tree/window
count and CPU extraction cost; a GPU name does not promise a faster full run.

On 2026-10-11, two haplotypes from the full 1KG chr22 THREADS ARG (88,734 trees)
had a maximum CPU/GPU posterior difference of 2.1e-25 and byte-identical final
summary tables. CLI inference times were CPU 3.32/3.17 s and GPU 3.38/3.03 s,
one run per haplotype, including startup and excluding ARG extraction. These
measurements show similar runtimes, not an established speedup. A separate
200 kb region pipeline produced nine identical normalized segment rows.
Historical benchmark scripts and reports were archived during cleanup.
