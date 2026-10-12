# Experimental PhyML CUDA backend

`./gu.sh phyml-GPU` runs the original PhyML tree search and model optimization
with a CUDA nucleotide pruning kernel. This is a local experimental build, not
an official GPU release. Source: https://github.com/stephaneguindon/phyml at
`42bb4b52f479938a3427dfe4357b861f0e306686` (3.3.20260528, GPL-3.0).

The adapter changes the partial-likelihood update in `src/lk.c`; the CUDA kernel
evaluates sites and rate categories using double precision without fused
multiply-add. Tip ambiguities, all-missing subtrees, zero-weight bootstrap
patterns and the original powers-of-two scaling are preserved. Branch/model
optimization, topology search, likelihood reductions and bootstrap sampling
remain the upstream CPU implementation. GU's VCF/LD/haplotype preparation also
remains CPU. The GU model is still HKY85+G4+I with 100 bootstrap replicates.

Build using the existing CUDA 12.8 installation:

```bash
bash f/phyml_gpu/build.sh
```

The original installed `phyml` and `phyml-mpi` binaries are not replaced.
The build is installed under `/mnt/d/software/gu/phyml-gpu`, with a source/binary
SHA-256 manifest. A dedicated source checkout is used. `PHYML_GPU_BIN` and
`PHYML_GPU_DEVICE` select binary/device. `PHYML_TREE_SEED` is propagated to both
CPU and GPU workers; the recorded comparisons used 42 and one worker.
GPU outputs use `<result-root>/phyml-GPU/<unit>/`, defaulting to
`/mnt/d/analysis/gu/1kg/phyml-GPU/`. Command lists and logs also use `phyml-GPU`.
Use `GU_RESULT_MODULE=phyml-GPU ./gu.sh final --dir-out <result-root>` to write
the corresponding summary under `<result-root>/final/phyml-GPU/`.
No CPU fallback is used for
missing CUDA, allocation errors, or unsupported models/scaling.

**This first implementation is slower on the tested small trees.** Host/device
transfers occur for each partial update. An 8-tip, 120-site synthetic alignment
took 0.09 s on CPU and 4.74 s on GPU; a 62-tip, 450-site 1KG-derived alignment
took 2.87 s and 131.92 s respectively (one run each, no bootstrap). These are
observations, not expected UKB speedups. GPU utilization alone is not evidence
of useful acceleration. Keeping likelihood buffers resident on the device and
batching tree operations would be a separate optimization effort.

Historical kernel validation covered 48 tip/internal, ambiguity, category, scaling and batch
size combinations. Historical tree comparisons also checked log likelihood,
unrooted splits, branch lengths, bootstrap counts and supports. Floating-point
changes can affect optimization/search, so byte equality is not a universal
guarantee. Historical benchmark scripts and reports were archived during cleanup.
