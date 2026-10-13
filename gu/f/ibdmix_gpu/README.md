# Experimental IBDmix CUDA caller

`/mnt/d/scripts/gu/gu.sh ibdmix-GPU` uses a local CUDA implementation of the
IBDmix segment stack. This is an experimental GU backend, not an upstream GPU
release. It was derived from PrincetonUniversity/IBDmix commit
`6a9d4bc67dcaf602e65cac0fa6191e1e7a1291ef` and is GPL-3.0 licensed (see LICENSE).
Upstream: https://github.com/PrincetonUniversity/IBDmix

The calling unit is a modern **individual's genotype**, compared with one archaic
individual. Autosomal phased haplotypes are not called separately. The selected
modern cohort supplies allele frequencies. Male non-PAR X is converted using the
existing GU X contract and independent non-PAR components are called separately.

CPU and GPU entry points default to exactly
`Altai,Chagyr,Denisova,Denisova25,Vindija`. New calls reject `.2013` references;
historical replication files remain readable.

## Run

```bash
cd /mnt/d/scripts/gu
bash f/ibdmix_gpu/build.sh

./gu.sh ibdmix-GPU --grch 37 --target ukb \
  --dir-target /mnt/f/gen/ukb/37/imp/chr \
  --dir-out /mnt/d/analysis/gu/ukb --postfix 003 \
  --keep-psam /mnt/f/gen/ukb/37/imp/chr1.psam,2000,3 --keep-males \
  --archaic-gen Altai,Chagyr,Denisova,Denisova25,Vindija
```

This example writes to `/mnt/d/analysis/gu/ukb003/ibdmix-GPU/<unit>/`.
The default 1KG destination is `/mnt/d/analysis/gu/1kg/ibdmix-GPU/<unit>/`.
`--dir-out` selects the shared result root; each module uses its exact CLI name
for the next directory, command lists and logs. GPU runs default to `--jobs 4`
chromosome workers, subject to the RAM budget. At the default `--memory-cap 32G`,
each UKB worker receives a 4096 MiB PLINK workspace, with 2048 MiB per worker
and 8192 MiB shared overhead reserved. An explicit `GU_PREP_MEMORY_MB` is honored
and can reduce concurrency. `--jobs 1` restores serial chromosome processing.
Preparation can overlap across chromosomes; CUDA calls share a per-user,
per-device lock and run one at a time to avoid adding their peak VRAM usage.
The lock coordinates GU IBDmix callers, not unrelated GPU programs. Each
chromosome still uses one population worker. Results retain method `ibdmix`; their
provenance records `cuda-segment-v1`, a distinct pipeline version and the
binary SHA-256. For a summary using these GPU calls, run
`GU_RESULT_MODULE=ibdmix-GPU ./gu.sh final --dir-out /mnt/d/analysis/gu/ukb003`.
It writes to `ukb003/final/ibdmix-GPU/`; use the same environment variable with
`./gu.sh shiny` to view it. CPU runs use `ukb003/ibdmix/` and `ukb003/final/`.

## What is accelerated

The GPU maintains a separate streaming segment stack for every selected modern
individual, accumulates LOD scores in double precision and emits the original
region statistics. The CPU computes population frequencies and the three
possible genotype LODs per site using the original upstream formulas. It also
handles parsing, masks, VCF conversion, compression, QC and final filtering.
Do not interpret caller speedups as whole-pipeline speedups.

Genotype preparation now uses the [shared HTSlib reader](../ibdmix_prepare/README.md):
one modern genotype pass supplies all five references, with QC during emission.
CPU/GPU callers receive the same canonical genotype text as before. Progress on
screen is limited to preparation and reference/chromosome completion; detailed
diagnostics are retained in log files.

UKB preparation also reports the current process, elapsed time and input file
offset every 30 seconds, including an estimate for the current read when it is
advancing. A file offset is not whole-chromosome completion: VCF export, audit,
shared genotype preparation and reference calls follow the PGEN subset pass.
The monitor is separate from the preparation implementation so adding progress
does not invalidate its scientific cache receipts. Already-running schedulers
keep their original concurrency; these defaults apply on the next launch.

Site batches preserve all pending stacks, so a batch boundary does not create
a segment boundary or change frequencies. Candidate storage grows before each
batch by appending pages of up to 2048 nodes per individual. Existing pages
remain allocated in place; growth neither copies the complete stack nor doubles
its capacity. Only the pages required for the next batch are added. This avoids
the old simultaneous allocation of a full old and new stack, which could fail
even when the live stack would fit in VRAM.
Growth is bounded by available VRAM and `IBDMIX_GPU_MAX_NODES`. Exceeding
the limit produces an explicit failure, never a truncated successful result
or a silent CPU fallback. Output order matches upstream's site/sample order.

## Configuration and limitations

The local build uses CUDA 12.8 from the existing `as3_mamba` environment and
`-arch=native`; rebuild on a different GPU. The caller itself is a native binary
and does not require PyTorch. Build overrides: `IBDMIX_NVCC`, `IBDMIX_RUNTIME`
(upstream source checkout), `IBDMIX_GPU_BIN` (output executable).

Runtime overrides, also recorded in generated chromosome commands:

| Environment variable | Default | Meaning |
| --- | --- | --- |
| `IBDMIX_GPU_DEVICE` | 0 | CUDA device ordinal |
| `IBDMIX_GPU_BATCH_SITES` | 2048 | Sites transferred per batch |
| `IBDMIX_GPU_NODES` | 2048 | Initial candidate nodes per individual |
| `IBDMIX_GPU_MAX_NODES` | 131072 | Upper bound for candidate nodes per individual |
| `IBDMIX_GPU_EVENTS` | 1048576 | Maximum emitted segments per batch |

The caller accepts one chromosome, canonical `generate_gt` input and the flags
used by GU, including sample selection, masks, error rates, MAC, LOD, inclusive
ends and `--more-stats`. It does not implement upstream `--write-snps` or
`--write-lods`; unsupported options fail explicitly. It rejects malformed
records, duplicate selected sample IDs, selection of the archaic as a modern
sample and unsorted positions.

VRAM consumption grows with cohort size and the longest retained candidate
stack. A test with 917 individuals is not validation for all approximately
500,000 UKB participants at once. Do not split a cohort into independent
frequency-estimation batches solely to fit GPU memory: that changes the model.
Use the original CPU entry point if a fixed cohort exceeds available memory.
`--memory-cap` controls system RAM, not GPU VRAM. Other GPU workloads also reduce
the VRAM available to this caller. Caller failures now include the reference,
population, exit status, detailed log path and the final error lines in the
main run log.

Rebuilding changes the GPU binary SHA-256 recorded in result provenance.
Existing results remain readable, but cannot be resumed with a different binary.
Use a new output postfix to retain the existing run when rerunning a full genome
with a rebuilt caller. `--replace-ibdmix TRUE` explicitly discards and recomputes
the selected existing results.

## Validation

Historical CPU/GPU checks compared complete outputs byte for byte, covering
archaic genotype codes, missingness, monomorphic sites, MAC thresholds, sample
subsets, masks, segment boundaries, storage growth and retained-suffix replay.
Small VCF fixtures exercised all five references through the existing pipeline.
Standalone test/benchmark scripts and historical logs are kept out of the source
tree; runtime backend and output checks remain in the production code.

## Local measurements, 2026-10-11

On the RTX 5090 Laptop GPU, 917 UKB individuals and the prepared Altai chr22
genotypes (746,710 sites) produced 68,502 raw segments, byte-identical to the
original caller in all three runs. CPU wall times were 28.21, 43.69 and 31.12 s;
GPU wall times were 13.56, 13.80 and 14.53 s. The median caller speedup was 2.26x.
These real-data timing runs used no exclusion mask; masks were exercised in
the synthetic equivalence and five-reference pipeline tests.

A male chrX region at X:10,000,000–12,000,000 also produced byte-identical output
(2,392 segments). That small input took 0.72 s on CPU and 0.81 s on GPU in one
run; GPU startup and transfer costs can outweigh benefits for small inputs.
This is not a whole-UKB or whole-chrX validation.

The paged-memory correction was also verified against the failing full chr3
Denisova call with the same 4,564 individuals and excluded-site mask. It completed
all 4,233,154 sites in 235.808 seconds of caller time, producing 1,071,267 raw
segments. The maximum stack held 83,968 nodes per individual (11.421 GiB of node
storage); the previous 65,536-to-131,072 reallocation would require 26.742 GiB
of simultaneous old and new node storage. This real-data check verifies completion;
CPU equivalence was tested separately in 28 synthetic comparisons and the
five-reference pipeline. Historical reports/logs were removed from the source tree.
