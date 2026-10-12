# Shared IBDmix genotype preparation

`f/ibdmix.sh` automatically builds and uses `prepare.cc` for CPU and CUDA runs.
HTSlib reads the modern VCF once per analysis unit; decoded genotypes are shared
by all missing archaic references. It reads VCF/BCF directly, avoiding the previous
`bcftools → text VCF → generate_gt` pipeline for every reference. All compressors
share four threads by default (`IBDMIX_PREPARE_THREADS`).

The merge follows upstream `generate_gt`: modern indels remain merge barriers,
archaic homozygous-reference records remain informative at shared sites, and
archaic-only non-reference sites are emitted only before modern EOF. Missing
alleles, haploid duplication, allele matching, ordering, sample identity and
chrX non-PAR filtering retain the existing contract. QC is accumulated while
writing, avoiding a second full genotype decompression. Files become reusable
only after all readers and compressors succeed.

Genotype outputs remain gzip-compatible text. Up to five compressed inputs are
kept under `/tmp` and removed after their calls; no large uncompressed modern VCF
is cached. Completed reference calls can still be resumed independently. The
pipeline version and preparation binary SHA-256 participate in provenance.

Build explicitly with `bash f/ibdmix_prepare/build.sh`. The build uses the GU
Conda environment's HTSlib, with overrides `IBDMIX_HTSLIB_PREFIX`, `CXX`, and
`IBDMIX_PREPARE_BIN`. Production runs rebuild only when sources change.

## Validation

In the GU environment, run:

```bash
python3 f/ibdmix_prepare/test_prepare.py
```

The eight tests compare all five references against the previous bcftools plus
upstream generator, including GT missingness/phasing/ploidy, allele mismatches,
indels, indexed regions, EOF and GRCh37/38 X non-PAR boundaries. They also compare
CPU/GPU raw and final results on autosomes and X, and resume a partially completed
chromosome without rewriting completed reference calls. Invalid input must fail
without committed genotype/QC outputs.

On 2026-10-12, 4,564 real UKB samples at chr1:10,000,000–15,000,000 produced
byte-identical genotype streams and matching QC for all five references
(about 90,600 rows each). Preparation including QC took **73.17 s** for the old
five-pass pipeline and **18.34 s** for the shared pass, **3.99×** faster in this
single timing comparison. This measures preparation only, not complete chromosome
calling. Details and hashes: `validation-20261012.json`.

Console output is limited to chromosome boundaries, shared preparation, reference
completion/reuse and failure summaries. Full diagnostics remain in the main log;
the native preparer's diagnostics are in `logs/<unit>.prepare.log`.
