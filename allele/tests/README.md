# Validation

On 2026-10-15, all **30 offline regression tests passed** in 57.0 seconds; see [last_run.txt](last_run.txt). Run them again with:

```bash
cd /mnt/d/scripts/allele
./allele.sh test
```

The installed conda environment and dependencies passed `./allele.sh inspect`. Reinstallation is not needed on this machine. Tests create temporary data and do not download files or call a real API.

The suite covers:

- The public `allele.sh all` command through four figures and six tables, using synthetic variants and an explicitly simulated Atlas client.
- Cohort VCF parsing, sample-column skipping while retaining INFO AC/AN/AF, indexed region boundaries, mixed indel/SNV indexing, zero counts, REF mismatch and empty input.
- Streaming across batch/file boundaries, merging split alternatives, rejecting unsorted input, avoiding partial output publication and refusing to overwrite existing prepare results.
- Exact weighted feature-count training versus expanded-row logistic regression, locus/block validation, primary-label exclusion, unavailable statistical fits and single-chromosome behavior.
- Atlas variant alignment, shuffled responses, interruption/resume, cache identity, nonfinite and zero vectors, known geometry values.
- Required gnomAD chromosome sets, unrelated-file exclusion, truncated BGZF and pending `.aria2` detection, and read-only input checks.
- Seven standalone downloader tests: chromosome lists, official checksum manifests, proxy overrides, resume flags, invalid arguments, failure propagation/cleanup, standalone copies and removal of download options from `install.sh`.

## Real 1KG data

Test data are under `/mnt/d/analysis/allele_1kg_validation_20261015`, separate from formal gnomAD output. All 22 local autosomal VCFs passed file/index readiness checks. Original chr1, chr2 and chr22 headers each contain 3,202 samples.

The complete original `chr22.vcf.gz` passed `prepare → select → stage1`: 1,066,557 VCF records yielded **914,117 loci**, **2,742,351 candidates**, and **34,353 selected loci**. There were 902,666 biallelic, 11,297 triallelic and 154 tetra-allelic loci. Extraction used about 240 MiB peak RSS in 5m20s; selection took 9.8s, and Stage 1 took 22.4s with about 933 MiB peak RSS. These timings include concurrent local I/O and are not whole-genome estimates. All 77,094 candidate rows in the previously tested chr22 1 Mb interval match the earlier extraction in every column. Results and timing logs are in `chr22_full_verified/` and `chr22_*_verified.log`.

The multi-chromosome integration test uses `chr1`, `chr2` and `chr22`, each at `20,000,001–21,000,000`. `bcftools view --no-update --drop-genotypes` retains the original cohort INFO values in small sites-only test VCFs. This exercises directory input in the same shape as a sites-summary dataset, without calling the data gnomAD.

The test produced **69,720 loci**, **209,160 candidates**, **2,199 selected loci**, four figures and six tables. Occupancy counts were 68,987 biallelic, 728 triallelic and five tetra-allelic sites. All seven pipeline stages passed; Stage 1 used chromosome-blocked cross-fitting (AUC 0.81836).

**Atlas vectors in this integration test are simulated.** `output_FAKE_ATLAS/` and `NOTICE.txt` label these artifacts accordingly. This verifies software wiring on real cohort variants, not actual molecular effects or biological associations. A separate live probe subsequently authenticated successfully and queried one public 1KG locus (three alternatives) across all six default scorers, returning 3 × 3,763 finite, correctly aligned values. Its record is `/mnt/d/analysis/allele_1kg_validation_20261015/live_api_check.json`. Full real-data Atlas integration and whole-genome throughput remain untested.

To repeat this manual integration test in a new directory:

```bash
/home/huangj/anaconda3/bin/conda run --no-capture-output -n allele \
  python /mnt/d/scripts/allele/tests/real_1kg_smoke.py \
  --out /mnt/d/analysis/allele_1kg_new_validation
```

The earlier 1 Mb chr22 pilot is retained in [real_1kg_pilot.json](real_1kg_pilot.json) for historical comparison; it used the previous Stage 1 fold assignment. See [real_1kg_validation.json](real_1kg_validation.json) for the current validation record, including full chr22 extraction/selection/Stage 1 and measured resource usage.

## Synthetic artifacts

By default fixtures are removed automatically. To retain the synthetic smoke run, choose a **new** directory:

```bash
ALLELE_TEST_OUT=/path/to/new/allele-smoke ./allele.sh test
```

The retained directory contains `input/`, `output/` and `pipeline.log`; `output/final/` has the synthetic figures/tables. Do not interpret them as research findings.
