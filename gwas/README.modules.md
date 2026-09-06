# WSL GWAS modules

`gwas.sh --help` includes the complete parameter list and runnable examples.
The entry point dispatches to `f/gwas_runner.py`; shell command plans and logs
are written beside each method's raw chromosome results. PLINK2, REGENIE and
SAIGE are independent modules. Internal REGENIE/SAIGE model fitting and testing
retain the software's own step terminology.

## Ready-to-run REGENIE commands

The phenotype file has been generated at `/mnt/d/data/ukb/phe/common/ukb.phe`.
To regenerate it:

```bash
bash /mnt/d/scripts/ukb/phe.sh --steps phe4gwas
```

First prepare the shared pruned BED and sparse GRM **once**. SAIGE is installed
and automatically selected by the script:

```bash
bash /mnt/d/scripts/gwas/gwas.sh prep_gwas --sparse-grm TRUE --threads 16
```

Then run all five traits:

```bash
bash /mnt/d/scripts/gwas/gwas.sh run_regenie --grch 37 --threads 16 \
  --pheno /mnt/d/data/ukb/phe/common/ukb.phe --pheno-name height,bald,bald12,cvd_cad.Yt2e,cvd_cad.t2e --type auto
```

The trait list and `--type auto` above are the defaults. REGENIE 4.1 was obtained
from its official release and is available at `/mnt/d/software/bin/regenie`.
An active environment's `regenie` takes precedence; use `--regenie PATH` to select
a particular executable. Time-to-event requires REGENIE >=4.

| Trait | Model | Complete phenotype/covariate cases |
|---|---|---:|
| height | Quantitative | 458073 |
| bald | Ordinal score treated as quantitative | 208282 |
| bald12 | Binary 0/1 | 114597 |
| cvd_cad.Yt2e | Binary 0/1 | 440615 |
| cvd_cad.t2e | Time-to-event, event column cvd_cad.Yt2e | 440615 |

These counts precede genotype overlap/QC. All methods use the same complete-case
numeric covariate design, including dummy-coded center. Constant/linearly
dependent covariates (e.g. sex in male-only bald data) are dropped. No rank inverse
normal transformation is applied, so quantitative BETA remains on the trait
scale. `bald` is not a proportional-odds model. Existing HPC BETA may differ if
its phenotype transformation, sample selection or covariates were different.

Default imputed data: `/mnt/h/ukbGen/<grch>/imp/chr<chr>.pgen/.psam/.pvar.zst`.
Shared pruning data remains in `/mnt/h/ukbGen/37/typ/typ.prune.{bed,bim,fam}`
for either build. Chromosomes 1–22 are the default. Raw association results stay
under `<output-dir>/<trait>/gwas/<method>/`; the canonical combined file is
`<output-dir>/<trait>/gwas/<trait>.<method>.gz`, with columns
`SNP CHR POS EA NEA EAF N BETA SE P` and a `.grch` sidecar. Very small P is capped
at 1e-300. Invalid result rows are counted in `.qc.txt`.

Use `--run FALSE` to generate plans only. Finished jobs are reused if commands,
input sizes/timestamps, and expected outputs match. Different builds/subsets
should use separate output roots. `--replace TRUE` explicitly permits rerunning.
Generated `run.cmd.sh` executes the corresponding JSON plan with failure checks.

## Sparse GRM and SAIGE

SAIGE 1.3.1 is installed in `/home/huangj/anaconda3/envs/saige` with R 4.4.3.
The script uses `f/saige_Rscript.sh` to isolate it from global R libraries.
RcppParallel is pinned to 5.1.9 for the binary's classic-TBB ABI; the environment
specification is in `f/saige.environment.yml`. No manual activation is required.
Rerun prep to add the sparse GRM; successfully completed pruning jobs resume
without recomputation when inputs and other settings are unchanged:

```bash
bash /mnt/d/scripts/gwas/gwas.sh prep_gwas --sparse-grm TRUE --threads 16
```

This produces stable names `typ.sparseGRM.mtx` and
`typ.sparseGRM.mtx.sampleIDs.txt`, linked to SAIGE's actual output names.
Prep excludes the build-37 MHC 25–35 Mb region, uses autosomal common markers,
and merges only the explicitly selected pruned BED files. Its default includes
GRM generation; `--sparse-grm FALSE` is the explicit REGENIE-only option.
Pruning uses `--hwe 1e-15 0`: the explicit zero retains a fixed threshold while
acknowledging PLINK2's large-sample HWE strictness check. A failed pruning job has
no success checkpoint and reruns; merely having an output directory is insufficient.

`run_saige` uses a sparse null model with LOCO disabled consistently in both
stages. PGEN is exported to dosage VCF with forced DS fields, preserving dosage
rather than silently replacing it with hard calls. These VCFs can be large.
PLINK2 and this SAIGE configuration skip T2E with an `UNSUPPORTED.txt` notice.
SAIGE sparse GRM generation and FALSE-to-TRUE resume have passed an actual
500-sample synthetic test. Full-UKB GRM generation has not been started here.

## Compare

```bash
bash /mnt/d/scripts/0data/gwas_post.sh compare \
  --gwas-files /mnt/d/data.BIG/gwas/main/common/bald/gwas/bald.gz,/mnt/d/data.BIG/gwas/self/common/bald/gwas/bald.regenie.gz \
  --mplot TRUE --compare-beta TRUE --compare-EAF TRUE \
  --output-dir /mnt/d/data.BIG/gwas/self/common/bald/qc/compare
```

Add any number of comma-separated GWAS files. Each follower is compared with
the first. `--p-threshold 5e-8 --significant first` is the default; `either` and
`both` select the union or intersection of significance. Position and alleles
are matched; simple SNP strand complements and swapped alleles are aligned.
Swapping changes BETA's sign and EAF to `1-EAF`. Palindromic SNPs and duplicate
allele/position keys are excluded from scatter plots. Input and paired QC TSVs,
harmonized data and PNGs are written to the output directory. Manhattan tracks
share genomic coordinates and y-axis limits; background points are thinned in
20kb/0.1-logP bins, with all significant points retained. GWAS builds must match;
conflicting `.grch` metadata is rejected. A legacy file without metadata must be
checked by the caller. Inputs are read one genome-wide file at a time; allow
sufficient memory for large summary-statistic tables.

## LDSC

```bash
bash /mnt/d/scripts/0data/gwas_post.sh ldsc \
  --gwas-files /path/A.gz,/path/B.gz --output-dir /path/ldsc-results \
  --merge-alleles /path/w_hm3.snplist \
  --ref-ld-chr /path/eur_w_ld_chr/ --w-ld-chr /path/eur_w_ld_chr/
```

The default software directory is `/mnt/d/software/ldsc`, run via conda environment
`ldsc`. `--python PATH` overrides conda. Munging, h2 and all pairwise rg are enabled
by default; `--run FALSE` only writes the commands. An explicit `--N` is required
if N is missing. The original GWAS is preserved; invalid P is removed and zero P
is capped in a separate munging input. Real LD-score reference files are required
for h2/rg and were not available for an end-to-end regression test here.

## Validation

`python3 /mnt/d/scripts/gwas/tests/test_modules.py` has passed:

- Four real PLINK2 analyses on synthetic data, then canonical result merging.
- Five real REGENIE 4.1 analyses, including Cox/T2E, model fitting and testing.
- Compressed PVAR handling, IID-only PSAM adaptation, actual pruning/BED merging.
- Effect-allele normalization for all three methods, signed BETA/EAF alignment,
  duplicate/palindromic exclusions, N-way plots and real LDSC munging.
- Real phenotype generation and complete-case checks (aggregate counts above).

Additional tests: `tests/test_ukb_hwe.py` checks the corrected HWE flags using
all 488377 typed UKB samples and 100 variants (no association analysis);
`tests/test_prep_resume.py` runs actual SAIGE GRM generation and asserts that
prune/merge output timestamps do not change when switching FALSE to TRUE.
The preferred CLI names are `--pheno` and `--pheno-name`; the previous
`--phenotype-file` and `--phenotypes` spellings remain compatible aliases.

Full UKB GWAS and full-data pruning/GRM have not been started. Runtime logs for
the synthetic tests are in `tests/test_modules.log`.

Method interfaces follow the primary documentation:
[REGENIE](https://rgcgithub.github.io/regenie/options/),
[SAIGE](https://saigegit.github.io/SAIGE-doc/docs/single_step1.html),
[LDSC](https://github.com/bulik/ldsc).
