# GRID v0.2 — software and data requirements

This release replaces the previous `grid/` workflow with a small dispatcher and
module scripts under `f/`. The public `jielab/scripts` implementation present at
the time of release is preserved as `grid.legacy.unmodified.sh`; installation
also creates a timestamped backup of the entire destination directory before
writing any file.

## Software

Required for the core PRS/PCA/evaluation workflow:

- Linux or WSL2; Bash 4+, GNU coreutils, `awk`, `gzip`, `sort`, `sha256sum`.
- PLINK 2 (`plink2`).
- Python 3.10 or 3.11 with NumPy, SciPy, pandas, h5py, scikit-learn and tskit.
- R 4.x with `data.table`, `ggplot2`, `patchwork`, `openxlsx`, `pROC`, and
  `Matrix`.
- The PRS-CSx source is already bundled under `csx/`; its Python dependencies
  are installed by `install/environment.grid.yml`.

Additional requirements for ARG inference:

- ARG-Needle core and the official `arg-needle-scripts` Python package.
- A C/C++ compiler toolchain, CMake, Boost, GSL and Eigen, as required to build
  ARG-Needle and its Python extension.
- `bgenix`, needed to inspect/index UK Biobank phased BGEN files.

The supplied installer creates the Conda environment and installs ARG-Needle:

```bash
bash install/install_grid_env.sh
conda activate grid-argneedle
```

`tslmm`, `tsinfer`, and `tsdate` are not required by the primary GRID workflow.
`tslmm` is retained only as an optional future supervised ARG-BLUP comparator;
it is not an ARG inference program.

## Required data

### 1. Multi-ancestry GWAS summary statistics

Default directory:

```text
/mnt/d/data.BIG/gwas/4grid/
```

Expected flat names, case-insensitive:

```text
height.AFR.gz  height.EAS.gz  height.EUR.gz  height.SAS.gz
ldl.AFR.gz     ldl.EAS.gz     ldl.EUR.gz     ldl.SAS.gz
t2dm.AFR.gz    t2dm.EAS.gz    t2dm.EUR.gz    t2dm.SAS.gz
```

Each file should provide chromosome, position or rsID, effect/non-effect
alleles, effect estimate (`BETA` or `OR`), P value, and preferably SE, EAF and
sample size. For T2DM, odds ratios are converted to log odds ratios. Supply
`--n-gwas AFR=...,EAS=...,EUR=...,SAS=...` when sample size cannot be recovered
reliably from the files.

For unbiased UKB evaluation, every discovery GWAS used in the primary analysis
must exclude UK Biobank. The code records inputs but cannot prove cohort overlap
from the file contents.

### 2. PRS-CSx LD reference panels

Default:

```text
/mnt/d/data.BIG/refLD/csx/
```

The usual AFR/EAS/EUR/SAS PRS-CSx LD-block directories and SNP-information file
must be present. GRID also extracts ancestry-specific SNP LD scores directly
from these HDF5 blocks; no second LD panel is required for the baseline
transportability model.

### 3. UK Biobank target genotypes for PCA and PRS scoring

Default:

```text
/mnt/d/data/ukb/gen/typ/chr1 ... chr22
```

Per-chromosome PLINK BED/BIM/FAM or PGEN/PVAR/PSAM data are accepted according
to the module help. These are used for PCA projection and score calculation;
they are separate from the phased BGEN files used by ARG-Needle.

### 4. UK Biobank phased genotype BGEN files for ARG-Needle

Default WSL path corresponding to `H:\ukbGen\37\hap`:

```text
/mnt/h/ukbGen/37/hap/
```

Expected Field 22438 naming for each chromosome:

```text
ukb22438_c1_b0_v2.bgen
ukb22438_c1_b0_v2_s*.sample
...
ukb22438_c22_b0_v2.bgen
ukb22438_c22_b0_v2_s*.sample
```

The ARG module verifies that selected chromosomes have the same sample IDs in
the same order. It exports a deterministic, ancestry-balanced pilot subset to
Oxford HAPS/SAMPLE rather than exporting all approximately biobank-scale
participants by default.

### 5. GRCh37 genetic maps

Required for ARG-Needle. Prepare one sex-averaged cumulative genetic map per
autosome, covering the variants exported from the phased BGENs. The map parser
accepts common HapMap/1000 Genomes-style tables containing chromosome, physical
position and cumulative cM; `build_argneedle_map.py` validates monotonic bp and
cM coordinates and writes the exact ARG-Needle input map.

Recommended directory:

```text
/mnt/d/data.BIG/refGen/maps/GRCh37/
```

Pass a different location with `--arg-map-dir`.

### 6. Phenotypes and covariates

Default:

```text
/mnt/d/data/ukb/phe/Rdata/phe.rds
```

The RDS must contain `eid`. It should contain height, LDL and T2DM variables plus
age, sex and available PCs/covariates. The evaluator tries conservative name
aliases; use `--phenotype-col` when a trait name is ambiguous, and
`--eval-covariates` to specify the exact adjustment set.

### 7. PCA/ancestry reference assets

The existing DiscoDivas assets are bundled and preserved:

```text
DiscoDivas/files/g1k_hm3_maf5_woamb_wolr.pca.weight
DiscoDivas/files/med.g1000.4pop.tsv
```

The existing 1000 Genomes reference paths used by the previous PCA/ancestry
workflow must remain available. The `ancestry` module provides global ancestry
labels/posteriors; it is not a local-ancestry inference method.

## Optional data

- External allele-age estimates, such as a build-matched GEVA table. This is a
  sensitivity/augmentation input only. The primary implementation derives
  branch/mutation ages from the inferred ARG when conversion to tskit succeeds.
- An explicit target/exclusion ID file, strongly recommended for removing any
  discovery-overlap participants and close relatives.

## Storage and compute

Use a native Linux/WSL ext4 directory for `--arg-scratch`; do not put
ARG-Needle temporary/intermediate files under `/mnt/c`, `/mnt/d`, or `/mnt/h`.
The chromosome-22, 20,000-person pilot is intentionally the default. Measure
its runtime, peak RAM and disk use before increasing sample size or chromosome
coverage. A full 1–22 UKB ARG is an HPC-scale job and is blocked unless
`--arg-full TRUE` is supplied explicitly.
