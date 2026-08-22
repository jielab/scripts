# GU: archaic introgression workflow

`gu` (“古”) is an integrated workflow for archaic introgression analysis. It combines
selected-locus evidence with three individual-level segment callers and converts their
outputs into a common, UKB-friendly data model.

The long-term endpoint is deliberately similar to UK Biobank ICD-10 long-format data:

```text
sample_id   segment_code            source_class     dosage
1234567     GU37-NEA-3-7A12F08C31  Neanderthal     1
1234567     GU37-GST-7-9912A1408C   Ghost/Unknown   1
2345678     GU37-DEN-12-4ED76000B2  Denisovan       2
```

Each participant can therefore carry **0-N introgressed-segment codes**, which can be
used in carrier, burden, survival, phenome-wide, proteomic or other association analyses.

Current GU version: see [`VERSION`](VERSION).

---

## 1. Methods and their intended roles

| GU method | Scale | Archaic reference | ARG | Main role |
|---|---|---:|---:|---|
| `loci_avcf` | selected GWAS loci | required | no | locus-level haplotype/phylogenetic evidence |
| `loci_asnp` | selected GWAS loci | published archaic SNP/haplotype map | no | fast screening / sensitivity analysis |
| `ibdmix` | individual, genome-wide | required | no | known-donor Neanderthal/Denisovan tracts |
| `trace` | individual, genome-wide | **not required** | **required** | deep/ghost/unknown archaic ancestry |
| `as3` | individual, genome-wide | archaic + African reference | no | scalable known-archaic classification |

The two `loci_*` methods are **not** treated as equivalent to whole-genome individual
callers. Their outputs are stored as locus evidence. `ibdmix`, `trace` and `as3` feed the
per-person segment database.

### Why keep several methods?

They answer different questions.

- **IBDmix** asks whether a modern segment is unusually identical-by-descent with a
  sequenced archaic genome.
- **TRACE** uses genealogy/ARG features and therefore can detect deep archaic ancestry
  even when the donor genome is not available.
- **AS3.1-mamba** is a modern deep-learning caller for known archaic classes.
- **loci_avcf** is useful when the question is a specific GWAS haplotype, such as the
  COVID-19 Neanderthal risk haplotype paradigm.
- **loci_asnp** is useful as a fast map-based sensitivity analysis, but should not be
  treated as a de novo genome-wide caller.

---

## 2. The most important TRACE rule

A finished ARG containing only 1000 Genomes samples **cannot directly call a new UKB
participant**.

TRACE extracts features for haplotype/sample nodes that are already present in the input
tree sequence. A reusable reference or anchor genealogy can help define the background,
but each target UKB haplotype must still be represented in the **final genealogy used by
TRACE**.

For UKB batching, GU therefore supports the following experimental design:

```text
fixed UKB anchor panel + target batch
                |
                v
             joint ARG
                |
                v
         TRACE on targets only
```

`TRACE_TARGET_SAMPLES` tells GU which people are targets. The anchors affect the ARG but
are not emitted as target TRACE calls.

Do **not** build hundreds of completely unrelated 1,000-person ARGs and assume their
TRACE posteriors are automatically exchangeable. Fixed anchors reduce this problem but
do not eliminate it; batch effects must be evaluated empirically.

---

# 3. Directory structure

The repository is self-contained under `gu/`:

```text
gu/
├── README.md
├── VERSION
├── gu.sh
├── install.sh
├── environment.yml
├── f/
│   ├── add_aa_from_1kg.sh
│   ├── arg.sh
│   ├── arg_tsinfer.py
│   ├── arg_vcf_prep.sh
│   ├── as3.sh
│   ├── as3_plan.py
│   ├── as3_prep.sh
│   ├── gen.clean.sh
│   ├── ibdmix.R
│   ├── ibdmix.sh
│   ├── ibdmix_legacy.sh
│   ├── loci.R
│   ├── loci.sh
│   ├── make_batches.py
│   ├── make_sample_map.py
│   ├── normalize_results.py
│   ├── sample_panel.py
│   ├── trace.sh
│   ├── trace_combine.py
│   ├── trace_report_legacy.py
│   ├── ukb.sh
│   └── vcf_gt_fix.py
└── shiny/
    ├── global.R
    ├── run_shiny.R
    ├── server.R
    └── ui.R
```

---

# 4. Installation

## 4.1 Linux / WSL assumptions

The scripts use Linux/WSL paths. Example Windows drive mappings are:

```text
D:\data.BIG       -> /mnt/d/data.BIG
H:\ukbGen\hap    -> /mnt/h/ukbGen/hap
H:\ukbGen\typ    -> /mnt/h/ukbGen/typ
```

A working Conda/Miniconda/Miniforge installation is required.

## 4.2 Clone

```bash
git clone https://github.com/jielab/pub.git
cd pub/gu
chmod +x gu.sh install.sh f/*.sh
```

If the directory was copied manually instead of cloned, simply enter the local `gu`
directory and run the same `chmod` command.

## 4.3 Install and check the single GU environment

```bash
./install.sh
```

The single environment installs the tools used by every current pipeline, including:

- Python 3.11;
- `bcftools`, `htslib`, `bedtools`, `samtools`;
- PLINK 2;
- `bgenix` for UKB BGEN inspection/indexing;
- `tskit`, `tsinfer`, `tsdate`, `bio2zarr`;
- TRACE;
- IBDmix;
- R + `data.table`, Shiny, Plotly, DT, DBI and RSQLite.

The environment pins `tsinfer>=0.5.1` and `tsdate>=0.2.7`, which are the current 2026
releases used by this code path.

## 4.4 AS3

`./install.sh` also installs ArchaicSeeker3.1-mamba into the same `gu` environment.
It invokes the upstream installer with `AS3_ENV_NAME=gu`, so no second Conda
environment or second environment YAML is created. The upstream repository is:
The official AS3 repository is:

```text
https://github.com/Shuhua-Group/ArchaicSeeker3.0
```

## 4.5 Optional ARG software

SINGER, Threads and ARG-Needle are **not required** for the default GU run. They are
important sensitivity/scaling candidates, but an ARG backend should not be considered
TRACE-valid merely because it can produce a tree.

Optional software can be installed separately when a non-default ARG backend is needed.

## 4.6 Configure paths

```bash
export GU_CHRS=22
```

The default development setup assumes:

```text
GU_REF_ROOT        /mnt/d/data.BIG/refGen
GU_ANALYSIS_ROOT   /mnt/d/analysis/gu
1KG GRCh37         /mnt/d/data.BIG/refGen/1kg/GRCH37
archaic GRCh37     /mnt/d/data.BIG/refGen/archaic/GRCH37
UKB phased BGEN    /mnt/h/ukbGen/hap
UKB typed PLINK    /mnt/h/ukbGen/typ
```

All defaults are built into `gu.sh`; one-off shell exports override them. An optional
local `gu.env` beside `gu.sh` is still sourced when present.

## 4.7 Check installation

```bash
./install.sh
```

The first call installs or repairs everything and performs a full health check. Later
calls perform the same check and print only `OK` when the YAML and installation are current.

---

# 5. Recommended validation order

Do not begin with chrX or the full UKB cohort.

Recommended sequence:

```text
1KG chr22
   -> IBDmix
   -> tsinfer/tsdate ARG
   -> TRACE
   -> AS3 (with appropriate official resources)
   -> normalize
   -> Shiny

then

UKB Field 22438 chr22, one anchor+target batch
   -> phased VCF
   -> AA transfer
   -> ARG
   -> TRACE target-only
   -> batch sensitivity

then

multiple UKB chr22 batches

then additional autosomes
```

chrX is intentionally left for a later dedicated analysis because PAR/non-PAR,
hemizygosity, sex-specific ploidy and X-specific evolutionary history require separate
handling.

---

# 6. 1000 Genomes test workflow

Set the development target:

```bash
export GU_TARGET=1kg
export GU_BUILD=37
export GU_TARGET_ROOT=/mnt/d/data.BIG/refGen/1kg/GRCH37
export GU_TARGET_VCF_DIR=$GU_TARGET_ROOT/vcf
export GU_SAMPLE_PANEL=$GU_TARGET_VCF_DIR/samples_v3.ALL.panel
export GU_ARCHAIC_ROOT=/mnt/d/data.BIG/refGen/archaic/GRCH37
export GU_CHRS=22
```

If the modern and archaic reference inputs have not yet been standardized:

```bash
./install.sh
```

## 6.1 `loci_avcf`

Example:

```bash
traits="bald bald12" ./gu.sh loci_avcf run
```

Main output root:

```text
$GU_ANALYSIS_ROOT/loci_avcf/
```

Use this for locus-specific modern/archaic haplotype and phylogenetic evidence.

## 6.2 `loci_asnp`

```bash
traits="bald bald12" ./gu.sh loci_asnp run
```

Main output root:

```text
$GU_ANALYSIS_ROOT/loci_asnp/
```

The published archaic-SNP/haplotype map must be present where configured by the locus
pipeline. Treat this method as prior-map evidence/sensitivity analysis.

## 6.3 IBDmix

Check inputs:

```bash
./gu.sh ibdmix check
```

Run:

```bash
./gu.sh ibdmix run
```

Rebuild reports without rerunning the expensive scan:

```bash
./gu.sh ibdmix report
```

Main results:

```text
$GU_ANALYSIS_ROOT/ibdmix/final/
$GU_ANALYSIS_ROOT/ibdmix/report/
```

The key object is the per-person interval table containing sample, chromosome,
start/end, archaic source/reference and IBDmix evidence.

## 6.4 Prepare an ARG for TRACE

The default engineering backend is `tsinfer + tsdate`:

```bash
export GU_ARG_BACKEND=tsinfer
export GU_ARG_DIR=$GU_TARGET_ROOT/arg
export GU_ARG_VCF_DIR=$GU_TARGET_ROOT/vcf.4arg
export GU_CHRS=22

./gu.sh arg prepare
./gu.sh arg build
./gu.sh arg check
```

Main outputs:

```text
$GU_ARG_DIR/chr22.trees
$GU_ARG_DIR/chr22.sample_map.tsv
$GU_ARG_DIR/chr22.arg_qc.json
```

### Ancestral allele handling

For 1KG, `arg_vcf_prep.sh` uses `INFO/AA` as the ancestral allele. Values such as
`A|||` are normalized to `A`. Sites with missing/invalid ancestral states are removed.
GU does **not** silently assume that `REF` or the major allele is ancestral.

### External ARGs

For sensitivity validation, an external full-chromosome `.trees/.tsz` ARG may be used:

```bash
export GU_ARG_BACKEND=external
export GU_ARG_DIR=/path/to/validated/trees
./gu.sh arg check
```

TRACE's official examples include SINGER and Relate ARG workflows. SINGER is especially
useful when multiple posterior ARG samples are needed.

**Important:** files such as `chr1_part1.trees` and `chr1_part2.trees` are genomic
coordinate chunks, not posterior samples of the same region. GU refuses to treat such
chunks as posterior ARG replicates.

## 6.5 TRACE

```bash
export TRACE_CHRS=22
./gu.sh trace run
```

Restart individual stages when required:

```bash
./gu.sh trace extract
./gu.sh trace infer
./gu.sh trace segments
./gu.sh trace report
```

Main results:

```text
$GU_ANALYSIS_ROOT/trace/final/trace_haplotype_segments.tsv.gz
$GU_ANALYSIS_ROOT/trace/report/
```

TRACE calls enter the unified database as `Ghost/Unknown` (`GST`) unless a later analysis
provides donor-specific evidence.

## 6.6 ArchaicSeeker3.1-mamba

Prepare the official build-matched resources, validate, then run:

```bash
./gu.sh as3 prep
./gu.sh as3 check
./gu.sh as3 run
```

Main output root:

```text
$GU_ANALYSIS_ROOT/as3/
```

The primary AS3 output recognized by GU is:

```text
introgression_prediction.bed
```

Current official AS3 coding is interpreted as:

```text
1 = Denisovan
2 = Neanderthal
3 = Mosaic
```

Do not assume a sparse GRCh37 UKB array-marker VCF is equivalent to the input
distribution used to validate the official AS3 model.

---

# 7. UKB Field 22438 phased haplotypes

## 7.1 What these local files are

The available files are:

```text
/mnt/h/ukbGen/hap/
  ukb22438_c1_b0_v2.bgen
  ukb22438_c1_b0_v2_s487162.sample
  ...
  ukb22438_c22_b0_v2.bgen
  ukb22438_c22_b0_v2_s487162.sample
```

These are UK Biobank **Field 22438, Haplotypes (WTCHG)**. They are not WGS haplotypes.
They are the phased genotype-array/phasing-marker release.

UKB/BGEN documentation states that:

- QC-passing genotypes were statistically phased using SHAPEIT3;
- phased calls were stored as hard-called phased genotypes;
- phasing chunks were ligated, so phase can be treated as consistent across each
  chromosome;
- the coordinates are GRCh37;
- the historical full-release phased dataset contained 487,409 post-QC individuals;
- the phased BGEN chromosome field is blank due to a historical processing issue, so
  the chromosome must be supplied/read from the single-chromosome file context/ID.

Useful references:

- UKB Field 22438: https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=22438
- BGEN in UK Biobank: https://enkre.net/cgi-bin/code/bgen/wiki?name=BGEN+in+the+UK+Biobank
- PLINK 2 BGEN input: https://www.cog-genomics.org/plink/2.0/input

The local `.sample` file is the authority for the samples actually present in the local
archive. GU does not infer the local sample count from today's UKB Showcase count.

## 7.2 Why use `/hap` instead of `/typ` for TRACE

`/mnt/h/ukbGen/typ` contains the older BED/BIM/FAM directly genotyped representation.
It is useful for selected-locus analyses and some ARG experiments, but phase is not
represented in ordinary PLINK1 BED/BIM/FAM.

`/mnt/h/ukbGen/hap` already contains phased haplotypes, so it is the preferred local
input for the current UKB ARG/TRACE experiment.

## 7.3 Required GRCh37 FASTA

For BGEN -> VCF conversion, GU does **not** assume BGEN allele 1 is the reference allele.
Set a build-matched GRCh37 FASTA:

```bash
export UKB_REF_FASTA=/path/to/human_g1k_v37.fasta
```

The file should be readable by PLINK 2; keeping a `.fai` index is also recommended.

GU imports UKB BGEN with `ref-unknown` and uses PLINK 2 `--ref-from-fa` to establish REF
where possible before any ancestral-allele annotation is transferred.

## 7.4 Inspect chr22 first

```bash
export GU_BUILD=37
export GU_CHRS=22
export UKB_HAP_ROOT=/mnt/h/ukbGen/hap

./gu.sh ukb inspect-hap
```

This:

1. finds the BGEN and matching Oxford `.sample` file;
2. creates `.bgen.bgi` when absent;
3. counts samples and variants;
4. hashes the sample IDs in BGEN order;
5. when multiple chromosomes are requested, verifies that sample IDs/order are the same.

Report:

```text
$GU_ANALYSIS_ROOT/ukb/hap.inspect.tsv
```

Inspect all autosomes only after chr22 works:

```bash
export GU_CHRS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
./gu.sh ukb inspect-hap
```

## 7.5 Create a GU UKB sample panel

If an ancestry-labelled UKB panel is not yet available:

```bash
export GU_SAMPLE_PANEL=$GU_ANALYSIS_ROOT/ukb/ukb.sample_panel.tsv
export GU_CHRS=22
./gu.sh ukb make-panel
```

Output columns:

```text
sample  pop  super_pop  sex
```

The automatically created panel uses:

```text
pop       = UKB
super_pop = ALL
```

and uses sex from the Oxford `.sample` file when available.

For the final large-scale study, replace `ALL` with a validated ancestry classification
before creating batches.

## 7.6 Create target batches and fixed UKB anchors

### First ALL-only experiment

```bash
export GU_BATCH_SIZE=1000
export GU_ANCHORS_PER_GROUP=1000
./gu.sh ukb batches
```

Output:

```text
$GU_ANALYSIS_ROOT/ukb/batches/
  anchors.samples.txt
  ALL.b0001.targets.txt
  ALL.b0001.joint.txt
  ALL.b0002.targets.txt
  ALL.b0002.joint.txt
  ...
  batch_manifest.tsv
```

For this ALL-only example:

```text
anchors.samples.txt     1000 fixed UKB anchor people
*.targets.txt           1000 target people per batch
*.joint.txt             anchors + targets
```

Anchors are removed from the target pool.

### An ancestry-labelled panel

If `super_pop` contains EUR/EAS/SAS/AFR/AMR and you set:

```bash
export GU_ANCHORS_PER_GROUP=200
```

GU selects 200 fixed anchors from each ancestry group (about 1,000 total with five
nonempty groups). Every batch uses the same complete anchor panel.

### Curated anchor list

```bash
export GU_ANCHOR_LIST=/path/to/fixed_ukb_anchor_ids.txt
export GU_ANCHORS_PER_GROUP=0
./gu.sh ukb batches
```

**The current `GU_ANCHOR_LIST` must contain UKB IDs present in Field 22438.** It is not a
1KG sample list. Adding 1KG haplotypes to the same ARG requires an explicit harmonized
multi-cohort ARG workflow and is not silently performed by GU.

## 7.7 Prepare one UKB batch: BGEN -> phased VCF -> INFO/AA

Start with the first batch only:

```bash
B=ALL.b0001

export GU_CHRS=22
export UKB_HAP_ROOT=/mnt/h/ukbGen/hap
export UKB_REF_FASTA=/path/to/human_g1k_v37.fasta

export UKB_KEEP=$GU_ANALYSIS_ROOT/ukb/batches/${B}.joint.txt
export UKB_VCF_OUT=$GU_ANALYSIS_ROOT/ukb/work/${B}/vcf
export UKB_ARG_VCF_OUT=$GU_ANALYSIS_ROOT/ukb/work/${B}/vcf.aa
export UKB_1KG_VCF_DIR=/mnt/d/data.BIG/refGen/1kg/GRCH37/vcf

./gu.sh ukb hap-arg-vcf
```

This performs two stages.

### Stage A: Field 22438 BGEN -> phased VCF

PLINK 2 imports the BGEN as `ref-unknown`, supplies chromosome 22 explicitly, subsets to
the requested UKB IDs, restricts to biallelic A/C/G/T SNPs, sets REF from the GRCh37
FASTA, and exports BGZF VCF.

GU then samples heterozygotes and aborts if unphased `0/1` calls are found where phased
`0|1` / `1|0` calls are expected.

Main output:

```text
$GU_ANALYSIS_ROOT/ukb/work/ALL.b0001/vcf/chr22.vcf.gz
```

### Stage B: transfer ancestral allele from 1KG

Field 22438 does not contain the `INFO/AA` tag needed by the current tsinfer preparation.
GU transfers AA from a **build-matched 1KG GRCh37 VCF only when all four fields match**:

```text
CHROM
POS
REF
ALT
```

There is no strand guess, allele-complement guess or liftover.

Main output:

```text
$GU_ANALYSIS_ROOT/ukb/work/ALL.b0001/vcf.aa/chr22.vcf.gz
$GU_ANALYSIS_ROOT/ukb/work/ALL.b0001/vcf.aa/chr22.aa_qc.tsv
```

If few sites receive AA, first check:

- GRCh37 vs GRCh38 mismatch;
- `22` vs `chr22` chromosome naming;
- wrong REF orientation/reference FASTA;
- nonmatching 1KG reference release.

## 7.8 Build the batch ARG

Point GU to this batch:

```bash
B=ALL.b0001

export GU_TARGET=ukb
export GU_BUILD=37
export GU_TARGET_ROOT=$GU_ANALYSIS_ROOT/ukb/work/${B}
export GU_TARGET_VCF_DIR=$GU_TARGET_ROOT/vcf.aa
export GU_ARG_VCF_DIR=$GU_TARGET_ROOT/vcf.4arg
export GU_ARG_DIR=$GU_TARGET_ROOT/arg
export GU_ARG_BACKEND=tsinfer
export GU_CHRS=22

./gu.sh arg prepare
./gu.sh arg build
./gu.sh arg check
```

Output:

```text
$GU_TARGET_ROOT/arg/chr22.trees
$GU_TARGET_ROOT/arg/chr22.sample_map.tsv
$GU_TARGET_ROOT/arg/chr22.arg_qc.json
```

The ARG contains the fixed anchors **and** target participants.

## 7.9 Run TRACE on targets only

```bash
B=ALL.b0001

export TRACE_CHRS=22
export TRACE_TARGET_SAMPLES=$GU_ANALYSIS_ROOT/ukb/batches/${B}.targets.txt
export TRACE_OUT=$GU_ANALYSIS_ROOT/ukb/trace/${B}

./gu.sh trace run
```

TRACE uses the complete joint ARG but only target nodes listed in
`TRACE_TARGET_SAMPLES` are sent through target inference/reporting.

Main result:

```text
$GU_ANALYSIS_ROOT/ukb/trace/ALL.b0001/final/trace_haplotype_segments.tsv.gz
```

Repeat with several batches before scaling further:

```bash
B=ALL.b0002
# repeat BGEN -> AA -> ARG -> TRACE with the same anchors
```

Compare between batches:

- ghost bp/person;
- number of tracts/person;
- tract length distribution;
- posterior distribution;
- genomic hotspots;
- duplicated/overlapping validation individuals if deliberately included;
- sensitivity to anchor composition.

Only after these checks should the analysis be expanded beyond chr22.

---

# 8. IBDmix with UKB Field 22438: sensitivity route

IBDmix does not require an ARG, but Field 22438 is much sparser than WGS. Therefore GU
supports an array-marker experiment but does **not** label it automatically as
production-equivalent to dense-sequence IBDmix.

Use a target-only phased VCF rather than repeating anchors as target calls:

```bash
B=ALL.b0001
export GU_CHRS=22
export UKB_KEEP=$GU_ANALYSIS_ROOT/ukb/batches/${B}.targets.txt
export UKB_VCF_OUT=$GU_ANALYSIS_ROOT/ukb/ibdmix_work/${B}/vcf
./gu.sh ukb hap-vcf
```

Then point the IBDmix wrapper at the batch root:

```bash
export GU_TARGET=ukb
export GU_BUILD=37
export GU_TARGET_ROOT=$GU_ANALYSIS_ROOT/ukb/ibdmix_work/${B}
export GU_TARGET_VCF_DIR=$GU_TARGET_ROOT/vcf
export IBDMIX_OUT=$GU_ANALYSIS_ROOT/ukb/ibdmix/${B}

./gu.sh ibdmix check
./gu.sh ibdmix run
```

The resulting UKB batch outputs are recognized recursively by `./gu.sh normalize`.

Before a cohort-wide run, benchmark array-marker IBDmix against dense 1KG calls after
artificially thinning the 1KG data to a comparable marker set.

---

# 9. AS3 and the local UKB phased BGEN data

The current official ArchaicSeeker3.1-mamba repository provides GRCh38/CHM13-oriented
resources and expects phased target VCFs. The local Field 22438 files are GRCh37 and
array-marker density.

Therefore:

- AS3 remains integrated into GU;
- it is appropriate for a build/resource combination that matches the official model;
- the local Field 22438 GRCh37 data should be treated as a sensitivity/engineering input,
  not silently as a model-equivalent production dataset.

Do not mix genome builds merely to make all five methods run.

---

# 10. Normalize all methods into the ICD-like data model

After any subset of methods has completed:

```bash
./gu.sh normalize
```

The normalizer reads both ordinary single-run outputs and UKB batch layouts such as:

```text
$GU_ANALYSIS_ROOT/ukb/trace/*/final/
$GU_ANALYSIS_ROOT/ukb/ibdmix/*/final/
$GU_ANALYSIS_ROOT/ukb/as3/**/
```

Main products:

```text
$GU_ANALYSIS_ROOT/Rshiny/
├── gu.sqlite
├── segments.tsv.gz
├── segment_catalog.tsv.gz
├── carriers.tsv.gz
├── sample_burden.tsv.gz
└── loci.tsv.gz
```

## 10.1 `segments.tsv.gz`

One method-level raw call per sample/haplotype/segment after harmonization:

```text
sample_id
method
source
source_class
chr
start
end
length_bp
haplotype
score
posterior
trait
locus_id
genome_build
batch_id
raw_file
segment_code
```

Raw tract coordinates are retained.

## 10.2 `segment_catalog.tsv.gz`

Exact boundaries differ between samples and callers. GU therefore clusters strongly
overlapping intervals within the same source class/chromosome and generates a canonical
segment code such as:

```text
GU37-NEA-3-7A12F08C31
GU37-DEN-12-4ED76000B2
GU37-GST-7-9912A1408C
```

Abbreviations:

```text
NEA = Neanderthal
DEN = Denisovan
MOS = Mosaic
GST = Ghost/Unknown
ARC = other archaic
```

For a formal frozen UKB codebook, first establish and validate the catalog on a discovery
set, then freeze/version it before final association analyses. The current automatic
catalog builder is intended for discovery and data harmonization, not as an immutable
international ontology.

## 10.3 `carriers.tsv.gz`

This is the table most analogous to UKB ICD-10 long format:

```text
sample_id
segment_code
source_class
dosage
n_methods
methods_support
max_score
max_posterior
```

One participant naturally has 0-N rows.

## 10.4 `sample_burden.tsv.gz`

Per-person summary by method/source class:

```text
sample_id
method
source_class
n_segments
total_bp
mean_length_bp
max_score
max_posterior
```

This table is convenient for global Neanderthal/Denisovan/ghost burden association.

## 10.5 `loci.tsv.gz`

Harmonized locus-level evidence from `loci_avcf` and `loci_asnp`.

## 10.6 `gu.sqlite`

SQLite is the primary backend for Shiny and large interactive queries. It avoids loading
all UKB segment rows into R memory.

---

# 11. Interactive visualization with R Shiny

After normalization:

```bash
./gu.sh shiny
```

Default local address:

```text
http://127.0.0.1:3838
```

Override when needed:

```bash
export GU_SHINY_HOST=0.0.0.0
export GU_SHINY_PORT=3838
./gu.sh shiny
```

The app reads `GU_SQLITE` (default `$GU_ANALYSIS_ROOT/Rshiny/gu.sqlite`) lazily.

Main views include:

1. **Overview** — genome-wide segment density/carrier burden;
2. **Individual** — all archaic segments for one 1KG/UKB sample, by chromosome,
   method, source and haplotype;
3. **Region / GWAS locus** — interactive regional interval view;
4. **Locus methods** — `loci_avcf` and `loci_asnp` evidence;
5. **Segment codes** — canonical GU codes and carriers;
6. **Method concordance** — interval-level overlap between whole-genome callers;
7. **Data** — searchable/downloadable normalized records.

This visualization is designed around **introgressed intervals, methods, source classes,
haplotypes, scores/posteriors and carrier frequencies**, rather than sequencing read
alignments, which is why it is more suitable for this project than IGV alone.

---

# 12. Useful command summary

## Installation

```bash
./install.sh
./install.sh --check
```

## 1KG chr22

```bash
export GU_CHRS=22
./gu.sh ibdmix check
./gu.sh ibdmix run
./gu.sh arg build
export TRACE_CHRS=22
./gu.sh trace run
./gu.sh normalize
./gu.sh shiny
```

## UKB Field 22438 chr22, first TRACE batch

```bash
export GU_BUILD=37
export GU_CHRS=22
export UKB_HAP_ROOT=/mnt/h/ukbGen/hap
export UKB_REF_FASTA=/path/to/human_g1k_v37.fasta
export GU_SAMPLE_PANEL=$GU_ANALYSIS_ROOT/ukb/ukb.sample_panel.tsv

./gu.sh ukb inspect-hap
./gu.sh ukb make-panel
GU_ANCHORS_PER_GROUP=1000 GU_BATCH_SIZE=1000 ./gu.sh ukb batches

B=ALL.b0001
export UKB_KEEP=$GU_ANALYSIS_ROOT/ukb/batches/${B}.joint.txt
export UKB_VCF_OUT=$GU_ANALYSIS_ROOT/ukb/work/${B}/vcf
export UKB_ARG_VCF_OUT=$GU_ANALYSIS_ROOT/ukb/work/${B}/vcf.aa
export UKB_1KG_VCF_DIR=/mnt/d/data.BIG/refGen/1kg/GRCH37/vcf
./gu.sh ukb hap-arg-vcf

export GU_TARGET=ukb
export GU_TARGET_ROOT=$GU_ANALYSIS_ROOT/ukb/work/${B}
export GU_TARGET_VCF_DIR=$GU_TARGET_ROOT/vcf.aa
export GU_ARG_VCF_DIR=$GU_TARGET_ROOT/vcf.4arg
export GU_ARG_DIR=$GU_TARGET_ROOT/arg
export GU_ARG_BACKEND=tsinfer
./gu.sh arg build

export TRACE_CHRS=22
export TRACE_TARGET_SAMPLES=$GU_ANALYSIS_ROOT/ukb/batches/${B}.targets.txt
export TRACE_OUT=$GU_ANALYSIS_ROOT/ukb/trace/${B}
./gu.sh trace run

./gu.sh normalize
./gu.sh shiny
```

---

# 13. Software/reference links

- TRACE: https://github.com/YulinZhang9806/trace
- IBDmix: https://github.com/PrincetonUniversity/IBDmix
- ArchaicSeeker3.1-mamba: https://github.com/Shuhua-Group/ArchaicSeeker3.0
- tsinfer: https://tskit.dev/tsinfer/
- tsdate: https://tskit.dev/tsdate/
- tskit: https://tskit.dev/
- SINGER: https://github.com/popgenmethods/SINGER
- BGEN/bgenix: https://enkre.net/cgi-bin/code/bgen/
- PLINK 2: https://www.cog-genomics.org/plink/2.0/
- UKB Field 22438: https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=22438

---

# 14. Scientific caveats

1. **ARG quality is part of TRACE calling accuracy.** A faster ARG is not automatically a
   better TRACE input.
2. **1KG ARG != universal ARG.** A UKB target must be represented/threaded into the final
   genealogy used for its TRACE call.
3. **Field 22438 is phased but sparse relative to WGS.** This is a major limitation for
   tract boundary resolution and ghost-ancestry sensitivity.
4. **Batching is a computational strategy, not a statistical theorem.** Validate anchor
   and batch effects before treating all batch posteriors as exchangeable.
5. **Genome build consistency is mandatory.** The local Field 22438 files are GRCh37;
   do not mix them with GRCh38 archaic/reference resources without explicit conversion.
6. **AA transfer is conservative.** GU only transfers `INFO/AA` on exact CHROM/POS/REF/ALT
   matches and deliberately refuses strand guessing.
7. **AS3 has model-specific input assumptions.** Use the official resource/model build
   and validate marker-density effects before UKB-scale inference.
8. **Canonical GU segment codes should be versioned/frozen before final association.**
   Discovery-time clustering can evolve as more calls are added.
