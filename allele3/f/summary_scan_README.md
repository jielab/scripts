# Unified population summary scanner

`3.allele3_scan.sh` streams population SNV summaries into a reusable GRCh38 store. The data source is selected with `--target`; the public actions are `scan`, `merge` and `verify`. Available columns determine count/frequency capability, independently of the cohort name. The detailed gnomAD engine is selected with `--mode detailed`, using the same generic actions.

## Commands

From the `allele3` directory:

```bash
# Inputs default to data/annot/gnomad/38/sites/chr*.vcf.bgz.
# ALL means the aggregate population; lowercase all requests every available group.
./3.allele3_scan.sh scan --target gnomad --jobs 4 --populations ALL,AFR,NFE --replace FALSE --grch 38

# Input defaults to data/annot/allofus/38/hg38_allofus.txt.gz.
./3.allele3_scan.sh scan --target allofus --min-af 0.0001 --replace FALSE --grch 38

# Include all available ancestry groups for coexistence versus pooling analyses.
./3.allele3_scan.sh scan --target allofus --populations all --replace FALSE --grch 38

# Explicit AF thresholds are distinct from gnomAD's default AC thresholds.
./3.allele3_scan.sh scan --target gnomad --threshold-kind AF --min-af 0.0001 --grch 38 \
  --af-thresholds 0,0.00001,0.0001,0.001

# Inspect the full Python interface without starting a scan.
./3.allele3_scan.sh --help-full --grch 38

# Optional exhaustive gzip verification; ordinary verification checks receipts,
# file sizes and column headers, without rereading enormous source VCFs.
./3.allele3_scan.sh verify --target allofus --read-all --grch 38

# Explicit historical deep analysis.
./3.allele3_scan.sh scan --target gnomad --mode detailed --replace FALSE --grch 38
```

The scanner itself uses only Python's standard library. Scheduler options and process-tree resource limits belong to the shared shell launcher. The shell respects the repository-wide `--replace` convention; explicitly use `--replace FALSE` for resume. Direct Python invocation defaults to resume.

## What an observation means

| Input | Positive ALT evidence | Counts and singleton status | REF evidence |
|---|---|---|---|
| gnomAD VCF | Exact `FILTER=PASS` and filtered `AC>0` | Preserve `AC`, per-ALT `AN` and reported `AF`; no `AC_raw` fallback | Reconstruct only from a coherent joint record with complete nonnegative counts and a common `AN` |
| Supplied All of Us frequency excerpt | `gvs_{population}_af>0` | No AC/AN columns are present, so counts and singleton status remain unknown | Not established by independent AF summaries |
| All of Us with standard count columns | Actual `gvs_{population}_ac>0` | Automatically retain `gvs_{population}_ac` and `_an` when supplied; exact singleton requires nonmissing AC=1 | A count-bearing header alone does not declare joint REF counts |
| Custom TSV | Explicit AC column if provided; otherwise AF | Determined by the declared column schema | Only with an explicit `joint_counts: true` declaration and coherent complete counts |

All alternatives at one coordinate are aggregated, whether stored as one comma-separated VCF ALT field or separate rows. Arrays are aligned by actual nucleotide in alphabetical order. REF is a reference-sequence identity; it is not an ancestral state, the major allele, or the first allele to arise. ALT ordering is not frequency rank or temporal order.

Identical duplicates are deduplicated without summing counts. Compatible duplicate fields can fill a missing value and receive a flag. Conflicting summaries become unknown and receive a flag. REF conflicts, unequal per-ALT AN, filtered records and non-SNV alleles at the same anchor prevent inappropriate joint reconstruction. Split records with equal AN still do not prove that the same individuals were called for every allele.

Zero AF is reported zero and can reflect rounding; an empty tab field, `.`, or `NA` remains missing. The All of Us parser preserves empty tab-delimited fields. It reads the supplied `#Chr Start End Ref Alt ...` header or a TSV export with standard VAT `contig position ref_allele alt_allele` names, discovers standard `gvs_{population}_af`, `_ac` and `_an` columns from that header, and uses 1-based coordinates. This is a text-summary scanner, not a Parquet reader. Every `gvs_max_*` field is excluded from population-specific and pooled-frequency evidence. Explicit `--schema` mappings override the defaults for any target. With AF and AN but no AC, AN is preserved and AC is not reconstructed from rounded AF. The accompanying ordinary `.idx` is not assumed to be a Tabix index.

### ALT richness and strict occupancy

- `n_alt_observed`: number of non-reference bases with positive source evidence in this population.
- `n_alt_qualified`: number passing the explicitly stated primary AC or AF threshold.
- `observed_bases` / `observed_state`: full A/C/G/T occupancy only when joint REF evidence is available; otherwise `NA` / `unknown`.
- `qualified_bases` / `qualified_state`: full occupancy after the primary threshold only when all threshold inputs and REF are known.
- `singleton_status`: whether at least one observed ALT has exact AC=1 (`present`, `absent`, or `unknown`); this describes ALT singletons and does not claim that a minority REF singleton is an ALT singleton.

Thus two observed ALTs are valid multi-ALT evidence even when the strict total is unknown. A joint record with two ALTs and REF count zero is strictly biallelic, not triallelic. Missingness can make observed/qualified ALT counts lower bounds; flags and aligned arrays are retained for downstream exclusion or sensitivity analyses.

Sensitivity thresholds are accumulated in window and population summaries without repeating hundreds of millions of site records. AC thresholds default to 1,2,5,10; AF thresholds default to positive AF at 0,1e-5,1e-4,1e-3. AF=0 never becomes observed at the zero threshold. gnomAD's optional AF analysis still requires filtered AC evidence.

### The user's chr1:10001 example

The supplied `T>A`, `T>C` and `T>G` rows have pooled AF 0.013301, 0.003448 and
0.001478. All three ALT bases are positively observed in that population.
The data support a four-base **REF-plus-observed-ALT candidate**. The existing
strict `observed_state` may remain unknown because the separate AF rows do not
provide joint REF counts, but the new evidence columns retain what is known:

| Field | Value for pooled chr1:10001 | Meaning |
| --- | --- | --- |
| `observed_alt_bases` | `A,C,G` | Positive source evidence |
| `record_bases` / `record_allele_count` | `A,C,G,T` / 4 | REF plus every reported SNV identity, possibly including zeros at other sites |
| `ref_plus_observed_alt_bases` / `ref_plus_observed_alt_count` | `A,C,G,T` / 4 | Explicit reference-based candidate convention |
| `observed_alleles_min` / `observed_alleles_max` | 3 / 4 | Reported evidence guarantees three observed bases; REF is unresolved |
| `possible_observed_states` | `triallelic,quadallelic` | Total observed-base states compatible with the source evidence |
| `singleton_status` | `unknown` | Exact AC=1 is unavailable from these rounded AFs |

The per-population conclusions differ. AFR has A and C positive, with reported
G AF=0; EUR has A/C/G positive; AMR and EAS have only A positive among these
reported alternatives. This is usable population information, not a reason to
exclude the site. `gvs_max_af` is not an extra population or the pooled AF.

The other appended fields are `zero_alt_bases`, `unknown_alt_bases`,
`unreported_alt_bases`, `qualified_alt_bases`, `record_state` and
`ref_plus_observed_alt_state`. An unreported base is unmeasured unless a coherent
joint record already establishes exact occupancy. Bounds concern evidence in
the supplied summaries; reported AF=0 does not certify absence below the
source's rounding or detection limit. REF+positive-ALT classification retains
REF as an explicit convention even when an actual joint record establishes
REF count zero; the strict occupancy and bound columns show that difference.

### Population coexistence, pooling and frequency-rank changes

With `--populations all`, the scanner makes one additional compact classification
while each site's rows are already in memory. It records separate raw-positive
and primary-threshold patterns:

- `within_population_multi_ALT`: at least one reported ancestry group contains
  two or more positively observed/qualifying ALT bases.
- `between_population_only`: at least two ancestry groups have complete evidence
  for the source's reported ALT set, no group has more than one such ALT, but
  their union contains at least two. This describes the selected summaries,
  not the absence of lower-frequency alleles below detection.
- `unresolved`: missing evidence, missing requested groups, an insufficient
  number of groups, or non-comparable thresholds prevents the above conclusion.
- `not_multi_ALT`: the required ancestry evidence is complete and its union
  contains fewer than two such ALT bases.

`ALL` is a pooled label and never serves as another independent ancestry group.
Summary scans default to `--populations all` (all source strata). Use uppercase `ALL` to select only the pooled stratum.
An AC-versus-AF threshold mismatch across groups prevents a between-only
threshold conclusion. Missing ALT rows in one population are not filled with zero.

For ALT pairs with positive evidence and available AF in at least two ancestry
groups, the classification also records whether their AF ranking reverses.
Tied values are labelled separately. These are descriptive ranks and do not
establish mutation birth order or selection. A frequency ranking can change
between groups even though the nucleotide identities remain the same.

Only sites with at least two primary-threshold ALTs in the pooled summary or
ancestry union are exported to the compact candidate table. Classification
counts include **all** sites, including unresolved and non-candidate sites.

## Output contract

Default store: `analysis/allele3/allele3_scan/{gnomad|allofus|custom}`. Override with `--out`.

| Path relative to the store | Content |
|---|---|
| `shards/chrN.sites.tsv.gz` | One row per coordinate and requested population, primary threshold only |
| `shards/chrN.controls.tsv.gz` | Sparse deterministic sample of sites with one qualified ALT |
| `shards/chrN.windows.tsv.gz` | Per-chromosome window counts over sensitivity thresholds |
| `qc/chrN.json` | Atomic completion receipt, source path/size/mtime, parameters, output sizes and issues |
| `qc/manifest.json` | Completion state, build, input paths, configuration, source accounting and shard catalog |
| `results/window_summary.tsv.gz` | Merged windows, numerically sorted by chromosome |
| `results/chromosome_summary.tsv` | Counts by chromosome, population and threshold |
| `results/population_summary.tsv` | Counts by population and threshold |
| `results/occupancy_summary.tsv` | Occupancy and ALT-richness summary counts |
| `results/singleton_summary.tsv` | Supported ALT singleton counts and unknown-status counts |
| `results/site_catalog.tsv` | Reusable site/control shard paths and row counts |
| `results/allele_evidence_summary.tsv` | Positive ALT richness, reference-based candidates, occupancy bounds and actual AF/AC/AN availability by population |
| `results/source_capabilities.tsv` | Effective input column mappings and the analyses they can support |
| `shards/chrN.population_patterns.tsv.gz` | Primary-threshold candidate sites with coexistence/pooling patterns and AF-rank information |
| `results/population_patterns.tsv` | All-site classification counts by chromosome and combined |
| `results/population_candidates.tsv.gz` | Streamed merge of the complete compact candidate tables |
| `reports/index.html` | Self-contained HTML summary, SVG plot and interpretation notes |

The additional population-pattern and evidence files are described below.

Site positions are 1-based. Windows are 0-based half-open. A trailing window is nominally `window_size` bases; no chromosome-length or callable-base denominator is invented. All counts are descriptive among source records. They cannot by themselves estimate enrichment per callable base or establish balancing selection, allele birth dates, or mutation rates.

Base site columns (retained for existing stores):

```text
chrom pos ref population alts alt_afs alt_acs alt_ans
n_alt_observed n_alt_qualified ref_af ref_ac AN
observed_bases qualified_bases observed_state qualified_state singleton_status
count_mode threshold_kind threshold ref_evidence flags control primary_af
```

New scans append the evidence columns below. Existing complete base-schema shards remain usable: `iter_sites` derives the same fields from the stored arrays, without changing those shards or rescanning the source.

The physical table uses tabs. ALT arrays and base lists use commas. Numeric missingness is `NA`. `primary_af` is the maximum known ALT AF and is not the major-allele frequency.

### Downstream streaming API

```python
from summary_scan import iter_sites

for site in iter_sites("analysis/allele3/allele3_scan/gnomad",
                       chromosomes="1,2,X", population="ALL"):
    # Chromosomes are numeric ordered; positions ascend; populations are
    # grouped at each coordinate. No chromosome-sized list is constructed.
    if site["n_alt_qualified"] >= 2:
        alt_identities = site["alts"]       # e.g. ["C", "G"]
        reported_afs = site["alt_afs"]     # e.g. [0.2, None]
        counts = site["base_ac"]           # actual nucleotide -> int or None

for control in iter_sites(store, population="ALL", controls_only=True):
    pass  # Reads sparse control files directly.
```

`population=None` yields every stored group. `alt_afs`, `alt_acs`, and `alt_ans` are lists aligned to `alts`. Numeric missing fields become `None`. `base_af`, `base_ac`, and `base_an` map nucleotide identities, including a possibly unknown REF value. The reader does not derive missing AF from AC/AN, reorder alleles by frequency, or reinterpret REF as the major allele. It additionally returns the typed evidence fields below for both new and existing shards.

Controls use a deterministic BLAKE2 hash of seed, chromosome, position and REF. Selection is independent of functional scores, AVI, HERV and genomic annotations. Sites need exactly one qualified ALT and known threshold evidence; conflicting/invalid summaries are excluded. Hash sampling creates an auditable starting pool, not automatically matched controls; downstream analyses must still account for frequency, mutation context, AN/callability, coverage and annotation-specific ascertainment.

## Bounded memory, sharding and resume

Records must be sorted by position within one contiguous block per chromosome. The scanner holds one coordinate's ALT summaries and one window's counters; it never builds a Python row list for the whole input. Separate chromosome VCFs can run with `--jobs N`. A single ordinary gzip/BGZF TSV is read sequentially once, regardless of its `.idx` sidecar. Numeric chromosome order is imposed at merge even when input blocks were lexicographically ordered.

`--replace FALSE` reuses completed shards only when source path/size/mtime and scientific parameters match. Completion uses atomic per-shard receipts and matching output sizes. Interrupted temporary files are recomputed. A whole-store match avoids source rereading. A partially completed ordinary gzip TSV must still be decompressed sequentially to reach remaining chromosomes; a Tabix-capable source or externally split files can avoid that cost.

For independent scheduler jobs writing disjoint chromosomes to the same store, use `--no-merge`; do one merge after every job has finished:

```bash
./3.allele3_scan.sh scan --target gnomad --chr 1 --no-merge --replace FALSE --wait TRUE --grch 38
./3.allele3_scan.sh scan --target gnomad --chr 2 --no-merge --replace FALSE --wait TRUE --grch 38
./3.allele3_scan.sh merge --target gnomad --wait TRUE --grch 38
./3.allele3_scan.sh verify --target gnomad --wait TRUE --grch 38
```

Chromosome jobs may run concurrently when their output chromosomes are disjoint. Do not launch two writers for the same chromosome. `--no-merge` deliberately does not update shared reports or `qc/manifest.json`. Merge rejects incompatible scientific parameters. The manifest lists actual completed chromosomes; absent chromosomes are not silently interpreted as variant-free.

Earlier completed shards remain readable and reusable. The new site-level
evidence can be derived directly from their existing columns. New aggregate
or population-pattern statistics absent from earlier completion receipts are
reported as `not_available_from_prior_receipts` and NA, rather than invented or
silently set to zero. A reuse operation can refresh the small reports without
rescanning source variants. A fresh scan produces every new aggregate while
streaming once; no automatic whole-store migration or full-file hash is required.

## Custom TSV schema

Example for 1-based columns and AF-only data:

```json
{
  "format": "tsv",
  "build": "GRCh38",
  "coordinate_system": "1-based",
  "chrom": "chromosome",
  "pos": "position",
  "ref": "reference",
  "alt": "alternative",
  "populations": {
    "ALL": {"af": "frequency"},
    "AFR": {"af": "frequency_afr"}
  },
  "joint_counts": false
}
```

Invoke with `./3.allele3_scan.sh scan --target custom --input cohort.tsv.gz --schema cohort.json --grch 38`. Optional schema keys are `end`, `filter_column`, and `pass_values` (defaults to `["PASS"]`). Each population can declare `ac`, `an`, and/or `af` column names. AC and AF vectors must align with comma-separated ALT bases; a scalar AN may broadcast. Counts stay unknown when absent. `joint_counts: true` requires AC and AN columns and asserts that all rows at a coordinate partition one common called-allele denominator. Only make this declaration when the source actually guarantees that property. `0-based-half-open` is also supported when explicitly declared. GRCh37/hg19 must be lifted and validated upstream; this workflow does not silently convert builds.

## Validation

```bash
python3 -m unittest discover -s f/tests -p test_summary_scan.py -v
```

The focused evidence tests also use the user's exact chr1:10001 A/C/G example, automatic standard count-column discovery, AF+AN without invented AC, explicit schemas on All of Us, and source-free reuse of earlier shards.

The tests exercise joint and split counts, unequal AN, REF absence, three/four observed bases, exact singletons, filtered/raw counts, mixed indels, duplicate conflicts, missing/zero AF, ancestry restriction, AF versus AC thresholds, custom coordinate conversion, sortedness failures, sparse controls, atomic resume, independent chromosome merges and parallel VCF workers. They are synthetic edge-case fixtures; no full-scale user cohort has been rerun in the development environment.
