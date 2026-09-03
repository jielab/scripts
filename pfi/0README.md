# Publication Fairness Index Pipeline (PFI)

This pipeline estimates whether a biomedical paper was published in a journal whose field-normalized prestige is higher or lower than what the paper's blinded scientific content would predict. The core analysis uses PubMed records published from 2000 onward, restricts the modeling population to papers with available PMC full text, samples 1,000,000 papers for prediction, and trains a blinded PubMedBERT/BiomedBERT multi-task model on a structured 1,000-paper review set.

The goal is not to label individual papers as "good" or "bad" in an absolute moral sense. The goal is to build a reproducible audit system for journal placement fairness after controlling for article content, article type, research domain, and publication year.

## Conceptual design

### Primary quantities

| Quantity | Meaning | Source |
|---|---|---|
| `actual_journal_score_0_100` | Field-normalized journal placement score | JCR/JCI-derived journal rank file |
| `predicted_paper_score` | Blinded content-derived paper score | PubMedBERT prediction from title/abstract/full text sections |
| `expected_journal_score_0_100` | Journal score expected from blinded content | PubMedBERT prediction trained against journal score in the review set |
| `fairness_index` | Actual minus expected journal score | `actual_journal_score_0_100 - expected_journal_score_0_100` |
| `fairness_z` | Domain/year/article-type standardized fairness residual | Robust z-score within strata |
| `fairness_category` | Five-level fairness classification | fair / modestly or strongly favored or unfavored |

### Five fairness categories

```text
fair                         abs(fairness_z) < 1
modestly_unfair_favored      1 <= fairness_z < 2
strongly_unfair_favored      fairness_z >= 2
modestly_unfair_unfavored    -2 < fairness_z <= -1
strongly_unfair_unfavored    fairness_z <= -2
```

`favored` means the actual journal placement is substantially higher than predicted from blinded content. `unfavored` means the actual journal placement is substantially lower than predicted from blinded content.

## Why full text only?

The default modeling and prediction population is restricted to `has_pmc_full_text=1`. PubMed-only and abstract-only records are kept in the master metadata table for coverage, but the default paper-score model uses full text so that methods, results, data availability, software, trial registration, and discussion claims can be evaluated.

## Why blinded text?

The model input must not contain journal names, authors, affiliations, funders, publishers, acknowledgments, or citation counts. Otherwise the model can simply learn journal prestige or author prestige instead of scientific content. The default blinded text contains:

```text
title
abstract
intro_text
methods_text
results_text
discussion_text
data_availability_text
software_availability_text
ethics_text
trial_registration_text
```

The review XLSX is also blinded. A separate key file stores PMID/PMCID/journal/author information for later audit but is not shown during scoring.

## Journal score construction

The preferred journal-rank input is a JCR-derived file with journal-year-category records. Because impact factors are not comparable across fields, the primary score should use category-normalized rank percentiles rather than raw impact factor. If a journal belongs to multiple JCR categories, the preferred score is:

1. match article domain to the closest JCR category/domain;
2. use the JIF percentile in that matched category;
3. if no matched category exists, use the mean percentile across categories;
4. keep `max_percentile` and `mean_percentile` for sensitivity analyses.

Expected input columns are documented in `conf/journal_rank_schema.csv`. Minimal required columns are:

```text
journal_key, jcr_year, jcr_category, journal_domain, jif_percentile, jif_quartile
```

Optional columns:

```text
journal_title, issn, eissn, impact_factor, journal_citation_indicator, total_cites
```

## Paper-score rubric

The structured review set uses 0-100 labels. The default dimensions are:

| Column | Definition |
|---|---|
| `new_results` | New empirical finding or important confirmation/refutation |
| `new_methods` | New method, design, assay, statistical framework, or analytical approach |
| `new_data` | New dataset, cohort, large-scale resource, trial, or deep phenotyping |
| `new_software` | Reusable software, algorithm, package, tool, database, or workflow |
| `clinical_public_health_value` | Potential clinical, translational, public-health, or policy relevance |
| `evidence_strength` | Whether the design and analysis support the claims |
| `reproducibility` | Data/code availability, transparent protocol, preregistration, validation |
| `claim_evidence_gap` | Penalty for overclaiming relative to evidence |
| `paper_score` | Overall paper-level score, 0-100 |

Suggested formula for internal consistency checking:

```text
paper_score = 0.18*new_results
            + 0.16*new_methods
            + 0.14*new_data
            + 0.12*new_software
            + 0.16*clinical_public_health_value
            + 0.14*evidence_strength
            + 0.10*reproducibility
            - 0.10*claim_evidence_gap
```

The formula is not forced; reviewers can edit `paper_score` manually after reading the blinded text. The formula-derived score is saved as `paper_score_formula` to identify inconsistent labels.

## Review-set design: 1,000 papers

A simple random sample of 1,000 papers is not enough because it will contain few landmark papers, few problematic papers, and too many middle cases. The default review set is therefore a stratified and anchored sample:

```text
400 field/year/journal-quartile balanced random papers
200 high-anchor papers
200 low/questionable-anchor papers
200 disagreement/gray-zone papers
```

### High anchors

High anchors are not automatically treated as perfect papers. They are used to make sure the review set contains high-end examples. Candidate sources:

- field/year top citation percentiles;
- papers later cited by guidelines or major consensus reports;
- major clinical trials, major methods papers, software/database papers;
- editor's choice/highlighted papers when the label is available;
- Nobel-associated papers or papers by Nobel laureates only as a calibration reference, not as proof of quality.

### Low/questionable anchors

Low anchors are also not automatically treated as bad science. They are used to improve coverage of the lower tail and problematic cases. Candidate sources:

- retracted papers or papers with expressions of concern;
- papers with known fabricated or unverifiable references;
- papers with major corrections;
- older papers with extremely low field/year-normalized citation uptake;
- papers with large claim-evidence-gap features.

Low citation alone must not be used as a low-quality label, especially for recent papers, niche fields, negative findings, replication studies, or local public-health studies.

### Expert-review practical solution

If independent experts are unavailable, use a defensible two-layer approach:

1. generate a blinded review file with model-free text features and anchor source flags;
2. conduct structured scoring with a fixed rubric;
3. score at least 150 overlapping papers twice to estimate reliability;
4. report intraclass correlation for continuous scores and weighted kappa for categorical bands;
5. treat citation/editorial/Nobel/retraction information as anchor metadata, not direct truth labels;
6. use sensitivity analyses excluding anchor papers and excluding papers with high reviewer disagreement.

This avoids the weak argument that the labels are arbitrary while still being feasible for a small research team.

## Folder layout

```text
/mnt/d/scripts/pfi                 code
/mnt/d/analysis/pfi/work           intermediate parquet files
/mnt/d/analysis/pfi/dat            generated CSV/parquet data tables
/mnt/d/analysis/pfi/ml             trained model files
/mnt/d/analysis/pfi/output         figures and reports
/mnt/d/analysis/pfi/log            logs
/mnt/g/pfi                         raw PubMed/PMC archives and large parquet archive
```

Shared paths and defaults are defined in:

```text
conf/pfi_paths.sh
```

## Environment

```bash
cd /mnt/d/scripts/pfi
./01_install_system_tools.sh
./01_create_conda_env.sh
conda activate pfi
```

The workflow entry point accepts explicit executable paths:

```bash
./pfi.sh --help
```

## Current defaults

The default run uses year 2000 onward, full-text articles, 1,000 review rows,
1,000 training rows, and a 1,000,000-row prediction sample. See
`./pfi.sh --help` for the corresponding command-line options.

## Main workflow

Use this order for a full run. Generated analysis files may be deleted and
rebuilt under `/mnt/d/analysis/pfi`; durable inputs such as
`journal_rank_real.csv` should live under `/mnt/d/scripts/pfi`.

### 0. Setup

```bash
cd /mnt/d/scripts/pfi
conda activate pfi

# Add --python and --rscript to any pfi.sh command when needed.
./pfi.sh --help
```

### 1. Build the base data

Run these once, or rerun from the first missing/failed step:

```bash
./pfi.sh download-pubmed
./pfi.sh parse-pubmed
./pfi.sh download-pmc
./pfi.sh parse-pmc
./pfi.sh build-master
./pfi.sh make-blinded
```

`07_make_blinded_text.sh` is resumable. If interrupted, rerun the same command.

### 2. Add the real journal-rank file

For scientific results, put the real JCR/JCI-derived file here:

```bash
/mnt/d/scripts/pfi/journal_rank_real.csv
```

Required columns are documented in `conf/journal_rank_schema.csv`. Minimal
columns:

```text
journal_key, jcr_year, jcr_category, journal_domain, jif_percentile, jif_quartile
```

For debugging only, a temporary proxy rank can be generated from
`articles_master` with `--allow-proxy-journal-rank TRUE`. Do not use proxy-rank
results for final interpretation.

### 3. Continue automatically after 06/07

This is the current shortest end-to-end path after `articles_master` and
`blinded_text` exist. It extracts features if needed, normalizes journal rank,
creates the review set, creates weak-supervision labels if no scored CSV exists,
trains the current text-model fallback, predicts, computes PFI, runs R figures,
and writes outlier-audit tables.

```bash
./pfi.sh after08 --journal-rank /mnt/d/scripts/pfi/journal_rank_real.csv --review-rows 1000 --max-train-rows 1000 --pred-sample-rows 1000000 --pred-limit 0
```

If you do not yet have real JCR/JCI data and only want to test the pipeline:

```bash
./pfi.sh after08 --review-rows 100 --max-train-rows 100 --pred-sample-rows 10000 --pred-limit 10000 --allow-proxy-journal-rank TRUE
```

### 4. Manual expert-label path

Use this path when you will fill the blinded review XLSX manually instead of
using weak-supervision pseudo labels.

```bash
./pfi.sh features
./pfi.sh normalize-rank --journal-rank /mnt/d/scripts/pfi/journal_rank_real.csv
./pfi.sh review-inputs
```

Fill:

```text
/mnt/d/analysis/pfi/dat/paper_score_review_1000_blinded.xlsx
```

Save the scored labels as:

```text
/mnt/d/analysis/pfi/dat/paper_score_review_1000_scored.csv
```

Then continue:

```bash
./pfi.sh prepare-training
./pfi.sh train
./pfi.sh predict
./pfi.sh compute
./pfi.sh analysis
./pfi.sh audit
```

### Useful controls

```bash
# Review sampling speed/coverage.
./pfi.sh review-inputs --review-max-scan-rows 200000

# Prediction size; 0 is no hard limit. Rebuild the sample when requested.
./pfi.sh predict --pred-sample-rows 1000000 --pred-limit 0 --pred-fresh-sample TRUE

# Optional LLM review batch after features exist.
./pfi.sh sample-llm --rows 1000
```

## Output contract

A completed run should produce:

```text
/mnt/d/analysis/pfi/dat/paper_score_review_1000_blinded.xlsx
/mnt/d/analysis/pfi/dat/paper_score_review_1000_scored.csv              # manual labels, if used
/mnt/d/analysis/pfi/dat/paper_score_review_1000_pseudo_scored.csv       # weak labels, if no manual labels
/mnt/d/analysis/pfi/dat/training_dataset.parquet
/mnt/d/analysis/pfi/ml/pubmedbert/
/mnt/d/analysis/pfi/dat/pubmedbert_prediction_fulltext_sample_1000000_seed1.parquet
/mnt/d/analysis/pfi/dat/pubmedbert_predictions.parquet
/mnt/d/analysis/pfi/dat/fairness_results.parquet
/mnt/d/analysis/pfi/dat/fairness_results.csv
/mnt/d/analysis/pfi/output/analysis/
```

## Interpretation cautions

1. This is a fairness audit, not a misconduct detector.
2. A strongly favored paper may be genuinely excellent in ways not captured by the model.
3. A strongly unfavored paper may be under-recognized, field-specific, negative, or local/public-health research.
4. Citation counts, Nobel links, editor choice labels, author prestige, and institution rank must be treated as audit variables or anchors, not as direct paper-quality truth.
5. The key validity check is whether the blinded model predicts paper-score labels and expected journal score with acceptable validation performance, and whether results are stable after excluding anchors and high-disagreement review cases.
