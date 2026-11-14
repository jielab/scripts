"""Write English provenance/use notes beside references and validation outputs."""
from __future__ import annotations
import argparse
from pathlib import Path
import shlex

GRCH38_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
STEM = "GCA_000001405.15_GRCh38_no_alt_analysis_set"


def write_reference_readme(directory, url, reads_root, analysis_root, tools_root, code_root):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    source = ("NCBI GRCh38 no-alt analysis set, assembly accession `GCA_000001405.15`."
              if url == GRCH38_URL else "User-specified GRCH38_URL. The local filename does not establish the assembly version.")
    text = f"""# Human validation reference and BWA-MEM2 index

## Source and files

{source}

Reference download: [{url}]({url}). Reference directory: `{directory}`.
Pinned BWA-MEM2 2.2.1 is installed under `{tools_root}/bwa-mem2/2.2.1/`,
with `bwa-mem2-launch` as the executable. The bundled software archive is checked
at installation. `reference.url.txt` records the actual reference URL.

| File | Purpose |
|---|---|
| `{STEM}.fna.gz` | Downloaded reference FASTA; not sample reads |
| `{STEM}.fa` | Uncompressed alignment reference |
| `*.fa.0123`, `*.fa.amb`, `*.fa.ann`, `*.fa.bwt.2bit.64`, `*.fa.pac` | Locally built BWA-MEM2 index |
| `*.fa.building.*` | Incomplete, unpublished index files |
| `reference.sha256`, `bwa-mem2-index.sha256` | Explicit reference/index build checksums |
| `BUILD.COMPLETE` | Written only after successful index construction and publication |
| `index.jobid` | Last submitted LSF job ID; not a success marker |
| `human-db.env` | Reference and executable locations |

Other pre-existing reference files are not automatically substituted. Select the
exact FASTA and its matching sidecars. A large file alone does not establish a
complete index. Full checksum verification is an explicit, optional compute-node
operation, not a prerequisite for each sample split.

## Interpretation and input

`metag.sh validate --sample DATASET/SAMPLE` checks all saved human candidates.
`--scope full` screens all original reads and compares ordinal/read-ID matched
Kraken calls, including unclassified, above-species and other-species calls.
Candidate-only validation cannot estimate missed human-compatible reads.

Default thresholds are MAPQ >=20, aligned fraction >=0.80 and identity >=0.95.
The default pair rule is any_mate. Viral decoys `chrEBV` and `NC_007605.1` are
excluded. GRCh38-only alignment supports human compatibility; it does not exclude
other primates or shared sequences. Validation preserves original Kraken labels
and does not enable production host rescue automatically.

No separate sample dataset is downloaded. Existing candidate FASTQs are read
from split outputs; full scope uses originals under `{reads_root}` or explicit
`--target-fastq` paths. Keep inputs until their required analyses finish.

Results normally use `{analysis_root}/qc/human/DATASET/SAMPLE/candidates/` or
`full/`. `--replace TRUE` replaces an analysis by default; `--replace FALSE`
reuses a complete result and rebuilds incomplete results. Read
`human-validation.tsv` first; JSON retains detailed provenance and thresholds.

```bash
cd {shlex.quote(str(code_root))}
bash 1.build-lib.sh human --config config.json --submit hpc
bash metag.sh validate --sample DATASET/SAMPLE --submit hpc
```

On a local Linux/WSL machine, select `--submit local` and configure `--root`.
This small independent index does not replace or rebuild the five-category
Kraken library.
"""
    (directory / "README.md").write_text(text, encoding="utf-8")


def write_validation_readme(out, r1, r2, reference, assignments=None, summary=None, *, scope="full", sources=None):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    selected = sources if sources is not None else [{"r1": r1, "r2": r2}]
    rows = "\n".join(f"| {row.get('bin', 'explicit')} | `{row.get('r1')}` | `{row.get('r2') or 'single-end'}` |" for row in selected)
    interpretation = ("All saved human candidates were examined. The denominator is candidate fragments; this does not estimate missed reads or the total human fraction."
                      if scope == "candidates" else "All selected original fragments were examined. When matching assignments are supplied, agreement includes supported human calls, potential misses and other-species conflicts.")
    text = f"""# Independent human validation

Scope: `{scope}`. {interpretation}

## Inputs

| Bin/source | R1 | R2 |
|---|---|---|
{rows or '| No candidates | — | — |'}

Reference: `{reference}`. Assignments: `{assignments or 'not supplied'}`.
Split summary: `{summary or 'not supplied'}`. Candidate scope uses the saved
summary for candidate/sample counts; it does not compare a subset digest against
the digest of the complete sample. Reference details are in
`{Path(reference).parent / 'README.md'}`.

## Results

- `human-validation.tsv`: counts, denominators and fractions.
- `human-validation.json`: detailed scope, reference, thresholds and interpretation.
- `human-read-evidence.tsv`: candidate-mode per-fragment evidence with both read IDs,
  contig, position, CIGAR, MAPQ, identity, aligned fraction, flags and pass/fail.
- `human-validation-crosstab.tsv`: full-mode Kraken/alignment cross-tabulation.
- `human-host-summary.json`: alignment policy, digest and fragment counts.
- `human-host-calls.bin`: one byte per fragment; bits 0/1 identify passing mates.
- `human-host.log`: aligner log; the full SAM stream is not retained.
- `complete.json`: completion and cheap file-size inventory, written last.

Human-only alignment establishes human compatibility. The any_mate rule needs one
passing mate; both_mates needs both. Proper-pair support is reported separately.
The final JSON records the actual thresholds. No sample reads are downloaded,
backed up here or relabelled by validation.

`--replace TRUE` rebuilds this analysis by default. `--replace FALSE` reuses a
complete current result and rebuilds an incomplete result. A README or partial
intermediate file alone is not a completion marker.
"""
    (out / "README.md").write_text(text, encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", required=True)
    parser.add_argument("--url", default=GRCH38_URL)
    parser.add_argument("--reads-root", required=True)
    parser.add_argument("--analysis-root", required=True)
    parser.add_argument("--tools-root", required=True)
    parser.add_argument("--code-root", required=True)
    args = parser.parse_args()
    write_reference_readme(args.directory, args.url, args.reads_root, args.analysis_root, args.tools_root, args.code_root)
