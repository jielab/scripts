# Native minimizer-level union query

This extension is a CPU classifier for the independently built reference
categories. It opens all completed category compact hash tables with `mmap`,
reads the query FASTQ stream once, and applies this rule to each minimizer:

```
union_taxid = LCA(category_1_lookup, ..., category_N_lookup)
```

The values are **internal IDs in the same taxonomy**, not reported NCBI taxids.
The original pinned Kraken classifier receives the union values and performs its
read-pair scoring, confidence calculation, minimum-hit-group check, and reporting
once. Shared evidence contributes once, and two conflicting category hits can
resolve to their common ancestor before a species decision is made.

## Install and command contract

The Python API is `native_union.install.install(cfg)`. It needs `tools_dir` and
optionally `setup_cores`. It compiles from the checked, bundled source archive
without network access. The independent installation lives at:

```
<tools_dir>/metag_union/bin/metag_classify
<tools_dir>/metag_union/install.json
<tools_dir>/metag_union/install.log
```

An interrupted compilation is reusable. A published executable is verified
against its receipt. A different extension/source identity needs a different
tools directory; an existing installation is never silently replaced.

The native interface is:

```bash
metag_classify \
  --metag-hashes /absolute/path/hash-list.txt \
  -H /absolute/path/first-category/hash.k2d \
  -t /absolute/path/first-category/taxo.k2d \
  -o /absolute/path/first-category/opts.k2d \
  -M -p 32 -T 0.1 -g 2 -Q 0 -S
```

The normal pipeline invokes this through the ordinal FASTQ adapter; users need
not run it directly. `-S` selects paired interleaved input. With no filenames,
native input is stdin. Standard output retains the native five-column format.
Reports (`-R`, `-K`, etc.) still use the original classifier's combined evidence.
`-H` is retained as a required compatibility flag, but the supplied list defines
the actual tables. Memory mapping is enforced. Protein databases, daemon mode,
and quick early-stop mode are rejected.

The list contains one absolute `hash.k2d` path per line, without headers, comments,
or blank lines. Before the first read is classified, the extension checks:

- Every listed table exists, has a valid native header and exact expected length.
- Every category's sibling `taxo.k2d` and `opts.k2d` are byte-identical to `-t` and
  `-o`. The new builder initializes options padding to make this reproducible.
- Category cell/value widths match, and the shared taxonomy fits those widths.
- The list is nonempty and contains no duplicate resolved table paths.

The parent catalog records every expected category, including explicitly empty
categories. An empty category has no hash entry. Missing, failed, or unfinished
categories must never be converted into empty-category records or omitted from
the catalog silently. This low-level executable only validates its supplied list;
the parent catalog validation enforces completeness before constructing that list.

## Performance and limits

The code preserves the upstream batched lookup and software prefetch per table.
Each worker reuses a small scratch vector for a read pair. It does not create
five intermediate classification files or scan/decompress FASTQ five times.
However, a minimizer generally needs one lookup in **each** category, so a
five-category query is not guaranteed to match a single monolithic index's query
speed. Shared page cache, storage latency, index residency, memory bandwidth,
thread count, and category occupancy all matter. Benchmark representative input
before opening hundreds of simultaneous query jobs on shared storage.

Kraken's compact hash tables use truncated fingerprints. Independent tables can
have different false-positive collisions and probing behavior from a monolithic
table; their union is not promised to produce identical results for every possible
read. It correctly combines the evidence actually stored in the category tables.
The pinned builder also has subblock-boundary behavior that can produce small
differences if build parameters change. Use a single fixed build parameter set,
and keep the recorded source and index identities with results.

This extension does not use upstream `k2` multi-database output merging. Merging
already thresholded final read calls cannot reproduce confidence and minimizer
hit-group decisions from the original unthresholded combined evidence.

## Local verification

Run the synthetic native tests in an isolated local directory or a small compute
job; no production data is read:

```bash
METAG_RUN_NATIVE_TESTS=1 python -m unittest discover \
  -s f/tests -p 'test_union*.py' -v
```

The tests build two categories plus a monolithic control from deterministic
synthetic sequences and the same full sequence map. They verify species reads,
shared sequences, conflicting paired reads, unclassified reads, confidence and
minimum-hit-group filtering, order independence, and metadata/truncation failures.
A separate C++ oracle compares every batched result against scalar category
lookups followed by explicit taxonomy LCA over more than ten thousand keys.
The integration test also publishes a real two-category catalog and runs the
ordinal adapter plus species splitter over 1,001 paired fragments. It checks
each bin's exact original FASTQ bytes, including repeated IDs, CRLF, plus-line
comments and a missing terminal newline; then verifies successful reuse and
rejection of a valid same-size replacement input.

The MIT license and copyright for the underlying Kraken source remain in the
vendored archive and extracted source. `classify` is left unchanged; the extension
is a separately named `metag_classify` executable.
