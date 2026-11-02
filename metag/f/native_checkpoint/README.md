# Resumable native index construction

This private extension builds a normal Kraken compact hash table. It retains the
**same in-memory table across prepared FASTA inputs**, then snapshots that table
at a complete-input boundary. Restart loads that saved mutable table and skips
the committed input prefix. It does not rebuild earlier references, merge five
independent labels, or treat a preparation receipt as a saved hash table.

`install(cfg)` in `install.py` independently extracts the checksum-verified,
pinned Kraken source and creates:

```
<tools_dir>/metag_checkpoint/bin/metag_build
<tools_dir>/metag_checkpoint/install.json
```

The original `kraken/bin/build_db`, masking executable, and official installation
receipts are unchanged. The installer changes only its private header copy:

1. `WriteTable` uses checked POSIX writes, handles short writes/EINTR, and checks
   both `fsync` and `close`. It creates a new file rather than truncating a saved one.
2. A private mutable table load initializes the OpenMP update locks. Upstream
   loading is intended for queries and does not initialize those locks.

The pinned deterministic minimizer scanning and LCA insertion algorithms are
unchanged. `-F`, protein mode, downsampling, and non-32-bit cells are forbidden.
`IndexOptions` is zero-initialized, including its padding, so separate category
builds can produce byte-identical `opts.k2d` files.

## Invocation contract

The standard native arguments are `-H HASH -t TAXONOMY -o OPTIONS -n NCBI_DIR
-m GLOBAL_SEQID_MAP -k K -l L -S SEED -p THREADS -B BLOCK -b SUBBLOCK -c CAPACITY`.
Use the **same full global sequence-to-taxid map** and taxonomy snapshot for every
category. Category-specific maps would change the internal taxon IDs and prevent
safe shared-taxonomy queries.

Additional arguments:

| Argument | Meaning |
|---|---|
| `--metag-inputs FILE` | UTF-8 TSV, no header: absolute FASTA path, byte size, lowercase SHA256, record count, base count |
| `--metag-context HEX` | Caller context, exactly 64 lowercase hex characters |
| `--metag-checkpoint-dir DIR` | Private persistent directory on the same filesystem as final hash |
| `--metag-checkpoint-seconds 14400` | Target interval between hash snapshots; evaluated at input boundaries |
| `--metag-max-seconds 396000` | Target soft-stop time for the invocation, including startup verification and hash loading |
| `--metag-stop-after-files N` | Pause after N additional input files in this invocation; test/controlled pause option |

Exit status **0** means every input is incorporated and the final hash was
published. Exit status **75** means a partial build has a valid durable checkpoint
and must be resubmitted. Every other nonzero status is an error. Do not convert
status 75 into successful build completion or allow downstream queries to begin.

`SIGUSR1`, `SIGUSR2`, or `SIGTERM` sets a soft-stop request. The current input must
finish before checkpointing. The wrapper must forward signals or use `exec` when
needed. A scheduler's eventual hard kill can arrive before a large checkpoint
finishes; signal handling does not create an unlimited grace period.

## Durability and storage

Each snapshot writes a new `stage-*` directory with a full `hash.k2d`, checked
native header/length, SHA256, machine receipt, and JSON receipt. It fsyncs payload
and metadata, renames the directory to immutable `cp-*`, then atomically publishes
and fsyncs `CURRENT`. Only after that commit does it remove older owned snapshots.
Incomplete stages are discarded on a verified restart. A kill during the first
generation's final publication gap can be recovered from that one fully checked
immutable generation even if `CURRENT` was not yet published.

`CURRENT` is one line:

```
cp-<timestamp>-<pid><TAB><SHA256 of receipt.tsv>
```

The canonical `receipt.tsv` binds the JSON receipt digest, table digest/size,
completed-input count, cumulative sequence/base counts, semantic parameters,
native executable, input manifest, mapping, and taxonomy identities. `receipt.json`
is a human/tool-friendly mirror. `progress.json` is updated atomically after each
file/checkpoint; `progress.csv` records progress over all invocations.

The latest valid snapshot is retained while a new snapshot is written. Plan for
**two complete hash payloads at peak**, plus prepared reference data and other
libraries. Filesystem free-space checks cannot see every project quota; checked
writes also detect quota exhaustion. There is deliberately no implicit cross-
filesystem multi-terabyte copy: final `-H` is published through an atomic hardlink,
so final DB and checkpoint directory must share a persistent filesystem. The
hardlink does not require a third payload copy.

SHA256 subprocesses are launched with `posix_spawnp`, safe argument arrays and
libc file actions; `POSIX_SPAWN_USEVFORK` is requested where GNU libc exposes it.
The builder never directly calls `fork()` or hand-written `vfork()`. This avoids
fork's address-space overcommit and page-table duplication problem while the
parent retains a multi-terabyte mutable hash table, including older Linux/glibc
environments where `posix_spawn` otherwise could choose `fork`.

## Resume identity and limits

Every restart checks the SHA256 of **all manifest FASTA inputs**, including the
completed prefix; then checks the saved hash before loading it. This costs an
additional reference/hash read on a restart but detects altered inputs and bad
snapshots. Inputs must remain immutable during the build. Changes in thread
count, block size, subblock size, k/l/seed, capacity, taxonomy, mapping, input
order/content, context, or native executable are rejected. Queue, walltime,
checkpoint interval, and soft-stop budget may change.

Checkpointing and soft-stop checks occur only after a whole prepared input file.
An interval of four hours is a target, not a hard maximum: processing one large
input and writing/checksumming a large table can exceed it. Bound prepared input
sizes and retain a generous walltime reserve. This extension cannot promise a
five-day total build or infer shared-storage throughput from a tiny fixture.
It makes completed checkpoint work reusable across separate jobs.

On hard interruption, work after the most recent committed snapshot may need
repeating. A SHA/identity mismatch fails closed; it never silently restarts at
zero. Existing uncheckpointed upstream `build_db` RAM state cannot be recovered
by installing this extension.

## Native verification

Compile with GCC >= 9, OpenMP and GNU coreutils `sha256sum`. Run the real native
tests with:

```bash
METAG_NATIVE_BUILD=/absolute/path/metag_build \
METAG_OFFICIAL_BUILD=/absolute/path/kraken/bin/build_db \
python3 -m unittest discover -s metag/f/tests -p 'test_checkpoint_native.py' -v
```

Tests compare uninterrupted/resumed and original deterministic hashes; deliver
real SIGKILL/SIGUSR2; force a partial hash write using `RLIMIT_FSIZE`; verify
preservation of the old CURRENT, category taxonomy/options equality, completed
input integrity, changed-parameter rejection, payload corruption rejection,
first-generation publication-gap recovery, and hardlink publication.
An additional real build/resume regression uses an `LD_PRELOAD` probe that makes
direct `fork`/`vfork` fail and verifies the GNU spawn flag; checkpoint hashing
must still complete with exactly the expected hash bytes.
