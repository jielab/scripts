#!/usr/bin/env python3
"""Species splitting with bounded transactional state and existing MetaCache DBs."""
from __future__ import annotations
import csv
import json
import os
from pathlib import Path
import signal
import shutil
import statistics
import subprocess
import sys
import time

from common import ROOT, DIR_RAW, DIR_SPLIT, FASTQ_RE, assert_reference_policy
from library import digest, file_state
from fastq_stream import prepare_reads, input_state
from sample_discovery import validate_dataset, validate_split_datasets
from sample_state import is_complete
from stream_cli import parser
from stream_runtime import (VERSION, sha256, durable_json, durable_unlink, locked, space_guard,
    safe_destination, library_snapshot, choose_backend, build_helper, helper_call, query_to_delta)


def valid_name(name):
    if not name or name in (".", "..") or "/" in name or "\\" in name or "\x00" in name:
        raise RuntimeError(f"Sample name must be one directory component: {name!r}")
    return name


def discover_pair_files(folder):
    return validate_dataset(folder)


def _selected_run(folder, stem):
    groups = {}
    for path in sorted(folder.rglob("*")):
        if not path.is_file(): continue
        m = FASTQ_RE.match(path.name)
        if not m or m.group(1) != stem: continue
        key = str(path.parent.resolve())
        mate = m.group(2) or m.group(3) or "single"
        row = groups.setdefault(key, {})
        if mate in row:
            raise RuntimeError(f"Ambiguous FASTQ copies: {row[mate]} and {path}")
        row[mate] = path
    if len(groups) > 1:
        raise RuntimeError(f"Sample {stem} occurs in several folders: {list(groups)}; use explicit --r1/--r2")
    if not groups: return None
    row = next(iter(groups.values()))
    if set(row) == {"1", "2"}: return stem, row["1"], row["2"]
    if set(row) == {"single"}: return stem, row["single"], None
    raise RuntimeError(f"Incomplete/mixed input for {stem}: {list(map(str,row.values()))}; use complete explicit inputs")


def discover(args):
    suffix = "__" + valid_name(args.output_suffix) if args.output_suffix else ""
    if args.command == "classify":
        name = valid_name(args.name)
        target = args.out_split or DIR_SPLIT / (name + suffix)
        result = [{"name": name+suffix, "dataset": "", "run": name, "r1": args.r1,
                   "r2": args.r2, "target": target}]
    elif args.sample:
        selector = args.sample.replace("\\", "/")
        selected_dataset = ""
        if "/" in selector:
            selected_dataset, stem = selector.split("/", 1)
        elif "__" in selector:
            selected_dataset, stem = selector.split("__", 1)
        else:
            stem = selector
        valid_name(stem)
        datasets = [x.strip() for x in args.dat_list.split(",") if x.strip()]
        if selected_dataset:
            valid_name(selected_dataset)
            if datasets and selected_dataset not in datasets:
                raise RuntimeError("--sample dataset is outside --dat-list")
            datasets = [selected_dataset]
        if not DIR_RAW.is_dir(): raise RuntimeError(f"Raw FASTQ root not found: {DIR_RAW}")
        datasets = datasets or sorted(p.name for p in DIR_RAW.iterdir() if p.is_dir() and not p.is_symlink())
        result = []
        for dataset in datasets:
            valid_name(dataset)
            folder = DIR_RAW / dataset
            if not folder.is_dir(): raise RuntimeError(f"Dataset not found: {folder}")
            found = _selected_run(folder, stem)
            if found:
                run, r1, r2 = found
                result.append({"name": f"{dataset}__{run}{suffix}", "dataset": dataset,
                    "run": run, "r1": r1, "r2": r2, "target": DIR_SPLIT/dataset/(run+suffix)})
        if len(result) != 1:
            raise RuntimeError(f"--sample {args.sample!r} matched {len(result)} samples; use dataset/run or explicit --r1/--r2")
    else:
        result = []
        for dataset, runs in validate_split_datasets(DIR_RAW, args.dat_list):
            valid_name(dataset)
            for run, r1, r2 in runs:
                result.append({"name": f"{dataset}__{run}{suffix}", "dataset": dataset,
                    "run": run, "r1": r1, "r2": r2, "target": DIR_SPLIT/dataset/(run+suffix)})
    for sample in result:
        valid_name(sample["name"])
        for key in ("r1", "r2"):
            if sample[key] is not None:
                path = Path(sample[key]).resolve(strict=True)
                if not path.is_file() or path.stat().st_size == 0:
                    raise RuntimeError(f"FASTQ missing/empty: {path}")
                if Path(str(path)+".aria2").exists():
                    raise RuntimeError(f"FASTQ still has an aria2 download marker: {path}.aria2")
                sample[key] = path
        if sample["r2"] is not None and sample["r1"].samefile(sample["r2"]):
            raise RuntimeError("R1 and R2 are the same file")
        sample["target"] = safe_destination(sample["target"])
    return result


def code_identity():
    files = ("stream_pipeline.py", "stream_runtime.py", "stream_cli.py", "fastq_stream.py",
             "stream_split.py", "compact_merge.cpp", "common.py", "split_metacache.py")
    return {name: sha256(ROOT/"f"/name) for name in files}


def context(sample, libraries, exe, backend, args):
    return {"version": VERSION, "name": sample["name"], "output": str(sample["target"]),
            "reads": [input_state(p) for p in (sample["r1"], sample["r2"]) if p],
            "libraries": libraries, "executable": {"path": str(exe), "sha256": sha256(exe), "backend": backend},
            "classifier": {"lowest": "species", "highest": "domain", "maxcand": 2, "hitmin": 5, "hitdiff": 100},
            "reference_evidence": "not_recorded", "query_scope": "all_completed_partitions",
            "sample_report": bool(args.screen), "query_cache": args.query_cache,
            "cache_gib": args.cache_gib, "code": code_identity()}


def recover_pending(work, identity, helper, progress):
    """Replay absolute replacements before querying another partition."""
    pending_path = work/"pending.json"
    if not pending_path.exists(): return progress
    pending = json.loads(pending_path.read_text())
    if pending.get("identity") != identity or pending.get("delta") != "pending.delta":
        raise RuntimeError(f"Incompatible pending transaction: {pending_path}")
    index = pending["index"]
    if index not in (len(progress["parts"]), len(progress["parts"])-1):
        raise RuntimeError("Pending transaction is outside the committed partition sequence")
    helper_call(helper, "apply", "--state", work/"candidates.bin", "--delta", work/"pending.delta")
    if index == len(progress["parts"]):
        progress["parts"].append(pending["result"])
        durable_json(work/"progress.json", progress)
    elif progress["parts"][index] != pending["result"]:
        raise RuntimeError("Pending and committed query records disagree")
    durable_unlink(pending_path)
    (work/"pending.delta").unlink(missing_ok=True)
    print(f"Recovered/committed partition {index+1}", flush=True)
    return progress


def clean_completed_work(work, identity):
    receipt = work/"run.json"
    if not receipt.is_file(): return
    record = json.loads(receipt.read_text())
    if record.get("identity") != identity or record.get("status") != "complete":
        raise RuntimeError("Refusing cleanup without a matching completed run receipt")
    reads = work/"reads"
    if reads.exists() and (reads.is_symlink() or
            (not record.get("cleanup_started") and not (reads/"prepared.json").is_file())):
        raise RuntimeError(f"Refusing cleanup of unowned prepared directory: {reads}")
    record["cleanup_started"] = True
    durable_json(receipt, record)
    for name in ("candidates.bin", "candidates.init.partial", "taxids.bin", "taxids.partial.bin", "pending.delta"):
        (work/name).unlink(missing_ok=True)
    if reads.exists():
        shutil.rmtree(reads)
    record["work_cleaned"] = True
    durable_json(receipt, record)


def run_sample(sample, libraries, exe, backend, helper, args):
    current = context(sample, libraries, exe, backend, args)
    identity = digest(current)
    target = sample["target"]
    work = safe_destination(args.work_root)/sample["name"]/identity[:20]
    if is_complete(target, identity):
        if not args.keep_work and (work/"run.json").is_file():
            with locked(work/"run.lock"):
                saved = json.loads((work/"run.json").read_text())
                if saved.get("identity") != identity:
                    raise RuntimeError("Completed output has a conflicting local run receipt")
                # A crash may occur after final output publication but before
                # the run receipt. The verified output receipt is authoritative.
                saved.update(status="complete", output=str(target))
                durable_json(work/"run.json", saved)
                clean_completed_work(work, identity)
        print(f"SKIP completed: {sample['name']} -> {target}", flush=True)
        return target
    if target.exists():
        raise RuntimeError(f"Existing output is incomplete or belongs to another version: {target}. "
                           "Use --output-suffix v2 (or an explicit new --out-split). Existing output is retained.")
    reserve = int(args.min_free_gib * 2**30)
    space_guard(work, reserve); space_guard(target, reserve)
    work.mkdir(parents=True, exist_ok=True)
    print(f"Sample: {sample['name']}\n  R1: {sample['r1']}\n  R2: {sample['r2']}\n  Work: {work}\n  Output: {target}", flush=True)
    with locked(work/"run.lock"), locked(target.parent/("."+target.name+".metag.lock")):
        if is_complete(target, identity): return target
        if target.exists(): raise RuntimeError(f"Output appeared while acquiring lock: {target}")
        receipt = work/"run.json"
        if receipt.exists():
            previous = json.loads(receipt.read_text())
            if previous.get("identity") != identity:
                raise RuntimeError(f"Work identity mismatch: {work}")
            if previous.get("status") == "complete" and (previous.get("work_cleaned") or previous.get("cleanup_started")):
                clean_completed_work(work, identity)
                for filename in ("progress.json", "export.json"):
                    (work/filename).unlink(missing_ok=True)
        durable_json(receipt, {"status": "running", "identity": identity, "context": current})
        prep = prepare_reads(sample["r1"], sample["r2"], work/"reads", cache_mode=args.query_cache,
                             min_free_bytes=reserve, max_cache_bytes=int(args.cache_gib*2**30))
        print(f"Validated {prep['fragments']:,} fragments; query cache: {prep['cache_mode']}", flush=True)
        query_states = {path: file_state(Path(path)) for path in
                        (prep["query_r1"], prep.get("query_r2"), prep["candidate_headers"]) if path}

        def verify_query_inputs():
            if any(file_state(Path(path)) != state for path, state in query_states.items()):
                raise RuntimeError("Prepared query input/header fingerprints changed during processing")
        state = work/"candidates.bin"
        progress_path = work/"progress.json"
        if progress_path.exists():
            progress = json.loads(progress_path.read_text())
            if progress.get("identity") != identity or progress.get("fragments") != prep["fragments"]:
                raise RuntimeError("Progress belongs to different inputs or libraries")
            if not state.is_file(): raise RuntimeError("Committed candidate state is missing")
        else:
            if (work/"pending.json").exists(): raise RuntimeError("Pending transaction has no base progress receipt")
            progress = {"identity": identity, "fragments": prep["fragments"], "parts": []}
            if not state.exists():
                space_guard(work, reserve, 64*prep["fragments"] + 4096)
                initial = work/"candidates.init.partial"
                initial.unlink(missing_ok=True)
                helper_call(helper, "init", "--state", initial, "--fragments", prep["fragments"])
                initial.replace(state)
            info = helper_call(helper, "info", "--state", state)
            if info.get("applied_parts") != 0:
                raise RuntimeError("Candidate state has no matching progress receipt")
            durable_json(progress_path, progress)
        progress = recover_pending(work, identity, helper, progress)
        info = helper_call(helper, "info", "--state", state)
        if info.get("dirty") or info.get("applied_parts") != len(progress["parts"]):
            raise RuntimeError("Candidate state/partition receipts disagree; retained for inspection")
        (work/"pending.delta").unlink(missing_ok=True)
        (work/"candidates.init.partial").unlink(missing_ok=True)
        parts = libraries["parts"]
        if len(progress["parts"]) > len(parts): raise RuntimeError("Too many committed partitions")
        for index, item in enumerate(progress["parts"]):
            if item.get("db") != parts[index]["db"]: raise RuntimeError("Committed partition order changed")
        logs = work/"logs"; logs.mkdir(exist_ok=True)
        for index in range(len(progress["parts"]), len(parts)):
            verify_query_inputs()
            row = parts[index]
            if any(file_state(Path(path)) != saved for path, saved in row["outputs"].items()):
                raise RuntimeError(f"Database changed during run: {row['db']}")
            space_guard(work, reserve, 32*prep["fragments"] + 4096)
            log = logs/f"part{index+1:04d}.log"
            print(f"[MetaCache {index+1}/{len(parts)}] {backend}: {Path(row['db']).name}; log={log}", flush=True)
            result = query_to_delta(exe, backend, helper, row["db"], prep, state,
                                    work/"pending.delta", log, args.batch_size, reserve)
            durable_json(work/"pending.json", {"identity": identity, "index": index,
                                               "delta": "pending.delta", "result": result})
            progress = recover_pending(work, identity, helper, progress)
            recent = [p["seconds"] for p in progress["parts"][-5:]]
            remaining = statistics.mean(recent)*(len(parts)-index-1)/3600
            print(f"  {result['seconds']:.1f}s; compact delta {result['delta_bytes']/2**20:.1f} MiB; "
                  f"rough remaining query time {remaining:.1f}h (recent partitions only)", flush=True)
        classifications = work/"taxids.bin"
        exported = work/"export.json"
        reuse = False
        if exported.exists():
            data = json.loads(exported.read_text())
            reuse = classifications.is_file() and data.get("identity") == identity and data.get("output") == file_state(classifications)
        if not reuse:
            temp = work/"taxids.partial.bin"; temp.unlink(missing_ok=True)
            space_guard(work, reserve, prep["fragments"]*8)
            summary = helper_call(helper, "export", "--state", state, "--taxonomy", libraries["taxonomy_path"],
                                  "--output", temp, "--hitmin", 5, "--hitdiff", 100)
            temp.replace(classifications)
            durable_json(exported, {"identity": identity, "output": file_state(classifications), "summary": summary})

        def verify_inputs():
            verify_query_inputs()
            if context(sample, library_snapshot(), exe, backend, args) != current:
                raise RuntimeError("Input reads, reference libraries, executable or code changed while processing")

        verify_inputs()
        gather = None
        if args.screen:
            from sourmash_screen import screen
            gather = screen(sample["r1"], sample["r2"], sample["name"], work/"sourmash")
        from stream_split import split
        print("Write species FASTQs from final full-reference classifications", flush=True)
        result = split(Path(prep["split_r1"]), Path(prep["split_r2"]) if prep.get("split_r2") else None,
            classifications, target, expected_fragments=prep["fragments"], input_identity=identity,
            before_publish=verify_inputs, gather=gather, expected_header_digest=prep["header_digest"],
            min_free_bytes=reserve, provenance={"identity": identity, "context": current,
                "fragments": prep["fragments"], "header_digest": prep["header_digest"],
                "query_cache_used": prep["cache_mode"], "partition_timing": progress["parts"]})
        durable_json(receipt, {"status": "complete", "identity": identity, "context": current,
                               "output": str(result), "fragments": prep["fragments"], "parts": len(parts)})
        if not args.keep_work:
            clean_completed_work(work, identity)
        return result


def benchmark(sample, libraries, exe, backend, helper, args):
    try:
        batches = list(dict.fromkeys(int(x) for x in args.benchmark_batches.split(",")))
    except ValueError as exc:
        raise RuntimeError("--benchmark-batches must be comma-separated positive integers") from exc
    if not batches or any(x < 512 for x in batches): raise RuntimeError("Benchmark batch sizes must be >=512")
    root = safe_destination(args.work_root)/"benchmarks"/(sample["name"]+"-"+time.strftime("%Y%m%d-%H%M%S")+f"-{os.getpid()}")
    root.mkdir(parents=True)
    reserve = int(args.min_free_gib*2**30)
    print(f"Benchmark auto-selected {sample['name']}\n  R1: {sample['r1']}\n  R2: {sample['r2']}\n  Results: {root}", flush=True)
    prep = prepare_reads(sample["r1"], sample["r2"], root/"reads", limit=args.benchmark_fragments,
                         cache_mode="none", min_free_bytes=reserve)
    all_parts = libraries["parts"]
    count = min(args.benchmark_parts, len(all_parts))
    indices = sorted({round(i*(len(all_parts)-1)/(count-1)) for i in range(count)}) if count > 1 else [0]
    results = []
    for batch in batches:
        for index in indices:
            folder = root/f"batch{batch}"/f"part{index+1:04d}"; folder.mkdir(parents=True)
            state = folder/"candidates.bin"
            helper_call(helper, "init", "--state", state, "--fragments", prep["fragments"])
            item = query_to_delta(exe, backend, helper, all_parts[index]["db"], prep, state,
                                  folder/"query.delta", folder/"query.log", batch, reserve)
            results.append({"batch_size": batch, "partition": index+1, "database": all_parts[index]["db"],
                            "fragments": prep["fragments"], "seconds": round(item["seconds"], 3),
                            "fragments_per_second": round(prep["fragments"]/item["seconds"], 2),
                            "delta_bytes": item["delta_bytes"]})
            durable_json(folder/"result.json", item)
            (folder/"query.delta").unlink(); state.unlink()
            print(f"Benchmark batch={batch}, part={index+1}: {item['seconds']:.2f}s", flush=True)
    with (root/"benchmark.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(results)
    durable_json(root/"benchmark.json", {"version": VERSION, "sample": sample["name"],
        "validated_scope": "prefix only; throughput diagnostic, not sensitivity validation",
        "requested_fragments": args.benchmark_fragments, "actual_fragments": prep["fragments"],
        "results": results, "includes_database_load": True, "species_output_created": False})
    print(f"PASS benchmark: {root/'benchmark.tsv'}", flush=True)
    return root


def main(argv=None):
    args = parser().parse_args(argv)
    assert_reference_policy()
    samples = discover(args)
    libraries = library_snapshot()
    safe_destination(args.work_root)
    print(f"MetaCache libraries: {len(libraries['parts'])} completed partitions; scope=full; reference evidence=not_recorded", flush=True)
    if args.plan:
        for sample in samples:
            print(f"{sample['name']}\n  R1: {sample['r1']}\n  R2: {sample['r2']}\n  Output: {sample['target']}")
        print(f"Samples: {len(samples)}; work root: {args.work_root}; query cache: {args.query_cache} <= {args.cache_gib:g} GiB")
        return 0
    exe, backend = choose_backend(args.backend)
    helper = build_helper()
    print(f"Backend: {backend} ({exe}); batch={args.batch_size}", flush=True)
    if args.benchmark:
        benchmark(samples[0], libraries, exe, backend, helper, args)
    else:
        for sample in samples:
            run_sample(sample, libraries, exe, backend, helper, args)
    return 0


def cli(argv=None):
    def stopped(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, stopped)
    try:
        return main(argv)
    except KeyboardInterrupt:
        print("Stopped. Committed partitions are retained; rerun the same command to resume.", file=sys.stderr)
        return 130
    except (RuntimeError, OSError, ValueError, KeyError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(cli())
