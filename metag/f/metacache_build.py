#!/usr/bin/env python3
from __future__ import annotations
import json, os, re, shutil, subprocess, time, fcntl
from pathlib import Path
from common import DIR_LIB, PART_COMPRESSED_GIB, atomic_json
from partition_refs import main as partition_refs
from local_metadata import TAXONOMY_FILES, prepare_taxonomy
from library import reference_identity, file_state
from build_limits import (run_limited, memory_gib, require_workflow_scope,
                          in_workflow_scope, run_workflow_limited)
from metacache_prepare import PreparationPool
from metacache_plan import plan_group

TOOLS = Path(os.environ.get("METAG_TOOLS", Path.home()/".cache/metag/tools"))
MC = TOOLS/"metacache"
CPU = MC/"metacache"
GPU = MC/"metacache_gpu"


def complete_part(db, identity):
    try:
        saved = json.loads(Path(str(db)+".complete.json").read_text())
        return (saved.get("identity") == identity and bool(saved.get("outputs")) and
                all(file_state(path) == state for path, state in saved["outputs"].items()))
    except (OSError, ValueError, KeyError):
        return False


def build_problem(db, log, returncode):
    if returncode:
        return f"exit status {returncode}"
    text = log.read_text(errors="replace")
    status = re.search(r"hashtable status[^\r\n]*(?:out of memory|not initialized|unknown error|probing length exceeded|invalid)[^\r\n]*", text, re.I)
    if status:
        detail = status.group(0).strip()
        if re.search(r"out of memory", detail, re.I):
            return (f"GPU VRAM allocation failed ({detail}); free GPU memory before resuming; "
                    "the builder sizes its initial hash table from free VRAM, so smaller input "
                    "partitions alone do not reduce this allocation")
        return f"GPU hash table reported an error ({detail})"
    if re.search(r"WARNING:|FAIL:|could not|not added to database|out of memory|capacity.*exceeded|targets remain unranked", text, re.I):
        return "upstream reported an incomplete build; zero exit status is insufficient"
    if any(float(value) >= 95 for value in re.findall(r"^features\s+\d+\s+\(([\d.]+)%\)", text, re.M)):
        return "GPU hash table reached its default load limit; input may be truncated"
    if "Total build time:" not in text or "Completed database writing." not in text:
        return "completion messages missing"
    outputs = [Path(str(db)+".meta"), *db.parent.glob(db.name+".cache*")]
    if len(outputs) < 2 or any(not p.is_file() or p.stat().st_size == 0 for p in outputs):
        return "database files missing/empty"
    return None


def archive_unfinished(db):
    paths = [p for p in db.parent.glob(db.name+".*")
             if p.suffix == ".meta" or ".cache" in p.name or
             p.name.endswith((".complete.json", ".build.log"))]
    if paths:
        folder = db.parent/"incomplete"/f"{db.name}-{time.time_ns()}"
        folder.mkdir(parents=True)
        for path in paths:
            path.rename(folder/path.name)

def backend():
    if GPU.is_file() and shutil.which("nvidia-smi"):
        return GPU, "gpu"
    if CPU.is_file():
        return CPU, "cpu"
    raise RuntimeError("MetaCache not installed; run ./install.sh")

def ensure_taxonomy():
    source = prepare_taxonomy()
    tax=DIR_LIB/"metacache"/"taxonomy"
    missing = [source/name for name in TAXONOMY_FILES if not (source/name).is_file()]
    if missing:
        raise RuntimeError("Missing local NCBI taxonomy:\n  " + "\n  ".join(map(str, missing)) +
                           "\nDownload taxdump.tar.gz as described in 0.README.md; extraction is automatic.")
    if source != tax:
        tax.mkdir(parents=True, exist_ok=True)
        for name in TAXONOMY_FILES:
            target = tax/name
            if target.is_file() and not target.is_symlink() and target.stat().st_size == (source/name).stat().st_size and target.stat().st_mtime_ns == (source/name).stat().st_mtime_ns:
                continue
            target.with_name(name + ".partial").unlink(missing_ok=True)
            shutil.copy2(source/name, target.with_name(name + ".partial"))
            target.with_name(name + ".partial").replace(target)
    print(f"NCBI taxonomy: local files in {source}", flush=True)
    return tax

def execute_job(job, tax, exe, memory_gb=16):
    require_workflow_scope(memory_gb)
    db=job["db"]
    dbroot=db.parent
    part_identity=job["identity"]
    used="reused"
    if complete_part(db, part_identity):
        print(f"SKIP existing {db}", flush=True)
    else:
        log=dbroot/(db.name+".build.log")
        archive_unfinished(db)
        print(f"[MetaCache build {db.name}] gpu: {len(job['inputs'])} input shards",flush=True)
        print(f"  Database: {db}\n  Build log: {log}",flush=True)
        command=[str(exe),"build",str(db),*[str(p) for p in job["inputs"]],
                 "-taxonomy",str(tax),"-kmerlen","16","-sketchlen","16",
                 "-winlen","127","-winstride","112","-parts","1"]
        with log.open("w") as lf:
            proc=run_limited(command,stdout=lf,stderr=subprocess.STDOUT,memory_gb=memory_gb)
        problem=build_problem(db,log,proc.returncode)
        if problem:
            raise RuntimeError(f"GPU MetaCache build incomplete: {db}: {problem}; see {log}. "
                               "Completed databases retained; CPU fallback is disabled.")
        outputs=[Path(str(db)+".meta"), *sorted(dbroot.glob(db.name+".cache*"))]
        atomic_json(Path(str(db)+".complete.json"), {
            "identity":part_identity,"coverage":job["coverage"],
            "outputs":{str(path):file_state(path) for path in outputs}})
        used="gpu"
    return used


def main(part_gib=PART_COMPRESSED_GIB, memory_gb=16, prepare_workers=2):
    require_workflow_scope(memory_gb)
    base = DIR_LIB/"metacache"
    base.mkdir(parents=True, exist_ok=True)
    with (base/"build.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("Another MetaCache workflow is running") from exc
        return build(part_gib, memory_gb, prepare_workers)


def build(part_gib=PART_COMPRESSED_GIB, memory_gb=16, prepare_workers=2):
    memory_gb = memory_gib(memory_gb)
    require_workflow_scope(memory_gb)
    exe, kind=backend()
    if kind != "gpu":
        raise RuntimeError("GPU builder required; CPU fallback is disabled")
    tax=ensure_taxonomy()
    if not (DIR_LIB/"reference"/"assembly_summary.txt").is_file():
        raise SystemExit("Reference mapping missing; run ./1.library.sh")

    base=DIR_LIB/"metacache"
    dbroot=base/"db"
    dbroot.mkdir(parents=True,exist_ok=True)
    receipt=dbroot/"build.json"
    identity = reference_identity() + f":part_gib={part_gib}"
    if receipt.is_file():
        saved=json.loads(receipt.read_text())
        if saved.get("status")=="complete" and saved.get("identity") == identity and all(
            complete_part(Path(row["db"]), row.get("identity", identity)) for row in saved.get("results", [])) and saved.get("results"):
            print(f"REUSE completed MetaCache DB: {receipt}")
            return 0
    atomic_json(receipt, {"status": "building", "identity": identity})

    # Custom conservative partitioning of the compressed reference files.
    # Upstream metacache-partition-genomes uses on-disk file size; with .fna.gz that
    # can underestimate final DB/VRAM substantially. We therefore use small 3-GiB
    # compressed partitions by default. Source /mnt/f/annot/metag is never modified.
    parts=partition_refs(part_gib)

    results=[]
    tasks=[]
    for i, part in enumerate(parts, 1):
        legacy=dbroot/f"part{i:04d}"
        if complete_part(legacy, identity):
            print(f"SKIP existing {legacy}", flush=True)
            results.append({"group":i,"first_shard":0,"input":str(part),"db":str(legacy),
                            "identity":identity,"backend":"reused"})
        else:
            tasks.append((i, part))

    print(f"MetaCache: {prepare_workers} input preparation processes; up to 4 x 512 MiB "
          f"inputs per GPU database; shared RAM cap {memory_gb} GiB", flush=True)
    with PreparationPool(tasks, base/"inputs", identity, workers=prepare_workers) as groups:
        for group, stream in groups:
            for job in plan_group(group, stream, dbroot, complete_part):
                db=job["db"]
                part_identity=job["identity"]
                used=execute_job(job, tax, exe, memory_gb)
                results.append({"group":group,"first_shard":job["coverage"][0]["shard"],
                                "inputs":[str(p) for p in job["inputs"]],"db":str(db),
                                "coverage":job["coverage"],"identity":part_identity,"backend":used})
    results.sort(key=lambda row:(row["group"],row["first_shard"]))
    for i,row in enumerate(results,1):row["part"]=i

    atomic_json(receipt, {
        "status":"complete","identity":identity,"partition_gib":part_gib,"parts":len(results),
        "taxonomy":str(tax),"results":results,
        "reference_policy":"KEEP source .fna.gz; derived bounded FASTA shards retained for resume",
    })
    print(f"PASS: {len(results)} MetaCache partitions")
    return 0

if __name__=="__main__":
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument("--part-gib",type=float,default=PART_COMPRESSED_GIB, help="compressed reference GiB per partition [3]")
    ap.add_argument("--memory-limit-gb", type=memory_gib, default=16)
    ap.add_argument("--prepare-workers",type=int,choices=range(1,5),default=2)
    a=ap.parse_args()
    import sys
    if not in_workflow_scope(a.memory_limit_gb):
        raise SystemExit(run_workflow_limited([sys.executable, __file__, *sys.argv[1:]],
                                             memory_gb=a.memory_limit_gb).returncode)
    raise SystemExit(main(a.part_gib, memory_gb=a.memory_limit_gb, prepare_workers=a.prepare_workers))
