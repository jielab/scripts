"""Small, dependency-free utilities for the independent metag HPC workflow."""
from __future__ import annotations

from contextlib import contextmanager
from datetime import date, datetime, timezone
import fcntl
import hashlib
import json
import yaml
import math
import os
from pathlib import Path
import shutil
import tempfile

CODE_ROOT = Path(__file__).resolve().parent.parent
PATH_KEYS = ("reference_raw", "metadata_dir", "library_dir", "reads_raw",
             "output_dir", "work_dir", "tools_dir")


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
                                     separators=(",", ":")).encode()).hexdigest()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb", buffering=8 << 20) as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def file_state(path):
    value = Path(path).stat()
    return {"size": value.st_size}


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp = tempfile.mkstemp(prefix="." + path.name + ".", suffix=".partial",
                                dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True, ensure_ascii=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp, path)
        directory = os.open(path.parent, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        if os.path.exists(temp):
            os.unlink(temp)


@contextmanager
def locked(path, shared=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Locks remain in place; unlinking one would create different lock inodes.
    with path.open("a+b") as handle:
        fcntl.flock(handle, fcntl.LOCK_SH if shared else fcntl.LOCK_EX)
        try:
            yield handle
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def check_space(path, min_gib):
    target = Path(path)
    while not target.exists() and target != target.parent:
        target = target.parent
    free = shutil.disk_usage(target).free
    if free < float(min_gib) * 2**30:
        raise RuntimeError(f"Free space below reserve at {target}: "
                           f"{free / 2**30:.2f} GiB available, {min_gib} GiB required. "
                           "Filesystem free space does not include every project quota.")
    return free


def software_sources(root):
    path = Path(root) / 'software.yml'
    data = yaml.safe_load(path.read_text(encoding='utf-8'))
    if not isinstance(data, dict) or not {'kraken', 'kunpeng'} <= data.keys():
        raise ValueError(f'Missing software source records in {path}')
    return data


def read_config_file(path):
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    return yaml.safe_load(text) if path.suffix in (".yml", ".yaml") else json.loads(text)


def load_config(config_path=None, root_override=None):
    defaults = read_config_file(CODE_ROOT / "hpc.yml")
    source = Path(config_path or os.environ.get("METAG_HPC_CONFIG", CODE_ROOT / "hpc.yml"))
    values = read_config_file(source)
    if not isinstance(values, dict):
        raise ValueError("HPC configuration must be a mapping")
    cfg = dict(defaults)
    cfg.update(values)
    # Resolve only unfrozen defaults; generated config.json keeps numeric budgets.
    if cfg.get("convert_memory_gib") == "default":
        budget = os.environ.get("MEMORY_DEFAULT_GIB", "160")
        cfg["convert_memory_gib"] = int(budget)
    if isinstance(cfg["convert_memory_gib"], bool) or not isinstance(cfg["convert_memory_gib"], (int, float)) or not math.isfinite(cfg["convert_memory_gib"]) or cfg["convert_memory_gib"] <= 0:
        raise ValueError("convert_memory_gib must be a positive finite number")
    if cfg.get("schema") != 1:
        raise ValueError("Unsupported HPC configuration schema")
    if cfg.get("database_layout", "single") not in ("single", "partitioned"):
        raise ValueError("database_layout must be single or partitioned")
    root_value = root_override or str(cfg["root"]).replace("{dir_data}", os.environ.get("DIR_DATA", "/work/sph-huang"))
    root = Path(root_value).expanduser().absolute()
    cfg["root"] = str(root)
    for key in PATH_KEYS:
        scratch = os.environ.get("SCRATCH_ANALYSIS_ROOT", "/scratch/" + "2026-10-31" + "/sph-huangj/analysis")
        value = Path(str(cfg[key]).replace("{scratch_analysis}", scratch).replace("{scratch_root}", str(Path(scratch).parent))).expanduser()
        cfg[key] = str(value if value.is_absolute() else root / value)
    for key in ("legacy_library_dir", "legacy_metadata_dir", "parallel_reference_dir"):
        if cfg.get(key):
            value = Path(str(cfg[key])).expanduser()
            cfg[key] = str(value if value.is_absolute() else root / value)
    # These names deliberately avoid the previous library and split outputs.
    cfg.update(db_dir=str(Path(cfg["library_dir"]) / "db"),
               kunpeng_db_dir=str(Path(cfg["library_dir"]) / "kunpeng"),
               prepared_dir=str(Path(cfg["library_dir"]) / "prepared"),
               taxonomy_dir=str(Path(cfg["metadata_dir"]) / "taxonomy"),
               reference_plan=str(Path(cfg["library_dir"]) / "reference_plan.json"),
               code_root=str(CODE_ROOT), config_path=str(source.absolute()))
    if cfg.get("database_layout", "single") == "partitioned":
        cfg["db_dir"] = str(Path(cfg.get("parallel_reference_dir") or Path(cfg["library_dir"]) / "parallel") / "index")
        if cfg.get("classifier_engine", "kraken") != "kraken" or cfg.get("build_kunpeng", False):
            raise ValueError("Partitioned per-minimizer union supports engine kraken only. Kun-peng conversion requires a separate completed single index and database_layout: single.")
    if cfg.get("reference_taxid_map"):
        mapping = Path(cfg["reference_taxid_map"]).expanduser()
        cfg["reference_taxid_map"] = str(mapping if mapping.is_absolute() else root / mapping)
    if cfg.get("integrity_seal"):
        seal = Path(cfg["integrity_seal"]).expanduser()
        cfg["integrity_seal"] = str(seal if seal.is_absolute() else root / seal)
    for key in ("prepare_threads", "prepare_cores", "prepare_max_files", "prepare_concurrency",
                "setup_cores", "build_cores", "build_threads", "classify_cores",
                "classify_threads", "classify_workers", "classify_concurrency",
                "sample_batch_size", "pilot_fragments", "minimum_hit_groups", "max_array_size"):
        if isinstance(cfg[key], bool) or not isinstance(cfg[key], int) or cfg[key] < 1:
            raise ValueError(f"{key} must be a positive integer")
    for key in ("prepare_target_gib", "prepare_memory_gib", "setup_memory_gib",
                "build_memory_gib", "classify_memory_gib", "worker_memory_gib",
                "prepare_hours", "setup_hours", "build_hours", "classify_hours"):
        if isinstance(cfg[key], bool) or not isinstance(cfg[key], (int, float)) or not math.isfinite(cfg[key]) or cfg[key] <= 0:
            raise ValueError(f"{key} must be a positive finite number")
    for key in ("min_free_gib", "system_reserve_gib"):
        if not isinstance(cfg[key], (int, float)) or not math.isfinite(cfg[key]) or cfg[key] < 0:
            raise ValueError(f"{key} must be finite and nonnegative")
    if not isinstance(cfg["masking"], bool):
        raise ValueError("masking must be true or false, not a string")
    if not isinstance(cfg.get("build_kunpeng", True), bool):
        raise ValueError("build_kunpeng must be true or false, not a string")
    if not isinstance(cfg["minimum_base_quality"], int) or cfg["minimum_base_quality"] < 0:
        raise ValueError("minimum_base_quality must be a nonnegative integer")
    confidence = cfg["confidence"]
    if not isinstance(confidence, (int, float)) or not math.isfinite(confidence) or not 0 <= confidence <= 1:
        raise ValueError("confidence must lie in [0,1]")
    if cfg["output_format"] not in ("tar", "directory"):
        raise ValueError("output_format must be tar or directory")
    if cfg["memory_accounting"] not in ("conservative", "shared"):
        raise ValueError("memory_accounting must be conservative or shared")
    if cfg.get("classifier_engine", "kraken") not in ENGINES:
        raise ValueError("classifier_engine must be kraken, kunpeng_full or kunpeng_memory")
    if cfg.get("classifier_engine", "kraken") != "kraken" and cfg["classify_threads"] < 3:
        raise ValueError("Pinned Kun-peng requires at least 3 classify threads")
    if cfg.get("kunpeng_install_mode", "binary") not in ("binary", "source"):
        raise ValueError("kunpeng_install_mode must be binary or source")
    if cfg.get("kunpeng_input_compression", "gzip") not in ("gzip", "none"):
        raise ValueError("kunpeng_input_compression must be gzip or none")
    if cfg.get("compare_placement", "matched") not in ("matched", "production"):
        raise ValueError("compare_placement must be matched or production")
    for key in ("kunpeng_hash_capacity", "kunpeng_buffer_bytes", "kunpeng_batch_size",
                "compare_threads", "compare_cores", "compare_repeats", "compare_fragments",
                "compare_concurrency_full", "compare_concurrency_memory"):
        if isinstance(cfg[key], bool) or not isinstance(cfg[key], int) or cfg[key] < 1:
            raise ValueError(f"{key} must be a positive integer")
    for engine in ("kunpeng_full", "kunpeng_memory"):
        for suffix in ("cores", "threads", "workers", "concurrency"):
            key = f"{engine}_{suffix}"
            if isinstance(cfg[key], bool) or not isinstance(cfg[key], int) or cfg[key] < 1:
                raise ValueError(f"{key} must be a positive integer")
        if cfg[f"{engine}_workers"] * (cfg[f"{engine}_threads"] + 2) > cfg[f"{engine}_cores"]:
            raise ValueError(f"{engine} workers/threads exceed its requested CPU cores")
        if cfg[f"{engine}_threads"] < 3:
            raise ValueError(f"Pinned Kun-peng requires {engine}_threads >= 3")
        for suffix in ("memory_gib", "hours"):
            key = f"{engine}_{suffix}"
            if isinstance(cfg[key], bool) or not isinstance(cfg[key], (int, float)) or not math.isfinite(cfg[key]) or cfg[key] <= 0:
                raise ValueError(f"{key} must be a positive finite number")
    if cfg["compare_threads"] + 2 > cfg["compare_cores"]:
        raise ValueError("compare_cores must leave room for FASTQ processing")
    for key in ("kunpeng_max_work_gib", "kunpeng_memory_reserve_gib"):
        if isinstance(cfg[key], bool) or not isinstance(cfg[key], (int, float)) or not math.isfinite(cfg[key]) or cfg[key] < 0:
            raise ValueError(f"{key} must be finite and nonnegative")
    for phase in ("prepare", "build"):
        if cfg[phase + "_threads"] > cfg[phase + "_cores"]:
            raise ValueError(f"{phase}_threads exceeds its LSF CPU request")
    if cfg["classify_threads"] + 2 > cfg["classify_cores"]:
        raise ValueError("classify_cores must leave room for the FASTQ producer/consumer")
    if cfg["classify_workers"] * (cfg["classify_threads"] + 2) > cfg["classify_cores"]:
        raise ValueError("classify_workers × (classify_threads + 2) exceeds classify_cores")
    for key in (*PATH_KEYS, "root", "reference_taxid_map"):
        if cfg.get(key) and any(c in str(cfg[key]) for c in "\n\r\0"):
            raise ValueError(f"Control character in {key}")
    # Mutating paths must not be the raw inputs, nor ancestors of either input.
    raws = [Path(cfg[k]).resolve() for k in ("reference_raw", "reads_raw")]
    writes = [Path(cfg[k]).resolve() for k in ("metadata_dir", "library_dir", "output_dir", "work_dir", "tools_dir")]
    for dst in writes:
        for raw in raws:
            if dst == raw or dst in raw.parents or raw in dst.parents:
                raise ValueError(f"Output/cache path overlaps raw input: {dst} / {raw}")
    if len(set(writes)) != len(writes):
        raise ValueError("Metadata, library, output, work and tools paths must be distinct")
    return cfg


def code_identity():
    paths = [p for p in sorted((CODE_ROOT / "f").rglob("*")) if p.is_file()
             and p.suffix in (".py", ".sh", ".cc", ".h", ".patch") and "__pycache__" not in p.parts]
    paths += sorted(CODE_ROOT.glob("*.sh"))
    paths.append(CODE_ROOT / "software.yml")
    writer = CODE_ROOT / "f" / "native_checkpoint" / "checked_writer.txt"
    if writer.is_file():
        paths.append(writer)
    paths += sorted((CODE_ROOT / "f" / "vendor").glob("*source.tar.gz"))
    bundled = CODE_ROOT / "f" / "vendor" / "kun_peng-linux-x86_64"
    if bundled.is_file():
        paths.append(bundled)
    files = {str(p.relative_to(CODE_ROOT)): sha256(p) for p in paths}
    return fingerprint(files)


def require_compute(allow_local=False):
    if not os.environ.get("LSB_JOBID") and not allow_local and not (os.environ.get("submit") == "local" and os.environ.get("WORKFLOW_WORKER") == "1"):
        raise RuntimeError("This is a compute step. Submit its generated .cmd with bsub; "
                           "do not run it on the login node. --allow-local is only for a deliberate local test.")


ENGINES = ("kraken", "kunpeng_full", "kunpeng_memory")


def engine_snapshot(cfg):
    if cfg.get("classifier_engine", "kraken") == "kraken":
        if cfg.get("database_layout", "single") == "partitioned":
            from partitioned_build import database_snapshot
        else:
            from hpc_reference import database_snapshot
    else:
        if cfg.get("database_layout", "single") == "partitioned":
            raise ValueError("Kun-peng cannot consume the partitioned union catalog; use a separately completed single index for that benchmark")
        from hpc_kunpeng_index import database_snapshot
    return database_snapshot(cfg)


def engine_resources(cfg, engine=None, *, benchmark=False):
    """Select explicit execution resources; this never changes reference content."""
    result = dict(cfg)
    engine = engine or result.get("classifier_engine", "kraken")
    if engine not in ENGINES:
        raise ValueError(f"Unknown classifier engine: {engine}")
    result["classifier_engine"] = engine
    if engine != "kraken":
        for suffix in ("queue", "cores", "threads", "workers", "memory_gib", "hours", "concurrency"):
            result[f"classify_{suffix}"] = cfg[f"{engine}_{suffix}"]
        if engine == "kunpeng_memory":
            result["system_reserve_gib"] = cfg["kunpeng_memory_reserve_gib"]
    if benchmark:
        if engine != "kraken" and cfg["compare_threads"] < 3:
            raise ValueError("Kun-peng comparisons require compare_threads >= 3")
        result["classify_threads"] = cfg["compare_threads"]
        result["classify_cores"] = cfg["compare_cores"]
        result["classify_workers"] = 1
        if cfg.get("compare_placement", "matched") == "matched":
            # Identical node pool, requested CPUs and RAM for all three modes.
            result["classify_queue"] = cfg["classify_queue"]
            result["classify_memory_gib"] = cfg["classify_memory_gib"]
            result["classify_hours"] = cfg["classify_hours"]
            result["system_reserve_gib"] = cfg["system_reserve_gib"]
    return result


def effective_workers(cfg, snapshot):
    """Reserve enough physical memory and conservatively allow LSF RSS summation."""
    engine = cfg.get("classifier_engine", "kraken")
    if engine != "kraken":
        shards = [p for p in Path(snapshot["db_dir"]).glob("hash_*.k2d")
                  if p.name != "hash_config.k2d"]
        if not shards:
            raise RuntimeError("Kun-peng database has no hash shards")
        index_gib = sum(p.stat().st_size for p in shards) / 2**30
        budget = float(cfg["classify_memory_gib"]) - float(cfg["system_reserve_gib"])
        if engine == "kunpeng_full":
            per_worker = index_gib * 1.15 + float(cfg["worker_memory_gib"])
            note = "private full index per native process; no mmap sharing assumed"
        else:
            per_worker = (max(p.stat().st_size for p in shards) / 2**30 *
                          max(1, int(cfg["kunpeng_batch_size"])) + float(cfg["worker_memory_gib"]) +
                          int(cfg["kunpeng_buffer_bytes"]) * cfg["classify_threads"] / 2**30)
            note = "shard lower-bound estimate; sample-dependent native buffers still require measurement"
        capacity = min(math.floor(budget / per_worker),
                       cfg["classify_cores"] // (cfg["classify_threads"] + 2))
        if capacity < 1:
            raise RuntimeError(f"{engine} estimated worker floor {per_worker:.2f} GiB exceeds requested budget")
        result = dict(cfg)
        result["requested_classify_workers"] = cfg["classify_workers"]
        result["classify_workers"] = min(cfg["classify_workers"], capacity)
        print(f"{engine}: index={index_gib:.2f} GiB; workers={result['classify_workers']}; {note}", flush=True)
        return result
    names = list(snapshot["hashes"]) + ["taxo.k2d", "opts.k2d"] if snapshot.get("layout") == "partitioned" else ["hash.k2d", "taxo.k2d", "opts.k2d"]
    index_bytes = sum(Path(snapshot["db_dir"], name).stat().st_size for name in names)
    index_gib = index_bytes / 2**30
    budget = float(cfg["classify_memory_gib"]) - float(cfg["system_reserve_gib"])
    overhead = float(cfg["worker_memory_gib"])
    if cfg["memory_accounting"] == "shared":
        capacity = math.floor((budget - index_gib) / overhead)
    else:
        capacity = math.floor(budget / (index_gib + overhead))
    capacity = min(capacity, cfg["classify_cores"] // (cfg["classify_threads"] + 2))
    if capacity < 1:
        raise RuntimeError(f"Complete index {index_gib:.1f} GiB plus worker/reserve exceeds "
                           f"the {cfg['classify_memory_gib']} GiB classification budget. "
                           "Request a suitable node; this workflow never downsamples to fit.")
    effective = min(cfg["classify_workers"], capacity)
    print(f"Shared mmap index={index_gib:.2f} GiB; workers requested={cfg['classify_workers']}, "
          f"effective={effective}; scheduler accounting={cfg['memory_accounting']}", flush=True)
    result = dict(cfg)
    result["requested_classify_workers"] = cfg["classify_workers"]
    result["classify_workers"] = effective
    return result
