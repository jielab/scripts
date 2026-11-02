"""Matched-input, isolated classifier comparisons; agreement is not accuracy.

Every case runs in its own generated LSF process. Native build/conversion costs
are separate from these query/split timings. No command changes system caches.
"""
from contextlib import ExitStack, contextmanager
from collections import Counter, defaultdict
import csv
import gzip
import io
import itertools
import json
import math
import os
from pathlib import Path
import re
import resource
import socket
import sys
import tarfile
import time

from hpc_common import atomic_json, file_state, fingerprint, locked, sha256, utcnow
from ordinal_query import input_states

ENGINES = ("kraken", "kunpeng_full", "kunpeng_memory")
PARAMETERS = ("confidence", "minimum_hit_groups", "minimum_base_quality", "classify_threads")
PROFILE_PARAMETERS = PARAMETERS[:3]
SCHEMA = "metag-method-comparison-v1"
STATUSES = ("assigned_species", "higher_rank", "unclassified")


class ComparisonCaseError(RuntimeError):
    def __init__(self, receipt):
        self.receipt = receipt
        super().__init__(receipt["error"])


def _positive_integer(value, name, zero=False):
    if isinstance(value, bool) or not isinstance(value, int) or value < (0 if zero else 1):
        raise ValueError(f"{name} must be a {'nonnegative' if zero else 'positive'} integer")
    return value


def _profile_parameters(values):
    if not isinstance(values, dict) or set(values) - set(PROFILE_PARAMETERS):
        raise ValueError(f"Profile parameters must be a mapping using {PROFILE_PARAMETERS}; resources belong to the frozen HPC configuration")
    values = dict(values)
    for key, value in values.items():
        if key == "confidence":
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError("Profile confidence must lie in [0,1]")
        else:
            _positive_integer(value, key, zero=key == "minimum_base_quality")
    return values


def make_cases(samples, methods=None, repeats=2, limit=1_000_000, profiles=None):
    """Freeze matched prefixes without copying or mutating the original FASTQs.

    A profile is {name, classifier_engine, parameters}; parameters are declared
    before seeing results. A shared numerical confidence does not imply shared
    calibration across engines.
    """
    _positive_integer(repeats, "repeats")
    _positive_integer(limit, "limit", zero=True)
    samples = list(samples)
    if not samples or len({s["sample_id"] for s in samples}) != len(samples):
        raise ValueError("Comparison needs nonempty, distinct samples")
    if profiles is not None and methods is not None:
        raise ValueError("Choose methods or profiles, not both")
    if profiles is None:
        methods = ENGINES if methods is None else methods
        if isinstance(methods, str):
            methods = [m.strip() for m in methods.split(",") if m.strip()]
        profiles = [{"name": method, "classifier_engine": method, "parameters": {}} for method in methods]
    profiles = list(profiles)
    if not profiles:
        raise ValueError("At least one comparison profile is required")
    normalized, names = [], set()
    for profile in profiles:
        if not isinstance(profile, dict) or set(profile) - {"name", "classifier_engine", "parameters"}:
            raise ValueError("Profiles accept name, classifier_engine and parameters only")
        name = profile.get("name", "")
        engine = profile.get("classifier_engine")
        if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,79}", name) or name in names:
            raise ValueError("Profile names must be distinct safe path components")
        if engine not in ENGINES:
            raise ValueError(f"Unsupported comparison engine: {engine}")
        names.add(name)
        normalized.append((name, engine, _profile_parameters(profile.get("parameters", {}))))
    cases = []
    # Alternate profile order between repeats; cluster/cache state remains unknown.
    for sample in samples:
        for repetition in range(1, repeats + 1):
            order = normalized if repetition % 2 else list(reversed(normalized))
            for name, engine, parameters in order:
                context = {"schema": SCHEMA, "sample": dict(sample), "classifier_engine": engine,
                           "profile": name, "parameters": parameters, "repetition": repetition, "limit": limit}
                identity = fingerprint(context)
                cases.append({**context, "case_id": f"{name}.r{repetition:02d}.{identity[:16]}",
                              "case_fingerprint": identity})
    return cases


def _validate_case(case):
    expected = {k: v for k, v in case.items() if k not in ("case_id", "case_fingerprint")}
    if case.get("schema") != SCHEMA or fingerprint(expected) != case.get("case_fingerprint"):
        raise RuntimeError("Comparison case fingerprint is invalid")
    if case.get("classifier_engine") not in ENGINES:
        raise RuntimeError("Invalid comparison engine")
    _positive_integer(case["repetition"], "repetition")
    _positive_integer(case["limit"], "limit", zero=True)
    _profile_parameters(case["parameters"])
    if case["case_id"] != f"{case['profile']}.r{case['repetition']:02d}.{case['case_fingerprint'][:16]}":
        raise RuntimeError("Comparison case identifier differs from its contents")


def _usage():
    result = {}
    for label, who in (("self", resource.RUSAGE_SELF), ("children", resource.RUSAGE_CHILDREN)):
        usage = resource.getrusage(who)
        result[label] = {key: getattr(usage, key) for key in
                         ("ru_utime", "ru_stime", "ru_maxrss", "ru_inblock", "ru_oublock")}
    return result


def _resources(before, after):
    multiplier = 1024 if sys.platform.startswith("linux") else 1 if sys.platform == "darwin" else None
    result = {"measurement": "getrusage in one case process; peak values are not a process-tree sum",
              "job_peak_rss_bytes": None, "platform": sys.platform,
              "maxrss_native_unit": "KiB" if sys.platform.startswith("linux") else "bytes" if multiplier == 1 else "platform-defined"}
    for label in ("self", "children"):
        for key in ("ru_utime", "ru_stime", "ru_inblock", "ru_oublock"):
            result[f"{label}_{key}"] = after[label][key] - before[label][key]
        result[f"{label}_ru_maxrss"] = after[label]["ru_maxrss"]
    result["self_peak_rss_bytes"] = after["self"]["ru_maxrss"] * multiplier if multiplier else None
    result["max_child_peak_rss_bytes"] = after["children"]["ru_maxrss"] * multiplier if multiplier else None
    result["cpu_seconds"] = sum(result[f"{label}_{key}"] for label in ("self", "children") for key in ("ru_utime", "ru_stime"))
    return result


def _snapshot(cfg):
    if cfg.get("classifier_engine", "kraken") == "kraken":
        from hpc_reference import database_snapshot
    else:
        from hpc_kunpeng_index import database_snapshot
    return database_snapshot(cfg)


def run_case(cfg, case, case_dir, snapshot=None):
    """Run one native classifier and publish an isolated case.json receipt.

    Root orchestration must call this once per fresh .cmd process so getrusage
    high-water marks are not inherited from another benchmark case.
    """
    from hpc_query import run_sample, _complete
    _validate_case(case)
    case_dir = Path(case_dir).absolute()
    work_root = Path(cfg["work_dir"]).resolve()
    if case_dir.resolve() != case_dir or work_root not in case_dir.parents:
        raise RuntimeError("Comparison case must be an unaliased descendant of the configured work directory")
    local = dict(cfg)
    local.update(case["parameters"])
    local.update(classifier_engine=case["classifier_engine"], classify_workers=1,
                 output_dir=str(case_dir / "output"), work_dir=str(case_dir / "work"))
    if int(local.get("classify_threads", 8)) + 2 > int(local.get("classify_cores", 32)):
        raise RuntimeError("Comparison threads exceed the CPU allocation including parser overhead")
    snapshot = _snapshot(local) if snapshot is None else snapshot
    source_identity = snapshot.get("source_database_identity") or (snapshot["identity"] if case["classifier_engine"] == "kraken" else None)
    if not source_identity:
        raise RuntimeError("Kun-peng comparison index must record its exact source database identity")
    parameters = {key: local.get(key, {"confidence": 0.1, "minimum_hit_groups": 2,
                                      "minimum_base_quality": 0, "classify_threads": 8}[key]) for key in PARAMETERS}
    allocation = {"placement": local.get("compare_placement", "matched"),
                  "queue": local.get("classify_queue"), "cores": local.get("classify_cores"),
                  "threads": local.get("classify_threads", 8),
                  "memory_gib": local.get("classify_memory_gib"), "hours": local.get("classify_hours")}
    context = {"case": case, "inputs": input_states(case["sample"]),
               "database_identity": snapshot["identity"], "source_database_identity": source_identity,
               "parameters": parameters, "allocation": allocation, "output_format": local.get("output_format", "tar")}
    owner = {"schema": SCHEMA, "context": context, "identity": fingerprint(context)}
    lock = case_dir.parent / ".case-locks" / (case["case_id"] + ".lock")
    with locked(lock):
        marker = case_dir / ".comparison-owner.json"
        if case_dir.exists():
            if case_dir.is_symlink() or not marker.is_file() or json.loads(marker.read_text()) != owner:
                raise RuntimeError(f"Refusing unrelated comparison directory: {case_dir}")
        else:
            case_dir.mkdir(parents=True)
            atomic_json(marker, owner)
        receipt_path = case_dir / "case.json"
        previous = json.loads(receipt_path.read_text()) if receipt_path.is_file() else {}
        if previous.get("status") == "complete":
            if previous.get("identity") != owner["identity"]:
                raise RuntimeError("Comparison receipt identity changed")
            if not _complete(previous["output"], previous["output_identity"]):
                raise RuntimeError("Completed comparison output is missing")
            return {**previous, "reused_existing_case": True}
        record = {"schema": SCHEMA, "identity": owner["identity"], "case": case,
                  "status": "running", "attempt": previous.get("attempt", 0) + 1,
                  "started_utc": utcnow(), "hostname": socket.gethostname(),
                  "lsf_job_id": os.environ.get("LSB_JOBID"), "lsf_job_index": os.environ.get("LSB_JOBINDEX"),
                  "source_database_identity": source_identity, "database_identity": snapshot["identity"],
                  "index_bytes": sum(item.get("state", {}).get("size", 0) for item in snapshot.get("files", {}).values()),
                  "parameters": parameters, "allocation": allocation, "input_states": context["inputs"],
                  "cache_state": "unknown; no system cache manipulation; repeats may run on different hosts",
                  "confidence_note": "Numerical thresholds are recorded; calibration is not assumed equivalent between engines",
                  "build_cost_note": "Query/split measurement excludes index preparation, build and conversion"}
        atomic_json(receipt_path, record)
        before, started = _usage(), time.monotonic()
        try:
            result = run_sample(local, case["sample"], snapshot=snapshot, limit=case["limit"])
            elapsed = time.monotonic() - started
            summary = result["summary"]
            eligible = result["status"] == "complete" and not summary.get("query_reused", False)
            record.update(status="complete", completed_utc=utcnow(), output=result["output"],
                          output_identity=result["identity"], result_status=result["status"], summary=summary,
                          pipeline_wall_seconds=elapsed, resources=_resources(before, _usage()),
                          timing_eligible=eligible,
                          timing_exclusion="" if eligible else "query/output reused or recovered; this attempt is not a fresh end-to-end measurement")
            provenance = json.loads((Path(result["output"]) / "provenance.json").read_text())
            record["native"] = provenance.get("context", {}).get("native", {})
            query = provenance.get("query", {})
            record["query_stage_seconds"] = query.get("stage_seconds") or {
                item["stage"]: item["seconds"] for item in query.get("stages", [])}
            for key in ("input_staging_seconds", "normalized_input_bytes", "native_seconds", "peak_work_bytes_sampled"):
                record[key] = query.get(key)
            # A tiny process may exit between /proc samples. Zero is unobserved,
            # never evidence that a native classifier needed zero memory.
            record["native_peak_rss_bytes_sampled"] = query.get("native_peak_rss_bytes_sampled") or None
            atomic_json(receipt_path, record)
            return record
        except BaseException as exc:
            record.update(status="failed", completed_utc=utcnow(),
                          pipeline_wall_seconds=time.monotonic() - started,
                          resources=_resources(before, _usage()), timing_eligible=False,
                          error=f"{type(exc).__name__}: {exc}")
            atomic_json(receipt_path, record)
            raise ComparisonCaseError(record) from exc


@contextmanager
def _assignment_reader(output):
    output = Path(output)
    directory = output / "species/read-assignments.tsv.gz"
    archive = output / "species.tar"
    if directory.is_file() == archive.is_file():
        raise RuntimeError(f"Expected exactly one published assignment source: {output}")
    with ExitStack() as stack:
        if archive.is_file():
            tar = stack.enter_context(tarfile.open(archive, "r:"))
            entries = [m for m in tar.getmembers() if m.name == "read-assignments.tsv.gz"]
            if len(entries) != 1 or not entries[0].isfile():
                raise RuntimeError("Archive does not contain one regular read-assignments.tsv.gz")
            raw = stack.enter_context(tar.extractfile(entries[0]))
            zipped = stack.enter_context(gzip.GzipFile(fileobj=raw, mode="rb"))
        else:
            zipped = stack.enter_context(gzip.open(directory, "rb"))
        text = stack.enter_context(io.TextIOWrapper(zipped, encoding="utf-8", newline=""))
        reader = csv.DictReader(text, delimiter="\t")
        required = {"ordinal", "r1_read_id", "r2_read_id", "assigned_taxid", "assigned_rank", "species_taxid", "status"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames) or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise RuntimeError("Assignment table header is invalid")
        yield reader


def _uint(value, field):
    if not isinstance(value, str) or not re.fullmatch(r"[0-9]+", value):
        raise RuntimeError(f"Invalid {field}: {value!r}")
    return int(value)


def _assignments(record):
    expected = record["summary"]["fragments"]
    count = 0
    with _assignment_reader(record["output"]) as rows:
        for row in rows:
            count += 1
            ordinal = _uint(row["ordinal"], "ordinal")
            taxid = _uint(row["assigned_taxid"], "assigned_taxid")
            species = _uint(row["species_taxid"], "species_taxid")
            status = row["status"]
            if ordinal != count or count > expected:
                raise RuntimeError("Assignment ordinals are duplicated, missing, reordered or excessive")
            valid = ((status == "assigned_species" and species == taxid > 0 and row["assigned_rank"] == "species") or
                     (status == "higher_rank" and species == 0 and taxid > 0) or
                     (status == "unclassified" and species == taxid == 0))
            if not valid or not row["r1_read_id"] or bool(row["r2_read_id"]) != bool(record["summary"]["paired"]):
                raise RuntimeError("Assignment status/taxid/mate fields contradict one another")
            yield {"ordinal": ordinal, "taxid": taxid, "species": species, "status": status,
                   "r1": row["r1_read_id"], "r2": row["r2_read_id"]}
    if count != expected:
        raise RuntimeError(f"Assignment row count {count} differs from fragments {expected}")


def _check_pair_inputs(left, right):
    a, b = left["summary"], right["summary"]
    if (a["input_digest"], a["fragments"], a["paired"], left["source_database_identity"]) != (
            b["input_digest"], b["fragments"], b["paired"], right["source_database_identity"]):
        raise RuntimeError("Cannot compare different input bytes/order/count or source reference indexes")


def _count_pair(one, two, counts, status_pairs):
    if one is None or two is None or (one["ordinal"], one["r1"], one["r2"]) != (two["ordinal"], two["r1"], two["r2"]):
        raise RuntimeError("Methods did not publish matching ordinal read pairs")
    if one["species"] and two["species"]:
        counts["same_species" if one["species"] == two["species"] else "different_species"] += 1
    elif one["species"]:
        counts["left_species_right_unresolved"] += 1
    elif two["species"]:
        counts["right_species_left_unresolved"] += 1
    else:
        counts["both_unresolved"] += 1
    if (one["taxid"], one["status"]) == (two["taxid"], two["status"]):
        counts["same_normalized_assignment"] += 1
    status_pairs[one["status"] + ":" + two["status"]] += 1


def _pair_result(left, right, counts, status_pairs):
    a = left["summary"]
    both = counts["same_species"] + counts["different_species"]
    return {"sample_id": left["case"]["sample"]["sample_id"],
            "left_case": left["case"]["case_id"], "right_case": right["case"]["case_id"],
            "left_profile": left["case"]["profile"], "right_profile": right["case"]["profile"],
            "left_repeat": left["case"]["repetition"], "right_repeat": right["case"]["repetition"],
            "fragments": a["fragments"], **{key: counts[key] for key in
                ("same_species", "different_species", "left_species_right_unresolved",
                 "right_species_left_unresolved", "both_unresolved", "same_normalized_assignment")},
            "species_agreement_among_both_species": counts["same_species"] / both if both else None,
            "species_agreement_over_all_fragments": counts["same_species"] / a["fragments"],
            "status_pair_counts": json.dumps(dict(status_pairs), sort_keys=True),
            "interpretation": "agreement only; neither method is treated as ground truth"}


def _pair(left, right):
    """Two-case reference implementation, also useful for focused verification."""
    _check_pair_inputs(left, right)
    counts, status_pairs = Counter(), Counter()
    for one, two in itertools.zip_longest(_assignments(left), _assignments(right)):
        _count_pair(one, two, counts, status_pairs)
    return _pair_result(left, right, counts, status_pairs)


def _all_pairs(records):
    """Read each compatible case once while accumulating all pairwise counts.

    Memory holds one row per case and a fixed set of counters per pair, never
    all fragments. An invalid case/pair removes only affected comparisons;
    unrelated valid comparisons can still be reported as before.
    """
    pairs = list(itertools.combinations(range(len(records)), 2))
    agreements, errors, compatible = {}, {}, set()

    def fail(pair, exc):
        left, right = (records[index] for index in pair)
        errors.setdefault(pair, {"sample_id": left["case"]["sample"]["sample_id"],
            "left_case": left["case"]["case_id"], "right_case": right["case"]["case_id"],
            "error": f"{type(exc).__name__}: {exc}"})

    for pair in pairs:
        try:
            _check_pair_inputs(*(records[index] for index in pair))
            compatible.add(pair)
        except Exception as exc:
            fail(pair, exc)

    # Input equality forms cohorts. Cases without any compatible partner need
    # no second scan: summarize_cases() already validates each one separately.
    pending = {index for pair in compatible for index in pair}
    while pending:
        first = min(pending)
        members = [first] + sorted(index for index in pending - {first}
                                   if (first, index) in compatible)
        pending.difference_update(members)
        counters = {pair: (Counter(), Counter())
                    for pair in itertools.combinations(members, 2) if pair in compatible}
        with ExitStack() as stack:
            iterators = {}
            for index in members:
                iterator = iter(_assignments(records[index]))
                iterators[index] = iterator
                close = getattr(iterator, "close", None)
                if close is not None:
                    stack.callback(close)
            active, finished = set(members), set()
            while active and counters:
                rows = {}
                for index in sorted(active):
                    try:
                        rows[index] = next(iterators[index])
                    except StopIteration:
                        active.remove(index)
                        finished.add(index)
                    except Exception as exc:
                        active.remove(index)
                        for pair in list(counters):
                            if index in pair:
                                fail(pair, exc)
                                del counters[pair]
                for pair, (counts, status_pairs) in list(counters.items()):
                    left, right = pair
                    try:
                        if left in finished and right in finished:
                            agreements[pair] = _pair_result(records[left], records[right], counts, status_pairs)
                            del counters[pair]
                        else:
                            _count_pair(rows.get(left), rows.get(right), counts, status_pairs)
                    except Exception as exc:
                        fail(pair, exc)
                        del counters[pair]
    # Preserve the previous report ordering even when incompatible inputs cause
    # several independent cohorts to be scanned.
    return ([agreements[pair] for pair in pairs if pair in agreements],
            [errors[pair] for pair in pairs if pair in errors])


def _truth_rows(path, sample_id, expected, taxonomy=None):
    """Truth must cover exactly the tested prefix for each tested sample.

    species_taxid=0 explicitly means a known negative, not an unknown truth.
    reference_membership (optional) is in_reference, held_out or unknown.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    count = 0
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"sample_id", "ordinal", "species_taxid"}.issubset(reader.fieldnames):
            raise RuntimeError("Truth TSV needs sample_id, ordinal, species_taxid columns")
        for row in reader:
            if row["sample_id"] != sample_id:
                continue
            ordinal = _uint(row["ordinal"], "truth ordinal")
            species = _uint(row["species_taxid"], "truth species_taxid")
            if ordinal != count + 1 or ordinal > expected:
                raise RuntimeError("Truth rows must exactly cover tested ordinals in ascending order")
            membership = row.get("reference_membership") or "unknown"
            if membership not in ("in_reference", "held_out", "unknown"):
                raise RuntimeError("Invalid reference_membership in truth TSV")
            if taxonomy is not None and species and taxonomy.describe(species)["species_taxid"] != species:
                raise RuntimeError("Truth taxid must name a species in the frozen taxonomy")
            count += 1
            yield species, membership
    if count != expected:
        raise RuntimeError(f"Truth covers {count}/{expected} tested fragments for {sample_id}")


def _truth(record, path, taxonomy=None):
    counts = defaultdict(Counter)
    sample_id = record["case"]["sample"]["sample_id"]
    truth = _truth_rows(path, sample_id, record["summary"]["fragments"], taxonomy)
    for assigned, known in itertools.zip_longest(_assignments(record), truth):
        if assigned is None or known is None:
            raise RuntimeError("Truth and assignment counts differ")
        expected, membership = known
        for stratum in ("all", membership):
            c = counts[stratum]
            c["fragments"] += 1
            c["truth_species" if expected else "truth_negative"] += 1
            if assigned["species"]:
                c["assigned_species"] += 1
                c["correct_species" if expected and assigned["species"] == expected else "incorrect_species"] += 1
                if not expected:
                    c["negative_assigned_species"] += 1
            else:
                c[assigned["status"]] += 1
                if not expected and assigned["status"] == "unclassified":
                    c["negative_unclassified"] += 1
    result = []
    keys = ("fragments", "truth_species", "truth_negative", "assigned_species", "correct_species",
            "incorrect_species", "higher_rank", "unclassified", "negative_assigned_species", "negative_unclassified")
    for stratum in sorted(counts):
        c = counts[stratum]
        result.append({"case_id": record["case"]["case_id"], "sample_id": sample_id,
                       "profile": record["case"]["profile"], "reference_membership": stratum,
                       **{key: c[key] for key in keys},
                       "species_precision": c["correct_species"] / c["assigned_species"] if c["assigned_species"] else None,
                       "species_recall": c["correct_species"] / c["truth_species"] if c["truth_species"] else None})
    return result


def _write_tsv(path, rows, fallback):
    fields = list(dict.fromkeys(key for row in rows for key in row)) or list(fallback)
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize_cases(case_receipt_paths, out_dir, truth_path=None, taxonomy_dir=None):
    """Report all planned paths, including killed jobs with no final receipt.

    Invalid complete results and mismatched inputs make the report incomplete,
    with explicit errors. They are never silently treated as valid comparisons.
    """
    from hpc_query import _complete
    paths = [Path(path).absolute() for path in case_receipt_paths]
    if not paths or len(set(paths)) != len(paths):
        raise ValueError("Supply distinct receipt paths for every planned case")
    out = Path(out_dir).absolute()
    if out.resolve() != out:
        raise RuntimeError("Comparison report directory cannot be a symlink or alias")
    owner = {"schema": SCHEMA, "receipt_paths": list(map(str, paths)),
             "truth_path": str(Path(truth_path).absolute()) if truth_path else None}
    marker = out / ".report-owner.json"
    if out.exists():
        if not marker.is_file() or json.loads(marker.read_text()) != owner:
            raise RuntimeError("Refusing unrelated comparison report directory")
    else:
        out.mkdir(parents=True)
        atomic_json(marker, owner)
    taxonomy = None
    if taxonomy_dir:
        from taxonomy import Taxonomy
        taxonomy = Taxonomy(taxonomy_dir)
    truth_state = file_state(truth_path) if truth_path else None
    truth_hash = sha256(truth_path) if truth_path else None
    records, rows, errors = [], [], []
    for path in paths:
        row = {"case_id": path.parent.name, "receipt": str(path), "status": "missing"}
        try:
            if not path.is_file():
                rows.append(row)
                continue
            record = json.loads(path.read_text())
            if record.get("schema") != SCHEMA:
                raise RuntimeError("Invalid case receipt schema")
            case = record["case"]
            _validate_case(case)
            row.update(case_id=case["case_id"], status=record["status"],
                       sample_id=case["sample"]["sample_id"], profile=case["profile"],
                       classifier_engine=case["classifier_engine"], repetition=case["repetition"],
                       limit=case["limit"], hostname=record.get("hostname"), lsf_job_id=record.get("lsf_job_id"),
                       pipeline_wall_seconds=record.get("pipeline_wall_seconds"),
                       timing_eligible=record.get("timing_eligible", False),
                       error=record.get("error", ""), parameters=json.dumps(record.get("parameters", {}), sort_keys=True),
                       source_database_identity=record["source_database_identity"], index_bytes=record.get("index_bytes"))
            resources = record.get("resources", {})
            allocation = record.get("allocation", {})
            for key in ("placement", "queue", "cores", "threads", "memory_gib", "hours"):
                row["requested_" + key] = allocation.get(key)
            for key in ("self_peak_rss_bytes", "max_child_peak_rss_bytes", "job_peak_rss_bytes", "cpu_seconds"):
                row[key] = resources.get(key)
            for key in ("native_peak_rss_bytes_sampled", "input_staging_seconds", "normalized_input_bytes",
                        "native_seconds", "peak_work_bytes_sampled"):
                row[key] = record.get(key)
            row["query_stage_seconds"] = json.dumps(record.get("query_stage_seconds", {}), sort_keys=True)
            if record["status"] == "complete":
                if not _complete(record["output"], record["output_identity"]):
                    raise RuntimeError("Published output completion marker is missing")
                summary = record["summary"]
                if not summary.get("count_conservation") or not summary.get("input_digest") or summary.get("fragments", 0) < 1:
                    raise RuntimeError("Case did not validate positive fragment counts and input identity")
                # Validate even a lone method so malformed assignments cannot hide.
                observed = Counter(assignment["status"] for assignment in _assignments(record))
                expected_statuses = {key: int(summary.get("statuses", {}).get(key, 0)) for key in STATUSES}
                if any(observed[key] != expected_statuses[key] for key in STATUSES):
                    raise RuntimeError("Assignment statuses differ from the published summary")
                row.update(fragments=summary["fragments"], reads=summary["reads"],
                           input_digest=summary["input_digest"], output=record["output"],
                           classification_seconds=summary["classification_seconds"], split_seconds=summary["split_seconds"],
                           archive_seconds=summary.get("archive_seconds"), artifact_bytes=summary.get("artifact_bytes"),
                           **expected_statuses,
                           fragments_per_second=summary["fragments"] / record["pipeline_wall_seconds"]
                           if record.get("timing_eligible") and record.get("pipeline_wall_seconds", 0) > 0 else None)
                records.append(record)
        except Exception as exc:
            row.update(status="invalid", timing_eligible=False, error=f"{type(exc).__name__}: {exc}")
            errors.append({"receipt": str(path), "error": row["error"]})
        rows.append(row)
    grouped = defaultdict(list)
    for record in records:
        grouped[record["case"]["sample"]["sample_id"]].append(record)
    agreements, truth_rows = [], []
    for sample_id, group in grouped.items():
        comparisons, pair_errors = _all_pairs(group)
        agreements.extend(comparisons)
        errors.extend(pair_errors)
        if truth_path:
            for record in group:
                try:
                    truth_rows.extend(_truth(record, truth_path, taxonomy))
                except Exception as exc:
                    errors.append({"case_id": record["case"]["case_id"], "error": f"{type(exc).__name__}: {exc}"})
    if truth_path and file_state(truth_path) != truth_state:
        raise RuntimeError("Truth file changed while computing report")
    _write_tsv(out / "benchmark.tsv", rows, ("case_id", "status"))
    _write_tsv(out / "agreement.tsv", agreements, ("sample_id", "left_case", "right_case", "same_species"))
    if truth_path:
        _write_tsv(out / "truth.tsv", truth_rows, ("case_id", "reference_membership", "correct_species", "incorrect_species"))
    status_counts = dict(Counter(row["status"] for row in rows))
    report = {"schema": SCHEMA, "created_utc": utcnow(), "planned_cases": len(paths),
              "status": "complete" if not errors and len(records) == len(paths) else "incomplete",
              "case_statuses": status_counts, "valid_completed_cases": len(records),
              "timing_eligible_cases": sum(bool(record.get("timing_eligible")) for record in records),
              "pairwise_comparisons": len(agreements), "errors": errors,
              "truth": {"path": str(Path(truth_path).absolute()), "sha256": truth_hash,
                        "zero_means": "known negative, not unknown truth"} if truth_path else None,
              "interpretation": "Agreement between classifiers is not accuracy; species precision/recall require supplied truth",
              "timing_note": "Only timing_eligible cases represent fresh end-to-end query plus splitting; build/conversion costs excluded",
              "memory_note": "Self peak and maximum child peak are separate; aggregate LSF job peak is not measured here",
              "cache_note": "Cache state unknown; no cache flush; independent jobs can differ in host and shared-storage contention"}
    atomic_json(out / "report.json", report)
    return report
