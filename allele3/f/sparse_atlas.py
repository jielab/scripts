"""Bounded, optional AlphaGenome Atlas lookup and same-scorer geometry.

The default gnomAD scan never imports this module or calls an API. Planning
requires only the standard library. Actual lookup uses the verified Atlas
GetDenseVariantScores RPC (precomputed data, not model inference). Geometry
uses NumPy and cached raw effects; it does not make network requests.
"""
from __future__ import annotations

from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from contextlib import contextmanager
import csv
import dataclasses
from datetime import datetime, timezone
import enum
import fcntl
import gzip
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import re
import shutil
import tempfile
import threading
import time
from typing import Any


SCHEMA = "allele3-sparse-atlas-v1"
ENDPOINT = "dns:///gdmscience.googleapis.com:443"
VERIFIED_SDK_COMMIT = "55943112b83439e0c2512110f1ae607d4b5b492a"
DEFAULT_SCORERS = ("ATAC", "DNASE", "CHIP_TF", "CHIP_HISTONE", "CAGE", "PROCAP")
BASES = "ACGT"
REQUIRED_COLUMNS = ("site_id", "chrom", "pos", "ref", "anchor_ref", "anchor_alt",
                    "extra_alt", "s", "posterior")
SOURCE_CAVEAT = (
    "Atlas precomputed lookup; the service exposes no dataset-release pin. "
    "This cache records endpoint, SDK version, query settings, timestamps, "
    "fetch sessions and metadata hashes, not an immutable server release. "
    "Resume can span sessions; geometry requires matching complete axes and "
    "intervals. Matching metadata cannot prove unchanged model weights."
)
GEOMETRY_CAVEAT = (
    "Raw effects are compared separately for each scorer with a shared "
    "GRCh38 reference. A REF effect is the zero vector by definition. "
    "Distances are descriptive molecular predictions, not fitness effects. "
    "Response breadth is a fraction of predicted response dimensions above "
    "an absolute raw-effect threshold, not a count of biological functions "
    "or a validated measure of pleiotropy. The direction inferred from s "
    "remains an estimate, including when posterior confidence is low."
)


def _utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _hash(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     ensure_ascii=False, allow_nan=False).encode()).hexdigest()


def _file_hash(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def _json_write(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(data, handle, indent=2, ensure_ascii=False, allow_nan=False)
            handle.write("\n")
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _read_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _positive_int(value: int, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{label} must be a positive integer")
    return value


def _names(values, defaults, label) -> list[str]:
    values = defaults if values is None else values
    if isinstance(values, str):
        values = values.split(",")
    result = sorted({str(value).strip() for value in values if str(value).strip()})
    if not result and label == "scorers":
        raise ValueError("At least one scorer must be requested")
    pattern = r"[A-Za-z][A-Za-z0-9_]*" if label == "scorers" else r"(?:UBERON|CL|CLO|EFO):[0-9]+"
    if any(not re.fullmatch(pattern, name) for name in result):
        raise ValueError(f"Invalid {label}: {result}")
    return result


def _chrom(value: str) -> str:
    value = str(value).strip()
    value = value[3:] if value.lower().startswith("chr") else value
    if value.isdigit():
        value = str(int(value))
    else:
        value = value.upper()
    if value == "MT":
        value = "M"
    if value not in {*(str(i) for i in range(1, 23)), "X", "Y", "M"}:
        raise ValueError(f"Sparse Atlas requires a primary human chromosome, got {value!r}")
    return "chr" + value


def _variant(site: dict, alt: str) -> dict:
    return dict(chrom=site["chrom"], pos=site["pos"], ref=site["ref"], alt=alt)


def _variant_id(variant: dict) -> str:
    return f"{variant['chrom']}:{variant['pos']}:{variant['ref']}>{variant['alt']}"


def _variants(sites: list[dict]) -> list[dict]:
    unique = {}
    for site in sites:
        for alt in BASES:
            if alt != site["ref"]:
                value = _variant(site, alt)
                unique[_variant_id(value)] = value
    return [unique[key] for key in sorted(unique)]


def _site(row: dict, line: int) -> dict:
    if None in row or any(row.get(name) is None for name in REQUIRED_COLUMNS):
        raise ValueError(f"Malformed candidate row {line}")
    result = dict(row)
    result["site_id"] = str(row["site_id"]).strip()
    if not result["site_id"]:
        raise ValueError(f"Empty site_id at row {line}")
    result["chrom"] = _chrom(row["chrom"])
    try:
        result["pos"] = int(row["pos"])
        result["s"] = float(row["s"])
        result["posterior"] = float(row["posterior"])
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid position, s or posterior at row {line}") from exc
    if result["pos"] < 1 or not math.isfinite(result["s"]) or not math.isfinite(result["posterior"]):
        raise ValueError(f"Position must be 1-based and s/posterior finite at row {line}")
    if not 0 <= result["posterior"] <= 1:
        raise ValueError(f"posterior must lie in [0,1] at row {line}")
    for name in ("ref", "anchor_ref", "anchor_alt", "extra_alt"):
        result[name] = str(row[name]).strip().upper()
        if result[name] not in BASES or len(result[name]) != 1:
            raise ValueError(f"{name} must be one forward-strand A/C/G/T at row {line}")
    if len({result["anchor_ref"], result["anchor_alt"], result["extra_alt"]}) != 3:
        raise ValueError(f"anchor_ref, anchor_alt and extra_alt must be three distinct bases at row {line}")
    for name in ("build", "assembly", "genome_build"):
        if name in row and str(row[name]).strip().lower() not in {"grch38", "hg38"}:
            raise ValueError(f"Sparse Atlas only supports GRCh38; row {line} declares {row[name]!r}")
    if "strand" in row and str(row["strand"]).strip() not in {"+", "1", "forward"}:
        raise ValueError(f"Candidate row {line} must already be mapped to the GRCh38 forward strand")
    return result


def build_plan(candidate_tsv, out_dir, max_sites=2000, max_variants=6000,
               scorers=None, ontology_terms=None) -> dict:
    """Write an entirely offline bounded plan; never truncate an over-budget list.

    Each row denotes one site and supplies the ancient pair plus the observed
    extra base after mapping to GRCh38 forward-strand coordinates. All three
    non-REF substitutions are planned, even when one is not observed in gnomAD.
    Additional candidate columns are retained. No API key or SDK is needed.
    """
    _positive_int(max_sites, "max_sites")
    _positive_int(max_variants, "max_variants")
    candidate = Path(candidate_tsv).expanduser().resolve()
    scorers = _names(scorers, DEFAULT_SCORERS, "scorers")
    terms = _names(ontology_terms, (), "ontology_terms")
    sites, ids, loci = [], set(), set()
    opener = gzip.open if candidate.suffix in {".gz", ".bgz"} else open
    with opener(candidate, "rt", encoding="utf-8-sig", newline="") as handle:
        def lines():
            for line in handle:
                if line.startswith("##"):
                    if re.match(r"##\s*(build|assembly|genome_build)\s*[=:]", line, re.I):
                        declared = re.split(r"[=:]", line, maxsplit=1)[1].strip().lower()
                        if declared not in {"grch38", "hg38"}:
                            raise ValueError(f"Sparse Atlas only supports GRCh38, not {declared}")
                    continue
                if line.strip():
                    yield line[1:] if line.startswith("#site_id\t") else line
        reader = csv.DictReader(lines(), delimiter="\t")
        fields = list(reader.fieldnames or [])
        missing = set(REQUIRED_COLUMNS) - set(fields)
        if missing or len(set(fields)) != len(fields):
            raise ValueError(f"Candidate TSV needs unique column names and required columns; missing {sorted(missing)}")
        if any(name.startswith("atlas_") for name in fields):
            raise ValueError("Candidate columns starting atlas_ are reserved for Atlas outputs")
        for line, row in enumerate(reader, 2):
            site = _site(row, line)
            if site["site_id"] in ids:
                raise ValueError(f"Duplicate site_id: {site['site_id']}")
            locus = (site["chrom"], site["pos"])
            if locus in loci:
                raise ValueError(f"More than one candidate row at {locus}; supply one row per site")
            ids.add(site["site_id"])
            loci.add(locus)
            sites.append(site)
            if len(sites) > max_sites:
                raise ValueError(f"Candidate site budget exceeded: more than max_sites={max_sites}; no requests made")
            if 3 * len(sites) > max_variants:
                raise ValueError(f"Variant budget exceeded: three non-REF variants/site exceed max_variants={max_variants}; no requests made")
    variants = _variants(sites)
    body = dict(schema=SCHEMA, build="GRCh38", coordinate_system="1-based forward strand",
                candidate_columns=fields, candidates_sha256=_file_hash(candidate),
                sites=sites, variants=variants, scorers=scorers, ontology_terms=terms)
    output = Path(out_dir).expanduser().resolve()
    path = output / "plan.json"
    manifest = dict(body, plan_id=_hash(body), created_utc=_utc(),
                    candidate_tsv=str(candidate), plan_path=str(path),
                    n_sites=len(sites), n_variants=len(variants),
                    max_sites=max_sites, max_variants=max_variants,
                    network_requests_made=0, variant_rpc_count=len(variants),
                    implicit_reference="REF raw effect = zero; REF is not queried",
                    source_caveat=SOURCE_CAVEAT)
    if path.exists():
        existing = _load_plan(path)
        if existing["plan_id"] != manifest["plan_id"]:
            raise ValueError(f"Different plan already exists at {path}; choose a new output directory")
        return existing
    _json_write(path, manifest)
    return manifest


def _load_plan(path) -> dict:
    data = _read_json(Path(path))
    keys = ("schema", "build", "coordinate_system", "candidate_columns", "candidates_sha256",
            "sites", "variants", "scorers", "ontology_terms")
    if data.get("schema") != SCHEMA or data.get("build") != "GRCh38":
        raise ValueError("Unsupported sparse Atlas plan schema/build")
    if data.get("plan_id") != _hash({key: data[key] for key in keys}):
        raise ValueError("Sparse Atlas plan content or fingerprint changed")
    sites = [_site(row, i + 2) for i, row in enumerate(data["sites"])]
    if len({(row["chrom"], row["pos"]) for row in sites}) != len(sites):
        raise ValueError("Sparse Atlas plan contains duplicate loci")
    if len({row["site_id"] for row in sites}) != len(sites):
        raise ValueError("Sparse Atlas plan contains duplicate site IDs")
    if data["variants"] != _variants(sites) or data["n_sites"] != len(sites) or data["n_variants"] != len(data["variants"]):
        raise ValueError("Sparse Atlas variant list/count is inconsistent with candidate sites")
    if data["scorers"] != _names(data["scorers"], (), "scorers") or data["ontology_terms"] != _names(data["ontology_terms"], (), "ontology_terms"):
        raise ValueError("Sparse Atlas plan has invalid query settings")
    return data


def _settings(plan: dict) -> dict:
    return dict(schema=SCHEMA, build=plan["build"], scorers=plan["scorers"],
                ontology_terms=plan["ontology_terms"])


def _empty_catalog(plan: dict) -> dict:
    """An empty plan is a successful offline result, with no claimed SDK use."""
    settings = _settings(plan)
    source = dict(endpoint=ENDPOINT, sdk_version="not_loaded_empty_plan",
                  protocol="AtlasService.GetDenseVariantScores", lookup="precomputed",
                  verified_source_commit=VERIFIED_SDK_COMMIT, server_release="not_exposed")
    return dict(schema=SCHEMA, settings=settings, source=source,
                namespace=_hash(dict(settings=settings, source=source)), created_utc=_utc(),
                caveat=SOURCE_CAVEAT)


def _entry_path(cache: Path, catalog: dict, variant: dict) -> Path:
    return cache / _hash(dict(namespace=catalog["namespace"], variant=variant))


def _cache_entry(cache: Path, catalog: dict, variant: dict) -> tuple[Path, dict] | None:
    path = _entry_path(cache, catalog, variant)
    marker = path / "COMPLETE.json"
    if not marker.exists():
        if path.exists():
            raise ValueError(f"Unrecognized incomplete Atlas cache: {path}")
        return None
    complete = _read_json(marker)
    if complete.get("namespace") != catalog["namespace"] or complete.get("variant") != variant:
        raise ValueError(f"Atlas cache identity mismatch at {path}")
    for name in ("response.pb", "arrays.npz", "metadata.json"):
        file = path / name
        if not file.is_file() or _file_hash(file) != complete.get("files", {}).get(name):
            raise ValueError(f"Missing or changed completed Atlas cache file: {file}")
    return path, _read_json(path / "metadata.json")


def _json_value(value):
    """Preserve typed SDK metadata in JSON, without pickled Python objects."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else {"nonfinite_float": str(value)}
    if isinstance(value, enum.Enum):
        return dict(enum=type(value).__name__, name=value.name, value=_json_value(value.value))
    if dataclasses.is_dataclass(value):
        return {key: _json_value(item) for key, item in dataclasses.asdict(value).items()}
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if hasattr(value, "item"):
        return _json_value(value.item())
    # pandas.NA/NaT can occur in optional metadata columns.
    if type(value).__name__ in {"NAType", "NaTType"}:
        return None
    raise TypeError(f"Cannot losslessly export SDK metadata type {type(value).__name__}")


def _table(frame) -> dict:
    return dict(index=[_json_value(x) for x in frame.index.tolist()],
                index_names=[_json_value(x) for x in frame.index.names],
                columns=[str(x) for x in frame.columns],
                dtypes=[str(x) for x in frame.dtypes],
                data=[[_json_value(x) for x in row]
                      for row in frame.itertuples(index=False, name=None)])


def _response_packet(response, requested: dict) -> dict:
    """Use the official converter and retain the original protobuf verbatim."""
    import numpy as np
    from alphagenome.atlas import atlas
    from google.protobuf.json_format import MessageToDict

    actual = dict(chrom=_chrom(response.variant.chromosome), pos=response.variant.position,
                  ref=response.variant.reference_bases, alt=response.variant.alternate_bases)
    if actual != requested:
        raise ValueError(f"Atlas response variant mismatch: requested {requested}, received {actual}")
    raw_proto = response.SerializeToString()
    structure = MessageToDict(response, preserving_proto_field_name=True)
    for score in structure.get("scores", []):
        score.pop("scores", None)
        score.pop("calibrated_scores", None)
    by_name = {}
    for score in response.scores:
        name = score.variant_scorer.name
        if name in by_name:
            raise ValueError(f"Duplicate scorer in Atlas response: {name}")
        by_name[name] = score
    converted = atlas.convert_variant_scores_to_anndata([response], None)
    arrays, scorers = {}, []
    for i, (name, adata) in enumerate(sorted(converted.items())):
        raw_key = f"raw_{i}"
        arrays[raw_key] = np.asarray(adata.X).copy()
        item = dict(name=name, is_signed=bool(by_name[name].variant_scorer.is_signed),
                    raw_key=raw_key, shape=list(adata.X.shape), obs=_table(adata.obs),
                    var=_table(adata.var), uns=_json_value(dict(adata.uns)))
        if "quantiles" in adata.layers:
            key = f"quantiles_{i}"
            arrays[key] = np.asarray(adata.layers["quantiles"]).copy()
            item["quantiles_key"] = key
        scorers.append(item)
    return dict(protobuf=raw_proto, arrays=arrays, scorers=scorers,
                interval=structure.get("interval"), response_structure=structure,
                variant=actual)


class _FetchError(RuntimeError):
    def __init__(self, code: str, detail: str, *, retryable=False, fatal=False):
        self.code, self.retryable, self.fatal = code, retryable, fatal
        super().__init__(f"{code}: {detail}")


class _AtlasTransport:
    """Exactly one verified RPC per fetch; application code controls retries."""
    endpoint = ENDPOINT

    def __init__(self, api_key: str, timeout: float):
        try:
            import grpc
            from alphagenome.atlas import atlas_utils
            from alphagenome.data import genome
            from alphagenome.protos import atlas_service_pb2, atlas_service_pb2_grpc
            from alphagenome.protos.v1 import dna_model_pb2
        except ImportError as exc:
            raise RuntimeError("Optional Atlas query requires alphagenome==0.10.0; planning/scanning do not") from exc
        self.sdk_version = importlib.metadata.version("alphagenome")
        self.grpc, self.atlas_utils, self.genome = grpc, atlas_utils, genome
        self.messages, self.dna_messages = atlas_service_pb2, dna_model_pb2
        self.metadata = (("x-goog-api-key", api_key),)
        self.channel = grpc.secure_channel(
            self.endpoint, grpc.ssl_channel_credentials(),
            options=(("grpc.enable_retries", 0),),
        )
        try:
            grpc.channel_ready_future(self.channel).result(timeout=timeout)
        except BaseException:
            self.channel.close()
            raise
        self.stub = atlas_service_pb2_grpc.AtlasServiceStub(self.channel)

    def fetch(self, variant: dict, scorers: list[str], ontology_terms: list[str], timeout: float) -> dict:
        request = self.messages.GetDenseVariantScoresRequest(
            variant=self.genome.Variant(chromosome=variant["chrom"], position=variant["pos"],
                                        reference_bases=variant["ref"], alternate_bases=variant["alt"]).to_proto(),
            organism=self.dna_messages.ORGANISM_HOMO_SAPIENS,
            filter=self.atlas_utils.build_filter(requested_scorers=scorers,
                                                 ontology_terms=ontology_terms or None),
        )
        try:
            response = self.stub.GetDenseVariantScores(request, metadata=self.metadata, timeout=timeout)
        except self.grpc.RpcError as exc:
            code = exc.code().name
            # Never expose authentication material or server text in logs.
            raise _FetchError(code, "Atlas read failed",
                              retryable=code in {"UNAVAILABLE", "DEADLINE_EXCEEDED", "RESOURCE_EXHAUSTED", "ABORTED"},
                              fatal=code in {"UNAUTHENTICATED", "PERMISSION_DENIED", "UNIMPLEMENTED", "INVALID_ARGUMENT"}) from None
        return _response_packet(response, variant)

    def close(self):
        self.channel.close()


def _create_transport(api_key: str, timeout: float):
    return _AtlasTransport(api_key, timeout)


def _api_key(api_key_file, api_key_env) -> str:
    if api_key_file:
        value = Path(api_key_file).expanduser().read_text(encoding="utf-8").strip()
    else:
        value = os.environ.get(api_key_env, "").strip()
    if not value or any(character.isspace() for character in value):
        raise ValueError(f"A single API key is required in --api-key-file or environment {api_key_env}")
    return value


class _Limiter:
    def __init__(self, rate: float):
        self.interval, self.next_time = 1 / rate, 0.0
        self.lock = threading.Lock()

    def acquire(self, stop: threading.Event):
        with self.lock:
            now = time.monotonic()
            begin = max(now, self.next_time)
            self.next_time = begin + self.interval
        if stop.wait(max(0.0, begin - now)):
            raise _FetchError("CANCELLED", "Stopped before a new request", fatal=True)


def _retry_delay(attempt: int) -> float:
    return min(30.0, 2.0 ** min(attempt, 5))


def _save_packet(cache: Path, catalog: dict, variant: dict, packet: dict, session: str) -> dict:
    import numpy as np

    if packet.get("variant") != variant or not isinstance(packet.get("protobuf"), bytes):
        raise ValueError("Atlas response packet has missing raw protobuf or mismatched variant")
    if not packet["protobuf"]:
        raise ValueError("Atlas response protobuf is empty")
    arrays = packet["arrays"]
    names = set()
    for scorer in packet["scorers"]:
        if scorer["name"] in names or not isinstance(scorer.get("is_signed"), bool):
            raise ValueError("Duplicate scorer or missing signedness in Atlas packet")
        names.add(scorer["name"])
        raw = np.asarray(arrays[scorer["raw_key"]])
        if raw.ndim != 2 or raw.dtype.kind != "f" or list(raw.shape) != scorer["shape"]:
            raise ValueError("Atlas raw score shape/dtype is inconsistent")
        if raw.shape != (len(scorer["obs"]["data"]), len(scorer["var"]["data"])):
            raise ValueError("Atlas raw scores and labels have different shapes")
    path = _entry_path(cache, catalog, variant)
    partial = Path(tempfile.mkdtemp(prefix=".partial-", dir=cache))
    metadata = dict(schema=SCHEMA, variant=variant, scorers=packet["scorers"],
                    interval=packet.get("interval"), response_structure=packet.get("response_structure"),
                    namespace=catalog["namespace"], source=catalog["source"],
                    fetch_session=session, fetched_utc=_utc())
    metadata["axis_metadata_sha256"] = _hash(metadata["scorers"])
    try:
        (partial / "response.pb").write_bytes(packet["protobuf"])
        np.savez_compressed(partial / "arrays.npz", **arrays)
        _json_write(partial / "metadata.json", metadata)
        _json_write(partial / "COMPLETE.json", dict(namespace=catalog["namespace"], variant=variant,
                    files={name: _file_hash(partial / name) for name in ("response.pb", "arrays.npz", "metadata.json")}))
        os.replace(partial, path)
    finally:
        if partial.exists():
            shutil.rmtree(partial)
    return metadata


@contextmanager
def _run_lock(output: Path):
    output.mkdir(parents=True, exist_ok=True)
    with (output / ".atlas.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        yield


def run_sparse(plan_path, out_dir, api_key_file=None, api_key_env="ALPHAGENOME_API_KEY",
               max_variants=6000, workers=2, retries=3, timeout=60, *,
               requests_per_second=2.0, max_consecutive_failures=3) -> dict:
    """Explicit opt-in execution; cache completed variants and stop bounded failures.

    The hard budget is unique variants in the entire plan, including cached
    variants. At most retries+1 RPC attempts are made per missing variant in
    this invocation. Authentication/configuration errors stop immediately.
    Partial completion raises RuntimeError after writing query.json; resume
    with the same plan/out_dir reuses every verified complete cache entry.
    """
    _positive_int(max_variants, "max_variants")
    _positive_int(workers, "workers")
    _positive_int(max_consecutive_failures, "max_consecutive_failures")
    if isinstance(retries, bool) or not isinstance(retries, int) or retries < 0:
        raise ValueError("retries must be a nonnegative integer")
    if not math.isfinite(timeout) or timeout <= 0 or not math.isfinite(requests_per_second) or requests_per_second <= 0:
        raise ValueError("timeout and requests_per_second must be finite and positive")
    plan = _load_plan(plan_path)
    if plan["n_variants"] > max_variants:
        raise ValueError(f"Variant budget exceeded before connection: {plan['n_variants']} > {max_variants}")
    output = Path(out_dir).expanduser().resolve()
    with _run_lock(output):
        return _run_locked(plan, output, api_key_file, api_key_env, workers, retries, timeout,
                           requests_per_second, max_consecutive_failures)


def _run_locked(plan, output, api_key_file, api_key_env, workers, retries, timeout, rate, max_failures):
    cache = output / "cache"
    cache.mkdir(exist_ok=True)
    query_path, catalog_path = output / "query.json", cache / "catalog.json"
    if query_path.exists() and _read_json(query_path).get("plan_id") != plan["plan_id"]:
        raise ValueError("Sparse Atlas output belongs to another plan; use a separate output directory")
    settings = _settings(plan)
    catalog = _read_json(catalog_path) if catalog_path.exists() else None
    if catalog and catalog.get("settings") != settings:
        raise ValueError("Atlas cache request settings differ from this plan")
    if catalog and catalog.get("namespace") != _hash(dict(settings=settings, source=catalog["source"])):
        raise ValueError("Atlas cache namespace fingerprint changed")
    if catalog is None and any(cache.iterdir()):
        raise ValueError("Unrecognized existing Atlas cache; refusing to adopt it")
    complete, missing = {}, []
    for variant in plan["variants"]:
        entry = _cache_entry(cache, catalog, variant) if catalog else None
        if entry:
            complete[_variant_id(variant)] = entry[1]["fetch_session"]
        else:
            missing.append(variant)
    session = _hash(dict(plan=plan["plan_id"], time=_utc()))[:20]
    result = dict(schema=SCHEMA, plan_id=plan["plan_id"], plan_path=str(Path(plan["plan_path"]).resolve()),
                  cache_dir=str(cache), query_json=str(query_path), started_utc=_utc(),
                  fetch_session=session, n_sites=plan["n_sites"], n_variants=plan["n_variants"],
                  reused_variants=len(complete), requested_variants=len(missing),
                  max_variant_rpc_attempts=len(missing) * (retries + 1), rpc_attempts=0,
                  workers=workers, retries=retries, timeout=timeout, requests_per_second=rate,
                  max_consecutive_failures=max_failures, failures={}, source_caveat=SOURCE_CAVEAT)
    if not missing:
        if catalog is None:
            catalog = _empty_catalog(plan)
            _json_write(catalog_path, catalog)
        result.update(status="complete", completed_variants=len(complete), pending_variants=0,
                      finished_utc=_utc(), source=catalog["source"],
                      cache_fetch_sessions=sorted(set(complete.values())), empty_plan=not plan["sites"])
        _json_write(query_path, result)
        return result
    transport = _create_transport(_api_key(api_key_file, api_key_env), timeout)
    source = dict(endpoint=transport.endpoint, sdk_version=transport.sdk_version,
                  protocol="AtlasService.GetDenseVariantScores", lookup="precomputed",
                  verified_source_commit=VERIFIED_SDK_COMMIT, server_release="not_exposed")
    try:
        if catalog is not None and catalog["source"] != source:
            raise ValueError("SDK/endpoint/protocol changed since cached queries; use a separate output directory")
        if catalog is None:
            catalog = dict(schema=SCHEMA, settings=settings, source=source,
                           namespace=_hash(dict(settings=settings, source=source)), created_utc=_utc(),
                           caveat=SOURCE_CAVEAT)
            _json_write(catalog_path, catalog)
        result.update(status="running", source=source)
        _json_write(query_path, result)
        stop, limiter, counter_lock = threading.Event(), _Limiter(rate), threading.Lock()

        def fetch_one(variant):
            for attempt in range(retries + 1):
                limiter.acquire(stop)
                with counter_lock:
                    result["rpc_attempts"] += 1
                try:
                    packet = transport.fetch(variant, plan["scorers"], plan["ontology_terms"], timeout)
                    return _save_packet(cache, catalog, variant, packet, session)
                except _FetchError as exc:
                    if exc.fatal:
                        stop.set()
                    if not exc.retryable or attempt == retries or stop.wait(_retry_delay(attempt)):
                        raise
                except BaseException:
                    stop.set()
                    raise
            raise AssertionError("unreachable retry loop")

        todo = iter(missing)
        consecutive = 0
        interrupted = None
        with ThreadPoolExecutor(max_workers=workers) as pool:
            pending = {}
            def schedule_one():
                if stop.is_set():
                    return
                variant = next(todo, None)
                if variant is not None:
                    pending[pool.submit(fetch_one, variant)] = variant
            for _ in range(min(workers, len(missing))):
                schedule_one()
            try:
                while pending:
                    done, _ = wait(pending, return_when=FIRST_COMPLETED)
                    for future in done:
                        variant = pending.pop(future)
                        key = _variant_id(variant)
                        try:
                            metadata = future.result()
                            complete[key] = metadata["fetch_session"]
                            consecutive = 0
                        except _FetchError as exc:
                            if exc.code != "CANCELLED":
                                result["failures"][key] = dict(code=exc.code, retryable=exc.retryable)
                                consecutive += 1
                            if exc.fatal or consecutive >= max_failures:
                                stop.set()
                        except Exception as exc:
                            # Exception class only: a transport error string could contain a key.
                            result["failures"][key] = dict(code=type(exc).__name__, retryable=False)
                            stop.set()
                        schedule_one()
                    result.update(completed_variants=len(complete), pending_variants=plan["n_variants"] - len(complete))
                    _json_write(query_path, result)
            except BaseException as exc:
                stop.set()
                interrupted = exc
        # A Ctrl-C may race with a successful in-flight response; recover it from disk.
        complete = {}
        for variant in plan["variants"]:
            entry = _cache_entry(cache, catalog, variant)
            if entry:
                complete[_variant_id(variant)] = entry[1]["fetch_session"]
        result.update(completed_variants=len(complete), pending_variants=plan["n_variants"] - len(complete),
                      cache_fetch_sessions=sorted(set(complete.values())), finished_utc=_utc(),
                      status="complete" if len(complete) == plan["n_variants"] else "partial")
        _json_write(query_path, result)
        if interrupted is not None:
            raise interrupted
        if result["status"] != "complete":
            raise RuntimeError(f"Sparse Atlas query incomplete: {len(complete)}/{plan['n_variants']} variants cached; see {query_path}; resume the same plan")
        return result
    finally:
        transport.close()


def _axis_rows(table: dict, ignored=()) -> list[str]:
    columns = table["columns"]
    if len(set(columns)) != len(columns):
        raise ValueError("duplicate axis metadata columns")
    rows = []
    for row in table["data"]:
        if len(row) != len(columns):
            raise ValueError("axis metadata row length mismatch")
        labels = {key: value for key, value in zip(columns, row) if key not in ignored}
        rows.append(json.dumps(labels, sort_keys=True, separators=(",", ":"), ensure_ascii=False))
    if len(set(rows)) != len(rows):
        raise ValueError("ambiguous duplicate biological axis labels")
    return rows


def _vector(entry: tuple[Path, dict], scorer_name: str):
    import numpy as np

    path, metadata = entry
    matches = [item for item in metadata["scorers"] if item["name"] == scorer_name]
    if len(matches) != 1:
        return None, "missing_scorer", "Requested scorer absent from this cached response"
    item = matches[0]
    if not item["is_signed"]:
        return None, "unsigned_scorer", "Unsigned scores do not define a directional response vector"
    try:
        # obs variant/interval identify the input, not a gene/junction axis.
        rows = _axis_rows(item["obs"], ignored=("variant", "interval"))
        columns = _axis_rows(item["var"])
    except ValueError as exc:
        return None, "ambiguous_axes", str(exc)
    with np.load(path / "arrays.npz", allow_pickle=False) as arrays:
        raw = np.asarray(arrays[item["raw_key"]], dtype=np.float64)
    if raw.shape != (len(rows), len(columns)) or not raw.size:
        return None, "empty_or_invalid_shape", "Raw array has empty or incompatible axes"
    if not np.isfinite(raw).all():
        return None, "nonfinite_raw", "At least one raw effect is not finite"
    features = [(row, column) for row in rows for column in columns]
    order = sorted(range(len(features)), key=features.__getitem__)
    return dict(values=raw.reshape(-1)[order], axes=[features[i] for i in order],
                interval=metadata.get("interval"), session=metadata["fetch_session"]), "ok", ""


def _direction_metrics(extra, start, end) -> dict:
    import numpy as np

    delta, displacement = end - start, extra - start
    squared = float(np.dot(delta, delta))
    if squared == 0:
        return dict(cosine=None, projection_fraction=None, projection_raw=None, orthogonal_l2=None)
    norm = math.sqrt(squared)
    projection_raw = float(np.dot(displacement, delta)) / norm
    fraction = projection_raw / norm
    length = float(np.linalg.norm(displacement))
    cosine = max(-1.0, min(1.0, projection_raw / length)) if length > 0 else None
    return dict(cosine=cosine, projection_fraction=fraction, projection_raw=projection_raw,
                orthogonal_l2=float(np.linalg.norm(displacement - fraction * delta)))


GEOMETRY_COLUMNS = [
    "atlas_scorer", "atlas_status", "atlas_reason", "atlas_n_features", "atlas_axis_sha256",
    "atlas_distance_to_anchor_ref", "atlas_distance_to_anchor_alt", "atlas_anchor_delta_l2",
    "atlas_signed_cosine_from_anchor_ref", "atlas_projection_fraction_from_anchor_ref",
    "atlas_projection_raw_from_anchor_ref", "atlas_orthogonal_residual_l2",
    "atlas_estimated_increasing_base", "atlas_estimated_decreasing_base", "atlas_direction_status",
    "atlas_distance_to_estimated_increasing", "atlas_distance_to_estimated_decreasing",
    "atlas_signed_cosine_toward_estimated_increasing", "atlas_projection_fraction_toward_estimated_increasing",
    "atlas_projection_raw_toward_estimated_increasing", "atlas_fetch_sessions",
]
BREADTH_COLUMNS = ["atlas_scorer", "atlas_role", "atlas_base", "atlas_threshold", "atlas_n_features",
                   "atlas_n_above_threshold", "atlas_predicted_response_breadth", "atlas_axis_sha256"]


def summarize_geometry(plan_path, cache_dir, out_dir, effect_thresholds=(0.1, 0.5, 1.0)) -> dict:
    """Export descriptive geometry from complete matched axes; never query Atlas.

    s>0 means the mapped ancient anchor_alt is the estimated increasing base;
    s<0 means anchor_ref. Posterior is retained, not used to reverse this sign.
    At s=0, selection-oriented labels/metrics are absent. Zero anchor deltas
    retain distances and breadth but leave direction metrics undefined.
    """
    import numpy as np

    plan = _load_plan(plan_path)
    thresholds = sorted({float(value) for value in effect_thresholds})
    if not thresholds or any(not math.isfinite(value) or value <= 0 for value in thresholds):
        raise ValueError("effect_thresholds must be finite and positive")
    cache = Path(cache_dir).expanduser().resolve()
    if not (cache / "catalog.json").exists() and (cache / "cache" / "catalog.json").exists():
        cache /= "cache"
    if not (cache / "catalog.json").exists() and plan["n_variants"]:
        raise ValueError(f"No Atlas cache catalog at {cache}; execute an explicit query first")
    catalog = _read_json(cache / "catalog.json") if (cache / "catalog.json").exists() else _empty_catalog(plan)
    if catalog.get("settings") != _settings(plan) or catalog.get("namespace") != _hash(dict(settings=catalog.get("settings"), source=catalog.get("source"))):
        raise ValueError("Atlas cache settings or namespace differ from the plan")
    output = Path(out_dir).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest_path = output / "geometry.json"
    if manifest_path.exists() and _read_json(manifest_path).get("plan_id") != plan["plan_id"]:
        raise ValueError("Geometry output belongs to a different plan")
    geometry_rows, breadth_rows, counts, all_sessions = [], [], {}, set()
    for site in plan["sites"]:
        entries = {base: _cache_entry(cache, catalog, _variant(site, base))
                   for base in (site["anchor_ref"], site["anchor_alt"], site["extra_alt"])
                   if base != site["ref"]}
        for scorer in plan["scorers"]:
            row = dict(site, atlas_scorer=scorer, atlas_status="ok", atlas_reason="")
            vectors = {}
            for base, entry in entries.items():
                if entry is None:
                    row.update(atlas_status="missing_variant", atlas_reason=f"No complete cache for REF>{base}")
                    break
                vector, status, reason = _vector(entry, scorer)
                if status != "ok":
                    row.update(atlas_status=status, atlas_reason=f"REF>{base}: {reason}")
                    break
                vectors[base] = vector
            if row["atlas_status"] == "ok":
                exemplar = next(iter(vectors.values()))
                if any(value["axes"] != exemplar["axes"] for value in vectors.values()):
                    row.update(atlas_status="axis_mismatch", atlas_reason="Complete obs/var biological labels differ between alleles")
                elif any(value["interval"] != exemplar["interval"] for value in vectors.values()):
                    row.update(atlas_status="interval_mismatch", atlas_reason="Atlas scoring intervals differ between alleles")
            if row["atlas_status"] == "ok":
                raw = {base: value["values"] for base, value in vectors.items()}
                raw[site["ref"]] = np.zeros_like(exemplar["values"])
                anc_ref, anc_alt, extra = (raw[site[key]] for key in ("anchor_ref", "anchor_alt", "extra_alt"))
                n_features = len(extra)
                axis_hash = _hash(exemplar["axes"])
                sessions = sorted({value["session"] for value in vectors.values()})
                all_sessions.update(sessions)
                row.update(atlas_n_features=n_features, atlas_axis_sha256=axis_hash,
                           atlas_distance_to_anchor_ref=float(np.linalg.norm(extra - anc_ref)),
                           atlas_distance_to_anchor_alt=float(np.linalg.norm(extra - anc_alt)),
                           atlas_anchor_delta_l2=float(np.linalg.norm(anc_alt - anc_ref)),
                           atlas_fetch_sessions=",".join(sessions))
                direction = _direction_metrics(extra, anc_ref, anc_alt)
                row.update(atlas_signed_cosine_from_anchor_ref=direction["cosine"],
                           atlas_projection_fraction_from_anchor_ref=direction["projection_fraction"],
                           atlas_projection_raw_from_anchor_ref=direction["projection_raw"],
                           atlas_orthogonal_residual_l2=direction["orthogonal_l2"])
                if site["s"] == 0:
                    row["atlas_direction_status"] = "undefined_s_zero"
                else:
                    increasing, decreasing = ((site["anchor_alt"], site["anchor_ref"]) if site["s"] > 0
                                              else (site["anchor_ref"], site["anchor_alt"]))
                    oriented = _direction_metrics(extra, raw[decreasing], raw[increasing])
                    row.update(atlas_estimated_increasing_base=increasing, atlas_estimated_decreasing_base=decreasing,
                               atlas_direction_status="estimated_from_s_sign_posterior_retained",
                               atlas_distance_to_estimated_increasing=float(np.linalg.norm(extra - raw[increasing])),
                               atlas_distance_to_estimated_decreasing=float(np.linalg.norm(extra - raw[decreasing])),
                               atlas_signed_cosine_toward_estimated_increasing=oriented["cosine"],
                               atlas_projection_fraction_toward_estimated_increasing=oriented["projection_fraction"],
                               atlas_projection_raw_toward_estimated_increasing=oriented["projection_raw"])
                if row["atlas_anchor_delta_l2"] == 0:
                    row.update(atlas_status="degenerate_anchor_delta", atlas_reason="Anchor effect vectors coincide; directional geometry is undefined")
                for role, field in (("anchor_ref", "anchor_ref"), ("anchor_alt", "anchor_alt"), ("extra", "extra_alt")):
                    base = site[field]
                    for threshold in thresholds:
                        above = int(np.count_nonzero(np.abs(raw[base]) > threshold))
                        breadth_rows.append(dict(site, atlas_scorer=scorer, atlas_role=role, atlas_base=base,
                                                 atlas_threshold=threshold, atlas_n_features=n_features,
                                                 atlas_n_above_threshold=above,
                                                 atlas_predicted_response_breadth=above / n_features,
                                                 atlas_axis_sha256=axis_hash))
            counts[row["atlas_status"]] = counts.get(row["atlas_status"], 0) + 1
            geometry_rows.append(row)
    paths = {"geometry": output / "geometry.tsv", "breadth": output / "breadth.tsv"}
    for key, rows, columns in (("geometry", geometry_rows, GEOMETRY_COLUMNS), ("breadth", breadth_rows, BREADTH_COLUMNS)):
        fd, temporary = tempfile.mkstemp(prefix=key + ".", suffix=".tmp", dir=output)
        try:
            with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=plan["candidate_columns"] + columns,
                                        delimiter="\t", lineterminator="\n", extrasaction="raise")
                writer.writeheader()
                writer.writerows(rows)
            os.replace(temporary, paths[key])
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    result = dict(schema=SCHEMA, plan_id=plan["plan_id"], created_utc=_utc(), cache_dir=str(cache),
                  source=catalog["source"], cache_namespace=catalog["namespace"], source_caveat=SOURCE_CAVEAT,
                  caveat=GEOMETRY_CAVEAT, effect_thresholds=thresholds, threshold_rule="abs(raw_effect) > threshold",
                  n_sites=plan["n_sites"], n_geometry_rows=len(geometry_rows), n_breadth_rows=len(breadth_rows),
                  status_counts=counts, cache_fetch_sessions=sorted(all_sessions),
                  outputs={key: str(path) for key, path in paths.items()}, manifest_path=str(manifest_path),
                  files={key: _file_hash(path) for key, path in paths.items()})
    _json_write(manifest_path, result)
    return result
