#!/usr/bin/env python3
from __future__ import annotations
_description_assoc_pipeline = 'Plan and run the independent allele3 association DAG.\n\nThe planner uses only the Python standard library. It never reads a genotype\nmatrix, imports a scientific package, or starts a model on a submission host.\nActual stages use the existing batch_runner.py worker and durable receipts.'
#!/usr/bin/env python3
"""Plan and run the independent allele3 association DAG.

The planner uses only the Python standard library. It never reads a genotype
matrix, imports a scientific package, or starts a model on a submission host.
Actual stages use the existing batch_runner.py worker and durable receipts.
"""

import argparse
import copy
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from contextlib import contextmanager

ROOT = Path(__file__).resolve().parents[1]
SHARED = ROOT.parent / "0f"
sys.path.insert(0, str(SHARED))
from workflow_layout import named_directory
DEFAULT_CONFIG = ROOT / "config" / "burden.json"
SCHEMA = "allele3-assoc-workflow-1"
JSON_LIMIT = 16 * 1024 * 1024
SMALL_HASH_LIMIT = 1024 * 1024
NAME_RE = re.compile(r"[A-Za-z][A-Za-z0-9_]{0,99}\Z")
MODULES = ('burden.py',)


def now():
    return datetime.now(timezone.utc).isoformat()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False, ensure_ascii=False).encode()).hexdigest()


def read_json(path):
    path = Path(path)
    if path.stat().st_size > JSON_LIMIT:
        raise ValueError(f"Configuration/receipt exceeds {JSON_LIMIT} bytes: {path}")
    with path.open() as handle:
        text = handle.read(JSON_LIMIT + 1)
    if len(text) > JSON_LIMIT:
        raise ValueError(f"Configuration/receipt grew while reading: {path}")
    result = json.loads(text)
    if not isinstance(result, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return result


def merge_config(base, override):
    result = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = merge_config(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def text(value, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be nonempty text")
    if any(ord(c) < 32 or ord(c) == 127 for c in value):
        raise ValueError(f"Control character in {label}")
    return value


def boolean(value):
    if isinstance(value, bool):
        return value
    if str(value).upper() in {"TRUE", "1", "YES"}:
        return True
    if str(value).upper() in {"FALSE", "0", "NO"}:
        return False
    raise argparse.ArgumentTypeError("Expected TRUE or FALSE")


def chromosomes(value):
    values = value if isinstance(value, list) else str(value).split(",")
    result = []
    for token in values:
        token = str(token).strip().removeprefix("chr")
        if re.fullmatch(r"\d+-\d+", token):
            start, stop = map(int, token.split("-"))
            if stop < start:
                raise ValueError("Reversed chromosome range")
            items = list(range(start, stop + 1))
        elif token.isdigit():
            items = [int(token)]
        else:
            raise ValueError("assoc chromosome scope must be autosomes 1-22, e.g. 1-22 or 1,2")
        for item in items:
            if item not in range(1, 23):
                raise ValueError("assoc currently requires autosomal diploid genotypes (chromosomes 1-22)")
            if str(item) in result:
                raise ValueError(f"Repeated chromosome in scope: {item}")
            result.append(str(item))
    if not result:
        raise ValueError("Empty chromosome scope")
    return result


def outcomes(value):
    values = value if isinstance(value, list) else str(value).split(",")
    result = [str(v).strip() for v in values]
    if not result or any(not NAME_RE.fullmatch(v) for v in result):
        raise ValueError("--Y requires comma-separated outcome keys (letters, numbers, underscores)")
    if len(set(result)) != len(result):
        raise ValueError("Repeated outcome in --Y")
    return result


def absolute(value, base):
    if value is None or value == "":
        return None
    value = os.path.expandvars(text(value, "path"))
    if "${" in value:
        raise ValueError(f"Unresolved environment/root token in path: {value}")
    path = Path(value).expanduser()
    return str((path if path.is_absolute() else Path(base) / path).absolute())


def executable(value, base):
    value = text(value, "scientific Python")
    if "/" not in value and "\\" not in value:
        return shutil.which(value) or value
    return absolute(value, base)


def shared_roots(mode):
    # Use the same authoritative table as 0f/0os.sh and LE8. No search for data.
    os.environ["submit"] = mode
    os.environ["SUBMIT"] = mode
    sys.path.insert(0, str(SHARED))
    from workflow_paths import resolved_table
    return resolved_table()


def expand_roots(value, roots):
    if isinstance(value, dict):
        return {k: expand_roots(v, roots) for k, v in value.items()}
    if isinstance(value, list):
        return [expand_roots(v, roots) for v in value]
    if isinstance(value, str):
        return re.sub(r"\$\{([A-Z][A-Z0-9_]*)\}",
                      lambda m: roots.get(m.group(1), os.environ.get(m.group(1), m.group(0))), value)
    return value


def resolve(a):
    sys.path.insert(0, str(SHARED))
    from workflow_mode import resolve_mode
    from ukb_followup import resolve_followup, validate_record
    mode = resolve_mode(a.submit)
    roots = shared_roots(mode)
    cfg_path = Path(a.config or DEFAULT_CONFIG).expanduser().absolute()
    fu = resolve_followup()
    default_config, supplied = read_json(DEFAULT_CONFIG), read_json(cfg_path)
    for value in (default_config, supplied):
        validate_record(value, "phenotypes", "end_date", fu)
    config = merge_config(default_config, supplied)
    config["phenotypes"]["end_date"] = fu["date"]
    config["followup_provenance"] = fu
    config = expand_roots(config, roots)
    for section in ("runtime", "catalog", "exposures", "phenotype_source", "phenotypes", "models"):
        if not isinstance(config.get(section), dict):
            raise ValueError(f"config.{section} must be an object")
    base = cfg_path.parent
    for section, keys in {
        "catalog": ("counts", "sets_bed"),
        "exposures": ("pgen_template", "vcf_template", "vcf", "samples", "upstream_qc", "dosage_sites"),
        "phenotype_source": ("cohort", "dir_ukb_phe", "unrelated_keep"),
    }.items():
        for key in keys:
            config[section][key] = absolute(config[section].get(key), base)
    cli_paths = {
        "counts": ("catalog", "counts"), "sets_bed": ("catalog", "sets_bed"),
        "pgen_template": ("exposures", "pgen_template"), "vcf_template": ("exposures", "vcf_template"),
        "vcf": ("exposures", "vcf"), "samples": ("exposures", "samples"),
        "upstream_qc": ("exposures", "upstream_qc"), "dosage_sites": ("exposures", "dosage_sites"),
        "cohort": ("phenotype_source", "cohort"), "dir_ukb_phe": ("phenotype_source", "dir_ukb_phe"),
        "unrelated_keep": ("phenotype_source", "unrelated_keep"),
    }
    if any(getattr(a, key) is not None for key in ("pgen_template", "vcf_template", "vcf")):
        for key in ("pgen_template", "vcf_template", "vcf"):
            config["exposures"][key] = None
    if a.dir_ukb_phe is not None:
        config["phenotype_source"]["cohort"] = None
    for option, (section, key) in cli_paths.items():
        if getattr(a, option) is not None:
            config[section][key] = absolute(getattr(a, option), Path.cwd())
    for option in ("source_cohort", "population", "source_scope"):
        if getattr(a, option) is not None:
            config["catalog"][option] = getattr(a, option)
    if a.build is not None:
        config["catalog"]["build"] = a.build
        config["exposures"]["build"] = a.build
    config["outcomes"] = outcomes(a.Y if a.Y is not None else config["outcomes"])
    config["exposures"]["chromosomes"] = chromosomes(
        a.chromosomes if a.chromosomes is not None else config["exposures"].get("chromosomes", "1-22"))
    config["data_label"] = a.data_label or config.get("data_label", "REAL")
    if config["data_label"] not in {"REAL", "SIMULATED"}:
        raise ValueError("data_label must be REAL or SIMULATED")
    runtime = config["runtime"]
    for key in ("cores", "memory_limit_gb", "hours", "queue"):
        if getattr(a, key) is not None:
            runtime[key] = getattr(a, key)
    chosen = a.python or runtime.get("python") or os.environ.get("ALLELE_ASSOC_PYTHON") or os.environ.get("ALLELE_PYTHON") or os.environ.get("PYTHON_BIN")
    if not chosen:
        chosen = str(Path(roots["DIR_WORK"]) / ".conda/envs/allele3-evolution/bin/python") if mode == "hpc" else sys.executable
    runtime["python"] = executable(chosen, base)
    out = Path(absolute(a.out, Path.cwd()) or str(Path(roots["ANALYSIS_ROOT"]) / "allele3/assoc"))
    from batch_runner import integer, positive
    runtime["cores"] = integer(runtime["cores"], "cores")
    runtime["memory_limit_gb"] = float(positive(runtime["memory_limit_gb"], "memory_limit_gb"))
    runtime["hours"] = float(positive(runtime["hours"], "hours"))
    jobs_n = integer(a.jobs_N if a.jobs_N is not None else os.environ.get("jobs_N", "8" if mode == "hpc" else "1"),
                     "jobs_N", maximum=10000 if mode == "hpc" else 8)
    if mode == "local" and runtime["memory_limit_gb"] > 32:
        raise ValueError("Local --memory-limit-gb must be <=32 GiB; larger analyses require an HPC allocation")
    job_memory = runtime.get("job_memory_gb")
    if job_memory is None:
        job_memory = runtime["memory_limit_gb"] / jobs_n if mode == "local" else runtime["memory_limit_gb"]
    runtime["job_memory_gb"] = float(positive(job_memory, "runtime.job_memory_gb"))
    if runtime["job_memory_gb"] > runtime["memory_limit_gb"]:
        raise ValueError("runtime.job_memory_gb exceeds --memory-limit-gb")
    wait = a.wait if a.wait is not None else boolean(os.environ.get("SUBMIT_WAIT", "FALSE"))
    return config, cfg_path, out, mode, jobs_n, wait, roots


def configuration_issues(config):
    errors = []
    catalog, exposure, source = (config[k] for k in ("catalog", "exposures", "phenotype_source"))
    for key in ("counts", "source_cohort", "population", "catalog_id"):
        if not catalog.get(key):
            errors.append(f"Configure catalog.{key} (no UKB WGS/counts path is inferred)")
    if sum(bool(exposure.get(k)) for k in ("pgen_template", "vcf_template", "vcf")) != 1:
        errors.append("Configure exactly one exposures.pgen_template, exposures.vcf_template, or exposures.vcf; supply an actual WGS source")
    if catalog.get("build") not in {"GRCh37", "GRCh38"} or exposure.get("build") != catalog.get("build"):
        errors.append("catalog.build and exposures.build must explicitly match (GRCh37 or GRCh38); no implicit liftOver")
    if exposure.get("evidence_kind") != "direct_gt":
        errors.append("Association exposures require evidence_kind=direct_gt; imputed/dosage evidence is not substituted for multi-allelic GT")
    if exposure.get("pgen_template") and not exposure.get("upstream_qc"):
        errors.append("PGEN requires exposures.upstream_qc documenting upstream genotype QC; PGEN does not contain per-call DP/GQ/AD")
    if not source.get("cohort") and not source.get("dir_ukb_phe"):
        errors.append("Configure phenotype_source.cohort or --dir-ukb-phe (DIR_UKB/Rdata/all.rds)")
    if config["models"].get("relatedness_policy", "require_unrelated") == "require_unrelated" and not source.get("unrelated_keep"):
        errors.append("models.relatedness_policy=require_unrelated requires phenotype_source.unrelated_keep / --unrelated-keep; no kinship adjustment is invented")
    return errors


def template_path(template, chromosome):
    value = template.replace("{chrom}", chromosome)
    if "{" in value or "}" in value:
        raise ValueError("Only {chrom} is supported in a genotype template: " + template)
    return value


def input_specs(config, cfg_path):
    """Only metadata/stat calls; do not open genotype streams or cohort tables."""
    specs = {}

    def add(path, label, required=True, hash_small=True):
        if not path:
            return
        path = str(Path(path).expanduser().absolute())
        old = specs.get(path)
        item = dict(path=path, label=label, required=bool(required), hash_small=bool(hash_small))
        if old:
            item["required"] = item["required"] or old["required"]
            item["hash_small"] = item["hash_small"] or old["hash_small"]
        specs[path] = item

    add(DEFAULT_CONFIG, "default configuration")
    add(cfg_path, "requested configuration")
    for path in [Path(__file__), ROOT / "allele3.sh", ROOT / "f/batch_runner.py",
                 ROOT / "f/analysis_integrity.py", SHARED / "qiming_resources.py",
                 SHARED / "0os.sh", SHARED / "workflow_env.sh", SHARED / "workflow_paths.py",
                 SHARED / "ukb_followup.py", SHARED / "phenotype.R"]:
        add(path, "workflow source")
    for filename in MODULES:
        add(ROOT / "f" / filename, "scientific stage source")
    for path in sorted((ROOT / "f").glob("burden*.py")):
        add(path, "association source")
    catalog, exposure, source = (config[k] for k in ("catalog", "exposures", "phenotype_source"))
    for key in ("counts", "sets_bed"):
        add(catalog.get(key), "catalog." + key)
    for key in ("samples", "upstream_qc", "dosage_sites"):
        add(exposure.get(key), "exposures." + key)
    for chrom in exposure["chromosomes"]:
        if exposure.get("pgen_template"):
            prefix = template_path(exposure["pgen_template"], chrom)
            add(prefix + ".pgen", "WGS PGEN chr" + chrom, hash_small=False)
            add(prefix + ".psam", "WGS PSAM chr" + chrom)
            plain, compressed = Path(prefix + ".pvar"), Path(prefix + ".pvar.zst")
            if plain.is_file() and compressed.is_file():
                raise ValueError("Ambiguous PVAR: both .pvar and .pvar.zst exist for " + prefix)
            # Freeze both possible paths, including absence of the unused one.
            add(plain, "WGS PVAR chr" + chrom, required=not compressed.is_file(), hash_small=False)
            add(compressed, "WGS compressed PVAR chr" + chrom, required=compressed.is_file(), hash_small=False)
        elif exposure.get("vcf_template") or exposure.get("vcf"):
            path = template_path(exposure.get("vcf_template") or exposure["vcf"], chrom)
            add(path, "WGS VCF chr" + chrom, hash_small=False)
            # Indices enable bounded fetch. Streaming is supported without one;
            # adding/removing an index still invalidates an old frozen plan.
            add(path + ".tbi", "optional VCF tabix index", required=False, hash_small=False)
            add(path + ".csi", "optional VCF CSI index", required=False, hash_small=False)
    if source.get("cohort"):
        add(source["cohort"], "portable phenotype cohort")
        phenotype_file = Path(source["cohort"])
    elif source.get("dir_ukb_phe"):
        phenotype_file = Path(source["dir_ukb_phe"]) / "Rdata/all.rds"
        add(phenotype_file, "LE8 phenotype source")
    else:
        phenotype_file = None
    if phenotype_file is not None and config["phenotypes"].get("bbc_columns"):
        # The bridge consults this exact sidecar only when requested BBC fields
        # are absent from the primary cohort. Freeze presence/absence without
        # opening either table on the planner host.
        add(phenotype_file.parent / "bbc.rds", "optional BBC phenotype sidecar", required=False)
    add(source.get("unrelated_keep"), "unrelated sample keep file")
    python = config["runtime"]["python"]
    add(shutil.which(python) or python, "scientific Python interpreter", hash_small=False)
    return list(specs.values())


def preflight(config, cfg_path, out, mode, jobs_n):
    errors = configuration_issues(config)
    specs = input_specs(config, cfg_path)
    checked = []
    for item in specs:
        path = Path(item["path"])
        state = "present" if path.is_file() else "missing_or_not_file"
        if state == "present" and (path.stat().st_size == 0 or not os.access(path, os.R_OK)):
            state = "empty_or_unreadable"
        if item["required"] and state != "present":
            errors.append(f"{item['label']}: {state}: {path}")
        checked.append(dict(path=str(path), label=item["label"], required=item["required"], status=state))
    python = config["runtime"]["python"]
    if not shutil.which(python):
        errors.append("Scientific Python is not executable on this host: " + python)
    from analysis_integrity import InputPath, check_output_paths
    protected = [InputPath(item["path"]) for item in specs]
    source = config["phenotype_source"]
    if not source.get("cohort") and source.get("dir_ukb_phe"):
        protected.append(InputPath(source["dir_ukb_phe"], "phenotype directory"))
    check_output_paths(out, protected)
    resources = dict(queue=str(config["runtime"]["queue"]), cores=config["runtime"]["cores"],
                     memory_gib=config["runtime"]["job_memory_gb"], hours=config["runtime"]["hours"])
    from batch_runner import normalize_resources
    normalized = normalize_resources({}, {"resources": resources, "lsf_memory_scope": "host"}, mode)
    warnings = ["Preflight checks configuration and file metadata only; genotype/cohort schemas and scientific packages are checked by compute stages.",
                "The phenotype follow-up end date must match the supplied cohort release."]
    if mode == "local":
        warnings.append("Local memory is a total admission budget, bounded by existing system/cgroup limits; this workflow does not create a new process-tree hard cap.")
        warnings.append("Local per-job admission memory defaults to total budget/jobs_N. Use jobs_N=1 or explicit runtime.job_memory_gb for memory-heavy cohort tables.")
    if config["models"].get("relatedness_policy") == "fixed_effect_exploratory":
        warnings.append("Explicit fixed_effect_exploratory analysis is not a mixed model or a kinship-adjusted confirmatory test.")
    if config["exposures"].get("vcf") and len(config["exposures"]["chromosomes"]) > 1:
        path = config["exposures"]["vcf"]
        if not Path(path + ".tbi").is_file() and not Path(path + ".csi").is_file():
            warnings.append("A single unindexed VCF will be streamed once per chromosome shard; use chromosome VCFs or a BGZF/indexed source for large WGS.")
    return dict(schema=SCHEMA, status="ready_metadata_only" if not errors else "blocked",
                mode=mode, out=str(out), outcomes=config["outcomes"], chromosomes=config["exposures"]["chromosomes"],
                jobs_N=jobs_n, resources=normalized, memory_limit_gb=config["runtime"]["memory_limit_gb"],
                job_memory_gb=config["runtime"]["job_memory_gb"], data_label=config["data_label"],
                errors=errors, warnings=warnings, files=checked), specs


def stage_specs(config, out, config_file):
    """Return the explicit scientific DAG, before the batch worker guard."""
    py = config["runtime"]["python"]
    stages = []
    cat_dir, exp_dir = out / "catalog", out / "exposures"
    catalog = cat_dir / "catalog.tsv.gz"
    sets = cat_dir / "site_sets.tsv.gz"
    phenotype = out / "phenotypes/phenotypes.tsv.gz"

    def push(name, module, args, deps, outputs):
        stages.append(dict(id=name, argv=[py, str(ROOT / "f" / "burden.py"), module, *[str(v) for v in args]],
                           deps=deps, outputs=[str(p) for p in outputs]))

    cat_args = []
    for key in ("counts", "population", "build", "source_cohort", "catalog_id", "source_scope", "min_ac", "min_call_rate", "sets_bed"):
        value = config["catalog"].get(key)
        if value is not None:
            cat_args.extend(["--" + key.replace("_", "-"), value])
    push("catalog", 'catalog', [*cat_args, "--out", cat_dir], [],
         [catalog, sets, cat_dir / "catalog_audit.tsv.gz", cat_dir / "manifest.json"])
    exposure = config["exposures"]
    exp_args = ["--catalog", catalog, "--site-sets", sets]
    for key in ("build", "cohort_id", "evidence_kind", "pgen_template", "vcf_template", "vcf", "samples", "upstream_qc",
                "min_dp", "min_gq", "min_allele_balance", "min_sample_call_rate", "max_missing_fraction", "dosage_sites", "max_dosage_rows", "max_aggregation_mb"):
        if exposure.get(key) is not None:
            exp_args.extend(["--" + key.replace("_", "-"), exposure[key]])
    shards = []
    for chromosome in exposure["chromosomes"]:
        stage_id = "exposure_chr" + chromosome.zfill(2)
        target = out / "exposure_shards" / ("chr" + chromosome)
        outputs = [target / "exposures.tsv.gz", target / "riskset.tsv.gz", target / "site_qc.tsv.gz",
                   target / "samples.tsv", target / "manifest.json"]
        if exposure.get("dosage_sites"):
            outputs.append(target / "locus_dosages.tsv.gz")
        push(stage_id, 'exposure', [*exp_args, "--chromosomes", chromosome, "--out", target], ["catalog"], outputs)
        shards.append((stage_id, target))
    merge_args = ["merge", "--catalog", catalog, "--site-sets", sets,
                  "--chromosomes", ",".join(exposure["chromosomes"])]
    for _, target in shards:
        merge_args.extend(["--shard-dir", target])
    if exposure.get("max_aggregation_mb") is not None:
        merge_args.extend(["--max-aggregation-mb", exposure["max_aggregation_mb"]])
    merge_outputs = [exp_dir / "exposures.tsv.gz", exp_dir / "riskset.tsv.gz", exp_dir / "site_qc.tsv.gz",
                     exp_dir / "samples.tsv", exp_dir / "manifest.json"]
    if exposure.get("dosage_sites"):
        merge_outputs.append(exp_dir / "locus_dosages.tsv.gz")
    push("exposures", 'exposure', [*merge_args, "--out", exp_dir], [name for name, _ in shards], merge_outputs)
    source = config["phenotype_source"]
    phe_args = ["--Y", ",".join(config["outcomes"]), "--config", config_file, "--data-label", config["data_label"]]
    phe_args.extend(["--cohort", source["cohort"]] if source.get("cohort") else ["--dir-ukb-phe", source["dir_ukb_phe"]])
    push("phenotypes", 'phenotype', [*phe_args, "--out", phenotype], [],
         [phenotype, Path(str(phenotype) + ".manifest.json"), Path(str(phenotype) + ".audit.tsv")])
    model_ids = []
    for outcome in config["outcomes"]:
        model_id = "model_" + outcome
        target = out / "models" / outcome
        args = ["--phenotypes", phenotype, "--exposures", exp_dir / "exposures.tsv.gz", "--Y", outcome,
                "--config", config_file, "--data-label", config["data_label"], "--out-dir", target]
        if source.get("unrelated_keep"):
            args.extend(["--unrelated-keep", source["unrelated_keep"]])
        if exposure.get("dosage_sites"):
            args.extend(["--locus-dosages", exp_dir / "locus_dosages.tsv.gz"])
        push(model_id, 'model', args, ["exposures", "phenotypes"],
             [target / "association_results.tsv", target / "cohort_flow.tsv", target / "model_manifest.json", target / "analysis_families.tsv"])
        model_ids.append(model_id)
    report_dir = out / "reports"
    push("report", 'report', ["--results-root", out / "models", "--Y", ",".join(config["outcomes"]), "--out-dir", report_dir], model_ids,
         [report_dir / filename for filename in ("association_results.tsv", "analysis_families.tsv", "forest.png", "forest.svg", "report.html", "manifest.json")])
    return stages


def immutable_json(path, value):
    from batch_runner import atomic_json
    path = Path(path)
    if path.exists():
        if read_json(path) != value:
            raise ValueError("Immutable request changed: " + str(path))
    else:
        atomic_json(path, value)


@contextmanager
def output_lock(out):
    out.mkdir(parents=True, exist_ok=True)
    with (out / ".assoc.workflow.lock").open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError("Another association planner owns this output directory") from None
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def ensure_previous_finished(out, run_key):
    marker = out / "runtime/burden/current.json"
    if not marker.is_file():
        return
    previous = read_json(marker)
    if previous.get("run_key") == run_key:
        return
    prior_dir = Path(previous["plan_dir"])
    state_path = prior_dir / "submission.json"
    if not state_path.exists():
        return  # A written-only plan started no payload.
    state = read_json(state_path)
    if state.get("mode") == "local" and state.get("status") in {"complete", "failed"}:
        return
    jobs = state.get("jobs", {})
    if not jobs:
        return
    frozen = read_json(prior_dir / "plan.json")
    complete = True
    for job in frozen["jobs"]:
        path = prior_dir / "runs" / (job["id"] + ".json")
        if not path.is_file():
            complete = False
            break
        row = read_json(path)
        submitted = jobs.get(job["id"], {})
        if row.get("status") != "succeeded" or (state.get("mode") == "hpc" and str(row.get("lsf_job_id")) != str(submitted.get("job_id"))):
            complete = False
            break
    if not complete:
        raise RuntimeError("This --out has a different unfinished/failed HPC or unresolved workflow. Use a fresh --out or resolve the existing jobs first: " + str(prior_dir))


def prepare_workflow(config, cfg_path, out, mode, jobs_n, specs, plan_dir=None):
    from batch_runner import atomic_json, normalize_plan, path_identity
    identities, declared = [], []
    for item in specs:
        path = Path(item["path"])
        full = item["hash_small"] and path.is_file() and path.stat().st_size <= SMALL_HASH_LIMIT
        identities.append(path_identity(path, full))
        declared.append(dict(path=str(path), sha256=full))
    runtime = config["runtime"]
    environment = {key: value for key, value in os.environ.items()
                   if key.startswith(("DIR_", "_OS_DEFAULT_", "OS_")) or key in {
                       "PATH", "SCRIPTS_ROOT", "ANALYSIS_ROOT", "DATE_SCRATCH", "PYTHONPATH", "LD_LIBRARY_PATH", "R_BIN"}}
    environment.update(PYTHONDONTWRITEBYTECODE="1", PYTHONUNBUFFERED="1", OMP_NUM_THREADS=str(runtime["cores"]))
    run_key = digest(dict(schema=SCHEMA, config=config, mode=mode, jobs_N=jobs_n, inputs=identities,
                          environment=environment, driver_python=sys.executable))
    environment["ALLELE_ASSOC_RUN_KEY"] = run_key
    ensure_previous_finished(out, run_key)
    request_root = out / "runtime/burden/requests"
    request_dir = named_directory(request_root, run_key, "run")
    config_file = request_dir / "config.json"
    immutable_json(config_file, config)
    stages = stage_specs(config, out, config_file)
    stage_root = out / "runtime/burden/stages" / request_dir.name
    request = dict(schema=SCHEMA, run_key=run_key, mode=mode, out=str(out), config_file=str(config_file),
                   config=config, input_snapshot=identities, stages=stages)
    request_file = request_dir / "request.json"
    immutable_json(request_file, request)
    declared.extend([dict(path=str(config_file), sha256=True), dict(path=str(request_file), sha256=True)])
    jobs = [dict(id=stage["id"], deps=stage["deps"], cwd=str(ROOT), env=environment,
                 argv=[runtime["python"], str(Path(__file__).resolve()),"stage", "--request", str(request_file), "--stage-id", stage["id"]],
                 outputs=[*stage["outputs"], str(stage_root / (stage["id"] + ".json"))]) for stage in stages]
    raw = dict(name="allele3.burden", driver_python=sys.executable, max_running=jobs_n,
               local_memory_gib=runtime["memory_limit_gb"], lsf_memory_scope="host",
               resources=dict(queue=str(runtime["queue"]), cores=runtime["cores"],
                              memory_gib=runtime["job_memory_gb"], hours=runtime["hours"]),
               inputs=declared, jobs=jobs, provenance=dict(schema=SCHEMA, run_key=run_key,
               data_label=config["data_label"], followup_provenance=config["followup_provenance"], source_config=str(cfg_path), resolved_config=str(config_file),
               request=str(request_file), source_scope=config["catalog"].get("source_scope"),
               memory_policy="local admission budget; HPC LSF reservation; no new hard cap asserted"))
    plan = normalize_plan(raw, mode)
    plan_dir = Path(plan_dir).expanduser().absolute() if plan_dir else named_directory(out / "jobs", plan["identity"], "run")
    from analysis_integrity import InputPath, check_output_paths
    check_output_paths(plan_dir, [*[InputPath(item["path"]) for item in specs],
                                  InputPath(config_file), InputPath(request_file)])
    if (plan_dir / "plan.json").exists():
        previous = read_json(plan_dir / "plan.json")
        if previous.get("identity") != plan["identity"]:
            raise ValueError("The requested --plan-dir belongs to a different frozen plan; choose a new directory: " + str(plan_dir))
    current = dict(schema=SCHEMA, run_key=run_key, plan_identity=plan["identity"], plan_dir=str(plan_dir),
                   request=str(request_file), mode=mode, status="prepared", updated_at=now())
    atomic_json(out / "runtime/burden/current.json", current)
    return plan, plan_dir, current


def write_review_commands(plan_dir, request_file):
    """Readable payloads; submitting must still use receipt-aware submit.sh."""
    request = read_json(request_file)
    lines = ["#!/usr/bin/env bash", "# Scientific command review only. All commands below are comments.",
             "# Execute submit.sh in this directory to retain allocations, dependencies and receipts.",
             "# The resolved request is " + shlex.quote(str(request_file)), ""]
    for stage in request["stages"]:
        lines.extend(["# Stage: " + stage["id"] + "; dependencies: " + (", ".join(stage["deps"]) or "none"),
                      "# " + shlex.join(stage["argv"]), ""])
    target = Path(plan_dir) / "payloads.review.sh"
    target.write_text("\n".join(lines) + "\n")


def verify_source_snapshot(request):
    from batch_runner import path_identity
    for expected in request["input_snapshot"]:
        actual = path_identity(expected["path"], "sha256" in expected)
        if actual != expected:
            raise ValueError("Association input/source changed since planning: " + expected["path"])


def run_stage(argv):
    p = argparse.ArgumentParser(description="Internal association stage; use the generated batch plan")
    p.add_argument("--request", required=True)
    p.add_argument("--stage-id", required=True)
    a = p.parse_args(argv)
    request = read_json(a.request)
    if request.get("schema") != SCHEMA:
        raise ValueError("Unsupported association request")
    out = Path(request["out"])
    stages = {item["id"]: item for item in request["stages"]}
    if a.stage_id not in stages:
        raise ValueError("Unknown association stage")
    if os.environ.get("ALLELE_BATCH_WORKER") != "1":
        raise RuntimeError("Scientific stages must run through batch_runner worker receipts; use the public assoc command")
    if request["mode"] == "hpc" and not re.fullmatch(r"[1-9][0-9]*", os.environ.get("LSB_JOBID", "")):
        raise RuntimeError("HPC scientific payload requires an LSF allocation")
    if request["mode"] == "local" and not os.environ.get("LSB_JOBID") and shutil.which("bsub"):
        raise RuntimeError("Local scientific payload refused on an LSF submission host without an allocation")
    current = read_json(out / "runtime/burden/current.json")
    if current.get("run_key") != request["run_key"]:
        raise RuntimeError("This plan no longer owns --out; an old .cmd cannot overwrite a different association request")
    stage = stages[a.stage_id]
    from batch_runner import atomic_json, path_identity, validate_output_file
    receipts = out / "runtime/burden/stages" / Path(a.request).parent.name
    receipts.mkdir(parents=True, exist_ok=True)
    receipt = receipts / (a.stage_id + ".json")
    row = dict(schema="allele3-assoc-stage-1", run_key=request["run_key"], stage=a.stage_id,
               status="running", started_at=now(), argv=stage["argv"], dependencies=stage["deps"],
               input_snapshot=request["input_snapshot"], lsf_job_id=os.environ.get("LSB_JOBID"))
    atomic_json(receipt, row)
    try:
        verify_source_snapshot(request)
        result = subprocess.run(stage["argv"], cwd=ROOT, check=False)
        if result.returncode:
            row.update(status="failed", returncode=result.returncode)
        else:
            verify_source_snapshot(request)
            output_ids = [path_identity(path) for path in stage["outputs"]]
            if any("size" not in item for item in output_ids):
                raise ValueError("Stage exited zero without its complete declared outputs")
            for item in output_ids:
                validate_output_file(item["path"])
            row.update(status="succeeded", returncode=0, outputs=output_ids)
    except (OSError, ValueError, RuntimeError) as exc:
        row.update(status="failed", returncode=2, error=str(exc))
    finally:
        row["finished_at"] = now()
        atomic_json(receipt, row)
    return row["returncode"]


def parser():
    p = argparse.ArgumentParser(description=_description_assoc_pipeline, formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Workflow: catalog -> chromosome exposure shards -> exposure merge;
phenotypes -> per-outcome models -> combined report. Failed dependencies block
their downstream jobs. Re-running an unchanged request validates and reuses
completed receipts; changed inputs receive a different frozen plan.

Usage Examples (optional native workflow; configure actual counts and WGS first):
  bash allele3.sh burden --config assoc.site.json --preflight
  bash allele3.sh burden --config assoc.site.json --Y cvd_cad,ra,depression --plan-only
  bash allele3.sh burden --config assoc.site.json --wait
  bash allele3.sh burden --config assoc.site.json --submit local --jobs-N 2
  bash allele3.sh burden demo --out /absolute/path/demo

--plan-only/--write-cmd never submits or executes a scientific payload.
Local memory is an admission budget, not a newly established process-tree cap.
Local job memory defaults to --memory-limit-gb / --jobs-N; config runtime.job_memory_gb
overrides that reservation. Use --jobs-N 1 for memory-heavy cohort tables.
HPC job memory defaults to --memory-limit-gb and uses the shared QiMing rounding.
PCA, HERV and genome-scan entry points remain independent of this workflow.
""")
    p.add_argument("--config", help="JSON override merged with config/burden.json; relative paths use this file's directory")
    p.add_argument("--out", "--dir-out", dest="out", help="Persistent output [ANALYSIS_ROOT/allele3/assoc]")
    p.add_argument("--Y", "--trait", dest="Y", help="Outcome keys, default cvd_cad,ra,depression")
    p.add_argument("--submit", choices=("auto", "local", "hpc"), help="Optional override; detect LSF automatically")
    p.add_argument("--plan-only", "--write-cmd", "--dry-run", action="store_true", dest="plan_only")
    p.add_argument("--preflight", action="store_true", help="Metadata-only validation; return nonzero for missing/unconfigured inputs")
    p.add_argument("--cores", type=int)
    p.add_argument("--memory-limit-gb", type=float, help="Local total admission budget (<=32 GiB); HPC default per-job memory")
    p.add_argument("--hours", type=float)
    p.add_argument("--queue")
    p.add_argument("--jobs-N", type=int, help="Concurrent local jobs (max 8) / active HPC lane ceiling")
    p.add_argument("--wait", nargs="?", const=True, type=boolean, help="Wait for submitted HPC stages [FALSE]")
    p.add_argument("--plan-dir", "--job-dir", dest="plan_dir", help="Optional exact frozen-plan directory; changed plans require a new directory")
    p.add_argument("--retry-failed", action="store_true", help="Explicit LSF EXIT retry with dependency repair, through batch_runner")
    p.add_argument("--python", help="Scientific interpreter; planner itself needs only Python stdlib")
    p.add_argument("--counts", help="Joint allele-count source; split per-ALT summaries are not silently accepted")
    p.add_argument("--source-cohort")
    p.add_argument("--population")
    p.add_argument("--source-scope")
    p.add_argument("--build", choices=("GRCh37", "GRCh38"), help="Verified common build for counts and genotype source")
    p.add_argument("--sets-bed")
    genotypes = p.add_mutually_exclusive_group()
    genotypes.add_argument("--pgen-template", help="Explicit WGS prefix, optionally containing {chrom}")
    genotypes.add_argument("--vcf-template", help="Explicit WGS VCF path containing {chrom}")
    genotypes.add_argument("--vcf", help="Explicit multi-allelic VCF; an index is recommended for large data")
    p.add_argument("--chromosomes", help="Autosome scope, default 1-22; comma lists and numeric ranges supported")
    p.add_argument("--samples", help="Fixed sample roster with sample_id, shared by all shards")
    p.add_argument("--upstream-qc", help="Documented upstream genotype QC; required for PGEN")
    p.add_argument("--dosage-sites", help="Prespecified bounded locus-level extraction list")
    phenotypes = p.add_mutually_exclusive_group()
    phenotypes.add_argument("--dir-ukb-phe", help="LE8 phenotype directory [DIR_UKB]; uses Rdata/all.rds")
    phenotypes.add_argument("--cohort", help="Portable, explicitly supplied phenotype cohort table")
    p.add_argument("--unrelated-keep", help="Explicit unrelated sample list for the primary fixed-effect analysis")
    p.add_argument("--data-label", choices=("REAL", "SIMULATED"))
    return p


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    try:
        if argv and argv[0] == "stage":
            return run_stage(argv[1:])
        if argv and argv[0] == "run":
            argv = argv[1:]
        a = parser().parse_args(argv)
        config, cfg_path, out, mode, jobs_n, should_wait, roots = resolve(a)
        check, specs = preflight(config, cfg_path, out, mode, jobs_n)
        if a.preflight:
            print(json.dumps(check, indent=2, ensure_ascii=False, allow_nan=False))
            return 0 if not check["errors"] else 2
        problems = configuration_issues(config)
        if problems:
            raise ValueError("Missing/invalid configuration:\n  - " + "\n  - ".join(problems))
        if check["errors"] and not a.plan_only:
            raise ValueError("Preflight failed; no scientific jobs submitted:\n  - " + "\n  - ".join(check["errors"]))
        if check["errors"]:
            print("[assoc] plan contains unresolved file checks; rerun preflight on the execution host before submission", file=sys.stderr)
        import batch_runner
        from batch_runner import atomic_json, dispatch_frozen
        with output_lock(out):
            plan, plan_dir, current = prepare_workflow(config, cfg_path, out, mode, jobs_n, specs, a.plan_dir)
            atomic_json(plan_dir / "preflight.json", check)
            write_review_commands(plan_dir, current["request"])
            # Waiting must not hold the output planning lock for hours. Worker
            # receipts and the current-run guard protect the scientific outputs.
            result = dispatch_frozen(plan, plan_dir, write_only=a.plan_only, jobs_n=jobs_n,
                                     wait=False, retry_failed=a.retry_failed)
            current.update(status=result.get("status", "submitted"), updated_at=now())
            atomic_json(out / "runtime/burden/current.json", current)
        if should_wait and mode == "hpc" and not a.plan_only:
            from batch_runner import wait_for_plan
            wait_for_plan(plan, plan_dir, result)
            current.update(status="complete", updated_at=now())
            atomic_json(out / "runtime/burden/current.json", current)
        print(json.dumps(dict(schema=SCHEMA, status=current["status"], mode=mode, out=str(out),
                              plan_dir=str(plan_dir), run_key=current["run_key"],
                              jobs=len(plan["jobs"]), outcomes=config["outcomes"],
                              unresolved_checks=len(check["errors"])), sort_keys=True))
        return 0
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as exc:
        print("allele3-assoc: " + str(exc), file=sys.stderr)
        return 2


def _entry_assoc_pipeline():
    raise SystemExit(main())



if __name__ == "__main__":
    raise SystemExit(main())
