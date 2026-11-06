#!/usr/bin/env python3
"""Lightweight stage/trait DAG planning; no numerical imports or input preparation."""
from __future__ import annotations

# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(_paths_root / "0f"))
from workflow_paths import configured_path
import argparse
import json
import math
import os
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "0f"))
from workflow_jobs import Job, load_plan, resources, run_plan, write_plan

ENTRIES = {"grid": "grid.sh", "pca": "1.pca.sh", "csx": "2.csx.sh", "disco": "3.disco.sh", "prsformer": "4.prsformer.sh", "install": "install.sh"}
IDENTITY = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_.-]*$")
COMMON_ENV = {"submit", "SUBMIT", "DATE_SCRATCH", "DIR_WORK", "DIR_DATA", "DIR_SCRATCH", "DIR_SOFTWARE",
              "ANALYSIS_ROOT", "SCRIPTS_ROOT", "OS_FILES_ROOT", "OS_SNP_LIST", "SHARED_PYTHON", "WORKFLOW_PYTHON"}


def _env():
    return {k: v for k, v in os.environ.items()
            if (k in COMMON_ENV or k.startswith(("GRID_", "PRSFORMER_")))
            and not re.search(r"TOKEN|PASSWORD|SECRET|CREDENTIAL", k, re.I)}


def option(tokens, names):
    for i, token in enumerate(tokens):
        if token in names:
            return tokens[i + 1] if i + 1 < len(tokens) else None
        for name in names:
            if token.startswith(name + "="):
                return token.split("=", 1)[1]
    return None


def subset(tokens, options, flags=()):
    """Project-specific prerequisite argv, never pass unrelated model options."""
    result = []
    i = 0
    while i < len(tokens):
        token = tokens[i]
        name, sep, value = token.partition("=")
        if name in flags:
            result.append(token)
        elif name in options:
            if sep:
                result += [name, value]
            elif i + 1 < len(tokens):
                result += [name, tokens[i + 1]]
                i += 1
            else:
                raise ValueError(f"Missing value: {name}")
        i += 1
    return result


def without(tokens, options):
    result, i = [], 0
    while i < len(tokens):
        token = tokens[i]
        if token.split("=", 1)[0] in options:
            i += 1 if "=" in token else 2
        else:
            result.append(token)
            i += 1
    return result


def score_subset(tokens):
    """Use the native scorer's chromosome order and permanent subdirectory."""
    text = option(tokens, {"--chrs", "--chr"}) or os.getenv("GRID_CHRS", "1-22")
    chromosomes = []
    for part in re.split(r"[,;\s]+", text):
        if not part:
            continue
        if re.fullmatch(r"\d+-\d+", part):
            first, last = map(int, part.split("-"))
            chromosomes.extend(range(first, last + 1))
        else:
            part = part.replace("chr", "").replace("CHR", "")
            if not part.isdigit():
                raise ValueError("Invalid autosome list: " + text)
            chromosomes.append(int(part))
    chromosomes = list(dict.fromkeys(chromosomes))
    if not chromosomes or any(c < 1 or c > 22 for c in chromosomes):
        raise ValueError("Invalid autosome list: " + text)
    return "" if chromosomes == list(range(1, 23)) else "chr" + ",".join(map(str, chromosomes))


def _profile(args):
    cfg = json.loads(Path(args.hpc_config).read_text())
    if cfg.get("schema_version") != 1:
        raise ValueError("Unsupported HPC resource profile")
    return cfg


def _resources(args, config, stage):
    r = {**config.get("defaults", {}), **config.get("stages", {}).get(stage, {})}
    local = os.getenv("submit", "local") == "local"
    if local:
        r.update(config.get("local", {"cores": 8, "memory_gib": 16, "gpus": 0}))
        if os.getenv("GRID_MEMORY_CAP_GB"):
            r["memory_gib"] = float(os.environ["GRID_MEMORY_CAP_GB"])
    gpu_capable = stage in {"weights", "score", "blockprepare", "blocknet", "prsformer", "prsformer2", "countprepare", "twins", "install"}
    wants_gpu = gpu_capable and (args.device == "cuda" or (args.device == "auto" and r.get("gpus", 0) > 0))
    if wants_gpu:
        if not local:
            r.update(config.get("gpu", {}))
        else:
            r["gpus"] = 1
        if args.gpu_queue:
            r["queue"] = args.gpu_queue
        if args.gpu_count is not None:
            if args.gpu_count < 1:
                raise ValueError("A CUDA stage requires --gpu-count >= 1")
            r["gpus"] = args.gpu_count
    else:
        r["gpus"] = 0
        if args.queue:
            r["queue"] = args.queue
    overrides = (("gpu_cores", "cores"), ("gpu_memory_gb", "memory_gib"), ("gpu_hours", "hours")) if wants_gpu else (
        ("job_cores", "cores"), ("job_memory_gb", "memory_gib"), ("job_hours", "hours"))
    for attr, key in overrides:
        value = getattr(args, attr, None)
        if value is not None:
            r[key] = value
    if str(r.get("queue")) == "short" and int(r.get("cores", 0)) < 40:
        r["cores"] = 40
    requested_threads = args.gpu_threads if wants_gpu else args.threads
    if requested_threads and requested_threads > int(r["cores"]) and stage not in {"architecture", "final", "evaluate", "targetfit"}:
        raise ValueError(f"Requested threads exceed allocated cores for {stage}")
    if not local and not (args.gpu_queue if wants_gpu else args.queue) and os.environ.get('HPC_CPU_QUEUE'):
        if wants_gpu:
            r['queue'] = os.getenv('GRID_GPU_QUEUE') or config.get('gpu', {}).get('queue') or os.environ['HPC_GPU_QUEUE']
        elif not config.get('stages', {}).get(stage, {}).get('queue'):
            memory = float(r['memory_gib']); cores = int(r['cores'])
            tier = 'HPC_CPU_QUEUE' if memory <= 180 and cores <= 40 else 'HPC_MEDIUM_QUEUE' if memory <= 480 and cores <= 64 else 'HPC_HIGHMEM_QUEUE' if memory <= 960 and cores <= 120 else 'HPC_HUGEMEM_QUEUE'
            r['queue'] = os.environ[tier]
    return resources(r)


RAW_OPTIONS = {"--pheno-file", "--ancestry-file", "--group-col", "--group-file", "--split-group-file", "--relationship-file", "--dir-gen", "--snpinfo", "--chrs", "--remove", "--height-col", "--ldl-col", "--t2dm-col"}
PREP_OPTIONS = {"--cohort", "--prepared", "--covariates", "--config", "--block-config"} | RAW_OPTIONS
FIT_OPTIONS = {"--cohort", "--prepared", "--family-col", "--scores", "--score-file", "--score-dir", "--validation-n", "--covariates", "--seed", "--prs-standardization", "--models", "--block-config"}
BLOCK_OPTIONS = {"--cohort", "--prepared", "--family-col", "--config", "--block-config", "--ld-block-dir", "--ld-block-files", "--blocks", "--score-dir", "--covariates", "--validation-n", "--seed", "--epochs", "--batch-size"}


def count_block_config(tokens):
    """Resolve the same outcome-free producer config used by count workers."""
    path = Path(option(tokens, {"--config"}) or ROOT / "config/count_abm.json").expanduser()
    if path.is_dir():
        path /= "count_abm.json"
    config = json.loads(path.read_text())
    if config.get("schema_version") != "grid.count_abm.config.1":
        raise ValueError("GRID count requires count_abm.json, not a historical fitted-model config")
    explicit = option(tokens, {"--block-config"})
    if explicit:
        return str(Path(explicit).expanduser().resolve())
    block = Path(config["block_config"]).expanduser()
    return str((block if block.is_absolute() else path.parent / block).resolve())


def build_jobs(args, extra):
    cfg = _profile(args)
    traits = args.traits.split(",")
    if not traits or any(not IDENTITY.fullmatch(t) for t in traits) or len(set(traits)) != len(traits):
        raise ValueError("--traits requires unique safe trait names")
    root = str(Path(args.analysis_root).expanduser().resolve())
    if args.entry == "grid" and args.module in {"prepare", "all"} and not option(extra, {"--input-manifest"}):
        extra = without(extra, {"--block-config"}) + ["--block-config", count_block_config(extra)]
    base_env = _env()
    base_env.update(GRID_OUTPUT_ROOT=root, GRID_ANALYSIS_ROOT=root, WORKFLOW_WORKER="1")
    for key, names in {
        "GRID_PHE_FILE": {"--pheno-file", "--phe-file", "--phe"},
        "GRID_TARGET_DIR": {"--dir-gen"}, "GRID_CSX_SNPINFO": {"--snpinfo", "--csx-snpinfo"},
        "GRID_ANCESTRY_FILE": {"--ancestry-file"}, "GRID_GROUP_FILE": {"--group-file", "--split-group-file"},
        "GRID_PCA_FILE": {"--pca-file"},
        "GRID_REL_FILE": {"--relationship-file"}, "GRID_SCORE_DIR": {"--score-dir", "--dir-score"},
        "GRID_GWAS_DIR": {"--dir-gwas"}, "GRID_REMOVE": {"--remove", "--exclude-ids"},
    }.items():
        value = option(extra, names)
        if value is not None:
            base_env[key] = value
    if base_env.get("GRID_PCA_FILE") and not base_env.get("GRID_ANCESTRY_FILE"):
        base_env["GRID_ANCESTRY_FILE"] = str(Path(base_env["GRID_PCA_FILE"]).parent / "ukb.ancestry.auto.tsv.gz")
    local = os.getenv("submit", "local") == "local"
    jobs = []
    def add(key, stage, entry, module, selected, depends=(), passed=None, device=False):
        r = _resources(args, cfg, stage)
        requested_threads = (args.gpu_threads or args.threads or (2 if entry == "prsformer" else 8)) if r["gpus"] else (args.threads or 8)
        stage_threads = min(requested_threads, r["cores"])
        env = dict(base_env)
        env.update(GRID_MEMORY_CAP_GB=str(math.ceil(r["memory_gib"])), GRID_THREADS=str(stage_threads),
                   GRID_JOBS=str(max(1, min(4, r["cores"] // stage_threads))))
        cpu_only = stage not in {"weights", "score", "blockprepare", "blocknet", "prsformer", "prsformer2", "countprepare", "twins", "install"}
        worker_device = "cpu" if cpu_only else (args.device if local else ("cuda" if r["gpus"] else "cpu"))
        env["GRID_CSX_BACKEND"] = worker_device
        env["WORKFLOW_DEVICE"] = worker_device
        command = ["bash", str(ROOT / ENTRIES[entry])]
        if entry not in {"pca", "install"}:
            command.append(module)
        if entry != "install" and module != "self-check":
            command += ["--traits", ",".join(selected)]
            output_flag = "--output-root" if entry in {"csx", "disco", "pca"} and module in {"all", "weights", "score"} else "--analysis-root"
            command += [output_flag, root]
        command += list(extra if passed is None else passed)
        if entry == "csx" and stage in {"weights", "score"}:
            children = option(command, {"--jobs"})
            if children and int(children) * int(env["GRID_THREADS"]) > r["cores"]:
                raise ValueError("CSx --jobs * --threads exceeds --job-cores; reduce within-job parallelism or request more cores")
        if device:
            command += ["--device", worker_device,
                        "--threads", str(stage_threads)]
        # Preserve the public caller's meaning of relative scientific paths.
        jobs.append(Job(key, command, list(depends), r, env, str(Path.cwd())))
        return key

    def bundle(trait, entry, start):
        """One allocation per trait; retain native stage commands and caches."""
        steps = jobs[start:]
        del jobs[start:]
        # Sequential steps need their peak CPU/RAM, not the sum. A GPU bundle
        # must also fit CPU preparation. The largest stage walltime (or the
        # explicit --job-hours override) is the entire bundle's time budget.
        allocation = dict(next((j.resources for j in steps if j.resources["gpus"]), steps[0].resources))
        for key in ("cores", "memory_gib", "hours", "gpus"):
            allocation[key] = max(j.resources[key] for j in steps)
        allocation = resources(allocation)
        env = dict(base_env)
        env.update(GRID_MEMORY_CAP_GB=str(math.ceil(allocation["memory_gib"])),
                   GRID_THREADS=str(min(args.threads or 8, allocation["cores"])))
        env["GRID_TRAIT_STEPS"] = json.dumps([
            {"key": j.key, "argv": j.argv,
             "env": {k: j.env[k] for k in ("GRID_MEMORY_CAP_GB", "GRID_THREADS", "GRID_JOBS",
                                           "GRID_CSX_BACKEND", "WORKFLOW_DEVICE")}}
            for j in steps], separators=(",", ":"))
        output_flag = "--analysis-root" if entry == "prsformer" else "--output-root"
        command = ["bash", str(ROOT / ENTRIES[entry]), "all", "--traits", trait,
                   output_flag, root, *extra]
        jobs.append(Job(f"{trait}.all", command, [], allocation, env, str(Path.cwd())))

    module = args.module
    checking = "--check" in extra and args.entry != "install"
    if checking:
        supported = ((args.entry == "pca" and module == "all")
                     or (args.entry == "csx" and module in {"all", "weights", "score"})
                     or (args.entry == "disco" and module in {"all", "score"})
                     or (args.entry == "prsformer" and module == "prsformer2"))
        if not supported:
            raise ValueError(f"{args.entry} {module} does not support --check; use --plan-only to inspect its DAG")
        # Native check paths validate inputs without preparing the shared cohort
        # or fitting target models. Do not append ordinary downstream workers.
        if args.entry in {"csx", "disco"}:
            stage = "weights" if args.entry == "csx" and module != "score" else "score"
            for trait in traits:
                add(f"{trait}.check", stage, args.entry, module, [trait])
            return jobs
    if args.entry in {"pca", "install"}:
        if module != "all":
            raise ValueError(f"{args.entry} accepts options without a positional module")
        add(args.entry, args.entry, args.entry, "all", traits)
    elif args.entry == "grid":
        if module in {"self-check", "cleanup"}:
            add(module, "final", "grid", module, traits)
            return jobs
        if module in {"assess", "abm"}:
            module = "evaluate"
        if module not in {"prepare", "twins", "evaluate", "architecture", "final", "all"}:
            raise ValueError("GRID fitted block models use 2.csx.sh blockridge/blocknet")
        if module == "architecture":
            add("count.cpu_architecture", "architecture", "grid", "architecture", traits,
                passed=without(extra, RAW_OPTIONS))
            return jobs
        preparation = []
        if module in {"prepare", "all"} and not option(extra, {"--input-manifest"}):
            preparation = [add("shared.cpu_prepare", "prepare", "csx", "prepare", traits,
                               passed=subset(extra, PREP_OPTIONS - {"--config"}))]
            # Shared genotype/block writers remain serialized on CPU nodes.
            for trait in traits:
                preparation = [add(f"{trait}.count.cpu_blockprepare", "countgeometry", "csx", "blockprepare",
                    [trait], preparation, passed=subset(extra, BLOCK_OPTIONS - {"--config"}), device=True)]
        end = []
        for trait in traits:
            parent = preparation
            passed = without(extra, RAW_OPTIONS)
            if module in {"prepare", "all"}:
                parent = [add(f"{trait}.count.cpu_geometry", "countgeometry", "grid", "geometry",
                              [trait], parent, passed=passed, device=True)]
                phase = "radius_twins" if module == "all" else "radius"
                kind = "gpu" if _resources(args, cfg, "countprepare")["gpus"] else "cpu"
                parent = [add(f"{trait}.count.{kind}_{phase}", "countprepare", "grid", phase,
                              [trait], parent, passed=passed, device=True)]
            elif module == "twins":
                kind = "gpu" if _resources(args, cfg, "twins")["gpus"] else "cpu"
                parent = [add(f"{trait}.count.{kind}_twins", "twins", "grid", "twins",
                              [trait], passed=passed, device=True)]
            if module in {"all", "evaluate"}:
                parent = [add(f"{trait}.count.cpu_evaluate", "evaluate", "grid", "evaluate",
                              [trait], parent, passed=passed, device=True)]
            end.extend(parent)
        if module == "all":
            end = [add("count.cpu_architecture", "architecture", "grid", "architecture", traits, end,
                       passed=without(extra, RAW_OPTIONS))]
        if module in {"all", "final"}:
            add("count.cpu_final", "final", "grid", "final", traits, end, passed=without(extra, RAW_OPTIONS))
    elif args.entry == "csx":
        if module not in {"all", "weights", "score", "prepare", "targetfit", "blockprepare", "blockridge", "blocknet"}:
            raise ValueError(f"Unknown CSx module: {module}")
        if module in {"all", "weights", "score"}:
            legacy = [t for t in extra]
            # target-only options belong to targetfit, never to upstream MCMC.
            clean = []
            i = 0
            while i < len(legacy):
                name = legacy[i].split("=", 1)[0]
                if name in {"--cohort", "--prepared", "--family-col", "--validation-n", "--covariates", "--prs-standardization", "--scores", "--score-file", "--block-config"}:
                    i += 1 if "=" in legacy[i] else 2
                else:
                    clean.append(legacy[i]); i += 1
            for trait in traits:
                start = len(jobs)
                parent = []
                if module in {"all", "weights"}:
                    parent = [add(f"{trait}.weights", "weights", "csx", "weights", [trait], passed=clean)]
                if module in {"all", "score"}:
                    parent = [add(f"{trait}.score", "score", "csx", "score", [trait], parent, passed=clean)]
                if module == "all":
                    if not option(extra, {"--cohort", "--prepared"}):
                        add(f"{trait}.prepare", "prepare", "csx", "prepare", [trait],
                            passed=subset(extra, PREP_OPTIONS - {"--score-dir"}))
                    fitted_options = subset(extra, FIT_OPTIONS)
                    chrom_subset = score_subset(extra)
                    if chrom_subset:
                        fitted_options += ["--score-subset", chrom_subset]
                    add(f"{trait}.targetfit", "targetfit", "csx", "targetfit", [trait], passed=fitted_options)
                    bundle(trait, "csx", start)
        elif module == "prepare":
            add("shared.prepare", "prepare", "csx", "prepare", traits, passed=extra)
        elif module == "targetfit":
            for trait in traits:
                add(f"{trait}.targetfit", "targetfit", "csx", "targetfit", [trait], passed=extra)
        else:
            deps = []
            if not option(extra, {"--cohort", "--prepared"}):
                deps = [add("shared.prepare", "prepare", "csx", "prepare", traits, passed=subset(extra, PREP_OPTIONS))]
            for trait in traits:
                deps = [add(f"{trait}.blockprepare", "blockprepare", "csx", "blockprepare", [trait], deps,
                            passed=subset(extra, BLOCK_OPTIONS, {"--replace"}), device=True)]
            if module != "blockprepare":
                for trait in traits:
                    add(f"{trait}.{module}", module, "csx", module, [trait], deps,
                        passed=without(extra, RAW_OPTIONS), device=True)
    elif args.entry == "disco":
        if module not in {"all", "score", "evaluate", "prepare"}:
            raise ValueError("Disco modules: all, score, evaluate, prepare")
        if module == "prepare":
            add("shared.prepare", "prepare", "disco", "prepare", traits, passed=subset(extra, PREP_OPTIONS))
        else:
            for trait in traits:
                start = len(jobs)
                if module == "all" and not option(extra, {"--cohort", "--prepared"}):
                    add(f"{trait}.prepare", "prepare", "disco", "prepare", [trait],
                        passed=subset(extra, PREP_OPTIONS - {"--score-dir"}))
                add(f"{trait}.{module}", "evaluate" if module == "evaluate" else "score", "disco", module, [trait])
                if module == "all":
                    bundle(trait, "disco", start)
    elif args.entry == "prsformer":
        if module not in {"prepare", "prsformer", "prsformer2", "all"}:
            raise ValueError("Unknown PRSformer module")
        if checking:
            for trait in traits:
                add(f"{trait}.{module}.check", module, "prsformer", module, [trait], device=True)
            return jobs
        # Every model invocation verifies/restores shared inputs on CPU first.
        # One producer serves all traits/methods; GPU workers never rebuild it.
        parent = [add("shared.cpu_prepare", "prepare", "prsformer", "prepare", traits,
                      passed=subset(extra, PREP_OPTIONS), device=True)]
        if module != "prepare":
            from uuid import uuid4
            stage_dir = str(Path(root) / "jobs/prsformer_stages" / uuid4().hex)
            passed = without(extra, RAW_OPTIONS)
            for method in (["prsformer", "prsformer2"] if module == "all" else [module]):
                for trait in traits:
                    kind = "gpu" if _resources(args, cfg, method)["gpus"] else "cpu"
                    key = add(f"{trait}.{method}.{kind}_train_predict", method, "prsformer", method,
                              [trait], parent, passed=passed + ["--workflow-stage", "model"], device=True)
                    jobs[-1].env["GRID_PRSFORMER_STAGE_DIR"] = stage_dir
                    add(f"{trait}.{method}.cpu_report", "evaluate", "prsformer", method,
                        [trait], [key], passed=passed + ["--workflow-stage", "report"], device=True)
                    jobs[-1].env["GRID_PRSFORMER_STAGE_DIR"] = stage_dir
    return jobs


def check_upstream_files(jobs):
    """Check external stage outputs before submission, without polling any jobs.

    Outputs of ancestors in this invocation may not exist yet. Their ordinary
    scheduler dependencies still order the workers. Scientific validation stays
    in those workers; this check only stats files and reads small path metadata.
    """
    from grid_paths import read_config

    by_key = {job.key: job for job in jobs}
    parents = {}
    def ancestors(key):
        if key not in parents:
            parents[key] = set(by_key[key].depends)
            for parent in by_key[key].depends:
                parents[key].update(ancestors(parent))
        return parents[key]

    required, produced = [], {}
    def path(job, value):
        value = Path(value).expanduser()
        return (value if value.is_absolute() else Path(job.cwd or Path.cwd()) / value).resolve()

    def need(job, value, hint):
        required.append((job.key, path(job, value), hint))

    def output(job, value):
        produced.setdefault(path(job, value), set()).add(job.key)

    for job in jobs:
        tokens = job.argv[2:]
        entry = Path(job.argv[1]).name
        module = tokens[0] if tokens and not tokens[0].startswith("-") else "all"
        env = {**os.environ, **job.env}
        root = path(job, env["GRID_OUTPUT_ROOT"])
        traits = (option(tokens, {"--traits", "--trait"}) or "height,ldl,t2dm").split(",")
        pca = option(tokens, {"--pca-file"}) or env.get("GRID_PCA_FILE") or configured_path("DIR_GEN", "ukb/pca_proj/ukb.discodivas.pca.tsv.gz")
        ancestry = option(tokens, {"--ancestry-file"}) or env.get("GRID_ANCESTRY_FILE") or str(Path(pca).parent / "ukb.ancestry.auto.tsv.gz")
        score_root = option(tokens, {"--dir-score", "--score-dir"}) or env.get("GRID_SCORE_DIR") or configured_path("DIR_GEN", "ukb/pgs")

        def projection():
            need(job, pca, "bash 1.pca.sh")
            if str(ancestry).lower() not in {"", "none"}:
                need(job, ancestry, "bash 1.pca.sh")

        def block_config():
            source = option(tokens, {"--block-config"})
            if not source:
                config_file = option(tokens, {"--config"})
                if entry == "2.csx.sh" and module.startswith("block"):
                    source = config_file or ROOT / "config/abm.blocks.json"
                else:
                    config = read_config(path(job, config_file or ROOT / "config/abm.json"), root)
                    source = config["block_config"]
            source = path(job, source)
            config = read_config(source, root)
            # read_config normally sees the worker environment. During planning
            # an explicit --dir-score lives in job.env, not os.environ yet.
            from workflow_paths import relocate_defaults
            raw = json.loads(source.read_text())
            if relocate_defaults(raw.get("score_root")) == configured_path("DIR_GEN", "ukb/pgs") and env.get("GRID_SCORE_DIR"):
                config["score_root"] = env["GRID_SCORE_DIR"]
            return config

        # Explicit standalone inputs are external even when a default producer
        # is in the plan. Honor relative paths and per-trait filename templates.
        external_cohort = option(tokens, {"--cohort", "--prepared"})
        for names in ({"--cohort"}, {"--prepared"}, {"--input-manifest"}):
            explicit = option(tokens, names)
            if explicit:
                for trait in traits:
                    need(job, explicit.format(trait=trait), "provide the requested input file")

        if module == "prepare" and entry in {"2.csx.sh", "3.disco.sh", "4.prsformer.sh"}:
            projection()
            cache = Path(block_config()["prsformer_cache"])
            for name in ("data.tsv.gz", "prepare.json", "variants.tsv.gz", "genotypes.npy"):
                output(job, cache / name)

        # Bundled preparation has no separate LSF job, but its external inputs
        # must still be validated before anything is submitted.
        if module == "all" and job.env.get("GRID_TRAIT_STEPS"):
            if entry == "4.prsformer.sh" or not external_cohort:
                projection()
            if entry in {"2.csx.sh", "3.disco.sh"}:
                explicit_scores = option(tokens, {"--scores", "--score-file"})
                if explicit_scores:
                    for trait in traits:
                        need(job, explicit_scores.format(trait=trait), "provide the requested score file")

        if entry == "2.csx.sh" and module in {"weights", "score", "all"}:
            subset_name = score_subset(tokens)
            for trait in traits:
                if module in {"score", "all"}:
                    output(job, Path(score_root) / trait / subset_name / "1csx.scores.rds")
                if module == "score" and not any(k.endswith(".weights") for k in ancestors(job.key)):
                    gwas_root = path(job, option(tokens, {"--dir-gwas"}) or env.get("GRID_GWAS_DIR") or configured_path("DIR_GWAS", "4grid/common"))
                    for pop in ("AFR", "EAS", "EUR", "SAS"):
                        tag = f"{trait.lower()}.{pop}"
                        choices = [gwas_root / tag / "gwas" / (tag + ".gz")]
                        if trait.lower() == "t2dm" and pop == "AFR":
                            choices.append(gwas_root / "t2dm.AFA/gwas/t2dm.AFA.gz")
                        if gwas_root.is_dir():
                            choices.extend(p for p in gwas_root.iterdir() if p.name.lower() == tag.lower() + ".gz")
                        source = next((p for p in choices if p.is_file()), choices[0])
                        weights = str(source).removesuffix(".gz") + ("." + subset_name if subset_name else "") + ".csx.gz"
                        need(job, weights, "bash 2.csx.sh weights")
                        need(job, weights + ".signature", "bash 2.csx.sh weights")

        if entry == "3.disco.sh" and module in {"all", "score"}:
            need(job, pca, "bash 1.pca.sh")
            for trait in traits:
                need(job, Path(score_root) / trait / score_subset(tokens) / "1csx.scores.rds", "bash 2.csx.sh")

        if (entry == "2.csx.sh" and module == "targetfit") or (entry == "3.disco.sh" and module == "evaluate"):
            filename = "1csx.scores.rds" if entry == "2.csx.sh" else "2disco.scores.rds"
            explicit = option(tokens, {"--scores", "--score-file"})
            for trait in traits:
                score = explicit.format(trait=trait) if explicit else Path(score_root) / trait / (option(tokens, {"--score-subset"}) or "") / filename
                need(job, score, "bash " + entry)
            if not external_cohort:
                cache_data = path(job, Path(block_config()["prsformer_cache"]) / "data.tsv.gz")
                if not cache_data.is_file() and "shared.prepare" not in ancestors(job.key):
                    projection()  # This standalone evaluator prepares missing inputs itself.

        if entry == "2.csx.sh" and module == "blockprepare":
            config = block_config()
            for trait in traits:
                score_dir = option(tokens, {"--score-dir"}) or config["score_root"]
                need(job, Path(score_dir) / trait / "1csx.scores.rds", "bash 2.csx.sh")
                provenance = path(job, Path(config["csx_analysis_root"]) / trait / "provenance.json")
                need(job, provenance, "bash 2.csx.sh")
                if provenance.is_file():
                    models = json.loads(provenance.read_text()).get("models", {})
                    for pop in ("AFR", "EAS", "EUR", "SAS"):
                        weights = models.get("csx." + pop, {}).get("weights_file")
                        if weights:
                            need(job, weights, "bash 2.csx.sh weights")

        if entry == "4.prsformer.sh" and module in {"prsformer", "prsformer2"}:
            cache = path(job, block_config()["prsformer_cache"])
            for name in ("prepare.json", "data.tsv.gz", "variants.tsv.gz"):
                need(job, cache / name, "bash 4.prsformer.sh prepare")
            stage = option(tokens, {"--workflow-stage"})
            if stage in {"model", "report"}:
                for trait in traits:
                    handoff = Path(env["GRID_PRSFORMER_STAGE_DIR"]) / f"{module}.{trait}.json"
                    if stage == "model":
                        output(job, handoff)
                    else:
                        need(job, handoff, "complete the preceding GPU model stage")

        if entry == "grid.sh" and module in {"geometry", "radius", "radius_twins", "prepare", "twins", "evaluate", "architecture"}:
            for trait in traits:
                folder = root / "count_abm" / trait
                if module in {"radius", "radius_twins"}:
                    need(job, folder / "GEOMETRY_COMPLETE.json", "bash grid.sh prepare")
                elif module == "twins":
                    need(job, folder / "PREPARE_COMPLETE.json", "bash grid.sh prepare")
                elif module == "evaluate":
                    need(job, folder / "TWINS_COMPLETE.json", "bash grid.sh twins")
                elif module == "architecture":
                    need(job, folder / "PREPARE_COMPLETE.json", "bash grid.sh prepare")
                if module != "architecture":
                    stages = {"radius": ["prepare"], "radius_twins": ["prepare", "twins"]}.get(module, [module])
                    for stage in stages:
                        output(job, folder / (stage.upper() + "_COMPLETE.json"))

    missing = {}
    for key, file, hint in required:
        if not (produced.get(file, set()) & ancestors(key)) and not file.is_file():
            missing.setdefault((file, hint), set()).add(key)
    if missing:
        lines = ["Missing prerequisite files; no jobs were submitted:"]
        for (file, hint), consumers in missing.items():
            lines.extend([f"  [MISSING] {file}", f"    Required by: {', '.join(sorted(consumers))}; first: {hint}"])
        lines.append("Run this command again after the required files are available.")
        raise ValueError("\n".join(lines))


def parser():
    p = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    p.add_argument("--entry", choices=ENTRIES, required=True)
    p.add_argument("--workflow-module", dest="module", default="all", help=argparse.SUPPRESS)
    p.add_argument("--traits", "--trait", default="height,ldl,t2dm")
    p.add_argument("--analysis-root", "--output-root", "--dir-analysis-root", "--dir-output-root",
                   default=os.getenv("GRID_OUTPUT_ROOT", configured_path('ANALYSIS_ROOT', 'grid')))
    p.add_argument("--hpc-config", default=os.getenv("GRID_HPC_CONFIG", str(ROOT / "config/hpc.json")))
    p.add_argument("--job-dir")
    p.add_argument("--resume-plan", type=Path)
    p.add_argument("--queue")
    p.add_argument("--gpu-queue")
    p.add_argument("--job-cores", type=int)
    p.add_argument("--job-memory-gb", type=float)
    p.add_argument("--job-hours", type=float)
    p.add_argument("--gpu-count", type=int)
    p.add_argument("--gpu-cores", type=int)
    p.add_argument("--gpu-memory-gb", type=float)
    p.add_argument("--gpu-hours", type=float)
    p.add_argument("--gpu-threads", type=int)
    p.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    p.add_argument("--threads", type=int)
    p.add_argument("--jobs-N", type=int, default=int(os.getenv("jobs_N", "8")))
    p.add_argument("--wait", nargs="?", const="TRUE", choices=["TRUE", "FALSE"], default=os.getenv("SUBMIT_WAIT", "FALSE"))
    p.add_argument("--plan-only", action="store_true")
    p.add_argument("--dry-run", nargs="?", const="TRUE", choices=["TRUE", "FALSE"], default="FALSE")
    return p


def main(argv=None):
    p = parser()
    tokens = list(sys.argv[1:] if argv is None else argv)
    # The public module is the first project token after --entry ENTRY. Do not
    # let argparse mistake an unknown scientific option's VALUE for a module.
    start = 2 if len(tokens) >= 2 and tokens[0] == "--entry" else 0
    if len(tokens) > start and not tokens[start].startswith("-"):
        module = tokens.pop(start)
        tokens += ["--workflow-module", module]
    args, extra = p.parse_known_args(tokens)
    try:
        mode = os.getenv("submit", "local")
        preview = args.plan_only or args.dry_run == "TRUE"
        if args.entry == "csx":
            legacy_stage = option(extra, {"--stage"})
            if legacy_stage:
                if args.module != "all" and args.module != legacy_stage:
                    raise ValueError("Positional CSx module and --stage disagree")
                args.module = legacy_stage
                extra = without(extra, {"--stage"})
            backend = option(extra, {"--backend"})
            if backend:
                if backend not in {"auto", "cpu", "cuda"} or (args.device != "auto" and args.device != backend):
                    raise ValueError("--backend must be auto|cpu|cuda and agree with --device")
                args.device = backend
                extra = without(extra, {"--backend"})
        if args.resume_plan:
            path = args.resume_plan.resolve()
            if not preview:
                _, plan = load_plan(path)
                check_upstream_files([Job(**job) for job in plan["jobs"]])
        else:
            jobs = build_jobs(args, extra)
            if not preview:
                check_upstream_files(jobs)
            path = write_plan(jobs, args.job_dir or Path(args.analysis_root) / "jobs", mode=mode, name=f"grid-{args.entry}")
        if not preview:
            import shutil
            if not shutil.which(os.getenv("WORKFLOW_PYTHON", "python3")):
                raise ValueError("WORKFLOW_PYTHON is unavailable on this submission host; set the shared worker interpreter before submitting")
        result = run_plan(path, submit=mode, jobs=args.jobs_N, wait=args.wait == "TRUE", dry_run=preview)
        print(json.dumps(result, indent=2))
        return 0
    except (OSError, RuntimeError, ValueError) as exc:
        p.exit(2, f"ERROR: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
