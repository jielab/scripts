#!/usr/bin/env python3
"""Shared input preparation and an outcome-blind small target-validation roster.

No fitted PRSformer, block model, or twin predictor is an input. The historical
directory ``prsformer/cache`` is retained as a storage location only.
"""
from __future__ import annotations

#!/usr/bin/env python3
# Resolve data roots through the shared 0os.sh table.
import sys as _paths_sys
from pathlib import Path as _PathsPath
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.append(str(_paths_root / "0f"))
from workflow_paths import configured_path

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
from grid_paths import GRID, analysis_default, read_config, relocate


def file_record(path):
    from file_state import file_record as record
    return record(path)


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    os.replace(temporary, path)


def identified(frame):
    frame = frame.copy()
    column = next((c for c in ("id", "eid", "IID", "#IID", "ID") if c in frame), None)
    if column is None:
        raise ValueError("Input table requires id/eid/IID")
    if column != "id":
        frame = frame.rename(columns={column: "id"})
    if frame.id.isna().any():
        raise ValueError("Missing participant identity")
    frame["id"] = frame.id.astype(str).str.replace(r"^([0-9]+)\.0$", r"\1", regex=True)
    if frame.id.duplicated().any():
        raise ValueError("Duplicate participant identity")
    return frame


def table(path):
    path = Path(path)
    if path.suffix.lower() == ".rds":
        from twin import read_table
        return identified(read_table(path))
    return identified(pd.read_csv(path, sep="\t", dtype={k: str for k in
        ("id", "eid", "IID", "#IID", "ID", "family_id", "group", "group_id")}, low_memory=False))


def load_cohort(analysis_root, trait, cohort=None, prepared=None, group_col=None):
    root = Path(analysis_root)
    source = Path(cohort) if cohort else root / "tmp/prsformer/cache/data.tsv.gz"
    metadata = {}
    if cohort or (not prepared and source.is_file()):
        frame = table(source)
        meta_path = source.parent / "prepare.json"
        if meta_path.is_file():
            metadata = json.loads(meta_path.read_text())
    else:
        source = Path(prepared) if prepared else root / "tmp/cache/abm/prepared.joblib"
        from twin import load_source
        frame, _, metadata = load_source(source, trait)
        frame = identified(frame)
    y_column = trait if trait in frame else "y"
    if y_column not in frame or "split" not in frame:
        raise ValueError(f"{source}: requires {trait}/y and frozen split columns")
    frame["y"] = pd.to_numeric(frame[y_column], errors="raise")
    frame["split"] = frame.split.astype(str).replace({"training": "train", "testing": "test", "validation": "train"})
    if not set(frame.split) <= {"train", "test", "unused"} or "test" not in set(frame.split):
        raise ValueError("Preserve a pre-existing train/test outer roster")
    group = group_col or next((c for c in ("group_id", "family_id", "group") if c in frame), None)
    if group is None:
        raise ValueError("Family/group IDs are required; provide a verified group column or prepare the grouped cohort")
    if group not in frame or frame[group].isna().any() or frame[group].astype(str).str.strip().eq("").any():
        raise ValueError("Missing family/group identity")
    frame["group_id"] = frame[group].astype(str)
    original_test = set(frame.loc[frame.split.eq("test"), "id"])
    from kinship import active, correct_grid, publication_metadata
    frame, correction_meta = correct_grid(frame, {}, root)
    if original_test != set(frame.loc[frame.split.eq("test"), "id"]):
        raise ValueError("The frozen test roster changed")
    if set(frame.loc[frame.split.eq("test"), "group_id"]) & set(frame.loc[frame.split.eq("train"), "group_id"]):
        raise ValueError("Family/group crosses development and untouched test")
    if np.isinf(frame.y.to_numpy(float)).any():
        raise ValueError("Infinite phenotype")
    if trait == "t2dm" and not set(frame.y.dropna().unique()) <= {0, 1}:
        raise ValueError("T2DM requires binary baseline phenotype 0/1")
    ancestry = next((c for c in ("ancestry", "target") if c in frame), None)
    frame["ancestry"] = frame[ancestry].fillna("UNASSIGNED").astype(str) if ancestry else "UNASSIGNED"
    meta = dict(source=str(source.resolve()), source_file_record=file_record(source),
                grouping_status="declared_group_identity", groups_verified=False,
                endpoint_definition="baseline_t2dm" if trait == "t2dm" else trait,
                original_test_N=len(original_test), **correction_meta)
    if active(root):
        meta.update(publication_metadata(root))
    else:
        evidence_path = root / "tmp/cache/shared_preparation/GROUPS_VERIFIED.json"
        if evidence_path.is_file():
            evidence = json.loads(evidence_path.read_text())
            for filename, digest in evidence["files"].items():
                if not Path(filename).is_file() or file_record(filename) != digest:
                    raise ValueError("Shared family source changed: " + filename)
            group_table = pd.read_csv(evidence["mapping"], sep="\t", dtype=str).set_index("eid")
            if not set(frame.id) <= set(group_table.index) or not np.array_equal(
                    group_table.loc[frame.id, "group"].to_numpy(str), frame.group_id.to_numpy(str)):
                raise ValueError("Shared family identities differ from connected-component evidence")
            meta.update(grouping_status="verified_from_kinship_components", groups_verified=True,
                        group_verification=evidence)
    meta["phenotype_definitions"] = metadata.get("phenotype_definitions", metadata.get("metadata", {}).get("phenotype_definitions", {}))
    return frame.reset_index(drop=True), meta


def small_validation_roster(frame, validation_n=5000, seed=20261102, tune_fraction=.2):
    """Read identities/groups/splits only; no phenotype or predictor affects roles."""
    if validation_n < 4 or not 0 < tune_fraction < 1:
        raise ValueError("validation_n >= 4 and 0 < tune_fraction < 1 required")
    roster = frame[["id", "group_id", "split"]].copy()
    if roster.id.duplicated().any() or roster[["id", "group_id", "split"]].isna().any().any():
        raise ValueError("Malformed cohort identities")
    dev = roster.split.eq("train")
    test = roster.split.eq("test")
    if set(roster.loc[dev, "group_id"]) & set(roster.loc[test, "group_id"]):
        raise ValueError("Family/group crosses development and untouched test")
    sizes = roster.loc[dev].groupby("group_id", sort=True).size().to_dict()
    key = lambda group, tag: hashlib.sha256(f"{seed}:{tag}:{group}".encode()).hexdigest()
    chosen, count = [], 0
    for group in sorted(sizes, key=lambda g: key(g, "small-validation")):
        if count + sizes[group] <= validation_n:
            chosen.append(group)
            count += sizes[group]
    if len(chosen) < 2:
        raise ValueError("Small validation budget must contain at least two whole development groups")
    ordered = sorted(chosen, key=lambda g: key(g, "block-tune"))
    n_tune = max(1, min(len(ordered) - 1, round(tune_fraction * len(ordered))))
    tune_groups = set(ordered[:n_tune])
    chosen = set(chosen)
    roster["role"] = np.where(test, "test", "development_unused")
    roster.loc[roster.group_id.isin(chosen) & dev, "role"] = "target_validation"
    roster.loc[roster.split.eq("unused"), "role"] = "excluded"
    roster["block_role"] = roster.role.replace({"target_validation": "build"})
    roster.loc[roster.group_id.isin(tune_groups) & dev, "block_role"] = "tune_model"
    for column in ("role", "block_role"):
        if roster.groupby("group_id")[column].nunique().max() != 1:
            raise ValueError("A family/group crosses fitting roles")
    return roster


def frozen_roster(frame, root, validation_n=5000, seed=20261102):
    """One shared small-validation budget for every benchmark, with immutable test."""
    from resume import run_lock
    out = Path(root) / "qc/benchmarks/_LOCAL_ONLY"
    out.mkdir(parents=True, exist_ok=True)
    path, receipt = out / "small_validation_roster.tsv.gz", out / "small_validation_roster.json"
    expected = small_validation_roster(frame, validation_n, seed)
    canon = expected.sort_values("id").reset_index(drop=True)
    digest = hashlib.sha256(canon.to_csv(index=False).encode()).hexdigest()
    contract = dict(schema="grid.small_validation.1", validation_n=validation_n, seed=seed,
                    assignment="outcome_blind_sha256_whole_groups", roster_sha256=digest,
                    test_N=int(expected.role.eq("test").sum()),
                    actual_validation_N=int(expected.role.eq("target_validation").sum()))
    with run_lock(out):
        reusable = False
        if os.environ.get("GRID_REPLACE", "TRUE") == "FALSE" and receipt.is_file() and path.is_file():
            try:
                saved = pd.read_csv(path, sep="\t", dtype=str).sort_values("id").reset_index(drop=True)
                reusable = saved.equals(canon.astype(str))
            except (OSError, ValueError):
                pass
        if not reusable:
            temporary = path.with_name(path.name + ".tmp")
            expected.to_csv(temporary, sep="\t", index=False, compression="gzip")
            os.replace(temporary, path)
            atomic_json(receipt, contract)
    return expected, dict(contract, path=str(path))


def preparation_options(args, config=None):
    """Pure genotype/phenotype preparation options, independent of CSx scores."""
    from prsformer import parser as preparation_parser
    config = config or {}
    root = Path(args.analysis_root)
    block_path = getattr(args, "block_config", None) or config.get("block_config", str(GRID / "config/abm.blocks.yml"))
    block = read_config(block_path, root)
    cache = Path(block["prsformer_cache"])
    options = preparation_parser().parse_args(["prepare", "--cache-dir", str(cache)])
    saved = cache / "prepare.json"
    if saved.is_file() and os.environ.get("GRID_REPLACE", "TRUE") == "FALSE":
        # Reuse recorded preparation settings when its JSON is complete.
        try:
            saved_options = json.loads(saved.read_text()).get("options", {})
        except (OSError, ValueError, TypeError):
            saved_options = {}
        for key, value in saved_options.items():
            if hasattr(options, key):
                setattr(options, key, value)
    else:
        for key, value in vars(options).copy().items():
            if isinstance(value, str):
                setattr(options, key, relocate(value, root))
        options.train_fraction, options.validation_fraction = .4, .1
        options.pheno_file = os.environ.get("GRID_PHE_FILE", configured_path("DIR_UKB", "Rdata", "all.rds"))
        options.dir_gen = os.environ.get("GRID_TARGET_DIR", relocate(configured_path('DIR_GEN', 'ukb/37/typ'), root))
        options.ancestry_file = os.environ.get("GRID_ANCESTRY_FILE", options.ancestry_file)
        options.snp_list = os.environ.get("GRID_CSX_SNPINFO", options.snp_list)
        options.remove = os.environ.get("GRID_REMOVE", options.remove)
        options.split_group_file = os.environ.get("GRID_GROUP_FILE", "none")
    for key in ("pheno_file", "ancestry_file", "dir_gen", "remove", "group_col", "chrs",
                "covariates", "height_col", "ldl_col", "t2dm_col"):
        value = getattr(args, key, None)
        if value is not None:
            setattr(options, key, value)
    if getattr(args, "snpinfo", None):
        options.snp_list = args.snpinfo
    if getattr(args, "group_file", None):
        options.split_group_file = args.group_file
    options.traits, options.cache_dir = "height,ldl,t2dm", str(cache)
    options.replace = os.environ.get("GRID_REPLACE", "TRUE") == "TRUE"
    return options, block


def preparation_command(args, config=None):
    options, _ = preparation_options(args, config)
    command = [sys.executable, str(GRID / "f/prsformer.py"), "prepare"]
    for key, value in vars(options).items():
        if key in {"command", "func", "replace", "check"} or value is None:
            continue
        flag = "--" + key.replace("_", "-")
        if isinstance(value, bool):
            if value:
                command.append(flag)
        else:
            command.extend([flag, str(value)])
    return command


def align_projected_cohort(frame, analysis_root, *, frozen=False):
    """Bind shared IDs to usable reference PCs before any family split is frozen."""
    path = Path(os.environ.get("GRID_PCA_FILE") or
                configured_path('DIR_GEN', 'ukb/pca_proj/ukb.discodivas.pca.tsv.gz'))
    if not path.is_file():
        raise FileNotFoundError(f"Shared preparation requires projected PCs: {path}; run 1.pca.sh first")
    pca = table(path)
    pcs = [f"PC{i}" for i in range(1, 11)]
    if not set(pcs) <= set(pca):
        raise ValueError("Shared preparation requires projected PC1-PC10")
    if not np.isfinite(pca[pcs].apply(pd.to_numeric, errors="raise").to_numpy(float)).all():
        raise ValueError("Projected PCs contain nonfinite values; rebuild the PCA input")
    ids = identified(frame).id
    keep = ids.isin(pca.id).to_numpy()
    missing = int((~keep).sum())
    if frozen:
        if missing:
            raise ValueError(f"Frozen shared cohort lacks projected PCs for {missing} participants; "
                             "use a new analysis root to rebuild a common roster")
        return frame
    if not keep.any():
        raise ValueError("No shared participants have projected PCs")
    selected = frame.loc[keep].copy()
    atomic_json(Path(analysis_root) / "tmp/cache/shared_preparation/PCA_COHORT.json", dict(
        schema="grid.shared_pca_cohort.1", source=str(path.resolve()), source_file_record=file_record(path),
        required_pcs=pcs, input_N=len(frame), retained_N=len(selected), missing_pca_N=missing,
        selected_ids_sha256=hashlib.sha256("\n".join(ids[keep]).encode()).hexdigest(),
        rule="intersect available projected PCs before freezing family-based splits",
        uses_score_or_outcome_values=False))
    print(f"[shared inputs] PCA coverage: {len(frame):,} -> {len(selected):,}; "
          f"missing projected PCs: {missing:,}", flush=True)
    return selected


def prepare_frozen_roster(options, analysis_root, relationship_file=None):
    """On a fresh run freeze IDs/families without reading any score or model."""
    import prsformer
    root, cache = Path(analysis_root), Path(options.cache_dir)
    # Restore the shared roster paths on every invocation, including replacement.
    # A completed genotype cache alone must not leave independent split defaults active.
    outer = root / "tmp/cache/abm/split.rds"
    legacy = root / "tmp/abm/grid.split.rds"
    if not outer.is_file() and not legacy.is_file():
        genotype_sources = prsformer.genotype_files(options.dir_gen, prsformer.chromosomes(options.chrs))
        if not options.split_group_file or options.split_group_file == "none":
            source = Path(relationship_file or os.environ.get("GRID_REL_FILE", relocate(configured_path('DIR_GEN', 'ukb/37/typ/ukb_rel.dat'))))
            if not source.is_file():
                raise FileNotFoundError(f"Fresh shared preparation needs UKB relationship data: {source}; "
                    "set --relationship-file/GRID_REL_FILE or provide --group-file (eid,group)")
            edges = pd.read_csv(source, sep=r"\s+", dtype={"ID1":str, "ID2":str})
            if not {"ID1", "ID2", "Kinship"} <= set(edges) or edges[["ID1", "ID2", "Kinship"]].isna().any().any():
                raise ValueError("Relationship file requires nonmissing ID1, ID2, Kinship")
            from kinship import components
            _, _, _, sample_path, mode = genotype_sources[0]
            ids = prsformer.sample_ids(sample_path, mode)
            mapping = Path(analysis_root) / "tmp/cache/shared_preparation/raw_family_components.tsv.gz"
            mapping.parent.mkdir(parents=True, exist_ok=True)
            families = pd.DataFrame(dict(eid=ids, group=components(edges, ids).to_numpy(str)))
            families.to_csv(mapping, sep="\t", index=False)
            atomic_json(mapping.parent / "GROUPS_VERIFIED.json", dict(schema="grid.shared_family_components.1",
                mapping=str(mapping), source=str(source.resolve()),
                component_rule="all supplied UKB relationship edges; no additional threshold",
                split_uses_outcome_values=False, N_genotype_ids=len(families), N_components=int(families.group.nunique()),
                files={str(source.resolve()):file_record(source), str(mapping):file_record(mapping)}))
            options.split_group_file = str(mapping)
        initial = prsformer.make_cohort(options, genotype_sources)
        if "group" not in initial:
            raise ValueError("The declared family file was not carried into the raw cohort")
        initial = align_projected_cohort(initial, root)
        cohort = initial[["eid", "target", "group"]].rename(columns={"target": "ancestry", "group": "family_id"})
        settings = SimpleNamespace(**{k:v for k,v in vars(options).items() if k not in {"func", "command"}})
        settings.cache_dir = str(root / "tmp/cache/abm")
        settings.prsformer_root = str(root)
        settings.out_root = str(root / "tmp/abm")
        settings.split_file = "none"
        # The historical 50/50 outer assignment and its family hash are retained;
        # make_cohort's temporary partitions above do not become published roles.
        prsformer.grid_data_helpers().freeze_split(cohort, settings)
    roster = cache / "roster"
    prsformer.shared_roster(SimpleNamespace(cache_dir=str(root / "tmp/cache/abm"), grid_split=str(legacy),
        roster_dir=str(roster), seed=options.seed, check=False))
    align_projected_cohort(table(roster / "prsformer.split.tsv.gz"), root, frozen=True)
    options.keep = str(roster / "common.keep")
    options.split_file = str(roster / "prsformer.split.tsv.gz")
    options.split_group_file = str(roster / "families.tsv.gz")
    options.require_shared_cohort = True


def prepared_for_trait(root, trait):
    """Per-trait caches prevent concurrent traits changing each other's identity."""
    current = Path(root) / "tmp/cache/shared_features" / trait / "prepared.joblib"
    if current.is_file():
        return current
    legacy = Path(root) / "tmp/cache/abm/prepared.joblib"
    if legacy.is_file():
        import joblib
        saved = joblib.load(legacy)  # trusted local user-produced cache only
        if saved.get("state") == "features_ready" and trait in saved.get("tables", {}):
            return legacy
    return current


def complete_score_features(args, config=None):
    """Join only requested scores onto frozen IDs, retaining missing outcomes.

    Unlike historical raw_tables, this does not intersect traits or create a
    new cohort. It requires no other trait's CSx score and no fitted model.
    Existing complete legacy caches are read-only; new caches are per trait.
    """
    import joblib
    from resume import run_lock
    from twin import load_source
    root = Path(args.analysis_root)
    config = config or {}
    block_path = getattr(args, "block_config", None) or config.get("block_config", GRID / "config/abm.blocks.yml")
    block = read_config(block_path, root)
    cache = Path(block["prsformer_cache"])
    source = cache / "data.tsv.gz"
    preparation = json.loads((cache / "prepare.json").read_text())
    covariate_text = preparation.get("options", {}).get("covariates", "age,sex,PC1,PC2")
    covariates = [] if covariate_text in {"none", ""} else covariate_text.split(",")
    requested = getattr(args, "traits", "height,ldl,t2dm")
    traits = requested.split(",") if isinstance(requested, str) else list(requested)
    results = {}
    for trait in traits:
        target = prepared_for_trait(root, trait)
        if target == root / "tmp/cache/abm/prepared.joblib":
            load_source(target, trait)  # validate every existing source stamp
            results[trait] = target
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        with run_lock(target.parent):
            frame, metadata = load_cohort(root, trait, cohort=source)
            score_path = Path(block["score_root"]) / trait / "1csx.scores.rds"
            scores = table(score_path).set_index("id")
            if not set(frame.id) <= set(scores.index):
                raise ValueError(f"{trait} scores do not cover every fixed cohort ID")
            features = []
            for pop in ("AFR", "EAS", "EUR", "SAS"):
                column = next((c for c in (f"csx.{pop}", f"CSX_{pop}", f"csx_{pop}") if c in scores), None)
                if column is None:
                    raise ValueError(f"Missing {trait} population score {pop}")
                name = f"csx.{pop}"
                frame[name] = pd.to_numeric(scores.loc[frame.id, column], errors="raise").to_numpy(float)
                features.append(name)
            if not np.isfinite(frame[features].to_numpy(float)).all():
                raise ValueError("Nonfinite score in fixed cohort")
            paths = [source, cache / "prepare.json", score_path]
            match_pcs = []
            pca_path = Path(os.environ.get("GRID_PCA_FILE", relocate(configured_path('DIR_GEN', 'ukb/pca_proj/ukb.discodivas.pca.tsv.gz'))))
            if pca_path.is_file():
                pca = table(pca_path).set_index("id")
                if not set(frame.id) <= set(pca.index):
                    raise ValueError("Projected PC table does not cover the fixed cohort")
                for index in range(1, int(config.get("distance_pcs", 10)) + 1):
                    name = f"PC{index}"
                    if name in pca:
                        column = "match_" + name
                        frame[column] = pd.to_numeric(pca.loc[frame.id, name], errors="raise").to_numpy(float)
                        match_pcs.append(column)
                paths.append(pca_path)
            if set(covariates) - set(frame):
                raise ValueError("Covariate absent from raw shared cohort")
            frame["family_id"] = frame.group_id
            frame = frame.rename(columns={"id": "eid"})
            cohort = frame[["eid", "family_id", "split", "ancestry"]].copy()
            evidence = root / "tmp/cache/shared_preparation/GROUPS_VERIFIED.json"
            if evidence.is_file():
                paths.append(evidence)
                paths.extend(Path(p) for p in json.loads(evidence.read_text())["files"])
            activation = root / "local/kinship/ACTIVE.json"
            if activation.is_file():
                paths.append(activation)
            records = []
            for path in dict.fromkeys(paths):
                stat = path.stat()
                records.append(dict(path=str(path.resolve()), size=stat.st_size,
                    mtime_ns=stat.st_mtime_ns, file_record=file_record(path)))
            identity = dict(schema="grid.shared_score_features.1", trait=trait, sources=records,
                cohort_ids_sha256=hashlib.sha256("\n".join(cohort.eid).encode()).hexdigest(),
                implementation_file_record=file_record(__file__))
            receipt_path = target.parent / "COMPLETE.json"
            if target.is_file() and receipt_path.is_file():
                receipt = json.loads(receipt_path.read_text())
                if receipt.get("identity") == identity and receipt.get("file_record") == file_record(target):
                    results[trait] = target
                    continue
            metadata.update(sources=records, labels_kept_missing_until_evaluation=True,
                abm_preparation="score_join_on_immutable_shared_ids_no_training_no_other_trait_scores",
                gwas_overlap_status="unverified", phenotype_definitions=preparation.get("phenotype_definitions", {}))
            value = dict(state="features_ready", cohort=cohort, tables={trait:frame}, metadata=metadata,
                feature_groups={trait:dict(csx=features, covariates=covariates, ancestry=match_pcs)},
                settings=dict(traits=[trait], covariates=",".join(covariates)), scoring_artifacts={})
            temporary = target.with_name(target.name + ".tmp")
            joblib.dump(value, temporary, compress=3)
            load_source(temporary, trait)
            os.replace(temporary, target)
            atomic_json(receipt_path, dict(identity=identity, file_record=file_record(target)))
            results[trait] = target
    return results


def ensure_inputs(args, config=None):
    """Prepare shared raw inputs once; lock excludes concurrent cache producers."""
    from kinship import analysis_guard
    from resume import run_lock
    root = Path(args.analysis_root)
    config = config or read_config(getattr(args, "config", None) or GRID / "config/benchmarks.yml", root)
    command = preparation_command(args, config)
    if getattr(args, "dry_run", False):
        import shlex
        print("[shared inputs; no model training] " + shlex.join(command))
        return
    lock_folder = root / "tmp/cache/shared_preparation"
    lock_folder.mkdir(parents=True, exist_ok=True)
    with analysis_guard(root), run_lock(lock_folder):
        options, block = preparation_options(args, config)
        prepare_frozen_roster(options, root, getattr(args, "relationship_file", None))
        # Directly call the reusable input producer; its train/predict functions
        # are never called. This stage needs no CSx scores on a fresh analysis.
        import prsformer
        from dosage_storage import dosage_lock, prsformer_cache_lock
        # Independent trait jobs may already be training from this cache. A
        # verified reuse needs only read locks, otherwise an unnecessary writer
        # lock would wait for those entire training jobs to finish.
        reusable = False
        cache = Path(options.cache_dir)
        with dosage_lock(block, wait=True), prsformer_cache_lock(cache):
            if (cache / "prepare.json").is_file() and (cache / "genotypes.npy").is_file():
                files = prsformer.genotype_files(options.dir_gen, prsformer.chromosomes(options.chrs))
                signature, _ = prsformer.identity(options, files)
                reusable = not options.replace and prsformer.validate_cache(cache, signature)
        if not reusable:
            with dosage_lock(block, exclusive=True, wait=True), prsformer_cache_lock(cache, exclusive=True):
                prsformer.prepare(options)
        else:
            print(f"[shared inputs] REUSE verified prepared cache: {cache}", flush=True)
        if getattr(args, "complete_features", True):
            # Count ABM may request its existing score feature table after CSx
            # scoring. Benchmark preparation itself does not need that table.
            complete_score_features(args, config)
        atomic_json(lock_folder / "PREPARED.json", dict(schema="grid.shared_inputs.1", status="completed",
            cache=str(Path(options.cache_dir)), fitted_models_required=False, score_features_completed=bool(getattr(args, "complete_features", True)),
            preparation_file_record=file_record(Path(options.cache_dir) / "prepare.json")))
    print("[shared inputs] READY: frozen cohort and genotype cache; no model trained", flush=True)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("module", nargs="?", choices=["prepare", "roster"], default="prepare")
    p.add_argument("--analysis-root", "--output-root", default=analysis_default())
    p.add_argument("--config", default=str(GRID / "config/benchmarks.yml"))
    p.add_argument("--block-config")
    p.add_argument("--traits", "--trait", default="height,ldl,t2dm")
    p.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--validation-n", type=int, default=5000)
    p.add_argument("--seed", type=int, default=20261102)
    p.add_argument("--cohort")
    p.add_argument("--prepared")
    p.add_argument("--complete-features", action="store_true", help="After CSx scoring, also prepare the count-ABM score table")
    for option in ("pheno-file", "ancestry-file", "dir-gen", "remove", "group-col", "chrs", "snpinfo",
                   "covariates", "height-col", "ldl-col", "t2dm-col"):
        p.add_argument("--" + option)
    p.add_argument("--group-file", "--split-group-file", help="Required for a fresh cohort: eid,group family components")
    p.add_argument("--relationship-file", help="UKB ID1,ID2,Kinship; default GRID_REL_FILE or data/gen/ukb/37/typ/ukb_rel.dat")
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args(argv)
    if a.threads < 1:
        p.error("--threads must be positive")
    if a.module == "prepare":
        ensure_inputs(a)
    else:
        for trait in a.traits.split(","):
            frame, _ = load_cohort(a.analysis_root, trait, a.cohort, a.prepared)
            _, receipt = frozen_roster(frame, a.analysis_root, a.validation_n, a.seed)
            print(json.dumps(receipt))



if __name__ == "__main__":
    main()
