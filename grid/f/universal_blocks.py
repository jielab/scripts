"""One PRS-CSx META coefficient per SNP, on the same BED for every target person."""

import gzip
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import joblib
import numpy as np
import pandas as pd

from benchmark2 import CalibratedModel
from block_models import BlockTransform, fit_glmnet, predict_frozen, feature_schema
from core2 import aligned, atomic_json, atomic_table, group_folds, read_table, stamp
from output_layout import project_path, trait_path

BASE_COLS = [
    "global_score",
    "global_score_squared",
    "age",
    "sex",
    "age_squared",
    "age_sex",
    "score_sex",
    "score_age",
]


def posterior(cfg, root, trait, mode):
    from driver2 import read_manifest_weights
    from infer import phi_key

    fixed = cfg["benchmark"]["input_phi"]
    if str(fixed).lower() == "auto":
        fixed = cfg["inference"]["phi_grid"][0]
    phi = "auto" if mode == "auto" else phi_key(fixed)
    table = read_manifest_weights(root, trait, cfg["chromosomes"], family="csx_meta", phi=phi)
    return phi, table.iloc[0]


def inputs_for(cfg, root, trait, method):
    from ancestry_regions import bed_paths
    from consensus3 import consensus_files
    from genomics import pfile

    mode = method.split(".")[-1]
    code = Path(__file__).parent
    files = [Path(__file__), project_path(root, "tmp/cohort/participants.tsv.gz")]
    settings = {"method": method}
    if ".prepare." in method:
        if cfg["regions"]["layout"] != "consensus":
            raise ValueError("The universal models require a shared GRCh38 BED.")
        phi, model = posterior(cfg, root, trait, mode)
        files.extend([trait_path(root, trait, "inference/weights/manifest.tsv"), *model.weight_files])
        files.extend(consensus_files(bed_paths(cfg)["shared"]))
        for chrom in cfg["chromosomes"]:
            prefix = Path(cfg["inputs"]["genotype_dir"]) / f"chr{chrom}"
            _, pvar = pfile(prefix)
            files.extend([Path(str(prefix) + ".pgen"), Path(str(prefix) + ".psam"), pvar])
        files.extend(
            code / k
            for k in ("ancestry_regions.py", "genomics.py", "scoring_parallel.py", "consensus3.py")
        )
        settings.update(phi=phi, scoring=cfg["scoring"], bed=str(bed_paths(cfg)["shared"]))
    else:
        files.extend(
            [
                trait_path(root, trait, f"universal.prepare.{mode}/complete.json"),
                project_path(root, f"tmp/cohort/{trait}.development.tsv.gz"),
            ]
        )
        files.extend(code / k for k in ("benchmark2.py", "block_models.py", "block_glmnet.R"))
        settings.update(
            block_models=cfg["block_models"],
            folds=cfg["benchmark"]["oof_folds"],
            seed=cfg["cohort"]["seed"],
        )
    return {"settings": settings, "files": [stamp(p) for p in files]}


def region_atlas(x, features, rows):
    """Accumulate TRAIN moments stably; a zero-score block stays in the catalog."""
    n = 0
    mean = np.zeros(x.shape[1])
    m2 = mean.copy()
    cross = mean.copy()
    total_mean = 0.0
    total_m2 = 0.0
    for start in range(0, len(rows), 512):
        data = np.asarray(x[rows[start : start + 512]], dtype="float64")
        if not np.isfinite(data).all():
            raise ValueError("Nonfinite regional score")
        totals = data.sum(axis=1)
        block_mean = data.mean(axis=0)
        batch_total = totals.mean()
        delta = block_mean - mean
        dt = batch_total - total_mean
        size = len(totals)
        factor = n * size / (n + size)
        data -= block_mean
        totals -= batch_total
        m2 += np.einsum("ij,ij->j", data, data) + delta**2 * factor
        cross += totals @ data + dt * delta * factor
        total_m2 += totals @ totals + dt**2 * factor
        mean += delta * size / (n + size)
        total_mean += dt * size / (n + size)
        n += size
    if n < 2 or total_m2 <= 0 or m2.sum() <= 0:
        raise ValueError("The TRAIN roster has insufficient regional variation")
    cols = ["feature_index", "block_id", "chrom38", "start38", "end38", "n_snps"]
    frame = features[cols].copy()
    frame["n_train"] = n
    frame["score_mean"] = mean
    frame["score_sd"] = np.sqrt(m2 / (n - 1))
    frame["variance_share"] = m2 / m2.sum()
    frame["corr_with_sum"] = np.divide(
        cross, np.sqrt(m2 * total_m2), out=np.full_like(m2, np.nan), where=m2 > 0
    )
    return frame


def training_frequencies(cfg, chrom, assigned, out, roster, parallel):
    """Recover the actual PLINK TRAIN allele frequencies used for missing dosages."""
    from genomics import ensure_frequencies, read_roster, write_keep

    serial_path = out / "tmp" / "META" / f"chr{chrom}" / "train.afreq"
    with TemporaryDirectory(dir=out, prefix="freq-") as tmp:
        path = serial_path
        if parallel or not path.is_file():
            # Parallel scoring removes its per-chromosome scratch after merging.
            tmp = Path(tmp)
            keep = tmp / "train.keep"
            write_keep(read_roster(roster), keep, True)
            extract = tmp / "ids.txt"
            assigned.TARGET_ID.to_csv(extract, index=False, header=False)
            path = ensure_frequencies(
                Path(cfg["inputs"]["genotype_dir"]) / f"chr{chrom}",
                keep,
                extract,
                tmp / "train",
                plink=cfg["runtime"]["plink"],
                threads=cfg["runtime"]["threads"],
                memory_mb=cfg["runtime"]["plink_memory_mb"],
            )
        freq = pd.read_csv(path, sep="\t", dtype={"ID": str}).set_index("ID", verify_integrity=True)
    if not assigned.TARGET_ID.isin(freq.index).all():
        raise ValueError(f"Missing TRAIN frequency record on chr{chrom}")
    freq = freq.loc[assigned.TARGET_ID]
    if not np.array_equal(freq.ALT.to_numpy(), assigned.A1.to_numpy()):
        raise ValueError(f"TRAIN frequency allele differs from the scoring allele on chr{chrom}")
    return pd.to_numeric(freq.ALT_FREQS, errors="raise").to_numpy()


def export_weights(cfg, row, trait, phi, features, out, roster):
    """Distribute one harmonized weight file for the shared BED, including gap SNPs."""
    from genomics import assign_blocks, chromosome_weight_sets, index_weight_files

    columns = [
        "CHR",
        "BP38",
        "SNP",
        "TARGET_ID",
        "A1",
        "A2",
        "BETA",
        "OFFSET",
        "TRAIN_ALT_FREQ",
        "block_id",
    ]
    paths = index_weight_files({"META": row.weight_files})
    output = out / "allelic_weights.tsv.gz"
    temporary = out / "allelic_weights.tsv.gz.partial"
    digest = hashlib.sha256()
    parallel = json.loads((out / "metadata.json").read_text()).get("chromosome_processes", 1) > 1
    n_total = n_assigned = 0
    with gzip.open(temporary, "wb", compresslevel=4) as stream:
        header = ("\t".join(columns) + "\n").encode()
        digest.update(header)
        stream.write(header)
        for chrom, blocks in features.groupby("chrom38", sort=True):
            weights = chromosome_weight_sets(paths, chrom)["META"].copy()
            indices = assign_blocks(weights.BP38, blocks)
            assigned = indices >= 0
            counts = np.bincount(indices[assigned], minlength=len(blocks))
            if not np.array_equal(counts, blocks.n_snps.to_numpy()):
                raise ValueError(f"Exported SNP counts disagree with scored blocks on chr{chrom}")
            weights["block_id"] = np.where(
                assigned, blocks.block_id.to_numpy()[indices.clip(min=0)], ""
            )
            weights["TRAIN_ALT_FREQ"] = np.nan
            if assigned.any():
                weights.loc[assigned, "TRAIN_ALT_FREQ"] = training_frequencies(
                    cfg, chrom, weights.loc[assigned], out, roster, parallel
                )
            n_total += len(weights)
            n_assigned += int(assigned.sum())
            for start in range(0, len(weights), 40000):
                chunk = (
                    weights.iloc[start : start + 40000][columns]
                    .to_csv(
                        sep="\t", index=False, header=False, lineterminator="\n", float_format="%.17g"
                    )
                    .encode()
                )
                stream.write(chunk)
                digest.update(chunk)
    temporary.replace(output)
    spec = dict(
        trait=trait,
        phi=phi,
        source=row.score_name,
        build="GRCh38",
        filename=output.name,
        uncompressed_sha256=digest.hexdigest(),
        n_variants=n_total,
        n_in_bed=n_assigned,
        dosage_allele="A1 (target ALT)",
        per_variant_score="dosage * BETA + OFFSET",
        missing_dosage="2 * TRAIN_ALT_FREQ for a covered diploid SNP",
        missing_block_id="outside shared BED; excluded from block scores and covered global sum",
    )
    atomic_json(spec, out / "weight_spec.json")


def prepare(cfg, root, trait, mode, out):
    from ancestry_regions import score

    phi, row = posterior(cfg, root, trait, mode)
    if cfg["regions"]["layout"] != "consensus":
        raise ValueError("Universal scores require a shared GRCh38 BED")
    roster = project_path(root, "tmp/cohort/participants.tsv.gz")
    local = dict(cfg, regions=dict(cfg["regions"], representation="meta"))
    score(local, roster, {"META": row.weight_files}, out, phi)
    people = aligned(read_table(roster, ["id", "role"]), read_table(out / "ids.tsv").id)
    x = np.load(out / "scores.npy", mmap_mode="r")
    catalog = pd.read_csv(out / "features.tsv", sep="\t")
    if x.shape != (len(people), len(catalog)) or not catalog["map"].eq("META").all():
        raise ValueError("Universal score matrix does not match its catalog")
    export_weights(cfg, row, trait, phi, catalog, out, roster)
    atlas = region_atlas(x, catalog, np.flatnonzero(people.role.eq("train")))
    atomic_table(atlas.assign(trait=trait, source=row.score_name, phi=phi), out / "atlas.tsv")
    return [
        "scores.npy",
        "features.tsv",
        "ids.tsv",
        "metadata.json",
        "sum_check.tsv",
        "region_counts.tsv",
        "covered_global.tsv.gz",
        "unassigned_variants.tsv.gz",
        "atlas.tsv",
        "allelic_weights.tsv.gz",
        "weight_spec.json",
    ]


def global_design(people, totals):
    # Common, prespecified global nonlinear terms in both arms. No ancestry labels, PCs or gates.
    age = pd.to_numeric(people.age, errors="raise")
    sex = pd.to_numeric(people.sex, errors="raise")
    z = np.asarray(totals, dtype=float)
    a = (age - 55.0) / 10.0
    return pd.DataFrame(
        dict(
            global_score=z,
            global_score_squared=z * z,
            age=a,
            sex=sex,
            age_squared=a * a,
            age_sex=a * sex,
            score_sex=z * sex,
            score_age=z * a,
        )
    )


def predict(bundle, blocks, ids, covariates, catalog, weight_spec):
    """Predict for new people using age, sex and aligned one-weight block scores."""
    if weight_spec["uncompressed_sha256"] != bundle["weight_spec"]["uncompressed_sha256"]:
        raise ValueError("Query weights differ from the trained model")
    if feature_schema(catalog) != bundle["feature_schema"]:
        raise ValueError("Query BED/catalog differs from the trained model")
    people = aligned(covariates, ids, "query covariates")
    totals = np.empty(len(people))
    if blocks.shape != (len(people), len(catalog)):
        raise ValueError("Query matrix and IDs/catalog disagree")
    for i in range(0, len(people), 512):
        totals[i : i + 512] = np.sum(blocks[i : i + 512], axis=1, dtype=np.float64)
    base = bundle["global_model"].predict(global_design(people, totals))
    pred = predict_frozen(bundle, blocks, base, ids, feature_manifest=catalog)
    pred[bundle["global_name"] + ".prediction"] = base
    return pred


def fit(cfg, root, trait, mode, out):
    name = f"universal.{mode}"
    global_name = f"universal.global.{mode}"
    folder = trait_path(root, trait, f"universal.prepare.{mode}")
    ids = read_table(folder / "ids.tsv").id
    co = aligned(
        read_table(
            project_path(root, "tmp/cohort/participants.tsv.gz"),
            ["id", "age", "sex", "role", "group_id"],
        ),
        ids,
    )
    if co.group_id.isna().any() or co.groupby("group_id").role.nunique().max() != 1:
        raise ValueError("Missing family or family crossing cohort roles")
    truth = aligned(
        read_table(project_path(root, f"tmp/cohort/{trait}.development.tsv.gz")),
        ids,
        allow_missing=True,
    )
    if truth.loc[co.role.eq("test"), "y"].notna().any():
        raise ValueError("Development file contains TEST labels")
    y = np.asarray(
        pd.to_numeric(truth.y.where(co.role.isin(["train", "tune"])), errors="raise"), dtype=float
    )
    tr = np.flatnonzero(co.role.eq("train") & np.isfinite(y))
    tu = np.flatnonzero(co.role.eq("tune") & np.isfinite(y))
    if len(tr) < 20 or len(tu) < 10:
        raise ValueError("Insufficient TRAIN/TUNE data")
    binary = trait == "t2dm"
    if binary and any(set(y[ix]) != {0, 1} for ix in (tr, tu)):
        raise ValueError("Both T2DM classes required")
    x = np.load(folder / "scores.npy", mmap_mode="r")
    catalog = pd.read_csv(folder / "features.tsv", sep="\t")
    totals = np.empty(len(co))
    for i in range(0, len(co), 512):
        totals[i : i + 512] = np.sum(x[i : i + 512], axis=1, dtype=np.float64)
    design = global_design(co, totals)

    def train(rows):
        return CalibratedModel(BASE_COLS, "binary" if binary else "continuous", 1.0).fit(
            design.iloc[rows], y[rows]
        )

    glm = train(tr)
    base = glm.predict(design)
    folds = group_folds(co.group_id.iloc[tr], cfg["benchmark"]["oof_folds"], cfg["cohort"]["seed"])
    oof = np.full(len(co), np.nan)
    for fold in np.unique(folds):
        hold = tr[folds == fold]
        oof[hold] = train(tr[folds != fold]).predict(design.iloc[hold])
    transform = BlockTransform().fit(x, tr)
    model = fit_glmnet(
        x,
        transform,
        y,
        base,
        oof,
        tr,
        tu,
        binary,
        out / "glmnet",
        cfg["block_models"],
        cfg["cohort"]["seed"],
        cfg["runtime"]["rscript"],
        cfg["runtime"]["threads"],
    )
    bundle = dict(
        transform=transform,
        models={name: model},
        binary=binary,
        global_model=glm,
        global_name=global_name,
        weight_spec=json.loads((folder / "weight_spec.json").read_text()),
        feature_schema=feature_schema(catalog),
        genetic_input="one META posterior on shared GRCh38 BED",
    )
    pred = predict_frozen(bundle, x, base, ids, feature_manifest=catalog)
    pred[global_name + ".prediction"] = base
    atomic_table(pred, out / "predictions.tsv.gz")
    joblib.dump(bundle, out / "model.joblib")
    atlas = pd.read_csv(folder / "atlas.tsv", sep="\t")
    atlas["correction_coef_per_score_unit"] = 0.0
    atlas.loc[transform.active, "correction_coef_per_score_unit"] = (
        model["blend"] * model["coef"][1:] / transform.scale
    )
    atomic_table(atlas, out / "block_coefficients.tsv")
    atomic_json(
        dict(
            trait=trait,
            method=name,
            global_control=global_name,
            weight_source=atlas.source.iloc[0],
            n_blocks=len(catalog),
            n_active=len(transform.active),
            global_features=BASE_COLS,
            train_n=len(tr),
            tune_n=len(tu),
            test_labels_used=False,
            **{"fit": model["fit"]},
        ),
        out / "model.json",
    )
    return ["predictions.tsv.gz", "model.joblib", "model.json", "block_coefficients.tsv"]


def predict_cli(argv):
    import argparse
    from core2 import identified

    cli = argparse.ArgumentParser(
        description="Apply a frozen single-weight GRID model to precomputed block scores"
    )
    cli.add_argument(
        "--model",
        required=True,
        type=Path,
        help="Frozen model.joblib from universal.auto or universal.fixed",
    )
    cli.add_argument(
        "--blocks",
        required=True,
        type=Path,
        help="Directory with scores.npy, ids.tsv, features.tsv, weight_spec.json",
    )
    cli.add_argument(
        "--covariates", required=True, type=Path, help="TSV with id, age and the training sex coding"
    )
    cli.add_argument("--out", required=True, type=Path, help="Output prediction TSV or TSV.GZ")
    args = cli.parse_args(argv)
    files = [
        args.blocks / name for name in ("scores.npy", "ids.tsv", "features.tsv", "weight_spec.json")
    ]
    if args.out.resolve() in {p.resolve() for p in [args.model, args.covariates, *files]}:
        cli.error("Choose an output path distinct from the model and inputs")
    scores = predict(
        joblib.load(args.model),
        np.load(files[0], mmap_mode="r", allow_pickle=False),
        identified(read_table(files[1]), "block rows").id,
        read_table(args.covariates, ["id", "age", "sex"]),
        pd.read_csv(files[2], sep="\t"),
        json.loads(files[3].read_text()),
    )
    atomic_table(scores, args.out)
    print(f"Scored {len(scores):,} participants: {args.out}")
