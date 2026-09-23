"""Train / freeze / evaluate / project with an untouched outer half-cohort."""
from pathlib import Path
import copy
import importlib.metadata
import json
import joblib
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, GroupShuffleSplit
from common import words, dump, log, fingerprints, VERSION, digest
from io_data import prepare, read_table, numeric, map_metabolites
from preprocess import outcomes, MolecularPreprocessor, MetadataDesign
from models import oof_quality, tune_logistic, RiskCalibrator
from reference import MolecularGeometry
from survival import CensoringKM, metrics, loss
from neural import token_partition, train_encoder, encode, device_for
from borrowing import ReferenceBank, select_references, tune_bank, reference_utility, SupportRule, stable_seed
from mosaic import ModuleMosaic


PRIMARY = "panome_transformer"


def split_people(p, a):
    if a.split_file:
        table = read_table(a.split_file, a.id_col)
        if "split" not in table or not set(table.split) <= {"build", "tune", "calibration", "test"}:
            raise ValueError("Split file requires build/tune/calibration/test labels")
        part = table.set_index(a.id_col).reindex(p[a.id_col]).split
        if part.isna().any():
            raise ValueError("Split file must cover all eligible participants")
        part = part.to_numpy(str)
    else:
        def divide(ix, fraction, seed):
            if a.group_col:
                one, two = next(GroupShuffleSplit(1, test_size=fraction, random_state=seed).split(
                    ix, groups=p.iloc[ix][a.group_col].to_numpy(str)))
                return ix[one], ix[two]
            labels = (p.iloc[ix].event.eq(1) & p.iloc[ix].time.le(a.horizon)).to_numpy(int)
            return train_test_split(ix, test_size=fraction, random_state=seed, stratify=labels)
        develop, test = divide(np.arange(len(p)), .5, a.seed)
        bt, cal = divide(develop, .2, a.seed+1)
        build, tune = divide(bt, .25, a.seed+2)
        part = np.full(len(p), "test", object)
        part[build], part[tune], part[cal] = "build", "tune", "calibration"
    if set(part) != {"build", "tune", "calibration", "test"}:
        raise ValueError("All four outer partitions are required")
    if a.group_col:
        if p[a.group_col].isna().any() or p[a.group_col].astype(str).str.strip().eq("").any():
            raise ValueError("Complete family component IDs required")
        if pd.DataFrame({"group":p[a.group_col].to_numpy(str), "split":part}).groupby("group").split.nunique().max()>1:
            raise ValueError("Family components cross split boundaries")
    return part


def calibration_halves(p, mask, idcol, groupcol, seed):
    # Outcome-blind, stable under row order. Families stay together.
    unit = p[groupcol if groupcol else idcol].astype(str).to_numpy()
    levels = sorted(set(unit[mask]), key=lambda s: stable_seed(s, seed+832))
    if len(levels)<2:
        raise ValueError("Need two independent calibration units")
    audit_units = set(levels[:len(levels)//2])
    audit = mask & np.array([s in audit_units for s in unit])
    return mask & ~audit, audit


def unknown_rows(p, designs):
    unknown = np.zeros(len(p), bool)
    for design in designs:
        if design is None:
            continue
        norm = design.normalized(p)
        for col, levels in design.levels.items():
            unknown |= (norm[col].notna() & ~norm[col].isin(levels)).to_numpy()
    return unknown


def runtime_manifest(a):
    deps = {}
    for name in ["numpy", "pandas", "scipy", "scikit-learn", "joblib", "lightgbm", "torch", "pyreadr"]:
        try: deps[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError: deps[name] = "not_installed"
    paths = [] if a.demo else [a.phe_file, a.omics_file]
    paths += [s for s in [a.split_file, a.module_file] if s]
    if a.biom == "met" and a.met_input == "raw" and not a.demo:
        paths.append(a.met_map)
    config = {k:v for k,v in vars(a).items() if k not in ["resume", "replace", "train_only"]}
    data = dict(version=VERSION, config=config, dependencies=deps,
        inputs=fingerprints(paths, a.full_input_hash),
        code=fingerprints(list(Path(__file__).parent.glob("*.py")) +
                          [Path(__file__).with_name("export_rds.R"),
                           Path(__file__).parent.parent/"panome.sh"], True))
    data["signature"] = digest(data)
    return data


def train(a, out):
    out = Path(out)
    manifest = runtime_manifest(a)
    if a.resume and (out/"manifest.json").exists():
        old = json.loads((out/"manifest.json").read_text())
        if manifest["signature"] != old["signature"]:
            raise ValueError("Resume requires identical code, data fingerprints, dependencies and configuration")
    dump(out/"manifest.json", manifest)
    prepared = out/"input"; prepared.mkdir(exist_ok=True)
    log("START", "prepare")
    p = prepare(a, prepared)
    raw = np.load(prepared/"raw.npy")
    features = (prepared/"features.txt").read_text().splitlines()
    p, audit = outcomes(p, a)
    eligible = p.eligible.to_numpy(bool)
    raw, p = raw[eligible], p.loc[eligible].reset_index(drop=True)
    p["split"] = split_people(p, a)
    p[[a.id_col, "split"]].to_csv(out/"split_before_qc.csv", index=False)
    build = p.split.eq("build").to_numpy()
    proto = MolecularPreprocessor(a.feature_missing, words(a.residualize), words(a.categorical), a.transform)
    scaled = proto.scale_transform(raw[build])
    keep = (np.mean(~np.isfinite(scaled), axis=0)<a.feature_missing) & (np.nanstd(scaled, axis=0)>1e-8)
    if keep.sum()<3:
        raise ValueError("Fewer than three usable assays")
    missing = np.mean(~np.isfinite(proto.scale_transform(raw[:, keep])), axis=1)
    qc = missing<=a.sample_missing
    pd.DataFrame({a.id_col:p[a.id_col], "split":p.split, "missing_fraction":missing,
                  "included":qc}).to_csv(out/"sample_qc.csv", index=False)
    raw, p, missing = raw[qc], p.loc[qc].reset_index(drop=True), missing[qc]
    build, tune, cal, test = [p.split.eq(s).to_numpy() for s in ["build", "tune", "calibration", "test"]]
    calfit, cala = calibration_halves(p, cal, a.id_col, a.group_col, a.seed)
    if min(build.sum(), tune.sum(), calfit.sum(), cala.sum(), test.sum())<30:
        raise ValueError("Too few people after QC / calibration-audit separation")
    p["role"] = p.split
    p.loc[calfit, "role"], p.loc[cala, "role"] = "calibration_fit", "calibration_audit"
    p[[a.id_col, "split", "role"]].to_csv(out/"split.csv", index=False)
    if a.shuffle_development_outcomes:
        rng = np.random.default_rng(a.seed+912)
        for role in ["build", "tune", "calibration_fit", "calibration_audit"]:
            ix = np.flatnonzero(p.role.eq(role))
            permuted = rng.permutation(ix)
            p.loc[ix,["time","event"]] = p.loc[permuted,["time","event"]].to_numpy()
        log("CONTROL", "development_outcomes", "Joint time/event pairs permuted within each development role; not a permutation p-value")
    audit.update(n_after_qc=len(p), sample_missing_exclusions=int((~qc).sum()),
        retained_features=int(keep.sum()), group_split=bool(a.group_col),
        role_counts=p.role.value_counts().to_dict(), test_Y_not_used_after_initial_stratification=True,
        development_outcomes_shuffled=a.shuffle_development_outcomes)
    dump(out/"cohort_audit.json", audit)
    prep = proto.fit(raw[build], p.loc[build], allowed=keep)
    x, observed = prep.transform(raw, p)
    clinical = MetadataDesign(words(a.covariates), words(a.categorical)).fit(p.loc[build])
    c = clinical.transform(p)
    if c.shape[1]==0:
        raise ValueError("At least one clinical comparator covariate required")
    tech = MetadataDesign(list(dict.fromkeys(words(a.covariates)+words(a.residualize))), words(a.categorical)).fit(p.loc[build])
    unknown = unknown_rows(p, [clinical, getattr(prep, "design", None)])
    km = CensoringKM().fit(p.loc[build], a.horizon, a.min_censor_survival)
    y, w = np.zeros(len(p), int), np.zeros(len(p))
    y[~test], w[~test] = km.labels_weights(p.loc[~test])
    for name, mask in [("build", build), ("tune", tune), ("calibration_fit", calfit), ("calibration_audit", cala)]:
        if np.sum(y[mask])<a.min_events or np.sum((w[mask]>0)&(y[mask]==0))<a.min_events:
            raise ValueError(f"Too few known cases/controls in {name}; change horizon/data or --min-events for a demo")
    ids = p[a.id_col].to_numpy(str)
    groups = p[a.group_col].to_numpy(str) if a.group_col else None
    bg = None if groups is None else groups[build]
    tg = None if groups is None else groups[tune]
    f_names = [features[j] for j in prep.keep]
    log("DONE", "prepare", f"N={len(p)}, assays={x.shape[1]}, build={build.sum()}, test={test.sum()}")
    quality, importance = oof_quality(raw[build], p.loc[build].reset_index(drop=True), features, a, out)
    geometry = MolecularGeometry().fit(x[build], importance[prep.keep], a.dimensions, a.seed)
    pca_z = geometry.unsupervised.transform(x)/geometry.scale_u/np.sqrt(len(geometry.scale_u))
    membership, token_names = token_partition(x[build], f_names, a.tokens, a.seed, a.module_file)
    pd.DataFrame({"feature":f_names, "token":[token_names[j] for j in membership]}).to_csv(out/"token_membership.csv", index=False)
    models, rawpred, model_tuning = {}, {}, []
    blocks = {"clinical":c, "elasticnet":np.c_[c,x], "protein_elasticnet":x,
        "clinical_pca":np.c_[c,pca_z], "clinical_technical":np.c_[tech.transform(p), 1-observed.mean(1)]}
    log("START", "classical_comparators")
    for name, block in blocks.items():
        model, rows = tune_logistic(block, y, w, build, tune, name, a, .5 if "elasticnet" in name else 0)
        models[name] = model; model_tuning.extend(rows)
        rawpred[name] = model.predict_proba(block)[:, 1]
    if a.tree != "none":
        best = None
        for leaves in [7, 15, 31]:
            if a.tree == "lightgbm":
                from lightgbm import LGBMClassifier
                obj = LGBMClassifier(n_estimators=a.tree_estimators, num_leaves=leaves,
                    learning_rate=.03, min_child_samples=40, reg_lambda=10, random_state=a.seed,
                    n_jobs=a.cores, verbosity=-1, deterministic=True, force_col_wise=True)
            else:
                from sklearn.ensemble import HistGradientBoostingClassifier
                obj = HistGradientBoostingClassifier(max_iter=a.tree_estimators, max_leaf_nodes=leaves,
                    learning_rate=.03, min_samples_leaf=40, l2_regularization=10,
                    early_stopping=False, random_state=a.seed)
            known = build & (w>0)
            obj.fit(blocks["elasticnet"][known], y[known], sample_weight=w[known]/w[known].mean())
            pred = obj.predict_proba(blocks["elasticnet"][tune])[:, 1]
            score = float(np.average(loss(y[tune], pred), weights=w[tune]))
            model_tuning.append(dict(model=a.tree, leaves=leaves, validation_logloss=score))
            if best is None or score<best[0]: best = score, obj
        models[a.tree] = best[1]
        rawpred[a.tree] = best[1].predict_proba(blocks["elasticnet"])[:, 1]
    pd.DataFrame(model_tuning).to_csv(out/"comparator_tuning.csv", index=False)
    # Independent models answer whether attention, pretraining, and retrieval help.
    encoders, embeddings = {}, {}
    for name in a.experiments:
        objective = "direct" if name == "direct" else "retrieval"
        architecture = name if name in ["mlp","uniform"] else "transformer"
        log("START", "neural_"+name)
        obj, ssl = train_encoder(x[build], observed[build], y[build], w[build], ids[build], bg,
            membership, a, out/"neural"/name, architecture, objective, pretrain=name!="no_pretrain",
            cross_mode="metric" if name=="metric" else "qkv")
        encoders[name] = obj
        z, direct, _ = encode(obj, x, observed, device_for(a.device), a.batch_size)
        obj.cpu(); embeddings[name] = z
        rawpred[name+"_direct_head"] = direct
        if name == "transformer" and ssl is not None:
            encoders["ssl"] = ssl
            embeddings["ssl"] = encode(ssl, x, observed, device_for(a.device), a.batch_size)[0]
            ssl.cpu()
        log("DONE", "neural_"+name)
    primary_z = embeddings["transformer"]
    embedding_audit = []
    for name, zz in embeddings.items():
        eigen = np.maximum(np.linalg.eigvalsh(np.cov(zz[build].T)), 0)
        embedding_audit.append(dict(model=name, mean_coordinate_SD=float(zz[build].std(0).mean()),
            effective_rank=float(eigen.sum()**2/max(np.sum(eigen**2),1e-30)),
            dimensions=zz.shape[1], warning="representation_collapse" if eigen.sum()<1e-6 else ""))
    pd.DataFrame(embedding_audit).to_csv(out/"embedding_diagnostics.csv",index=False)
    xb, ob, yb, wb, bid = x[build], observed[build], y[build], w[build], ids[build]
    known_indices = np.flatnonzero(wb>0)
    primary_anchor = select_references(primary_z[build], quality, a.panel_size, "reliable", a.seed)
    qualified = len(primary_anchor)==a.panel_size and len(np.unique(yb[primary_anchor]))==2
    fallback = len(primary_anchor)<2
    if fallback:
        primary_anchor = select_references(primary_z[build], quality, a.panel_size, "diversity", a.seed)
    banks, specs, tuning, direct_map = {}, {}, [], {}
    def add_bank(name, encoder_name, anchor, kernel="attention", literal=False, balanced=False,
                 z_override=None, full=False):
        if len(anchor)<2:
            log("SKIP", name, "fewer than 2 reference candidates"); return
        zz = embeddings[encoder_name] if z_override is None else z_override
        encoder = encoders.get(encoder_name) if kernel in ["attention","equal_heads"] else None
        bank = ReferenceBank(zz[build], xb, ob, bid, bg, yb, wb, anchor, encoder,
                             inclusion_correction=balanced, seed=a.seed)
        ks = [1] if literal else a.all_k if full else a.match_k
        ks = sorted(set(min(k, len(anchor)) for k in ks))
        # Scalar temperature candidates matter: unlike v4, the bandwidth does not
        # force every chosen donor into a narrow, nearly uniform weight interval.
        chosen, rows = tune_bank(bank, zz[tune], ids[tune], tg, y[tune], w[tune], ks,
            [1.] if literal else a.temperatures, [0.] if literal else a.prior_grid,
            name, kernel, literal)
        specs[name] = dict(encoder=encoder_name, **chosen,
            space="full_proteome" if z_override is x else "pca" if z_override is pca_z else "neural")
        banks[name] = bank; tuning.extend(rows)
        log("TUNED",name,f"references={len(bank.ids)}; k={chosen['k']}; temperature={chosen['temperature']}; prior={chosen['strength']}")
        og = None if groups is None else groups[~build]
        pred, _ = bank.match(zz[~build], ids[~build], og, **chosen)
        rawpred[name] = np.full(len(p), np.nan); rawpred[name][~build] = pred
        pd.DataFrame({a.id_col:bank.ids, "horizon_label":bank.y,
            "fit_gain":quality.iloc[anchor].fit_gain.to_numpy(),
            "OOF_probability":quality.iloc[anchor].OOF_probability.to_numpy(),
            "class_sampling_correction":bank.correction}).to_csv(out/f"panel_{name}.csv", index=False)
    log("START", "reference_panels")
    add_bank(PRIMARY, "transformer", primary_anchor)
    for mode in ["topfit", "stratified_fit", "random", "diversity"]:
        anchor = select_references(primary_z[build], quality, a.panel_size, mode, a.seed)
        if mode in ["topfit", "stratified_fit"]:
            add_bank("copy1_"+mode, "transformer", anchor, literal=True)
        add_bank("transformer_"+mode+"_panel", "transformer", anchor)
    topfit = select_references(primary_z[build], quality, a.panel_size, "topfit", a.seed)
    add_bank("copy1_fullproteome", "none", topfit, kernel="euclidean", literal=True, z_override=x)
    add_bank("copy1_reliable_panel", "transformer", primary_anchor, literal=True)
    balanced = select_references(primary_z[build], quality, a.panel_size, "reliable", a.seed, balanced=True)
    add_bank("transformer_balanced_panel", "transformer", balanced, balanced=True)
    for size in a.panel_sizes:
        if size != a.panel_size:
            add_bank(f"transformer_reliable{size}", "transformer",
                select_references(primary_z[build], quality, size, "reliable", a.seed))
    add_bank("transformer_fullbank", "transformer", known_indices, full=True)
    full = banks["transformer_fullbank"]
    utility, utility_table = reference_utility(full, primary_z[tune], ids[tune], tg, y[tune], w[tune])
    utility_table.to_csv(out/"tuning_reference_utility.csv", index=False)
    all_utility = np.zeros(build.sum()); all_utility[full.indices] = utility
    utility_anchor = select_references(primary_z[build], quality, a.panel_size, "utility", a.seed, all_utility)
    add_bank("transformer_utility_panel", "transformer", utility_anchor)
    add_bank("equal_heads_same_panel", "transformer", primary_anchor, kernel="equal_heads")
    add_bank("embedding_knn_same_panel", "transformer", primary_anchor, kernel="euclidean")
    add_bank("pca_same_panel", "none", primary_anchor, kernel="euclidean", z_override=pca_z)
    add_bank("pca_fullbank", "none", known_indices, kernel="euclidean", z_override=pca_z, full=True)
    for name in ["mlp", "no_pretrain", "ssl", "uniform", "metric"]:
        if name in encoders:
            add_bank(name+"_same_panel", name, primary_anchor, kernel="euclidean" if name=="ssl" else "attention")
            if name == "mlp":
                add_bank("mlp_fullbank", name, known_indices, full=True)
    primary = banks[PRIMARY]; primary_spec = specs[PRIMARY]
    match_kw = {k:v for k,v in primary_spec.items() if k not in ["encoder", "space"]}
    og = None if groups is None else groups[~build]
    for name, change in [("random_donors_same_panel", {"random_candidates":True}),
                         ("permuted_values_same_panel", {"permuted":True})]:
        specs[name] = dict(primary_spec, **change); banks[name] = primary
        rawpred[name] = np.full(len(p), np.nan)
        rawpred[name][~build] = primary.match(primary_z[~build], ids[~build], og, **match_kw, **change)[0]
    primary_raw, detail = primary.match(primary_z[~build], ids[~build], og, **match_kw)
    info = primary.describe(x[~build], observed[~build], primary_raw, detail)
    info.insert(0, a.id_col, ids[~build]); info["role"] = p.role.to_numpy()[~build]
    log("START", "module_mosaic")
    mosaic = ModuleMosaic().fit(xb, ob, bid, bg, yb, wb, primary_anchor, membership, token_names, a.seed)
    module_risk, module_reference = mosaic.components(x[~build], ids[~build], og)
    mosaic.tune(module_risk[tune[~build]], y[tune], w[tune])
    for name, risk in [("mosaic_equal", module_risk.mean(1)), ("mosaic_weighted", module_risk@mosaic.weights)]:
        rawpred[name] = np.full(len(p), np.nan); rawpred[name][~build] = risk
    pd.DataFrame({"token":token_names, "mixture_weight":mosaic.weights}).to_csv(out/"mosaic_weights.csv", index=False)
    pd.DataFrame(tuning).to_csv(out/"panel_tuning.csv", index=False)
    dump(out/"panel_specs.json", specs)
    # Tune set selects candidate hyperparameters. Calibration and auditing are disjoint.
    calibrators, predictions = {}, {}
    for name, risk in rawpred.items():
        if name.startswith("copy1_"):
            predictions[name] = risk.copy()
        else:
            calibrators[name] = RiskCalibrator().fit(risk[calfit], y[calfit], w[calfit])
            predictions[name] = calibrators[name].predict(risk)
    combinations = {"panome_clinical":(PRIMARY, "clinical"),
                    "panome_elasticnet_hybrid":(PRIMARY, "elasticnet"),
                    "mosaic_clinical":("mosaic_weighted", "clinical")}
    for name, (source, covariate) in combinations.items():
        obj = RiskCalibrator().fit(rawpred[source][calfit], y[calfit], w[calfit], rawpred[covariate][calfit])
        calibrators[name] = obj; predictions[name] = obj.predict(rawpred[source], rawpred[covariate])
    dump(out/"calibration_parameters.json", {name:dict(coefficients=obj.coef,
         logit_mean=obj.mean, logit_scale=obj.sd) for name,obj in calibrators.items()})
    rule = SupportRule().fit(info.loc[tune[~build]], predictions[PRIMARY][tune], y[tune], w[tune],
                             a.accept_quantile, a.min_match_ess)
    supported, reasons = rule.apply(info, missing[~build], unknown[~build], a.sample_missing)
    info["supported_match"], info["rejection_reason"] = supported, reasons
    info["estimated_squared_error"] = rule.error_score(info)
    info["missing_fraction"], info["unknown_category"] = missing[~build], unknown[~build]
    dev_rows = []
    for subset, mask in [("tune", tune), ("calibration_audit", cala)]:
        for name, risk in predictions.items():
            dev_rows.append(dict(subset=subset, model=name, **metrics(y[mask], w[mask], risk[mask])))
    dev_rows = pd.DataFrame(dev_rows); dev_rows.to_csv(out/"development_metrics.csv", index=False)
    primary_met = metrics(y[cala], w[cala], predictions[PRIMARY][cala])
    constant = float(np.average(y[calfit], weights=w[calfit]))
    null = metrics(y[cala], w[cala], np.full(cala.sum(), constant))
    ready_reasons = []
    if not qualified: ready_reasons.append("requested_reliable_panel_incomplete_or_one_class")
    if supported[cala[~build]].mean()<a.min_coverage: ready_reasons.append("low_audit_coverage")
    if primary_met["Brier_IPCW"]>=null["Brier_IPCW"]: ready_reasons.append("no_independent_audit_Brier_gain")
    if primary_met["LogLoss_IPCW"]>=null["LogLoss_IPCW"]: ready_reasons.append("no_independent_audit_logloss_gain")
    if not np.isfinite(primary_met["AUC_IPCW"]) or primary_met["AUC_IPCW"]<=.5:
        ready_reasons.append("no_independent_audit_discrimination")
    readiness = dict(status="provisional_pass" if not ready_reasons else "not_ready", reasons=ready_reasons,
        reliable_candidates=int(quality.reliable_candidate.sum()), requested_size=a.panel_size,
        actual_size=len(primary.ids), primary_model=PRIMARY, fallback_diversity=fallback,
        audit_metrics=primary_met, audit_null=null, audit_coverage=float(supported[cala[~build]].mean()),
        interpretation="independent audit point-estimate screen; not an individual accuracy guarantee or cohort validity test")
    dump(out/"reference_readiness.json", readiness)
    dump(out/"support_thresholds.json", dict(thresholds=rule.thresholds, min_ess=rule.min_ess,
        error_threshold=rule.error_threshold, source="tune; error model uses tune Y; audit is separate"))
    for name, value in predictions.items():
        info[name] = value[~build]
    raw_frame = pd.DataFrame({name+"_raw":value[~build] for name,value in rawpred.items()})
    info = pd.concat([info.reset_index(drop=True),raw_frame],axis=1)
    info["raw_borrowed_risk"] = primary_raw
    info["reference_ready"] = readiness["status"]=="provisional_pass"
    info["prediction_released"] = info.supported_match & info.reference_ready
    info["released_net_risk"] = np.where(info.prediction_released, info[PRIMARY], np.nan)
    info.to_csv(out/"development_and_test_individuals.csv.gz", index=False, compression="gzip")
    coverage_rules = {}
    for quantile in [.5,.75,.9,.95,.99,1.]:
        alternate = copy.copy(rule)
        alternate.thresholds = {key:float(info.loc[tune[~build],key].quantile(quantile)) for key in rule.thresholds}
        alternate.error_threshold = float(np.quantile(rule.error_score(info.loc[tune[~build]]),quantile))
        coverage_rules[quantile] = alternate
    profile_build = np.column_stack([xb[:,membership==j].mean(1) for j in range(len(token_names))])
    profile_center, profile_scale = profile_build.mean(0), np.maximum(profile_build.std(0),1e-6)
    bundle = dict(version=VERSION, config=vars(a), features=features, retained_features=f_names,
        prep=prep, clinical=clinical, technical=tech, geometry=geometry, models=models,
        encoders=encoders, banks=banks, specs=specs, calibrators=calibrators, combinations=combinations,
        rule=rule, readiness=readiness, censoring=km, membership=membership, token_names=token_names,
        mosaic=mosaic, constant_risk=constant, primary=PRIMARY, coverage_rules=coverage_rules,
        profile_center=profile_center, profile_scale=profile_scale)
    joblib.dump(bundle, out/"model_bundle.joblib", compress=3)
    dump(out/"MODEL_FROZEN.json", dict(version=VERSION, primary=PRIMARY,
        test_Y_used_for_fitting_or_selection=False,
        test_Y_used_for_initial_stratification=not bool(a.split_file),
        independent_calibration_audit=True, artifact=fingerprints([out/"model_bundle.joblib"], True)))
    testrel = test[~build]
    info.loc[testrel].reset_index(drop=True).to_csv(out/"test_individuals.csv", index=False)
    p.loc[test, [a.id_col, "time", "event"]+([a.group_col] if a.group_col else [])].to_csv(out/"test_outcomes.csv", index=False)
    p.loc[test, [a.id_col]+[v for v in ["age", "sex", "center"] if v in p]].to_csv(out/"test_demographics.csv", index=False)
    dump(out/"prediction_columns.json", list(predictions))
    # Assay means within train-defined tokens describe profiles; they are not pathways.
    profile = (np.column_stack([x[test][:, membership==j].mean(1) for j in range(len(token_names))])-profile_center)/profile_scale
    profiles = pd.DataFrame(profile, columns=token_names); profiles.insert(0, a.id_col, ids[test])
    profiles.to_csv(out/"molecular_profiles.csv", index=False)
    from evidence import write_evidence, masked_validation, same_risk_pairs, module_perturbations
    from attention_audit import audit_attention
    testdetail = {k:v[testrel] for k,v in detail.items()}
    write_evidence(out, a, ids[test], x[test], observed[test], primary, testdetail,
        predictions[PRIMARY][test], f_names, membership, token_names, calibrators[PRIMARY], match_kw["strength"])
    same_risk_pairs(ids[test], predictions["elasticnet"][test], profile, token_names,
                    a.pair_caliper, a.max_pairs).to_csv(out/"same_risk_pairs.csv", index=False)
    pd.DataFrame({a.id_col:np.repeat(ids[test], len(token_names)), "token":np.tile(token_names, test.sum()),
        "reference_eid":module_reference[testrel].ravel(), "raw_module_risk":module_risk[testrel].ravel(),
        "mixture_weight":np.tile(mosaic.weights, test.sum())}).to_csv(out/"test_mosaic_matches.csv.gz", index=False, compression="gzip")
    if a.explanation_samples:
        log("START", "masked_validation")
        masked_validation(bundle, x[test], observed[test], ids[test], None if groups is None else groups[test], out, a)
        module_perturbations(bundle, x[test], observed[test], ids[test], None if groups is None else groups[test], out, a)
    if a.attention_samples:
        log("START", "attention_audit")
        audit_attention(bundle,x[test],observed[test],ids[test],None if groups is None else groups[test],out,a)
    log("DONE", "model_frozen", readiness["status"])
    return bundle


def predict_bundle(bundle, x, observed, p, device="cpu", details=True):
    cfg = bundle["config"]
    ids = p[cfg["id_col"]].to_numpy(str)
    groups = p[cfg["group_col"]].to_numpy(str) if cfg["group_col"] else None
    if groups is not None and p[cfg["group_col"]].isna().any():
        raise ValueError("Missing query family component")
    c = bundle["clinical"].transform(p)
    geometry = bundle["geometry"]
    pca = geometry.unsupervised.transform(x)/geometry.scale_u/np.sqrt(len(geometry.scale_u))
    blocks = {"clinical":c, "elasticnet":np.c_[c,x], "protein_elasticnet":x,
        "clinical_pca":np.c_[c,pca], "clinical_technical":np.c_[bundle["technical"].transform(p), 1-observed.mean(1)]}
    raw = {name:obj.predict_proba(blocks.get(name, blocks["elasticnet"]))[:, 1]
           for name,obj in bundle["models"].items()}
    embeddings = {}
    for name,obj in bundle["encoders"].items():
        embeddings[name], head, _ = encode(obj, x, observed, device, cfg["batch_size"])
        obj.cpu()
        if name!="ssl": raw[name+"_direct_head"] = head
    detail = None
    for name, bank in bundle["banks"].items():
        spec = bundle["specs"][name]
        z = x if spec["space"]=="full_proteome" else pca if spec["space"]=="pca" else embeddings[spec["encoder"]]
        kw = {k:v for k,v in spec.items() if k not in ["encoder", "space"]}
        raw[name], result = bank.match(z, ids, groups, **kw)
        if name==PRIMARY: detail = result
    component, ref = bundle["mosaic"].components(x, ids, groups)
    raw["mosaic_equal"], raw["mosaic_weighted"] = component.mean(1), component@bundle["mosaic"].weights
    predictions = {name:bundle["calibrators"][name].predict(v) if name in bundle["calibrators"] else v
                   for name,v in raw.items()}
    for name,(source,covariate) in bundle["combinations"].items():
        predictions[name] = bundle["calibrators"][name].predict(raw[source], raw[covariate])
    info = bundle["banks"][PRIMARY].describe(x, observed, raw[PRIMARY], detail)
    unknown = unknown_rows(p, [bundle["clinical"], getattr(bundle["prep"], "design", None)])
    supported, reason = bundle["rule"].apply(info, 1-observed.mean(1), unknown, cfg["sample_missing"])
    info.insert(0, cfg["id_col"], ids)
    for name,v in predictions.items(): info[name] = v
    info["supported_match"], info["rejection_reason"] = supported, reason
    info["reference_ready"] = bundle["readiness"]["status"]=="provisional_pass"
    info["prediction_released"] = supported & info.reference_ready
    info["released_net_risk"] = np.where(info.prediction_released, info[PRIMARY], np.nan)
    info["raw_borrowed_risk"] = raw[PRIMARY]
    info["estimated_squared_error"] = bundle["rule"].error_score(info)
    info["missing_fraction"], info["unknown_category"] = 1-observed.mean(1), unknown
    return info, detail


def project(run_dir, phe_path, omics_path, output, r_bin="Rscript", met_input="named", device="cpu"):
    bundle = joblib.load(Path(run_dir)/"model_bundle.joblib")
    cfg = bundle["config"]; idcol = cfg["id_col"]
    cols = list(dict.fromkeys([idcol]+words(cfg["covariates"])+words(cfg["residualize"])+([cfg["group_col"]] if cfg["group_col"] else [])))
    p = read_table(phe_path, idcol, cols, r_bin)
    omics = read_table(omics_path, idcol, r_bin=r_bin)
    mapping_audit = None
    if cfg["biom"]=="prot": omics.columns = [c if c==idcol else str(c).upper() for c in omics.columns]
    elif met_input=="raw":
        omics, mapping_audit = map_metabolites(omics, cfg["met_map"], idcol,
                                               nonnegative=cfg["transform"] == "log1p")
    if omics.columns.duplicated().any(): raise ValueError("Assay names collide after normalization")
    if not set(omics[idcol])<=set(p[idcol]): raise ValueError("Query IDs lack baseline metadata")
    p = p.set_index(idcol).loc[omics[idcol]].reset_index()
    # Entirely absent assays can be imputed, but their missingness affects support.
    absent = set(bundle["features"])-set(omics)
    omics = omics.reindex(columns=[idcol]+bundle["features"])
    raw = numeric(omics, bundle["features"], "projection")
    x, observed = bundle["prep"].transform(raw, p)
    info, detail = predict_bundle(bundle, x, observed, p, device_for(device))
    output = Path(output); output.parent.mkdir(parents=True, exist_ok=True)
    info.to_csv(output, index=False)
    if mapping_audit is not None:
        mapping_audit.to_csv(output.with_name(output.stem+"_metabolite_mapping.csv"), index=False)
    from evidence import write_match_table
    bank = bundle["banks"][PRIMARY]
    write_match_table(output.with_name(output.stem+"_references.csv.gz"), cfg["id_col"],
                      p[idcol].to_numpy(str), bank, detail)
    dump(output.with_name(output.stem+"_audit.json"), dict(absent_assays=sorted(absent),
        n=len(p), outcome_columns_required=False, query_rows_attend_only_to_build_references=True,
        released=int(info.prediction_released.sum())))
    log("DONE", "project", f"N={len(p)}; released={info.prediction_released.sum()}")


def evaluate(out):
    from evaluation import evaluate_frozen
    return evaluate_frozen(Path(out))
