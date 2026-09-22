"""An outer 50% test set and disjoint build/tune/calibration development sets."""
from pathlib import Path
import importlib.metadata
import json
import joblib
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, GroupShuffleSplit
from sklearn.ensemble import HistGradientBoostingClassifier
from common import words, dump, log, fingerprints, VERSION
from io_data import prepare, read_table, numeric, map_metabolites
from preprocess import outcomes, MolecularPreprocessor, MetadataDesign
from models import oof_quality, tune_logistic, RiskCalibrator
from reference import MolecularGeometry, ReferencePanel, AcceptanceRule, choose_panel, make_profiles
from survival import CensoringKM, metrics, loss, bootstrap_contrasts
from report import matched_pairs, figures, write_report
from explain import MosaicReference, individual_evidence, masked_copy_stability


def split_people(p, a):
    if a.split_file:
        s = read_table(a.split_file, a.id_col)
        if "split" not in s or not set(s.split) <= {"build", "tune", "calibration", "test"}:
            raise ValueError("Split file requires build/tune/calibration/test labels")
        part = s.set_index(a.id_col).reindex(p[a.id_col]).split
        if part.isna().any():
            raise ValueError("Split file does not cover every eligible ID")
        part = part.to_numpy(str)
    else:
        idx = np.arange(len(p))
        y = (p.event.eq(1) & p.time.le(a.horizon)).to_numpy(int)
        def divide(ix, fraction, seed):
            if a.group_col:
                groups = p.iloc[ix][a.group_col].to_numpy(str)
                first, second = next(GroupShuffleSplit(1, test_size=fraction, random_state=seed).split(ix, groups=groups))
                return ix[first], ix[second]
            return train_test_split(ix, test_size=fraction, random_state=seed, stratify=y[ix])
        development, test = divide(idx, .5, a.seed)
        build_tune, calibration = divide(development, .2, a.seed+1)
        build, tune = divide(build_tune, .25, a.seed+2)
        part = np.full(len(p), "test", object)
        part[build], part[tune], part[calibration] = "build", "tune", "calibration"
    if a.group_col:
        if p[a.group_col].isna().any() or p[a.group_col].astype(str).str.strip().eq("").any():
            raise ValueError("Family components must be complete before splitting")
        if pd.DataFrame({"group": p[a.group_col].to_numpy(str), "split": part}).groupby("group").split.nunique().max() > 1:
            raise ValueError("Related participants cross split boundaries")
    if set(part) != {"build", "tune", "calibration", "test"}:
        raise ValueError("All four partitions must be nonempty")
    return part


def unknown_rows(p, designs):
    unknown = np.zeros(len(p), bool)
    for design in designs:
        if design is not None:
            norm = design.normalized(p)
            for col, levels in design.levels.items():
                unknown |= (norm[col].notna() & ~norm[col].isin(levels)).to_numpy()
    return unknown


def runtime_manifest(a):
    deps = {}
    for name in ["numpy", "pandas", "scipy", "scikit-learn", "joblib", "lightgbm"]:
        try:
            deps[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            deps[name] = "not_installed"
    inputs = [] if a.demo else [a.phe_file, a.omics_file]
    inputs += [s for s in [a.split_file, a.module_file] if s]
    if a.biom == "met" and a.met_input == "raw" and not a.demo:
        inputs.append(a.met_map)
    return dict(version=VERSION, config=vars(a), dependencies=deps,
                inputs=fingerprints(inputs, a.full_input_hash),
                code=fingerprints(list(Path(__file__).parent.glob("*.py")) +
                                  [Path(__file__).with_name("export_rds.R"),
                                   Path(__file__).parent.parent/"panome.sh"], True))


def train(a, out):
    dump(out/"manifest.json", runtime_manifest(a))
    log("START", "prepare")
    prepared = out/"input"
    prepared.mkdir(exist_ok=True)
    # Keep original metadata types/values; CSV is an audit export, not model input.
    p = prepare(a, prepared)
    raw = np.load(prepared/"raw.npy")
    features = (prepared/"features.txt").read_text().splitlines()
    p, audit = outcomes(p, a)
    eligible = p.eligible.to_numpy(bool)
    raw, p = raw[eligible], p.loc[eligible].reset_index(drop=True)
    p["split"] = split_people(p, a)
    p[[a.id_col, "split"]].to_csv(out/"split_before_qc.csv", index=False)
    build = p.split.eq("build").to_numpy()
    initial = MolecularPreprocessor(a.feature_missing, words(a.residualize), words(a.categorical), a.transform)
    transformed = initial.scale_transform(raw[build])
    keep = (np.mean(~np.isfinite(transformed), axis=0) < a.feature_missing) & (np.nanstd(transformed, axis=0) > 1e-8)
    if keep.sum() < 3:
        raise ValueError("Fewer than three usable assays")
    miss = np.mean(~np.isfinite(raw[:, keep]), axis=1)
    qc = miss <= a.sample_missing
    exclusion = p[[a.id_col, "split", "event", "time"]].copy()
    exclusion["missing_fraction"], exclusion["included"] = miss, qc
    exclusion.to_csv(out/"sample_qc.csv", index=False)
    raw, p, miss = raw[qc], p.loc[qc].reset_index(drop=True), miss[qc]
    part = p.split.to_numpy()
    build, tune, cal, test = [part == s for s in ["build", "tune", "calibration", "test"]]
    if min(build.sum(), tune.sum(), cal.sum(), test.sum()) < 30:
        raise ValueError("Partitions too small after sample QC")
    audit.update(n_after_qc=len(p), sample_missing_exclusions=int((~qc).sum()),
                 retained_features=int(keep.sum()), split_summary={s: dict(n=int((part==s).sum()),
                 events=int(p.loc[part==s, "event"].sum())) for s in np.unique(part)},
                 group_split=bool(a.group_col))
    dump(out/"cohort_audit.json", audit)
    p[[a.id_col, "split"]].to_csv(out/"split.csv", index=False)
    prep = MolecularPreprocessor(a.feature_missing, words(a.residualize), words(a.categorical), a.transform).fit(
        raw[build], p.loc[build], allowed=keep)
    x, observed = prep.transform(raw, p)
    clinical = MetadataDesign(words(a.covariates), words(a.categorical)).fit(p.loc[build])
    c = clinical.transform(p)
    if c.shape[1] == 0:
        raise ValueError("At least one clinical covariate is required for the comparative analysis")
    unknown = unknown_rows(p, [clinical, getattr(prep, "design", None)])
    dump(out/"category_audit.json", {"clinical": clinical.unknown_categories(p.loc[~build]),
                                    "residual": prep.design.unknown_categories(p.loc[~build]) if prep.adjust is not None else {}})
    km = CensoringKM().fit(p.loc[build], a.horizon, a.min_censor_survival)
    # Test outcomes are deliberately not converted into training labels/weights.
    y, w = np.zeros(len(p), int), np.zeros(len(p))
    y[~test], w[~test] = km.labels_weights(p.loc[~test])
    for name, mask in [("build", build), ("tune", tune), ("calibration", cal)]:
        if np.sum(y[mask]) < a.min_events or np.sum((w[mask] > 0) & (y[mask] == 0)) < a.min_events:
            raise ValueError(f"Too few known cases/controls in {name}; shorten horizon or increase data")
    log("DONE", "prepare", f"N={len(p)}; proteins={x.shape[1]}; build={build.sum()}; test={test.sum()}")
    log("START", "reference_quality")
    quality, importance = oof_quality(raw[build], p.loc[build].reset_index(drop=True), features, a, out)
    geometry = MolecularGeometry().fit(x[build], importance[prep.keep], a.dimensions, a.seed)
    z = geometry.transform(x)
    f_names = [features[j] for j in prep.keep]
    profile_weights, profile_names, profile_table = make_profiles(x[build], geometry, f_names, a.seed, a.module_file)
    profile_table.to_csv(out/"molecular_blocks.csv", index=False)
    log("DONE", "reference_quality", f"reliable={quality.reliable_candidate.sum()}/{build.sum()}")
    ids = p[a.id_col].to_numpy(str)
    groups = p[a.group_col].to_numpy(str) if a.group_col else None
    bg = groups[build] if groups is not None else None
    log("START", "reference_panels")
    all_idx = choose_panel(z[build], quality, build.sum(), "all", a.seed)
    full_panel = ReferencePanel().fit(z[build], x[build], ids[build], bg, y[build], w[build],
                                     all_idx, a.donor_neighbors, a.prior_strength)
    tuning, panel_models, panel_specs, candidates = [], {}, {}, []
    def select_model(name, mode, size, ks, risk_mode="local"):
        anchor = choose_panel(z[build], quality, size, mode, a.seed)
        if len(anchor) < 2:
            tuning.append(dict(model=name, panel_size=size, actual_size=len(anchor), k=0,
                               status="insufficient_candidates", AUC_IPCW=np.nan, Brier_IPCW=np.nan))
            return None
        panel = full_panel.subset(anchor)
        best = None
        for k in ks:
            if k >= len(anchor):
                continue
            pred = panel.predict(z[tune], x[tune], observed[tune], ids[tune],
                                 groups[tune] if groups is not None else None, k, risk_mode, False)
            met = metrics(y[tune], w[tune], pred)
            tuning.append(dict(model=name, panel_size=size, actual_size=len(anchor), k=k,
                               status="completed", **met))
            score = np.average(loss(y[tune], pred), weights=w[tune])
            if best is None or score < best[0]:
                best = (score, int(k), met)
        if best is not None:
            panel_models[name] = panel
            panel_specs[name] = dict(k=best[1], mode=risk_mode, requested_size=int(size), actual_size=len(anchor),
                                     selection=mode, validation=best[2])
            return best[0]
        return None
    # Literal proposal: globally best OOF-fit people, a single copied 0/1 outcome.
    select_model("copy1_topfit", "topfit", a.panel_size, [1], "label")
    select_model("copyk_topfit", "topfit", a.panel_size, a.match_k, "label")
    select_model("random_panel", "random", a.panel_size, a.match_k)
    select_model("diversity_panel", "diversity", a.panel_size, a.match_k)
    select_model("all_reference", "all", len(all_idx), a.all_k, "label")
    for size in a.panel_sizes:
        name = f"reliable_{size}"
        value = select_model(name, "reliable", size, a.match_k)
        if value is not None:
            candidates.append((value, name))
    primary_name = f"reliable_{a.panel_size}"
    qualified_panel = primary_name in panel_models
    # The 100-person claim is prespecified. Larger panels are sensitivity models,
    # not permission to quietly replace the primary panel after seeing results.
    selected = primary_name if qualified_panel else "diversity_panel"
    panel_models["panome"] = panel_models[selected]
    panel_specs["panome"] = panel_specs[selected].copy()
    panel_specs["panome"]["selected_candidate"] = selected
    # Full-proteome matching implements the user's literal X-to-X proposal.
    for source, name in [("copy1_topfit", "copy1_fullproteome"), ("copyk_topfit", "copyk_fullproteome")]:
        if source in panel_models:
            panel_models[name] = full_panel.subset(panel_models[source].anchor_indices)
            panel_models[name].z = panel_models[name].x.copy()
            panel_specs[name] = dict(panel_specs[source], space="full_proteome")
    # Same anchor identities, a label-free PCA metric: isolates metric learning.
    du = geometry.unsupervised.n_components_
    unsup_z = z[:, :du]
    unsup = ReferencePanel().fit(unsup_z[build], x[build], ids[build], bg, y[build], w[build],
        panel_models["panome"].anchor_indices, a.donor_neighbors, a.prior_strength)
    panel_models["unsupervised_matching"] = unsup
    panel_specs["unsupervised_matching"] = dict(panel_specs["panome"], space="unsupervised")
    pd.DataFrame(tuning).to_csv(out/"panel_tuning.csv", index=False)
    dump(out/"panel_specs.json", panel_specs)
    # Candidate panels remain saved for sensitivity; only final Panome enters deployment.
    log("DONE", "reference_panels", f"selected={selected}; K={panel_specs['panome']['k']}")
    log("START", "comparators")
    models, model_rows = {}, []
    blocks = {"clinical": c, "elasticnet": np.c_[c, x], "protein_elasticnet": x,
              "clinical_pca": np.c_[c, geometry.unsupervised.transform(x)],
              "clinical_technical": np.c_[c, 1-observed.mean(1)]}
    # Technical control adds plate/batch metadata and missing-assay fraction.
    tech_cols = list(dict.fromkeys(words(a.covariates)+words(a.residualize)))
    technical = MetadataDesign(tech_cols, words(a.categorical)).fit(p.loc[build])
    blocks["clinical_technical"] = np.c_[technical.transform(p), 1-observed.mean(1)]
    for name, block in blocks.items():
        models[name], rows = tune_logistic(block, y, w, build, tune, name, a,
                                          .5 if "elasticnet" in name else 0.)
        model_rows.extend(rows)
        log("DONE", name)
    tree_rows = []
    if a.tree != "none":
        best = None
        for leaves in [7, 15, 31]:
            if a.tree == "lightgbm":
                try:
                    from lightgbm import LGBMClassifier
                except ImportError as exc:
                    raise RuntimeError("--tree lightgbm requires LightGBM; install requirements.txt or explicitly use --tree hist") from exc
                model = LGBMClassifier(n_estimators=a.tree_estimators, num_leaves=leaves, learning_rate=.03,
                                       min_child_samples=50, reg_lambda=10, n_jobs=a.cores,
                                       random_state=a.seed, verbosity=-1, deterministic=True, force_col_wise=True)
            else:
                model = HistGradientBoostingClassifier(max_iter=a.tree_estimators, max_leaf_nodes=leaves,
                    learning_rate=.03, min_samples_leaf=50, l2_regularization=10,
                    early_stopping=False, random_state=a.seed)
            ix = build & (w > 0)
            model.fit(blocks["elasticnet"][ix], y[ix], sample_weight=w[ix]/w[ix].mean())
            pred = model.predict_proba(blocks["elasticnet"][tune])[:, 1]
            score = np.average(loss(y[tune], pred), weights=w[tune])
            tree_rows.append(dict(model=a.tree, leaves=leaves, validation_logloss=score))
            if best is None or score < best[0]:
                best = (score, model)
        models[a.tree], blocks[a.tree] = best[1], blocks["elasticnet"]
    pd.DataFrame(model_rows+tree_rows).to_csv(out/"comparator_tuning.csv", index=False)
    log("DONE", "comparators")
    log("START", "calibration_and_support")
    rawpred = {name: model.predict_proba(blocks[name])[:, 1] for name, model in models.items()}
    infos = {}
    matching_diagnostics = []
    for name, panel in panel_models.items():
        spec = panel_specs[name]
        space = x if spec.get("space") == "full_proteome" else unsup_z if spec.get("space") == "unsupervised" else z
        # No build prediction is needed; self/relative exclusions are still enforced.
        pred, info, jj, ww = panel.predict(space[~build], x[~build], observed[~build], ids[~build],
            groups[~build] if groups is not None else None, spec["k"], spec["mode"])
        fullpred = np.full(len(p), np.nan); fullpred[~build] = pred
        rawpred[name] = fullpred
        diag = info.copy()
        diag.insert(0, a.id_col, ids[~build]); diag["split"] = part[~build]
        diag["model"] = name
        diag["panel_size"] = len(panel.ids)
        diag["label_case_fraction"] = float(panel.labels.mean())
        matching_diagnostics.append(diag)
        if name == "panome":
            infos[name], match_ids, match_weights = info, jj, ww
    log("START", "mosaic_matching")
    all_profiles = x @ profile_weights
    mosaic = MosaicReference().fit(all_profiles[build], x[build], ids[build], bg,
        y[build], w[build], panel_models["panome"].anchor_indices, profile_names, a)
    mosaic_values, mosaic_ids, mosaic_distances = mosaic.block_predictions(
        all_profiles[~build], ids[~build], groups[~build] if groups is not None else None)
    mosaic.tune(mosaic_values[tune[~build]], y[tune], w[tune]).to_csv(out/"mosaic_module_weights.csv", index=False)
    for name, value in [("mosaic_equal", mosaic_values.mean(1)), ("mosaic_weighted", mosaic_values @ mosaic.weights)]:
        rawpred[name] = np.full(len(p), np.nan); rawpred[name][~build] = value
    log("DONE", "mosaic_matching", f"blocks={len(profile_names)}")
    pd.concat(matching_diagnostics, ignore_index=True).to_csv(out/"matching_diagnostics.csv.gz", index=False, compression="gzip")
    calibrators, predictions = {}, {}
    for name, pred in rawpred.items():
        if name in ["copy1_topfit", "copy1_fullproteome"]:
            # Preserve literal COPY Y; calibration would hide its 0/1 behavior.
            predictions[name] = pred.copy()
            continue
        calibrators[name] = RiskCalibrator().fit(pred[cal], y[cal], w[cal])
        predictions[name] = calibrators[name].predict(pred)
    calibrators["panome_clinical"] = RiskCalibrator().fit(rawpred["panome"][cal], y[cal], w[cal],
                                                        rawpred["clinical"][cal])
    predictions["panome_clinical"] = calibrators["panome_clinical"].predict(rawpred["panome"], rawpred["clinical"])
    for name, source, covariate in [("mosaic_clinical", "mosaic_weighted", "clinical"),
                                    ("panome_elasticnet_hybrid", "panome", "elasticnet")]:
        calibrators[name] = RiskCalibrator().fit(rawpred[source][cal], y[cal], w[cal], rawpred[covariate][cal])
        predictions[name] = calibrators[name].predict(rawpred[source], rawpred[covariate])
    dump(out/"calibration_parameters.json", {name: dict(coefficients=obj.coef, logit_mean=obj.mean,
          logit_scale=obj.sd) for name,obj in calibrators.items()})
    relative_tune = tune[~build]
    rule = AcceptanceRule().fit(infos["panome"].loc[relative_tune], a.accept_quantile, a.min_match_ess)
    supported, reasons = rule.apply(infos["panome"], miss[~build] > a.sample_missing, unknown[~build])
    provisional_coverage = float(supported[relative_tune].mean())
    met = panel_specs["panome"]["validation"]
    null = metrics(y[tune], w[tune], np.full(tune.sum(), full_panel.prior))
    readiness_reasons = []
    if not qualified_panel:
        readiness_reasons.append("no_qualifying_reliable_panel; diversity fallback is exploratory")
    if int(quality.reliable_candidate.sum()) < a.panel_size:
        readiness_reasons.append("fewer_than_requested_reliable_people")
    if len(panel_models["panome"].z) < panel_specs["panome"]["requested_size"]:
        readiness_reasons.append("selected_panel_incomplete")
    prediction_warnings = []
    if not geometry.has_supervised_signal:
        prediction_warnings.append("no_stable_linear_metric_signal; PCA metric fallback")
    if not np.isfinite(met["AUC_IPCW"]) or met["AUC_IPCW"] <= .5:
        prediction_warnings.append("no_positive_tuning_discrimination")
    if met["Brier_IPCW"] >= null["Brier_IPCW"]:
        prediction_warnings.append("no_raw_Brier_gain_over_constant_risk")
    if provisional_coverage < a.min_coverage:
        readiness_reasons.append("insufficient_tuning_coverage")
    if quality.loc[quality.reliable_candidate, "horizon_label"].nunique() < 2:
        readiness_reasons.append("reliable_candidates_do_not_include_both_outcomes")
    readiness = dict(status="provisional_pass" if not readiness_reasons else "not_ready",
                     reasons=readiness_reasons, reliable_people=int(quality.reliable_candidate.sum()),
                     prediction_warnings=prediction_warnings,
                     selected=selected, validation_metrics=met, validation_null=null,
                     validation_supported_coverage=provisional_coverage,
                     interpretation="tuning-set diagnostic; selection optimism remains; not a cohort validity test")
    dump(out/"reference_readiness.json", readiness)
    dump(out/"support_thresholds.json", dict(thresholds=rule.thresholds, min_ess=rule.min_ess,
                                            quantile=rule.quantile, source="tune_X_only"))
    panel = panel_models["panome"]
    panel_table = quality.iloc[panel.anchor_indices].copy()
    panel_table["local_risk"], panel_table["donor_ESS"] = panel.local_risk, panel.donor_ess
    panel_table["donor_events"] = panel.donor_events
    panel_table.to_csv(out/"reference_panel.csv", index=False)
    for name, obj in panel_models.items():
        table = quality.iloc[obj.anchor_indices].copy()
        table["local_risk"] = obj.local_risk
        table.to_csv(out/f"panel_{name}.csv", index=False)
    # All learned objects and decision rules are now frozen before test Y is evaluated.
    bundle = dict(version=VERSION, config=vars(a), features=features, retained_features=f_names,
                  prep=prep, clinical=clinical, technical=technical, geometry=geometry, models=models,
                  panels=panel_models, panel_specs=panel_specs, calibrators=calibrators,
                  mosaic=mosaic,
                  rule=rule, readiness=readiness, censoring=km,
                  profile_weights=profile_weights, profile_names=profile_names)
    joblib.dump(bundle, out/"model_bundle.joblib", compress=3)
    dump(out/"MODEL_FROZEN.json", dict(version=VERSION, selection=selected,
                                      test_Y_used_for_fitting_or_model_selection=False,
                                      test_Y_used_for_initial_split_stratification=not bool(a.split_file),
                                      artifacts=fingerprints([out/"model_bundle.joblib"], True)))
    test_in_other = test[~build]
    test_info = infos["panome"].loc[test_in_other].reset_index(drop=True)
    test_info.insert(0, a.id_col, ids[test])
    test_info["supported_match"] = supported[test_in_other]
    test_info["rejection_reason"] = reasons[test_in_other]
    test_info["panel_ready"] = readiness["status"] == "provisional_pass"
    test_info["prediction_released"] = test_info.supported_match & test_info.panel_ready
    test_info["missing_fraction"] = miss[test]
    test_info["unknown_category"] = unknown[test]
    for name, pred in predictions.items():
        test_info[name] = pred[test]
        if name in rawpred:
            test_info[name+"_raw"] = rawpred[name][test]
    test_info["released_net_risk"] = np.where(test_info.prediction_released, predictions["panome_clinical"][test], np.nan)
    test_info.to_csv(out/"test_individuals.csv", index=False)
    test_outcomes = p.loc[test, [a.id_col, "time", "event"]+([a.group_col] if a.group_col else [])].copy()
    test_outcomes.to_csv(out/"test_outcomes.csv", index=False)
    dump(out/"prediction_columns.json", list(predictions))
    jj, ww = match_ids[test_in_other], match_weights[test_in_other]
    pd.DataFrame({a.id_col: np.repeat(ids[test], jj.shape[1]), "rank": np.tile(np.arange(1, jj.shape[1]+1), test.sum()),
                  "reference_eid": panel.ids[jj.ravel()], "copy_weight": ww.ravel(),
                  "reference_local_risk": panel.local_risk[jj.ravel()]}).to_csv(out/"test_reference_matches.csv", index=False)
    profile = x[test] @ profile_weights
    profiles = pd.DataFrame(profile, columns=profile_names); profiles.insert(0, a.id_col, ids[test])
    profiles.to_csv(out/"molecular_profiles.csv", index=False)
    matched_pairs(ids[test], predictions["elasticnet"][test], profile, profile_names,
                  a.pair_caliper, a.max_pairs).to_csv(out/"same_risk_pairs.csv", index=False)
    log("START", "individual_explanations")
    individual_evidence(out, a.id_col, ids[test], x[test], panel, jj, ww, f_names,
        profile, profile_names, mosaic_values[test_in_other], mosaic_ids[test_in_other],
        mosaic_distances[test_in_other], mosaic.weights)
    stability_results = []
    for name in words(a.mask_models):
        if name not in panel_models:
            raise ValueError(f"Unknown/unsupported --mask-models entry: {name}")
        table = masked_copy_stability(raw[test], p.loc[test].reset_index(drop=True), prep, geometry,
            panel_models[name], panel_specs[name], None, ids[test], groups[test] if groups is not None else None, a, out)
        if len(table):
            table.insert(0,"model",name); stability_results.append(table)
    if stability_results:
        pd.concat(stability_results,ignore_index=True).to_csv(out/"masked_copy_stability.csv",index=False)
    log("DONE", "individual_explanations")
    # Save only calibration-free, tune-X thresholds for alternate coverage points.
    thresholds = {str(q): AcceptanceRule().fit(infos["panome"].loc[relative_tune], q, a.min_match_ess)
                  for q in [.5, .75, .9, .95, .99, 1.]}
    joblib.dump(thresholds, out/"coverage_rules.joblib")
    p.loc[test, [a.id_col]+[v for v in ["age", "sex", "center"] if v in p]].to_csv(out/"test_demographics.csv", index=False)
    log("DONE", "calibration_and_support", f"panel={readiness['status']}; tune_coverage={provisional_coverage:.1%}")
    return bundle


def evaluate(out):
    bundle = joblib.load(out/"model_bundle.joblib")
    cfg = bundle["config"]
    info = pd.read_csv(out/"test_individuals.csv", dtype={cfg["id_col"]: str})
    p = pd.read_csv(out/"test_outcomes.csv", dtype={cfg["id_col"]: str})
    if not info[cfg["id_col"]].equals(p[cfg["id_col"]]):
        raise ValueError("Evaluation IDs are not aligned")
    y, w = bundle["censoring"].labels_weights(p)
    names = json.loads((out/"prediction_columns.json").read_text())
    predictions = {name: info[name].to_numpy() for name in names}
    supported = info.supported_match.to_numpy(bool)
    masks = dict(all=np.ones(len(p), bool), supported=supported, rejected=~supported)
    rows = []
    for subset, mask in masks.items():
        for name, pred in predictions.items():
            rows.append(dict(model=name, subset=subset, n=int(mask.sum()), cases=int(y[mask].sum()),
                             known=int((w[mask]>0).sum()), coverage=float(mask.mean()),
                             **metrics(y[mask], w[mask], pred[mask])))
    pd.DataFrame(rows).to_csv(out/"test_metrics.csv", index=False)
    diag = pd.read_csv(out/"matching_diagnostics.csv.gz")
    diag = diag[diag.split == "test"]
    explanation = diag.groupby("model").agg(panel_size=("panel_size", "first"),
        label_case_fraction=("label_case_fraction", "first"),
        mean_reference_ESS=("reference_ESS", "mean"), mean_donor_ESS=("donor_ESS", "mean"),
        mean_matching_distance=("nearest_distance", "mean"),
        mean_profile_reconstruction_RMSE=("reconstruction_RMSE", "mean"),
        used_reference_people=("nearest_reference", "nunique")).reset_index()
    comparison = pd.DataFrame(rows).query('subset == "all"').merge(explanation, how="left", on="model")
    mosaic_table = pd.read_csv(out/"test_mosaic_matches.csv.gz")
    mosaic_mask = comparison.model.str.startswith("mosaic")
    comparison.loc[mosaic_mask,"n_modules"] = len(bundle["profile_names"])
    comparison.loc[mosaic_mask,"mean_distinct_module_references"] = mosaic_table.groupby(cfg["id_col"]).reference_eid.nunique().mean()
    if (out/"masked_copy_stability.csv").exists():
        stability = pd.read_csv(out/"masked_copy_stability.csv")
        if "model" in stability:
            stability = stability.groupby("model")[["masked_RMSE","match_Jaccard","absolute_raw_risk_change"]].mean().reset_index()
            comparison = comparison.merge(stability,how="left",on="model")
    comparison.to_csv(out/"approach_comparison.csv", index=False)
    groups = p[cfg["group_col"]].to_numpy(str) if cfg["group_col"] else None
    bootstrap_contrasts(y, w, predictions, {k: masks[k] for k in ["all", "supported"]},
                        groups, cfg["bootstrap"], cfg["seed"]+200).to_csv(out/"paired_contrasts.csv", index=False)
    curve = []
    for q, rule in joblib.load(out/"coverage_rules.joblib").items():
        mask, _ = rule.apply(info, info.missing_fraction.to_numpy()>cfg["sample_missing"], info.unknown_category.to_numpy(bool))
        for name, pred in predictions.items():
            curve.append(dict(quantile=float(q), coverage=float(mask.mean()), n=int(mask.sum()), model=name,
                              **metrics(y[mask], w[mask], pred[mask])))
    pd.DataFrame(curve).to_csv(out/"coverage_curve.csv", index=False)
    calibration = []
    for name, pred in predictions.items():
        bins = pd.qcut(pd.Series(pred), 5, labels=False, duplicates="drop").to_numpy()
        for b in np.unique(bins[np.isfinite(bins)]):
            ix = bins == b
            calibration.append(dict(model=name, bin=int(b)+1, n=int(ix.sum()),
                                    predicted=float(pred[ix].mean()), observed_IPCW=float(np.average(y[ix], weights=w[ix]))))
    pd.DataFrame(calibration).to_csv(out/"test_calibration.csv", index=False)
    demographic = pd.read_csv(out/"test_demographics.csv")
    subgroup = []
    for col in [c for c in ["sex", "center"] if c in demographic]:
        for level, sub in demographic.groupby(col, dropna=False):
            ix = sub.index.to_numpy()
            subgroup.append(dict(variable=col, level=level, n=len(ix), supported_fraction=float(supported[ix].mean()),
                                 **metrics(y[ix], w[ix], predictions["panome_clinical"][ix])))
    pd.DataFrame(subgroup).to_csv(out/"subgroup_coverage.csv", index=False)
    figures(out)
    audit = json.loads((out/"cohort_audit.json").read_text())
    write_report(out, audit, bundle["readiness"], bundle["panel_specs"]["panome"], cfg["horizon"])
    dump(out/"DONE.json", dict(version=VERSION, n_test=len(p), cases_by_horizon=int(y.sum()),
                               supported_fraction=float(supported.mean()),
                               panel_status=bundle["readiness"]["status"]))
    log("DONE", "test_evaluation", f"N={len(p)}; cases_by_horizon={y.sum()}; supported={supported.mean():.1%}")


def project(run_dir, phe_path, omics_path, output, r_bin="Rscript", met_input="named"):
    bundle = joblib.load(Path(run_dir)/"model_bundle.joblib")
    cfg = bundle["config"]
    idcol = cfg["id_col"]
    cols = list(dict.fromkeys([idcol]+words(cfg["covariates"])+words(cfg["residualize"])+([cfg["group_col"]] if cfg["group_col"] else [])))
    p = read_table(phe_path, idcol, cols, r_bin)
    omics = read_table(omics_path, idcol, r_bin=r_bin)
    if cfg["biom"] == "prot":
        omics.columns = [c if c == idcol else str(c).upper() for c in omics.columns]
    elif met_input == "raw":
        omics, _ = map_metabolites(omics, cfg["met_map"], idcol)
    if omics.columns.duplicated().any():
        raise ValueError("Assay names collide after normalization")
    if not set(omics[idcol]) <= set(p[idcol]):
        raise ValueError("Some query IDs lack baseline covariates")
    p = p.set_index(idcol).loc[omics[idcol]].reset_index()
    # Excluded training assays may be omitted, retained assays must exist.
    required = set(bundle["retained_features"])
    if not required <= set(omics):
        raise ValueError("Missing trained assays: "+", ".join(sorted(required-set(omics))[:20]))
    for name in set(bundle["features"])-set(omics):
        omics[name] = np.nan
    raw = numeric(omics, bundle["features"], "projection")
    prep = bundle["prep"]
    x, observed = prep.transform(raw, p)
    z = bundle["geometry"].transform(x)
    c = bundle["clinical"].transform(p)
    ids = p[idcol].to_numpy(str)
    groups = p[cfg["group_col"]].to_numpy(str) if cfg["group_col"] else None
    if cfg["group_col"] and p[cfg["group_col"]].isna().any():
        raise ValueError("Missing query family groups")
    panel, spec = bundle["panels"]["panome"], bundle["panel_specs"]["panome"]
    rawrisk, info, jj, ww = panel.predict(z, x, observed, ids, groups, spec["k"], spec["mode"])
    cp = bundle["models"]["clinical"].predict_proba(c)[:, 1]
    risk = bundle["calibrators"]["panome_clinical"].predict(rawrisk, cp)
    blocks = {"clinical": c, "elasticnet": np.c_[c,x], "protein_elasticnet": x,
              "clinical_pca": np.c_[c,bundle["geometry"].unsupervised.transform(x)],
              "clinical_technical": np.c_[bundle["technical"].transform(p),1-observed.mean(1)]}
    raw_predictions = {}
    for name, model in bundle["models"].items():
        block = blocks[name] if name in blocks else blocks["elasticnet"]
        raw_predictions[name] = model.predict_proba(block)[:,1]
    for name, ref in bundle["panels"].items():
        params = bundle["panel_specs"][name]
        space = x if params.get("space") == "full_proteome" else z[:, :bundle["geometry"].unsupervised.n_components_] if params.get("space") == "unsupervised" else z
        raw_predictions[name] = ref.predict(space, x, observed, ids, groups, params["k"], params["mode"], False)
    module_values, module_ids, module_distances = bundle["mosaic"].block_predictions(x @ bundle["profile_weights"], ids, groups)
    raw_predictions["mosaic_equal"] = module_values.mean(1)
    raw_predictions["mosaic_weighted"] = module_values @ bundle["mosaic"].weights
    for name, values in raw_predictions.items():
        info[name] = bundle["calibrators"][name].predict(values) if name in bundle["calibrators"] else values
    for name, source, covariate in [("panome_clinical", "panome", "clinical"),
                                    ("mosaic_clinical", "mosaic_weighted", "clinical"),
                                    ("panome_elasticnet_hybrid", "panome", "elasticnet")]:
        info[name] = bundle["calibrators"][name].predict(raw_predictions[source],raw_predictions[covariate])
    unknown = unknown_rows(p, [bundle["clinical"], getattr(prep, "design", None)])
    missing = 1-observed.mean(1)
    supported, reason = bundle["rule"].apply(info, missing>cfg["sample_missing"], unknown)
    release = supported & (bundle["readiness"]["status"] == "provisional_pass")
    info.insert(0, idcol, ids)
    info["exploratory_net_risk"] = risk
    info["prediction_released"], info["supported_match"] = release, supported
    info["released_net_risk"] = np.where(release, risk, np.nan)
    info["rejection_reason"] = reason
    info["panel_status"] = bundle["readiness"]["status"]
    info["missing_fraction"], info["unknown_category"] = missing, unknown
    info.to_csv(output, index=False)
    module_names = bundle["profile_names"]
    pd.DataFrame({idcol: np.repeat(ids,len(module_names)),"module":np.tile(module_names,len(ids)),
                  "reference_eid":module_ids.ravel(),"module_distance":module_distances.ravel(),
                  "module_local_risk":module_values.ravel(),
                  "mixture_weight":np.tile(bundle["mosaic"].weights,len(ids))}).to_csv(
                      Path(output).with_name(Path(output).stem+"_mosaic.csv.gz"),index=False,compression="gzip")
    log("DONE", "project", f"N={len(p)}; released={release.sum()}")
