"""Aggregate manuscript contrasts from frozen person-level predictions."""
import numpy as np
import pandas as pd
from scipy.linalg import qr
from scipy.optimize import minimize
from scipy.special import expit

from core2 import atomic_table
from evaluate2 import aligned, boolean_mask, interval

PRIMARY = "csx.population"
SELECT = "selected_support_0p20"
CORE = ("abm", PRIMARY, "disco")
SETS = ("all", SELECT)


def populations(ancestry):
    return [("ALL", np.ones(len(ancestry), bool))] + [
        (name, ancestry == name) for name in sorted(set(ancestry))]


def comparisons(ancestry, methods=()):
    pairs = [("abm", PRIMARY), ("abm", "disco"), ("disco", PRIMARY),
             ("disco", "csx.input"), (PRIMARY, "cs.best"),
             (PRIMARY, "cs.EUR"), (PRIMARY, "cs.linear_combination"),
             ("abm", "abm.global_residual")]
    if ancestry in ("AFR", "EAS", "SAS"):
        pairs.append((PRIMARY, "cs." + ancestry))
    return pairs


def summarize_draws(row, point, replicates):
    low, high, n, status = interval(replicates, len(replicates))
    return dict(row, difference=point, ci_low=low, ci_high=high, valid_draws=n,
                ci_status=status, bootstrap_replicates=len(replicates))


def family_contrasts(y, predictions, groups, ancestry, selected, trait, reps, seed, only_pairs=None):
    """Resample families once across methods and nested sets. Recompute R² denominators."""
    binary = trait == "t2dm"
    rows = []
    for i, (pop, pop_mask) in enumerate(populations(ancestry)):
        take = pop_mask & np.isfinite(y)
        yy = y[take]
        pairs = comparisons(pop, predictions)
        if "direct.blocknet" in predictions:
            from ldl_report import ABLATIONS
            pairs.extend([(m,r) for m,r in ABLATIONS if m in predictions and r in predictions])
            pairs.extend([("direct.blocknet",PRIMARY),("abm","direct.blocknet")])
        for model in ("block.glmnet","block.lightgbm","abm.deep"):
            if model in predictions:
                pairs.extend([(model,PRIMARY),("abm",model)])
        if only_pairs is not None:
            pairs=only_pairs
        methods = sorted({m for pair in pairs for m in pair} | {"covariate"})
        available = [m for m in methods if m in predictions]
        if yy.size:
            centered = yy - yy.mean()
            values = np.column_stack([np.ones(len(yy)), centered, centered**2] +
                                     [(predictions[m][take] - yy)**2 for m in available])
        else:
            values = np.zeros((0, 3 + len(available)))
        _, group = np.unique(groups[take], return_inverse=True)
        ng = len(np.unique(group))
        masks = [np.ones(len(yy), bool), selected[take]]
        totals = np.zeros((ng, 2, values.shape[1]))
        for s, mask in enumerate(masks):
            for col in range(values.shape[1]):
                totals[:, s, col] = np.bincount(group[mask], weights=values[mask, col], minlength=ng)
        # Each replicate is small after reducing family summaries to two sets.
        raw = np.full((reps + 1, 2, values.shape[1]), np.nan)
        raw[0] = totals.sum(axis=0)
        if ng >= 2:
            rng = np.random.default_rng(seed + i)
            for b in range(1, reps + 1):
                weights = np.bincount(rng.integers(ng, size=ng), minlength=ng)
                raw[b] = np.einsum("g,gsc->sc", weights, totals)
        n = raw[:, :, 0]
        with np.errstate(invalid="ignore", divide="ignore"):
            sst = raw[:, :, 2] - raw[:, :, 1]**2 / n
        lookup = {m: raw[:, :, j + 3] for j, m in enumerate(available)}
        no = np.full(n.shape, np.nan)

        def divide(a, denom):
            return np.divide(a, denom, out=no.copy(), where=denom > 0)

        for a, ref in pairs:
            diff = lookup.get(a, no) - lookup.get(ref, no)
            metrics = {"Brier" if binary else "MSE": divide(diff, n)}
            metrics["relative_loss_reduction_pct"] = 100 * divide(-diff, lookup.get(ref, no))
            if not binary:
                metrics["predictive_R2"] = divide(-diff, sst)
                metrics["partial_R2"] = divide(-diff, lookup.get("covariate", no))
            for s, selection in enumerate(SETS):
                nfamilies = len(np.unique(group[masks[s]]))
                for metric, vector in metrics.items():
                    point = vector[0, s]
                    draws = vector[1:, s] if nfamilies >= 2 else np.full(reps, np.nan)
                    row = dict(trait=trait, ancestry=pop, selection=selection, comparison="method_minus_reference",
                               method=a, reference=ref, metric=metric, n=int(masks[s].sum()), families=nfamilies,
                               n_cases=int(yy[masks[s]].sum()) if binary else np.nan,
                               availability="available" if a in lookup and ref in lookup else "method_unavailable")
                    rows.append(summarize_draws(row, point, draws))
        if not binary and only_pairs is None:
            for m in CORE:
                for metric, denom in [("predictive_R2", sst), ("partial_R2", lookup.get("covariate", no))]:
                    val = 1 - divide(lookup.get(m, no), denom)
                    delta = val[:, 1] - val[:, 0]
                    nf = len(np.unique(group[masks[1]]))
                    row = dict(trait=trait, ancestry=pop, selection=SELECT, comparison="selected_minus_all",
                               method=m, reference=m, metric=metric, n=int(masks[1].sum()), families=nf,
                               n_cases=np.nan, availability="available" if m in lookup else "method_unavailable")
                    rows.append(summarize_draws(row, delta[0], delta[1:] if nf >= 2 else np.full(reps, np.nan)))
    return pd.DataFrame(rows)


def twin_summary(test, abm, y, selected, trait, predictions):
    support = pd.to_numeric(abm.n_support, errors="raise").to_numpy(float)
    ess = pd.to_numeric(abm.support_group_ess, errors="raise").to_numpy(float)
    supported = boolean_mask(abm.supported)
    rows = []
    for pop, mask in populations(test.ancestry.to_numpy()):
        for label, membership in zip((*SETS, "rejected_support_0p20"),
                                      (np.ones(len(test), bool), selected, ~selected)):
            ix = mask & membership
            good = ix & np.isfinite(y)
            row = dict(trait=trait, ancestry=pop, selection=label, roster_n=int(mask.sum()),
                       n=int(ix.sum()), evaluable_n=int(good.sum()), coverage=ix.sum()/mask.sum(),
                       observed_outcome_mean=np.mean(y[good]) if good.any() else np.nan,
                       n_cases=int(y[good].sum()) if trait=="t2dm" else np.nan,
                       n_zero_neighbors=int((ix & (support == 0)).sum()),
                       supported_fraction=supported[ix].mean() if ix.any() else np.nan)
            fields = {"n_support": support, "family_ess": ess}
            for name in ("n_support_groups","support_fraction","n_donors_used"):
                if name in abm:
                    fields[name] = pd.to_numeric(abm[name],errors="raise").to_numpy(float)
            for name, field in fields.items():
                for p in (.1, .5, .9):
                    row[f"{name}_p{int(p*100)}"] = np.quantile(field[ix], p) if ix.any() else np.nan
            change = predictions["abm"][ix] - predictions[PRIMARY][ix]
            row["fraction_abm_changed_csx"] = np.mean(abs(change)>1e-12) if ix.any() else np.nan
            row["mean_absolute_prediction_change"] = np.mean(abs(change)) if ix.any() else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def score_association(y, q, z, groups, binary):
    """Descriptive TEST association conditional on a full-rank covariate design."""
    base = dict(status="not_estimable", partial_r2=np.nan, incremental_r2=np.nan,
                log_or_per_sd=np.nan, log_or_low=np.nan, log_or_high=np.nan)
    if len(y) <= q.shape[1] + 2 or np.std(z) <= 1e-12 or np.std(y) == 0:
        return base
    z = (z - z.mean()) / z.std()
    residual = z - q @ (q.T @ z)
    if residual @ residual < 1e-10:
        return base
    if not binary:
        ry = y - q @ (q.T @ y)
        sst = np.sum((y - y.mean())**2)
        if ry @ ry <= 1e-12 * sst:
            return base
        improvement = (ry @ residual)**2 / (residual @ residual)
        return dict(base, status="available", partial_r2=improvement/(ry @ ry), incremental_r2=improvement/sst)
    # Orthonormal columns span the same nuisance design. Rescale for optimization.
    design = np.column_stack([q * np.sqrt(len(y)), z])
    def objective(beta):
        eta = design @ beta
        return np.mean(np.logaddexp(0, eta)-y*eta), design.T @ (expit(eta)-y)/len(y)
    fit = minimize(objective, np.zeros(design.shape[1]), jac=True, method="BFGS",
                   options={"gtol":1e-7, "maxiter":200})
    if not fit.success:
        return base
    prob = expit(design @ fit.x)
    info = design.T @ ((prob * (1-prob))[:, None] * design)
    if not np.isfinite(info).all() or np.linalg.cond(info) > 1e12:
        return base
    vector = np.linalg.solve(info, np.eye(len(info))[:, -1])
    influence = np.bincount(groups, weights=(y-prob)*(design @ vector))
    count = len(influence)
    if count < 2:
        return base
    se = np.sqrt((influence @ influence) * count/(count-1) * (len(y)-1)/(len(y)-len(info)))
    beta = fit.x[-1]
    return dict(base, status="low_case_control_count_lt20" if min(y.sum(), len(y)-y.sum()) < 20 else "available",
                log_or_per_sd=beta, log_or_low=beta-1.96*se, log_or_high=beta+1.96*se)


def genetic_terms(path, test, benchmark, y, covariates, trait):
    """Evaluate exported anchor genetic components and interpolation on the same complete cases."""
    rows = []
    columns = ["trait","ancestry","score","n","n_cases","n_complete_covariates",
               "matched_anchor","status","partial_r2","incremental_r2","log_or_per_sd","log_or_low","log_or_high"]
    if not path.is_file() or "disco.genetic_score" not in benchmark:
        return pd.DataFrame([dict(trait=trait, status="interpolation_unavailable")], columns=columns)
    raw = pd.read_csv(path, sep="\t", dtype={"id":str}, usecols=lambda c: c in {"id","disco"} or c.startswith("component_"))
    parts = aligned(raw, test.id, "Disco interpolation components")
    if not np.allclose(parts.disco, pd.to_numeric(benchmark["disco.genetic_score"]), rtol=1e-7, atol=1e-9):
        raise ValueError(f"{trait}: full-model interpolation differs from the frozen TEST genetic score")
    scores = parts.drop(columns="id").apply(pd.to_numeric)
    cov = test[covariates].apply(pd.to_numeric).to_numpy(float)
    for pop, mask in populations(test.ancestry.to_numpy()):
        use = mask & np.isfinite(y) & np.isfinite(cov).all(1) & np.isfinite(scores).all(axis=1).to_numpy()
        yy = y[use]
        X = cov[use].copy()
        if len(yy):
            scale = X.std(0)
            X = (X - X.mean(0)) / np.where(scale > 0, scale, 1)
            X = np.column_stack([np.ones(len(yy)), X])
            q, r, _ = qr(X, mode="economic", pivoting=True)
            q = q[:, :np.linalg.matrix_rank(r)]
        else:
            q = np.empty((0, 0))
        _, groups = np.unique(test.group_id[use], return_inverse=True)
        for name in scores:
            stat = score_association(yy, q, scores[name].to_numpy()[use], groups, trait == "t2dm")
            rows.append(dict(trait=trait, ancestry=pop, score=name, n=len(yy),
                             n_cases=int(yy.sum()) if trait=="t2dm" else np.nan,
                             n_complete_covariates=int((mask & np.isfinite(cov).all(1)).sum()),
                             matched_anchor=("component_"+pop in scores) if pop!="ALL" else None, **stat))
    return pd.DataFrame(rows, columns=columns)


def manuscript_evaluation(trait, out, cohort, outcomes, benchmark, abm, components, covariates, reps, seed):
    test = cohort.loc[cohort.role.eq("test")].reset_index(drop=True)
    truth = aligned(outcomes, test.id, "test outcomes")
    base = aligned(benchmark, test.id, "benchmark predictions")
    local = aligned(abm, test.id, "ABM support")
    selected = boolean_mask(local[SELECT])
    y = pd.to_numeric(truth.y).to_numpy(float)
    predictions = {c[:-11]: pd.to_numeric(frame[c]).to_numpy(float) for frame in (base, local)
                   for c in frame if c.endswith(".prediction")}
    paired = family_contrasts(y, predictions, test.group_id.to_numpy(), test.ancestry.to_numpy(),
                              selected, trait, reps, seed)
    atomic_table(paired, out / f"{trait}.contrasts.tsv")
    atomic_table(twin_summary(test, local, y, selected, trait, predictions), out / f"{trait}.support.tsv")
    atomic_table(genetic_terms(components, test, base, y, covariates, trait), out / f"{trait}.anchor_associations.tsv")
