"""State profiles, individual attribution and descriptive matched-risk pairs."""
import json
import joblib
import numpy as np
import pandas as pd
import torch
from representation import load_encoder
from common import dump, log, bh

def state_association(p, clinical, graph, a, out):
    if graph["soft"].shape[1]<2:
        dump(out/"state_association_status.json",dict(status="unavailable",reason="No multistate solution"))
        return
    te = p.split.to_numpy()=="test"
    columns = {f"cov_{j}":clinical[te,j] for j in range(clinical.shape[1])
               if np.std(clinical[te,j])>1e-8}
    d = pd.DataFrame(columns)
    # Fixed S1 reference. If S1 is absent in test, do not silently redefine it.
    if not np.any(graph["state"][te]==0):
        dump(out/"state_association_status.json",dict(status="unavailable",reason="Reference state S1 absent from test"))
        return
    for state in range(1,graph["soft"].shape[1]):
        indicator = (graph["state"][te]==state).astype(float)
        if 0<indicator.sum()<len(indicator):
            d[f"state_S{state+1}"] = indicator
    try:
        if a.outcome_type=="survival":
            from lifelines import CoxPHFitter
            from lifelines.statistics import proportional_hazard_test
            d["time"],d["event"] = p.loc[te,"time"].to_numpy(),p.loc[te,"event"].to_numpy()
            fitted = CoxPHFitter().fit(d,"time","event")
            result = fitted.summary.reset_index().rename(columns={"covariate":"term"})
            ph = proportional_hazard_test(fitted,d,time_transform="rank").summary.reset_index()
            ph.to_csv(out/"test_state_PH.csv",index=False)
        else:
            import statsmodels.api as sm
            fitted = sm.OLS(p.loc[te,"target"].to_numpy(),sm.add_constant(d)).fit(cov_type="HC3")
            ci = fitted.conf_int()
            result = pd.DataFrame(dict(term=fitted.params.index,beta=fitted.params.values,
                se=fitted.bse.values,p=fitted.pvalues.values,lower=ci.iloc[:,0].values,upper=ci.iloc[:,1].values))
        if "term" not in result:
            result = result.rename(columns={result.columns[0]:"term"})
        result = result[result.term.astype(str).str.startswith("state_")].copy()
        result["fdr_states"] = bh(result.p)
        result.to_csv(out/"test_state_association.csv",index=False)
        dump(out/"state_association_status.json",dict(status="completed",interpretation="Exploratory, discovery-defined state associations"))
    except (ValueError,ArithmeticError,np.linalg.LinAlgError) as exc:
        dump(out/"state_association_status.json",dict(status="failed",reason=str(exc)))
    except Exception as exc:
        # Lifelines convergence exceptions are not all ValueError subclasses.
        if "Convergence" not in type(exc).__name__:
            raise
        dump(out/"state_association_status.json",dict(status="failed",reason=str(exc)))

def attribute_ae(x, observed, p, features, a, root, out):
    path = root/"s5_predict/clinical_ae.joblib"
    if not path.exists() or a.attribution_samples==0:
        return
    model = load_encoder(root/"s3_representation/ae.pt")
    fitted = joblib.load(path)
    coefficients = np.zeros(len(fitted.active))
    coefficients[fitted.active] = np.asarray(fitted.model.coef_).ravel()
    beta = torch.as_tensor((coefficients/fitted.scaler.scale_)[-a.latent:],dtype=torch.float32)
    test = np.flatnonzero(p.split.to_numpy()=="test")
    rng = np.random.default_rng(a.seed)
    selected = np.sort(rng.choice(test,min(a.attribution_samples,len(test)),replace=False))
    rows,diagnostics = [],[]
    for start in range(0,len(selected),32):
        ix = selected[start:start+32]
        target = torch.as_tensor(np.array(x[ix]))
        mask = torch.as_tensor(np.array(observed[ix]),dtype=torch.float32)
        derivative = torch.zeros_like(target)
        for step in range(a.ig_steps):
            point = (target*((step+.5)/a.ig_steps)).requires_grad_(True)
            value = (model.encode(point,mask)*beta).sum()
            derivative += torch.autograd.grad(value,point)[0].detach()/a.ig_steps
        contributions = (target*derivative).numpy()
        with torch.no_grad():
            delta = ((model.encode(target,mask)-model.encode(torch.zeros_like(target),mask))*beta).sum(1).numpy()
        for local,index in enumerate(ix):
            for j in np.argsort(-np.abs(contributions[local]))[:20]:
                rows.append(dict(person_id=p.iloc[index][a.id_col],feature=features[j],
                    contribution=float(contributions[local,j]),observed=bool(observed[index,j]),
                    model="clinical_ae",scale="log_hazard" if a.outcome_type=="survival" else "target_units"))
            diagnostics.append(dict(person_id=p.iloc[index][a.id_col],prediction_delta=float(delta[local]),
                attribution_sum=float(contributions[local].sum()),
                integration_error=float(contributions[local].sum()-delta[local])))
    pd.DataFrame(rows).to_csv(out/"person_ae_attributions.csv",index=False)
    pd.DataFrame(diagnostics).to_csv(out/"attribution_completeness.csv",index=False)

def matched_pairs(x,p,z,graph,features,a,root,out):
    test = pd.read_csv(root/"s5_predict/test_predictions.csv",dtype={a.id_col:str})
    chosen = a.match_model
    if a.outcome_type=="survival":
        key = f"{chosen}_net_risk_{a.primary_horizon:g}y"
        caliper = a.pair_caliper
    else:
        key = f"{chosen}_prediction"
        train_sd = float(p.loc[p.split=="train","target"].std())
        caliper = a.pair_caliper*train_sd
    if key not in test:
        dump(out/"matched_pairs_status.json",dict(status="unavailable",reason=f"Prespecified reference missing: {key}"))
        return
    # Pair only by locked predictions and molecular state, never test event labels.
    order = np.argsort(test[key].to_numpy(),kind="stable")
    used,rows,feature_rows = set(),[],[]
    lookup = {str(v):i for i,v in enumerate(p[a.id_col])}
    for pos,i in enumerate(order):
        if i in used:
            continue
        state = test.iloc[i].state
        candidate = None
        for j in order[pos+1:]:
            if test.iloc[j][key]-test.iloc[i][key]>caliper:
                break
            if j not in used and test.iloc[j].state != state:
                candidate = j
                break
        if candidate is None:
            continue
        j = candidate
        used.update([i,j])
        one,two = test.iloc[i],test.iloc[j]
        ai,bi = lookup[str(one[a.id_col])],lookup[str(two[a.id_col])]
        number = len(rows)+1
        rows.append(dict(pair=number,person_1=one[a.id_col],person_2=two[a.id_col],
            state_1=int(one.state),state_2=int(two.state),prediction_1=float(one[key]),
            prediction_2=float(two[key]),difference=float(two[key]-one[key]),
            reference_model=chosen,scale="net_risk" if a.outcome_type=="survival" else "target_units",
            euclidean_AE_distance=float(np.linalg.norm(z[ai]-z[bi]))))
        delta = x[ai]-x[bi]
        for k in np.argsort(-np.abs(delta))[:20]:
            feature_rows.append(dict(pair=number,feature=features[k],
                person_1_z=float(x[ai,k]),person_2_z=float(x[bi,k]),difference=float(delta[k])))
        if len(rows)>=a.max_pairs:
            break
    pd.DataFrame(rows,columns=["pair","person_1","person_2","state_1","state_2","prediction_1",
        "prediction_2","difference","reference_model","scale","euclidean_AE_distance"]).to_csv(
        out/"matched_risk_pairs.csv",index=False)
    pd.DataFrame(feature_rows,columns=["pair","feature","person_1_z","person_2_z","difference"]).to_csv(
        out/"matched_pair_features.csv",index=False)
    dump(out/"matched_pairs_status.json",dict(status="completed",n_pairs=len(rows),caliper=caliper,
        interpretation="Descriptive examples selected on predicted risk and state; not evidence of distinct causal mechanisms."))

def interpret(x, observed,p,z,graph,atlas,features,clinical,a,root,out):
    tr = p.split.to_numpy()=="train"
    profile = []
    for state in range(atlas.n_states):
        rows = np.flatnonzero(tr)[atlas.labels==state]
        masked = np.where(observed[rows],x[rows],np.nan)
        means = np.nanmean(masked,axis=0)
        for j,feature in enumerate(features):
            profile.append(dict(state=state+1,feature=feature,mean_standardized=float(means[j]),
                                 n_observed=int(observed[rows,j].sum())))
    pd.DataFrame(profile).to_csv(out/"training_state_profiles.csv",index=False)
    persons = p[[a.id_col,"split"]].copy()
    persons["state"],persons["novel"] = graph["state"]+1,graph["novelty"]
    persons["neighbor_distance"] = graph["neighbor_distance"]
    persons["effective_neighbors"] = graph["effective_neighbors"]
    persons["observed_fraction"] = observed.mean(1)
    for j in range(z.shape[1]):
        persons[f"AE{j+1}"] = z[:,j]
    for j in range(graph["soft"].shape[1]):
        persons[f"state_weight_{j+1}"] = graph["soft"][:,j]
    for j in range(graph["diffusion"].shape[1]):
        persons[f"diffusion{j+1}"] = graph["diffusion"][:,j]
    persons.to_csv(out/"person_molecular_states.csv",index=False)
    counts = persons.groupby(["split","state"]).size().reset_index(name="n")
    counts.to_csv(out/"state_counts.csv",index=False)
    # Independent metadata checks are descriptive; no circular protein cluster P values.
    state_association(p,clinical,graph,a,out)
    attribute_ae(x,observed,p,features,a,root,out)
    matched_pairs(x,p,z,graph,features,a,root,out)
    attention_path = root/"s5_predict/person_attention.npy"
    if attention_path.exists():
        attention = np.load(attention_path,mmap_mode="r")
        rows = []
        for i in np.flatnonzero(p.split.to_numpy()=="test"):
            for rank,j in enumerate(np.argsort(-attention[i])[:5],1):
                ref = graph["neighbors"][i,j]
                rows.append(dict(person_id=p.iloc[i][a.id_col],rank=rank,
                    reference_person_id=atlas.train_ids[ref],attention_weight=float(attention[i,j])))
        pd.DataFrame(rows).to_csv(out/"person_reference_attention.csv",index=False)
    audit = json.loads((root/"s2_preprocess/cohort_audit.json").read_text())
    status = pd.read_csv(root/"s5_predict/model_status.csv")
    successful = status.loc[status.status=="completed","model"].tolist()
    text = f"""# Panome {a.trait}: {a.biom}

This run uses {'synthetic software-test data' if a.demo else 'the configured real input files'}.
The outcome is {a.outcome_type}. After eligibility and molecular QC there are {len(p)} people and {len(features)} features.
Split counts are {p.split.value_counts().to_dict()}. The clinical covariates are {a.covariates}; age uses a training-fitted spline.
Molecular residualization uses {a.residualize or 'none'}. This covariate set is not a validated clinical risk calculator.

The train-only graph contains {atlas.n_states} states and {atlas.summary['n_components']} connected components.
Partition status is {atlas.summary['state_status']}. Mean edge-dropout ARI is {atlas.summary['perturbation_ARI_mean']}.
Low stability is a limitation; one-state fallback does not establish a discrete subtype.
State weights summarize reference-neighborhood membership, not biological pathway percentages.

Completed prediction models: {', '.join(successful)}.
Read s5_predict/model_status.csv for unavailable or failed models, test_metrics.csv for held-out performance,
paired_contrasts.csv for paired increments, and metric_limitations.json for unsupported horizons.
The Transformer comparison shares one frozen, label-free module encoder. Only the supervised prediction head changes.
The random-reference control has the same row-attention architecture. It tests whether molecular neighbors matter.
Neural survival heads use a piecewise-exponential proportional-hazards likelihood with observed exposure time.
Classical survival models retain right censoring. Death is censored; probabilities are net risks.

Individual outputs retain molecular coordinates, state weights, distance and missingness.
AE integrated gradients describe only the clinical_ae model on the processed input scale, holding metadata and missingness fixed.
They do not explain the row-attention model. Reference attention weights describe model weighting, not causal influence.
Matched-risk pairs are descriptive and selected without test outcomes; molecular separation alone is not independent biological evidence.

Landmark analyses use frozen baseline scores among people still event-free and observed at the landmark.
They assess residual-follow-up discrimination, not causal direction. No protein is automatically excluded or downweighted.
Bootstrap intervals are conditional on the fitted models; family groups are resampled together when supplied.
Full retraining across seeds/splits, external validation, independent phenotypes and competing-death risks remain separate research requirements.
The baseline atlas does not infer disease trajectories, interventions or molecular ancestry.
"""
    (out/"REPORT.md").write_text(text,encoding="utf-8")
