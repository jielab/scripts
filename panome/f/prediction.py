"""Locked classical and individual-conditioned prediction comparisons."""
import warnings
import joblib
import numpy as np
import pandas as pd
from scipy.stats import norm, t as student_t
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge, ElasticNet
from sklearn.ensemble import HistGradientBoostingRegressor
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sksurv.ensemble import GradientBoostingSurvivalAnalysis
from common import words, bh, dump, log
from preprocess import MetadataDesign
from evaluation import survival_y, primary_metric, evaluate, paired_uncertainty, landmarks
from neural import NeuralModel, random_neighbors

def marginal_scan(x, p, features, kind):
    """Marginal prediction screen, not a fully adjusted etiological PWAS."""
    if kind == "quantitative":
        y = p.target.to_numpy(float, copy=True)
        y -= y.mean()
        z = x.astype(float) - x.mean(0)
        xx = (z*z).sum(0)
        beta = (z.T@y)/np.maximum(xx, 1e-12)
        rss = np.maximum((y*y).sum()-beta*beta*xx, 0)
        se = np.sqrt(rss/max(len(y)-2, 1)/np.maximum(xx, 1e-12))
        pv = 2*student_t.sf(np.abs(beta/np.maximum(se, 1e-12)), max(len(y)-2,1))
        result = pd.DataFrame(dict(feature=features, beta=beta, se=se, p=pv,
                                   converged=np.isfinite(pv)))
    else:
        y = survival_y(p)
        order = np.argsort(y["time"])
        times, event = y["time"][order], y["event"][order]
        starts = np.r_[0, np.flatnonzero(np.diff(times))+1]
        deaths = np.add.reduceat(event.astype(float), starts)
        rows = []
        for start in range(0, x.shape[1], 64):
            z = x[order, start:start+64].astype(float)
            beta = np.zeros(z.shape[1])
            event_sum = z[event].sum(0)
            for _ in range(60):
                eta = np.clip(z*beta, -50, 50)
                weight = np.exp(eta-eta.max(0))
                s0 = np.cumsum(weight[::-1], axis=0)[::-1][starts]
                s1 = np.cumsum((weight*z)[::-1], axis=0)[::-1][starts]
                s2 = np.cumsum((weight*z*z)[::-1], axis=0)[::-1][starts]
                mu = s1/np.maximum(s0, 1e-100)
                grad = event_sum-(deaths[:,None]*mu).sum(0)
                info = (deaths[:,None]*np.maximum(s2/np.maximum(s0,1e-100)-mu*mu,0)).sum(0)
                step = np.clip(grad/np.maximum(info,1e-9), -.5,.5)
                beta += step
                if np.max(np.abs(step))<1e-7:
                    break
            se = 1/np.sqrt(np.maximum(info,1e-9))
            for j in range(z.shape[1]):
                rows.append(dict(feature=features[start+j], beta=beta[j], se=se[j],
                    p=float(2*norm.sf(abs(beta[j]/se[j]))),
                    converged=bool(abs(step[j])<1e-5 and np.isfinite(beta[j]) and abs(beta[j])<15)))
        result = pd.DataFrame(rows)
    result["fdr"] = bh(result.p.where(result.converged, np.nan))
    return result

def make_blocks(c, x, z, pca, graph, score, varying=True):
    soft = graph["soft"]
    hard = np.eye(soft.shape[1])[graph["state"]][:,1:]
    blocks = {
        "clinical": c,
        "clinical_pwas": np.c_[c, score],
        "clinical_pca": np.c_[c, pca],
        "clinical_ae": np.c_[c, z],
        "clinical_elasticnet": np.c_[c, x],
        "panome_neighborhood": np.c_[c, z, graph["neighborhood"]],
        "clinical_tree": np.c_[c, x],
    }
    if soft.shape[1]>1:
        blocks["clinical_state"] = np.c_[c, hard]
        blocks["clinical_soft_state"] = np.c_[c, soft[:,1:]]
        blocks["panome_state"] = np.c_[c, z, soft[:,1:]]
        if varying:
            # Reference-state z coefficient plus state-dependent deviations.
            interactions = np.concatenate([z*soft[:,j:j+1] for j in range(1,soft.shape[1])], axis=1)
            blocks["panome_varying"] = np.c_[c, z, soft[:,1:], interactions]
    if graph["diffusion"].shape[1]:
        blocks["clinical_diffusion"] = np.c_[c, graph["diffusion"]]
    return blocks

class BlockModel:
    def fit(self, block, p, kind, name, a):
        self.kind, self.name = kind, name
        tr, va = p.split.to_numpy()=="train", p.split.to_numpy()=="validation"
        self.scaler = StandardScaler().fit(block[tr])
        scaled = self.scaler.transform(block)
        self.active = np.std(scaled[tr],axis=0)>1e-8
        scaled = scaled[:,self.active]
        if scaled.shape[1] == 0:
            raise ValueError("No variable predictors in this model")
        target = survival_y(p) if kind=="survival" else p.target.to_numpy(float)
        tuning = []
        if name=="clinical_tree":
            if kind=="survival":
                fitted = GradientBoostingSurvivalAnalysis(n_estimators=a.tree_estimators,
                    learning_rate=.05, max_depth=2, min_samples_leaf=20, max_features="sqrt",
                    random_state=a.seed).fit(scaled[tr],target[tr])
            else:
                fitted = HistGradientBoostingRegressor(max_iter=a.tree_estimators,
                    max_leaf_nodes=15, min_samples_leaf=20, l2_regularization=1.,
                    early_stopping=False, random_state=a.seed).fit(scaled[tr],target[tr])
            candidates = [(primary_metric(p.loc[va], fitted.predict(scaled[va]),kind), None,fitted)]
        elif name=="clinical_elasticnet" and kind=="survival":
            fitted = CoxnetSurvivalAnalysis(l1_ratio=.5, n_alphas=a.coxnet_alphas,
                alpha_min_ratio=.01, fit_baseline_model=True, max_iter=100000).fit(scaled[tr],target[tr])
            candidates = [(primary_metric(p.loc[va],fitted.predict(scaled[va],alpha=float(alpha)),kind),
                           float(alpha),fitted) for alpha in fitted.alphas_]
        else:
            grid = np.logspace(-3,1,a.coxnet_alphas) if name=="clinical_elasticnet" else [.1,1.,10.,100.,1000.]
            candidates = []
            for alpha in grid:
                try:
                    if kind=="survival":
                        fitted = CoxPHSurvivalAnalysis(alpha=alpha, ties="breslow", n_iter=200).fit(scaled[tr],target[tr])
                    elif name=="clinical_elasticnet":
                        fitted = ElasticNet(alpha=alpha, l1_ratio=.5,max_iter=30000,
                            random_state=a.seed).fit(scaled[tr],target[tr])
                    else:
                        fitted = Ridge(alpha=alpha).fit(scaled[tr],target[tr])
                    value = primary_metric(p.loc[va],fitted.predict(scaled[va]),kind)
                    if np.isfinite(value):
                        candidates.append((value, float(alpha),fitted))
                except (ValueError, ArithmeticError, np.linalg.LinAlgError):
                    continue
        if not candidates:
            raise ValueError(f"{name}: all training/validation fits failed")
        self.validation_value, self.alpha, self.model = max(candidates,key=lambda row:row[0])
        for value,alpha,_ in candidates:
            tuning.append(dict(model=name, alpha=alpha, validation_metric=value,
                               selection_metric="Harrell_C" if kind=="survival" else "R2"))
        return tuning

    def matrix(self, block):
        return self.scaler.transform(block)[:,self.active]

    def predict(self, block):
        kw = {"alpha":self.alpha} if self.kind=="survival" and self.name=="clinical_elasticnet" else {}
        return self.model.predict(self.matrix(block),**kw)

    def survival_at(self, block, horizon):
        kw = {"alpha":self.alpha} if self.name=="clinical_elasticnet" else {}
        functions = self.model.predict_survival_function(self.matrix(block),**kw)
        return np.array([fn(horizon) for fn in functions])

def run_predictions(x, mask, p, z, pca, transformer, graph, atlas, features, a, out):
    tr, te = p.split.to_numpy()=="train", p.split.to_numpy()=="test"
    clinical = MetadataDesign(words(a.covariates),words(a.categorical)).fit(p.loc[tr])
    c = clinical.transform(p)
    joblib.dump(clinical,out/"clinical_preprocessor.joblib")
    dump(out/"clinical_category_audit.json",clinical.unknown_categories(p.loc[~tr]))
    pwas = marginal_scan(x[tr],p.loc[tr],features,a.outcome_type)
    pwas.to_csv(out/"training_marginal_pwas.csv",index=False)
    selected = pwas.loc[pwas.converged & np.isfinite(pwas.p)].sort_values("p").head(a.pwas_top)
    if selected.empty:
        raise ValueError("No usable training marginal features")
    indices = np.array([features.index(v) for v in selected.feature])
    pw = dict(indices=indices,beta=selected.beta.to_numpy())
    joblib.dump(pw,out/"pwas_score.joblib")
    score = x[:,indices]@pw["beta"]
    support = (graph["soft"][tr]*p.loc[tr,"event"].to_numpy()[:,None]).sum(0) if a.outcome_type=="survival" else graph["soft"][tr].sum(0)
    threshold = a.min_expert_events if a.outcome_type=="survival" else max(30,a.latent+5)
    varying = graph["soft"].shape[1]>1 and np.all(support>=threshold)
    dump(out/"expert_support.json",dict(effective_support=support, required_per_state=threshold,
                                      varying_coefficients_enabled=varying))
    blocks = make_blocks(c,x,z,pca,graph,score,varying)
    if not a.tree:
        blocks.pop("clinical_tree")
    metrics, calibration, limitations, tuning, statuses = [],[],[],[],[]
    test = p.loc[te,[a.id_col]+(["time","event"] if a.outcome_type=="survival" else ["target"])].copy()
    test["state"],test["novel"] = graph["state"][te]+1,graph["novelty"][te]
    predictions = {}
    train_score_sd = {}
    models_index = []
    if not varying:
        statuses.append(dict(model="panome_varying",status="unavailable",
            reason="No multistate partition or insufficient effective events/people per expert"))

    def record(name, all_prediction, survival_callback, artifact, model_type, inputs):
        if not np.isfinite(all_prediction).all():
            raise ValueError("Nonfinite predictions")
        prediction = all_prediction[te]
        rows,cal,limits,risks = evaluate(p.loc[tr],p.loc[te],prediction,survival_callback,
                                        a.outcome_type,a.horizons,name)
        metrics.extend(rows);calibration.extend(cal);limitations.extend(limits)
        predictions[name] = prediction
        train_score_sd[name] = float(np.std(all_prediction[tr]))
        test[f"{name}_prediction"] = prediction
        for horizon,risk in risks.items():
            test[f"{name}_net_risk_{horizon:g}y"] = risk
        statuses.append(dict(model=name,status="completed",reason=""))
        models_index.append(dict(name=name,type=model_type,artifact=artifact,inputs=inputs))

    for name,block in blocks.items():
        log("START",name)
        try:
            fitted = BlockModel()
            tuning.extend(fitted.fit(block,p,a.outcome_type,name,a))
            all_prediction = fitted.predict(block)
            callback = (lambda horizon,fit=fitted,b=block: fit.survival_at(b[te],horizon)) if a.outcome_type=="survival" else None
            record(name,all_prediction,callback,f"{name}.joblib","classical",name)
            joblib.dump(fitted,out/f"{name}.joblib")
            log("DONE",name,f"validation={fitted.validation_value:.4f}")
        except (ValueError,ArithmeticError,np.linalg.LinAlgError) as exc:
            if name=="clinical":
                raise
            statuses.append(dict(model=name,status="failed",reason=str(exc)))
            log("NOTE",name,str(exc))
    if a.neural:
        # The same frozen transformer representation is shared by all three attention controls.
        scaler = StandardScaler().fit(transformer[tr])
        u = scaler.transform(transformer).astype("float32")
        joblib.dump(scaler,out/"transformer_scaler.joblib")
        bank = u[tr].copy()
        np.save(out/"reference_bank.npy",bank)
        ids = p[a.id_col].to_numpy(str)
        groups = p[a.group_col].to_numpy(str) if a.group_col else None
        random = random_neighbors(atlas,ids,groups,a.seed)
        jobs = [("neural_mlp","mlp",x,None,"molecular"),
                ("transformer_column","column",u,None,"transformer"),
                ("panome_transformer","row",u,graph["neighbors"],"neighbors"),
                ("transformer_random","row",u,random,"random")]
        for name,mode,inputs,neighbors,input_name in jobs:
            log("START",name)
            try:
                model = NeuralModel(mode,a.outcome_type,a.token_dim,a.attention_heads)
                model.fit(inputs,c,p,neighbors,bank,a,out,name)
                predicted = model.predict(inputs,c,neighbors,bank,batch=a.batch_size,
                                          return_attention=name=="panome_transformer")
                if name=="panome_transformer":
                    all_prediction,attention = predicted
                    # Every query's attention is over a fixed training-only reference neighborhood.
                    np.save(out/"person_attention.npy",attention)
                else:
                    all_prediction = predicted
                callback = (lambda horizon,m=model,s=all_prediction[te]: m.survival(s,horizon)) if a.outcome_type=="survival" else None
                record(name,all_prediction,callback,f"{name}.pt","neural",input_name)
                tuning.append(dict(model=name,alpha=np.nan,validation_metric=model.best_validation_loss,
                    selection_metric="piecewise_exponential_NLL" if a.outcome_type=="survival" else "standardized_MSE"))
                log("DONE",name,f"validation_loss={model.best_validation_loss:.4f}")
            except (ValueError,ArithmeticError) as exc:
                statuses.append(dict(model=name,status="failed",reason=str(exc)))
                log("NOTE",name,str(exc))
    pd.DataFrame(statuses).to_csv(out/"model_status.csv",index=False)
    pd.DataFrame(metrics).to_csv(out/"test_metrics.csv",index=False)
    pd.DataFrame(calibration,columns=["model","horizon","group","n","predicted","observed","at_risk"]).to_csv(
        out/"test_calibration.csv",index=False)
    pd.DataFrame(tuning).to_csv(out/"validation_tuning.csv",index=False)
    test.to_csv(out/"test_predictions.csv",index=False)
    dump(out/"metric_limitations.json",limitations)
    dump(out/"models.json",models_index)
    dump(out/"training_score_sd.json",train_score_sd)
    paired_uncertainty(p.loc[te],predictions,a.outcome_type,
        p.loc[te,a.group_col].to_numpy(str) if a.group_col else None,
        a.bootstrap,a.seed,out)
    landmarks(p.loc[tr],p.loc[te],predictions,a,out)
    return c
