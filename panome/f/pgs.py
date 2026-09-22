"""Optional PGS overlay. It does not change atlas fitting or primary prediction."""
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score
from io_data import read_table, numeric
from common import dump, log

def pgs_overlay(x, observed, p, clinical, features, a, out):
    if not a.pgs:
        dump(out/"pgs_status.json",dict(status="not_requested"))
        return
    scores = read_table(a.pgs_file,a.id_col,r_bin=a.r_bin)
    canonical = {}
    for col in scores:
        key = col.upper() if a.biom=="prot" else col
        if key in canonical:
            raise ValueError("PGS column names collide after harmonization")
        canonical[key] = col
    aligned = scores.set_index(a.id_col).reindex(p[a.id_col])
    train = p.split.to_numpy()=="train"
    test = p.split.to_numpy()=="test"
    mapping,rows = [],[]
    residual_sum = np.zeros(len(p))
    residual_count = np.zeros(len(p),dtype=int)
    for j,feature in enumerate(features):
        expected = feature+".pgs"
        key = expected.upper() if a.biom=="prot" else expected
        if key not in canonical:
            continue
        column = canonical[key]
        raw = aligned[column]
        g = pd.to_numeric(raw,errors="coerce").to_numpy(float)
        if (raw.notna().to_numpy() & np.isnan(g)).any():
            raise ValueError(f"Non-numeric PGS: {column}")
        available = np.isfinite(g) & observed[:,j]
        tr,te = train & available,test & available
        mapping.append(dict(feature=feature,score_column=column,n_train=int(tr.sum()),n_test=int(te.sum())))
        if tr.sum()<max(50,clinical.shape[1]+10) or te.sum()<20 or np.std(g[tr])<1e-8:
            continue
        mean,scale = np.mean(g[tr]),np.std(g[tr])
        gs = (g-mean)/scale
        c = np.c_[np.ones(len(p)),clinical]
        base_beta = np.linalg.lstsq(c[tr],x[tr,j],rcond=None)[0]
        full = np.c_[c,gs]
        beta = np.linalg.lstsq(full[tr],x[tr,j],rcond=None)[0]
        prediction = full@beta
        base = c@base_beta
        complete = np.isfinite(prediction)&available
        residual_sum[complete] += np.abs(x[complete,j]-prediction[complete])
        residual_count[complete] += 1
        rows.append(dict(feature=feature,score_column=column,n_train=int(tr.sum()),
            n_test=int(te.sum()),beta_pgs_per_train_sd=float(beta[-1]),
            test_R2_clinical=float(r2_score(x[te,j],base[te])),
            test_R2_clinical_pgs=float(r2_score(x[te,j],prediction[te])),
            test_delta_R2=float(r2_score(x[te,j],prediction[te])-r2_score(x[te,j],base[te]))))
    pd.DataFrame(mapping,columns=["feature","score_column","n_train","n_test"]).to_csv(
        out/"pgs_exact_mapping.csv",index=False)
    pd.DataFrame(rows,columns=["feature","score_column","n_train","n_test","beta_pgs_per_train_sd",
        "test_R2_clinical","test_R2_clinical_pgs","test_delta_R2"]).to_csv(out/"pgs_reconstruction.csv",index=False)
    result = p[[a.id_col,"split"]].copy()
    result["n_observed_matched_features"] = residual_count
    result["mean_absolute_unexplained_component"] = np.divide(residual_sum,residual_count,
        out=np.full(len(p),np.nan),where=residual_count>0)
    result.to_csv(out/"person_pgs_discordance.csv",index=False)
    dump(out/"pgs_status.json",dict(status="completed" if rows else "no_estimable_matches",
        source=a.pgs_source, discovery_overlap=a.pgs_overlap, source_file=a.pgs_file,
        evaluated_features=len(rows),
        interpretation="PGS predicts an adult molecular measurement conditional on clinical covariates. "
        "It is not a biomarker measured at birth, a causal component, or a genetic/environmental partition. "
        "GWAS overlap is not removed by this train/test split; known/unknown overlap limits interpretation."))
    log("DONE","PGS overlay",f"{len(rows)} matched features; discovery_overlap={a.pgs_overlap}")
