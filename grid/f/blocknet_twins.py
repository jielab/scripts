"""Prospective trait precision screening using a frozen BlockNet and independent calibration families."""
import hashlib

import joblib
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingRegressor

from abm2 import check_roster
from core2 import aligned, atomic_json, atomic_table, identified, read_table
from output_layout import project_path, trait_path

POLICIES = ("regional_count", "latent_count", "risk_only", "risk_count", "risk_local")
PREFIX = "selected_twin_"


def tie_key(value):
    return int.from_bytes(hashlib.blake2b(str(value).encode(),digest_size=8,person=b'GRID-TWIN').digest(),'big')


def freeze_top(scores, reference, ids, fraction=.20):
    """Exact top fraction in reference, stable threshold for subsequent participants including ties."""
    keys=np.array([tie_key(i) for i in ids],dtype=np.uint64)
    r=np.flatnonzero(reference & np.isfinite(scores))
    if not len(r):
        return np.zeros(len(scores),bool),None
    keep=min(len(r),max(1,int(np.ceil(reference.sum()*fraction))))
    last=r[np.lexsort((keys[r],-scores[r]))[keep-1]]
    value=scores[last]
    selected=(scores>value)|((scores==value)&(keys<=keys[last]))
    return selected,dict(score=float(value),tie_key=int(keys[last]),n_reference=int(reference.sum()),n_eligible=len(r))


def count_twins(z, cohort, donor, thresholds, y, workers, radius_k):
    """Exact TRAIN counts within ancestry; query and donor families are disjoint by cohort contract."""
    groups=cohort.group_id.to_numpy(str)
    pop=cohort.ancestry.to_numpy(str)
    query=~cohort.role.eq("train").to_numpy()
    count=np.full(len(z),np.nan); families=count.copy(); spread=count.copy(); kth=count.copy()
    radii={}
    for ancestry in sorted(set(pop[query])):
        bank=np.flatnonzero(donor & (pop==ancestry))
        ref=np.flatnonzero(thresholds & (pop==ancestry))
        rows=np.flatnonzero(query & (pop==ancestry))
        if len(bank)<radius_k or len(ref)<20:
            radii[ancestry]=dict(status="insufficient_reference",donors=len(bank),reference=len(ref))
            continue
        tree=cKDTree(z[bank])
        reference_dist=tree.query(z[ref],k=[radius_k],workers=workers)[0].ravel()
        radius=max(float(np.median(reference_dist)),1e-10)
        radii[ancestry]=dict(status="available",radius=radius,donors=len(bank),reference=len(ref))
        for first in range(0,len(rows),256):
            batch=rows[first:first+256]
            kth[batch]=tree.query(z[batch],k=[radius_k],workers=workers)[0].ravel()/radius
            neighbors=tree.query_ball_point(z[batch],radius,workers=workers)
            for i,local in zip(batch,neighbors):
                matches=bank[np.asarray(local,int)]
                count[i]=len(matches)
                families[i]=len(np.unique(groups[matches]))
                if len(matches)>=2:
                    spread[i]=np.std(y[matches],ddof=1)
    return count,families,kth,spread,radii


def fit(cfg,root,trait):
    opt=cfg["ldl_twins"]
    folder=trait_path(root,trait,"blocknet.twins");folder.mkdir(parents=True,exist_ok=True)
    cohort=identified(read_table(project_path(root,"tmp/cohort/participants.tsv.gz")))
    check_roster(cohort)
    direct=trait_path(root,trait,"direct.blocknet")
    pred=identified(read_table(direct/"predictions.tsv.gz"))
    index=pd.Index(pred.id).get_indexer(cohort.id)
    embedding=np.load(direct/"embedding.npy",mmap_mode="r")
    if np.any(index<0) or len(embedding)!=len(pred) or len(pred)!=len(cohort):
        raise ValueError("BlockNet embeddings, frozen prediction IDs and cohort differ")
    prediction=pred["direct.blocknet.prediction"].to_numpy(float)[index]
    develop=identified(read_table(project_path(root,f"tmp/cohort/{trait}.development.tsv.gz")))
    cal=cohort.role.eq("calibration").to_numpy()
    fit_half=np.asarray([bool(tie_key("cal/"+g)%2) for g in cohort.group_id])
    train=cohort.role.eq("train").to_numpy()
    permitted=train|(cal&fit_half)
    y=np.full(len(cohort),np.nan)
    y[permitted]=pd.to_numeric(aligned(develop,cohort.id[permitted]).y,errors="raise")
    donor=train&np.isfinite(y)
    if donor.sum()<3:
        raise ValueError("Too few TRAIN donors with observed outcomes")
    if not 0<float(opt["coverage"])<1 or int(opt["dimensions"])<1 or int(opt["radius_k"])<1:
        raise ValueError("Twin dimensions and radius_k must be positive; coverage must be between 0 and 1")
    dims=min(int(opt["dimensions"]),int(donor.sum()-1),embedding.shape[1])
    pca=PCA(n_components=dims,svd_solver="full").fit(embedding[index[donor]])
    z=pca.transform(embedding[index])
    threshold=cal&~fit_half
    n,families,kth,spread,radii=count_twins(z,cohort,donor,threshold,y,cfg["runtime"]["threads"],int(opt["radius_k"]))
    bank_ok=np.isfinite(n)
    abm=identified(read_table(trait_path(root,trait,"abm/abm.predictions.tsv.gz")))
    region_count=np.full(len(y),np.nan);query=~train
    region_count[query]=pd.to_numeric(aligned(abm,cohort.id[query]).n_support,errors="raise")
    raw=cohort[["age","sex","PC1","PC2","PC3","PC4"]].apply(pd.to_numeric,errors="raise").to_numpy(float)
    basic=np.column_stack([prediction,raw])
    counts=np.column_stack([basic,np.log1p(n),np.log1p(families)])
    local=np.column_stack([counts,kth,spread])
    extended={"risk_only":basic,"risk_count":counts,"risk_local":local}
    fit_rows=cal&fit_half&bank_ok&np.isfinite(y)
    scores={"regional_count":region_count,"latent_count":n}
    heads={}
    result=cohort[["id","role","ancestry"]].copy()
    result["twin_count"]=n;result["twin_families"]=families
    result["twin_kth_distance_ratio"]=kth;result["twin_train_sd"]=spread
    for name,matrix in extended.items():
        if fit_rows.sum()<2:
            scores[name]=np.full(len(y),np.nan)
            result[name+".twin_expected_loss"]=scores[name]
            continue
        head=HistGradientBoostingRegressor(max_iter=80,max_leaf_nodes=7,min_samples_leaf=max(5,min(150,int(fit_rows.sum()/30))),
            l2_regularization=10.,early_stopping=False,random_state=cfg["cohort"]["seed"])
        head.fit(matrix[fit_rows],np.square(y[fit_rows]-prediction[fit_rows]))
        loss=np.full(len(y),np.nan)
        loss[bank_ok]=np.maximum(0,head.predict(matrix[bank_ok]))
        result[name+".twin_expected_loss"]=loss
        scores[name]=-loss;heads[name]=head
    cuts={}
    for name,values in scores.items():
        selected=np.zeros(len(y),bool);cuts[name]={}
        for ancestry,info in radii.items():
            if info["status"]!="available":continue
            member=bank_ok&cohort.ancestry.eq(ancestry).to_numpy()
            eligible=member&np.isfinite(values)
            if name in ("regional_count","latent_count"):
                eligible &= values>0
            # A zero count cannot be advertised as a matched twin. Coverage is always reported.
            pool=threshold&member
            scores_for_cut=np.where(eligible,values,-np.inf)
            mask,record=freeze_top(scores_for_cut,pool,cohort.id.to_numpy(),opt["coverage"])
            selected|=mask&eligible
            cuts[name][ancestry]=record
        result[PREFIX+name]=selected
    # Use the independent policy-calibration predictions, never TEST outcomes, to define risk cells.
    reference_predictions=prediction[threshold]
    edges=(np.unique(np.quantile(reference_predictions,np.arange(.1,1.,.1)))
           if reference_predictions.size else np.array([],dtype=float))
    result["twin_risk_bin"]=np.searchsorted(edges,prediction,side="right")
    atomic_table(result,folder/"predictions.tsv.gz")
    joblib.dump(dict(pca=pca,error_models=heads,donor_z=z[donor],donor_id=cohort.id.to_numpy()[donor],
                     donor_family=cohort.group_id.to_numpy()[donor],donor_ancestry=cohort.ancestry.to_numpy()[donor],
                     donor_outcome=y[donor]),folder/"model.joblib",compress=3)
    atomic_json(dict(trait=trait,radii=radii,thresholds=cuts,risk_edges=edges.tolist(),dimensions=dims,radius_k=int(opt["radius_k"]),
                     coverage=float(opt["coverage"]),error_fit_n=int(fit_rows.sum()),threshold_n=int(threshold.sum()),
                     head_features={"risk_only":["BlockNet prediction","age","sex","PC1","PC2","PC3","PC4"],
                                    "risk_count":"risk_only + log donor count + log donor family count",
                                    "risk_local":"risk_count + relative kth distance + donor outcome SD"},
                     error_fit="calibration families (hash partition 1), fixed direct.blocknet squared error",
                     threshold_role="disjoint calibration families (partition 0), features only",
                     count_definition="TRAIN donors with observed outcome in the same ancestry, within the frozen embedding radius",
                     test_outcomes_used=False),folder/"model.json")
