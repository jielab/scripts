#!/usr/bin/env python3
"""Synthetic selective-training regression suite from the 2026-10-04 delivery.
Imports the actual integrated c1.abm.py; temporary fits stay in /tmp.
Run: python validation/check_c1.py
"""
from __future__ import annotations
import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
import sys
import tempfile
import types
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
spec=importlib.util.spec_from_file_location("c1_abm",Path(__file__).resolve().parents[1]/"f/c1.abm.py")
m=importlib.util.module_from_spec(spec);sys.modules["c1_abm"]=m;spec.loader.exec_module(m)


def args(**changes):
    p=argparse.ArgumentParser();m.s6_options(p);d=vars(p.parse_args([]))
    d.update(seed=2026,trait="synthetic",biom="prot",id_col="eid",group_col="family",covariates="age,sex",
      residualize="",categorical="sex",transform="none",feature_missing=.2,sample_missing=.2,
      horizon=5,min_censor_survival=.05,teacher_c=.1,max_iter=2000,tree="hist",tree_estimators=20,
      cores=2,c_grid=[.01,.1],bootstrap=20,split_file="",module_file="",resume=False,
      selective_folds=3,selective_gate_trees=20,selective_components=6,selective_neighbors=30,
      selective_min_events=5,selective_coverages="0.4,0.6,1.0",shuffle_development_outcomes=False)
    d.update(changes)
    return m.s6_validate(types.SimpleNamespace(**d))


def data(n=2400,seed=315):
    rng=np.random.default_rng(seed)
    x=rng.normal(size=(n,12)).astype("float32")
    age=rng.uniform(40,70,n);sex=rng.integers(0,2,n)
    rate=.025*np.exp(.035*(age-55)+.5*sex+.9*x[:,0]-.6*x[:,1])
    t=rng.exponential(1/rate);c=rng.uniform(6,10,n)
    p=pd.DataFrame(dict(eid=["p"+str(i) for i in range(n)],family=["f"+str(i//2) for i in range(n)],
         age=age,sex=sex,time=np.minimum(t,c),event=(t<=c).astype(int)))
    return x,p,["X"+str(i) for i in range(x.shape[1])]


def run():
    passed=[]
    def done(name):passed.append(name);print("PASS "+name,flush=True)
    a=args();x,p,features=data()
    groups=m.s6_groups(p,a)
    fold=m.s6_group_folds(groups,3,90)
    for g in np.unique(groups):assert len(set(fold[groups==g]))==1
    shuffled=np.random.default_rng(2).permutation(len(p))
    assert np.array_equal(fold[shuffled],m.s6_group_folds(groups[shuffled],3,90))
    done("outcome-blind, family-preserving and row-order-stable folds")
    try:args(covariates="age,event")
    except ValueError:pass
    else:raise AssertionError("Outcome predictor was accepted")
    done("outcome/follow-up covariates rejected")
    ids=np.array([f"a{i}" for i in range(100)])
    score=np.zeros(100);support=np.linspace(0,1,100)
    rule=m.S6CoverageRule().fit(score,support,ids,.6,.99,7)
    assert rule.apply(score,support,ids).sum()==60
    full=m.S6CoverageRule().fit(score,support,ids,1.,.1,7)
    assert full.apply(score,support,ids).all()
    repeat=rule.apply(score[shuffled[:100]%100],support[shuffled[:100]%100],ids[shuffled[:100]%100])
    assert np.array_equal(repeat,rule.apply(score,support,ids)[shuffled[:100]%100])
    done("constant-score ties and true 100-percent coverage endpoint")
    # Pure probability borrowing must fall back instead of using a family donor.
    small_a=args(selective_neighbors=5)
    bank=m.S6Neighborhood().fit(x[:100],np.tile([0,1],50),np.ones(100),np.array(["same"]*100),small_a,3,features)
    prob,diag=bank.predict(x[100:105],np.full(5,.2),np.array(["same"]*5))
    assert np.allclose(prob,.2) and np.all(diag[:,2]==0)
    done("same-family donors excluded; empty neighborhood uses prior")
    prep=m.S6Preprocessor(.2,[],["sex"],"none").fit(x[:1000],p.iloc[:1000])
    mean=prep.mean.copy();median=prep.median.copy()
    prep.transform(x[1000:1010]*1e6,p.iloc[1000:1010])
    assert np.array_equal(mean,prep.mean) and np.array_equal(median,prep.median)
    done("test transformations cannot change fitted preprocessing")
    # Nontrivial censoring and sparse plate correction exercise UKB-specific paths.
    mixed=p.copy();rr=np.random.default_rng(55)
    early=rr.choice(len(p),len(p)//4,replace=False)
    mixed.loc[early,"time"]=rr.uniform(.1,3.,len(early));mixed.loc[early,"event"]=0
    censor=m.S6Censoring().fit(mixed,a.horizon,a.min_censor_survival)
    yy,ww=censor.labels_weights(mixed)
    assert (ww[early]==0).all() and (yy[early]==0).all()
    assert ww.max()>1 and np.isfinite(ww).all()
    assert censor.at([a.horizon])[0]<1
    done("pre-horizon censoring is not mislabeled as an observed control")
    technical=p.iloc[:600].copy();technical["plate"]=np.arange(600)%60
    z=x[:600].copy();z[np.arange(60),np.arange(60)%12]=np.nan
    corrected=m.S6Preprocessor(.2,["plate"],["sex","plate"],"none").fit(z,technical)
    zz,observed=corrected.transform(z,technical)
    assert np.isfinite(zz).all() and (~observed).sum()==60
    assert m.sparse.issparse(corrected.design.transform(technical))
    done("sparse technical residualization and missing molecular measurements")
    # Uno-style implementation checked against its explicit comparable-pair formula.
    km=m.S6Censoring().fit(p,5,.05)
    rng=np.random.default_rng(22);t=rng.integers(1,9,45).astype(float);e=rng.integers(0,2,45);risk=rng.integers(0,5,45)/5
    num=den=0.
    for i in range(len(t)):
        if e[i] and t[i]<=km.horizon:
            wi=1/km.at([t[i]],left=True)[0]**2
            for j in range(len(t)):
                if t[j]>t[i]:
                    den+=wi;num+=wi*(float(risk[i]>risk[j])+.5*float(risk[i]==risk[j]))
    assert abs(m.s6_uno(t,e,risk,km)-num/den)<1e-12
    assert abs(m.s6_uno(t,e,np.ones(len(t)),km)-.5)<1e-12
    done("IPCW concordance equals naive pair calculation; ties score 0.5")
    with tempfile.TemporaryDirectory(prefix="c1-selective-check-") as tmp:
        tmp=Path(tmp)
        module_file=tmp/"modules.tsv"
        pd.DataFrame({"feature":features[:6],"module":["A"]*3+["B"]*3}).to_csv(module_file,sep="\t",index=False)
        module_a=args(module_file=str(module_file))
        bank=m.S6Neighborhood().fit(x[:100],np.tile([0,1],50),np.ones(100),groups[:100],module_a,4,features)
        assert bank.transform(x[:3]).shape[1]==8 and bank.module_names==["A","B"]
        done("predeclared molecular modules included in geometry")
        # Flip the outcomes of one entire recipient fold. Its own raw gate scores,
        # hard masks and soft weights must be bit-for-bit unchanged.
        cf=args(tree="none",selective_gate_trees=15)
        base_x,base_p=x[:1800],p.iloc[:1800].copy()
        held=m.s6_group_folds(m.s6_groups(base_p,cf),cf.selective_folds,cf.seed+1200)==0
        _,v1=m.s6_crossfit_gate(base_x,base_p,features,cf)
        changed=base_p.copy();changed.loc[held,"event"]=1-changed.loc[held,"event"]
        _,v2=m.s6_crossfit_gate(base_x,changed,features,cf)
        for key in ["omics_gain","local_gain","training_selected","training_soft"]:
            assert np.array_equal(v1[key][held],v2[key][held]),key
        done("own-label poisoning cannot change honest fold-local training weights")
        keep,soft=m.s6_train_weights(v1,base_p.eid.to_numpy(),cf)
        rnd,rw=m.s6_matched_random(keep,soft,v1["OOF_clinical"],base_p.eid.to_numpy(),2)
        assert rnd.sum()==keep.sum() and np.allclose(np.sort(rw),np.sort(soft))
        assert np.min(soft)>=cf.selective_weight_floor and np.max(soft)<=1
        done("matched random sample size/weight distribution; bounded soft weights")
        # End-to-end: all fitted arms, calibration, frozen prediction, metrics, PNGs.
        p=p.copy();p["split"]=m.s6_split(p,a)
        out=tmp/"first";b=m.s6_fit(x,p,features,a,out);result=m.s6_evaluate(out)
        table=pd.read_csv(out/"test_individuals.csv")
        for _,sub in result.groupby("subset"):
            assert sub.n.nunique()==1 and sub.cases.nunique()==1
        assert {"AUC_IPCW","Uno_C_horizon","Brier_IPCW"}<=set(result)
        assert (out/"Fig_selective_training.png").is_file()
        done("end-to-end survival fitting, calibration, same-mask comparisons and figures")
        # Real inference does not need outcome/follow-up columns.
        test=np.flatnonzero(p.split.eq("test"))
        projection=m.s6_predict(b,x[test],p.iloc[test].drop(columns=["time","event"]))
        assert np.allclose(projection[b["primary"]],table[b["primary"]],atol=1e-12)
        done("serialized model reproduces outcome-free inference")
        # Refit with a permuted TEST outcome file and unchanged frozen split.
        p2=p.copy();perm=np.random.default_rng(200).permutation(test)
        p2.loc[test,["time","event"]]=p2.loc[perm,["time","event"]].to_numpy()
        b2=m.s6_fit(x,p2,features,a,tmp/"poisoned_test")
        t2=pd.read_csv(tmp/"poisoned_test"/"test_individuals.csv")
        for col in list(b["model"].calibrators)+["omics_gain","local_gain","supported_match","prediction_released"]:
            assert np.allclose(table[col],t2[col],equal_nan=True,atol=1e-12),col
        done("test-label permutation leaves every prediction and selection mask unchanged")
        before=hashlib.sha256((out/"model_bundle.joblib").read_bytes()).hexdigest()
        def forbidden_fit(*_a,**_k):raise AssertionError("Evaluation attempted model fitting")
        with patch.object(m.S6Calibrator,"fit",forbidden_fit),patch.object(m.S6Preprocessor,"fit",forbidden_fit),patch.object(m.S6GainGate,"fit",forbidden_fit):
            m.s6_evaluate(out)
        assert before==hashlib.sha256((out/"model_bundle.joblib").read_bytes()).hexdigest()
        done("evaluation is fit-free and cannot mutate frozen model")
        if importlib.util.find_spec("lightgbm"):
            la=args(tree="lightgbm")
            learner=m.s6_fit_tree(m.s6_make_tree(la,4),x[:300],(x[:300,0]>0).astype(int),np.ones(300))
            assert np.all(np.isfinite(learner.predict_proba(x[300:310])))
            done("installed LightGBM branch trains and predicts")
        else:print("SKIP LightGBM not installed")
    assert m.reference_parser().parse_args([]).abm_design=="selective"
    done("actual integrated host import and selective CLI")
    print(f"\n{len(passed)} checks passed. Synthetic data only; no UKB training and no R model execution.")
    return passed


if __name__=="__main__":
    with threadpool_limits(2):run()
