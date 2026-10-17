"""GRID genetic-profile twin predictor. No query outcome is accepted by project().

Own build outcomes score *out-of-fold* donor consistency. This is supervised
reference selection, not evidence that discarded observations are erroneous.
Continuous traits and baseline binary traits are supported; not survival data.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import hashlib
import json
import math
import warnings
from typing import Any
import numpy as np
from scipy.special import expit, logit
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold, StratifiedGroupKFold
from sklearn.neighbors import KDTree
from sklearn.decomposition import PCA
from sklearn.exceptions import ConvergenceWarning

VERSION = "GRID-TWIN-2026-10-14"
ROLE_NAMES = ("build", "tune_model", "tune_gate", "calibration_fit", "audit")


def hash_number(value, seed=0):
    return int.from_bytes(hashlib.sha256(f"{seed}:{value}".encode()).digest()[:8], "big")


def id_hash(ids):
    return hashlib.sha256("\n".join(sorted(map(str, ids))).encode()).hexdigest()


def state_hash(obj):
    """Canonical fitted-state hash, invariant to pickle aliasing and tree counters."""
    def encode(value):
        h=hashlib.sha256()
        def put(x):h.update(x if isinstance(x,bytes) else str(x).encode())
        if value is None:put("none")
        elif isinstance(value,np.generic):return encode(value.item())
        elif isinstance(value,np.ndarray):
            put("ndarray");put(value.dtype.str);put(repr(value.shape))
            if value.dtype.names:
                # Structured tree nodes may contain nonsemantic padding bytes.
                for name in value.dtype.names:put(name);put(encode(value[name]))
            elif value.dtype.hasobject:
                for x in value.ravel(order="C"):put(encode(x))
            else:put(np.ascontiguousarray(value).tobytes())
        elif isinstance(value,(str,bytes,int,float,bool)):
            put(type(value).__name__);put(repr(value))
        elif isinstance(value,dict):
            put("dict")
            for key in sorted(value,key=lambda x:(type(x).__name__,repr(x))):put(encode(key));put(encode(value[key]))
        elif isinstance(value,(list,tuple)):
            put(type(value).__name__)
            for x in value:put(encode(x))
        elif isinstance(value,(set,frozenset)):
            for x in sorted(encode(x) for x in value):put(x)
        elif isinstance(value,type):put(value.__module__+"."+value.__qualname__)
        elif isinstance(value,np.dtype):put(repr(value.descr))
        elif hasattr(value,"__dict__"):
            put(type(value).__module__+"."+type(value).__qualname__)
            fields=dict(vars(value))
            # The KDTree is deterministically reconstructed from ReferenceBank.z;
            # its query counters are not part of the fitted scientific model.
            if isinstance(value,ReferenceBank):fields.pop("tree",None)
            put(encode(fields))
        elif hasattr(value,"__getstate__"):
            put(type(value).__module__+"."+type(value).__qualname__);put(encode(value.__getstate__()))
        else:raise TypeError("Cannot fingerprint fitted-state type: "+str(type(value)))
        return h.digest()
    return encode(obj).hex()


@dataclass(frozen=True)
class Config:
    trait: str = "height"
    trait_type: str = "continuous"
    seed: int = 20260904
    folds: int = 5
    repeats: int = 2
    role_fractions: tuple = (.60, .15, .15, .05, .05)
    keep_fraction: float = .99
    pruning: str = "profile_concordance"
    soft_floor: float = .10
    primary: str = "ABM_clean"
    selector: str = "support_quality"
    primary_coverage: float = .50
    coverages: tuple = (.20, .40, .50, .60, .80, 1.)
    ridge_grid: tuple = (1., 10., 100., 1000.)
    logistic_c_grid: tuple = (.01, .1, 1.)
    pilot_ridge: float = 100.
    pilot_c: float = .1
    max_iter: int = 3000
    max_genetic_missing: float = .20
    max_covariate_missing: float = .20
    max_feature_missing: float = .30
    profile_kind: str = "genetic_summary"
    max_profile_dims: int = 32
    k: int = 64
    min_twins: int = 20
    min_group_ess: float = 15.
    radius_quantile: float = .90
    query_batch: int = 256
    gate_iterations: int = 100
    min_role_n: int = 30
    min_class_n: int = 10
    bootstrap: int = 1000
    grouping_status: str = "individual_unverified"
    upstream_gwas_overlap: str = "unknown"
    label_budget_note: str = "build fit; tune selection; separate probability calibration and audit"

    def validate(self):
        if self.trait_type not in {"continuous", "binary"}:
            raise ValueError("Only continuous or baseline binary outcomes are supported")
        if self.profile_kind not in {"genetic_summary", "variant_dosage", "provider_defined"}:
            raise ValueError("Invalid explicitly declared genetic profile kind")
        if len(self.role_fractions) != 5 or min(self.role_fractions) <= 0 or not np.isclose(sum(self.role_fractions), 1):
            raise ValueError("Five positive role fractions summing to one required")
        if not 0 < self.keep_fraction <= 1 or not 0 < self.soft_floor <= 1:
            raise ValueError("Invalid donor retention/weight")
        if self.pruning not in {"profile_concordance", "class_risk_stratified"}:
            raise ValueError("Unknown pruning mode")
        if self.selector not in {"support_quality", "expected_error", "support_only"}:
            raise ValueError("Unknown selector")
        if self.primary not in {"ABM_clean", "ABM_soft", "ABM_full"}:
            raise ValueError("Primary must be a declared ABM training arm")
        if self.primary_coverage not in self.coverages or any(not 0 < q <= 1 for q in self.coverages):
            raise ValueError("Primary coverage must be predeclared")
        if len(set(self.coverages)) != len(self.coverages):
            raise ValueError("Duplicate coverage levels")
        for value in (self.max_genetic_missing, self.max_covariate_missing):
            if not 0 <= value < 1: raise ValueError("Invalid sample missingness threshold")
        if not 0 < self.max_feature_missing < 1 or not 0 < self.radius_quantile < 1:
            raise ValueError("Invalid feature/radius quantile threshold")
        for name in ("folds", "repeats", "k", "min_twins", "query_batch", "min_role_n", "min_class_n", "max_iter", "max_profile_dims", "gate_iterations", "bootstrap"):
            if int(getattr(self,name)) != getattr(self,name) or getattr(self,name) < 1:
                raise ValueError(f"Invalid positive integer {name}")
        if self.folds < 2 or self.min_twins > self.k or not 1 <= self.min_group_ess <= self.k:
            raise ValueError("Invalid folds/neighbour support")
        if self.min_class_n < 2: raise ValueError("Need at least two examples of each class")
        if min(*self.ridge_grid, *self.logistic_c_grid, self.pilot_ridge, self.pilot_c) <= 0:
            raise ValueError("Model penalties must be positive")


def validate_ids(ids, groups=None):
    ids = np.asarray(ids).astype(str)
    groups = ids.copy() if groups is None else np.asarray(groups).astype(str)
    bad = {"", "nan", "none", "na", "<na>"}
    if len(ids) != len(groups) or len(set(ids)) != len(ids):
        raise ValueError("Unique IDs and aligned group labels required")
    if any(v.strip().lower() in bad for v in np.r_[ids,groups]):
        raise ValueError("Blank IDs/group labels are not allowed")
    return ids,groups


def assign_roles(groups, config):
    """Freeze role membership by group ID, not outcomes or measured fit quality."""
    edges = np.cumsum(config.role_fractions)
    values = np.array([hash_number(g,config.seed+19)/2**64 for g in groups])
    ix = np.minimum(np.searchsorted(edges,values,side="right"),4)
    return np.asarray(ROLE_NAMES)[ix]


def assert_roles(groups, roles):
    seen = {}
    for g,r in zip(groups,roles):
        if g in seen and seen[g] != r: raise ValueError("Reference group crosses roles")
        seen[g]=r
    if set(roles) != set(ROLE_NAMES): raise ValueError("All five development roles are required; never pass test rows")


def assert_labels(y, config, minimum=None):
    y=np.asarray(y,float)
    if not np.isfinite(y).all() or len(y)<(minimum or config.min_role_n):
        raise ValueError("Insufficient observed development outcomes")
    if config.trait_type == "binary":
        if not np.isin(y,[0,1]).all(): raise ValueError("Binary y must be 0/1, not incident survival times")
        if min(np.sum(y==v) for v in (0,1)) < config.min_class_n:
            raise ValueError("Too few cases/controls after declared QC/pruning; do not silently restore donors or change the split")
    elif np.var(y) < 1e-12: raise ValueError("Constant development outcome")


class NumericProfile:
    """Frozen preprocessing with genetic-only QC, no clinical dummy dilution."""
    def fit(self, X, genetic_indices, covariate_indices, config):
        X=np.asarray(X,float)
        self.width=X.shape[1];self.config=config
        self.genetic=np.asarray(genetic_indices,int);self.clinical=np.asarray(covariate_indices,int)
        if not len(self.genetic) or len(set(self.genetic))!=len(self.genetic) or len(set(self.clinical))!=len(self.clinical) or set(np.r_[self.genetic,self.clinical]) != set(range(self.width)):
            raise ValueError("Explicit genetic/covariate identities covering the schema are required")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",RuntimeWarning)
            variance=np.nanvar(np.where(np.isfinite(X),X,np.nan),axis=0)
        keep=(np.mean(~np.isfinite(X),axis=0)<config.max_feature_missing)&(variance>1e-12)
        self.genetic_qc=self.genetic[keep[self.genetic]]
        if not len(self.genetic_qc): raise ValueError("No usable variable genetic profile")
        self.keep=np.flatnonzero(keep)
        qc=self.qc(X)
        if qc.sum()<config.min_role_n: raise ValueError("Too few genetic-QC eligible training subjects")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",RuntimeWarning)
            var=np.nanvar(X[qc][:,self.keep],axis=0)
        self.keep=self.keep[np.isfinite(var)&(var>1e-12)]
        self.profile_keep=self.genetic_qc[np.isin(self.genetic_qc,self.keep)]
        if not len(self.profile_keep): raise ValueError("Genetic profile collapsed after sample QC")
        self.profile_positions=np.array([list(self.keep).index(j) for j in self.profile_keep],int)
        self.median=np.nanmedian(X[qc][:,self.keep],axis=0)
        if not np.isfinite(self.median).all(): raise ValueError("Entirely missing fitted predictor")
        z=np.where(np.isfinite(X[qc][:,self.keep]),X[qc][:,self.keep],self.median)
        self.mean=z.mean(0);self.scale=np.maximum(z.std(0),1e-8)
        return self
    def missing(self,X):
        X=np.asarray(X,float)
        g=np.mean(~np.isfinite(X[:,self.genetic_qc]),axis=1)
        c=np.mean(~np.isfinite(X[:,self.clinical]),axis=1) if len(self.clinical) else np.zeros(len(X))
        return g,c
    def qc(self,X):
        X=np.asarray(X,float)
        if X.ndim!=2 or X.shape[1]!=self.width:raise ValueError("Numeric input dimensions differ from frozen schema")
        g,c=self.missing(X)
        ok=(g<=self.config.max_genetic_missing)&(c<=self.config.max_covariate_missing)&~np.isinf(X).any(1)
        if self.config.profile_kind=="variant_dosage":
            a=X[:,self.genetic];bad=np.isfinite(a)&((a<0)|(a>2))
            ok &= ~bad.any(1)
        return ok
    def transform(self,X):
        X=np.asarray(X,float)
        if X.ndim!=2 or X.shape[1]!=self.width:raise ValueError("Input schema mismatch")
        z=X[:,self.keep]
        return (np.where(np.isfinite(z),z,self.median)-self.mean)/self.scale


class Predictor:
    def fit(self,z,y,weights,config,penalty):
        ok=np.asarray(weights)>0;assert_labels(np.asarray(y)[ok],config)
        self.binary=config.trait_type=="binary"
        m=LogisticRegression(C=float(penalty),solver="lbfgs",max_iter=config.max_iter,random_state=config.seed) if self.binary else Ridge(alpha=float(penalty))
        with warnings.catch_warnings():
            warnings.simplefilter("error",ConvergenceWarning)
            m.fit(np.asarray(z)[ok],np.asarray(y)[ok],sample_weight=np.asarray(weights)[ok]/np.mean(np.asarray(weights)[ok]))
        self.model=m;return self
    def predict(self,z):
        return self.model.predict_proba(z)[:,1] if self.binary else self.model.predict(z)
    @property
    def coefficients(self):return np.ravel(self.model.coef_)


class Calibration:
    def fit(self,p,y,binary):
        self.binary=binary
        if binary:
            self.model=LogisticRegression(C=1e6,solver="lbfgs",max_iter=2000).fit(logit(np.clip(p,1e-7,1-1e-7))[:,None],y)
        else:
            # Identity in outcome units. Do not use audit/test Y to recalibrate.
            self.model=None
        return self
    def predict(self,p):
        return self.model.predict_proba(logit(np.clip(p,1e-7,1-1e-7))[:,None])[:,1] if self.binary else np.asarray(p)


def loss(y,p,binary):
    if not binary:return (np.asarray(y)-p)**2
    p=np.clip(p,1e-7,1-1e-7)
    return -(np.asarray(y)*np.log(p)+(1-np.asarray(y))*np.log1p(-p))


def oof_quality(X,y,ids,groups,gidx,cidx,config):
    n=len(y);pred=np.full((config.repeats,n),np.nan);error=pred.copy();quality_qc=np.zeros((config.repeats,n),bool);audit=[]
    for repeat in range(config.repeats):
        if config.trait_type=="binary":
            split=StratifiedGroupKFold(config.folds,shuffle=True,random_state=config.seed+repeat)
            folds=split.split(X,y,groups)
        else:
            # Relabeling groups creates reproducible shuffled GroupKFold in older sklearn too.
            ug=np.unique(groups);rng=np.random.default_rng(config.seed+repeat)
            mapping=dict(zip(ug,rng.permutation(len(ug))))
            shuffled=np.array([mapping[g] for g in groups])
            folds=GroupKFold(config.folds).split(X,groups=shuffled)
        for fold,(tr,va) in enumerate(folds):
            if set(groups[tr])&set(groups[va]):raise AssertionError("OOF group leakage")
            prep=NumericProfile().fit(X[tr],gidx,cidx,config)
            tq=prep.qc(X[tr]);vq=prep.qc(X[va]);z=prep.transform(X[tr])
            penalty=config.pilot_c if config.trait_type=="binary" else config.pilot_ridge
            m=Predictor().fit(z,y[tr],tq.astype(float),config,penalty)
            pv=m.predict(prep.transform(X[va]));pred[repeat,va]=pv
            scale=1. if config.trait_type=="binary" else max(np.var(y[tr][tq]),1e-8)
            error[repeat,va]=np.where(vq,loss(y[va],pv,config.trait_type=="binary")/scale,np.nan)
            quality_qc[repeat,va]=vq
            audit.append(dict(repeat=repeat,fold=fold,training_N=len(tr),training_QC_excluded=int((~tq).sum()),recipient_N=len(va),recipient_QC_excluded=int((~vq).sum()),recipient_id_hash=id_hash(ids[va]),recipient_groups_absent_from_teacher=True))
    if not np.isfinite(pred).all():raise ValueError("Incomplete OOF teacher predictions")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning);err=np.nanmedian(error,axis=0)
    ok=quality_qc.all(0);err[~ok]=np.nan
    return dict(prediction=pred.mean(0),prediction_sd=pred.std(0),loss=err,technical_qc=ok,audit=audit)


def prune(y,quality,ids,final_qc,config):
    valid=quality['technical_qc']&final_qc&np.isfinite(quality['loss'])
    if not valid.any():raise ValueError("No valid OOF donors")
    edges=np.unique(np.quantile(quality['prediction'][valid],np.linspace(0,1,6)))[1:-1]
    bins=np.searchsorted(edges,quality['prediction'])
    controls=bins+(6*np.asarray(y,int) if config.trait_type=="binary" else 0)
    strata=controls if config.pruning=="class_risk_stratified" else np.zeros(len(y),int)
    keep=np.zeros(len(y),bool);soft=np.zeros(len(y));random_keep=keep.copy()
    ties=np.array([hash_number(i,config.seed+127) for i in ids],np.uint64)
    for s in np.unique(strata[valid]):
        ix=np.flatnonzero(valid&(strata==s));order=ix[np.lexsort((ties[ix],quality['loss'][ix]))]
        keep[order[:math.ceil(config.keep_fraction*len(order))]]=True
        ranks=(np.arange(len(order))+.5)/len(order)
        soft[order]=config.soft_floor+(1-config.soft_floor)*(1-ranks)
    for s in np.unique(controls[valid]):
        ix=np.flatnonzero(valid&(controls==s));ix=ix[np.argsort(ties[ix],kind='stable')]
        random_keep[ix[:int(keep[ix].sum())]]=True
    assert_labels(np.asarray(y)[keep],config)
    return dict(eligible=valid,keep=keep,soft=soft,random_keep=random_keep)


class Geometry:
    """Vector similarity, with a positive weight for every retained genetic feature.
    Optional high-dimensional projection retains a reconstruction-residual coordinate.
    This is not a full-SNP claim when the supplied features are aggregate scores.
    """
    def fit(self,z,coefficients,config):
        self.weights=np.sqrt(.25+.75*np.abs(coefficients)/max(np.mean(np.abs(coefficients)),1e-8))
        a=z*self.weights;self.pca=None;self.residual_scale=None
        if a.shape[1]>config.max_profile_dims:
            self.pca=PCA(n_components=min(config.max_profile_dims,a.shape[1],len(a)-1),svd_solver='randomized',random_state=config.seed).fit(a)
            residual=a-self.pca.inverse_transform(self.pca.transform(a))
            self.residual_scale=max(np.sqrt(np.mean(residual**2)),1e-8)
        self.dimensions=a.shape[1] if self.pca is None else self.pca.n_components_+1
        return self
    def transform(self,z):
        a=np.asarray(z)*self.weights
        if self.pca is not None:
            v=self.pca.transform(a);err=np.linalg.norm(a-self.pca.inverse_transform(v),axis=1)/self.residual_scale
            a=np.c_[v,err]
        return a/np.sqrt(self.dimensions)


class ReferenceBank:
    def fit(self,z,ids,groups,quality,config):
        self.z=np.asarray(z,float);self.ids=np.asarray(ids,str);self.groups=np.asarray(groups,str);self.quality=np.asarray(quality,float);self.config=config
        if len(self.z)<config.k or len(np.unique(self.groups))<2:raise ValueError("Insufficient technically valid reference donors")
        self.tree=KDTree(self.z,leaf_size=40);self.radius=None;return self
    def calibrate(self,z,ids,groups):
        kth=[];k=self.config.k
        for begin in range(0,len(z),self.config.query_batch):
            stop=min(begin+self.config.query_batch,len(z));count=min(len(self.ids),k+64)
            dist,ix=self.tree.query(z[begin:stop],k=count)
            for j,(ds,nn) in enumerate(zip(dist,ix)):
                good=(self.ids[nn]!=ids[begin+j])&(self.groups[nn]!=groups[begin+j])
                if good.sum()<k and count<len(self.ids):
                    ds,nn=self.tree.query(z[begin+j:begin+j+1],k=len(self.ids));ds,nn=ds[0],nn[0]
                    good=(self.ids[nn]!=ids[begin+j])&(self.groups[nn]!=groups[begin+j])
                if good.sum()>=k:kth.append(ds[good][k-1])
        if len(kth)<10:raise ValueError("Too few independent tune queries for a reference radius")
        self.radius=float(np.quantile(kth,self.config.radius_quantile))
        if self.radius<=1e-12:raise ValueError("Collapsed reference geometry")
        return self.radius
    def query(self,z,ids,groups):
        if self.radius is None:raise ValueError("Reference radius is not frozen")
        names=('eligible_match_count','matched_count','reference_groups','effective_reference_groups','nearest_radius_ratio','radius_ratio','donor_oof_loss')
        out={n:np.zeros(len(z),float) for n in names}
        out['nearest_radius_ratio'].fill(100);out['radius_ratio'].fill(100);out['donor_oof_loss'].fill(100)
        for start in range(0,len(z),self.config.query_batch):
            nn,dd=self.tree.query_radius(z[start:start+self.config.query_batch],self.radius,return_distance=True)
            for row,(idx,distance) in enumerate(zip(nn,dd)):
                i=start+row;ok=(self.ids[idx]!=ids[i])&(self.groups[idx]!=groups[i]);idx=idx[ok];distance=distance[ok]
                out['eligible_match_count'][i]=len(idx)
                if not len(idx):continue
                rank=np.lexsort((self.ids[idx],distance))[:self.config.k];idx=idx[rank];distance=distance[rank]
                w=np.exp(-.5*(distance/self.radius)**2);w/=w.sum()
                ug,inv=np.unique(self.groups[idx],return_inverse=True);gw=np.bincount(inv,weights=w)
                vals=(len(idx),len(ug),1/np.sum(gw**2),distance[0]/self.radius,distance[-1]/self.radius,w@self.quality[idx])
                for name,value in zip(names[1:],vals):out[name][i]=value
        out['matching_supported']=(out['matched_count']>=self.config.min_twins)&(out['effective_reference_groups']>=self.config.min_group_ess)
        return out


class FrozenRule:
    def fit(self,scores,allowed,ids,coverage,seed):
        self.seed=seed;self.coverage=coverage
        ix=np.flatnonzero(allowed&np.isfinite(scores));want=math.ceil(coverage*len(ids));count=min(want,len(ix))
        self.support_limited=count<want;self.threshold=np.inf;self.tie=0
        if count:
            ties=np.array([hash_number(s,seed) for s in ids],np.uint64);order=ix[np.lexsort((ties[ix],-scores[ix]))];i=order[count-1]
            self.threshold=float(scores[i]);self.tie=int(ties[i])
            if count==len(ix):self.threshold=-np.inf
        return self
    def apply(self,scores,allowed,ids):
        ties=np.array([hash_number(s,self.seed) for s in ids],np.uint64)
        return allowed&np.isfinite(scores)&((scores>self.threshold)|((scores==self.threshold)&(ties<=self.tie)))


class TwinModel:
    def __init__(self,config=Config()):config.validate();self.config=config;self.implementation_version=VERSION
    def fit(self,X,y,ids,groups,roles,features,genetic_indices,covariate_indices,csx_indices=(),disco_indices=()):
        cfg=self.config;X=np.asarray(X,float);y=np.asarray(y,float);ids,groups=validate_ids(ids,groups);roles=np.asarray(roles,str)
        if len(X)!=len(ids) or len(y)!=len(ids) or len(roles)!=len(ids):raise ValueError("Development arrays are not aligned")
        assert_roles(groups,roles);assert_labels(y,cfg)
        self.features=list(map(str,features))
        if len(set(self.features))!=len(self.features) or len(self.features)!=X.shape[1]:raise ValueError("Invalid named feature schema")
        self.ids=ids.copy();self.groups=groups.copy();self.roles=roles.copy()
        part={r:np.flatnonzero(roles==r) for r in ROLE_NAMES};b,t,g,cal,aud=(part[r] for r in ROLE_NAMES)
        for r in ROLE_NAMES:
            if len(part[r])<cfg.min_role_n:raise ValueError(f"Too few participants in frozen {r} role")
        self.prep=NumericProfile().fit(X[b],genetic_indices,covariate_indices,cfg);qc=self.prep.qc(X);z=self.prep.transform(X)
        for r in ROLE_NAMES[:-1]:assert_labels(y[part[r]][qc[part[r]]],cfg)
        quality=oof_quality(X[b],y[b],ids[b],groups[b],genetic_indices,covariate_indices,cfg)
        selection=prune(y[b],quality,ids[b],qc[b],cfg)
        weights={'full':qc[b].astype(float),'clean':selection['keep'].astype(float),'soft':selection['soft'],'random':selection['random_keep'].astype(float)}
        self.predictors={};self.columns={};self.tuning=[]
        parameters=cfg.logistic_c_grid if cfg.trait_type=='binary' else cfg.ridge_grid
        cpos=np.flatnonzero(np.isin(self.prep.keep,covariate_indices))
        specs={f'ABM_{arm}':(np.arange(z.shape[1]),w) for arm,w in weights.items()}
        if len(cpos):specs['Clinical']=(cpos,weights['full'])
        for name,indices in [('CSx_reference',csx_indices),('Disco_reference',disco_indices)]:
            take=np.flatnonzero(np.isin(self.prep.keep,list(covariate_indices)+list(indices)))
            if len(set(indices)&set(self.prep.keep)) and len(take):specs[name]=(take,weights['full'])
        for name,(columns,w) in specs.items():
            best=None
            for value in parameters:
                m=Predictor().fit(z[b][:,columns],y[b],w,cfg,value);p=m.predict(z[t][:,columns]);score=float(np.mean(loss(y[t][qc[t]],p[qc[t]],cfg.trait_type=='binary')))
                self.tuning.append(dict(model=name,penalty=value,tune_loss=score,training_N=int(np.sum(w>0)),tune_N=int(qc[t].sum())))
                if best is None or score<best[0]:best=(score,m)
            self.predictors[name]=best[1];self.columns[name]=columns
        self.fallback='CSx_reference' if 'CSx_reference' in self.predictors else 'ABM_full'
        self.calibration={}
        for name,m in self.predictors.items():
            p=m.predict(z[cal][:,self.columns[name]])
            self.calibration[name]=Calibration().fit(p[qc[cal]],y[cal][qc[cal]],cfg.trait_type=='binary')
        self.geometry=Geometry().fit(z[b][qc[b]][:,self.prep.profile_positions],self.predictors['ABM_full'].coefficients[self.prep.profile_positions],cfg)
        zz=self.geometry.transform(z[:,self.prep.profile_positions]);self.banks={}
        for arm,choose in [('full',selection['eligible']),('clean',selection['keep'])]:
            self.banks[arm]=ReferenceBank().fit(zz[b][choose],ids[b][choose],groups[b][choose],quality['loss'][choose],cfg)
        radius=self.banks['full'].calibrate(zz[t][qc[t]],ids[t][qc[t]],groups[t][qc[t]])
        self.banks['clean'].radius=radius # same caliper: filtering does not silently expand the matching domain
        pg,dg=self._raw(X[g],ids[g],groups[g]);gg=self._gate_inputs(pg,dg)
        target=(y[g]-pg[cfg.primary])**2;good=dg['technical_qc']
        self.gate=HistGradientBoostingRegressor(max_iter=cfg.gate_iterations,max_leaf_nodes=5,min_samples_leaf=max(10,min(100,int(good.sum()/10))),l2_regularization=10,early_stopping=False,random_state=cfg.seed+34).fit(gg[good],target[good])
        self.prediction_center=float(np.mean(y[b][qc[b]]));self.training_variance=float(np.var(y[b][qc[b]]))
        pc,dc=self._raw(X[cal],ids[cal],groups[cal]);scores=self._scores(pc,dc,ids[cal]);self.rules={}
        for name,score in scores.items():
            allowed=dc['technical_qc'] if name in {'baseline_confidence','lowrisk','random'} else dc['technical_qc']&dc['matching_supported']
            for q in cfg.coverages:self.rules[(name,q)]=FrozenRule().fit(score,allowed,ids[cal],q,cfg.seed+57)
        self.risk_edges=np.unique(np.quantile(pc[self.fallback][dc['technical_qc']],np.linspace(0,1,11)))[1:-1]
        self.donors=dict(eid=ids[b],group=groups[b],technical_qc=qc[b]&quality['technical_qc'],donor_eligible=selection['eligible'],retained=selection['keep'],random_retained=selection['random_keep'],soft_weight=selection['soft'],OOF_prediction=quality['prediction'],OOF_prediction_sd=quality['prediction_sd'],OOF_loss=quality['loss'],technical_excluded=~(qc[b]&quality['technical_qc']),profile_discordance_excluded=selection['eligible']&~selection['keep'])
        self.training_summary=dict(build_N=len(b),build_QC_N=int(qc[b].sum()),reference_eligible_N=int(selection['eligible'].sum()),retained_N=int(selection['keep'].sum()),retained_cases=int(y[b][selection['keep']].sum()) if cfg.trait_type=='binary' else None,reference_radius=radius,grouping_status=cfg.grouping_status,profile_kind=cfg.profile_kind,variant_dosage_QC_checked=cfg.profile_kind=='variant_dosage',raw_SNP_combinations_exhaustively_checked=False,profile_dimensions=len(self.prep.profile_keep),matching_dimensions=self.geometry.dimensions,pruning_is_not_proof_of_noise=True)
        self.role_summary=[dict(role=r,N=len(ix),QC_N=int(qc[ix].sum()),roster_hash=id_hash(ids[ix]),labels_used_for_fitting=r!='audit') for r,ix in part.items()]
        self.oof_audit=quality['audit'];self.frozen=True
        return self
    def _raw(self,X,ids,groups):
        z=self.prep.transform(X);qc=self.prep.qc(X)
        p={name:self.calibration[name].predict(m.predict(z[:,self.columns[name]])) for name,m in self.predictors.items()}
        zz=self.geometry.transform(z[:,self.prep.profile_positions]);local=self.banks['clean'].query(zz[qc],ids[qc],groups[qc]);d={name:np.zeros(len(X),dtype=a.dtype) for name,a in local.items()}
        for name,a in local.items():d[name][qc]=a
        for name in ('nearest_radius_ratio','radius_ratio'):d[name][~qc]=100.
        d['donor_oof_loss'][~qc]=100.;d['technical_qc']=qc
        d['genetic_missing_fraction'],d['covariate_missing_fraction']=self.prep.missing(X)
        proximity=np.exp(-np.minimum(d['nearest_radius_ratio'],100)**2)
        support=np.clip(d['effective_reference_groups']/self.config.k,0,1)
        concordance=np.exp(-np.clip(d['donor_oof_loss'],0,100));radius=np.exp(-np.minimum(d['radius_ratio'],100)**2)
        d.update(proximity_score=proximity,reference_support_score=support,donor_concordance_score=concordance,profile_radius_score=radius)
        d['pre_prediction_confidence']=np.where(qc&d['matching_supported'],100*(proximity*support*concordance*radius)**.25,0.)
        return p,d
    def _gate_inputs(self,p,d):
        return np.c_[p[self.config.primary],p[self.fallback],d['pre_prediction_confidence'],np.log1p(d['eligible_match_count']),np.log1p(d['effective_reference_groups']),d['nearest_radius_ratio'],d['radius_ratio'],d['donor_oof_loss']]
    def _scores(self,p,d,ids):
        expected=np.maximum(self.gate.predict(self._gate_inputs(p,d)),0)
        base=p[self.fallback]
        return dict(support_quality=d['pre_prediction_confidence']/100,expected_error=-expected,support_only=-d['radius_ratio']+.01*np.log1p(d['effective_reference_groups']),baseline_confidence=np.abs(base-(.5 if self.config.trait_type=='binary' else self.prediction_center)),lowrisk=-base,random=np.array([hash_number(i,self.config.seed+128)/2**64 for i in ids]))
    def project(self,X,ids,groups,features,*,audit=False):
        """Accepts only predictor values and identity. Never y, outcome, time or event."""
        if getattr(self,'implementation_version',None)!=VERSION:raise ValueError('Saved model implementation differs; preserve the original snapshot and explicitly refit')
        X=np.asarray(X,float);ids,groups=validate_ids(ids,groups)
        if list(features)!=self.features:raise ValueError("Feature names/order do not match the frozen model")
        if len(X)!=len(ids):raise ValueError("Unaligned query")
        if audit:
            allowed=self.roles=='audit'
            if not set(ids)<=set(self.ids[allowed]) or not set(groups)<=set(self.groups[allowed]):raise ValueError("Audit access restricted to frozen audit IDs")
        elif set(ids)&set(self.ids) or set(groups)&set(self.groups):raise ValueError("Query overlaps development individuals or reference groups")
        p,d=self._raw(X,ids,groups);scores=self._scores(p,d,ids);masks={}
        for (name,q),rule in self.rules.items():
            ok=d['technical_qc'] if name in {'baseline_confidence','lowrisk','random'} else d['technical_qc']&d['matching_supported']
            masks[f'{name}:{q:g}']=rule.apply(scores[name],ok,ids)
        selected=masks[f'{self.config.selector}:{self.config.primary_coverage:g}']
        policy=np.where(selected,p[self.config.primary],p[self.fallback]);policy[~d['technical_qc']]=np.nan
        for name in p:p[name]=np.where(d['technical_qc'],p[name],np.nan)
        p['ABM_policy']=policy
        reason=np.where(~d['technical_qc'],'genetic_or_covariate_QC_failed',np.where(~d['matching_supported'],'no_adequate_genetic_twins',np.where(selected,'selected','below_frozen_confidence_cutoff')))
        return dict(eid=ids,group=groups,predictions=p,diagnostics=d,selected=selected,reason=reason,masks=masks,expected_squared_error=-scores['expected_error'],frozen_hash=state_hash(self),evaluation_scope='development_audit' if audit else 'exploratory_heldout',confidence_semantics='reference support, not probability of an individually correct outcome',version=VERSION)
