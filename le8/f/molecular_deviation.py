"""Genetically conditioned molecular deviations and locked CAD validation.

Targets: contemporaneous objective plaque (primary), coronary calcium (separate),
known-CAD state (exploratory). A future incident diagnosis is never relabelled as
prevalent occult disease. No clinical diagnostic claims are hard-coded.
"""
from __future__ import annotations
from dataclasses import dataclass,asdict
from pathlib import Path
import json, math, warnings
import numpy as np
import pandas as pd
from scipy.special import logit
from sklearn.linear_model import Ridge,LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score,average_precision_score,roc_curve
from sklearn.exceptions import ConvergenceWarning
from abm_twin_core import NumericPrep,Config,Calibrator,logloss_rows,model_hash,json_write
from research_io import write_tables

VERSION='MOLECULAR-DEVIATION-2026-10-14'


@dataclass(frozen=True)
class DeviationConfig:
    seed:int=12345
    folds:int=5
    reference_min_n:int=100
    min_class_n:int=20
    sample_missing:float=.2
    ridge_alpha:float=1.
    c_grid:tuple=(.01,.1,1.)
    max_iter:int=3000
    bootstrap:int=500
    target_specificity:float=.95
    clinical_release:bool=False
    task:str='objective_plaque'
    trait:str='cvd_cad'
    layer:str='met'
    max_objective_lag_days:int=90
    gwas_overlap:str='unknown'
    omic_transform:str='identity'
    endpoint_sampling:str='ascertainment_sample'
    grouping_unit:str='family'
    def validate(self):
        if self.grouping_unit not in {'family','individual_unverified'}:raise ValueError('Invalid grouping unit')
        if self.task not in ('objective_plaque','coronary_calcium','known_cad_state'):
            raise ValueError('Unsupported target; incident events cannot stand in for present disease')
        if self.clinical_release:raise ValueError('Research software cannot authorize clinical release')
        if self.folds<2 or self.reference_min_n<10 or self.min_class_n<2 or self.bootstrap<1:raise ValueError('Invalid validation parameters')
        if not 0<self.target_specificity<1 or not 0<=self.sample_missing<1:raise ValueError('Invalid operating point/QC')
        if self.omic_transform not in ('identity','log1p'):raise ValueError('Choose identity or a prespecified log1p transform')


def objective_endpoint(d,config):
    """Returns endpoint eligibility, labels and exclusion reasons. Does not split data.
    objective_cad is a prespecified, independently adjudicated binary target.
    CT coronary calcium is NOT treated as an angiographic stenosis measure.
    """
    config.validate();n=len(d);reason=np.full(n,'eligible',object)
    if config.task=='known_cad_state':
        if 'known_cad_at_index' not in d:raise ValueError('Need known_cad_at_index')
        y=pd.to_numeric(d.known_cad_at_index,errors='coerce').to_numpy(float)
        ok=np.isin(y,[0,1]);reason[~ok]='unknown_baseline_status'
        return ok,y,reason
    required={'known_cad_at_index','index_date','objective_date','objective_cad','objective_modality'}
    if not required<=set(d):
        return np.zeros(n,bool),np.full(n,np.nan),np.full(n,'objective_endpoint_unavailable',object)
    known=pd.to_numeric(d.known_cad_at_index,errors='coerce').to_numpy(float)
    y=pd.to_numeric(d.objective_cad,errors='coerce').to_numpy(float)
    index=pd.to_datetime(d.index_date,errors='coerce');date=pd.to_datetime(d.objective_date,errors='coerce')
    delta=(date-index).dt.total_seconds().to_numpy()/86400
    allowed={'CCTA','INVASIVE_CORONARY_ANGIOGRAPHY'} if config.task=='objective_plaque' else {'CAC_CT'}
    modality=d.objective_modality.fillna('').astype(str).str.upper().to_numpy()
    ok=(known==0)&np.isin(y,[0,1])&np.isfinite(delta)&(np.abs(delta)<=config.max_objective_lag_days)&np.isin(modality,list(allowed))
    reason[known!=0]='known_or_unknown_CAD_history'
    reason[~np.isin(y,[0,1])]='objective_result_missing'
    reason[~np.isfinite(delta)]='objective_date_missing'
    reason[np.isfinite(delta)&(np.abs(delta)>config.max_objective_lag_days)]='not_contemporaneous'
    reason[~np.isin(modality,list(allowed))]='wrong_objective_modality'
    future=np.isfinite(delta)&(delta!=0)  # exclude intervening changes on either side of index
    clear=np.ones(n,bool)
    for col in ('intervening_cad','intervening_treatment'):
        val=pd.to_numeric(d[col],errors='coerce').to_numpy() if col in d else np.full(n,np.nan)
        clear &= val==0
    ok &= ~future|clear
    reason[future&~clear]='intervening_event_or_treatment_not_excluded'
    reason[ok]='eligible'
    return ok,y,reason


def bh(p):
    p=np.asarray(p,float);q=np.full(len(p),np.nan);ix=np.flatnonzero(np.isfinite(p))
    order=ix[np.argsort(p[ix])]
    if len(order):q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(order)+1))[::-1])[::-1])
    return q


class CovPrep:
    def fit(self,C):
        C=np.asarray(C,float)
        if C.ndim!=2 or not C.shape[1]:raise ValueError('Supply explicit numeric clinical/ancestry covariates')
        if np.isinf(C).any():raise ValueError('Infinite covariate')
        self.med=np.nanmedian(C,axis=0)
        if not np.isfinite(self.med).all():raise ValueError('Entirely missing covariate in reference training')
        x=np.where(np.isfinite(C),C,self.med);self.mean=x.mean(0);self.sd=np.maximum(x.std(0),1e-8)
        return self
    def transform(self,C):
        C=np.asarray(C,float)
        if C.ndim!=2 or C.shape[1]!=len(self.med):raise ValueError('Covariate schema mismatch')
        return (np.where(np.isfinite(C),C,self.med)-self.mean)/self.sd


class MolecularNormalizer:
    """Per-assay M~C+own-PGS versus M~C. The remainder is not pure environment."""
    def fit(self,M,G,C,reference,config):
        self.config=config;self.cov=CovPrep().fit(C[reference]);c=self.cov.transform(C)
        self.full=[];self.base=[];self.gmean=[];self.gsd=[];self.mmean=[];self.msd=[];self.nref=[]
        for j in range(M.shape[1]):
            ok=reference&np.isfinite(M[:,j])&np.isfinite(G[:,j])
            if ok.sum()<config.reference_min_n or np.std(M[ok,j])<1e-10 or len(np.unique(G[ok,j]))<2:
                self.full.append(None);self.base.append(None);self.gmean.append(np.nan);self.gsd.append(np.nan);self.mmean.append(np.nan);self.msd.append(np.nan);self.nref.append(int(ok.sum()));continue
            gm,gs=G[ok,j].mean(),G[ok,j].std();mm,ms=M[ok,j].mean(),M[ok,j].std()
            z=(M[ok,j]-mm)/ms;gx=(G[ok,j]-gm)/gs
            self.full.append(Ridge(alpha=config.ridge_alpha).fit(np.c_[c[ok],gx],z))
            self.base.append(Ridge(alpha=config.ridge_alpha).fit(c[ok],z))
            self.gmean.append(gm);self.gsd.append(gs);self.mmean.append(mm);self.msd.append(ms);self.nref.append(int(ok.sum()))
        for nm in ('gmean','gsd','mmean','msd'):setattr(self,nm,np.asarray(getattr(self,nm)))
        self.valid=np.array([m is not None for m in self.full]);return self
    def predict(self,G,C):
        c=self.cov.transform(C);a=np.full_like(G,np.nan,dtype=float);b=a.copy()
        for j in np.flatnonzero(self.valid):
            good=np.isfinite(G[:,j]);g=(G[good,j]-self.gmean[j])/self.gsd[j]
            a[good,j]=self.mmean[j]+self.msd[j]*self.full[j].predict(np.c_[c[good],g])
            b[good,j]=self.mmean[j]+self.msd[j]*self.base[j].predict(c[good])
        return a,b


class DeviationStudy:
    def __init__(self,config=DeviationConfig()):config.validate();self.config=config
    def _transform(self,M):
        M=np.asarray(M,float)
        if self.config.omic_transform=='log1p':
            if np.any(M[np.isfinite(M)]<=-1):raise ValueError('log1p domain violation')
            M=np.log1p(M)
        return M
    def fit_reference(self,M,G,C,ids,groups,reference,feature_names,covariate_names):
        """Only training-role input. Cross-fits every reference normalizer."""
        M=self._transform(M);G=np.asarray(G,float);C=np.asarray(C,float);reference=np.asarray(reference,bool)
        ids=np.asarray(ids,str);groups=np.asarray(groups,str)
        if M.shape!=G.shape or len(M)!=len(C) or len(M)!=len(ids):raise ValueError('Array alignment mismatch')
        if len(set(ids))!=len(ids) or any(v in ('','nan','None') for v in groups):raise ValueError('Invalid participant/family IDs')
        if np.isinf(M).any() or np.isinf(G).any():raise ValueError('Infinite molecular data')
        self.features=list(map(str,feature_names));self.covariates=list(map(str,covariate_names))
        if len(self.features)!=M.shape[1] or len(set(self.features))!=len(self.features) or len(self.covariates)!=C.shape[1]:raise ValueError('Named schema mismatch')
        self.normalizer=MolecularNormalizer().fit(M,G,C,reference,self.config)
        oof=np.full_like(M,np.nan);oof0=oof.copy();audit=[]
        split=GroupKFold(n_splits=self.config.folds)
        for f,(tr,va) in enumerate(split.split(M,groups=groups)):
            assert not set(groups[tr])&set(groups[va])
            norm=MolecularNormalizer().fit(M[tr],G[tr],C[tr],reference[tr],self.config)
            oof[va],oof0[va]=norm.predict(G[va],C[va]);audit.append(dict(fold=f,train_N=len(tr),validation_N=len(va),family_disjoint=True))
        self.reference_audit=audit;self.reference_table=[];self.resid_sd=[];self.base_sd=[]
        self.reference_train_ids=ids.copy();self.reference_train_groups=groups.copy()
        for j,name in enumerate(self.features):
            ok=reference&np.isfinite(M[:,j])&np.isfinite(oof[:,j])&np.isfinite(oof0[:,j])
            err=M[ok,j]-oof[ok,j];err0=M[ok,j]-oof0[ok,j]
            mse0=np.sum(err0**2);sst=np.sum((M[ok,j]-np.mean(M[ok,j]))**2) if ok.any() else 0
            self.resid_sd.append(float(np.sqrt(np.mean(err**2))) if len(err)>1 else np.nan)
            self.base_sd.append(float(np.sqrt(np.mean(err0**2))) if len(err0)>1 else np.nan)
            r2=1-np.sum(err**2)/mse0 if mse0>0 else np.nan
            self.reference_table.append(dict(feature=name,reference_N=self.normalizer.nref[j],oof_N=int(ok.sum()),
                own_PGS_partial_R2=r2,prediction_R2=1-np.sum(err**2)/sst if sst>0 else np.nan,
                calibration_slope=float(self.normalizer.full[j].coef_[-1]) if self.normalizer.valid[j] else np.nan,
                status='ok' if self.normalizer.valid[j] and np.isfinite(r2) else 'unavailable',weak_capture=not np.isfinite(r2) or r2<.005,
                gwas_overlap=self.config.gwas_overlap,interpretation='association calibration; not environmental or causal decomposition'))
        self.resid_sd=np.asarray(self.resid_sd);self.base_sd=np.asarray(self.base_sd)
        self.keep=np.flatnonzero(self.normalizer.valid&np.isfinite(self.resid_sd)&(self.resid_sd>1e-8)&np.isfinite(self.base_sd)&(self.base_sd>1e-8))
        if len(self.keep)<2:raise ValueError('Fewer than two calibrated molecular scores')
        self.train_oof_Z=(M-oof)/self.resid_sd;self.train_oof_Z0=(M-oof0)/self.base_sd
        self.classifiers={};self.classification_status='objective_endpoint_unavailable'
        return self
    def deviations(self,M,G,C):
        M=self._transform(M);expected,expected0=self.normalizer.predict(np.asarray(G,float),np.asarray(C,float))
        return dict(expected=expected,Z=(M-expected)/self.resid_sd,Z0=(M-expected0)/self.base_sd,M=M)
    def _blocks(self,M,G,C,Z=None,Z0=None):
        r=self.deviations(M,G,C)
        if Z is not None:r['Z']=Z
        if Z0 is not None:r['Z0']=Z0
        k=self.keep;c=np.asarray(C,float);m=r['M'][:,k];g=np.asarray(G,float)[:,k];z=r['Z'][:,k];z0=r['Z0'][:,k]
        return {'Clinical':c,'Clinical_G':np.c_[c,g],'Clinical_M':np.c_[c,m],
                'Clinical_M_G':np.c_[c,m,g],'Clinical_Z0':np.c_[c,z0],
                'Clinical_Z':np.c_[c,z],'Clinical_Z_G':np.c_[c,z,g]},r
    def fit_disease_models(self,development):
        """Development roles are train/tune/calibrate. No test labels accepted."""
        d=development;roles=np.asarray(d['role'],str);ids=np.asarray(d['eid'],str);groups=np.asarray(d['group'],str)
        if set(roles)!= {'train','tune','calibrate'}:raise ValueError('Need three development roles, no test')
        order=np.argsort(groups,kind='stable')
        if np.any((groups[order][1:]==groups[order][:-1]) & (roles[order][1:]!=roles[order][:-1])):
            raise ValueError('Family crosses roles')
        parts={r:np.flatnonzero(roles==r) for r in set(roles)};tr,tu,ca=[parts[r] for r in ('train','tune','calibrate')]
        if not np.array_equal(ids[tr],self.reference_train_ids):raise ValueError('Reference fit training roster differs')
        y=np.asarray(d['y'],float);valid=np.asarray(d['endpoint_eligible'],bool)&np.isin(y,[0,1])
        self.development_ids=ids;self.development_groups=groups
        blocks,r=self._blocks(d['M'],d['G'],d['C'])
        # Avoid normalizer overfit only for classifier training. Other roles are projected.
        kz=self.keep;blocks['Clinical_Z'][tr]=np.c_[d['C'][tr],self.train_oof_Z[:,kz]]
        blocks['Clinical_Z0'][tr]=np.c_[d['C'][tr],self.train_oof_Z0[:,kz]]
        blocks['Clinical_Z_G'][tr]=np.c_[d['C'][tr],self.train_oof_Z[:,kz],d['G'][tr][:,kz]]
        missing=np.maximum(np.mean(~np.isfinite(d['M'][:,kz]),axis=1),np.mean(~np.isfinite(d['G'][:,kz]),axis=1))
        valid &= (missing<=self.config.sample_missing)&np.all(np.isfinite(np.asarray(d['C'],float)),axis=1)
        self.development_qc=dict(total_N=len(y),endpoint_QC_N=int(valid.sum()))
        if not valid.any():return self
        for role,ix in parts.items():
            counts=[int(np.sum(valid[ix]&(y[ix]==v))) for v in (0,1)]
            if min(counts)<self.config.min_class_n:
                self.classification_status='insufficient_'+role+'_cases_or_controls';return self
        self.model_tuning=[];self.operating_thresholds={}
        for name,x in blocks.items():
            prep=CovPrep().fit(x[tr][valid[tr]]);a=prep.transform(x);best=None
            for c in self.config.c_grid:
                m=LogisticRegression(solver='saga',penalty='elasticnet',l1_ratio=.5,C=c,max_iter=self.config.max_iter,random_state=self.config.seed)
                with warnings.catch_warnings():
                    warnings.simplefilter('error',ConvergenceWarning)
                    m.fit(a[tr][valid[tr]],y[tr][valid[tr]].astype(int))
                p=m.predict_proba(a[tu][valid[tu]])[:,1];score=float(np.mean(logloss_rows(y[tu][valid[tu]],p)))
                self.model_tuning.append(dict(model=name,C=c,tune_logloss=score))
                if best is None or score<best[0]:best=score,m
            m=best[1];pc=m.predict_proba(a[ca][valid[ca]])[:,1];yc=y[ca][valid[ca]].astype(int)
            cal=Calibrator().fit(pc,yc,np.ones(len(yc)))
            pp=cal.predict(pc);threshold=float(np.quantile(pp[yc==0],self.config.target_specificity,method='higher'))
            self.classifiers[name]=(prep,m,cal);self.operating_thresholds[name]=threshold
        self.classification_status='fitted_locked_research'
        return self
    def project(self,M,G,C,ids,groups,features,covariates):
        """No disease labels/dates in this interface. Projection cannot retune."""
        ids=np.asarray(ids,str);groups=np.asarray(groups,str)
        if list(map(str,features))!=self.features or list(map(str,covariates))!=self.covariates:raise ValueError('Feature names/order mismatch')
        if len(set(ids))!=len(ids) or len(ids)!=len(M):raise ValueError('Projection ID mismatch')
        dev_ids=set(getattr(self,'development_ids',self.reference_train_ids));dev_groups=set(getattr(self,'development_groups',self.reference_train_groups))
        if set(ids)&dev_ids or set(groups)&dev_groups:raise ValueError('Projection overlaps development IDs/families')
        blocks,r=self._blocks(M,G,C);p={}
        for name,(prep,m,cal) in self.classifiers.items():p[name]=cal.predict(m.predict_proba(prep.transform(blocks[name]))[:,1])
        k=self.keep;missing=np.maximum(np.mean(~np.isfinite(M[:,k]),axis=1),np.mean(~np.isfinite(G[:,k]),axis=1))
        qc=(missing<=self.config.sample_missing)&np.all(np.isfinite(np.asarray(C,float)),axis=1)
        return dict(eid=ids,group=groups,probability=p,Z=r['Z'],Z0=r['Z0'],measured=r['M'],expected=r['expected'],covariate_values=np.asarray(C,float),qc=qc,
                    frozen_hash=model_hash(self),features=self.features,classification_status=self.classification_status)


def binary_metrics(y,p,threshold=.5):
    y=np.asarray(y,int);p=np.asarray(p,float)
    if not len(y):return dict(N=0,cases=0,controls=0,status='empty',AUC=np.nan,AUPRC=np.nan,Brier=np.nan,LogLoss=np.nan)
    both=len(np.unique(y))==2;positive=p>=threshold;tp=np.sum(positive&(y==1));fp=np.sum(positive&(y==0))
    return dict(N=len(y),cases=int(y.sum()),controls=int(np.sum(y==0)),status='ok' if both else 'one_class',
                AUC=float(roc_auc_score(y,p)) if both else np.nan,AUPRC=float(average_precision_score(y,p)) if both else np.nan,
                Brier=float(np.mean((y-p)**2)),LogLoss=float(np.mean(logloss_rows(y,p))),
                threshold=threshold,sensitivity=float(tp/max(1,np.sum(y==1))),specificity=float(1-fp/max(1,np.sum(y==0))),
                PPV=float(tp/max(1,tp+fp)))


def feature_state_associations(pred,y,eligible):
    """Evaluation-only adjusted mean deviation difference (NOT an odds ratio).
    Family-cluster sandwich standard errors; FDR across the assayed family.
    """
    from scipy.stats import t as student_t
    C=pred['covariate_values'];rows=[]
    for j,feature in enumerate(pred['features']):
        good=eligible&np.isfinite(pred['Z'][:,j])&np.all(np.isfinite(C),axis=1)
        row=dict(feature=feature,N=int(good.sum()),beta=np.nan,se=np.nan,p=np.nan,lo=np.nan,hi=np.nan,
                 effect_measure='adjusted mean Z difference, objective case minus control',status='insufficient_information')
        if good.sum()<30 or len(np.unique(y[good]))<2:rows.append(row);continue
        cov=C[good];cov=cov[:,np.std(cov,axis=0)>1e-10];x=np.c_[np.ones(good.sum()),y[good],cov]
        if np.linalg.matrix_rank(x)<x.shape[1]:row['status']='rank_deficient_covariates';rows.append(row);continue
        inverse=np.linalg.inv(x.T@x);b=inverse@x.T@pred['Z'][good,j];r=pred['Z'][good,j]-x@b
        groups,gi=np.unique(pred['group'][good],return_inverse=True);ng=len(groups)
        if ng<5:rows.append(row);continue
        score=r*(x@inverse[:,1]);score=np.bincount(gi,weights=score)
        var=float(score@score)*ng/(ng-1)*(len(r)-1)/(len(r)-x.shape[1]);se=np.sqrt(max(0,var))
        crit=student_t.ppf(.975,ng-1)
        row.update(beta=float(b[1]),se=float(se),p=float(2*student_t.sf(abs(b[1]/se),ng-1)) if se>0 else np.nan,
                   lo=float(b[1]-crit*se),hi=float(b[1]+crit*se),family_N=ng,status='estimated_cluster_robust')
        rows.append(row)
    d=pd.DataFrame(rows)
    if len(d):d['FDR']=bh(d.p)
    return d


def evaluate_deviation(model,pred,truth,output):
    """All cohort eligibility follows a prespecified endpoint; no best-test-model selection."""
    if pred['frozen_hash']!=model_hash(model):raise ValueError('Frozen model changed')
    if not np.array_equal(pred['eid'],np.asarray(truth['eid'],str)):raise ValueError('Truth ID/order mismatch')
    y=np.asarray(truth['y'],float);eligible=np.asarray(truth['endpoint_eligible'],bool)&np.isin(y,[0,1])&pred['qc']
    cfg=model.config;rows=[];roc=[];cal=[];dca=[];contrasts=[];dist=[];paired=[]
    for name,p in pred['probability'].items():
        rows.append(dict(model=name,task=cfg.task,**binary_metrics(y[eligible],p[eligible],model.operating_thresholds[name])))
        if len(np.unique(y[eligible]))==2:
            f,t,_=roc_curve(y[eligible],p[eligible]);roc.extend(dict(model=name,FPR=float(a),TPR=float(b)) for a,b in zip(f,t))
        for b in range(10):
            take=eligible&(p>=b/10)&(p<((b+1)/10) if b<9 else p<=1)
            if take.sum()>=10:cal.append(dict(model=name,bin=b,N=int(take.sum()),predicted=float(p[take].mean()),observed=float(y[take].mean())))
        for pt in np.linspace(.01,.5,50):
            pos=p[eligible]>=pt;yy=y[eligible];n=len(yy)
            if n:
                net=(np.sum(pos&(yy==1))-np.sum(pos&(yy==0))*pt/(1-pt))/n
                dca.append(dict(model=name,threshold=float(pt),net_benefit=float(net),treat_all=float(yy.mean()-(1-yy.mean())*pt/(1-pt)),treat_none=0.))
    pairs=[('incremental_PGS','Clinical_M_G','Clinical_M'),('deviation_vs_raw','Clinical_Z','Clinical_M'),
           ('genetic_vs_nongenetic_deviation','Clinical_Z','Clinical_Z0'),('residual_parameterization','Clinical_Z_G','Clinical_M_G')]
    names=pred['probability'];rng=np.random.default_rng(cfg.seed+816)
    u,inverse=np.unique(pred['group'],return_inverse=True)
    order=np.argsort(inverse,kind='stable')
    positions=np.split(order,np.cumsum(np.bincount(inverse))[:-1])
    samples={label:{m:[] for m in ('AUC','AUPRC','Brier','LogLoss')} for label,a,b in pairs if a in names and b in names}
    for _ in range(cfg.bootstrap if samples and eligible.any() else 0):
        ix=np.concatenate([positions[j] for j in rng.integers(0,len(u),len(u))]);ix=ix[eligible[ix]]
        if len(np.unique(y[ix]))<2:continue
        met={k:binary_metrics(y[ix],p[ix]) for k,p in names.items()}
        for label,a,b in pairs:
            if label in samples:
                for metric in samples[label]:samples[label][metric].append(met[a][metric]-met[b][metric])
    lookup={r['model']:r for r in rows}
    for label,a,b in pairs:
        if label not in samples:continue
        for metric,v in samples[label].items():
            valid=len(v)>=max(20,int(.8*cfg.bootstrap));bounds=np.quantile(v,[.025,.975]) if valid else (np.nan,np.nan)
            contrasts.append(dict(contrast=label,model=a,reference=b,metric=metric,delta=lookup[a][metric]-lookup[b][metric],
                lower=bounds[0],upper=bounds[1],valid_boot=len(v),status='estimated' if valid else 'insufficient_bootstrap',
                uncertainty=cfg.grouping_unit+' resampling; conditional on frozen development fitting',clinical_release=False))
    # Feature-wise held-out distributions are generated only when an objective/known-state endpoint exists.
    for j,f in enumerate(model.features):
        if j not in model.keep:continue
        for state in (0,1):
            take=eligible&(y==state)&np.isfinite(pred['Z'][:,j]);z=pred['Z'][take,j]
            if len(z)<10:continue
            counts,edges=np.histogram(z,bins=np.r_[-np.inf,np.arange(-5,5.25,.25),np.inf])
            for b,c in enumerate(counts):
                # Small cells are suppressed, not rounded up or changed to zero.
                dist.append(dict(feature=f,state=state,bin=b,lower=edges[b],upper=edges[b+1],N=int(c) if c>=10 else np.nan,
                                 status='visible' if c>=10 else 'suppressed_lt10'))
        take=eligible&np.isfinite(pred['Z'][:,j])&np.isfinite(pred['measured'][:,j])
        for state in (0,1):
            v=take&(y==state)
            if v.sum()>=10:paired.append(dict(feature=f,state=state,N=int(v.sum()),mean_measured=float(pred['measured'][v,j].mean()),
                 mean_expected=float(pred['expected'][v,j].mean()),mean_deviation=float(pred['Z'][v,j].mean()),sd_deviation=float(pred['Z'][v,j].std())))
    status=dict(module='molecular_deviation',task=cfg.task,classification_status=model.classification_status,
                validation_status='evaluated_internal_holdout' if eligible.any() and len(names) else 'not_validated_for_diagnosis',
                target_N=int(eligible.sum()),original_N=len(y),clinical_release=False,gwas_overlap=cfg.gwas_overlap,
                residuals_do_not_create_information=True,external_validation=False,endpoint_sampling=cfg.endpoint_sampling,
                PPV_and_net_benefit_scope='observed ascertainment sample; not automatically population screening')
    tables={'reference_calibration':pd.DataFrame(model.reference_table),'model_metrics':pd.DataFrame(rows),
            'paired_contrasts':pd.DataFrame(contrasts),'roc_points':pd.DataFrame(roc),'probability_calibration':pd.DataFrame(cal),
            'decision_curve':pd.DataFrame(dca),'feature_distribution':pd.DataFrame(dist),'feature_state_summary':pd.DataFrame(paired),
            'feature_state_association':feature_state_associations(pred,y,eligible),'status':pd.DataFrame([status])}
    for table in tables.values():table['grouping_unit']=cfg.grouping_unit
    write_tables(output,tables,'Molecular deviation',cfg.layer,cfg.trait,run_id=pred['frozen_hash'],workbook='molecular_deviation.results.xlsx')
    json_write(Path(output)/'status.json',status)
    return tables,status
