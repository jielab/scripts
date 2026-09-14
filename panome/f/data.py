"""Cohort and training-only preprocessing; no endpoint enters representation fitting."""
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.linear_model import Ridge
from sklearn.model_selection import train_test_split, GroupShuffleSplit


def csv_list(s):
    return [x.strip() for x in s.split(',') if x.strip()]


def outcomes(p, a):
    p = p.copy()
    def date(col):
        raw = p[col]
        # Numeric epoch days are ambiguous; the R adapter emits ISO dates.
        if pd.api.types.is_numeric_dtype(raw) and raw.notna().any():
            raise ValueError(f'{col}: dates must be ISO strings, not numeric')
        z = pd.to_datetime(raw, errors='coerce')
        if (raw.notna() & z.isna()).any():
            raise ValueError(f'Unparseable date in {col}')
        return z
    baseline = date(a.baseline_col)
    diagnosis = date(a.diagnosis_col)
    end = pd.Timestamp(a.end_date)
    censor = pd.concat([date(a.death_col), date(a.lost_col), pd.Series(end,index=p.index)],axis=1).min(axis=1)
    prevalent = diagnosis.notna() & (diagnosis <= baseline)
    dirty = pd.Series(False,index=p.index)
    if a.disease_evidence_col:
        evidence = pd.to_numeric(p[a.disease_evidence_col],errors='raise')
        dirty = evidence.gt(0) & diagnosis.isna()
    valid = baseline.notna() & (censor > baseline) & ~dirty
    event = valid & ~prevalent & diagnosis.notna() & (diagnosis > baseline) & (diagnosis <= censor)
    stop = censor.where(~event,diagnosis)
    p['time'] = (stop-baseline).dt.days / 365.25
    p['event'] = event.astype(int)
    healthy = ~prevalent
    for col in csv_list(a.healthy_date_cols):
        d = date(col)
        healthy &= ~(d.notna() & (d <= baseline))
    p['eligible'] = valid & ~prevalent & healthy & (p.time > 0)
    audit = dict(joined=len(p), prevalent=int(prevalent.sum()), same_day=int((diagnosis==baseline).sum()),
                 invalid_followup=int((~valid).sum()), unknown_diagnosis_with_evidence=int(dirty.sum()),
                 diagnosis_after_censor=int((diagnosis>censor).sum()),
                 excluded_other_baseline_disease=int((~healthy & ~prevalent).sum()),
                 eligible=int(p.eligible.sum()), incident=int(p.loc[p.eligible,'event'].sum()))
    return p, audit


def splits(p, a):
    idx = np.arange(len(p))
    if a.group_col:
        if p[a.group_col].isna().any():
            raise ValueError('Split-group IDs must be complete')
        groups = p[a.group_col].astype(str).values
        trv, te = next(GroupShuffleSplit(n_splits=1,test_size=.2,random_state=a.seed).split(idx,groups=groups))
        tr0, va0 = next(GroupShuffleSplit(n_splits=1,test_size=.25,random_state=a.seed+1).split(trv,groups=groups[trv]))
        tr, va = trv[tr0], trv[va0]
    else:
        trv, te = train_test_split(idx,test_size=.2,random_state=a.seed,stratify=p.event)
        tr, va = train_test_split(trv,test_size=.25,random_state=a.seed+1,stratify=p.event.iloc[trv])
    part = np.full(len(p),'test',dtype=object)
    part[tr]='train'; part[va]='validation'
    return part


def design(cols, categorical):
    cat = [x for x in cols if x in categorical]
    num = [x for x in cols if x not in categorical]
    return ColumnTransformer([
        ('num',make_pipeline(SimpleImputer(strategy='median',keep_empty_features=True),StandardScaler()),num),
        ('cat',make_pipeline(SimpleImputer(strategy='most_frequent',keep_empty_features=True),
                            OneHotEncoder(handle_unknown='ignore',sparse_output=False,drop='first')),cat)],
        remainder='drop')


class MolecularPreprocessor:
    def __init__(self, max_missing=.2, residual_cols=(), categorical=()):
        self.max_missing=max_missing; self.residual_cols=list(residual_cols); self.categorical=list(categorical)

    def fit(self, x, p):
        self.keep = np.flatnonzero((np.mean(~np.isfinite(x),axis=0)<self.max_missing) & (np.nanstd(x,axis=0)>1e-8))
        if len(self.keep)<3: raise ValueError('Fewer than three features pass training QC')
        z=x[:,self.keep]
        self.median=np.nanmedian(np.where(np.isfinite(z),z,np.nan),axis=0)
        z=np.where(np.isfinite(z),z,self.median)
        self.lower,self.upper=np.quantile(z,[.005,.995],axis=0)
        z=np.clip(z,self.lower,self.upper)
        self.adjust=None
        if self.residual_cols:
            if p[self.residual_cols].isna().all().any(): raise ValueError('All-missing residualization covariate')
            self.design=design(self.residual_cols,self.categorical)
            c=self.design.fit_transform(p)
            self.adjust=Ridge(alpha=1.,solver='cholesky').fit(c,z)
            z=z-self.adjust.predict(c)
        self.mean=z.mean(axis=0); self.scale=z.std(axis=0)
        self.scale[self.scale<1e-6]=1
        return self

    def transform(self,x,p):
        z=x[:,self.keep]
        mask=np.isfinite(z)
        z=np.clip(np.where(mask,z,self.median),self.lower,self.upper)
        if self.adjust is not None: z=z-self.adjust.predict(self.design.transform(p))
        z=np.clip((z-self.mean)/self.scale,-10,10).astype('float32')
        return z,mask
