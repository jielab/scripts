"""Training-only QC, clinical coding and molecular residualization."""
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler, SplineTransformer
from sklearn.linear_model import Ridge
from common import words

def parsed_date(raw, name):
    if pd.api.types.is_numeric_dtype(raw) and raw.notna().any():
        raise ValueError(f"{name}: numeric dates are ambiguous; use ISO dates or R Date")
    converted = pd.to_datetime(raw, errors="coerce")
    bad = raw.notna() & raw.astype(str).str.strip().ne("") & converted.isna()
    if bad.any():
        raise ValueError(f"Invalid dates in {name}")
    return converted

def outcomes(p, a):
    p = p.copy()
    if a.outcome_type == "quantitative":
        raw = p[a.target_col]
        value = pd.to_numeric(raw, errors="coerce")
        if (raw.notna() & value.isna()).any():
            raise ValueError(f"Non-numeric quantitative outcome: {a.target_col}")
        p["target"] = value
        p["eligible"] = np.isfinite(value)
        return p, dict(joined=len(p), eligible=int(p.eligible.sum()),
                       missing_target=int((~p.eligible).sum()), outcome_type="quantitative")
    baseline = parsed_date(p[a.baseline_col], a.baseline_col)
    diagnosis = parsed_date(p[a.diagnosis_col], a.diagnosis_col)
    censor = pd.concat([parsed_date(p[a.death_col], a.death_col),
                        parsed_date(p[a.lost_col], a.lost_col),
                        pd.Series(pd.Timestamp(a.end_date), index=p.index)], axis=1).min(axis=1)
    prevalent = diagnosis.notna() & (diagnosis <= baseline)
    unknown = pd.Series(False, index=p.index)
    if a.disease_evidence_col:
        evidence = pd.to_numeric(p[a.disease_evidence_col], errors="raise")
        unknown = evidence.gt(0) & diagnosis.isna()
    valid = baseline.notna() & (censor > baseline)
    healthy = pd.Series(True, index=p.index)
    for col in words(a.healthy_date_cols):
        d = parsed_date(p[col], col)
        healthy &= ~(d.notna() & (d <= baseline))
    event = valid & ~prevalent & diagnosis.notna() & (diagnosis > baseline) & (diagnosis <= censor)
    p["time"] = (censor.where(~event, diagnosis) - baseline).dt.days / 365.25
    p["event"] = event.astype(int)
    p["eligible"] = valid & ~prevalent & ~unknown & healthy
    return p, dict(joined=len(p), eligible=int(p.eligible.sum()),
        prevalent=int(prevalent.sum()), same_day=int((diagnosis == baseline).sum()),
        invalid_followup=int((~valid).sum()), unknown_diagnosis=int(unknown.sum()),
        other_baseline_disease=int((~healthy & ~prevalent).sum()),
        diagnosis_after_censor=int((diagnosis > censor).sum()),
        incident=int(p.loc[p.eligible, "event"].sum()), outcome_type="survival")

class MetadataDesign:
    """Normalize categorical types and fit all encodings on training rows only."""
    def __init__(self, cols, categorical, age_spline=True):
        self.cols = list(cols)
        self.categorical = [c for c in self.cols if c in categorical]
        self.age_spline = age_spline

    def normalized(self, p):
        p = p[self.cols].copy()
        for col in self.cols:
            if col in self.categorical:
                p[col] = p[col].map(lambda v: np.nan if pd.isna(v) else
                    str(int(v)) if isinstance(v, (float, np.floating)) and np.isfinite(v) and v.is_integer()
                    else str(v)).astype(object)
            else:
                raw = p[col]
                p[col] = pd.to_numeric(raw, errors="coerce").replace([np.inf, -np.inf], np.nan)
                if (raw.notna() & p[col].isna()).any():
                    raise ValueError(f"Non-numeric covariate {col}; declare it categorical if appropriate")
        return p

    def fit(self, p):
        p = self.normalized(p)
        if p.isna().all().any():
            raise ValueError("All-missing training covariate: " + ",".join(p.columns[p.isna().all()]))
        numeric = [c for c in self.cols if c not in self.categorical]
        spline = ["age"] if self.age_spline and "age" in numeric else []
        numeric = [c for c in numeric if c not in spline]
        pieces = []
        if numeric:
            pieces.append(("num", make_pipeline(SimpleImputer(strategy="median"),
                                               StandardScaler()), numeric))
        if spline:
            pieces.append(("age", make_pipeline(SimpleImputer(strategy="median"),
                SplineTransformer(n_knots=4, degree=3, include_bias=False),
                StandardScaler()), spline))
        if self.categorical:
            pieces.append(("cat", make_pipeline(SimpleImputer(strategy="most_frequent"),
                OneHotEncoder(handle_unknown="ignore", drop="first", sparse_output=False)),
                self.categorical))
        self.transformer = ColumnTransformer(pieces)
        self.transformer.fit(p)
        self.levels = {c: set(p[c].dropna()) for c in self.categorical}
        return self

    def transform(self, p):
        if not self.cols:
            return np.empty((len(p), 0), dtype="float32")
        return np.asarray(self.transformer.transform(self.normalized(p)), dtype="float32")

    def unknown_categories(self, p):
        p = self.normalized(p)
        return {c: int((p[c].notna() & ~p[c].isin(self.levels[c])).sum())
                for c in self.categorical}

class MolecularPreprocessor:
    def __init__(self, max_missing=.2, residual_cols=(), categorical=(), transform="none"):
        self.max_missing = max_missing
        self.residual_cols = list(residual_cols)
        self.categorical = list(categorical)
        self.transform_name = transform

    def scale_transform(self, x):
        x = np.asarray(x, dtype="float32")
        x = np.where(np.isfinite(x), x, np.nan)
        if self.transform_name == "log1p":
            if np.any(x < 0):
                raise ValueError("log1p requires nonnegative measured abundances; use --transform none for pretransformed data")
            x = np.log1p(x)
        return x

    def fit(self, x, p, allowed=None):
        x = self.scale_transform(x)
        keep = (np.mean(~np.isfinite(x), axis=0) < self.max_missing) & (np.nanstd(x, axis=0) > 1e-8)
        if allowed is not None:
            # Freeze the feature set used for sample QC; do not move its denominator.
            keep = np.asarray(allowed, dtype=bool).copy()
        self.keep = np.flatnonzero(keep)
        if len(self.keep) < 3:
            raise ValueError("Fewer than three molecular features pass training QC")
        z = x[:, self.keep]
        if np.any(np.all(~np.isfinite(z), axis=0)):
            raise ValueError("A retained assay has no observed training values after sample QC")
        self.lower, self.upper = np.nanquantile(z, [.005, .995], axis=0)
        z = np.clip(z, self.lower, self.upper)
        self.median = np.nanmedian(z, axis=0)
        z = np.where(np.isfinite(z), z, self.median)
        self.adjust = None
        if self.residual_cols:
            self.design = MetadataDesign(self.residual_cols, self.categorical).fit(p)
            c = self.design.transform(p)
            self.adjust = Ridge(alpha=1, solver="cholesky").fit(c, z)
            z = z - self.adjust.predict(c)
        self.mean, self.scale = z.mean(axis=0), z.std(axis=0)
        self.scale[self.scale < 1e-6] = 1
        return self

    def transform(self, x, p):
        z = self.scale_transform(x)[:, self.keep]
        observed = np.isfinite(z)
        z = np.clip(np.where(observed, z, self.median), self.lower, self.upper)
        if self.adjust is not None:
            z -= self.adjust.predict(self.design.transform(p))
        return np.clip((z - self.mean) / self.scale, -10, 10).astype("float32"), observed
