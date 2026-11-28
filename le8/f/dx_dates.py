"""Build five baseline endpoints from the cohort, enriched with linked AOD evidence."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

import numpy as np
import pandas as pd

from dx_layout import DATES as DEFINITIONS

BASE_ENDPOINTS = ("icd", "ref", "icd1", "ref1")
SPECIAL = ("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07")


def frame_hash(frame):
    return hashlib.sha256(pd.util.hash_pandas_object(frame,index=False).to_numpy().tobytes()).hexdigest()


def union_endpoint(c, patients, trait, baseline, birth, censor, valid, source_offsets=None):
    # A source can have a saved and a reconstructed FOD. Scan both. A routine
    # repeat visit after an earlier FOD is not a new source-specific first date.
    cols = {"icd10":("icd10","computed_icd10","official_icd10"),
            "ref":("ref","saved_ref","computed_ref"), "icd9":("icd9","computed_icd9","official_icd9"),
            "gp":("gp","computed_gp","gp2","gp3","computed_gp2","computed_gp3","official_gp2","official_gp3"),
            "self_report":("srd","computed_srd","official_srd"), "death":("death","computed_death","official_death"),
            "procedure":("opcs4","computed_opcs4")}
    idx = c.index
    first = pd.Series(pd.NaT,index=idx,dtype="datetime64[ns]")
    early_first = first.copy()
    latest_source_first = first.copy()
    near_before = first.copy()
    first_sources = pd.Series("",index=idx)
    early_sources = first_sources.copy()
    unresolved = pd.Series(False,index=idx)
    cutoff = baseline + pd.Timedelta(days=365)
    has_source = False
    for name, fields in cols.items():
        values = [patients[k] for k in fields if k in patients]
        orig = f"fod_{name}_{trait}"
        if name in ("icd10","ref") and orig in c:
            values.append(c[orig])
        earliest = pd.Series(pd.NaT,index=idx,dtype="datetime64[ns]")
        early = earliest.copy()
        bad = pd.Series(False,index=idx)
        for value in values:
            dt = pd.to_datetime(value,errors="coerce")
            invalid = dt.isin(pd.to_datetime(SPECIAL)) | dt.lt(birth) | (value.notna() & dt.isna())
            bad |= invalid
            dt = dt.mask(invalid)
            has_source |= bool(dt.notna().any())
            earliest = earliest.mask(dt.notna() & (earliest.isna() | dt.lt(earliest)),dt)
            before = dt.between(baseline-pd.Timedelta(days=365),baseline) & dt.le(censor)
            near_before = near_before.mask(before & (near_before.isna() | dt.gt(near_before)),dt)
            hit = dt.gt(baseline) & dt.le(cutoff) & dt.le(censor)
            early = early.mask(hit & (early.isna() | dt.lt(early)),dt)
        for k in fields:
            for prefix in ("cnt_", "cnt_undated_"):
                col = prefix+k
                if col in patients:
                    evidence = pd.to_numeric(patients[col],errors="coerce").gt(0)
                    bad |= evidence if prefix == "cnt_undated_" else evidence & earliest.isna()
        if "invalid_"+orig in c:
            bad |= c["invalid_"+orig].fillna(False)
        count = f"cnt_{name}_{trait}"
        if count in c:
            bad |= pd.to_numeric(c[count],errors="coerce").gt(0) & earliest.isna()
        if name=="ref" and "ref_status" in patients:
            bad |= patients.ref_status.isin(["undated_cache_evidence","special_or_invalid_date","valid_with_special_component"])
        if source_offsets is not None:
            label="srd" if name=="self_report" else name
            source_offsets["src_"+label+"_days"]=(earliest.where(valid & earliest.le(censor))-baseline).dt.days
            if name=="icd10":source_offsets["src_icd10_undated"]=bad
            if name=="gp":
                # A locally matched dated primary-care record. REF and composite/cached GP never confirm raw GP.
                fields=[k for k in ("computed_gp2","computed_gp3") if k in patients]
                direct=pd.concat([pd.to_datetime(patients[k],errors="coerce") for k in fields],axis=1) if fields else pd.DataFrame(index=idx)
                for k in direct:
                    dt=direct[k]
                    direct[k]=dt.mask(dt.isin(pd.to_datetime(SPECIAL)) | dt.lt(birth) | dt.gt(censor))
                gp=direct.min(axis=1) if fields else pd.Series(pd.NaT,index=idx,dtype="datetime64[ns]")
                source_offsets["src_raw_gp_days"]=(gp.where(valid)-baseline).dt.days
        unresolved |= bad
        # Source disagreement for the fixed 10-year study uses evidence through its horizon.
        seen = earliest.notna() & earliest.le(censor) & earliest.le(baseline+pd.Timedelta(days=3652.5))
        latest_source_first = latest_source_first.mask(seen & (latest_source_first.isna() | earliest.gt(latest_source_first)),earliest)
        earlier = earliest.notna() & (first.isna() | earliest.lt(first))
        tied = earliest.notna() & earliest.eq(first)
        first_sources.loc[earlier] = name
        first_sources.loc[tied] += "|"+name
        first = first.mask(earlier,earliest)
        hit = early.notna()
        early_sources.loc[hit] = (early_sources.loc[hit]+"|"+name).str.strip("|")
        early_first = early_first.mask(hit & (early_first.isna() | early.lt(early_first)),early)
    pre = first.le(baseline)
    state = pd.Series(np.select([~valid, np.full(len(c),not has_source),early_first.notna(),pre,
                       unresolved,censor.le(cutoff),first.gt(cutoff) & first.le(censor)],
                       ["invalid_followup","unavailable","excluded_year1","prevalent","uncertain","short_followup","incident"],
                       default="control"),index=idx)
    clean = state.isin(["incident","control"])
    date = first.where(state.isin(["incident","prevalent"]))
    # Retain valid post-censor source dates for auditing only. Never use those
    # dates as model inputs or put a post-censor diagnosis on an incident case.
    audit = pd.DataFrame(dict(any1=date,any1_state=state,any1_t2e=((first.where(state.eq("incident"),censor)-baseline).dt.days/365.25).where(clean),
        any1_Yt2e=state.eq("incident").astype(float).where(clean),any_first=first,any_first_source=first_sources,
        any_early_sources=early_sources,any_prior=pre,any_undated=unresolved,
        any_early_days=(early_first-baseline).dt.days,any_first_days=(first-baseline).dt.days,
        any_followup_days=(censor-baseline).dt.days),index=idx)
    audit["source_first_span_days"]=(latest_source_first-first).dt.days
    audit["any_near_before_days"]=(near_before-baseline).dt.days
    return audit


def derive(cohort, patients, trait, cutoff):
    for df in (cohort, patients):
        if df.eid.isna().any() or df.eid.duplicated().any():
            raise ValueError("Missing or duplicate participant ID")
    c = cohort.set_index("eid")
    if len(patients) and not patients.eid.isin(c.index).any():
        raise ValueError("AOD and phenotype cohort have no participant IDs in common")
    p = patients.set_index("eid").reindex(c.index)
    def dates(frame, col):
        return pd.to_datetime(frame.get(col, pd.Series(pd.NaT, index=c.index)), errors="coerce")
    baseline, birth = dates(c, "date_attend"), dates(c, "birth_date")
    censor = pd.concat([dates(c, "date_death"), dates(c, "date_lost"),
                        pd.Series(pd.Timestamp(cutoff), index=c.index)], axis=1).min(axis=1)
    valid = baseline.notna() & censor.ge(baseline) & (birth.isna() | baseline.ge(birth))
    def clean(frame, col):
        raw = dates(frame, col)
        invalid = raw.isin(pd.to_datetime(SPECIAL)) | raw.lt(birth)
        return raw.mask(invalid), invalid
    out = pd.DataFrame(index=c.index)
    summary, cross = [], []
    missing = pd.Series(False, index=c.index)
    case_dates, unknown, available = {}, {}, {}
    for key, source in (("icd", "icd10"), ("ref", "ref")):
        field = f"fod_{source}_{trait}"
        case_dates[key], bad = clean(c, field)
        count = pd.to_numeric(c.get(f"cnt_{source}_{trait}", pd.Series(0, index=c.index)), errors="coerce")
        unknown[key] = bad | c.get("invalid_" + field, missing) | (count.gt(0) & case_dates[key].isna())
        available[key] = field in c and case_dates[key].notna().any()
    for key in BASE_ENDPOINTS:
        src = key.rstrip("1")
        d = case_dates[src]
        early = d.gt(baseline) & d.le(baseline + pd.Timedelta(days=365)) & d.le(censor)
        excluded = early if key.endswith("1") else missing
        state = np.select([np.full(len(c), not available[src]), ~valid, unknown[src], excluded,
                           d.le(baseline), d.gt(baseline) & d.le(censor)],
                          ["unavailable", "invalid_followup", "uncertain", "excluded_year1", "prevalent", "incident"],
                          default="control")
        state = pd.Series(state, index=c.index)
        yin = state.isin(["incident", "control"])
        out[key] = d.where(yin | state.eq("prevalent"))
        out[key + "_state"] = state
        out[key + "_Yt2e"] = state.eq("incident").astype(float).where(yin)
        out[key + "_t2e"] = ((d.where(state.eq("incident"), censor) - baseline).dt.days / 365.25).where(yin)
        counts = state.value_counts()
        row = dict(trait=trait, dx=key, cohort_N=len(c), after_censor=int((valid & d.gt(censor)).sum()),
                   same_day=int((state.eq("prevalent") & d.eq(baseline)).sum()))
        for name in ("control", "incident", "prevalent", "excluded_year1", "uncertain", "invalid_followup", "unavailable"):
            row[name] = int(counts.get(name, 0))
        summary.append(row)
    base_endpoint_hash = frame_hash(out.reset_index())
    summary = [dict(row, endpoint_sha256=base_endpoint_hash) for row in summary]
    source_offsets={}
    union = union_endpoint(c, p, trait, baseline, birth, censor, valid,source_offsets)
    out = out.join(union)
    states = union.any1_state.value_counts()
    summary.append(dict(trait=trait, dx="any1", cohort_N=len(c),
        after_censor=int((valid & union.any_first.gt(censor)).sum()),
        same_day=int((union.any1_state.eq("prevalent") & union.any_first.eq(baseline)).sum()),
        **{state:int(states.get(state,0)) for state in ("control","incident","prevalent","excluded_year1",
            "uncertain","short_followup","invalid_followup","unavailable")},
        endpoint_sha256=frame_hash(union.drop(columns=["source_first_span_days","any_near_before_days"]).reset_index())))
    if source_offsets:out=out.join(pd.DataFrame(source_offsets,index=c.index))
    # Recorded offsets for the paired timing study. Post-censor records are audit-only.
    for key in ("icd","ref"):
        dt = case_dates[key]
        out[key+"_days"] = (dt.where(valid & dt.le(censor)) - baseline).dt.days
    for col in ("age","sex","drug.lipid","baseline_lipid_i0","bmi","sbp","dbp","bb_TC","bb_HDL","bb_HBA1C"):
        if col in c:
            out[col] = c[col]
    for i, k in enumerate(DEFINITIONS):
        for j in DEFINITIONS[i+1:]:
            tab = out.groupby([k+"_state", j+"_state"], observed=True).size()
            cross.extend(dict(trait=trait, dx_a=k, dx_b=j, state_a=a, state_b=b, N=int(n))
                         for (a,b),n in tab.items())
    return out.reset_index(), pd.DataFrame(summary), pd.DataFrame(cross)


def main(work, directory, snapshot, trait):
    work, directory, snapshot = map(Path, (work, directory, snapshot))
    datasets = [pd.read_csv(work / name, sep="\t", dtype={"eid": str}, low_memory=False)
                for name in ("cohort.tsv.gz", "patients.tsv.gz")]
    data, counts, cross = derive(*datasets, trait, (work/"followup.txt").read_text().strip())
    snapshot.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=directory) as stage:
        stage = Path(stage)
        f = stage/"phenotypes_LOCAL_ONLY.tsv.gz"
        data.to_csv(f, sep="\t", index=False, date_format="%Y-%m-%d", compression="gzip")
        fingerprint = hashlib.sha256(pd.util.hash_pandas_object(data, index=False).to_numpy().tobytes()).hexdigest()
        counts["short_followup"] = counts["short_followup"].fillna(0).astype(int)
        counts.to_csv(stage/"dx_counts.tsv", sep="\t", index=False)
        cross.to_csv(stage/"dx_transitions.tsv", sep="\t", index=False)
        hs = data.groupby(["any1_state","any_prior","any_undated","any_first_source","any_early_sources"],
                          dropna=False).size().reset_index(name="N")
        hs.to_csv(stage/"any1_sources.tsv",sep="\t",index=False)
        discordant=data.loc[data.icd_state.eq("incident") & data.ref_state.eq("prevalent")]
        if discordant.any1_state.isin(["incident","control"]).any():
            raise ValueError("A documented baseline REF case entered the Any1 risk set")
        audit=discordant.groupby(["icd1_state","any1_state"],dropna=False).size().reset_index(name="N")
        audit["source_N"]=len(discordant)
        audit["pct_of_discordant"]=100*audit.N/max(1,len(discordant))
        audit.assign(trait=trait,dx="any1").to_csv(stage/"any1_discordance.tsv",sep="\t",index=False)
        html = ('<!doctype html><meta charset="utf-8"><title>'+trait+' diagnosis dates</title><h1>'+trait+
                ' · Diagnosis dates</h1><p><a href="date-check.html">AOD disease course and source concordance</a></p>'
                '<p>Any1 excludes a first date from any eligible source on days 1–365, even if a different source records a baseline history. '
                'Other baseline histories count as prevalent. Follow-up must extend beyond day 365 for the incident risk set.</p>'+
                counts.drop(columns="endpoint_sha256").to_html(index=False)+
                '<h2>Matched transitions</h2>'+cross.to_html(index=False)+'<h2>Any1 source evidence</h2>'+hs.to_html(index=False))
        (stage/"date_check.html").write_text(html)
        # Each scheduled batch reads a frozen hard link, even during another preparation.
        f.chmod(0o600)
        with (snapshot.parent/"prepare.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            manifest = snapshot.parent/"manifest.json"
            metadata = {"endpoint_sha256":fingerprint, "counts":counts.to_dict("records")}
            if manifest.is_file() and json.loads(manifest.read_text())["endpoint_sha256"] != fingerprint:
                raise ValueError("Inputs changed during this batch; start a fresh plan.")
            if not manifest.is_file():
                temp_manifest = stage/"manifest.new"
                temp_manifest.write_text(json.dumps(metadata))
                os.replace(temp_manifest, manifest)
            if not snapshot.is_file():
                os.link(f, snapshot)
            (stage/"dx_manifest.json").write_text(json.dumps(metadata))
        # Publish only the canonical results/reports/runtime paths.
        for path in stage.iterdir():
            bucket={"date_check.html":"main/reports","phenotypes_LOCAL_ONLY.tsv.gz":"tmp",
                    "dx_manifest.json":"qc"}.get(path.name,"main/results")
            dest=directory/bucket/path.name
            dest.parent.mkdir(parents=True,exist_ok=True)
            os.replace(path,dest)


if __name__ == "__main__":
    main(*sys.argv[1:])
