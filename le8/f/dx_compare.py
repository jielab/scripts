"""Aggregate existing diagnosis-date runs; never open participant-level data."""
import json
import os
from pathlib import Path
import sys
import tempfile

import pandas as pd
from dx_layout import PRIMARY_DX, scope_for_read
from c1_report import reliability_summary_row, reliability_tables


TIMING = ("timing_groups","timing_anchor_validation","timing_scores","timing_medication","timing_biomarkers","timing_paired","timing_status")


INACTIVE = {"RUN","RUNNING","EXIT","FAILED","INTERRUPTED","SKIPPED","BLOCKED"}


def record(path):
    return json.loads(path.read_text()) if path.is_file() else {}


def execution_status(root, phase):
    state = phase.get("status","")
    job_id = str(phase.get("job_id",""))
    if state not in {"RUN","RUNNING"} or not job_id.isdigit():
        return state
    # A killed worker may leave RUNNING on disk. Only its own LSF accounting
    # footer establishes a terminal state. Absence of accounting changes nothing.
    for log in (root/"logs/scheduler").glob(f"*/**/*.{job_id}.out"):
        with log.open(errors="replace") as handle:
            accounting = False
            for line in handle:
                if line.startswith(f"Subject: Job {job_id}:"):
                    accounting = True
                if accounting and line.startswith("TERM_OWNER:"):
                    return "INTERRUPTED"
                if accounting and line.startswith("Exited with exit code "):
                    return "FAILED"
    return state


def report(root, traits, layers, dates):
    buckets = {k: [] for k in ("availability", "cohort_counts", "assay_counts", "transitions", "performance",
                               "coverage", "model_contrasts", "ABM_support", "selector_tests", "model_status", "glmnet_path", "any1_discordance", "any1_sources", *TIMING, "ham_sources", "historical_borrowing", "LDL_interaction", "LDL_slopes","glmnet_neighbors","glmnet_neighbor_status")}
    for trait in traits:
        aod = root/trait/"aod/main/results"
        counts = {}
        for name, file in (("cohort_counts", "dx_counts.tsv"), ("transitions", "dx_transitions.tsv")):
            if (aod/file).is_file():
                f = pd.read_csv(aod/file, sep="\t")
                if name == "transitions":
                    f = f[f.dx_a.isin(dates) & f.dx_b.isin(dates)]
                else:
                    f = f[f.dx.isin(dates)]
                    counts = f.set_index("dx").endpoint_sha256.to_dict()
                buckets[name].append(f)
        if "any1" in dates and (aod/"any1_sources.tsv").is_file():
            buckets["any1_sources"].append(pd.read_csv(aod/"any1_sources.tsv",sep="\t").assign(trait=trait,dx="any1"))
        if "any1" in dates and (aod/"any1_discordance.tsv").is_file():
            buckets["any1_discordance"].append(pd.read_csv(aod/"any1_discordance.tsv",sep="\t"))
        for layer in layers:
            for dx in dates:
                base = scope_for_read(root/trait/layer,dx)
                scope = base/"c1_correlate"
                lineage = record(base/"qc/dx.json")
                receipt = record(scope/"qc/C1_FINAL_INDEX.json")
                phases = [record(base/"logs"/(k+".status.json")) for k in ("analysis","prepare","primary","evaluation")]
                phases = [p for p in phases if p and (not p.get("steps") or
                          any(".c1." in s["key"] or "c2_cause" in s["key"] for s in p.get("steps",[])))]
                latest = max(phases, key=lambda p: p.get("started",""), default={})
                state = execution_status(root,latest)
                ok = bool(counts.get(dx) and lineage.get("dx") == dx and lineage.get("trait") == trait
                          and counts[dx] == lineage.get("endpoint_sha256"))
                label = dict(trait=trait, layer=layer, dx=dx)
                status = dict(label, execution=state or "no_scientific_status",
                              prediction=receipt.get("status","missing"), input_matches_aod=ok if dx != "default" else None,
                              abm_deep_status=("included" if "abm_deep" in receipt.get("methods", {}) else "optional_not_available"),
                              abm_deep_reason=receipt.get("optional_skipped_methods", {}).get("abm_deep", {}).get("reason", ""))
                buckets["availability"].append(pd.DataFrame([status]))
                if dx=="any1" and counts.get(dx) and (scope/"timing.xlsx").is_file():
                    valid_dates = set()
                    for ref in dates:
                        ref_base = scope_for_read(root/trait/layer,ref)
                        info = record(ref_base/"qc/dx.json")
                        updates = [record(ref_base/"logs"/(phase+".status.json")) for phase in
                                   ("analysis","prepare","primary","evaluation")]
                        ref_state = max(updates,key=lambda r:r.get("started","")).get("status")
                        if (counts.get(ref) and info.get("endpoint_sha256")==counts[ref] and
                                ref_state not in INACTIVE and not info.get("unavailable")):
                            valid_dates.add(ref)
                    with pd.ExcelFile(scope/"timing.xlsx") as timing:
                        for name in TIMING:
                            if name not in timing.sheet_names:
                                continue
                            local = timing.parse(name)
                            required = {"endpoint_sha256","trait","layer","dx"}
                            if (local.empty or not required <= set(local) or not
                                    all(local[k].eq(v).all() for k,v in {**label,"endpoint_sha256":counts[dx]}.items())):
                                continue
                            # Anchors and controlled training depend on the AOD, not C1 fits.
                            if name in ("timing_paired","timing_status") and "analysis" in local:
                                saved = local.analysis.eq("saved_C1_common_test")
                                refs = local.get("reference",pd.Series("",index=local.index)).astype(str).str.split(".").str[0]
                                current = refs.isin(valid_dates) & ("any1" in valid_dates)
                                if "reference_endpoint_sha256" in local:
                                    current &= local.reference_endpoint_sha256.eq(refs.map(counts))
                                local = local.loc[~saved | current]
                            if not local.empty:
                                buckets[name].append(local)
                if (not ok and dx != "default") or state in INACTIVE or lineage.get("unavailable"):
                    print(f"[LE8] SKIP {trait}.{layer}.{dx}: execution={state or 'unknown'}, input_matches_aod={ok}, unavailable={bool(lineage.get('unavailable'))}", flush=True)
                    continue
                assay = record(scope/"qc/yinyang.json")
                if assay and assay.get("dx") == dx and assay.get("endpoint_sha256") == counts.get(dx):
                    buckets["assay_counts"].append(pd.DataFrame([{**assay,**label}]))
                elif dx == "default" and (scope/"c1.cohort.xlsx").is_file():
                    book = scope/"c1.cohort.xlsx"
                    with pd.ExcelFile(book) as workbook:
                        if "results" in workbook.sheet_names:
                            cohort = workbook.parse("results")
                            buckets["assay_counts"].append(cohort.assign(**label))
                glm = receipt.get("methods",{}).get("glmnet") or receipt.get("unavailable_methods",{}).get("glmnet") or {}
                if glm.get("path"):
                    version = scope/"tmp/models/glmnet"/Path(glm["path"]).name
                    marker = record(version/"COMPLETE.json") or record(version/"UNAVAILABLE.json")
                    signature = glm.get("stored_signature") or glm.get("signature")
                    fit = record(version/"glmnet_fit.json") if signature and marker.get("signature")==signature else {}
                    if fit:
                        fields = ("status","solver_jerr","failed_lambda_index","lambda_requested","lambda_count",
                                  "selected_index","lambda","maxit_used","coordinate_passes","solver_seconds","fit_attempts","tune_logloss","reason")
                        row = {k:fit.get(k) for k in fields}
                        row["path_status"] = fit["path_status"]
                        buckets["glmnet_path"].append(pd.DataFrame([{**label,**row}]))
                entry = receipt.get("methods",{}).get("abm_prevalent",{})
                if isinstance(entry,dict) and entry.get("path"):
                    models = scope/"tmp/models/abm_prevalent"
                    saved = Path(entry["path"])
                    # Use only the version selected by the final index, also after copying/migration.
                    version = models/saved.name
                    finished = record(version/"COMPLETE.json") if saved.name not in ("",".","..") else {}
                    signature = entry.get("stored_signature") or entry.get("signature")
                    if signature and finished.get("signature") == signature:
                        meta = record(version/"MODEL_FROZEN.json")
                        if meta.get("method") == "abm_prevalent":
                            audit = meta.get("source_transform",{})
                            row = {k:meta.get(k) for k in ("alpha","history_use_status","fit_build_N",
                                   "fit_build_horizon_cases","max_calibrated_probability_change","historical_rows_in_incident_loss")}
                            row.update({k:audit.get(k) for k in ("source_offered_N","source_used_N","source_families","status","reason")})
                            buckets["historical_borrowing"].append(pd.DataFrame([{**row,**label}]))
                if receipt.get("evaluation") and (scope/"qc/evaluation/COMPLETE.json").is_file():
                    for name, filename in (("performance","performance"), ("coverage","coverage"), ("model_contrasts","contrasts"), ("model_status","evaluation_status"),
                            ("glmnet_neighbors","glmnet_neighbors"),("glmnet_neighbor_status","glmnet_neighbor_status")):
                        path = scope/"qc/evaluation"/(filename+".tsv")
                        if path.is_file():
                            try:
                                rows = pd.read_csv(path, sep="\t")
                            except pd.errors.EmptyDataError:
                                continue
                            buckets[name].append(rows.assign(**label))
                    tables, meta = reliability_tables(scope)
                    if tables:
                        for table in ("common_riskmatched_selection","reliability_increment_50"):
                            tests = tables.get(table, pd.DataFrame())
                            if {"primary_risk_model","model","reference"} <= set(tests):
                                tests = tests[tests.primary_risk_model.eq(True) & (
                                    (tests.model.eq("raw_count_support") & tests.reference.eq("raw_count_support_riskmatched")) |
                                    (tests.model.eq("P1_incident_support") & tests.reference.eq("P0_risk_only")) |
                                    (tests.model.eq("P2_incident_plus_prevalent") & tests.reference.eq("P1_incident_support")))]
                                if not tests.empty:
                                    buckets["selector_tests"].append(tests.assign(**label))
                        row = reliability_summary_row(trait,tables,meta)
                        buckets["ABM_support"].append(pd.DataFrame([{**row,**label}]))
                if dx in (PRIMARY_DX,"default") and layer == "met" and (base/"c2_cause").is_dir():
                    data = record(base/"c2_cause/c2.ldl_interaction.json")
                    for name, key in (("LDL_interaction","tests"),("LDL_slopes","slopes")):
                        rows = pd.DataFrame(data.get(key) or [])
                        if key == "tests" and rows.empty:
                            rows = pd.DataFrame([{"status":data.get("status","not_available")}])
                        if len(rows):
                            buckets[name].append(rows.assign(**label, definition_role="primary_ICD", calibration_partial_R2=data.get("calibration_partial_R2")))
    return {k:pd.concat(frames, ignore_index=True) if frames else pd.DataFrame() for k,frames in buckets.items()}


def main(root, traits, layers, dates):
    root = Path(root)
    from file_lock import file_lock
    with file_lock(root/"qc/dx_comparison.lock"):
        return write_comparison(root,traits,layers,dates)


def write_comparison(root, traits, layers, dates):
    archive = root/"main/results/dx_comparison.xlsx"
    if archive.is_file():
        previous = pd.read_excel(archive,sheet_name="availability",usecols=["trait","layer","dx"])
        selections = [list(dict.fromkeys([*old.dropna(),*new.split(",")])) for old,new in
                      ((previous.trait,traits),(previous.layer,layers),(previous.dx,dates))]
        # Same validity rules as the paper reader, including rejection of stale fits
        # when a newer failed/unfinished scope exists. Old archive-only scopes survive.
        from paper_summary import load_inputs
        tables = load_inputs(root,*selections)
        tables.pop("module_files",None)
    else:
        tables = report(root, traits.split(","), layers.split(","), dates.split(","))
    out = root/"main/results"
    reports = root/"main/reports"
    reports.mkdir(parents=True,exist_ok=True)
    out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out) as tmp:
        tmp = Path(tmp)
        book = tmp/"dx_comparison.xlsx"
        with pd.ExcelWriter(book,engine="openpyxl") as writer:
            for k, df in tables.items():
                df.to_excel(writer,sheet_name=k,index=False)
        body = ['<!doctype html><meta charset="utf-8"><title>Diagnosis-date comparison</title>',
                '<h1>Diagnosis-date comparison</h1><p><a href="../results/dx_comparison.xlsx">Full tables</a></p>',
                '<p>Counts use the full phenotype cohort. Performance uses the test population of each definition. '
                'Cross-definition changes are descriptive. The five definitions share participants. '
                'Historical borrowing is evaluated against its matched no-history prediction model. P2 is unavailable in the CVD primary selector protocol.</p>']
        for name in ("cohort_counts","assay_counts","availability","performance","ABM_support","LDL_interaction","LDL_slopes"):
            df = tables[name]
            if name == "performance" and "subset" in df:
                df = df.loc[df.subset.eq("all_common")]
                cols = ["trait","layer","dx","model","horizon","evaluated_N","cases","controls","C_index","Cindex","AUC","Brier","Harrell_C","AUC_IPCW","Brier_IPCW","LogLoss_IPCW","coverage_all"]
                df = df[[c for c in cols if c in df]]
            if name in {"cohort_counts","assay_counts"}:
                df = df.drop(columns="endpoint_sha256",errors="ignore")
            body.extend([f"<h2>{name}</h2>", df.to_html(index=False,na_rep="") if not df.empty else "<p>No eligible results.</p>"])
        (tmp/"dx_comparison.html").write_text("\n".join(body),encoding="utf-8")
        for trait in traits.split(","):
            for layer in layers.split(","):
                directory = root/trait/layer/"main/reports"
                directory.mkdir(parents=True,exist_ok=True)
                local = [f'<!doctype html><meta charset="utf-8"><title>{trait} / {layer}</title>',
                         f'<h1>{trait} / {layer}</h1>',
                         '<p><a href="../../../../main/reports/dx_comparison.html">All selected scopes</a></p>']
                for key in ("availability","assay_counts","ABM_support"):
                    frame = tables[key]
                    if not frame.empty:
                        frame = frame[frame.trait.eq(trait) & frame.layer.eq(layer)]
                        if not frame.empty:
                            local += [f"<h2>{key}</h2>",frame.drop(columns="endpoint_sha256",errors="ignore").to_html(index=False,na_rep="")]
                with tempfile.NamedTemporaryFile("w",dir=directory,delete=False,encoding="utf-8") as handle:
                    handle.write("\n".join(local))
                os.replace(handle.name,directory/"index.html")
        for f in tmp.iterdir():
            os.replace(f, (reports if f.suffix == ".html" else out)/f.name)


if __name__ == "__main__":
    main(*sys.argv[1:])
