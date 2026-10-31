#!/usr/bin/env python3
"""Read-only manuscript summary for GRID and upstream benchmark scripts.

No individual-level table is exported. No model is fitted here.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import pandas as pd
from final_outputs import aggregate, read_aggregate, current_abm_run
from final_compare import current_tables

PRIMARY = 0.5


def read(path):
    return read_aggregate(path)


def current_blocknet(root: Path, trait: str):
    return current_abm_run(root / "abm/FINAL" / trait / "blocknet.CURRENT.json")


def at_primary(path: Path):
    d = read(path)
    col = "coverage.nominal_development_coverage"
    if col not in d:
        raise ValueError(f"{path}: missing {col}")
    z = d[(pd.to_numeric(d[col], errors="coerce") - PRIMARY).abs() < 1e-12].copy()
    if len(z) != 1:
        raise ValueError(f"{path}: expected one primary-coverage row, got {len(z)}")
    return z.iloc[0]


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--analysis-root", default="/mnt/d/analysis/grid")
    p.add_argument("--traits", default="height,ldl,t2dm")
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args(argv)
    root = Path(a.analysis_root)
    out = root / "final"
    traits = [x.strip() for x in a.traits.split(",")]
    if not traits or len(set(traits)) != len(traits) or not set(traits) <= {"height", "ldl", "t2dm"}:
        p.error("Use unique traits from height,ldl,t2dm")

    if a.dry_run:
        print("[GRID final] read-only paper summary; no fitting")
        return

    # Revalidate the current producers even when invoked on its own; stale final
    # tables must not resurrect an old benchmark or a removed CURRENT pointer.
    verified_tables, verified_status = current_tables(root, traits)
    perf, pred, benefit, contrasts, key, status = [], [], [], [], [], []

    for trait in traits:
        run = current_blocknet(root, trait)
        if run is None:
            status.append(dict(component="GRID", trait=trait, status="not_completed",
                               source_script="grid.sh"))
            continue

        d = read(run / "main_performance.tsv")
        d["source_script"] = d["method"].map(lambda m:
            "2.csx.sh + GRID validation" if m == "csx-linear"
            else "grid.sh")
        d["scientific_role"] = d["method"].map(lambda m:
            "PRS-CSx aggregate-score baseline" if m == "csx-linear"
            else "GRID LD-block candidate" if m == "csx-blocknet"
            else "GRID selective policy" if str(m).startswith("ABM(")
            else "benchmark")
        d["method_class"] = d.method.map(lambda m: "baseline" if m == "csx-linear" else
            "block_model" if m == "csx-blocknet" else "policy" if str(m).startswith("ABM(") else "benchmark")
        d["current_run"] = run.name
        perf.append(d)

        missing = [name for name in ("predictability_curve.tsv", "benefit_curve.tsv")
                   if not (run / name).is_file()]
        if missing:
            status.append(dict(component="GRID assessment", trait=trait, status="not_completed",
                               source_script="grid.sh", run_id=run.name,
                               detail="Missing aggregate reports: " + ", ".join(missing)))
        pr = at_primary(run / "predictability_curve.tsv") if not missing else pd.Series(dtype=object)
        be = at_primary(run / "benefit_curve.tsv") if not missing else pd.Series(dtype=object)
        pred.append(dict(
            trait=trait, run_id=run.name, nominal_coverage=PRIMARY,
            actual_coverage=pr.get("coverage.actual_selected_coverage_query"),
            selected_base_loss=pr.get("matched.selected_base_loss.estimate"),
            control_base_loss=pr.get("matched.control_base_loss.estimate"),
            selected_minus_control_loss=pr.get("matched.predictability_loss_difference.estimate"),
            ci_low=pr.get("matched.predictability_loss_difference.ci_low"),
            ci_high=pr.get("matched.predictability_loss_difference.ci_high"),
            n_pairs=pr.get("matching.n_pairs"),
            matching_fraction=pr.get("matching.matching_fraction_of_selected"),
            grouping_status=pr.get("grouping_status"),
        ))
        benefit.append(dict(
            trait=trait, run_id=run.name, nominal_coverage=PRIMARY,
            actual_coverage=be.get("coverage.actual_selected_coverage_query"),
            selected_gain=be.get("matched.selected_gain_vs_base.estimate"),
            control_gain=be.get("matched.control_gain_vs_base.estimate"),
            gain_enrichment=be.get("matched.gain_enrichment.estimate"),
            ci_low=be.get("matched.gain_enrichment.ci_low"),
            ci_high=be.get("matched.gain_enrichment.ci_high"),
            n_pairs=be.get("matching.n_pairs"),
            matching_fraction=be.get("matching.matching_fraction_of_selected"),
            grouping_status=be.get("grouping_status"),
        ))

        if (run / "count_incremental_contrasts.tsv").is_file():
            c = read(run / "count_incremental_contrasts.tsv")
            c["trait"], c["run_id"] = trait, run.name
            contrasts.append(c)

        base = d[d.method.eq("csx-linear")]
        block = d[d.method.eq("csx-blocknet")]
        if len(base) == len(block) == 1:
            if trait == "t2dm":
                key += [
                    dict(trait=trait, result="blocknet_vs_csx",
                         metric="brier", baseline=float(base.iloc[0].brier),
                         GRID=float(block.iloc[0].brier),
                         delta=float(base.iloc[0].brier-block.iloc[0].brier),
                         direction="positive delta favours GRID"),
                    dict(trait=trait, result="blocknet_vs_csx",
                         metric="auprc", baseline=float(base.iloc[0].auprc),
                         GRID=float(block.iloc[0].auprc),
                         delta=float(block.iloc[0].auprc-base.iloc[0].auprc),
                         direction="positive delta favours GRID"),
                ]
            else:
                key.append(dict(trait=trait, result="blocknet_vs_csx",
                    metric="prediction_r2", baseline=float(base.iloc[0].prediction_r2),
                    GRID=float(block.iloc[0].prediction_r2),
                    delta=float(block.iloc[0].prediction_r2-base.iloc[0].prediction_r2),
                    direction="positive delta favours GRID"))
        for result, prefix, curve, metric in [
            ("predictability_selected_vs_control", "matched.predictability_loss_difference", pr, "baseline_prediction_loss_difference"),
            ("gain_enrichment", "matched.gain_enrichment", be, "candidate_gain_difference"),
            ("policy_vs_base", "whole_cohort.policy_gain_vs_base", be, "policy_loss_reduction_vs_base"),
            ("policy_vs_candidate", "whole_cohort.policy_gain_vs_candidate", be, "policy_loss_reduction_vs_candidate")]:
            estimate = curve.get(prefix + ".estimate")
            key.append(dict(trait=trait, result=result, metric=metric, delta=estimate,
                ci_low=curve.get(prefix + ".ci_low"), ci_high=curve.get(prefix + ".ci_high"),
                status="available" if pd.notna(estimate) else "unavailable",
                direction="negative favours selected" if result.startswith("predictability") else "positive favours GRID"))
        if (run / "count_incremental_contrasts.tsv").is_file():
            for row in read(run / "count_incremental_contrasts.tsv").query("contrast == 'D_vs_C'").to_dict("records"):
                key.append(dict(trait=trait, result="count_incremental_D_vs_C", metric="assessment_target_mse_reduction",
                    objective=row.get("objective"), baseline=row.get("comparator_target_mse"), GRID=row.get("augmented_target_mse"),
                    delta=row.get("mse_reduction"), ci_low=row.get("ci_low"), ci_high=row.get("ci_high"),
                    status=row.get("status"), direction="positive favours adding count"))
        else:
            key.append(dict(trait=trait, result="count_incremental_D_vs_C", metric="assessment_target_mse_reduction", status="unavailable"))
        status.append(dict(component="GRID", trait=trait, status="completed",
                           source_script="grid.sh", run_id=run.name))

    roles = {"PRSformer": "independent-trait adaptation",
             "PRSformer2 exploratory": "exploratory adaptation with separate calibration",
             "DiscoDivas": "same-roster development-only outcome adaptation"}
    for item in verified_status:
        method = item.get("method")
        if method in roles:
            status.append(dict(component=method, trait=item["trait"], status=item["status"],
                source_script="3.disco.sh" if method.startswith("Disco") else "4.prsformer.sh",
                scientific_role=roles[method], detail=item.get("reason")))
    for table in verified_tables:
        ext = table[table.method.isin(roles)].copy()
        if len(ext):
            ext["source_script"] = ext.method.map(lambda m: "3.disco.sh" if m.startswith("Disco") else "4.prsformer.sh")
            ext["scientific_role"] = ext.method.map(roles)
            ext["method_class"], ext["current_run"] = "benchmark", ext.run_id
            perf.append(ext)
    from kinship import active
    correction = active(root)
    if correction:
        status.append(dict(component="kinship correction", trait="all", status="active",
                           source_script="grid.sh + 4.prsformer.sh",
                           detail="All test IDs retained; development relatives excluded; original matrix reused"))
    status.append(dict(component="PCA", trait="all", status="prerequisite_QC",
                       source_script="1.pca.sh", scientific_role="prerequisite/QC, not a competing predictor"))

    from kinship_qc import public_qc
    qc = public_qc(root)
    # Every invocation replaces all six tables, including header-only outputs,
    # so no previous run's article data survives a missing CURRENT result.
    from regional_dimensionality import paper_table
    outputs = {
        "paper_regional_dimensionality": paper_table(root, traits),
        "paper_main_performance": (pd.concat(perf, ignore_index=True) if perf else
                                   pd.DataFrame(columns=["trait", "method", "source_script", "scientific_role", "current_run"])),
        "paper_predictability": pd.DataFrame(pred, columns=[
            "trait", "run_id", "nominal_coverage", "actual_coverage", "selected_base_loss", "control_base_loss",
            "selected_minus_control_loss", "ci_low", "ci_high", "n_pairs", "matching_fraction", "grouping_status"]),
        "paper_benefit": pd.DataFrame(benefit, columns=[
            "trait", "run_id", "nominal_coverage", "actual_coverage", "selected_gain", "control_gain",
            "gain_enrichment", "ci_low", "ci_high", "n_pairs", "matching_fraction", "grouping_status"]),
        "paper_count_contrasts": (pd.concat(contrasts, ignore_index=True) if contrasts else
                                  pd.DataFrame(columns=["trait", "run_id", "objective", "contrast", "status"])),
        "paper_key_results": pd.DataFrame(key, columns=[
            "trait", "result", "metric", "objective", "baseline", "GRID", "delta", "ci_low", "ci_high", "status", "direction"]),
        "paper_method_status": pd.DataFrame(status),
        "kinship_qc": qc,
    }
    for table in outputs.values():
        aggregate(table)
    out.mkdir(parents=True, exist_ok=True)
    for name, table in outputs.items():
        temporary = out / (name + ".tsv.tmp")
        table.to_csv(temporary, sep="\t", index=False)
        temporary.replace(out / (name + ".tsv"))
    print(f"[GRID final] article tables written to {out}")


if __name__ == "__main__":
    main()
