#!/usr/bin/env python3
"""Render measured HERV analyses with explicit estimands and availability."""
from __future__ import annotations
import argparse
import html
import json
import os
from pathlib import Path
from analysis_integrity import validate_manifest, identity
from herv_modern import rows
from herv_primary import select_primary_rows


def table_html(data, columns, limit=30):
    data = list(data)
    if not data:
        return "<p class='muted'>No estimable rows available.</p>"
    head = "".join(f"<th>{html.escape(c.replace('_',' '))}</th>" for c in columns)
    body = "".join("<tr>"+"".join(f"<td>{html.escape(str(r.get(c,'')))}</td>" for c in columns)+"</tr>" for r in data[:limit])
    return f"<div class='table'><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>"


def load_json(path):
    return json.loads(Path(path).read_text()) if Path(path).is_file() else {}


def recorded_output(base, name, result):
    """Only files from a validated current manifest may be shown."""
    return name in result.get('output_inventory', {}) and (base / name).is_file()


def make_report(out, modern, analysis_root, pca_temporal_dir=None):
    out, modern, root = Path(out), Path(modern), Path(analysis_root)
    out.mkdir(parents=True, exist_ok=True)
    sections = ["<h1>ERV1 / ERVK and multi-ALT variation</h1><p class='lead'>Measured results, analysis scope and completion status.</p>"]
    sections.append("<h2>1. Primary family enrichment</h2>")
    primary_dir = out / "primary"
    primary_available = False
    pipeline_state = load_json(out / "pipeline_status.json")
    primary_state = next((r for r in pipeline_state.get("stages", []) if isinstance(r, dict) and r.get("stage") == "primary"), {})
    try:
        if primary_state.get("status") in {"failed", "not_available", "blocked", "skipped"}:
            raise ValueError("Primary stage did not complete in this invocation")
        primary = validate_manifest(primary_dir, "primary_manifest.json", require_evidence=True)
        for expected in primary.get("source_counter_inputs", []):
            if identity(expected["path"]) != expected:
                raise ValueError("Primary counter inputs changed; recompute the primary stage")
        selected = select_primary_rows(rows(primary_dir / "enrichment_summary.tsv"))
        primary_available = True
        sections.append(table_html(selected, ["label", "threshold", "observed_multi", "expected_multi", "oe", "ci_low", "ci_high", "fdr", "status"]))
        if recorded_output(primary_dir, "primary_enrichment.png", primary):
            sections.append("<img src='primary/primary_enrichment.png' alt='Primary ERV family enrichment'>")
        sections.append("<p><a href='primary/enrichment_summary.tsv'>Complete primary result table</a> · <a href='primary/primary_manifest.json'>Analysis definition and input provenance</a></p>")
        details = {k: primary.get(k) for k in ("estimand", "target_scope", "background_scope", "minimum_training_sites", "bootstrap_replicates", "seed", "missing_rule")}
        sections.append("<pre>" + html.escape(json.dumps(details, ensure_ascii=False, indent=2)) + "</pre>")
    except (OSError, ValueError, KeyError) as exc:
        sections.append("<p class='muted'>Primary estimate unavailable: " + html.escape(str(exc)) + ". Other completed analyses remain reportable.</p>")
    sections.append("<p>AC5/AC10 requires two distinct nonreference ALT bases to meet the specified count threshold. Reference-base observation is established separately by joint genotypes. This SNV-denominator statistic is not a mutation rate per callable genomic base.</p>")
    status = load_json(out/"pipeline_status.json")
    stages = status.get("stages", [])
    if isinstance(stages, dict):
        stages = [dict(stage=k, **v) for k,v in stages.items()]
    stage_status = {r.get("stage"): r for r in stages}
    if stages:
        sections += ["<h2>Execution status</h2>", table_html(stages, ["stage","status","result_status","reason"], 30)]
    evidence = [
        ("archaic", "archaic_enrichment.tsv", "Archaic tract enrichment in GU cohorts",
         "Finalized IBDmix cohorts use each chromosome scope's analyzed sample panel including zero-tract people. Autosomes retain LOCO and EUR LD-block bootstrap, with BH FDR across selected cohorts. X/Y are separate descriptive O/E comparisons against matched loci outside each label on the same chromosome, including sites without LD blocks. Existing counters do not resolve PAR/non-PAR or sex strata: X/Y have no P, CI or FDR and are not significance tests. Cohorts may overlap; carrier prevalence is not archaic allele frequency.", "archaic_enrichment.png"),
        ("spatial", "primary_spatial_partitions.tsv", "2. Spatial distribution within annotated ERVs",
         "The primary estimator is decomposed into strand-aware aligned-base head, body and tail partitions, with short elements and ambiguous overlaps retained. The end length and normalized-bin count are recorded in the spatial manifest. Partitions use expected-count weights and the same target domain as their total. External flanks are separate. Joined repeat records do not necessarily represent complete proviruses or TSS.", "primary_spatial_partitions.png"),
        ("summary", "common_support_enrichment.tsv", "5. Common-support sensitivity analysis",
         "Unsupported background strata are removed from both observed and expected counts. This sensitivity estimand has a different support rule from the primary full-scope estimator.", "family_enrichment.png"),
        ("windows", "frozen_peak_validation.tsv", "6. Family-specific observed-boundary profiles",
         "Strand-aware windows use the width, step and flank recorded in the windows manifest. The central window-to-distal-flank ratio tests positional concentration. Peak location is selected on discovery chromosomes, then evaluated on held-out chromosomes. These alignment edges are not transcription start sites.", "boundary_profiles.png"),
        ("genotypes", "cohort_retention.tsv", "7. Multi-person genotype concordance",
         "Actual A/C/G/T counts establish whether three observed bases coexist in sampled populations. Results are conditional on the modern candidate set and are not a new genome-wide enrichment estimate. Sample-independent replication requires verified non-overlap with gnomAD.", "genotype_concordance.png"),
        ("ancient", "fixed_panel_family_trajectories.tsv", "3. Full ancient cohort: modern-defined alleles over the past 10,000 years",
         "The full gnomAD HERV catalog, with no default 500-site cap, is followed across every original PSAM sample. Primary, secondary and additional ALT bases are audited separately. Primary time models use officially unrelated ancient individuals aged 0–10,000 BP; modern EUR samples remain a separate comparison group. Exact represented nucleotides are counted by their native PGEN allele indices. A reference base cannot be inferred as two minus one ALT when split records have unknown other-ALT encoding. Unmeasurable bases remain unavailable. Shared PCs and available region, data-source and coverage covariates support descriptive two-sided temporal associations, not selection coefficients or proof of historical tri-allelic occupancy.", "family_allele_time.png"),
        ("selection", "published_selection_summary.tsv", "4. Published exact-allele selection evidence",
         "Published Akbari 2026 statistics are attached to their exact original REF/ALT contrasts. The modern catalog is fixed before this join, and no coefficient is transferred to an unrepresented third base. Raw S, SE, calibrated X, P_X, posterior and FDR are preserved. Signed effects are not reoriented under the default raw_only contract. Counts are descriptive base identities, not independent loci or a new enrichment test; the ancient data overlap the source study.", None),
        ("ld", "association.tsv", "8. HERV and population LD structure",
         "Common physical windows and the intersection of the four ancestry maps supply a shared denominator. Chromosome bootstrap preserves paired regional comparisons. Association is not evidence that HERVs drive LD, especially without sequence and recombination covariates.", "ld_association.png"),
    ]
    evidence_order = {name:i for i,name in enumerate(("spatial", "ancient", "selection", "archaic", "summary", "windows", "genotypes", "ld"))}
    evidence.sort(key=lambda item:evidence_order[item[0]])
    for directory, filename, title, description, figure in evidence:
        stage_name = "time" if directory == "ancient" else directory
        current = stage_status.get(stage_name, {})
        manifest_name = {"summary":"summary_manifest.json", "spatial":"spatial_manifest.json", "windows":"windows_manifest.json",
                         "genotypes":"genotypes_manifest.json", "ancient":"manifest.json", "selection":"manifest.json", "ld":"manifest.json", "archaic":"manifest.json"}[directory]
        # Read the canonical stage output only.
        candidates = [out/"time"] if directory == "ancient" else [out/directory]
        recorded = Path(current["result_manifest"]) if current.get("result_manifest") else None
        base = next((p for p in candidates if recorded is not None and (p/manifest_name).resolve() == recorded.resolve()), None)
        if base is None:
            base = next((p for p in candidates if (p/manifest_name).is_file() or (p/filename).is_file()), candidates[0])
        sections += [f"<h2>{title}</h2><p>{description}</p>"]
        path = base/filename
        if current.get("status") in {"failed", "not_available", "error", "blocked", "skipped"}:
            sections.append("<p class='muted'>This stage did not complete in the current run: "
                            + html.escape(str(current.get("status"))) + ". "
                            + html.escape(str(current.get("reason", "See execution status and stage log."))) + "</p>")
            if path.is_file():
                sections.append("<p class='muted'>Previous output exists but is not displayed as current evidence. "
                                "Its inputs or parameters may differ from this run.</p>")
            continue
        try:
            result = validate_manifest(base, manifest_name, require_evidence=True)
        except (OSError, ValueError) as exc:
            sections.append("<p class='muted'>Current evidence unavailable: " + html.escape(str(exc)) + "</p>")
            continue
        complete_states = {"complete", "complete_descriptive_imputed_panel", "complete_descriptive_all_allele_panel", "complete_adjustment_unavailable"}
        # A direct module rerun can invalidate its manifest after an earlier
        # successful pipeline run. Either record of incompletion blocks display.
        result_state = next((s for s in (result.get("status"), current.get("result_status")) if s and s not in complete_states), None)
        if recorded is not None and not (base/manifest_name).is_file():
            result_state = "missing_result_manifest"
        if result_state:
            sections.append("<p class='muted'>The stage audit did not produce a usable effect estimate: "
                            + html.escape(str(result_state)) + ". Previous measurement files, if present, are not displayed.</p>")
            continue
        if recorded_output(base, filename, result):
            if current.get("status") not in {"complete", "reused"}:
                sections.append("<p class='muted'>Previously saved output. This stage was not validated in this invocation; "
                                "the table below is not a new measurement from the current run.</p>")
            data = list(rows(path))
            if directory == "ancient":
                sections.append("<p>Genetic adjustment: " + html.escape(str(result.get('genetic_adjustment', 'unverified'))) + ". Raw family curves are not PC-adjusted.</p>")
                sections.append("<p>Actual full-catalog execution counts: <code>" + html.escape(json.dumps(
                    {k:result.get(k) for k in ("n_PSAM_samples", "sample_age_scopes", "catalog_sites", "exact_HQ_ALT_contrasts", "n_fixed_contrasts", "n_tested_primary_strict_AC5", "maximum_genotype_chunk_rows") if k in result}, ensure_ascii=False)) + "</code>. See the complete sample and allele-availability audits below.</p>")
            if directory == "summary":
                data = [r for r in data if r.get("label", "").startswith("family:") and r.get("sensitivity") == "exclude_combined"]
            if directory == 'archaic':
                data = [r for r in data if r.get('threshold') == '10' and r.get('sensitivity') == 'exclude_combined' and r.get('label', '').endswith(':any')]
                if recorded_output(base, 'datasets.tsv', result):
                    sections.append(f"<p><a href='{str((base/'datasets.tsv').relative_to(out))}'>Cohort inputs, sample counts, and unavailable cohorts</a></p>")
            if directory == "spatial":
                data = [r for r in data if r.get("scope") == "eur_ld_partition" and r.get("sensitivity") in {"exclude_combined", "strict"}]
            if directory == "ld":
                data = data[:20]
            relative = path.relative_to(out)
            sections.append(f"<p><a href='{html.escape(str(relative))}'>{html.escape(filename)}</a></p>")
            if data:
                columns = (["family","threshold","partition","observed_multi","expected_multi","oe","ci_low","ci_high","expected_weight","excess_contribution","status"]
                           if directory == "spatial" else ["dataset", "label", "threshold", "sensitivity", "n_eligible", "observed_multi", "expected_multi", "oe", "status"]
                           if directory == "archaic" else list(data[0])[:10])
                sections.append(table_html(data, columns, 30 if directory == "spatial" else 20))
            else:
                sections.append("<p class='muted'>The output contains no estimable rows; this is not a negative biological result.</p>")
            if figure and recorded_output(base, figure, result):
                sections.append(f"<img src='{html.escape(str((base/figure).relative_to(out)))}' alt='{html.escape(title)}'>")
            if directory == "selection":
                for name in ("published_selection_alleles.tsv.gz", "availability.tsv", "manifest.json"):
                    if recorded_output(base, name, result) or name == manifest_name and (base/name).is_file():
                        sections.append(f"<p><a href='{str((base/name).relative_to(out))}'>{name}</a></p>")
            if directory == "ancient":
                sections.append('<h3>Same-person HERV matched-panel statistics</h3><p>Incomplete personal means are observed-pair conditional. Formal personal models require the complete frozen panel; period comparisons retain fixed pair and composition weights. Person and chromosome uncertainty have separate scopes. These quantities are not ancient gnomAD O/E.</p>')
                for name in ("allele_time_heatmap.png", "allele_panel_availability.png", "F6_joint_measurability.png"):
                    if recorded_output(base, name, result):
                        sections.append(f"<img src='{str((base/name).relative_to(out))}' alt='Full ancient allele trajectories and representation'>")
                for name in ("allele_trajectories.tsv.gz", "measurement_plans.json", "adjusted_allele_trends.tsv.gz", "all_allele_audit.tsv.gz", "sample_time_region_coverage.tsv", "catalog_coverage.tsv", "matched_locus_role_trajectories.tsv", "model_covariates.json", "joint_measurability_funnel.tsv", "joint_person_statistics.tsv.gz", "joint_fixed_pair_panel.tsv", "joint_fixed_panel_trajectories.tsv", "joint_standardized_trajectories.tsv", "joint_person_models.tsv", "joint_missingness_models.tsv", "standardized_family_trajectories.tsv", "regional_fixed_panel_family_trajectories.tsv", "joint_manifest.json"):
                    if recorded_output(base, name, result):
                        sections.append(f"<p><a href='{str((base/name).relative_to(out))}'>{name}</a></p>")
            if directory == "spatial":
                if recorded_output(base, "primary_whole_element_profiles.png", result):
                    sections.append(f"<img src='{str((base/'primary_whole_element_profiles.png').relative_to(out))}' alt='Primary expected-count definition along aligned repeat bases'>")
                for name in ("primary_spatial_identity.tsv", "primary_counter_reconciliation.tsv", "family_counter_reconciliation.tsv", "spatial_annotation_audit.tsv"):
                    if recorded_output(base, name, result):
                        sections.append(f"<p><a href='{str((base/name).relative_to(out))}'>{name}</a></p>")
                sections.append("<p>The common-support sensitivity applies the recorded training minimum and removes unsupported observed and expected counts together.</p>")
                for name in ("spatial_partitions.tsv", "whole_element_profiles.tsv", "spatial_identity.tsv", "spatial_contrasts.tsv"):
                    if recorded_output(base, name, result):
                        sections.append(f"<p><a href='{str((base/name).relative_to(out))}'>{name}</a></p>")
        else:
            sections.append("<p class='muted'>Not run / inputs unavailable. No effect estimate is inferred.</p>")
    if pca_temporal_dir:
        sections.append('<h2>Genome-wide PCA temporal structure (separate estimand)</h2>')
        try:
            temporal = validate_manifest(pca_temporal_dir, require_evidence=True)
            sections.append('<p>Fixed shared coordinates; raw, fixed-composition and adjusted PC responses. This is not HERV selection evidence.</p>')
            sections.append('<pre>' + html.escape(json.dumps({k: temporal.get(k) for k in ('status', 'sample_counts', 'coordinate_contract', 'limitations')}, ensure_ascii=False, indent=2)) + '</pre>')
            for name in ('temporal_models.tsv', 'standardized_trajectories.tsv', 'report.html'):
                if recorded_output(Path(pca_temporal_dir), name, temporal):
                    link = os.path.relpath(Path(pca_temporal_dir) / name, out)
                    sections.append(f'<p><a href="{html.escape(link)}">PCA {name}</a></p>')
        except (OSError, ValueError) as exc:
            sections.append('<p>PCA temporal results unavailable: ' + html.escape(str(exc)) + '</p>')
    cohort = root/"final/herv_cohort/herv_cohort_burden.tsv"
    if not cohort.is_file() and os.environ.get("ALLELE_COHORT_TABLE"):
        cohort = Path(os.environ["ALLELE_COHORT_TABLE"])
    if cohort.is_file():
        n = len({r["sample_id"] for r in rows(cohort)})
        sections.append(f"<h2>Measured cohort context</h2><p>The existing cohort table contains {n:,} people. Its measure is nonreference allele-copy burden within shared reference HERV-related intervals, including indels. It is not a count of HERV insertions or personal multi-allelic loci.</p>")
    sections.append("<h2>Interpretation</h2><p>Regional association and descriptive temporal trends do not establish mutation origin or selection. Independent genotype support, repeat-aware mapping, and an appropriate population design are required for those questions.</p>")
    css = "body{font:16px/1.6 system-ui,sans-serif;color:#20333c;max-width:1160px;margin:40px auto;padding:0 24px;background:#fbfcfc}h1{font-size:34px}h2{margin-top:36px;color:#176b87}.lead{font-size:20px}.muted{color:#687a83}img{max-width:100%;display:block;margin:24px auto;background:white}.table{overflow:auto}table{border-collapse:collapse;font-size:13px;width:100%;background:white}th,td{padding:9px;border-bottom:1px solid #dae4e8;text-align:left;white-space:nowrap}th{background:#eaf2f5}a{color:#146d97}code{font-size:13px}"
    text = "<!doctype html><html lang='en'><meta charset='utf-8'><title>HERV integrated evidence report</title><style>"+css+"</style><body>"+"\n".join(sections)+"</body></html>"
    temporary = out / "report.html.partial"
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, out / "report.html")
    return dict(status="complete", report=str(out/"report.html"), primary_available=primary_available)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=["report"])
    p.add_argument("--out", required=True)
    p.add_argument("--modern-dir", required=True)
    p.add_argument("--analysis-root", required=True)
    p.add_argument("--pca-temporal-dir")
    a = p.parse_args(argv)
    result = make_report(a.out, a.modern_dir, a.analysis_root, a.pca_temporal_dir)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
