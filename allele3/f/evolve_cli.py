#!/usr/bin/env python3
"""Ancient-DNA, fixed PCA/time and optional sparse Atlas workflow."""
from __future__ import annotations

import argparse
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from ancient_join import (index_target_variants, join_selection, select_candidates,
                          table_rows, _stage_identity, completed_stage, finish_stage)
from core import VERSION, dump_json, fingerprint

PREFIX = "AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes"
DEFAULT_ANCIENT = "/mnt/f/gen/ancient"


def parser():
    p = argparse.ArgumentParser(description="allele3 evolution: local ancient trajectories, PCA × time, bounded optional Atlas",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("command", choices=["evolve", "prepare-ancient", "join-selection", "select-candidates",
                   "panel-index", "trajectories", "atlas-plan", "atlas-query", "atlas-summarize",
                   "compare-functions", "pca-plan", "pca-time", "report-evolution", "audit-ancient", "audit-counts"])
    p.add_argument("--out", "--dir-out", default="/mnt/d/analysis/allele3/evolution")
    p.add_argument("--ancient-dir", default=DEFAULT_ANCIENT)
    p.add_argument("--modern-out", default="/mnt/d/analysis/allele3/gnomad_local")
    p.add_argument("--pgen-prefix", help="Default: ANCIENT_DIR/AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes")
    p.add_argument("--selection")
    p.add_argument("--join-strategy", choices=["selection-first", "modern-first"], default="selection-first",
                   help="Selection-first indexes the complete published SNV coordinate universe, then validates actual modern intersections")
    p.add_argument("--counts-tsv", help="Optional legacy check counts.tsv for descriptive audit; no replacement for per-site VCF data")
    p.add_argument("--no-plots", action="store_true", help="Skip the optional counts-audit plot")
    p.add_argument("--hq-variants")
    p.add_argument("--supplement-xlsx")
    p.add_argument("--qpadm", help="Optional one verified 3way or 4way TSV; no basis mixing")
    p.add_argument("--ancestry-model", choices=["3way", "4way"])
    p.add_argument("--metadata", help="Canonical metadata from prepare-ancient; otherwise prepare automatically for evolve")
    p.add_argument("--fasta38", default=os.environ.get("ALLELE_FASTA", "/mnt/f/gen/1hgp/38/GRCH38.fasta.gz"))
    p.add_argument("--fasta37", help="GRCh37 / hs37d5 FASTA with existing .fai; required for join")
    p.add_argument("--chain38to37", help="Local hg38ToHg19.over.chain.gz")
    p.add_argument("--chain37to38", help="Local hg19ToHg38.over.chain.gz")
    p.add_argument("--intersection", help="Existing join/selection_intersection.tsv.gz")
    p.add_argument("--candidates", help="Candidate TSV; default OUT/candidates/candidates.tsv")
    p.add_argument("--selected-variants", help="Indexed targets TSV; default OUT/panel/selected_variants.tsv")
    p.add_argument("--max-sites", type=int, default=2000, help="Total strong+weak candidate sites; Atlas requests <= 3 times this")
    p.add_argument("--max-variants", type=int, default=9000, help="Hard budget for UNIQUE Atlas variants, not retry attempts")
    p.add_argument("--strong-posterior", type=float, default=.99)
    p.add_argument("--weak-posterior", type=float, default=.5)
    p.add_argument("--min-extra-ac", type=int, default=5)
    p.add_argument("--block-bp", type=int, default=1000000)
    p.add_argument("--include-mhc", action="store_true", help="Sensitivity analysis: include chr6:25–34 Mb GRCh37")
    p.add_argument("--allow-selection-warnings", action="store_true")
    p.add_argument("--seed", type=int, default=20261021)
    p.add_argument("--bins", default="0,2000,4000,6000,8000,10000")
    p.add_argument("--min-samples", type=int, default=10)
    p.add_argument("--weights", help="Fixed stratum/weight TSV; unsupported bins stay NA")
    p.add_argument("--count-kind", choices=["hardcall", "dosage"], default="hardcall")
    p.add_argument("--skip-trajectories", action="store_true", help="Plan/join without opening PGEN; leaves explicit missing trajectory stage")
    p.add_argument("--with-atlas", action="store_true", help="Explicitly execute the bounded Atlas plan after local analysis")
    p.add_argument("--atlas-plan", help="Existing plan.json; default OUT/atlas_plan/plan.json")
    p.add_argument("--atlas-cache", help="Existing run directory or cache directory for geometry")
    p.add_argument("--scorers", help="Comma-separated scorer enum names; default six signed accessibility/TF/transcription scorers")
    p.add_argument("--ontology-terms", help="Optional comma-separated ontology IDs for a smaller hypothesis-driven panel")
    p.add_argument("--api-key-file", help="Secret file; contents never written into manifests")
    p.add_argument("--api-key-env", default="ALPHAGENOME_API_KEY")
    p.add_argument("--workers", type=int, default=2)
    p.add_argument("--requests-per-second", type=float, default=2.)
    p.add_argument("--retries", type=int, default=3)
    p.add_argument("--timeout", type=float, default=60.)
    p.add_argument("--geometry", help="Existing geometry TSV for compare-functions")
    p.add_argument("--with-pca", action="store_true", help="Execute generated local PLINK2 PCA/projection script and render 3D view")
    p.add_argument("--eigenvec", help="Existing eigenvec or common-basis .sscore for all plotted samples")
    p.add_argument("--pca-basis", choices=["joint_balanced", "modern_projection"], default="joint_balanced")
    p.add_argument("--pca-basis-label", help="Provenance label for an externally provided shared PCA coordinate table")
    p.add_argument("--pca-maf", type=float, default=.05)
    p.add_argument("--pca-geno", type=float, default=.1)
    p.add_argument("--pca-max-per-group", type=int, default=200)
    p.add_argument("--pca-components", type=int, default=10)
    p.add_argument("--pca-target-flank-bp", type=int, default=1000000)
    p.add_argument("--dir-pgen-cache", help="Optional verified execution copies of PGEN/PVAR/PSAM on a fast local filesystem; outputs remain under --out")
    p.add_argument("--pca-threads", type=int, default=4)
    p.add_argument("--pca-memory-mb", type=int, default=8192)
    p.add_argument("--pca-exclude-bed", help="Optional GRCh37 BED for long-range LD or other background exclusions")
    p.add_argument("--pca-max-variants", type=int, default=12, help="Maximum candidate genotypes embedded in 3D HTML")
    p.add_argument("--plink2", default="plink2")
    p.add_argument("--html", help="Output HTML for pca-time; default OUT/pca/pca_time.html")
    p.add_argument("--simulated", action="store_true", help="Prominent synthetic-data label for examples only")
    p.add_argument("--with-universes", action="store_true")
    p.add_argument("--u1-sites", type=int, default=20000)
    p.add_argument("--u3-sites", type=int, default=3000)
    p.add_argument("--controls-per-case", type=int, default=4)
    p.add_argument("--ld-block-dir", default="/mnt/f/LD/block")
    p.add_argument("--with-ukb", action="store_true")
    p.add_argument("--ukb-mfi-dir", default="/mnt/f/gen/ukb/37/imp_mfi")
    p.add_argument("--metag-census", help="Completed metag genome census JSON, reused read-only")
    p.add_argument("--version", action="version", version=VERSION)
    return p


def paths(a):
    ancient, out = Path(a.ancient_dir).expanduser().resolve(), Path(a.out).expanduser().resolve()
    prefix = Path(a.pgen_prefix).expanduser().resolve() if a.pgen_prefix else ancient / PREFIX
    return dict(out=out, prefix=prefix, pgen=Path(str(prefix) + ".pgen"), psam=Path(str(prefix) + ".psam"),
                pvar=Path(str(prefix) + ".pvar"),
                selection=Path(a.selection) if a.selection else ancient / "Selection_Summary_Statistics_01OCT2025.tsv.gz",
                hq=Path(a.hq_variants) if a.hq_variants else ancient / "HighQualityVariants.tsv",
                supplement=Path(a.supplement_xlsx) if a.supplement_xlsx else ancient / "41586_2026_10358_MOESM3_ESM.xlsx",
                metadata=Path(a.metadata) if a.metadata else out / "metadata" / "sample_metadata.tsv",
                candidates=Path(a.candidates) if a.candidates else out / ("universes" if a.with_universes else "candidates") / "candidates.tsv",
                selected=Path(a.selected_variants) if a.selected_variants else out / "panel" / "selected_variants.tsv",
                intersection=Path(a.intersection) if a.intersection else out / "join" / "selection_intersection.tsv.gz",
                atlas_plan=Path(a.atlas_plan) if a.atlas_plan else out / "atlas_plan" / "plan.json")


def _require(value, name):
    if not value:
        raise ValueError(f"{name} is required for reciprocal build harmonization")
    return value


def _prepare(a, d):
    from ancient_metadata import prepare_metadata
    out = d["out"] / "metadata"
    inputs = dict(psam=d["psam"], supplement=d["supplement"])
    if a.qpadm:
        inputs["qpadm"] = a.qpadm
    identity = _stage_identity("metadata", inputs, dict(ancestry_model=a.ancestry_model))
    if old := completed_stage(out, identity):
        return old
    result = prepare_metadata(d["psam"], d["supplement"], out, qpadm_path=a.qpadm,
                              ancestry_model=a.ancestry_model)
    return finish_stage(out, identity, result, ["sample_metadata.tsv", "metadata_qc.tsv", "metadata_qc.json"])


def _join(a, d):
    return join_selection(a.modern_out, d["selection"], d["out"] / "join", a.fasta38,
                          _require(a.fasta37, "--fasta37"), _require(a.chain38to37, "--chain38to37"),
                          _require(a.chain37to38, "--chain37to38"), strategy=a.join_strategy)


def _counts(a, d):
    from counts_audit import audit_counts
    if not a.counts_tsv:
        raise ValueError("audit-counts requires --counts-tsv")
    return audit_counts(a.counts_tsv, d["out"] / "counts_audit", make_plots=not a.no_plots)


def _select(a, d):
    if a.with_universes:
        from universes import build_universes
        return build_universes(d['intersection'],d['out']/'universes',ld_block_dir=a.ld_block_dir,
            u1_sites=a.u1_sites,strong_posterior=a.strong_posterior,weak_posterior=a.weak_posterior,
            min_extra_ac=a.min_extra_ac,controls_per_case=a.controls_per_case,seed=a.seed,
            exclude_mhc=not a.include_mhc,allow_filter_warnings=a.allow_selection_warnings)
    return select_candidates(d["intersection"], d["out"] / "candidates", max_sites=a.max_sites,
                             strong_posterior=a.strong_posterior, weak_posterior=a.weak_posterior,
                             min_extra_ac=a.min_extra_ac, block_bp=a.block_bp,
                             exclude_mhc=not a.include_mhc, allow_filter_warnings=a.allow_selection_warnings, seed=a.seed)


def _panel(a, d):
    return index_target_variants(d["candidates"], d["pvar"], d["hq"], d["out"] / "panel")


def _trajectories(a, d):
    from trajectory import aggregate_trajectories, iter_pgen_variants
    out = d["out"] / "trajectories"
    inputs = {key: d[key] for key in ("pgen", "psam", "selected", "metadata")}
    if a.weights:
        inputs["weights"] = a.weights
    identity = _stage_identity("trajectories", inputs,
                               dict(bins=a.bins, min_samples=a.min_samples, count_kind=a.count_kind))
    if old := completed_stage(out, identity):
        return old
    samples, variants = iter_pgen_variants(d["pgen"], d["psam"], d["selected"], count_kind=a.count_kind)
    result = aggregate_trajectories(variants, samples, d["metadata"], out,
                                    bins=tuple(float(x) for x in a.bins.split(",")), min_samples=a.min_samples,
                                    weights_path=a.weights)
    return finish_stage(out, identity, result, ["trajectory_bins.tsv", "trajectory_standardized.tsv",
                        "sample_bin_composition.tsv", "sample_inclusion.tsv", "allele_availability.tsv", "trajectory_manifest.json"])


def _atlas_plan(a, d):
    from sparse_atlas import build_plan
    source=d['out']/'universes'/'U3.tsv' if a.with_universes else d['candidates']
    return build_plan(source, d["out"] / "atlas_plan", max_sites=a.u3_sites if a.with_universes else a.max_sites,
                      max_variants=a.max_variants, scorers=a.scorers,
                      ontology_terms=a.ontology_terms)


def _atlas_query(a, d):
    from sparse_atlas import run_sparse
    return run_sparse(d["atlas_plan"], d["out"] / "atlas", api_key_file=a.api_key_file,
                      api_key_env=a.api_key_env, max_variants=a.max_variants, workers=a.workers,
                      retries=a.retries, timeout=a.timeout, requests_per_second=a.requests_per_second)


def _atlas_summarize(a, d):
    from sparse_atlas import summarize_geometry
    return summarize_geometry(d["atlas_plan"], a.atlas_cache or d["out"] / "atlas", d["out"] / "geometry")


def _pca_plan(a, d):
    from pca_time import prepare_pca_plan
    prefix = d["prefix"]
    if a.dir_pgen_cache:
        from input_cache import cache_pgen_prefix
        prefix, cache_record = cache_pgen_prefix(prefix, a.dir_pgen_cache)
        dump_json(d["out"] / "pca_input_cache.json", cache_record)
    pca_targets=d['selected']
    if a.with_universes:
        # U1 covers the genome: expanding 20k sites by 1Mb would erase all PCA markers.
        # Exclude U2 strong loci and flanks; U1 sites are only sampling/background.
        from core import tsv_writer
        pca_targets=d['out']/'pca_exclusion_targets.tsv'
        target_rows={(r['chrom37'],int(r['pos37'])) for r in table_rows(d['out']/'universes'/'U2_all_strong.tsv')}
        content='chrom\tpos\tbuild\n'+''.join(f'{c}\t{p}\tGRCh37\n' for c,p in sorted(target_rows))
        if pca_targets.exists():
            if pca_targets.read_text()!=content: raise ValueError('PCA strong-target exclusions changed')
        else: pca_targets.write_text(content)
    return prepare_pca_plan(prefix, d["metadata"], d["hq"], d["out"] / "pca", basis=a.pca_basis,
                            maf=a.pca_maf, geno=a.pca_geno, max_per_group=a.pca_max_per_group,
                            seed=a.seed, exclude_bed=a.pca_exclude_bed, target_variants=pca_targets, plink2=a.plink2,
                            n_components=a.pca_components, target_flank_bp=a.pca_target_flank_bp,
                            threads=a.pca_threads, memory_mb=a.pca_memory_mb)


def _pca_time(a, d, eigenvec=None):
    from pca_time import render_pca_time
    source = eigenvec or a.eigenvec
    if not source:
        raise ValueError("pca-time needs --eigenvec from one fixed PCA basis for all epochs")
    label = a.pca_basis if eigenvec else (a.pca_basis_label or "用户提供的共享 PCA 基底；请核对其训练样本与标记")
    kwargs = dict(basis_label=label, max_variants=a.pca_max_variants, simulated=a.simulated, count_kind=a.count_kind, candidate_path=d.get("candidates"))
    if d["selected"].exists() and d["pgen"].exists() and d["psam"].exists():
        from core import tsv_writer
        rows = list(table_rows(d["selected"]))
        by_id = {r["variant_id"]: r for r in rows}
        links = d["out"] / "panel" / "site_variant_links.tsv"
        ordered = []
        if links.exists() and not a.selected_variants:
            for link in table_rows(links):
                r = by_id.get(link["variant_id"])
                if r and r not in ordered:
                    ordered.append(r)
                if len(ordered) >= a.pca_max_variants:
                    break
        else:
            ordered = rows[:a.pca_max_variants]
        if ordered:
            bounded = d["out"] / "pca" / "overlay_targets.tsv"
            h, w = tsv_writer(bounded, list(ordered[0]))
            with h:
                w.writerows(ordered)
            dump_json(d["out"] / "pca" / "overlay_selection.json",
                      dict(rule="first candidate sites in pre-functional selection order, anchor then extra; bounded visualization only",
                           total_available=len(rows), displayed=len(ordered), max_variants=a.pca_max_variants,
                           source=fingerprint(d["selected"])))
            kwargs.update(pgen_path=d["pgen"], psam_path=d["psam"], selected_variants_path=bounded)
    return render_pca_time(source, d["metadata"], a.html or d["out"] / "pca" / "pca_time.html", **kwargs)


def run_evolution(a, d):
    from evolve_report import compare_geometry, write_report
    started = time.monotonic()
    result = dict(version=VERSION, status="running", stages={}, api_requested=bool(a.with_atlas),
                  pca_execution_requested=bool(a.with_pca), modern_out=str(Path(a.modern_out).resolve()))
    path = d["out"] / "evolution.json"

    def stage(name, fn):
        print(f"[evolution] {name}", flush=True)
        value = fn()
        result["stages"][name] = value
        dump_json(path, result)
        return value

    try:
        if a.counts_tsv:
            stage("counts_audit", lambda: _counts(a, d))
        if a.metadata:
            result["stages"]["metadata"] = dict(source=fingerprint(d["metadata"]), adapter="user_supplied_canonical")
        else:
            stage("metadata", lambda: _prepare(a, d))
        stage("join", lambda: _join(a, d))
        stage("candidates", lambda: _select(a, d))
        if a.with_ukb:
            if not a.with_universes: raise ValueError('--with-ukb requires --with-universes so the full U0 is scanned first')
            from ukb_mfi import scan_mfi_dir
            stage('ukb_mfi',lambda: scan_mfi_dir(a.ukb_mfi_dir,d['out']/'universes'/'U0.tsv.gz',d['out']/'ukb'/'U0_mfi.tsv'))
        stage("panel", lambda: _panel(a, d))
        if not a.skip_trajectories:
            stage("trajectories", lambda: _trajectories(a, d))
        if not a.with_universes: stage("atlas_plan", lambda: _atlas_plan(a, d))
        pca = stage("pca_plan", lambda: _pca_plan(a, d))
        if a.with_pca:
            # Executable paths/arguments are written with shlex.quote by planner.
            # No shell interpolation of secrets or uncontrolled command strings.
            script = pca["script_path"]
            stage("pca_execution", lambda: _execute_pca(script))
            stage("pca_visualization", lambda: _pca_time(a, d, pca["projected_path"]))
        elif a.eigenvec:
            stage("pca_visualization", lambda: _pca_time(a, d))
        if a.with_universes:
            from final_figures import write_local_review
            stage('local_review',lambda: write_local_review(a.modern_out,d['out'],a.metag_census,make_plots=not a.no_plots))
            from universes import select_u3
            stage('U3_selection',lambda: select_u3(d['out']/'universes',d['out']/'universes'/'U3.tsv',a.u3_sites,a.seed))
            stage('atlas_plan',lambda: _atlas_plan(a,d))
        if a.with_atlas:
            stage("atlas_query", lambda: _atlas_query(a, d))
            geometry = stage("atlas_geometry", lambda: _atlas_summarize(a, d))
            stage("functional_comparison", lambda: compare_geometry(geometry["outputs"]["geometry"], d["out"] / "comparison", seed=a.seed))
            if a.with_universes:
                stage('final_figures',lambda: write_local_review(a.modern_out,d['out'],a.metag_census,make_plots=not a.no_plots))
        result.update(status="complete", elapsed_seconds=time.monotonic() - started)
    except BaseException as exc:
        result.update(status="incomplete", error_type=type(exc).__name__, elapsed_seconds=time.monotonic() - started)
        # Do not persist potentially credential-bearing exception text.
        dump_json(path, result)
        write_report(d["out"])
        raise
    result["report"] = write_report(d["out"])
    dump_json(path, result)
    return result


def _execute_pca(script):
    process = subprocess.run(["bash", str(script)], check=False)
    if process.returncode:
        raise RuntimeError(f"PLINK2 PCA script failed with exit {process.returncode}; inspect its local logs")
    return dict(status="complete", script=str(script))


def audit(a, d):
    result = dict(files={})
    for key in ("psam", "pvar", "pgen", "selection", "hq", "supplement"):
        p = d[key]
        result["files"][key] = fingerprint(p) if p.exists() else dict(path=str(p), missing=True)
    if d["psam"].exists():
        result["psam_data_rows"] = sum(1 for _ in table_rows(d["psam"], ("IID",)))
    for key in ("pvar", "selection", "hq"):
        if d[key].exists():
            rows = table_rows(d[key])
            first = next(rows, None)
            rows.close()
            result[key + "_first_record"] = first
    result["note"] = "PVAR full record count is obtained by panel-index; wc -l includes headers. PGEN fields are not inferred from filename."
    dump_json(d["out"] / "ancient_audit.json", result)
    return result


def main(argv=None):
    a = parser().parse_args(argv)
    d = paths(a)
    d["out"].mkdir(parents=True, exist_ok=True)
    try:
        with (d["out"] / ".evolution.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            actions = {"evolve": run_evolution, "prepare-ancient": _prepare, "join-selection": _join,
                       "select-candidates": _select, "panel-index": _panel, "trajectories": _trajectories,
                       "atlas-plan": _atlas_plan, "atlas-query": _atlas_query, "atlas-summarize": _atlas_summarize,
                       "pca-plan": _pca_plan, "pca-time": _pca_time, "audit-ancient": audit,
                       "audit-counts": _counts}
            if a.command == "report-evolution":
                from evolve_report import write_report
                result = dict(report=write_report(d["out"]))
            elif a.command == "compare-functions":
                from evolve_report import compare_geometry
                if not a.geometry:
                    raise ValueError("compare-functions requires --geometry")
                result = compare_geometry(a.geometry, d["out"] / "comparison", seed=a.seed)
            else:
                result = actions[a.command](a, d)
        # Large manifests stay in files; avoid printing secrets or 6,000 variants.
        print(json.dumps({k: result[k] for k in ("status", "manifest_path", "metadata_path", "candidate_path",
                                                 "plan_path", "report", "report_path", "qc", "pvar_records") if k in result},
                         ensure_ascii=False, indent=2))
        return 0
    except (ValueError, OSError, RuntimeError) as exc:
        print(f"allele3 evolution: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
