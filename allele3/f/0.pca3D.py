#!/usr/bin/env python3
from __future__ import annotations

# The shared utility filename is sorted first; use a stable import name in Python.
import importlib.util as _common_import
import sys as _common_sys
from pathlib import Path as _CommonPath
if "allele_common" not in _common_sys.modules:
    _common_spec = _common_import.spec_from_file_location(
        "allele_common", _CommonPath(__file__).with_name("0.common.py"))
    _common_module = _common_import.module_from_spec(_common_spec)
    _common_sys.modules[_common_spec.name] = _common_module
    try:
        _common_spec.loader.exec_module(_common_module)
    except BaseException:
        _common_sys.modules.pop(_common_spec.name, None)
        raise

# Standalone PCA planning, rendering, dosage export and carrier centroids.





# ========================================================================
# variant centroid
# ========================================================================

# Dosage export and carrier centroids in one shared PCA coordinate system.

import argparse, csv, gzip, json, math, sys
from collections import defaultdict, Counter
from pathlib import Path
from evolve import _coordinates
from evolve import iter_pgen_variants


def op(path):return gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)


def export_dosages(pgen,psam,selected,out,count_kind='hardcall'):
    """Stream original PVAR-indexed alleles, including explicit unavailable rows."""
    ids,variants=iter_pgen_variants(pgen,psam,selected,count_kind=count_kind)
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True);n=0
    with gzip.open(out,'wt',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['sample_id','variant_id','dosage','status','count_kind'])
        try:
            for r in variants:
                if r['dosages'] is None:
                    w.writerow(['.',r['variant_id'],'NA',r['status'],count_kind]);continue
                for iid,d in zip(ids,r['dosages']):
                    w.writerow([iid,r['variant_id'],'NA' if d is None else d,'missing_call' if d is None else 'ok',count_kind]);n+=1
        finally:variants.close()
    return dict(output=str(out),rows=n,count_kind=count_kind)


def run_centroids(a):
    if a.bin_years<=0:raise ValueError('Positive --bin-years required')
    coords,coordinate_kind=_coordinates(a.coordinates);meta={}
    with open(a.metadata) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            iid=r.get('sample_id') or r.get('IID')
            if iid in meta:raise ValueError(f'Duplicate metadata sample: {iid}')
            if not iid or r.get('include','1')!='1':continue
            try:age=float(r.get('age_bp') or r.get('age'))
            except (TypeError,ValueError):continue
            if not math.isfinite(age) or age<0:continue
            meta[iid]=(r.get('cohort') or ('modern' if age==0 else 'ancient'),age)
    agg=defaultdict(lambda:[0,0.,0.,0.]);unavailable=defaultdict(set);measured=set();missing=Counter()
    current=None;seen_samples=set();finished=set()
    with op(a.dosages) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            iid=r.get('sample_id') or r.get('IID');vid=r.get('variant_id') or r.get('site_id')
            if not vid:raise ValueError('Missing dosage variant ID')
            if current!=vid:
                if vid in finished:raise ValueError('Dosage table must group rows by variant_id')
                if current is not None:finished.add(current)
                current=vid;seen_samples=set()
            if iid in seen_samples:raise ValueError('Duplicate dosage sample/variant')
            seen_samples.add(iid)
            if r.get('status','ok') not in {'ok','assayed_pass','measured'}:
                unavailable[vid].add(r.get('status','unavailable'));continue
            if iid not in coords or iid not in meta:missing[vid]+=1;continue
            try:d=float(r['dosage'])
            except (KeyError,ValueError):raise ValueError('Measured dosage must be numeric')
            if not math.isfinite(d) or not 0<=d<=2:raise ValueError('Invalid diploid dosage')
            cohort,age=meta[iid];k=-1 if cohort=='modern' else int(age//a.bin_years)
            p1,p2=coords[iid];z=agg[vid,cohort,k];z[0]+=1;z[1]+=d;z[2]+=d*p1;z[3]+=d*p2;measured.add(vid)
    rows=[]
    for (vid,cohort,k),(n,ac,p1,p2) in sorted(agg.items()):
        rows.append(dict(variant_id=vid,cohort=cohort,start_bp=0 if k<0 else k*a.bin_years,end_bp=0 if k<0 else (k+1)*a.bin_years,n=n,allele_count=ac,af=ac/(2*n),pc1=p1/ac if ac>0 else 'NA',pc2=p2/ac if ac>0 else 'NA',status='ok' if ac>0 else 'observed_zero'))
    for vid in sorted((set(unavailable)|set(missing))-measured):
        rows.append(dict(variant_id=vid,cohort='NA',start_bp='NA',end_bp='NA',n=0,allele_count='NA',af='NA',pc1='NA',pc2='NA',status='not_projectable:'+','.join(sorted(unavailable[vid] or {'missing_coordinates_or_age'}))))
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    fields=['variant_id','cohort','start_bp','end_bp','n','allele_count','af','pc1','pc2','status']
    with open(out/'allele_centroids.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
    valid=[r for r in rows if r['status']=='ok']
    if valid:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig=plt.figure(figsize=(9,7));ax=fig.add_subplot(projection='3d')
        for vid in sorted({r['variant_id'] for r in valid})[:12]:
            selected=sorted([r for r in valid if r['variant_id']==vid],key=lambda r:r['start_bp'])
            ax.plot([r['pc1'] for r in selected],[r['pc2'] for r in selected],[(r['start_bp']+r['end_bp'])/2 for r in selected],marker='o',label=vid)
        ax.set(xlabel='PC1 carrier centroid',ylabel='PC2 carrier centroid',zlabel='Age BP');ax.legend(fontsize=6);fig.tight_layout()
        fig.savefig(out/'allele_centroid_trajectories.png',dpi=180)
        plt.close(fig)
    result=dict(status='complete',rows=len(rows),unavailable_variants=len((set(unavailable)|set(missing))-measured),coordinate_kind=coordinate_kind,bin_years=a.bin_years,rule='No measured individual dosage -> no PCA projection. Modern and ancient cohorts remain separate. Missing ages remain unavailable.')
    (out/'manifest.json').write_text(json.dumps(result,indent=2)+'\n');return result


def main_centroid(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--coordinates',required=True);p.add_argument('--metadata',required=True)
    p.add_argument('--dosages');p.add_argument('--pgen-prefix');p.add_argument('--selected-variants');p.add_argument('--count-kind',choices=['hardcall','dosage'],default='hardcall');p.add_argument('--out',required=True);p.add_argument('--bin-years',type=int,default=2000)
    a=p.parse_args(argv)
    if not a.dosages:
        if not a.pgen_prefix or not a.selected_variants:p.error('Supply --dosages or both --pgen-prefix and --selected-variants')
        a.dosages=str(Path(a.out)/'selected_variant_dosages.tsv.gz');export_dosages(a.pgen_prefix+'.pgen',a.pgen_prefix+'.psam',a.selected_variants,a.dosages,a.count_kind)
    print(json.dumps(run_centroids(a),indent=2))


# ========================================================================
# pca3d cli
# ========================================================================

# Independent GRCh37 PCA workflow and offline PC1/PC2/age visualization.
# 
# Reuses calculation libraries without importing the allele3 workflow dispatcher.


import argparse
from contextlib import nullcontext
import csv
import fcntl
import json
from pathlib import Path
import subprocess
import sys

from evolve import _atomic_text, _fingerprint, _pca_rows as _rows, _evolution_sha as _sha, prepare_pca_plan, render_pca_time

PREFIX = "AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes"


def parser():
    p = argparse.ArgumentParser(description=__doc__, allow_abbrev=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("command", choices=("run", "plan", "render", "atlas", "prepare-basis", "project-plan", "project"), nargs="?", default="run")
    p.add_argument("--out", "--dir-out", type=Path, default=Path("/mnt/f/gen/ancient/allele3/pca3d"))
    p.add_argument("--replace", choices=("TRUE", "FALSE"), default="TRUE", help="For run: recalculate from scratch, or resume verified results")
    p.add_argument("--ancient-dir", "--dir-ancient", type=Path, default=Path("/mnt/f/gen/ancient"))
    p.add_argument("--pgen-prefix", type=Path)
    p.add_argument("--metadata", type=Path, help="Canonical TSV; required for custom panels / render")
    p.add_argument("--hq-variants", type=Path, help="GRCh37 background variant TSV")
    p.add_argument("--supplement-xlsx", type=Path, help="Original ancient dataset metadata adapter")
    p.add_argument("--qpadm", type=Path)
    p.add_argument("--ancestry-model", choices=("3way", "4way"))
    p.add_argument("--dir-pgen-cache", type=Path)
    p.add_argument("--eigenvec", type=Path, help="Existing common-basis .sscore / .eigenvec for render")
    p.add_argument("--html", type=Path)
    p.add_argument("--source-html", type=Path, help="Existing complete allele3-data HTML for lossless atlas re-rendering")
    p.add_argument("--selected-variants", type=Path, help="Optional indexed PVAR targets for allele overlays")
    p.add_argument("--target-variants", type=Path, help="Optional separate GRCh37 background exclusions")
    p.add_argument("--pca-basis", choices=("joint_balanced", "modern_projection", "ancient_balanced"), default="joint_balanced")
    p.add_argument("--pca-basis-label")
    p.add_argument("--pca-basis-min-age-bp", type=float)
    p.add_argument("--pca-basis-max-age-bp", type=float, help="Ancient training limit; ancient_balanced defaults to 10,000 BP")
    p.add_argument("--pca-maf", type=float, default=.05)
    p.add_argument("--pca-geno", type=float, default=.1)
    p.add_argument("--pca-components", type=int, default=10)
    p.add_argument("--pca-max-per-group", type=int, default=200)
    p.add_argument("--pca-time-bin-years", type=int, default=2000)
    p.add_argument("--pca-layer-years", type=int, default=1000, help="Display-only layer width; independent of PCA training time bins")
    p.add_argument("--pca-time-view", choices=("layered", "continuous"), default="layered")
    p.add_argument("--pca-color-by", choices=("region", "layer", "age", "cohort", "source"), default="region")
    p.add_argument("--pca-time-scale", choices=("linear", "log1p"), default="linear", help="Atlas display only; never changes original dates or genotype PCA")
    p.add_argument("--pca-display-max-age-bp", type=float, default=10000, help="Main atlas age limit; older samples remain in the full-range supplement")
    p.add_argument("--skip-atlas", action="store_true", help="Skip additional publication atlas after run/render")
    p.add_argument("--pca-target-flank-bp", type=int, default=1000000)
    p.add_argument("--pca-exclude-bed", type=Path)
    p.add_argument("--pca-threads", type=int, default=4)
    p.add_argument("--pca-memory-mb", type=int, default=8192)
    p.add_argument("--seed", type=int, default=20261021)
    p.add_argument("--plink2", default="plink2")
    p.add_argument("--pca-max-variants", type=int, choices=range(1, 13), default=12)
    p.add_argument("--count-kind", choices=("hardcall", "dosage"), default="hardcall")
    p.add_argument("--title")
    p.add_argument("--simulated", action="store_true")
    p.add_argument("--basis-dir", type=Path, help="Completed frozen PCA directory for external projection")
    p.add_argument("--reuse-basis-dir", type=Path, help="Verify a completed basis before rendering its existing coordinates")
    p.add_argument("--target-glob", default="/mnt/f/gen/ukb/37/imp/chr[0-9]*.pgen")
    p.add_argument("--project-out", type=Path)
    p.add_argument("--reference-metadata", type=Path, help="External projection reference metadata")
    p.add_argument("--target-metadata", type=Path, help="Optional external target metadata snapshot")
    p.add_argument("--provenance-manifest", type=Path, help="Historical reference-metadata receipt")
    p.add_argument("--target-info", type=Path, help="Exact GRCh37 target variant quality table; see PCA guide")
    p.add_argument("--min-info", type=float, default=.8)
    p.add_argument("--min-overlap", type=float, default=.9)
    p.add_argument("--min-call-rate", type=float, default=.9)
    p.add_argument("--allow-unknown-info", action="store_true", help="Explicit exploratory projection without verified imputation INFO; audited")
    p.add_argument("--allow-palindromic", action="store_true", help="Explicitly retain A/T and C/G markers after external strand verification")
    return p


def metadata_path(a, prefix):
    if a.metadata:
        return a.metadata.expanduser().resolve()
    if a.command == "render":
        raise ValueError("render requires --metadata; no panel or ancient data are needed without overlays")
    if a.pgen_prefix and not a.supplement_xlsx:
        raise ValueError("Custom --pgen-prefix requires --metadata (or an explicit compatible --supplement-xlsx)")
    from evolve import prepare_metadata
    directory = a.out / "metadata"
    metadata = directory / "sample_metadata.tsv"
    receipt = directory / "pca3d_metadata.json"
    supplement = a.supplement_xlsx or a.ancient_dir / "41586_2026_10358_MOESM3_ESM.xlsx"
    psam = Path(str(prefix) + ".psam")
    inputs = {"psam": _fingerprint(psam), "supplement": _fingerprint(supplement),
              "qpadm": _fingerprint(a.qpadm) if a.qpadm else None, "ancestry_model": a.ancestry_model}
    if receipt.exists():
        saved = json.loads(receipt.read_text())
        if saved["inputs"] != inputs or any(not Path(path).is_file() or _sha(path) != digest
                for path, digest in saved["outputs"].items()):
            raise ValueError("Metadata inputs/output changed; use a new --out or supply verified --metadata")
        return metadata
    if directory.exists() and any(directory.iterdir()):
        raise ValueError("Existing metadata has no pca3d receipt; pass --metadata explicitly or use a new --out")
    result = prepare_metadata(psam, supplement, directory, qpadm_path=a.qpadm, ancestry_model=a.ancestry_model)
    outputs = {result[key]: _sha(result[key]) for key in ("metadata_path", "qc_path", "manifest_path")}
    _atomic_text(receipt, json.dumps(dict(inputs=inputs, outputs=outputs), indent=2) + "\n")
    return metadata


def overlay_options(a, prefix):
    if not a.selected_variants:
        return {}
    required = {"variant_id", "pvar_index", "chrom", "pos", "ref", "alt", "status"}
    selected, total, seen = [], 0, set()
    for row in _rows(a.selected_variants, required):
        total += 1
        if row["variant_id"] in seen:
            raise ValueError("Duplicate overlay variant_id: " + row["variant_id"])
        seen.add(row["variant_id"])
        if len(selected) < a.pca_max_variants:
            selected.append(row)
    if not selected:
        raise ValueError("--selected-variants is empty")
    folder = a.out / "pca"
    folder.mkdir(parents=True, exist_ok=True)
    bounded = folder / "overlay_targets.tsv"
    if bounded.resolve() == a.selected_variants.expanduser().resolve():
        raise ValueError("Overlay input must not be the generated overlay_targets.tsv; use a separate input file")
    with bounded.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(selected[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected)
    _atomic_text(folder / "overlay_selection.json", json.dumps(dict(
        rule="first rows in explicit input order; bounded visualization only", total_available=total,
        displayed=len(selected), max_variants=a.pca_max_variants, source=_fingerprint(a.selected_variants)), indent=2) + "\n")
    return dict(pgen_path=Path(str(prefix) + ".pgen"), psam_path=Path(str(prefix) + ".psam"),
                selected_variants_path=bounded)


def workflow(a):
    if a.pca_layer_years < 1:
        raise ValueError('--pca-layer-years must be positive')
    if a.command == 'prepare-basis':
        if not a.basis_dir or not a.project_out:
            raise ValueError('prepare-basis requires original --basis-dir and separate --project-out')
        subprocess.run([sys.executable,str(Path(__file__).with_name('pca_basis.py')),
            '--basis-dir',str(a.basis_dir),'--out',str(a.project_out),'--plink2',a.plink2,
            '--threads',str(a.pca_threads),'--memory-mb',str(a.pca_memory_mb)],check=True)
        return dict(status='complete',operation='separate_unique_coordinate_projection_basis',output=str(a.project_out))
    if a.command in {"project-plan", "project"}:
        if not a.basis_dir:
            raise ValueError("External projection requires --basis-dir with a completed fixed PCA basis")
        destination = a.project_out or a.out / "ukb_projection"
        command = [sys.executable, str(Path(__file__).with_name("pca_project.py")),
            "plan" if a.command == "project-plan" else "run", "--basis-dir", str(a.basis_dir),
            "--target-glob", a.target_glob, "--out", str(destination), "--plink2", a.plink2,
            "--threads", str(a.pca_threads), "--memory-mb", str(a.pca_memory_mb),
            "--min-info", str(a.min_info), "--min-overlap", str(a.min_overlap),
            "--min-call-rate", str(a.min_call_rate)]
        reference_metadata = a.reference_metadata or a.metadata
        if reference_metadata: command += ["--reference-metadata", str(reference_metadata)]
        if a.target_metadata: command += ["--target-metadata", str(a.target_metadata)]
        if a.provenance_manifest: command += ["--provenance-manifest", str(a.provenance_manifest)]
        if a.target_info: command += ["--target-info", str(a.target_info)]
        if a.allow_unknown_info: command += ["--allow-unknown-info"]
        if a.allow_palindromic: command += ["--allow-palindromic"]
        subprocess.run(command, check=True)
        return dict(status="planned" if a.command == "project-plan" else "complete",
                    operation="external_fixed_basis_projection", output=str(destination))
    if a.command == "atlas":
        command = [sys.executable, str(Path(__file__).with_name("pca_atlas.py")), "--out", str(a.out),
            "--time-scale", a.pca_time_scale, "--max-age-bp", str(a.pca_display_max_age_bp),
            "--layer-years", str(a.pca_layer_years)]
        if a.source_html:
            if a.eigenvec or a.metadata:
                raise ValueError("atlas uses either --source-html or --eigenvec plus --metadata")
            if a.pca_basis_label:
                raise ValueError("Source HTML retains its original basis label; --pca-basis-label is for shared coordinate input")
            command += ["--source-html", str(a.source_html)]
        elif a.eigenvec and a.metadata:
            command += ["--coordinates", str(a.eigenvec), "--metadata", str(a.metadata)]
            if a.pca_basis_label:
                command += ["--basis-label", a.pca_basis_label]
        else:
            raise ValueError("atlas needs --source-html or both --eigenvec and --metadata")
        if a.simulated:
            command.append("--simulated")
        subprocess.run(command, check=True)
        return dict(status="complete", operation="display_only_existing_coordinates", output=str(a.out))
    if a.source_html:
        raise ValueError("Use atlas for --source-html re-rendering")
    if a.command == "render" and not a.eigenvec:
        raise ValueError("render requires --eigenvec with shared PCA coordinates")
    if a.reuse_basis_dir:
        if a.command != "render" or a.eigenvec.resolve() != (a.reuse_basis_dir / "all.sscore").resolve():
            raise ValueError("Completed-basis reuse requires render with that basis's all.sscore")
        from pca_project import validate_completed_coordinates
        validate_completed_coordinates(a.reuse_basis_dir)
        reused_plan = json.loads((a.reuse_basis_dir / "pca_plan.json").read_text())
        saved = reused_plan["inputs"]["metadata"]
        if not a.metadata:
            raise ValueError("Completed-basis rendering requires its verified metadata")
        from pca_contract import historical_identity
        provenance = historical_identity(a.metadata, recorded=saved)
        if provenance['historical_metadata_verification'] != 'verified':
            raise ValueError('Historical metadata lacks a trusted full SHA256 receipt; default verified render cannot proceed')
        _atomic_text(a.out / 'metadata_provenance.json', json.dumps(provenance, indent=2) + '\n')
        a.pca_basis_label = a.pca_basis_label or reused_plan["settings"]["basis"]
    if a.command != "render" and a.eigenvec:
        raise ValueError("Use render for existing --eigenvec coordinates")
    prefix = (a.pgen_prefix or a.ancient_dir / PREFIX).expanduser().resolve()
    metadata = metadata_path(a, prefix)
    plan = None
    if a.command != "render":
        if a.pgen_prefix and not a.hq_variants:
            raise ValueError("Custom panels require --hq-variants; background QC cannot be inferred from a PGEN")
        if a.dir_pgen_cache:
            from allele_common import cache_pgen_prefix
            prefix, cache_record = cache_pgen_prefix(prefix, a.dir_pgen_cache)
            _atomic_text(a.out / "pca_input_cache.json", json.dumps(cache_record, indent=2) + "\n")
        plan = prepare_pca_plan(prefix, metadata, a.hq_variants or a.ancient_dir / "HighQualityVariants.tsv",
            a.out / "pca", basis=a.pca_basis, maf=a.pca_maf, geno=a.pca_geno,
            n_components=a.pca_components, max_per_group=a.pca_max_per_group,
            time_bin_years=a.pca_time_bin_years, seed=a.seed, exclude_bed=a.pca_exclude_bed,
            target_variants=a.target_variants or a.selected_variants, target_flank_bp=a.pca_target_flank_bp,
            plink2=a.plink2, threads=a.pca_threads, memory_mb=a.pca_memory_mb,
            basis_min_age_bp=0 if a.pca_basis == "ancient_balanced" and a.pca_basis_min_age_bp is None else a.pca_basis_min_age_bp,
            basis_max_age_bp=10000 if a.pca_basis == "ancient_balanced" and a.pca_basis_max_age_bp is None else a.pca_basis_max_age_bp)
        if a.command == "plan":
            return dict(status="planned", manifest_path=plan["manifest_path"], script_path=plan["script_path"])
        subprocess.run(["bash", plan["script_path"]], check=True)
    source = plan["projected_path"] if plan else a.eigenvec
    label = a.pca_basis_label or (a.pca_basis if plan else "User-supplied shared PCA basis")
    overlay = overlay_options(a, prefix)
    result = render_pca_time(source, metadata, a.html or a.out / "pca" / "pca_time.html",
        basis_label=label, max_variants=a.pca_max_variants, simulated=a.simulated, title=a.title,
        count_kind=a.count_kind, layer_years=a.pca_layer_years, time_view=a.pca_time_view,
        color_by=a.pca_color_by, **overlay)
    if not a.skip_atlas:
        subprocess.run([sys.executable, str(Path(__file__).with_name("pca_atlas.py")),
            "--source-html", result["html_path"], "--out", str(a.out / "pca" / "atlas"),
            "--time-scale", a.pca_time_scale, "--max-age-bp", str(a.pca_display_max_age_bp),
            "--layer-years", str(a.pca_layer_years)], check=True)
    if overlay:
        centroid = run_centroids
        from types import SimpleNamespace
        dosages=a.out/'pca'/'selected_variant_dosages.tsv.gz'
        export_dosages(overlay['pgen_path'],overlay['psam_path'],overlay['selected_variants_path'],dosages,a.count_kind)
        centroid(SimpleNamespace(coordinates=source,metadata=metadata,dosages=dosages,out=a.out/'pca'/'centroids',bin_years=a.pca_time_bin_years))
    return dict(status="complete", **{k: result[k] for k in ("html_path", "manifest_path", "samples", "variants")})


def main_pca(argv=None):
    a = parser().parse_args(argv)
    a.out = a.out.expanduser().resolve()
    a.ancient_dir = a.ancient_dir.expanduser().resolve()
    try:
        context = nullcontext()
        if a.command == 'run':
            if a.eigenvec or a.reuse_basis_dir or a.source_html:
                raise ValueError('Use render/atlas for existing coordinates or HTML; run recalculates from raw genotypes')
            from fresh_run import output_run
            prefix = a.pgen_prefix or a.ancient_dir / PREFIX
            protected = [Path(str(prefix) + suffix) for suffix in ('.pgen', '.pvar', '.pvar.zst', '.psam')]
            protected += [a.metadata, a.hq_variants or a.ancient_dir / 'HighQualityVariants.tsv',
                          a.supplement_xlsx or a.ancient_dir / '41586_2026_10358_MOESM3_ESM.xlsx',
                          a.qpadm, a.target_variants, a.selected_variants, a.pca_exclude_bed]
            context = output_run(a.out, a.replace == 'TRUE', protected)
        with context:
            a.out.mkdir(parents=True, exist_ok=True)
            with (a.out / ".pca3d.lock").open("a") as lock:
                try:
                    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    raise ValueError("Another pca3d job is already running in this output directory") from None
                result = workflow(a)
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"pca3d: {exc}", file=sys.stderr)
        return 2

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'centroid':
        raise SystemExit(main_centroid(sys.argv[2:]))
    raise SystemExit(main_pca())
