#!/usr/bin/env python3
"""Publish manuscript tables and figures from the allele3 analysis products."""
import argparse
from collections import defaultdict
from contextlib import contextmanager
import fcntl
import html
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
from urllib.parse import quote

from paper_sources import discover, envdir
from research_common import html_table, number, write_tsv

MARKER = 'allele3-manuscript-2'
TITLES = {
    'T4b_paired_near_far_by_alt_count':'Proximal–distal differences by number of qualifying ALTs',
    'T4c_masked_near_far':'Sensitivity excluding masked focal sites and controls',
    'T5b_direct_splicing':'Within-site geometry using the direct AG_SPLICE prediction',
    'T8_matched_cohort_enrichment':'Shared observed-site support at AF ≥ 0.0001, by cohort and population',
    'T8b_matched_alt_similarity':'Shared-site AF-ranked ALT distances by population and score',
    'T8c_mean_adjusted_gradient':'Gradient tails after matching mean PHRED AVI',
    'T8d_matched_distance':'Paired proximal–distal contrast on shared cohort support',
    'T8e_sensitivity_eligibility':'Denominators and prediction availability for additional analyses',
    'S27_study_status':'Prespecified analyses, unavailable inputs and applied masks',
    'S28_study_inputs':'Input paths and availability for the downstream workflow',
    'T7b_carrier_genotypes':'Individual genotype and archaic-tract overlap, summarized by allele role',
    'T1_scan_counts':'Observed sites by cohort, population and qualification threshold',
    'T1b_base_evidence':'Evidence for three or four observed bases',
    'T1c_analysis_denominators':'Eligible sites and cases in each reference-specific test',
    'T2b_gradient_heterogeneity':'Three-ALT versus two-ALT O/E and high-versus-low gradient interaction',
    'T2_region_enrichment':'Region enrichment by outcome, cohort and reference',
    'T3_distance_profile':'AVI at successive distances from matched cases and controls',
    'T4_paired_near_far':'Paired near-minus-far AVI contrast',
    'T5_within_site':'AF-ranked, within-site AVI contrasts',
    'T6_herv_family':'HERV family enrichment in the detailed scan',
    'T6b_common_support':'HERV family and pooled estimates on supported background strata',
    'T6c_boundary_validation':'Discovery and held-out chromosome estimates at frozen repeat boundaries',
    'T7_archaic_counters':'Archaic tract enrichment from the detailed-scan workflow',
    'S1_pca_counts':'Samples displayed on the PC1, PC2 and age axes',
    'S2_site_capabilities':'Source capabilities by population',
    'S3_herv_scan':'Personal-genome HERV family stability',
    'S4_avi_cutoffs':'AVI window definitions and coverage',
    'S5_hotspot_windows':'Window tests and nominations by scale',
    'S6_hotspot_loci':'Merged nominated loci and their lead windows',
    'S7_density':'Density score coverage by scale, threshold and channel',
    'S8_allele_time':'Fixed-panel ancient allele trajectories',
    'S9_space':'PC1 × PC2 × time sample coverage',
    'S10_simulation':'Simulated persistence across replicates',
    'S11_test_availability':'Upstream test denominators and exclusion counts',
    'S12_selection':'Overlap with the published ancient selection catalogue',
    'S13_joint_time':'Measurability of the jointly measured ancient panel',
    'S14_population_patterns':'Within-population and pooled multi-ALT evidence',
    'S16_herv_spatial':'HERV enrichment by aligned element partition',
    'S17_ld':'HERV and population LD association',
    'S18_genotypes':'Genotype retention and measurable A/C/G/T states',
    'S19_fixed_panel_time':'Ancient case-control fixed-panel time summaries',
    'S22_region_catalogues':'Number of region candidates by evidence tier',
    'S24_herv_time_models':'HERV paired time-association estimates with chromosome sensitivity',
    'S25_herv_partitions':'Prespecified head/body and tail/body comparisons',
    'S26_fixed_pairs':'Frozen HERV case-control pairs by membership and allele role',
    'S23_ancient_coverage':'Representation of the modern catalogue in ancient genotypes',
    'S20_window_contrasts':'Matched-window annotation contrasts',
    'S21_pca_stability':'Sensitivity to the choice of PCA basis',
    'S15_archaic_carriers':'Within-person archaic tract and allele co-occurrence',
}


def band(value):
    try:
        lo, hi = map(int, value.split(','))
        if 0 <= lo < hi:
            return lo, hi
    except ValueError:
        pass
    raise argparse.ArgumentTypeError('Use two ordered bp limits, e.g. 0,1000')


def path_arg(value):
    return Path(value).expanduser().resolve()


def options(argv):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis-root', type=path_arg, default=envdir('ALLELE_ANALYSIS_ROOT', Path(os.environ.get('ANALYSIS_ROOT','.'))/'allele3'),
                   help='Project analysis directory containing main/ and 3.allele3_scan/')
    p.add_argument('--result-root', type=path_arg, help='Exact directory containing the AVI, HERV, archaic and hotspot modules.')
    p.add_argument('--out', type=path_arg, help='Paper directory [paper/ inside the project]')
    p.add_argument('--scan-root', type=path_arg, help='3.allele3_scan directory containing gnomad/ and allofus/')
    p.add_argument('--pca-root', type=path_arg, help='Permanent pca3d directory')
    p.add_argument('--herv-root', type=path_arg, help='Completed 1KG HERV scan or archive')
    p.add_argument('--avi-root', type=path_arg, help='Completed 4.avi_scan directory')
    p.add_argument('--targets', default='gnomad,allofus', help='Comma-separated gnomad,allofus,custom')
    p.add_argument('--near', type=band, default=(0,1000), help='Near band inner,outer bp [0,1000]')
    p.add_argument('--far', type=band, default=(100000,1000000), help='Far band inner,outer bp [100000,1000000]')
    p.add_argument('--bootstrap', type=int, default=1000, help='Chromosome bootstrap draws for T4 only [1000]')
    p.add_argument('--seed', type=int, default=20261115)
    p.add_argument('--grch', choices=['38'], default='38')
    a = p.parse_args(argv)
    a.targets = list(dict.fromkeys(t.strip() for t in a.targets.split(',')))
    if set(a.targets)-{'gnomad','allofus','custom'}:
        p.error('Choose --targets from gnomad,allofus,custom')
    if a.near[1] > a.far[0] or a.bootstrap < 20 or a.seed < 0:
        p.error('Use disjoint ordered bands, at least 20 bootstrap draws and a nonnegative seed')
    a.out = a.out or a.analysis_root/'paper'
    return a


@contextmanager
def publication(out):
    # Only the manuscript-owned main directory is replaced. Legacy scientific
    # inputs may coexist directly under paper/ on the first run.
    for name in ('jobs','logs','qc','tmp'):
        (out/name).mkdir(parents=True, exist_ok=True)
    dest = out/'main'
    backup = out/'tmp'/'previous_main'
    with (out/'qc'/'paper.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        def owned(path):
            return (path/'reports'/'.allele3-paper').is_file()
        if backup.exists():
            if not owned(backup):
                raise ValueError(f'Unexpected backup at {backup}')
            if not dest.exists():
                backup.rename(dest)
            elif owned(dest):
                shutil.rmtree(backup)
            else:
                raise ValueError(f'Inspect {dest} and {backup} before publishing')
        if dest.exists() and (not dest.is_dir() or (any(dest.iterdir()) and not owned(dest))):
            raise ValueError(f'Paper output contains an existing scientific main directory: {dest}')
        with tempfile.TemporaryDirectory(prefix='build-', dir=out/'tmp') as temp:
            stage = Path(temp)
            for folder in ('main/results','main/reports','qc'):
                (stage/folder).mkdir(parents=True)
            yield stage
            (stage/'main/reports/.allele3-paper').write_text(MARKER+'\n')
            if dest.exists():
                if owned(dest):
                    dest.rename(backup)
                else:
                    dest.rmdir()
            try:
                (stage/'main').rename(dest)
            except BaseException:
                if backup.exists():
                    backup.rename(dest)
                raise
            for audit in sorted((stage/'qc').iterdir()):
                audit.replace(out/'qc'/audit.name)
            if backup.exists():
                shutil.rmtree(backup)


class Report:
    def __init__(self, a, stage, sources):
        self.a, self.stage, self.sources = a, stage, sources
        self.tables = defaultdict(list)
        self.coverage, self.figures = [], {}

    def check(self, source, label, path, status, n=None, detail=''):
        self.coverage.append(dict(product=label, source=source.key, status=status, rows=n,
                                  detail=detail or source.problem, path=str(path or source.root)))

    def table(self, key, source, names, required=(), transform=None):
        from paper_tables import rows
        path = source.get(*names)
        if not path:
            stage = Path(names[0]).name.split('.')[0]
            detail = source.problem or ('Upstream stage unavailable: '+stage if stage in source.excluded_stages else 'Missing '+names[0])
            if source.profile_disabled and 'avi_distance' in names[0]:
                detail = 'Distance profiles disabled in upstream request'
            self.check(source,key,None,'unavailable',detail=detail)
            return
        try:
            before = path.stat()
            records = list(transform(path) if transform else rows(path,required))
        except (ValueError, KeyError, OSError, EOFError) as ex:
            self.check(source,key,path,'unreadable',detail=str(ex))
            return
        after = path.stat()
        if (before.st_size,before.st_mtime_ns) != (after.st_size,after.st_mtime_ns):
            raise RuntimeError(f'Input changed during the read: {path}. Rerun after the producer finishes.')
        self.check(source,key,path,'present' if records else 'header_only',len(records))
        for row in records:
            self.tables[key].append(dict(row, source=source.key, target=source.target or row.get('target','NA'),
                                         receipt_status=source.status, reference=source.reference or row.get('reference','NA'),
                                         population=row.get('population',source.population),
                                         **({'gu_cohort':source.gu_cohort} if source.gu_cohort else {})))

    def calculate(self):
        from paper_tables import near_far, window_summary, simulation_summary, carrier_summary, time_models, fixed_pair_counts
        def stratified_pairs(path):
            from paper_tables import rows
            import csv
            from research_common import open_text
            with open_text(path) as handle:
                if 'qualifying_alts' not in (csv.DictReader(handle,delimiter='\t').fieldnames or ()):
                    return
            for k in (2,3):
                for row in near_far(path,self.a.near,self.a.far,self.a.bootstrap,self.a.seed,k):
                    yield dict(row,qualifying_alts=k)
        def direct_splicing(path):
            from paper_tables import rows
            return (r for r in rows(path) if r['scale']=='SPLICE')
        self.table('T7b_carrier_genotypes',self.sources['archaic_carriers'],['results/role_by_tract.tsv'])
        inventory=self.sources['study_inventory']
        self.table('S27_study_status',inventory,['study_status.tsv'])
        self.table('S28_study_inputs',inventory,['study_inputs.tsv'])
        for key,src in self.sources.items():
            if not key.startswith('study/'):continue
            if key.endswith('_mean_adjusted'):
                self.table('T8c_mean_adjusted_gradient',src,['results/enrichment.tsv'])
            elif src.reference=='alleles':
                self.table('T8b_matched_alt_similarity',src,['results/within_site_tests.tsv'])
                self.table('T8e_sensitivity_eligibility',src,['results/geometry_status.tsv'])
            else:
                self.table('T8_matched_cohort_enrichment',src,['results/enrichment.tsv'])
                self.table('T8e_sensitivity_eligibility',src,['results/status_counts.tsv'])
                if src.reference=='alphagenome':
                    for group,suffix in [('all',''),('exclude_mask','_exclude_mask')]:
                        self.table('T8d_matched_distance',src,['results/avi_distance_pairs'+suffix+'.tsv.gz'],
                                   transform=lambda p,g=group:(dict(r,sensitivity=g) for r in near_far(p,self.a.near,self.a.far,self.a.bootstrap,self.a.seed)))
        for target in self.a.targets:
            scan = self.sources['scan/'+target]
            for table,name in [('T1_scan_counts','population_summary'),('T1b_base_evidence','allele_evidence_summary'),
                               ('S2_site_capabilities','source_capabilities'),('S14_population_patterns','population_patterns')]:
                self.table(table,scan,[f'results/{name}.tsv'])
            region_sources=[s for s in self.sources.values() if s.target==target and
                            (s.key in ('alphagenome/'+target,'herv/'+target) or s.gu_cohort)]
            for s in region_sources:
                ref=s.reference
                self.table('T2_region_enrichment',s,['results/enrichment.tsv'],('region','outcome','sensitivity','oe'))
                self.table('T1c_analysis_denominators',s,['results/status_counts.tsv'],('measure','n'))
                highlights = s.meta.get('region_highlights', {})
                if not s.problem and highlights:
                    for key, value in highlights.items():
                        if isinstance(value,(int,float)):
                            self.tables['S22_region_catalogues'].append(dict(source=s.key,target=target,reference=ref,
                                                                          population=s.population,gu_cohort=s.gu_cohort,measure=key,value=value))
                if ref=='alphagenome':
                    self.table('T2b_gradient_heterogeneity',s,['results/gradient_heterogeneity.tsv'],('contrast','ratio','fdr'))
            avi = self.sources['alphagenome/'+target]
            self.table('T3_distance_profile',avi,['results/avi_distance_summary.tsv'],
                       ('metric','inner_bp','outer_bp','case_minus_control'))
            self.table('T4_paired_near_far',avi,['results/avi_distance_pairs.tsv.gz'],
                       transform=lambda p:near_far(p,self.a.near,self.a.far,self.a.bootstrap,self.a.seed))
            self.table('T4c_masked_near_far',avi,['results/avi_distance_pairs_exclude_mask.tsv.gz'],
                       transform=lambda p:near_far(p,self.a.near,self.a.far,self.a.bootstrap,self.a.seed))
            self.table('T4b_paired_near_far_by_alt_count',avi,['results/avi_distance_pairs.tsv.gz'],transform=stratified_pairs)
            allele = self.sources['allele_test/'+target]
            self.table('T5b_direct_splicing',allele,['results/within_site_tests.tsv'],transform=direct_splicing)
            self.table('T5_within_site',allele,['results/within_site_tests.tsv'],
                       transform=lambda p:(r for r in __import__('paper_tables').rows(p,('metric','scale','sensitivity','estimate')) if r['scale'] in ('RAW','PHRED')))
            for name in ('matching_status','geometry_status'):
                self.table('S11_test_availability',allele,[f'results/{name}.tsv'],
                           transform=lambda p,n=name:tagged(p,'audit',n))
        for key,module,names in [
            ('T6_herv_family','herv_enrichment',['results/primary.enrichment_summary.tsv']),
            ('T6b_common_support','herv_enrichment',['results/summary.common_support_enrichment.tsv']),
            ('T6c_boundary_validation','herv_enrichment',['results/windows.frozen_peak_validation.tsv']),
            ('T7_archaic_counters','herv_enrichment',['results/archaic.archaic_enrichment.tsv']),
            ('S16_herv_spatial','herv_enrichment',['results/spatial.primary_spatial_partitions.tsv']),
            ('S17_ld','herv_enrichment',['results/ld.association.tsv']),
            ('S18_genotypes','herv_enrichment',['results/genotypes.cohort_retention.tsv']),
            ('S19_fixed_panel_time','allele_time',['results/fixed_panel_time_summary.tsv']),
            ('S25_herv_partitions','herv_enrichment',['results/spatial.spatial_contrasts.tsv']),
            ('S23_ancient_coverage','herv_enrichment',['results/time.catalog_coverage.tsv']),
            ('S23_ancient_coverage','allele_time',['results/catalog_coverage.tsv']),
            ('S20_window_contrasts','hotspots',['results/characterize.annotation_contrasts.tsv']),
            ('S21_pca_stability','pca_study',['stability.tsv']),
            ('S3_herv_scan','herv_scan',['results/family_stability.tsv']),
            ('S4_avi_cutoffs','avi_scan',['results/cutoffs.tsv']),
            ('S6_hotspot_loci','hotspots',['results/windows.loci.tsv']),
            ('S8_allele_time','allele_time',['results/fixed_panel_family_trajectories.tsv']),
            ('S8_allele_time','herv_enrichment',['results/time.fixed_panel_family_trajectories.tsv']),
            ('S9_space','space',['results/coordinate_cells.tsv']),
            ('S12_selection','herv_enrichment',['results/selection.published_selection_summary.tsv']),
            ('S13_joint_time','herv_enrichment',['results/time.joint_measurability_funnel.tsv'])]:
            self.table(key,self.sources[module],names)
        for module,filename,key in [('hotspots','windows.all_windows','S5_hotspot_windows'),
                                    ('density','density_windows','S7_density')]:
            self.table(key,self.sources[module],[f'results/{filename}.tsv.gz'],
                       transform=lambda p,d=module=='density':window_summary(p,d))
        self.table('S10_simulation',self.sources['base_simulation'],['results/summary.tsv'],transform=simulation_summary)
        if 'occupancy/archaic_alleles' in self.sources:
            self.table('S15_archaic_carriers',self.sources['occupancy/archaic_alleles'],
                       ['results/allele_tract_summary.tsv','allele_tract_summary.tsv'],transform=carrier_summary)
        herv = self.sources['herv_enrichment']
        self.table('S24_herv_time_models',herv,['results/time.joint_person_models.tsv'],transform=time_models)
        self.table('S26_fixed_pairs',herv,['results/time.joint_fixed_pair_panel.tsv'],transform=fixed_pair_counts)
        self.pca()

    def pca(self):
        from paper_tables import rows
        from paper_figures import pca_figure
        source = self.sources['pca']
        path = source.get('tables/atlas/original_points.tsv.gz','pca/atlas/original_points.tsv.gz',
                          'atlas/original_points.tsv.gz','original_points.tsv.gz')
        if path:
            points = rows(path,['pc1','pc2','age_bp','cohort'])
        else:
            coordinates = source.get('coordinates/all.sscore','pca/all.sscore')
            metadata = source.get('tables/sample_metadata.tsv','metadata/sample_metadata.tsv')
            if coordinates and metadata:
                from pca_plot import read_coordinate_source
                path = coordinates
                try:
                    payload,_ = read_coordinate_source(coordinates,metadata)
                except (ValueError,KeyError,OSError) as ex:
                    self.check(source,'S1_pca_counts',path,'unreadable',detail=str(ex))
                    return
                points = payload['points']
            else:
                # The space table is an explicitly different sample frame.
                source = self.sources['space']
                path = source.get('metadata/sample_audit.tsv.gz','results/sample_audit.tsv.gz')
                if not path:
                    self.check(self.sources['pca'],'S1_pca_counts',None,'unavailable',
                               detail='No published PC points, usable coordinate/metadata pair or space sample audit')
                    return
                points = (dict(r,pc1=r['PC1'],pc2=r['PC2']) for r in rows(path,['PC1','PC2','age_bp','cohort','eligible'])
                          if str(r['eligible']).lower() in ('1','true'))
        try:
            result = pca_figure(self.stage/'main/reports',points)
        except (ValueError,KeyError,OSError,EOFError) as ex:
            self.check(source,'S1_pca_counts',path,'unreadable',detail=str(ex))
            return
        self.check(source,'S1_pca_counts',path,'present' if result else 'no_eligible_points',len(result))
        if result:
            self.tables['S1_pca_counts'] = [dict(r,source=source.key) for r in result]
            self.figures['S1_pca_time'] = f'PC1 and PC2 are genetic coordinates. The vertical axis uses recorded sample age. Sample frame: {source.key}.'

    def native_figures(self):
        for key,module,names,title in [
            ('S2_herv_families','herv_enrichment',['reports/primary.primary_enrichment'],'Detailed-scan HERV family enrichment'),
            ('S12_archaic_enrichment','herv_enrichment',['reports/archaic.archaic_enrichment'],'Archaic tract enrichment from the detailed scan'),
            ('S3_herv_spatial','herv_enrichment',['reports/spatial.primary_spatial_partitions'],'Spatial partitions of annotated HERV elements'),
            ('S4_density','density',['reports/density_scores'],'Genome-window density scores'),
            ('S5_simulation','base_simulation',['reports/base_simulation'],'Simulated four-base dynamics'),
            ('S6_allele_time','allele_time',['reports/family_allele_time','reports/ancient_fixed_panel_time'],'Ancient fixed-panel allele measurements'),
            ('S7_space','space',['reports/allele_space_examples'],'Allele distributions in the sampled PC/time space'),
            ('S8_herv_stability','herv_scan',['figures/family_stability'],'Personal-genome HERV stability in 1KG'),
            ('S9_avi_genome','avi_scan',['figures/avi_landscape'],'Genome-wide AVI landscape'),
            ('S10_pca_stability','pca_study',['figures/pca_time_stability'],'PCA basis sensitivity across time'),
            ('S11_herv_time','herv_enrichment',['reports/time.family_allele_time'],'HERV fixed-panel allele trajectories'),
        ]:
            source = self.sources[module]
            for name in names:
                found = False
                for ext in ('png','svg','pdf'):
                    path = source.get(name+'.'+ext)
                    if path:
                        shutil.copyfile(path,self.stage/'main/reports'/(key+'.'+ext))
                        found = True
                        break
                if found:
                    self.figures[key] = title + ' (upstream figure).'
                    break

    def write(self):
        from paper_figures import figures
        self.figures = {**figures(self.stage/'main/reports',self.tables),**self.figures}
        self.native_figures()
        for key, data in self.tables.items():
            if data:
                fields = list(dict.fromkeys(k for r in data for k in r))
                write_tsv(self.stage/'main/results'/(key+'.tsv'),fields,data)
        write_tsv(self.stage/'qc/coverage.tsv',['product','source','status','rows','detail','path'],self.coverage)


def tagged(path, key, value):
    from paper_tables import rows
    return (dict(r,**{key:value}) for r in rows(path))


def link(path, label, out):
    href = quote(os.path.relpath(path,out), safe='/')
    return '<a href="'+html.escape(href,quote=True)+'">'+html.escape(label)+'</a>'


def status_overview(rep):
    evidence = []
    def describe(name, data, effect):
        finite = sum(number(r.get(effect)) is not None for r in data)
        return dict(question=name, rows=len(data), finite_estimates=finite,
                    with_95pct_interval=sum(number(r.get('ci_low')) is not None and number(r.get('ci_high')) is not None for r in data),
                    availability='estimated' if finite else 'no estimable row in supplied products')
    for target in rep.a.targets:
        subset = lambda key: [r for r in rep.tables.get(key,[]) if r['target']==target]
        for ref in ('alphagenome','herv','archaic'):
            result = [r for r in subset('T2_region_enrichment') if r['reference']==ref and r['outcome']=='multi_ALT']
            evidence.append(dict(target=target,**describe(ref+' region enrichment',result,'oe')))
        result = [r for r in subset('T4_paired_near_far') if r['metric'] in ('sd_site_mean','mean_abs_gradient')]
        evidence.append(dict(target=target,**describe('Paired near/far AVI variability',result,'estimate')))
        result = [r for r in subset('T5_within_site') if r['sensitivity']=='all_available' and r['metric']=='normalized_weak_pair_distance_minus_other_pairs']
        evidence.append(dict(target=target,**describe('Lower-AF pair similarity, three observed ALTs',result,'estimate')))
    for name,key in [('Detailed-scan HERV families','T6_herv_family'),('Detailed-scan archaic tracts','T7_archaic_counters')]:
        evidence.append(dict(target='gnomad (detailed)', **describe(name, rep.tables.get(key,[]),'oe')))
    return evidence


def finish(rep, main, notes):
    a, stage = rep.a, rep.stage
    reports = a.out/'main/reports'
    records, modules = [], []
    for s in rep.sources.values():
        files = sorted(s.files.items())
        modules.append(dict(module=s.key,status=s.status,issue=s.problem or ', '.join(sorted(s.excluded_stages)),
                            n_indexed=len(files),path=str(s.root),receipt=str(s.receipt or '')))
        for name,path in files:
            records.append(dict(module=s.key,product=name,path=str(path),bytes=path.stat().st_size,
                                status=s.problem or ('stage unavailable' if Path(name).name.split('.')[0] in s.excluded_stages else s.status)))
    write_tsv(stage/'qc/sources.tsv',['module','product','path','bytes','status'],records)
    write_tsv(stage/'qc/modules.tsv',['module','status','issue','n_indexed','path','receipt'],modules)
    evidence = status_overview(rep)
    write_tsv(stage/'main/results/hypotheses.tsv',list(evidence[0]),evidence)
    settings = {}
    for s in rep.sources.values():
        settings[s.key] = dict(path=str(s.root),receipt=str(s.receipt or ''),status=s.status,
                               parameters={k:s.meta[k] for k in ('population','parameters','context','genome_build',
                               'build','target','ref','background','threshold_kind','score_parameters','eligibility',
                               'reference','request','chromosomes') if k in s.meta})
    run = dict(schema=MARKER, analysis_root=str(a.analysis_root), result_root=str(main), targets=a.targets,
               near_bp=a.near, far_bp=a.far, bootstrap=a.bootstrap, seed=a.seed, notes=notes, inputs=settings)
    (stage/'qc/run.json').write_text(json.dumps(run,indent=2)+'\n')
    document = ['''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>Multi-allelic variation | Manuscript results</title><style>
body{font:16px/1.48 system-ui;color:#163740;max-width:1200px;margin:auto;padding:2em}h1,h2{color:#123f49}
h2{border-top:1px solid #cbdcda;padding-top:.7em}img{max-width:min(850px,100%);height:auto}a{color:#07676c}
table{border-collapse:collapse;font-size:13px}td,th{padding:.4em .65em;border-bottom:1px solid #dce6e6;text-align:left}
td{overflow-wrap:anywhere}summary{cursor:pointer}figure{margin:1em 0 2em}.scroll{overflow-x:auto}
</style><h1>Multi-allelic variation in the human genome</h1>''']
    document += ['<p>Analysis products from '+html.escape(str(main))+'.</p>']
    for message in notes:
        document.append('<p>'+html.escape(message)+'</p>')
    document.append('<nav><a href="#main">Main estimates</a> · <a href="#tables">Tables</a> · <a href="#figures">Figures</a> · <a href="#inputs">Source products and availability</a></nav>')
    document.append('<h2>Evidence by question and target</h2>'+html_table(evidence,list(evidence[0]),len(evidence)))
    document.append('<h2 id="main">Main estimates</h2>')
    for title,key,condition,columns in [
        ('AVI gradient tails by number of qualifying ALTs','T2_region_enrichment',
         lambda r:r['region'] in ('AVI_mean_abs_gradient_high','AVI_mean_abs_gradient_low') and r['sensitivity']=='all' and r['outcome'] in ('multi_ALT','two_ALT','three_ALT'),
         'target population region outcome observed expected oe ci_low ci_high p fdr'),
        ('Paired block test of the gradient pattern','T2b_gradient_heterogeneity',
         lambda r:r['sensitivity']=='all',
         'target population contrast ratio ci_low ci_high p fdr n_blocks'),
        ('Matched near-versus-far comparison','T4_paired_near_far', lambda r:r['metric'] in ('sd_site_mean','mean_abs_gradient','mean_alt_range'),
         'target population metric n_cases n_chromosomes case_change control_change estimate ci_low ci_high p fdr status'),
        ('Similarity between the two lower-AF alleles','T5_within_site',
         lambda r:r['metric']=='normalized_weak_pair_distance_minus_other_pairs',
         'target population scale sensitivity n_sites n_blocks estimate ci_low ci_high p fdr'),
        ('Two observed ALTs versus the unreported substitution','T5_within_site',
         lambda r:r['metric']=='normalized_extra_ALT_distance_to_top_minus_unreported',
         'target population scale sensitivity n_sites estimate ci_low ci_high fdr'),
        ('Frozen repeat-edge windows on held-out chromosomes','T6c_boundary_validation',
         lambda r:r.get('status')=='held_out_chromosome_confirmation',
         'family boundary start_bp end_bp validation_threshold validation_local_ratio validation_ci_low validation_ci_high validation_fdr'),
        ('HERV families and archaic tracts','T2_region_enrichment',lambda r:r['reference'] in ('herv','archaic') and r['outcome']=='multi_ALT' and r['sensitivity']=='all',
         'target population reference gu_cohort region observed expected oe ci_low ci_high p fdr'),
        ('Detailed HERV counters','T6_herv_family',lambda r:r.get('label') in ('family:ERV1','family:ERVK') and r.get('sensitivity')=='exclude_combined' and r.get('threshold') in ('5','10'),
         'label threshold sensitivity n_eligible observed_multi expected_multi oe ci_low ci_high bootstrap_p fdr status'),
    ]:
        data = [r for r in rep.tables.get(key,()) if condition(r)]
        document.append('<h3>'+title+'</h3>'+(('<div class="scroll">'+html_table(data,columns.split(),len(data))+'</div>') if data else '<p>No matching estimate in supplied results.</p>'))
    document.append(f'''<h2>Definitions and interpretation</h2><p>T2 reports O/E among eligible observed variant sites.
Two or three qualifying nonreference ALTs and proven three- or four-base occupancy are separate outcomes. T1b retains measured and unresolved evidence.
AF-only summaries have limited singleton and REF evidence. Repeated thresholds and related cohorts should be interpreted separately.</p>
<p>T3 reports the case-minus-control difference in each distance band. T4 compares {a.near[0]:,}–{a.near[1]:,} bp with {a.far[0]:,}–{a.far[1]:,} bp.
Its estimate is (case near − case far) − (control near − control far), restricted to cases with the same controls and weights in both bands.
T4 pools sites with at least two qualifying ALTs. T4b separates sites with exactly two and three qualifying ALTs when the upstream pair table records that count. Separate two-ALT and three-ALT region estimates are in T2 when those outcomes were run. The source excludes the focal tile. SD averages local tile SDs. The adjacent-position gradient measures spatial change, while ALT range measures variation among substitutions at one position.
Positive T4 values indicate greater proximal change in cases relative to controls. T4 uses {a.bootstrap:,} whole-chromosome bootstrap draws and BH across estimable metrics within each target.
T2b reports three-ALT O/E divided by two-ALT O/E in each gradient tail, and the ratio between those two contrasts. New scientific runs resample identical genomic blocks for all four O/E estimates. Legacy aggregate rows supply no paired-block uncertainty for this interaction.
Other intervals, P and FDR values keep their upstream definitions. All of Us is displayed separately when available.</p>
<p>T5 ranks the three qualifying observed ALTs by decreasing AF in the selected population. With dᵢⱼ = |AVIᵢ − AVIⱼ|, the main contrast is
[d₂₃ − (d₁₂ + d₁₃)/2] / (max AVI − min AVI). Negative values favor closer scores for the lower-AF pair.
AF ties and zero ranges are excluded from this contrast. REF has no ALT impact score.
These are distances between scalar, reference-relative predictions. AF ranks do not date mutations. AVI incorporates gnomAD-frequency information.
Archaic interval overlap and simulated persistence on their own do not establish allele ancestry or selection.</p>''')
    document.append('<h2 id="tables">Manuscript tables</h2><ul>')
    document.append('<p>T8 series: the two cohorts share eligible GRCh38 positions and REF, positive evidence in each source, the AF cutoff, population label and mask. Each cohort retains its own observed ALT counts and AF ordering. Shared positions are an observed-variant denominator. For ancestry-specific matches, the same label does not guarantee identical sampling. T5b uses the direct splicing prediction wherever three scores are available; T8b includes the same validation on shared sites. Missing predictions and equal scores are counted separately.</p>')
    for key,title in TITLES.items():
        if rep.tables.get(key):
            document.append('<li>'+link(a.out/'main/results'/(key+'.tsv'),key+' · '+title,reports)+f' ({len(rep.tables[key]):,} rows)</li>')
    document.append('</ul><h2 id="figures">Manuscript figures</h2>')
    for key,caption in rep.figures.items():
        links = [link(reports/p.name,p.suffix[1:].upper(),reports) for p in sorted((stage/'main/reports').glob(key+'.*'))]
        preview = next((p for ext in ('png','svg') if (p:=stage/'main/reports'/(key+'.'+ext)).is_file()),None)
        img = f'<img loading="lazy" alt="{html.escape(key)}" src="{preview.name}">' if preview else ''
        document.append(f'<figure>{img}<figcaption><b>{html.escape(key)}</b> · {html.escape(caption)} · {" / ".join(links)}</figcaption></figure>')
    document.append('<h2 id="inputs">Source products and availability</h2><p>')
    for name in ('coverage.tsv','hypotheses.tsv','modules.tsv','sources.tsv','run.json'):
        path = a.out/('main/results' if name=='hypotheses.tsv' else 'qc')/name
        document.append(link(path,name,reports)+' &nbsp; ')
    document.append('</p><p>Files-only products have no completion receipt. Non-complete research tests are excluded from estimates. Source links open the original files and require the analysis tree. Full site, window and allele catalogues are indexed in place.</p>')
    document.append('<p>QC coverage records missing, header-only and unreadable inputs separately. The sources.tsv inventory records each original file path.</p>')
    document.append('</html>')
    (stage/'main/reports/index.html').write_text('\n'.join(document),encoding='utf-8')


def main(argv=None):
    a = options(argv)
    if a.analysis_root.is_relative_to(a.out) or a.out == a.analysis_root/'main':
        raise ValueError('Choose the project paper/ directory or a dedicated manuscript directory')
    with publication(a.out) as stage:
        main_root,sources,notes = discover(a)
        for source in sources.values():
            if a.out.is_relative_to(source.root) or source.root.is_relative_to(a.out/'main'):
                raise ValueError('Paper output overlaps scientific input: '+str(source.root))
        report = Report(a,stage,sources)
        report.calculate()
        report.write()
        count = sum(bool(x) for x in report.tables.values())
        if not count and not report.figures:
            raise ValueError('No scientific tables or figures found. Checked ' + ', '.join(str(s.root) for s in sources.values()))
        finish(report,main_root,notes)
    print(f'{count} manuscript tables, {len(report.figures)} figures. {a.out / "main/reports/index.html"}')
    return 0


if __name__=='__main__':
    try:
        sys.exit(main())
    except (ValueError,OSError,RuntimeError) as ex:
        sys.exit('allele3.paper: '+str(ex))
