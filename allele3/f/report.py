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

# Scientific figures, evolution reports and the combined results report.





# ========================================================================
# Final figures
# ========================================================================

# Static scientific exports and an explicit availability inventory for Fig1–9.

from pathlib import Path
from collections import defaultdict, Counter
import csv, json, math, re
from allele_common import dump_json, fingerprint, open_text, tsv_writer

def _plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt

def _save(fig,path,formats=('png','pdf')):
    paths=[]
    for ext in formats:
        p=Path(str(path)+'.'+ext); fig.savefig(p,dpi=180,bbox_inches='tight'); paths.append(str(p))
    _plot().close(fig); return paths

def figure_rows(path):
    with open_text(path) as f: yield from csv.DictReader(f,delimiter='\t')

def pca_display(points, layer_years=1000, time_view='layered', color_by='region'):
    """Display-only time bins; original ages and the shared PCs stay intact."""
    import colorsys
    if not isinstance(layer_years, int) or isinstance(layer_years, bool) or layer_years < 1:
        raise ValueError('pca-layer-years must be a positive integer')
    if time_view not in {'layered', 'continuous'}:
        raise ValueError('pca-time-view must be layered or continuous')
    if color_by not in {'region', 'layer', 'age', 'cohort', 'source'}:
        raise ValueError('Unsupported pca-color-by')
    palette=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9','#6F4C9B']
    def colors(values):
        result={}
        for i,value in enumerate(sorted(set(values))):
            rgb=colorsys.hsv_to_rgb((i*.381966)%1,.65,.72)
            result[value]=palette[i] if i<len(palette) else '#'+''.join(f'{round(c*255):02x}' for c in rgb)
        return result
    counts=Counter(-1 if p['cohort']=='modern' else int(p['age_bp']//layer_years) for p in points)
    layers=[dict(index=k,start_bp=0 if k==-1 else k*layer_years,
                 end_bp=0 if k==-1 else (k+1)*layer_years,
                 display_bp=0 if k==-1 else (k+.5)*layer_years,n=n,
                 label='Modern' if k==-1 else f'{k*layer_years:,}–{(k+1)*layer_years:,} BP')
            for k,n in sorted(counts.items())]
    return dict(layer_years=layer_years,time_view=time_view,color_by=color_by,layers=layers,
                region_colors=colors(p['region'] for p in points),
                source_colors=colors(p['source'] for p in points),
                bin_rule='Ancient: [start, end), shown at bin midpoint; modern: separate plane at 0 BP. Original ages retained for hover and filtering.')

def pca_centroids(points, bin_years=2000):
    """Keep modern samples out of every ancient centroid and trajectory."""
    import numpy as np
    groups = defaultdict(list)
    for point in points:
        cohort = point['cohort']
        epoch = -1 if cohort == 'modern' else int(point['age_bp'] // bin_years)
        groups[(point['region'], cohort, epoch)].append(point)
    rows = []
    for (region, cohort, epoch), values in sorted(groups.items()):
        x = np.asarray([p['pc1'] for p in values]); y = np.asarray([p['pc2'] for p in values]); n = len(values)
        rows.append(dict(region=region, cohort=cohort, start_bp=0 if epoch < 0 else epoch*bin_years,
            end_bp=0 if epoch < 0 else (epoch+1)*bin_years, n=n, pc1=float(x.mean()), pc2=float(y.mean()),
            pc1_se=float(x.std(ddof=1)/math.sqrt(n)) if n >= 2 else 'NA',
            pc2_se=float(y.std(ddof=1)/math.sqrt(n)) if n >= 2 else 'NA'))
    return rows


def static_pca(payload,directory,candidate_path=None):
    import numpy as np
    plt=_plot(); d=Path(directory); d.mkdir(parents=True,exist_ok=True)
    points=payload['points']; x=np.array([r['pc1'] for r in points]); y=np.array([r['pc2'] for r in points]); age=np.array([r['age_bp'] for r in points])
    display=payload.get('display') or pca_display(points)
    width=display['layer_years']; layers=display['layers']; layered=display['time_view']=='layered'
    regions=sorted({r['region'] for r in points}); colors=display['region_colors']
    bins=np.array([-1 if p['cohort']=='modern' else int(p['age_bp']//width) for p in points])
    z=np.where(bins<0,0,(bins+.5)*width)/1000 if layered else age/1000
    def bounds(values):
        lo,hi=float(min(values)),float(max(values)); pad=max((hi-lo)*.06,abs(lo)*.01,.00001)
        return lo-pad,hi+pad
    xr,yr=bounds(x),bounds(y)
    zmax=max(width,max(age,default=0),max((r['end_bp'] for r in layers),default=width))/1000
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    def decorate(ax, planes=True):
        ax.set(xlabel='PC1 (shared score)',ylabel='PC2 (shared score)',
               zlabel='Time layer (kyr BP)' if layered else 'Age (kyr BP)',xlim=xr,ylim=yr,zlim=(-.03*zmax,1.04*zmax))
        ax.set_box_aspect((1.15,1,1.65)); ax.view_init(elev=17,azim=-58)
        for axis in (ax.xaxis,ax.yaxis,ax.zaxis):
            axis.pane.fill=False
            axis._axinfo['grid']['color']=(.7,.76,.8,.18)
        if planes and layered:
            for layer in layers:
                zz=layer['display_bp']/1000
                verts=[(xr[0],yr[0],zz),(xr[1],yr[0],zz),(xr[1],yr[1],zz),(xr[0],yr[1],zz)]
                plane=Poly3DCollection([verts],facecolors='#75AFC6',edgecolors='#6A91A3',linewidths=.4,alpha=.07)
                ax.add_collection3d(plane)
            if len(layers)<=25:
                ax.set_zticks([r['display_bp']/1000 for r in layers],
                    ['Modern' if r['index']==-1 else f"{r['start_bp']/1000:g}–{r['end_bp']/1000:g}" for r in layers],fontsize=7)
        ax.tick_params(axis='x',labelsize=8);ax.tick_params(axis='y',labelsize=8)
    fig=plt.figure(figsize=(12,12)); ax=fig.add_subplot(111,projection='3d')
    mode=display['color_by']
    numeric=mode in {'age','layer'}
    numeric_values=age if mode=='age' else np.where(bins<0,0,(bins+.5)*width)
    for cohort,marker in [('ancient','o'),('modern','D')]:
        categories=regions if mode=='region' else sorted(display['source_colors']) if mode=='source' else [cohort]
        for category in categories:
            idx=[i for i,p in enumerate(points) if p['cohort']==cohort and (mode not in {'region','source'} or p[mode]==category)]
            if not idx:continue
            c=numeric_values[idx] if numeric else colors[category] if mode=='region' else display['source_colors'][category] if mode=='source' else '#0072B2' if cohort=='ancient' else '#D55E00'
            scat=ax.scatter(x[idx],y[idx],z[idx],s=5 if cohort=='ancient' else 12,alpha=.75,marker=marker,
                c=c,depthshade=False,linewidths=0,label=f'{category}'+(' · modern' if cohort=='modern' and mode in {'region','source'} else ''),
                **(dict(cmap='viridis',vmin=0,vmax=zmax*1000) if numeric else {}))
    decorate(ax)
    if numeric:fig.colorbar(scat,ax=ax,shrink=.5,pad=.08,label='Age (BP)' if mode=='age' else 'Layer midpoint (BP)')
    elif len(ax.get_legend_handles_labels()[1])<=20:ax.legend(fontsize=8,loc='upper left',bbox_to_anchor=(.9,1),frameon=False)
    fig.suptitle(f'Ancient DNA · shared PCA through time\n{width:,}-year layers' if layered else 'Ancient DNA · shared PCA through time',fontsize=17,y=.94)
    fig.text(.12,.07,f"{len(points):,} individuals · color: {mode} · modern: diamonds\n"+
             ('Ancient points shown at bin midpoints; original dates retained in HTML.\n' if layered else '')+
             'One shared PCA basis; older samples are higher.',fontsize=10,color='#496174')
    if payload.get('simulated'):fig.text(.12,.04,'SIMULATED DATA — not observed ancient DNA',color='#aa3300',weight='bold')
    outputs=_save(fig,d/'pca_time_overview')
    cols=4; nrows=math.ceil(len(layers)/cols)
    fig,axes=plt.subplots(nrows,cols,figsize=(14,3.2*nrows),squeeze=False,sharex=True,sharey=True)
    for j,ax in enumerate(axes.flat):
        if j>=len(layers): ax.set_visible(False); continue
        layer=layers[j]; mask=bins==layer['index']
        ax.scatter(x[mask],y[mask],s=5,c=[colors[p['region']] for i,p in enumerate(points) if mask[i]],alpha=.65,linewidths=0)
        ax.set(title=f"{layer['label']}; n={sum(mask):,}",xlabel='PC1',ylabel='PC2',xlim=xr,ylim=yr)
        ax.spines[['top','right']].set_visible(False)
    fig.suptitle(f'Fixed {width:,}-year slices · shared axes · color: region · modern separate',y=1.01)
    fig.tight_layout();outputs+=_save(fig,d/'pca_time_slices')
    h,w=tsv_writer(d/'pca_time_layers.tsv',['index','start_bp','end_bp','display_bp','n','label'])
    with h:w.writerows(layers)
    centroid=pca_centroids(points)
    fig,ax=plt.subplots(figsize=(9,7))
    for region in regions:
        trajectory=[]
        for row in centroid:
            if row['region']!=region: continue
            ancient=row['cohort']=='ancient'
            if ancient: trajectory.append((row['start_bp'],row['pc1'],row['pc2']))
            ax.errorbar(row['pc1'],row['pc2'],xerr=1.96*row['pc1_se'] if row['pc1_se']!='NA' else None,
                yerr=1.96*row['pc2_se'] if row['pc2_se']!='NA' else None,color=colors[region],alpha=.55,
                fmt='o' if ancient else 'D',ms=3 if ancient else 6)
        trajectory.sort(reverse=True)
        ax.plot([t[1] for t in trajectory],[t[2] for t in trajectory],color=colors[region],label=region,alpha=.75)
    ax.set(xlabel='PC1',ylabel='PC2',title='Ancient regional centroids; modern diamonds separate\nDescriptive 95% mean intervals')
    if len(regions)<=16: ax.legend(fontsize=7)
    outputs+=_save(fig,d/'pca_time_centroids')
    if centroid:
        h,w=tsv_writer(d/'pca_time_centroids.tsv',list(centroid[0]))
        with h: w.writerows(centroid)
    variants=payload['variants']; by_key={(r['chrom'],int(r['pos']),r['ref'],r['alt']):r for r in variants}
    loci=[]
    if candidate_path and Path(candidate_path).exists():
        seen=set(); complement=str.maketrans('ACGT','TGCA')
        for r in figure_rows(candidate_path):
            if r['site_id'] in seen: continue
            key=(r['chrom37'],int(r['pos37']),r['ref37'],r['alt37'])
            if key not in by_key: continue
            seen.add(r['site_id']); anchor=by_key[key]; extra=r['extra_alt']
            if r['strand37to38']=='-': extra=extra.translate(complement)
            extra=by_key.get((*key[:3],extra))
            vals=anchor['values']
            increasing=None if vals is None else [None if v is None else v if float(r['s'])>0 else 2-v for v in vals]
            decreasing=None if increasing is None else [None if v is None else 2-v for v in increasing]
            loci.append((r['site_id'],[('Anchor increasing',increasing),('Anchor decreasing',decreasing),('Extra allele',extra['values'] if extra else None)]))
    if not loci:
        loci=[(r['variant_id'],[(f'ALT {r["alt"]} copies',r['values'])]) for r in variants]
    for name,panels in loci:
        fig=plt.figure(figsize=(6*len(panels),6))
        for j,(label,vals) in enumerate(panels,1):
            ax=fig.add_subplot(1,len(panels),j,projection='3d')
            ax.scatter(x,y,z,c='#bcbcbc',s=4,alpha=.3)
            if vals is not None:
                idx=[i for i,v in enumerate(vals) if v is not None]
                scat=ax.scatter(x[idx],y[idx],z[idx],c=[vals[i] for i in idx],vmin=0,vmax=2,cmap='viridis',s=8)
                fig.colorbar(scat,ax=ax,shrink=.5,label='Copies of the named allele')
            ax.set(title=label+(' (unavailable)' if vals is None else ''),xlabel='PC1',ylabel='PC2',zlabel='kyr BP')
            decorate(ax)
        fig.suptitle(name+'; grey = unavailable, 0/1/2 = copies')
        outputs+=_save(fig,d/('pca_time_locus_'+re.sub(r'[^A-Za-z0-9_.-]','_',name)),formats=('png',))
    return dict(outputs=outputs,centroid_rows=len(centroid),locus_figures=len(loci),slice_width_bp=width,
                centroid_bin_years=2000,time_view=display['time_view'],color_by=mode,layer_table=str(d/'pca_time_layers.tsv'))

def reuse_metag_census(path,out):
    data=json.loads(Path(path).read_text())
    if data.get('assembly')!='hg38' or not data.get('inputs',{}).get('repeats_sha256'):
        raise ValueError('Expected versioned hg38 metag census with input checksums')
    inputs=data['inputs']
    for key in ('genome','repeats'):
        if not Path(inputs[key]).is_file(): raise ValueError('Shared metag census source missing: '+inputs[key])
    # Retain metag's exact class policy; it differs from allele3's ERV-related + Gypsy.
    dump_json(out,dict(source_report=fingerprint(path),source_inputs={k:fingerprint(inputs[k]) for k in ('genome','repeats')},
        source_sha256={k:inputs[k+'_sha256'] for k in ('genome','repeats')},report=data,
        reuse='read-only metag result; no index/genome download; no reassignment of its family policy',
        allele3_policy=dict(erv_core=['ERV1','ERVK','ERVL'],erv_related=['ERV1','ERVK','ERVL','ERVL-MaLR','Gypsy'],ltr_all='all LTR including solo LTR')))
    return data

def write_local_review(modern_dir,evolution_dir,metag_census=None,make_plots=True):
    out=Path(evolution_dir); summary=Path(modern_dir)/'summary'; figures=out/'figures'; figures.mkdir(parents=True,exist_ok=True)
    inventory=[]; plt=_plot() if make_plots else None
    def bar(number,title,labels,values,ylabel):
        if not make_plots: inventory.append(dict(figure=number,status='plots_disabled',outputs=[])); return
        if not values: inventory.append(dict(figure=number,status='unavailable_data',outputs=[])); return
        fig,ax=plt.subplots(figsize=(max(7,min(15,len(labels)*.5)),5)); ax.bar(range(len(values)),values)
        ax.set_xticks(range(len(labels)),labels,rotation=60,ha='right',fontsize=8); ax.set(title=title,ylabel=ylabel)
        inventory.append(dict(figure=number,status='measured',outputs=_save(fig,figures/f'Fig{number}')))
    spectrum=summary/'architecture.global_strata.tsv.gz'
    if spectrum.exists() and make_plots:
        counts=defaultdict(lambda:[0,0,0])
        for r in figure_rows(spectrum):
            if r['sensitivity']=='all' and r['population_mode']=='all':
                for i,k in enumerate(('n_sites','n_multi','n_three')): counts[int(r['threshold'])][i]+=int(r[k])
        tiers=sorted(counts);fig,axes=plt.subplots(1,2,figsize=(11,4))
        base=counts.get(1,[0,0,0]);axes[0].bar(['1 ALT','2 ALTs','3 ALTs'],[base[0]-base[1],base[1]-base[2],base[2]])
        axes[0].set(title='Modern eligible SNV loci',ylabel='Loci',yscale='log')
        for j,label in ((1,'Multi-ALT'),(2,'Three ALT')):
            axes[1].plot(tiers,[counts[t][j] for t in tiers],marker='o',label=label)
        axes[1].set(title='AC threshold persistence; fixed eligible denominator',xlabel='Minimum AC per qualifying ALT',ylabel='Loci');axes[1].legend()
        inventory.append(dict(figure=1,status='measured',outputs=_save(fig,figures/'Fig1')))
    else: bar(1,'',[],[],'')
    ld=summary/'ld_block_ranked.tsv'
    rs=[r for r in figure_rows(ld) if r['ancestry']=='EUR' and r['threshold']=='1' and r['sensitivity']=='all' and r['population_mode']=='all'][:20] if ld.exists() else []
    bar(2,'EUR LD blocks: top descriptive O/E', [r['block_id'] for r in rs],[float(r['oe']) for r in rs],'Observed / LOCO expected')
    consensus=summary/'ld_hotspot_consensus.tsv'
    if consensus.exists() and make_plots:
        overlaps=list(figure_rows(consensus))
        if overlaps:
            fig,ax=plt.subplots(figsize=(12,5));chroms=sorted({r['chrom'] for r in overlaps},key=lambda c:(not c.isdigit(),int(c) if c.isdigit() else c));ys={c:i for i,c in enumerate(chroms)}
            for r in overlaps:
                ax.plot([int(r['start'])/1e6,int(r['end'])/1e6],[ys[r['chrom']]]*2,lw=4,alpha=.8,color=plt.get_cmap('viridis')((int(r['n_ancestries'])-2)/2))
            ax.set_yticks(range(len(chroms)),chroms);ax.set(xlabel='GRCh38 position (Mb)',ylabel='Chromosome',title='Physical overlap of ancestry-specific top 1% LD blocks; 2–4 supporting ancestries')
            inventory[-1]['outputs']+=_save(fig,figures/'Fig2_consensus')
    erv=summary/'erv_enrichment.tsv'; totals=defaultdict(lambda:[0.,0.,0])
    if erv.exists():
        for r in figure_rows(erv):
            if r['threshold']!='1' or r['sensitivity']!='all' or r['population_mode']!='all': continue
            label=r['label']
            if label.endswith(('upstream','downstream')): continue
            z=totals[label]; z[0]+=int(r['observed_multi'])
            if r['expected_multi']=='NA': z[2]+=1
            else: z[1]+=float(r['expected_multi'])
    valid={k:o/e for k,(o,e,miss) in totals.items() if e and not miss}
    if valid and make_plots:
        fig,axes=plt.subplots(2,2,figsize=(14,10))
        for ax,labels,title in [(axes[0,0],[k for k in valid if k.startswith('inside_') or k.startswith('negative_control:')],'Bodies and negative controls'),
                                (axes[0,1],sorted([k for k in valid if k.startswith('family:')],key=lambda k:-valid[k])[:15],'LTR families')]:
            ax.bar(range(len(labels)),[valid[k] for k in labels]);ax.set_xticks(range(len(labels)),labels,rotation=60,ha='right',fontsize=7);ax.axhline(1,color='grey',ls='--');ax.set(title=title,ylabel='Observed / expected')
        bins=['0-100','100-200','200-500','500-1000','1000-2000','2000-5000','5000-10000']
        for ax,bt in ((axes[1,0],'5p'),(axes[1,1],'3p')):
            for tier in ('erv_core','erv_related','ltr_all'):
                ax.plot(range(7),[valid.get(f'{tier}:{bt}_boundary:{b}',float('nan')) for b in bins],marker='o',label=tier)
            ax.axhline(1,color='grey',ls='--');ax.set_xticks(range(7),bins,rotation=35,ha='right');ax.set(title=bt+' observed joined boundary',xlabel='Distance (bp)',ylabel='Observed / expected');ax.legend(fontsize=8)
        fig.tight_layout();inventory.append(dict(figure=3,status='measured',outputs=_save(fig,figures/'Fig3')))
    else: bar(3,'',[],[],'')
    pops=summary/'population_classes.tsv'; counts=Counter()
    if pops.exists():
        for r in figure_rows(pops):
            if r['threshold']=='1': counts[r['classification']]+=int(r['count'])
    bar(4,'Population partitioning of multi-ALT loci',list(counts),list(counts.values()),'Loci')
    candidates=out/'universes'/'U2_all_strong.tsv'; cases=list(figure_rows(candidates)) if candidates.exists() else []
    if cases and make_plots:
        fig,ax=plt.subplots(figsize=(8,5)); ax.scatter([float(r['s']) for r in cases],[float(r['extra_af']) for r in cases],s=10,alpha=.5)
        ax.set(xlabel='Published anchor selection coefficient',ylabel='Modern extra ALT AF',yscale='log',title='All eligible strong contrasts; extra ALT selection is unestimated')
        inventory.append(dict(figure=5,status='measured',outputs=_save(fig,figures/'Fig5')))
        trajectory=out/'trajectories'/'trajectory_bins.tsv'
        if trajectory.exists():
            chosen=[];seen=set();series=defaultdict(list);keys=set();complement=str.maketrans('ACGT','TGCA')
            for r in cases:
                if r['site_id'] in seen:continue
                seen.add(r['site_id']);chosen.append(r)
                extra=r['extra_alt'].translate(complement) if r['strand37to38']=='-' else r['extra_alt']
                keys.update([(r['chrom37'],r['pos37'],r['ref37'],a) for a in (r['alt37'],extra)])
                if len(chosen)==6:break
            for r in figure_rows(trajectory):
                key=(r['chrom'],r['pos'],r['ref'],r['alt'])
                if key in keys and r['cohort']=='ancient' and r['group_type']=='raw_all' and r['bin_status']=='ok':
                    try:af=float(r['af']);age=(float(r['start_bp'])+float(r['end_bp']))/2000
                    except ValueError:continue
                    if math.isfinite(af) and math.isfinite(age):series[key].append((age,af))
            fig,axes=plt.subplots(2,3,figsize=(14,8),squeeze=False)
            max_age=max([t for values in series.values() for t,f in values],default=10)
            for ax,r in zip(axes.flat,chosen):
                anchor=(r['chrom37'],r['pos37'],r['ref37'],r['alt37']);extra=r['extra_alt'].translate(complement) if r['strand37to38']=='-' else r['extra_alt']
                a=sorted(series[anchor]);e=sorted(series[(*anchor[:3],extra)])
                if a:
                    inc=[f if float(r['s'])>0 else 1-f for t,f in a]
                    ax.plot([t for t,f in a],inc,'o-',label='Anchor increasing')
                    ax.plot([t for t,f in a],[1-f for f in inc],'o--',label='Anchor decreasing')
                else:ax.text(.02,.14,'Anchor: unavailable',color='grey',transform=ax.transAxes)
                if e:ax.plot(*zip(*e),'s-',label='Extra allele')
                else:ax.text(.02,.04,'Extra: unavailable',color='grey',transform=ax.transAxes)
                ax.set(title=r['site_id'],xlabel='Age (kyr BP)',ylabel='Observed allele frequency',ylim=(0,1),xlim=(max_age+.5,0))
                if a or e:ax.legend(fontsize=7)
            for ax in list(axes.flat)[len(chosen):]:ax.set_visible(False)
            fig.suptitle('Ancient strong-anchor trajectories: raw frequencies; cohort composition is unadjusted')
            fig.tight_layout();inventory[-1]['outputs']+=_save(fig,figures/'Fig5_trajectories')
    else: inventory.append(dict(figure=5,status='unavailable_data' if not cases else 'plots_disabled',outputs=[]))
    pca=out/'pca'/'pca_time.html.json'; pcam=json.loads(pca.read_text()) if pca.exists() else {}
    inventory.append(dict(figure=6,status='measured' if pcam.get('static_figures') else 'unavailable_PCA',outputs=pcam.get('static_figures',{}).get('outputs',[])))
    ukb=out/'ukb'/'U0_mfi.tsv'; points=[]; n=0
    if ukb.exists():
        import random
        rng=random.Random(20261022)
        for r in figure_rows(ukb):
            if r['ukb_status']!='present' or r['gnomad_nfe_af']=='NA' or r['ukb_alt_af']=='NA' or r['ukb_info']=='NA' or float(r['ukb_info'])<.8: continue
            n+=1; pt=(float(r['gnomad_nfe_af']),float(r['ukb_alt_af']))
            if len(points)<10000: points.append(pt)
            else:
                j=rng.randrange(n)
                if j<len(points): points[j]=pt
    if points and make_plots:
        fig,ax=plt.subplots(figsize=(6,6)); ax.scatter(*zip(*points),s=3,alpha=.3); ax.plot([0,1],[0,1],'--',color='grey')
        ax.set(xlabel='gnomAD NFE AF',ylabel='UKB ALT AF',title=f'Exact GRCh37 alleles; INFO >= 0.8; n={n:,}')
        inventory.append(dict(figure=7,status='measured',outputs=_save(fig,figures/'Fig7')))
    else: inventory.append(dict(figure=7,status='unavailable_exact_frequency_pairs' if not points else 'plots_disabled',outputs=[]))
    geometry=out/'geometry'/'geometry.tsv'
    if geometry.exists() and make_plots:
        panels=defaultdict(list)
        for r in figure_rows(geometry):
            try: val=float(r['atlas_signed_cosine_toward_estimated_increasing'])
            except (ValueError,TypeError): continue
            if r['atlas_status']=='ok' and math.isfinite(val): panels[r['atlas_scorer']].append(val)
        if panels:
            fig,ax=plt.subplots(figsize=(10,5));ax.boxplot(list(panels.values()),tick_labels=list(panels));ax.tick_params(axis='x',rotation=35)
            ax.set(ylabel='Signed cosine toward increasing anchor',title='U3 functional geometry; comparable feature axes within scorer')
            inventory.append(dict(figure=8,status='measured',outputs=_save(fig,figures/'Fig8')))
        else: inventory.append(dict(figure=8,status='unavailable_valid_geometry',outputs=[]))
    else: inventory.append(dict(figure=8,status='plots_disabled' if geometry.exists() else 'pending_U3_AlphaGenome; no fabricated functional results',outputs=[]))
    if make_plots:
        fig,ax=plt.subplots(figsize=(13,3)); ax.axis('off')
        labels=['Modern observation','AC persistence','Population partition','LD / ERV architecture','Ancient anchor','UKB replication','U3 function']
        for i,label in enumerate(labels):
            x=.07+i*.145; ax.text(x,.55,label.replace(' ','\n',1),ha='center',va='center',bbox=dict(boxstyle='round',facecolor='#e5edf4'),transform=ax.transAxes)
            if i<6: ax.annotate('',xy=(x+.102,.55),xytext=(x+.055,.55),xycoords='axes fraction',arrowprops=dict(arrowstyle='->'))
        ax.text(.5,.1,'Evidence workflow, not inferred allele age or a causal life cycle',ha='center',transform=ax.transAxes)
        inventory.append(dict(figure=9,status='schematic',outputs=_save(fig,figures/'Fig9')))
    if metag_census: reuse_metag_census(metag_census,out/'shared_metag_census.json')
    evidence=integrated_evidence(modern_dir,evolution_dir)
    result=dict(figures=inventory,integrated_evidence=evidence,review_inputs=dict(modern=str(summary),evolution=str(out)),
        interpretation='Unavailable layers remain NA. LD/ERV are descriptive conditional enrichment. Missing MFI is not extinction. AlphaGenome is selected only after local results.')
    dump_json(figures/'figure_manifest.json',result)
    (out/'LOCAL_REVIEW.md').write_text('# FINAL local results review\n\n'+ '\n'.join(f'- Fig {r["figure"]}: {r["status"]}' for r in inventory)+'\n\nNo thresholds were relaxed to produce a signal. Inspect LOCO coverage, unmatched U2 controls, exact UKB orientation, annotation availability and sensitivity results before biological interpretation.\n')
    return result

def integrated_evidence(modern_dir,evolution_dir):
    """Transparent evidence completeness/criteria, never a causal or age score."""
    from scan import BlockIndex, read_blocks
    from scan import parse_joined_rmsk, ERVOverlapSweep
    out=Path(evolution_dir);candidate=out/'universes'/'candidates.tsv'
    if not candidate.exists(): return dict(status='unavailable_candidates')
    data=[];lookup=defaultdict(list)
    for r in figure_rows(candidate):
        x=dict(site_id=r['site_id'],selection_id=r['selection_id'],extra_alt=r['extra_alt'],chrom=r['chrom'],pos=int(r['pos']),
            group=r['group'],persistence_ac_ge5=int(float(r['extra_ac'])>=5),anchor_strong=int(float(r['posterior'])>=.99),
            ukb_exact_info_ge08='NA',ld_block38='NA',ld_oe='NA',inside_erv_core='NA',functional_geometry_measured='NA')
        data.append(x);lookup[(r['site_id'],r['selection_id'],r['extra_alt'])].append(x)
    ukb=out/'ukb'/'U0_mfi.tsv'
    if ukb.exists():
        for r in figure_rows(ukb):
            if r['target_role']!='extra': continue
            for x in lookup.get((r.get('site_id'),r.get('selection_id'),r.get('extra_alt')),[]):
                if r['ukb_status']=='unavailable_file' or (r['ukb_status']=='present' and r['ukb_info']=='NA'): continue
                x['ukb_exact_info_ge08']=int(r['ukb_status']=='present' and float(r.get('ukb_info') or 0)>=.8) if r['ukb_status']=='present' else 0
    run=Path(modern_dir)/'run.json';options=json.loads(run.read_text()).get('options',{}) if run.exists() else {}
    if options.get('ld_block_dir'):
        blocks,_=read_blocks(Path(options['ld_block_dir'])/'EUR.38.bed','EUR','38');idx=BlockIndex(blocks)
        en=Path(modern_dir)/'summary'/'ld_block_enrichment.tsv';oe={}
        if en.exists():
            oe={r['block_id']:r['oe'] for r in figure_rows(en) if r['ancestry']=='EUR' and r['sensitivity']=='all' and r['threshold']=='1' and r['population_mode']=='all'}
        for x in data:
            b=idx.locate(x['chrom'],x['pos'])
            if b: x['ld_block38']=b.block_id;x['ld_oe']=oe.get(b.block_id,'NA')
            else:x['ld_block38']='outside_partition'
    annotation_cache=options.get('architecture_annotation_cache')
    if annotation_cache:
        grouped=defaultdict(list)
        for x in data:grouped[x['chrom']].append(x)
        for c,group in grouped.items():
            path=Path(annotation_cache)/f'chr{c}.tsv.gz'
            if not path.exists():continue
            sweep=ERVOverlapSweep(parse_joined_rmsk(path,chromosome=c))
            for x in sorted(group,key=lambda r:r['pos']):x['inside_erv_core']=int(any(e.tier=='erv_core' for e in sweep.locate(c,x['pos'])))
    geo=out/'geometry'/'geometry.tsv'
    if geo.exists():
        measured={r['site_id'] for r in figure_rows(geo) if r['atlas_status']=='ok'}
        for x in data:
            if x['site_id'] in measured:x['functional_geometry_measured']=1
    components=['persistence_ac_ge5','anchor_strong','ukb_exact_info_ge08','functional_geometry_measured']
    for x in data:
        available=[x[k] for k in components if x[k]!='NA']
        x['criteria_observed']=len(available);x['criteria_positive']=sum(available)
        x['descriptive_evidence_fraction']=sum(available)/len(available) if available else 'NA'
    dest=out/'figures'/'integrated_evidence.tsv'
    h,w=tsv_writer(dest,list(data[0]) if data else ['status'])
    with h:w.writerows(data)
    result=dict(status='measured',rows=len(data),output=str(dest),components=components,
        interpretation='Equal-weight descriptive criteria fraction among observed components; not calibrated probability, allele age, selection on extra allele or causal evidence. LD/ERV context stays separate from directional evidence.',
        missing='NA; both numerator and denominator report observed layers')
    dump_json(out/'figures'/'integrated_evidence.json',result);return result


# ========================================================================
# Evolve report
# ========================================================================

# Exploratory, matched molecular comparisons and a local evolution report.


from collections import Counter, defaultdict
import html
import json
import math
from pathlib import Path

from evolve import number, table_rows
from allele_common import dump_json, fingerprint, tsv_writer


def _bh(values):
    """Benjamini-Hochberg within the explicitly supplied family."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    result, bound = [None] * len(values), 1.
    for rank in range(len(order), 0, -1):
        i = order[rank - 1]
        bound = min(bound, values[i] * len(order) / rank)
        result[i] = bound
    return result


def compare_geometry(geometry_path, out_dir, *, draws=10000, seed=20261021, min_pairs=10):
    """Within-scorer normalized distance contrast, paired sign flips and bootstrap.

    This is observational, conditional on candidate selection. Physical block
    thinning does not guarantee LD independence; uncertainty is exploratory.
    Never regard weak selection evidence as proof of neutrality.
    """
    import numpy as np
    if draws < 100 or min_pairs < 2:
        raise ValueError("At least 100 resampling draws and min_pairs >= 2 are required")
    out = Path(out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    grouped, qc = defaultdict(dict), Counter()
    fields = ["site_id", "group", "match_set", "scorer", "status", "closeness_to_increasing",
              "n_features", "axis_sha256", "locus_block"]
    handle, writer = tsv_writer(out / "site_closeness.tsv", fields)
    with handle:
        for row in table_rows(geometry_path, ("site_id", "group", "match_set", "atlas_scorer", "atlas_status")):
            qc[row["atlas_status"]] += 1
            up = number(row.get("atlas_distance_to_estimated_increasing"))
            down = number(row.get("atlas_distance_to_estimated_decreasing"))
            valid = row["atlas_status"] == "ok" and up is not None and down is not None and up + down > 0
            score = (down - up) / (down + up) if valid else None
            writer.writerow(dict(site_id=row["site_id"], group=row["group"], match_set=row["match_set"],
                                 scorer=row["atlas_scorer"], status=row["atlas_status"],
                                 closeness_to_increasing=score if score is not None else ".",
                                 n_features=row.get("atlas_n_features", "."), axis_sha256=row.get("atlas_axis_sha256", "."),
                                 locus_block=row.get("locus_block", ".")))
            if not valid or row["match_set"] in ("", ".") or row["group"] not in ("strong", "weak_signal"):
                continue
            key = (row["atlas_scorer"], row["match_set"])
            bucket=grouped[key].setdefault(row['group'],[])
            if any(x[2]==row['site_id'] for x in bucket) or (row['group']=='strong' and bucket):
                raise ValueError(f"Repeated case/site in match set/scorer: {key}")
            bucket.append((score,row.get('atlas_n_features'),row['site_id'],row.get('atlas_axis_sha256'),row.get('locus_block','.')))
    by_scorer = defaultdict(list)
    handle, writer = tsv_writer(out / "paired_closeness.tsv", ["scorer", "axis_sha256", "match_set", "strong_site", "weak_site", "strong_closeness", "weak_closeness", "difference", "status", "n_controls"])
    with handle:
        for (scorer, pair), groups in sorted(grouped.items()):
            if set(groups) != {"strong", "weak_signal"}:
                qc["incomplete_geometry_pairs"] += 1
                continue
            case, controls = groups['strong'][0], groups['weak_signal']
            status=('feature_count_mismatch' if any(case[1]!=c[1] for c in controls) else
                    'feature_axis_mismatch' if not case[3] or case[3] in ('.','NA') or any(case[3]!=c[3] for c in controls) else 'ok')
            control=(sum(c[0] for c in controls)/len(controls),case[1],','.join(c[2] for c in controls),case[3])
            difference = case[0] - control[0]
            writer.writerow(dict(scorer=scorer, axis_sha256=case[3] or ".", match_set=pair, strong_site=case[2], weak_site=control[2],
                                 strong_closeness=case[0], weak_closeness=control[0], difference=difference, status=status,n_controls=len(controls)))
            if status == "ok":
                by_scorer[(scorer, case[3])].append((difference,case[4] if case[4] not in (".","NA","") else pair))
            else:
                qc[status] += 1
    results, tested = [], []
    rng = np.random.default_rng(seed)
    for (scorer, axis), values in sorted(by_scorer.items()):
        arr = np.asarray([x[0] for x in values], dtype=float)
        clusters=defaultdict(list)
        for value,block in values: clusters[block].append(value)
        groups=list(clusters.values())
        row = dict(scorer=scorer, axis_sha256=axis, n_pairs=len(arr),n_clusters=len(groups), mean_difference=float(arr.mean()),
                   bootstrap_ci_low=".", bootstrap_ci_high=".", sign_flip_p=".", bh_q=".",
                   status="too_few_pairs", inference="exploratory_matched_observational")
        if len(arr) >= min_pairs and len(groups)>=2:
            means, permuted = [], []
            for start in range(0, draws, 100):
                size = min(100, draws - start)
                totals=np.asarray([sum(g) for g in groups]);counts=np.asarray([len(g) for g in groups])
                chosen=rng.integers(0,len(groups),size=(size,len(groups)))
                means.extend((totals[chosen].sum(axis=1)/counts[chosen].sum(axis=1)).tolist())
                signs=rng.choice((-1.,1.),size=(size,len(groups)))
                permuted.extend(((signs*totals).sum(axis=1)/len(arr)).tolist())
            p = (1 + sum(abs(x) >= abs(row["mean_difference"]) for x in permuted)) / (draws + 1)
            lo, hi = np.quantile(means, [.025, .975])
            row.update(bootstrap_ci_low=float(lo), bootstrap_ci_high=float(hi), sign_flip_p=p, status="ok")
            tested.append((len(results), p))
        results.append(row)
    for (index, _), q in zip(tested, _bh([p for _, p in tested])):
        results[index]["bh_q"] = q
    fields = ["scorer", "axis_sha256", "n_pairs", "n_clusters", "mean_difference", "bootstrap_ci_low", "bootstrap_ci_high",
              "sign_flip_p", "bh_q", "status", "inference"]
    handle, writer = tsv_writer(out / "matched_functional_summary.tsv", fields)
    with handle:
        writer.writerows(results)
    result = dict(input=fingerprint(geometry_path), draws=draws, seed=seed, min_pairs=min_pairs,
                  metric="(distance_to_decreasing-distance_to_increasing)/(distance_to_decreasing+distance_to_increasing)",
                  qc=dict(qc), resampling="Whole case LD/proximity clusters; each case compared with the mean of all available exact-feature-axis controls",
                  multiple_testing="BH over tested scorer x exact feature-panel comparisons in this output only",
                  limitation="No causal or balancing-selection claim; approximate physical blocks do not establish independent LD loci. "
                             "Exact biological feature axes are matched, but correlated tracks and heterogeneous biological effects remain.",
                  outputs=[str(out / name) for name in ("site_closeness.tsv", "paired_closeness.tsv", "matched_functional_summary.tsv")])
    dump_json(out / "comparison_manifest.json", result)
    return result


def write_report(out_dir):
    """Small self-contained HTML index of actual local stage outputs."""
    out = Path(out_dir).resolve()
    cards, sections = [], []
    manifests = [("counts_audit", "原始 check 计数审计", "manifest.json"),
                 ("metadata", "样本信息", "metadata_qc.json"), ("join", "现代—古代选择统计连接", "manifest.json"),
                 ("candidates", "候选与匹配", "manifest.json"), ("universes", "U0–U3 独立分析集", "manifest.json"),
                 ("ukb", "UKB MFI 精确匹配", "U0_mfi.tsv.manifest.json"), ("panel", "古代面板覆盖", "manifest.json"),
                 ("trajectories", "频率轨迹", "trajectory_manifest.json"),
                 ("atlas", "稀疏功能查询", "query.json"), ("geometry", "分子效应几何", "geometry.json")]
    for folder, title, filename in manifests:
        path = out / folder / filename
        if not path.exists():
            continue
        data = json.loads(path.read_text())
        values = data.get("qc", data.get("target_counts", data.get("counts", data.get("status_counts", {}))))
        if values:
            rows = "".join(f"<tr><th>{html.escape(str(k))}</th><td>{html.escape(str(v))}</td></tr>" for k, v in values.items())
            cards.append(f"<section><h2>{title}</h2><table>{rows}</table></section>")
    for relative, title in [("figures/figure_manifest.json", "Fig1–Fig9 实际产物与缺失状态"),
                             ("figures/integrated_evidence.tsv", "观测证据分量与描述性分数"),
                             ("universes/U2_all_strong.tsv", "全部合格强信号"),
                             ("universes/U2_clustered_bootstrap.tsv", "按 LD block 聚类的 bootstrap"),
                             ("shared_metag_census.json", "共享 metag 全基因组 ERV/LTR 覆盖率"),
                             ("pca/pca_time.html", "PC1–PC2–时间：交互式三维图"),
                             ("counts_audit/counts_review.html", "旧 check 计数：按染色体的描述性汇总"),
                             ("trajectories/trajectory_bins.tsv", "按年代与地区分层的轨迹表"),
                             ("trajectories/trajectory_standardized.tsv", "固定构成标准化轨迹"),
                             ("panel/site_variant_links.tsv", "古代面板是否覆盖额外碱基"),
                             ("candidates/candidates.tsv", "候选位点与匹配对照"),
                             ("comparison/matched_functional_summary.tsv", "探索性功能比较"),
                             ("pca/commands.sh", "PCA 计算脚本")]:
        if (out / relative).exists():
            sections.append(f'<li><a href="{relative}">{title}</a></li>')
    text = """<!doctype html><html lang="zh"><meta charset="utf-8"><title>allele3 evolution</title>
<style>body{font:16px system-ui,sans-serif;background:#f4f7fb;color:#172b46;max-width:1150px;margin:32px auto;padding:0 24px}
h1{font-size:32px}h2{font-size:20px}section{background:white;border:1px solid #dde4ef;border-radius:12px;padding:20px;margin:16px 0}
table{border-collapse:collapse;width:100%}th,td{text-align:left;padding:7px;border-bottom:1px solid #edf1f7}a{color:#135fb2}.note{border-left:4px solid #e4aa42;padding:12px 18px;background:#fff8e8}</style>
<h1>allele3：现代分布、古代轨迹与分子效应</h1>
<p>此报告只汇总本地已完成步骤。空缺的 Atlas 或 PCA 结果表示该步骤尚未执行，不代表零效应或没有结构。</p>
<div class="note">REF 不是祖先状态；未列入古代双等位插补面板的第三碱基属于未覆盖。
本流程的选择统计来自已发表的精确等位基因对比，不把它转赋给第三碱基；频率曲线与 PCA 不能单独证明选择。</div>
<h2>结果入口</h2><ul>""" + "".join(sections) + "</ul>" + "".join(cards) + "</html>"
    (out / "report.html").write_text(text, encoding="utf-8")
    return str(out / "report.html")


# ========================================================================
# Report
# ========================================================================

# Six-figure focused report with explicit missing evidence and no invented claims.

import argparse, csv, json, shutil, math
from pathlib import Path


def report_rows(path):
    if not path or not Path(path).exists():return []
    with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))


def make_report(a):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True);figs=out/'figures';figs.mkdir(exist_ok=True)
    population=Path(a.population_dir);modern=Path(a.modern_dir) if a.modern_dir else population.parent/'gnomad_local';summary=modern/'summary'
    evolution=Path(a.evolution_dir) if a.evolution_dir else population.parent/'evolution'
    kg=Path(a.replication_dir) if a.replication_dir else population.parent/'population_1kg'
    pop=report_rows(population/'population_balanced_summary.tsv');effects=report_rows(kg/'effect_direction_agreement.tsv')
    enrichment=report_rows(a.herv_enrichment or summary/'enrichment_summary.tsv')
    if a.archaic_enrichment and Path(a.archaic_enrichment)!=Path(a.herv_enrichment or summary/'enrichment_summary.tsv'):enrichment+=report_rows(a.archaic_enrichment)
    inventory=[]
    def save(fig,number,title):
        paths=[]
        for ext in ('png','pdf'):
            p=figs/f'Fig{number}.{ext}';fig.savefig(p,dpi=180,bbox_inches='tight');paths.append(str(p))
        plt.close(fig);inventory.append(dict(figure=number,title=title,status='measured',outputs=paths))
    # Original QC summary is kept visible; AC1 is explicitly background only.
    qc=report_rows(summary/'occupancy_thresholds.tsv')
    if qc:
        fig,ax=plt.subplots(figsize=(8,4))
        ax.bar([r['threshold'] for r in qc],[int(r['n_multi']) for r in qc],color=['#aaa' if r['threshold'] in {'1','2'} else '#35678c' for r in qc])
        ax.set(xlabel='Per-ALT AC threshold (1/2: QC; 5: primary; 10: strict)',ylabel='Observed loci with ≥2 qualifying ALTs',title='Persistence of multi-ALT occupancy');ax.set_yscale('log')
        save(fig,1,'Persistence and QC')
    else:inventory.append(dict(figure=1,status='unavailable',outputs=[]))
    fig,axes=plt.subplots(1,2,figsize=(13,5));classes=['within_group_coexistence','population_partitioned','unresolved_after_balance'];colors=['#35678c','#d78a39','#aaa']
    for ax,records,key,title in [(axes[0],pop,'classification','gnomAD equal AN'),(axes[1],report_rows(kg/'population_1kg_subsample_summary.tsv'),'classification','1KG repeated equal-sample sensitivity')]:
        combos=sorted({(int(r.get('global_ac_threshold',r.get('threshold',0))),int(r.get('target_an',0))) for r in records})
        bottom=[0.]*len(combos)
        for c,color in zip(classes,colors):
            vals=[]
            for t,n in combos:
                r=next((r for r in records if int(r.get('global_ac_threshold',r.get('threshold',0)))==t and int(r.get('target_an',0))==n and r[key]==c),None)
                vals.append(float(r.get('fraction',r.get('mean',0)) or 0) if r else 0)
            ax.bar(range(len(combos)),vals,bottom=bottom,label=c,color=color);bottom=[a+b for a,b in zip(bottom,vals)]
        ax.set_xticks(range(len(combos)),[f'AC≥{t}\nAN={n}' if n else f'AC≥{t}' for t,n in combos]);ax.set_ylim(0,1);ax.set_title(title)
        if not records:ax.text(.5,.5,'Unavailable',ha='center',transform=ax.transAxes)
    axes[0].legend(fontsize=7);fig.tight_layout();save(fig,2,'Balanced population architecture and individual-level replication')
    regions=report_rows(summary/'ld_hotspot_consensus_regions.tsv')
    if regions:
        fig,ax=plt.subplots(figsize=(10,max(4,len(regions)*.22)))
        for i,r in enumerate(regions):ax.plot([int(r['start'])/1e6,int(r['end'])/1e6],[i,i],lw=4)
        ax.set_yticks(range(len(regions)),[f"chr{r['chrom']} AC≥{r['threshold']}" for r in regions]);ax.set_xlabel('Position (Mb)');save(fig,3,'Merged four-ancestry LD consensus regions')
    else:inventory.append(dict(figure=3,status='no_four_way_consensus',outputs=[]))
    for number,archaic in ((4,False),(5,True)):
        data=[r for r in enrichment if ('::archaic:' in r['label'])==archaic and r['sensitivity']=='all' and r['status']=='complete' and r['ci_low'] not in {'','None','NA'}]
        # Display all main labels; detailed sensitivity/FDR tests stay in the table.
        if data:
            fig,ax=plt.subplots(figsize=(11,max(5,len(data)*.17)))
            for i,r in enumerate(data):
                point=float(r['oe']);lo=float(r['ci_low']);hi=float(r['ci_high']);color='#35678c' if r['threshold']=='5' else '#d78a39'
                ax.plot([lo,hi],[i,i],color=color);ax.plot(point,i,'o',color=color,ms=3)
            ax.axvline(1,color='black',ls='--',lw=1);ax.set_yticks(range(len(data)),[f"{r['label']} AC≥{r['threshold']}" for r in data],fontsize=5);ax.set_xlabel('Observed / LOCO expected; EUR LD-block bootstrap 95% CI');save(fig,number,'Archaic carrier prevalence' if archaic else 'HERV body, family, boundary and repeat controls')
        else:inventory.append(dict(figure=number,status='no_supported_enrichment',outputs=[]))
    centroid=evolution/'pca'/'centroids'/'allele_centroid_trajectories.png'
    if centroid.exists():
        outputs=[]
        for ext in ('png','pdf'):
            source=centroid.with_suffix('.'+ext)
            if source.exists():dest=figs/f'Fig6.{ext}';shutil.copyfile(source,dest);outputs.append(str(dest))
        inventory.append(dict(figure=6,status='measured',outputs=outputs,title='PC1–PC2–time carrier centroids'))
    else:inventory.append(dict(figure=6,status='no_projectable_carrier_centroids',outputs=[]))
    robust=[r for r in enrichment if r['sensitivity']=='exclude_combined' and r['status']=='complete' and r['fdr'] not in {'','None','NA'} and float(r['fdr'])<.05 and (float(r['ci_low'])>1 or float(r['ci_high'])<1)]
    paired={r['label'] for r in robust if r['threshold']=='5'} & {r['label'] for r in robust if r['threshold']=='10'}
    agreement=bool(effects) and all(r['direction_agreement']=='True' for r in effects) and {r['threshold'] for r in effects}>={'5','10'}
    decision='effect_size_review_required' if paired or agreement else 'primary_claim_not_established'
    lines=['# allele3 focused results','','AC≥5 为主分析；AC≥10 为严格复现；AC1 仅作 QC/emergence。','','## Population']
    lines += [f"- AC≥{r['global_ac_threshold']}, target AN={r['target_an']}, {r['classification']}: {r['n']} (fraction={r['fraction'] or 'NA'})" for r in pop]
    lines += ['','## 1KG replication']+[f'- {r}' for r in effects]
    lines += ['','1KG 使用分层抽取的发现位点；其分类比例是该位点集合的条件结果，不能当作全基因组比例。缺失 exact ALT 不作为零。','','## Enrichment',f'- 完成 {len(enrichment)} 个 label/AC/sensitivity 组合；完整 difficult-region 控制下 AC5/10 均通过的 label 数：{len(paired)}。','- 所有 P 值来自中心化 LD-block bootstrap，并进行 BH FDR；不使用 naive per-site z 作为主 P。','- 未提供 low-mappability 等必要注释时，不把未完成的严格敏感性检验写成通过。','','## Scientific decision',f'- `{decision}`','- 只有大效应、AC10 一致、独立复现/完整敏感性检验通过，才可支持主论文假设；“大效应”的科学判读仍需结合结果，不能靠显著性替代。','- 未满足条件时不把 ancestry hypothesis 放进标题；不把未测定与阴性结果混为一谈。','','## Main figures (maximum six)']
    for item in inventory:
        lines.append(f"- Fig{item['figure']}: {item['status']} "+' '.join(f'[{Path(p).name}]({Path(p).relative_to(out)})' for p in item['outputs']))
    lines += ['','旧版 evolution/figures 中的额外图作为支持材料；AlphaGenome 默认不执行。']
    report=out/'FOCUSED_REPORT.md';report.write_text('\n'.join(lines)+'\n')
    result=dict(status='complete',version='4.0.0-focused-20261024',report=str(report),decision=decision,main_claim_ready=False,figures=inventory,robust_ac5_ac10_labels=sorted(paired),population_direction_agreement=agreement)
    (out/'report_manifest.json').write_text(json.dumps(result,indent=2)+'\n');return result


def main_report(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--population-dir',required=True);p.add_argument('--modern-dir');p.add_argument('--evolution-dir');p.add_argument('--replication-dir');p.add_argument('--herv-enrichment');p.add_argument('--archaic-enrichment');p.add_argument('--out',required=True)
    print(json.dumps(make_report(p.parse_args(argv)),indent=2))


