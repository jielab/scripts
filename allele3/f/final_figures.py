"""Static scientific exports and an explicit availability inventory for Fig1–9."""
from __future__ import annotations
from pathlib import Path
from collections import defaultdict, Counter
import csv, json, math, re
from core import dump_json, fingerprint, open_text, tsv_writer

def _plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt

def _save(fig,path):
    paths=[]
    for ext in ('png','pdf'):
        p=Path(str(path)+'.'+ext); fig.savefig(p,dpi=180,bbox_inches='tight'); paths.append(str(p))
    _plot().close(fig); return paths

def rows(path):
    with open_text(path) as f: yield from csv.DictReader(f,delimiter='\t')

def static_pca(payload,directory,candidate_path=None):
    import numpy as np
    plt=_plot(); d=Path(directory); d.mkdir(parents=True,exist_ok=True)
    points=payload['points']; x=np.array([r['pc1'] for r in points]); y=np.array([r['pc2'] for r in points]); age=np.array([r['age_bp'] for r in points])
    regions=sorted({r['region'] for r in points}); colors={r:plt.get_cmap('tab20')(i%20) for i,r in enumerate(regions)}
    fig=plt.figure(figsize=(10,8)); ax=fig.add_subplot(111,projection='3d')
    for r in regions:
        idx=[i for i,p in enumerate(points) if p['region']==r]
        ax.scatter(x[idx],y[idx],age[idx]/1000,s=3,alpha=.45,c=[colors[r]],label=r)
    ax.set(xlabel='PC1 (shared score)',ylabel='PC2 (shared score)',zlabel='Age (kyr BP)',title='Fixed PCA basis, all individuals')
    if len(regions)<=16: ax.legend(fontsize=7,loc='upper left',bbox_to_anchor=(1,1))
    outputs=_save(fig,d/'pca_time_overview')
    maxbin=max(1,int(max(age,default=0)//2000)+1); cols=4; nrows=math.ceil(maxbin/cols)
    fig,axes=plt.subplots(nrows,cols,figsize=(14,3.3*nrows),squeeze=False,sharex=True,sharey=True)
    for k,ax in enumerate(axes.flat):
        if k>=maxbin: ax.set_visible(False); continue
        mask=(age>=k*2000)&(age<(k+1)*2000)
        ax.scatter(x[mask],y[mask],s=4,c=[colors[p['region']] for i,p in enumerate(points) if mask[i]],alpha=.5)
        ax.set(title=f'{k*2}–{(k+1)*2} kyr BP; n={sum(mask)}',xlabel='PC1',ylabel='PC2')
    fig.suptitle('Fixed 2-kyr slices; identical coordinate axes'); outputs+=_save(fig,d/'pca_time_slices')
    groups=defaultdict(list)
    for i,p in enumerate(points): groups[(p['region'],int(p['age_bp']//2000))].append(i)
    centroid=[]
    fig,ax=plt.subplots(figsize=(9,7))
    for region in regions:
        trajectory=[]
        for (r,k),idx in sorted(groups.items()):
            if r!=region: continue
            a=x[idx]; b=y[idx]; n=len(idx)
            se1=float(a.std(ddof=1)/math.sqrt(n)) if n>=2 else None
            se2=float(b.std(ddof=1)/math.sqrt(n)) if n>=2 else None
            centroid.append(dict(region=r,start_bp=k*2000,end_bp=(k+1)*2000,n=n,pc1=float(a.mean()),pc2=float(b.mean()),
                pc1_se=se1 if se1 is not None else 'NA',pc2_se=se2 if se2 is not None else 'NA'))
            trajectory.append((k,a.mean(),b.mean()))
            ax.errorbar(a.mean(),b.mean(),xerr=1.96*se1 if se1 is not None else None,yerr=1.96*se2 if se2 is not None else None,color=colors[r],alpha=.55,fmt='o',ms=3)
        trajectory.sort(reverse=True)
        ax.plot([t[1] for t in trajectory],[t[2] for t in trajectory],color=colors[region],label=region,alpha=.75)
    ax.set(xlabel='PC1',ylabel='PC2',title='Region centroids: older to younger; descriptive 95% mean intervals')
    if len(regions)<=16: ax.legend(fontsize=7)
    outputs+=_save(fig,d/'pca_time_centroids')
    if centroid:
        h,w=tsv_writer(d/'pca_time_centroids.tsv',list(centroid[0]))
        with h: w.writerows(centroid)
    variants=payload['variants']; by_key={(r['chrom'],int(r['pos']),r['ref'],r['alt']):r for r in variants}
    loci=[]
    if candidate_path and Path(candidate_path).exists():
        seen=set(); complement=str.maketrans('ACGT','TGCA')
        for r in rows(candidate_path):
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
            ax.scatter(x,y,age/1000,c='#bcbcbc',s=4,alpha=.3)
            if vals is not None:
                idx=[i for i,v in enumerate(vals) if v is not None]
                scat=ax.scatter(x[idx],y[idx],age[idx]/1000,c=[vals[i] for i in idx],vmin=0,vmax=2,cmap='viridis',s=8)
                fig.colorbar(scat,ax=ax,shrink=.5,label='Copies of the named allele')
            ax.set(title=label+(' (unavailable)' if vals is None else ''),xlabel='PC1',ylabel='PC2',zlabel='kyr BP')
        fig.suptitle(name+'; grey = unavailable, 0/1/2 = copies')
        outputs+=_save(fig,d/('pca_time_locus_'+re.sub(r'[^A-Za-z0-9_.-]','_',name)))
    return dict(outputs=outputs,centroid_rows=len(centroid),locus_figures=len(loci),slice_width_bp=2000)

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
        for r in rows(spectrum):
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
    rs=[r for r in rows(ld) if r['ancestry']=='EUR' and r['threshold']=='1' and r['sensitivity']=='all' and r['population_mode']=='all'][:20] if ld.exists() else []
    bar(2,'EUR LD blocks: top descriptive O/E', [r['block_id'] for r in rs],[float(r['oe']) for r in rs],'Observed / LOCO expected')
    consensus=summary/'ld_hotspot_consensus.tsv'
    if consensus.exists() and make_plots:
        overlaps=list(rows(consensus))
        if overlaps:
            fig,ax=plt.subplots(figsize=(12,5));chroms=sorted({r['chrom'] for r in overlaps},key=lambda c:(not c.isdigit(),int(c) if c.isdigit() else c));ys={c:i for i,c in enumerate(chroms)}
            for r in overlaps:
                ax.plot([int(r['start'])/1e6,int(r['end'])/1e6],[ys[r['chrom']]]*2,lw=4,alpha=.8,color=plt.get_cmap('viridis')((int(r['n_ancestries'])-2)/2))
            ax.set_yticks(range(len(chroms)),chroms);ax.set(xlabel='GRCh38 position (Mb)',ylabel='Chromosome',title='Physical overlap of ancestry-specific top 1% LD blocks; 2–4 supporting ancestries')
            inventory[-1]['outputs']+=_save(fig,figures/'Fig2_consensus')
    erv=summary/'erv_enrichment.tsv'; totals=defaultdict(lambda:[0.,0.,0])
    if erv.exists():
        for r in rows(erv):
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
        for r in rows(pops):
            if r['threshold']=='1': counts[r['classification']]+=int(r['count'])
    bar(4,'Population partitioning of multi-ALT loci',list(counts),list(counts.values()),'Loci')
    candidates=out/'universes'/'U2_all_strong.tsv'; cases=list(rows(candidates)) if candidates.exists() else []
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
            for r in rows(trajectory):
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
        for r in rows(ukb):
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
        for r in rows(geometry):
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
    from ldblocks import BlockIndex,read_blocks
    from erv import parse_joined_rmsk,ERVOverlapSweep
    out=Path(evolution_dir);candidate=out/'universes'/'candidates.tsv'
    if not candidate.exists(): return dict(status='unavailable_candidates')
    data=[];lookup=defaultdict(list)
    for r in rows(candidate):
        x=dict(site_id=r['site_id'],selection_id=r['selection_id'],extra_alt=r['extra_alt'],chrom=r['chrom'],pos=int(r['pos']),
            group=r['group'],persistence_ac_ge5=int(float(r['extra_ac'])>=5),anchor_strong=int(float(r['posterior'])>=.99),
            ukb_exact_info_ge08='NA',ld_block38='NA',ld_oe='NA',inside_erv_core='NA',functional_geometry_measured='NA')
        data.append(x);lookup[(r['site_id'],r['selection_id'],r['extra_alt'])].append(x)
    ukb=out/'ukb'/'U0_mfi.tsv'
    if ukb.exists():
        for r in rows(ukb):
            if r['target_role']!='extra': continue
            for x in lookup.get((r.get('site_id'),r.get('selection_id'),r.get('extra_alt')),[]):
                if r['ukb_status']=='unavailable_file' or (r['ukb_status']=='present' and r['ukb_info']=='NA'): continue
                x['ukb_exact_info_ge08']=int(r['ukb_status']=='present' and float(r.get('ukb_info') or 0)>=.8) if r['ukb_status']=='present' else 0
    run=Path(modern_dir)/'run.json';options=json.loads(run.read_text()).get('options',{}) if run.exists() else {}
    if options.get('ld_block_dir'):
        blocks,_=read_blocks(Path(options['ld_block_dir'])/'EUR.38.bed','EUR','38');idx=BlockIndex(blocks)
        en=Path(modern_dir)/'summary'/'ld_block_enrichment.tsv';oe={}
        if en.exists():
            oe={r['block_id']:r['oe'] for r in rows(en) if r['ancestry']=='EUR' and r['sensitivity']=='all' and r['threshold']=='1' and r['population_mode']=='all'}
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
        measured={r['site_id'] for r in rows(geo) if r['atlas_status']=='ok'}
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
