"""Manuscript figures, drawn only for endpoints with finite estimates."""
from collections import Counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from research_common import number

COLORS = ('#116768', '#ad5331', '#585486')
VARIABLES = {'sd_site_mean':'Between-position SD', 'mean_abs_gradient':'Adjacent-position gradient',
             'mean_alt_range':'Within-position ALT range'}
GEOMETRY = {'normalized_weak_pair_distance_minus_other_pairs':'Lower-AF pair vs both top-AF pairs',
            'normalized_least_ALT_distance_to_middle_minus_top':'Lowest to middle vs highest AF',
            'normalized_extra_ALT_distance_to_top_minus_unreported':'Second vs unreported ALT, relative to top'}


def save(fig, root, name):
    fig.savefig(root / f'{name}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def point(ax, x, y, lo, hi, color, horizontal):
    if lo is not None and hi is not None:
        if horizontal:
            ax.hlines(y, lo, hi, color=color, linewidth=1.6)
        else:
            ax.vlines(x, lo, hi, color=color, linewidth=1.6)
    ax.plot(x, y, 'o', color=color, ms=4)


def series(rows):
    return sorted({r['target'] for r in rows})


def legend(ax, groups):
    ax.legend([Line2D([],[],ls='',marker='o',color=COLORS[i]) for i in range(len(groups))],
              groups, frameon=False, fontsize=9)


def forest(root, filename, rows, labels, field, xlabel, null):
    rows = [r for r in rows if number(r.get(field)) is not None and r['key'] in labels]
    if not rows:
        return False
    keys = [k for k in labels if any(r['key']==k for r in rows)]
    groups = series(rows)
    fig, ax = plt.subplots(figsize=(8, 1.65+0.48*len(keys)), layout='constrained')
    for r in rows:
        i = groups.index(r['target'])
        y = keys.index(r['key']) + (i-(len(groups)-1)/2)*.19
        point(ax, number(r[field]), y, number(r.get('ci_low')), number(r.get('ci_high')), COLORS[i], True)
    ax.set_yticks(range(len(keys)), [labels[k] for k in keys])
    ax.axvline(null, color='.55', ls='--', lw=1)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.grid(axis='x', alpha=.16)
    legend(ax, groups)
    save(fig, root, filename)
    return True


def figures(root, tables):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    output = {}
    regions = [r for r in tables.get('T2_region_enrichment',[]) if r['sensitivity']=='all' and r['outcome']=='multi_ALT']
    outcomes = {'multi_ALT':'At least two ALTs','two_ALT':'Exactly two ALTs','three_ALT':'Exactly three ALTs'}
    labels = {(outcome,tail):f'{title} · {tail} gradient' for outcome,title in outcomes.items() for tail in ('low','high')}
    gradient = [dict(r,key=(r['outcome'],r['region'].rsplit('_',1)[1])) for r in tables.get('T2_region_enrichment',[])
                if r['reference']=='alphagenome' and r['sensitivity']=='all' and r['outcome'] in outcomes
                and r['region'] in ('AVI_mean_abs_gradient_low','AVI_mean_abs_gradient_high')]
    if forest(root, 'F1_avi_windows', gradient, labels, 'oe','Observed / expected sites (95% CI)',1):
        output['F1_avi_windows'] = 'Gradient tails in scored windows. Outcomes count qualifying nonreference ALTs. Each O/E has its own fitted expectation. The combined group includes the two component groups.'
    distance = [r for r in tables.get('T3_distance_profile',[]) if r['metric'] in VARIABLES and number(r['case_minus_control']) is not None]
    if distance:
        groups = series(distance)
        fig, axs = plt.subplots(1,3,figsize=(12.4,3.9),layout='constrained')
        bands = sorted({(int(r['inner_bp']),int(r['outer_bp'])) for r in distance})
        for (metric,title),ax in zip(VARIABLES.items(),axs):
            for r in distance:
                if r['metric'] != metric:
                    continue
                i = groups.index(r['target'])
                x = bands.index((int(r['inner_bp']),int(r['outer_bp']))) + (i-(len(groups)-1)/2)*.12
                point(ax,x,number(r['case_minus_control']),number(r.get('ci_low')),number(r.get('ci_high')),COLORS[i],False)
            ax.set_xticks(range(len(bands)), [f'{lo/1000:g}–{hi/1000:g}' for lo,hi in bands],rotation=35,ha='right')
            ax.axhline(0,color='.6',ls='--',lw=.9)
            ax.set(xlabel='Distance band (kb)',title=title)
        axs[0].set_ylabel('Multi-ALT case minus matched controls')
        legend(axs[-1],groups)
        save(fig,root,'F2_distance_profile')
        output['F2_distance_profile'] = 'Matched case minus control within each distance band. Bars show upstream chromosome-bootstrap 95% intervals. Available cases can differ by band.'
    # Metrics have different units, so each paired comparison has its own axis.
    paired = [r for r in tables.get('T4_paired_near_far',[]) if r['metric'] in VARIABLES and number(r['estimate']) is not None]
    if paired:
        fig, axs = plt.subplots(1,3,figsize=(11.2,3.2),layout='constrained')
        groups = series(paired)
        for ax,(metric,title) in zip(axs,VARIABLES.items()):
            for r in paired:
                if r['metric']==metric:
                    i=groups.index(r['target'])
                    point(ax,i,r['estimate'],r.get('ci_low'),r.get('ci_high'),COLORS[i],False)
            ax.axhline(0,color='.6',ls='--',lw=.9)
            ax.set_xticks(range(len(groups)),groups)
            ax.set(xlim=(-.55,len(groups)-.45),title=title)
        axs[0].set_ylabel('Case near−far minus control near−far')
        save(fig,root,'F3_paired_near_far')
        output['F3_paired_near_far'] = 'Near-versus-far contrast within the same focal sites and the same weighted controls. Positive estimates indicate a larger near-to-far difference in cases. T4 reports both component changes and coverage.'
    sensitivities = {'all_available':'All eligible','equal_reported_AN':'Equal AN'}
    allele = [dict(r,key=(r['metric'],r['sensitivity'])) for r in tables.get('T5_within_site',[])
              if r['scale']=='PHRED' and r['sensitivity'] in sensitivities and r['metric'] in GEOMETRY]
    labels = {(m,s):f'{label} · {group}' for m,label in GEOMETRY.items() for s,group in sensitivities.items()}
    if forest(root,'F4_af_geometry',allele,labels,'estimate','Normalized PHRED distance contrast (95% CI)',0):
        output['F4_af_geometry'] = 'AF ranks use the selected population. Equal AN restricts to ALTs with the same reported allele number. The last contrast compares the second observed ALT and the unreported substitution, each against the top observed ALT. RAW results are in T5.'
    family = [dict(r,key=r['label'],target='AC ≥ '+r['threshold']) for r in tables.get('T6b_common_support',[])
              if r.get('sensitivity')=='exclude_combined' and r.get('threshold') in ('5','10')]
    family_labels = {k:k.removeprefix('family:').replace('_',' ') for k in
                     ('family:ERV1','family:ERVK','family:ERVL','family:ERVL-MaLR','inside_erv_core','inside_any_ltr')}
    if forest(root,'F5b_herv_families',family,family_labels,'oe','Observed / expected on common support (95% CI)',1):
        output['F5b_herv_families'] = 'HERV families and pooled repeat annotations after combined exclusions. T6b reports excluded counts and the fraction with supported backgrounds.'
    peaks = [dict(r,key=(r['family'],r['boundary']),target='Held-out chromosomes',
                  ci_low=r['validation_ci_low'],ci_high=r['validation_ci_high'])
             for r in tables.get('T6c_boundary_validation',[]) if r.get('validation_threshold')=='5'
             and number(r.get('validation_local_ratio')) is not None]
    labels = {r['key']:f"{r['family']} {r['boundary']} · {r['start_bp']} to {r['end_bp']} bp" for r in peaks}
    if forest(root,'F5c_herv_boundaries',peaks,labels,'validation_local_ratio','Local enrichment ratio in frozen window (95% CI)',1):
        output['F5c_herv_boundaries'] = 'AC ≥ 5. Windows were chosen on discovery chromosomes and evaluated on held-out chromosomes. Distances follow aligned repeat orientation. T6c retains exploratory thresholds and status.'
    overlaps = [dict(r,key=(r['reference'],r['region'])) for r in regions if r['reference']=='herv' or (r['reference']=='archaic' and r.get('gu_cohort')=='1kg')]
    labels = {k:f'{k[0]} · {k[1]}' for k in sorted({r['key'] for r in overlaps})}
    if forest(root,'F5_region_overlap',overlaps,labels,'oe','Observed / expected multi-ALT sites',1):
        output['F5_region_overlap'] = 'Conditional O/E in annotated HERV families and 1KG IBDmix tracts. See T2 for each YGB batch separately. Each reference retains its original analysis and multiplicity adjustment. A tract overlap alone does not identify which allele is introgressed.'
    matched = []
    for row in tables.get('T8_matched_cohort_enrichment',[]):
        region = row.get('region','')
        if row['sensitivity']!='all' or row['outcome']!='multi_ALT':continue
        if region in ('AVI_mean_abs_gradient_low','AVI_mean_abs_gradient_high') or row['reference']=='herv':
            matched.append(dict(row,key=row['population']+' · '+region))
    labels = {s:s for s in sorted({r['key'] for r in matched})}
    if forest(root,'F8_matched_cohorts',matched,labels,'oe','Observed / expected (95% CI)',1):
        output['F8_matched_cohorts']='gnomAD and All of Us on shared eligible positions, by population. Points are cohort-specific estimates; intervals use genomic blocks.'
    splice=[dict(r,key='direct') for r in tables.get('T5b_direct_splicing',[])
            if r['metric']=='normalized_weak_pair_distance_minus_other_pairs' and r['sensitivity']=='all_available']
    if forest(root,'F8b_direct_splicing',splice,{'direct':'Lower-AF pair vs remaining pairs'},'estimate','Normalized splicing-score distance contrast (95% CI)',0):
        output['F8b_direct_splicing']='Direct AG_SPLICE score among three qualifying ALT sites with ranked AF and complete, variable splicing predictions.'
    stratified=[dict(r,key=str(r['qualifying_alts'])+' ALTs · '+VARIABLES[r['metric']])
                for r in tables.get('T4b_paired_near_far_by_alt_count',[]) if r['metric'] in VARIABLES]
    labels={s:s for s in sorted({r['key'] for r in stratified})}
    if forest(root,'F3b_distance_by_alt_count',stratified,labels,'estimate','Matched near-minus-far difference (metric-specific units)',0):
        output['F3b_distance_by_alt_count']='Separate paired contrasts for two and three qualifying ALTs. Each metric retains its scoring units.'
    return output


def pca_figure(root, points):
    cleaned = []
    for row in points:
        x, y, age = (number(row.get(k)) for k in ('pc1','pc2','age_bp'))
        cohort = row.get('cohort')
        if None not in (x,y,age) and age >= 0 and cohort in ('ancient','modern'):
            cleaned.append((cohort,x,y,age))
    if not cleaned:
        return []
    fig = plt.figure(figsize=(7.4,5.7),layout='constrained')
    ax = fig.add_subplot(projection='3d')
    counts = Counter()
    for cohort, color in [('ancient',COLORS[0]),('modern',COLORS[1])]:
        subset = [(x,y,a/1000) for c,x,y,a in cleaned if c==cohort]
        if subset:
            xs,ys,ages = zip(*subset)
            ax.scatter(xs,ys,ages,s=4,alpha=.32,color=color,rasterized=True,label=f'{cohort} (n={len(subset):,})')
        for age in (a for c,_,_,a in cleaned if c==cohort):
            counts[(cohort,int(age//1000))] += 1
    ax.set(xlabel='PC1',ylabel='PC2',zlabel='Sample age (ka BP)')
    ax.legend(frameon=False)
    save(fig,root,'S1_pca_time')
    return [dict(cohort=c,age_from_bp=e*1000,age_to_bp_exclusive=(e+1)*1000,n_plotted=n)
            for (c,e),n in sorted(counts.items())]
