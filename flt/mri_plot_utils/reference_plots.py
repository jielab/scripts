#!/usr/bin/env python3
"""Plot genuine regression outputs; never fill unmeasured modalities with synthetic beta.

python s4_4_mri.plot.py --results /path/01_all_linear_associations.csv --out figures
Optional: --tract-manifest verified_tract_masks.csv for anatomical panel b.
Default anatomy uses the verified local UKB atlas. --publication-only assembles the MRI Figure 4 layout.
"""
import argparse
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(p) for p in (ROOT/'.tools/enigma_packages',
    ROOT/'.tools/ENIGMA', ROOT/'.tools/mri_plot_packages') if p.exists()]
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
# Compatibility with this project's ENIGMA version and NumPy 2 / Matplotlib 3.10.
if not hasattr(np, 'string_'): np.string_ = np.bytes_
if not hasattr(np, 'unicode_'): np.unicode_ = np.str_
if not hasattr(matplotlib.cm, 'get_cmap'):
    matplotlib.cm.get_cmap = lambda name=None, lut=None: matplotlib.colormaps.get_cmap(name or 'viridis').resampled(lut or 256)

# FreeSurfer aparc order, omitting unknown and corpuscallosum.
DK = '''bankssts caudalanteriorcingulate caudalmiddlefrontal cuneus entorhinal
fusiform inferiorparietal inferiortemporal isthmuscingulate lateraloccipital
lateralorbitofrontal lingual medialorbitofrontal middletemporal parahippocampal
paracentral parsopercularis parsorbitalis parstriangularis pericalcarine postcentral
posteriorcingulate precentral precuneus rostralanteriorcingulate rostralmiddlefrontal
superiorfrontal superiorparietal superiortemporal supramarginal frontalpole
temporalpole transversetemporal insula'''.split()
KEYS = [h+'_'+r for h in ('lh','rh') for r in DK]
FAMILIES = ['ThicknessDK','FA27','MD27']
plt.rcParams.update({'font.family':'Arial','font.size':10,'pdf.fonttype':42,
                     'axes.spines.top':False,'axes.spines.right':False})

def stars(q):
    return '***' if q<.001 else '**' if q<.01 else '*' if q<.05 else ''

def save(fig,out,name):
    for ext in ('png','pdf'):
        fig.savefig(out/(name+'.'+ext),dpi=220,bbox_inches='tight',facecolor='white')
    plt.close(fig)

def short_label(s):
    return (str(s).replace('Volume of grey matter in ','').replace('Volume of ','')
            .replace(' (left)',' (L)').replace(' (right)',' (R)'))


def panel_a(r,out):
    from enigmatoolbox.plotting import plot_cortical
    from enigmatoolbox.utils.parcellation import parcel_to_surface
    d=r[(r.exposure==EXPOSURE)&(r.family=='ThicknessDK')].copy().set_index('atlas_key')
    if d.index.duplicated().any():raise ValueError('Duplicate DK labels')
    unknown=set(d.index)-set(KEYS)
    if unknown:raise ValueError(f'Unknown DK labels: {unknown}')
    d=d.reindex(KEYS)
    good=d.status.eq('ok') & np.isfinite(d.beta)
    if not good.any():raise ValueError('No successful thickness results')
    b=d.beta.where(good).to_numpy(float)
    limit=max(np.nanmax(np.abs(b)),1e-6)
    audit=d[['field_id','beta','FDR_family','status']].copy()
    audit['display_status']=np.where(~good,'missing_or_failed',np.where(d.FDR_family<.05,'FDR_significant','not_significant'))
    audit.to_csv(out/'panel_a_mapping.csv',index_label='atlas_key')
    renders=out/'cache/surface_renders';renders.mkdir(parents=True,exist_ok=True)
    for tag,values in [('all',b),('FDR05',np.where(good & (d.FDR_family<.05),b,np.nan))]:
        if not np.isfinite(values).any():
            print('No significant thickness parcels; FDR surface omitted.');continue
        surf=parcel_to_surface(values,'aparc_fsa5',fill=np.nan)
        plot_cortical(array_name=surf,surface_name='fsa5',cmap='RdBu_r',
                      color_bar=True,color_range=(-limit,limit),nan_color=(.87,.87,.87,1),
                      size=(1600,500),screenshot=True,filename=str(renders/f'a_thickness_{tag}.png'),
                      transparent_bg=False,interactive=False,offscreen=True)
        fig=plt.figure(figsize=(16,5.7))
        ax=fig.add_axes([0,.09,1,.84]);ax.imshow(plt.imread(renders/f'a_thickness_{tag}.png'));ax.axis('off')
        fig.text(.025,.96,'a  FLT and cortical thickness | '+('all measured regions' if tag=='all' else 'family FDR < 0.05'),fontsize=16,va='top')
        fig.text(.86,.88,'MRI SD per FLT point',fontsize=10)
        fig.text(.24,.83,'Left hemisphere',ha='center',fontsize=11)
        fig.text(.68,.83,'Right hemisphere',ha='center',fontsize=11)
        fig.text(.025,.035,'covs.all. Grey: '+('unmeasured/failed parcels.' if tag=='all' else 'unmeasured or FDR >= 0.05.')+'  Bilateral temporal poles unmeasured.',fontsize=10)
        save(fig,out,f'a_thickness_{tag}')
    # Reference-specific annotation selection omitted; mapping CSV preserves every parcel.



def tract_forest(r,out):
    fig,axs=plt.subplots(1,2,figsize=(13,10),sharey=True)
    for ax,fam in zip(axs,['FA27','MD27']):
        d=r[(r.exposure==EXPOSURE)&(r.family==fam)].copy().sort_values('field_id')
        y=np.arange(len(d));ok=d.status.eq('ok')
        for j,(_,v) in enumerate(d.iterrows()):
            if not ok.iloc[j]:continue
            color='#b2182b' if v.beta>=0 else '#2166ac'
            ax.errorbar(v.beta,j,xerr=[[v.beta-v.LCI],[v.HCI-v.beta]],fmt='o',
                        color=color,alpha=1 if v.FDR_family<.05 else .3,capsize=2)
            ax.annotate(stars(v.FDR_family),(v.HCI,j),xytext=(3,0),textcoords='offset points',va='center')
        ax.set_yticks(y,d.atlas_key);ax.axvline(0,color='.6',lw=1)
        ax.set_title(fam.replace('27','')+(' | no FDR-significant associations' if not (d.FDR_family<.05).any() else ''));ax.set_xlabel('MRI SD per FLT point (95% CI)')
    axs[0].invert_yaxis()
    fig.suptitle('b  White-matter tract associations | coefficient view',x=.05,ha='left')
    fig.text(.05,.015,'Opaque points: family FDR < 0.05. Faded points: not significant. Missing/failed tests have no point.',fontsize=9)
    fig.tight_layout(rect=(0,.04,1,.96));save(fig,out,'b_tract_coefficients')



def supplementary_forests(r,out):
    """Every volume and WMH estimate; categories change display only, not FDR."""
    d=r[r.exposure.eq(EXPOSURE) & r.family.isin(['Volume139','WMH'])]
    for (fam,category),group in d.groupby(['family','category'],sort=False):
        group=group.sort_values('field_id')
        for page,start in enumerate(range(0,len(group),32),1):
            chunk=group.iloc[start:start+32]
            fig,ax=plt.subplots(figsize=(12,max(3.5,len(chunk)*.29+1.6)))
            for j,(_,v) in enumerate(chunk.iterrows()):
                if v.status!='ok':continue
                ax.errorbar(v.beta,j,xerr=[[v.beta-v.LCI],[v.HCI-v.beta]],fmt='o',
                    color='#b2182b' if v.beta>=0 else '#2166ac',
                    alpha=1 if v.FDR_family<.05 else .3,capsize=2)
                ax.annotate(stars(v.FDR_family),(v.HCI,j),xytext=(4,0),textcoords='offset points',va='center')
            ax.set_yticks(np.arange(len(chunk)),[short_label(s) for s in chunk.label],fontsize=8)
            ax.invert_yaxis();ax.axvline(0,color='.6',lw=1)
            ax.set_xlabel('MRI SD per FLT point (95% CI)')
            ax.set_title(f'FLT | {fam}: {category} | page {page}',loc='left')
            fig.text(.015,.015,'BH across all 139 volumes (WMH: one test). Opaque: q<0.05; faded: not significant.',fontsize=9)
            fig.tight_layout(rect=(0,.045,1,1))
            safe=''.join(c if c.isalnum() else '_' for c in str(category))
            save(fig,out,f's_FLT_{fam}_{safe}_{page:02d}')




