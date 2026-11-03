#!/usr/bin/env python3
"""mechanism analysis steps. Run with dose.sh; imports do not execute analysis."""

# Shared PNG byte limit; PDF/SVG export settings are preserved.
import sys as _png_sys
from pathlib import Path as _PngPath
_png_sys.path.append(str(_PngPath(__file__).resolve().parents[2] / "0f"))
from png_output import save_png_figure

import importlib.util
import sys
from pathlib import Path

if "dose_common" not in sys.modules:
    spec = importlib.util.spec_from_file_location("dose_common", Path(__file__).with_name("0.common.py"))
    common = importlib.util.module_from_spec(spec)
    sys.modules["dose_common"] = common
    spec.loader.exec_module(common)
else:
    common = sys.modules["dose_common"]


# ===========================================================================
# mechanism / mri_prepare
def run_mri_prepare():
    """Stream selected fields without decompressing the whole UKB table to disk.

    Called by dose.sh before the MRI regression step. Keep participant data locally.
    No imputation or numeric transformation is performed.
    """
    import argparse
    import csv
    import gzip
    import json
    import os
    from pathlib import Path
    import sys
    output_dir = common.output_dir; ukb_dir = common.ukb_dir; mri_dictionary = common.mri_dictionary


    extract = common.extract_columns

    def main():
        p = argparse.ArgumentParser(description=__doc__)
        p.add_argument('--project', default=str(output_dir()))
        p.add_argument('--ukb', default=str(ukb_dir()))
        p.add_argument('--dictionary', default=str(mri_dictionary()))
        p.add_argument('--mode', choices=['visit','mri','all'], default='all')
        a = p.parse_args()
        out = Path(a.project)/'ukb/05mechanism/04mri/inputs'
        if a.mode in ('visit','all'):
            extract(Path(a.ukb)/'rap/raw/pheno.tab.gz',out/'mri_visit_covariates.csv',
                    ['eid','p21003_i2','p54_i2','p53_i2'],
                    ['eid','imaging_age','imaging_site','imaging_date'])
        if a.mode in ('mri','all'):
            with open(a.dictionary, encoding='utf-8-sig', newline='') as f:
                fields = ['p'+row['field_id']+'_i2' for row in csv.DictReader(f)]
            raw = Path(a.ukb)/'rap/img.tab.gz'
            with gzip.open(raw,'rt') as f:
                available = set(f.readline().strip().split('\t'))
            fields = list(dict.fromkeys(['eid']+[f for f in fields if f in available]+['p26521_i2']))
            extract(raw,out/'mri_selected.tsv',fields)


    main()


# ===========================================================================
# mechanism / mri_plot
def run_mri_plot():
    """Plot genuine regression outputs; never fill unmeasured modalities with synthetic beta.

    ./dose.sh mechanism --steps mri_plot -- --results /path/01_all_linear_associations.csv --out figures
    Optional: --tract-manifest verified_tract_masks.csv for anatomical panel b.
    Default anatomy uses the verified local UKB atlas. --publication-only assembles the MRI Figure 4 layout.
    """
    import argparse
    from pathlib import Path
    output_dir = common.output_dir
    import sys
    ROOT = output_dir()
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
            save_png_figure(fig, out/(name+'.'+ext),dpi=220,bbox_inches='tight',facecolor='white')
        plt.close(fig)

    def short_label(s):
        return (str(s).replace('Volume of grey matter in ','').replace('Volume of ','')
                .replace(' (left)',' (L)').replace(' (right)',' (R)'))

    def observed_volume_preview(r,out):
        d=r[(r.exposure=='dose') & (r.family=='Volume139')].copy()
        d=d[d.status=='ok']
        if d.empty:return
        hits=d[d.FDR_family<.05].sort_values('P')
        fig,ax=plt.subplots(figsize=(11,max(3.5,len(hits)*.55+1.5)))
        if hits.empty:
            ax.text(.5,.5,'No volume association passes FDR < 0.05',ha='center',transform=ax.transAxes)
        else:
            y=np.arange(len(hits))
            ax.errorbar(hits.beta,y,xerr=np.array([hits.beta-hits.LCI,hits.HCI-hits.beta]),
                        fmt='o',color='#b2182b',capsize=4,lw=2)
            ax.set_yticks(y,[short_label(s) for s in hits.label]);ax.invert_yaxis()
            for yi,(_,v) in zip(y,hits.iterrows()):
                ax.annotate(f' q={v.FDR_family:.4f}',(v.HCI,yi),xytext=(5,0),textcoords='offset points',va='center',fontsize=9)
            ax.axvline(0,color='.65',ls='--',lw=1)
            ax.set_xlim(min(0,float(hits.LCI.min()))-.003,float(hits.HCI.max())+.008)
        ax.set_xlabel('Standardized beta (95% CI)')
        ax.set_title('DOSE and brain volumes | observed MRI outcomes',loc='left',pad=18)
        fig.text(.01,.01,'FDR across 139 volumes. These are volume associations, not MODERN-style thickness/FA/MD panels.',fontsize=9,color='.35')
        fig.tight_layout(rect=(0,.04,1,1));save(fig,out,'DOSE_volume_FDR05')

    def panel_a(r,out):
        common.prepare_enigma()
        from enigmatoolbox.plotting import plot_cortical
        from enigmatoolbox.utils.parcellation import parcel_to_surface
        d=r[(r.exposure=='dose')&(r.family=='ThicknessDK')].copy().set_index('atlas_key')
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
            fig.text(.025,.96,'a  DOSE and cortical thickness | '+('all measured regions' if tag=='all' else 'family FDR < 0.05'),fontsize=16,va='top')
            fig.text(.86,.88,'Standardized beta',fontsize=10)
            fig.text(.24,.83,'Left hemisphere',ha='center',fontsize=11)
            fig.text(.68,.83,'Right hemisphere',ha='center',fontsize=11)
            fig.text(.025,.035,'M2 + baseline age squared + baseline centre. Grey: '+('unmeasured/failed parcels.' if tag=='all' else 'unmeasured or FDR >= 0.05.')+'  Bilateral temporal poles unmeasured.',fontsize=10)
            save(fig,out,f'a_thickness_{tag}')
        annotate_cortical_panels(d,out)


    def annotate_cortical_panels(d,out):
        """Locate leaders inside visible DK parcels using unlit binary atlas-mask renders.

        Identical surfaces, cameras, viewport and screenshot dimensions are used for
        beta and atlas renders. Only label layout is manual; anatomical endpoints
        are computed from the rendered parcel masks, never guessed from the image.
        """
        from enigmatoolbox.plotting import plot_cortical
        from enigmatoolbox.utils.parcellation import parcel_to_surface
        from scipy.ndimage import distance_transform_edt, label
        renders=out/'cache/surface_renders'
        masks=renders/'annotation_masks';masks.mkdir(exist_ok=True)
        viewports=None
        # View order as displayed by the ENIGMA fsa5 renderer: LH outer/inner, RH inner/outer.
        labels=[
            (0,'lh_caudalmiddlefrontal','Caudal middle\nfrontal gyrus',.24,65),
            (0,'lh_precentral','Precentral gyrus',.72,40),
            (0,'lh_postcentral','Postcentral\ngyrus',.90,375),
            (0,'lh_middletemporal','Middle temporal\ngyrus',.40,445),
            (1,'lh_superiorfrontal','Superior frontal\ngyrus',.72,65),
            (1,'lh_precuneus','Precuneus',.25,65),
            (1,'lh_entorhinal','Entorhinal cortex',.56,445),
            (2,'rh_superiorfrontal','Superior frontal\ngyrus',.25,65),
            (2,'rh_isthmuscingulate','Isthmus of\ncingulate cortex',.75,445),
            (3,'rh_precentral','Precentral gyrus',.30,65),
            (3,'rh_parsopercularis','Pars opercularis',.80,65),
            (3,'rh_superiortemporal','Superior temporal\ngyrus',.55,445),
        ]
        audit=[]
        for view,key,text,tx,ty in labels:
            if key not in d.index or d.loc[key,'status']!='ok' or d.loc[key,'FDR_family']>=.05:
                continue
            values=np.zeros(68);values[KEYS.index(key)]=1
            p=plot_cortical(array_name=parcel_to_surface(values,'aparc_fsa5',fill=np.nan),
                surface_name='fsa5',cmap='Reds',color_bar=True,color_range=(0,1),
                nan_color=(.87,.87,.87,1),size=(1600,500),interactive=False,
                offscreen=True,return_plotter=True,actor__lighting=False,
                mapper__interpolateScalarsBeforeMapping=False)
            path=masks/(key+'.png');p.screenshot(str(path),transparent_bg=False)
            if viewports is None:
                viewports=sorted([ren.VTKObject.GetViewport() for group in p.renderers.values()
                    for ren in group if ren.VTKObject.GetActors().GetNumberOfItems()],key=lambda v:v[0])
                if len(viewports)!=4:raise ValueError('Expected four cortical viewports')
            p.close()
            rgb=plt.imread(path)[:,:,:3];height,width=rgb.shape[:2]
            xmin,_,xmax,_=viewports[view];lo=int(round(xmin*width));hi=int(round(xmax*width))
            crop=rgb[:,lo:hi]
            mask=(crop[:,:,0]>.3)&(crop[:,:,1]<.2)&(crop[:,:,2]<.25)
            components,n=label(mask)
            if not n:raise ValueError(f'No visible atlas pixels for {key} in view {view}')
            counts=np.bincount(components.ravel());counts[0]=0
            mask=components==counts.argmax()
            if mask.sum()<15:raise ValueError(f'Insufficient visible parcel area for {key}')
            dist=distance_transform_edt(mask)
            y,x=np.unravel_index(np.argmax(dist),dist.shape);x+=lo
            audit.append(dict(atlas_key=key,field_id=d.loc[key,'field_id'],view=view,
                label=text.replace('\n',' '),display_label=text,target_x=int(x),target_y=int(y),
                text_x=lo+tx*(hi-lo),text_y=ty,visible_pixels=int(mask.sum()),
                beta=d.loc[key,'beta'],FDR_family=d.loc[key,'FDR_family']))
        pd.DataFrame(audit).to_csv(out/'panel_a_annotation_coordinates.csv',index=False)
        for tag in ('all','FDR05'):
            source=renders/f'a_thickness_{tag}.png'
            if not source.exists():continue
            fig=plt.figure(figsize=(20,7))
            ax=fig.add_axes([.01,.10,.98,.80]);ax.imshow(plt.imread(source));ax.axis('off')
            for item in audit:
                ax.annotate(item['display_label'],xy=(item['target_x'],item['target_y']),
                    xytext=(item['text_x'],item['text_y']),ha='center',va='center',fontsize=10,
                    arrowprops=dict(arrowstyle='-',color='#333333',lw=1.05,shrinkA=3,shrinkB=0),
                    bbox=dict(facecolor='white',edgecolor='none',pad=.7),annotation_clip=False)
            for view,caption in enumerate(['Left lateral','Left medial','Right medial','Right lateral']):
                xmin,_,xmax,_=viewports[view]
                ax.text((xmin+xmax)*width/2,495,caption,ha='center',fontsize=10,color='.3')
            fig.text(.03,.95,'a  DOSE and cortical thickness | '+('all measured regions' if tag=='all' else 'family FDR < 0.05'),fontsize=17)
            fig.text(.885,.865,'Standardized beta',fontsize=10)
            fig.text(.03,.045,'Labels identify selected FDR-significant DK regions; leader endpoints lie within visible atlas parcels.',fontsize=10)
            fig.text(.03,.018,'M2 + baseline age squared + baseline centre. Grey: '+('unmeasured regions.' if tag=='all' else 'unmeasured or FDR >= 0.05.')+' Bilateral temporal poles unmeasured.',fontsize=10)
            save(fig,out,f'a_thickness_{tag}_annotated')

    def tract_forest(r,out):
        fig,axs=plt.subplots(1,2,figsize=(13,10),sharey=True)
        for ax,fam in zip(axs,['FA27','MD27']):
            d=r[(r.exposure=='dose')&(r.family==fam)].copy().sort_values('field_id')
            y=np.arange(len(d));ok=d.status.eq('ok')
            for j,(_,v) in enumerate(d.iterrows()):
                if not ok.iloc[j]:continue
                color='#b2182b' if v.beta>=0 else '#2166ac'
                ax.errorbar(v.beta,j,xerr=[[v.beta-v.LCI],[v.HCI-v.beta]],fmt='o',
                            color=color,alpha=1 if v.FDR_family<.05 else .3,capsize=2)
                ax.annotate(stars(v.FDR_family),(v.HCI,j),xytext=(3,0),textcoords='offset points',va='center')
            ax.set_yticks(y,d.atlas_key);ax.axvline(0,color='.6',lw=1)
            ax.set_title(fam.replace('27','')+(' | no FDR-significant associations' if not (d.FDR_family<.05).any() else ''));ax.set_xlabel('Standardized beta (95% CI)')
        axs[0].invert_yaxis()
        fig.suptitle('b  White-matter tract associations | coefficient view',x=.05,ha='left')
        fig.text(.05,.015,'Opaque points: family FDR < 0.05. Faded points: not significant. Missing/failed tests have no point.',fontsize=9)
        fig.tight_layout(rect=(0,.04,1,.96));save(fig,out,'b_tract_coefficients')

    def tract_anatomy(r,out,manifest):
        """Render actual supplied MNI tract masks, not invented streamlines.
        All masks must be registered to the same MNI world coordinates.
        Probabilities must be 0..1; binary masks also accepted; thresholds explicit.
        """
        import nibabel as nib
        from skimage.measure import marching_cubes
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        m=pd.read_csv(manifest).set_index('atlas_key')
        if m.index.duplicated().any():raise ValueError('Duplicate tract atlas_key')
        selected=r[(r.exposure=='dose') & r.family.isin(['FA27','MD27']) &
                   r.status.eq('ok') & (r.FDR_family<.05)]
        if selected.empty:print('No FDR-significant tracts; anatomical b omitted.');return
        needed=set(selected.atlas_key)
        if not needed.issubset(m.index):raise ValueError('Missing tract masks: '+str(needed-set(m.index)))
        meshes={}
        for key in needed:
            row=m.loc[key]
            if pd.isna(row.mask_path) or not str(row.mask_path).strip():raise ValueError(f'No mask path: {key}')
            p=Path(str(row.mask_path));p=p if p.is_absolute() else Path(manifest).parent/p
            im=nib.load(p);a=im.get_fdata()
            if a.ndim!=3 or not np.isfinite(a).all() or a.min()<0 or a.max()>1:
                raise ValueError(f'{key}: require finite 3D binary / 0..1 probability mask')
            th=float(row.threshold)
            if not 0<th<1:raise ValueError('Threshold must lie in (0,1)')
            vv,ff,_,_=marching_cubes(a,level=th,step_size=2)
            vv=nib.affines.apply_affine(im.affine,vv);meshes[key]=(vv,ff)
        lim=max(float(selected.beta.abs().max()),1e-6);norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim)
        cmap=plt.get_cmap('RdBu_r');fig=plt.figure(figsize=(13,5))
        allv=np.concatenate([v for v,f in meshes.values()]);lo=allv.min(0);hi=allv.max(0)
        for i,fam in enumerate(['FA27','MD27']):
            ax=fig.add_subplot(1,2,i+1,projection='3d')
            for _,v in selected[selected.family==fam].iterrows():
                verts,faces=meshes[v.atlas_key]
                ax.add_collection3d(Poly3DCollection(verts[faces],facecolor=cmap(norm(v.beta)),edgecolor='none',alpha=.85))
            ax.set(xlim=(lo[0],hi[0]),ylim=(lo[1],hi[1]),zlim=(lo[2],hi[2]))
            ax.set_box_aspect(np.maximum(hi-lo,1));ax.view_init(elev=15,azim=180);ax.set_axis_off();ax.set_title(fam.replace('27',''))
        fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=fig.axes,shrink=.6,label='Standardized beta')
        fig.suptitle('b  White-matter associations | FDR < 0.05')
        fig.text(.03,.02,'Standard atlas tract masks in MNI space; not participant tractography.',fontsize=9)
        save(fig,out,'b_tract_anatomy')
        selected.to_csv(out/'panel_b_displayed_results.csv',index=False)

    def panel_c(r,out,selection='all',families=None,prefix='c_components'):
        families=FAMILIES if families is None else families
        c=r[r.exposure_type.eq('component') & r.family.isin(families)].copy()
        if c.empty:raise ValueError('Component results are required for panel c')
        total=r[r.exposure.eq('dose') & r.family.isin(families)]
        ids=[]
        for fam in families:
            d=total[total.family==fam].sort_values('field_id')
            if selection=='dose-fdr':d=d[d.FDR_family<.05]
            ids.extend(d.field_id.tolist())
        if not ids:
            (out/(prefix+'_'+selection+'_empty.txt')).write_text('No DOSE association passes family FDR < 0.05; no selected component columns.',encoding='utf-8')
            return
        c=c[c.field_id.isin(ids)].copy()
        c.loc[~c.status.eq('ok'),['beta','FDR_family']]=np.nan
        ex=list(c.exposure.drop_duplicates())
        # Multiple manageable sheets, retaining every preplanned feature by default.
        for page,start in enumerate(range(0,len(ids),40),1):
            cols=ids[start:start+40]
            b=c.pivot(index='exposure',columns='field_id',values='beta').reindex(index=ex,columns=cols)
            q=c.pivot(index='exposure',columns='field_id',values='FDR_family').reindex(index=ex,columns=cols)
            finite=c.beta[np.isfinite(c.beta)]
            lim=max(float(finite.abs().max()),1e-6) if len(finite) else .01
            fig,ax=plt.subplots(figsize=(max(10,len(cols)*.34+3),max(3,len(ex)*.4+2.5)))
            cm=plt.get_cmap('RdYlBu_r').copy();cm.set_bad('#d7d7d7')
            image=ax.imshow(b.to_numpy(),cmap=cm,vmin=-lim,vmax=lim,aspect='auto')
            labels=total.set_index('field_id').loc[cols]
            ax.set_xticks(np.arange(len(cols)),[f'{v.family}: '+(short_label(v.label) if v.family in ('Volume139','WMH') else str(v.atlas_key)) for _,v in labels.iterrows()],rotation=55,ha='right',fontsize=8)
            ax.set_yticks(np.arange(len(ex)),[s.removesuffix('.dose').replace('_',' ') for s in ex])
            for i in range(len(ex)):
                for j in range(len(cols)):ax.text(j,i,stars(q.iloc[i,j]),ha='center',va='center',fontsize=8)
            ax.set_title(f'{"c  " if prefix=="c_components" else ""}Component score associations | {selection} features | page {page}',loc='left',pad=15)
            fig.colorbar(image,ax=ax,shrink=.75,label='Standardized beta')
            fig.text(.01,.005,'* q<0.05; ** q<0.01; *** q<0.001. BH within exposure and modality. Grey: unavailable/failed.',fontsize=8)
            fig.tight_layout(rect=(0,.04,1,1));save(fig,out,f'{prefix}_{selection}_{page:02d}')
        c.to_csv(out/f'{prefix}_{selection}_source_results.csv',index=False)

    def supplementary_forests(r,out):
        """Every volume and WMH estimate; categories change display only, not FDR."""
        d=r[r.exposure.eq('dose') & r.family.isin(['Volume139','WMH'])]
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
                ax.set_xlabel('Standardized beta (95% CI)')
                ax.set_title(f'DOSE | {fam}: {category} | page {page}',loc='left')
                fig.text(.015,.015,'BH across all 139 volumes (WMH: one test). Opaque: q<0.05; faded: not significant.',fontsize=9)
                fig.tight_layout(rect=(0,.045,1,1))
                safe=''.join(c if c.isalnum() else '_' for c in str(category))
                save(fig,out,f's_DOSE_{fam}_{safe}_{page:02d}')

    def organize_outputs(out):
        """Keep publication figures apart from supplementary plots and audit data."""
        for folder in ('main','supplementary','audit'):(out/folder).mkdir(exist_ok=True)
        for p in list(out.iterdir()):
            if not p.is_file():continue
            if p.suffix in ('.png','.pdf'):
                # Annotated exports supersede duplicate unannotated composed panels.
                if p.stem in ('a_thickness_all','a_thickness_FDR05'):
                    if (out/(p.stem+'_annotated'+p.suffix)).exists():p.unlink();continue
                main=p.stem=='Figure4_MRI'
                dest=out/('main' if main else 'supplementary')/p.name
            else:dest=out/'audit'/p.name
            p.replace(dest)
        for p in list((out/'main').iterdir()):
            if p.is_file() and p.stem!='Figure4_MRI':p.replace(out/'supplementary'/p.name)

    def figure_manifest_and_combine(r,out):
        total=r[r.exposure.eq('dose')].copy()
        selected=total.family.isin(FAMILIES)&total.status.eq('ok')&(total.FDR_family<.05)
        coverage=total[['field_id','family','category','label','beta','FDR_family']].copy()
        coverage['main_a']=total.family.eq('ThicknessDK')&(total.FDR_family<.05)
        coverage['main_b']=total.family.isin(['FA27','MD27'])
        coverage['main_c_component_column']=selected
        coverage['volume_anatomical_panel']=np.where(total.family.eq('Volume139'),total.category,'')
        coverage['volume_anatomical_coloured']=total.family.eq('Volume139') & total.status.eq('ok') & (total.FDR_family<.05)
        coverage['supplement']=np.where(total.family.isin(['Volume139','WMH']),
            'total forest and component heatmap','complete component heatmap')
        coverage.to_csv(out/'figure_coverage.csv',index=False)
        total[selected].to_csv(out/'main_selected_regions.csv',index=False)
        plan=r.groupby(['exposure','family'],sort=False).size().reset_index(name='planned_tests')
        plan['method']='BH';plan['scope']='one score x one full MRI family; both hemispheres'
        plan.to_csv(out/'FDR_plan.csv',index=False)
        sources=[out/'a_thickness_FDR05_annotated.png',out/'b_tract_coefficients.png',
                 out/'c_components_dose-fdr_01.png']
        for bname in ('b_tract_coefficients','b_tract_anatomy'):
            sources[1]=out/(bname+'.png')
            if not all(p.exists() for p in sources):continue
            imgs=[plt.imread(p) for p in sources]
            heights=[im.shape[0]/im.shape[1]*16 for im in imgs]
            fig,axs=plt.subplots(3,1,figsize=(16,sum(heights)+.4),gridspec_kw={'height_ratios':heights})
            for ax,im in zip(axs,imgs):ax.imshow(im);ax.axis('off')
            fig.subplots_adjust(left=0,right=1,bottom=0,top=1,hspace=.025)
            save(fig,out,'Figure_abc_'+('anatomical' if bname=='b_tract_anatomy' else 'coefficient_view'))
        (out/'figure_caption.txt').write_text(
            'Associations of DOSE and component scores with MRI. M2 plus baseline age squared and centre; '
            'ICV additionally adjusted for volumes and WMH. Exposure and outcome directly standardized within each model. '
            'BH performed separately for each score within 66 thickness, 27 FA, 27 MD, 139 volume and 1 WMH measures. '
            'a: significant thickness parcels with selected anatomical labels; temporal poles unmeasured. '
            'b: UKB group-average tract surfaces in the anatomical version; FA/MD forest plots in the coefficient version. '
            'Surfaces are reference anatomy, not participant-specific streamlines or voxelwise significance. '
            'Group-average resources: University of Oxford / UK Biobank, CC BY 4.0. '
            'c: component coefficients for metrics selected by total-score family FDR < 0.05; '
            'stars use unchanged component-specific full-family FDR. '
            'Positive beta red, negative beta blue. Combined PDF embeds raster panels; individual b/c PDFs retain vector text.',encoding='utf-8')

    def main():
        p=argparse.ArgumentParser(description=__doc__)
        default_out=ROOT/'ukb/05mechanism/04mri'
        p.add_argument('--results',default=str(default_out/'01_all_linear_associations.csv'))
        p.add_argument('--out',default=str(default_out/'figures'))
        p.add_argument('--tract-manifest')
        p.add_argument('--publication-only',action='store_true',help='Assemble Figure4 MRI from current results and verified render caches')
        p.add_argument('--skip-tract-3d',action='store_true',help='Explicitly omit default UKB group-average tract rendering')
        p.add_argument('--skip-volume-3d',action='store_true',help='Explicitly omit volume atlas surfaces')
        p.add_argument('--skip-surface',action='store_true',help='Run b forest and c without ENIGMA')
        args=p.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
        r=pd.read_csv(args.results,dtype={'field_id':str})
        required={'exposure','exposure_type','field_id','family','atlas_key','label','status','beta','LCI','HCI','FDR_family'}
        if not required.issubset(r.columns):raise ValueError('Missing result columns: '+str(required-set(r.columns)))
        if 'family' not in r:r['family']='Volume139'
        if 'FDR_family' not in r:r['FDR_family']=r['FDR_within_exposure_139']
        if r.duplicated(['exposure','field_id']).any():raise ValueError('Duplicate exposure/field pairs')
        if args.publication_only:
            generate_publication = common.mri_publication().generate
            if not (out/'cache/tract_renders/FA27.png').exists():
                generate = common.mri_tracts().generate
                generate(r,out,common.mri_atlas_dir()/'ukb_tracts27')
            generate_publication(r,out);organize_outputs(out);return
        total=r[r.exposure.eq('dose') & r.status.eq('ok')].copy()
        total['FDR_significant']=total.FDR_family<.05
        total.groupby('family').agg(N_regions=('beta','size'),beta_min=('beta','min'),
            beta_max=('beta','max'),FDR05=('FDR_significant','sum')).to_csv(out/'DOSE_beta_ranges.csv')
        observed_volume_preview(r,out)
        missing=set(FAMILIES)-set(r.family)
        if missing:
            msg='a-c not generated: missing actual regression results for '+', '.join(sorted(missing))
            (out/'MISSING_INPUTS.txt').write_text(msg+'\nAlso need component rows for c.\n',encoding='utf-8')
            print(msg);return
        errors=[];notes=[]
        if not args.skip_volume_3d:
            try:
                generate_volume = common.mri_volume().generate
                generate_volume(r,out,common.mri_atlas_dir())
                notes.append('Volume anatomy: Harvard-Oxford cortex/subcortex and MNI-FNIRT cerebellum; reference anatomy, retaining BH across 139 volumes.')
            except Exception as e:errors.append('Volume anatomy: '+str(e))
        if not args.skip_surface:
            try:panel_a(r,out)
            except Exception as e:errors.append('a: '+str(e))
        tract_forest(r,out)
        if args.tract_manifest:
            try:tract_anatomy(r,out,args.tract_manifest)
            except Exception as e:errors.append('b anatomy: '+str(e))
        elif not args.skip_tract_3d:
            try:
                generate = common.mri_tracts().generate
                generate(r,out,common.mri_atlas_dir()/'ukb_tracts27')
                notes.append('b anatomy: UKB group-average 27-tract atlas; intensity threshold 1000; smoothing sigma 0.7 voxels. Not individual streamlines.')
            except Exception as e:errors.append('b anatomy: '+str(e))
        else:notes.append('b anatomy explicitly skipped; coefficient view generated')
        try:
            panel_c(r,out,'dose-fdr')
            panel_c(r,out,'all')
            panel_c(r,out,'all',['Volume139'],'s_volume_components')
            panel_c(r,out,'all',['WMH'],'s_WMH_components')
            supplementary_forests(r,out)
        except Exception as e:errors.append('c: '+str(e))
        (out/'plot_status.txt').write_text('\n'.join(errors+notes) if errors or notes else 'All requested panels generated.',encoding='utf-8')
        for e in errors+notes:print(e)
        if errors:raise RuntimeError('Required plot failed; see plot_status.txt')
        figure_manifest_and_combine(r,out)
        if not args.skip_surface and not args.skip_tract_3d and not args.tract_manifest:
            generate_publication = common.mri_publication().generate
            generate_publication(r,out)
        organize_outputs(out)

    main()


# ===========================================================================
# mechanism / mri_check
def run_mri_check():
    """Check completed MRI output tables and required plot artifacts."""
    from pathlib import Path
    import sys
    output_dir = common.output_dir
    import json
    ROOT = output_dir()
    sys.path.insert(0,str(ROOT/'.tools/mri_plot_packages'))
    import numpy as np
    import pandas as pd


    def bh(p):
        p=np.asarray(p,float)
        order=np.argsort(p)
        ranked=p[order]
        q=np.minimum.accumulate((ranked*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
        out=np.empty_like(q);out[order]=np.minimum(q,1)
        return out


    def main():
        out=ROOT/'ukb/05mechanism/04mri'
        r=pd.read_csv(out/'01_all_linear_associations.csv')
        assert len(r)==2600 and r.exposure.nunique()==10
        assert not r.duplicated(['exposure','field_id']).any()
        assert r.status.eq('ok').all(), 'Investigate failed models before accepting this run'
        assert r.age_center_source.eq('baseline').all()
        assert r['transform'].eq('none').all()
        assert np.isfinite(r[['beta','SE','P','LCI','HCI','N']]).all().all()
        assert ((r.LCI<r.beta)&(r.beta<r.HCI)&(r.SE>0)&(r.N>=100)).all()
        assert np.allclose(r.beta/r.exposure_SD,r.beta_per_score_point)
        for (_,fam),g in r.groupby(['exposure','family']):
            assert len(g)=={'Volume139':139,'ThicknessDK':66,'FA27':27,'MD27':27,'WMH':1}[fam]
            assert np.allclose(bh(g.P),g.FDR_family,atol=1e-12)
        for _,g in r.groupby('exposure'):
            assert np.allclose(bh(g.P),g.FDR_exposure_all,atol=1e-12)
        for _,g in r[r.exposure_type=='component'].groupby('family'):
            assert np.allclose(bh(g.P),g.FDR_components_family,atol=1e-12)
        assert np.allclose(bh(r.P),r.FDR_all_tests,atol=1e-12)
        counts=pd.read_csv(out/'00_observed_MRI_counts.csv').set_index('field_id')
        assert (r.N.to_numpy()<=counts.loc[r.field_id,'N_observed'].to_numpy()).all()
        cov=pd.read_csv(out/'00_covariates_used.csv')
        assert not cov.variable.str.contains('alcohol',case=False).any()
        assert {'baseline_age','baseline_age_squared','baseline_center'} <= set(cov.variable)
        assert r.adjust_icv.eq(r.family.isin(['Volume139','WMH'])).all()
        figs=out/'figures'
        for name in ['supplementary/a_thickness_FDR05_annotated.png','supplementary/b_tract_coefficients.png',
                     'supplementary/c_components_dose-fdr_01.png','supplementary/a_thickness_all_annotated.png',
                     'supplementary/c_components_all_01.png','supplementary/c_components_all_02.png',
                     'supplementary/c_components_all_03.png','supplementary/DOSE_volume_FDR05.png',
                     'supplementary/s_WMH_components_all_01.png']:
            assert (figs/name).stat().st_size>1000, name
        for name,n in [('c_components_all',1080),('s_volume_components_all',1251),('s_WMH_components_all',9)]:
            source=pd.read_csv(figs/'audit'/f'{name}_source_results.csv')
            assert len(source)==n
            merged=source.merge(r,on=['exposure','field_id'],suffixes=('_plot','_full'),validate='one_to_one')
            assert len(merged)==n and np.allclose(merged.FDR_family_plot,merged.FDR_family_full)
        selected=pd.read_csv(figs/'audit/c_components_dose-fdr_source_results.csv')
        targets=r[(r.exposure=='dose') & r.family.isin(['ThicknessDK','FA27','MD27']) & (r.FDR_family<.05)]
        assert set(selected.field_id)==set(targets.field_id) and len(selected)==len(targets)*9
        merged=selected.merge(r,on=['exposure','field_id'],suffixes=('_plot','_full'),validate='one_to_one')
        assert np.allclose(merged.FDR_family_plot,merged.FDR_family_full)
        annotations=pd.read_csv(figs/'audit/panel_a_annotation_coordinates.csv')
        assert (annotations.FDR_family<.05).all() and (annotations.visible_pixels>=15).all()
        mapping=pd.read_csv(figs/'audit/panel_a_mapping.csv')
        assert len(mapping)==68
        assert set(mapping.loc[mapping.display_status=='missing_or_failed','atlas_key'])=={'lh_temporalpole','rh_temporalpole'}
        tract=pd.read_csv(figs/'audit/displayed_tract_results.csv')
        expected_tract=r[r.exposure.eq('dose') & r.family.isin(['FA27','MD27']) & r.FDR_family.lt(.05)]
        assert set(zip(tract.family,tract.atlas_key))==set(zip(expected_tract.family,expected_tract.atlas_key))
        merged=tract.merge(expected_tract,on=['family','atlas_key'],suffixes=('_plot','_full'),validate='one_to_one')
        assert np.allclose(merged.beta_plot,merged.beta_full) and np.allclose(merged.FDR_family_plot,merged.FDR_family_full)
        atlas=pd.read_csv(figs/'audit/atlas_mapping_audit.csv')
        assert len(atlas)==27 and set(atlas.nifti_index_zero_based)==set(range(1,28))
        assert (atlas.loc[atlas.atlas_key.str.endswith('_l'),'centroid_x']<0).all()
        assert (atlas.loc[atlas.atlas_key.str.endswith('_r'),'centroid_x']>0).all()
        for name in ['supplementary/b_tract_anatomy.png','supplementary/Figure_abc_anatomical.pdf',
                     'supplementary/FA_three_views.png','supplementary/FA_tract_groups.png']:
            assert (figs/name).stat().st_size>1000,name
        summary={'status':'PASS','models':len(r),'exposures':int(r.exposure.nunique()),
                 'outcomes':int(r.field_id.nunique()),'N_min':int(r.N.min()),'N_max':int(r.N.max()),
                 'age_center_source':'baseline','BH_checks':'all four correction columns independently reproduced',
                 'missing_surface_parcels':['lh_temporalpole','rh_temporalpole']}
        summary['tract_3d']='PASS: 27 atlas indices and hemispheres checked; displayed coefficients match current results'
        vol=pd.read_csv(figs/'audit/volume_mapping_and_results.csv')
        actual=r[r.exposure.eq('dose') & r.family.eq('Volume139')]
        assert len(vol)==139 and vol.field_id.nunique()==139 and (vol.n_voxels>0).all()
        merged=vol.merge(actual,on='field_id',suffixes=('_plot','_full'),validate='one_to_one')
        assert len(merged)==139
        for name in ['beta','P','FDR_family']:
            assert np.allclose(merged[name+'_plot'],merged[name+'_full'])
        assert (vol.loc[vol.side.eq('Left'),'centroid_x']<0).all()
        assert (vol.loc[vol.side.eq('Right'),'centroid_x']>0).all()
        coverage=pd.read_csv(figs/'audit/cortical_projection_coverage.csv')
        assert len(coverage)==96 and coverage.field_id.nunique()==96 and (coverage.surface_faces>0).all()
        assert set(coverage.field_id)==set(vol.loc[vol.category.eq('Cortical'),'field_id'])
        for ext in ('png','pdf'):assert (figs/'supplementary'/('volume_brain_maps.'+ext)).stat().st_size>1000
        summary['volume_3d']='PASS: 139 unique atlas mappings; 96 cortical fields projected; coefficients and FDR unchanged'
        for ext in ('png','pdf'):assert (figs/'main'/('Figure4_MRI.'+ext)).stat().st_size>1000
        focus=pd.read_csv(figs/'audit/Figure4_subcortical_source.csv')
        assert len(focus)==4 and set(focus.label)=={'Left hippocampus','Right hippocampus','Left thalamus','Right thalamus'}
        selected=r[r.exposure.eq('dose') & r.field_id.isin(focus.field_id)]
        merged=focus.merge(selected,on='field_id',suffixes=('_plot','_full'),validate='one_to_one')
        for col in ['beta','FDR_family']:assert np.allclose(merged[col+'_plot'],merged[col+'_full'])
        point=pd.read_csv(figs/'audit/Figure4_WMH_point_source.csv')
        assert len(point)==1
        wmh=r[r.exposure.eq('dose') & r.family.eq('WMH')].iloc[0]
        for col in ['N','beta','LCI','HCI','P','FDR_family']:
            assert np.isclose(point.iloc[0][col],wmh[col],rtol=1e-12,atol=1e-14)
        axis=pd.read_csv(figs/'audit/Figure4_WMH_axis_audit.csv').iloc[0]
        assert axis.xmin < min(0.,wmh.LCI) and axis.xmax > max(0.,wmh.HCI)
        assert np.allclose(axis[['beta','LCI','HCI']].to_numpy(dtype=float),
                           wmh[['beta','LCI','HCI']].to_numpy(dtype=float))
        sizes=pd.read_csv(figs/'audit/Figure4_brain_display_sizes.csv')
        assert len(sizes)==6 and (sizes.brain_width_inches>0).all()
        assert np.allclose(sizes.loc[sizes.panel.isin(['b','c']),'brain_width_inches'],
                           sizes.loc[sizes.panel.eq('a'),'brain_width_inches'].median())
        for ext in ('png','pdf'):
            assert (figs/'supplementary'/('d_WMH_point_CI.'+ext)).stat().st_size>1000
        assert ((focus.x>0)&(focus.x<1400)&(focus.y>0)&(focus.y<1100)).all()
        summary['publication_figure']='PASS: four labelled regions and WMH point/CI match original results; b/c brain widths match median a width'
        (out/'verification.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
        print(json.dumps(summary,indent=2))


    main()


STEPS = {
    "mri_prepare": run_mri_prepare,
    "mri_plot": run_mri_plot,
    "mri_check": run_mri_check,
}

if __name__ == "__main__":
    common.run_cli(STEPS)
