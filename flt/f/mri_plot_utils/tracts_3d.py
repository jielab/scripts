#!/usr/bin/env python3
"""True YGB group-average tract surfaces coloured by observed association beta.
Requires numpy pandas scipy nibabel scikit-image matplotlib vtk (9.4+ for headless Linux).
"""
from pathlib import Path
import argparse, json, hashlib, urllib.request
import numpy as np
import pandas as pd
import nibabel as nib
from scipy.ndimage import gaussian_filter, binary_fill_holes, label
from skimage.measure import marching_cubes
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
KEYS='ar_l ar_r atr_l atr_r cgc_l cgc_r cgh_l cgh_r cst_l cst_r fma fmi ifo_l ifo_r ilf_l ilf_r mcp ml_l ml_r ptr_l ptr_r slf_l slf_r str_l str_r unc_l unc_r'.split()
BASE='https://www.fmrib.ox.ac.uk/%75%6B%62iobank/group_means/%75%6B%62/'
NAMES=['TractsSummed_Tracts1-27.nii.gz','T1_T2FLAIR_fMRI.nii.gz']
MASK_URL='https://raw.githubusercontent.com/Washington-University/HCPpipelines/master/global/templates/MNI152_T1_1mm_brain_mask.nii.gz'

def mesh(a,affine,level,step=2):
    vertices,faces,_,_=marching_cubes(a,level=level,step_size=step)
    vertices=nib.affines.apply_affine(affine,vertices)
    pts=vtk.vtkPoints();pts.SetData(numpy_to_vtk(vertices,deep=True))
    cells=vtk.vtkCellArray()
    cells.ImportLegacyFormat(numpy_to_vtk(np.c_[np.full(len(faces),3),faces].ravel().astype(np.int64),deep=True,array_type=vtk.VTK_ID_TYPE))
    poly=vtk.vtkPolyData();poly.SetPoints(pts);poly.SetPolys(cells)
    smooth=vtk.vtkWindowedSincPolyDataFilter();smooth.SetInputData(poly);smooth.SetNumberOfIterations(15)
    smooth.SetPassBand(.12);smooth.BoundarySmoothingOff();smooth.NormalizeCoordinatesOn();smooth.Update()
    normals=vtk.vtkPolyDataNormals();normals.SetInputConnection(smooth.GetOutputPort());normals.SplittingOff();normals.Update()
    result=vtk.vtkPolyData();result.DeepCopy(normals.GetOutput());return result

def actor(poly,color,opacity=1):
    mapper=vtk.vtkPolyDataMapper();mapper.SetInputData(poly);mapper.ScalarVisibilityOff()
    a=vtk.vtkActor();a.SetMapper(mapper);p=a.GetProperty();p.SetColor(*color[:3]);p.SetOpacity(opacity)
    p.SetAmbient(.25);p.SetDiffuse(.75);p.SetSpecular(.16);p.SetSpecularPower(20)
    return a

BETA_LIMIT=.02
def render(brain,meshes,rows,direction,size=(1200,1000)):
    win=vtk.vtkRenderWindow();win.SetOffScreenRendering(1);win.SetSize(*size);win.SetMultiSamples(0);win.SetAlphaBitPlanes(1)
    ren=vtk.vtkRenderer();ren.SetBackground(1,1,1);win.AddRenderer(ren)
    ren.SetUseDepthPeeling(1);ren.SetMaximumNumberOfPeels(100);ren.SetOcclusionRatio(.1)
    ren.AddActor(actor(brain,(.76,.79,.82),.09))
    for _,r in rows.iterrows():
        color=TRACT_CMAP(TRACT_NORM(r.beta))
        ren.AddActor(actor(meshes[r.atlas_key],color,.98))
    cam=ren.GetActiveCamera();center=np.array([0,-18,10]);direction=np.asarray(direction,float)
    cam.SetFocalPoint(*center);cam.SetPosition(*(center+direction/np.linalg.norm(direction)*450))
    cam.SetViewUp(*( (0,1,0) if abs(direction[2])>abs(direction[0])+abs(direction[1]) else (0,0,1)))
    cam.ParallelProjectionOn();cam.SetParallelScale(105);ren.ResetCameraClippingRange()
    win.Render();f=vtk.vtkWindowToImageFilter();f.SetInput(win);f.SetInputBufferTypeToRGB();f.ReadFrontBufferOff();f.Update()
    im=f.GetOutput();w,h,_=im.GetDimensions();array=vtk_to_numpy(im.GetPointData().GetScalars()).reshape(h,w,3)[::-1].copy()
    win.Finalize();return array

def save(fig,out,name):
    for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext),dpi=300,bbox_inches='tight',facecolor='white')
    plt.close(fig)

def generate(r,out,atlas,threshold=1000,shared_limit=None,render_only=False):
    from types import SimpleNamespace
    args=SimpleNamespace(threshold=threshold)
    out=Path(out);atlas=Path(atlas)
    maskfile=atlas/'MNI152_T1_1mm_brain_mask.nii.gz'
    r=r[r.exposure.eq(EXPOSURE) & r.family.isin(['FA27','MD27'])].copy()
    expected=[line.split()[0] for line in (atlas/'structureList.source.txt').read_text().splitlines() if line.strip()]
    if expected!=KEYS:raise ValueError('AutoPtx source order differs from renderer')
    if r.duplicated(['family','atlas_key']).any():raise ValueError('Duplicate tract results')
    for fam in ['FA27','MD27']:
        h=r[r.family.eq(fam)]
        if len(h)!=27 or set(h.atlas_key)!=set(KEYS):raise ValueError('Require all 27 named tracts')
        from scipy.stats import false_discovery_control
        if not np.allclose(false_discovery_control(h.P),h.FDR_family,atol=1e-10):raise ValueError('FDR mismatch')
    selected=r[r.status.eq('ok') & r.FDR_family.lt(.05)]
    global BETA_LIMIT, TRACT_CMAP, TRACT_NORM
    BETA_LIMIT=max(.02,float(selected.beta.abs().max())) if len(selected) else .02
    if shared_limit is not None:
        if shared_limit <= 0 or (len(selected) and selected.beta.abs().max()>shared_limit): raise ValueError("Invalid shared colour range")
        BETA_LIMIT=shared_limit
    mixed=shared_limit is not None or (selected.beta<0).any()
    TRACT_CMAP=plt.get_cmap('RdBu_r' if mixed else 'Reds')
    TRACT_NORM=TwoSlopeNorm(vmin=-BETA_LIMIT,vcenter=0,vmax=BETA_LIMIT) if mixed else Normalize(0,BETA_LIMIT)
    selected.to_csv(out/'displayed_tract_results.csv',index=False)
    im=nib.load(atlas/NAMES[0]);bg=nib.load(atlas/NAMES[1])
    if im.shape!=(182,218,182,28) or not np.allclose(im.affine,bg.affine):raise ValueError('Unexpected atlas geometry')
    # Volume 0 is the sum, volumes 1..27 follow official AutoPtx order.
    mapping=[];meshes={};manifest=[]
    masks=out/'cache/selected_tract_masks';masks.mkdir(parents=True,exist_ok=True)
    for i,key in enumerate(KEYS,1):
        a=np.asarray(im.dataobj[...,i],dtype=np.float32)
        xyz=nib.affines.apply_affine(im.affine,np.argwhere(a>a.max()*.15)).mean(0)
        if key.endswith('_l') and xyz[0]>=0 or key.endswith('_r') and xyz[0]<=0:raise ValueError('Hemisphere check failed')
        mapping.append(dict(atlas_key=key,nifti_index_zero_based=i,centroid_x=xyz[0],centroid_y=xyz[1],centroid_z=xyz[2],intensity_max=float(a.max())))
        if key in set(selected.atlas_key):
            if not 0<args.threshold<a.max():raise ValueError('Invalid isosurface threshold')
            meshes[key]=mesh(gaussian_filter(a,.7),im.affine,args.threshold)
            scale=float(a.max());path=masks/(key+'.nii.gz');nib.save(nib.Nifti1Image(a/scale,im.affine),path)
            manifest.append(dict(atlas_key=key,mask_path=str(path.relative_to(out)),threshold=args.threshold/scale,space='MNI152',source=BASE+NAMES[0]))
    pd.DataFrame(mapping).to_csv(out/'atlas_mapping_audit.csv',index=False)
    pd.DataFrame(manifest).to_csv(out/'verified_tract_manifest.csv',index=False)
    t1=np.asarray(bg.dataobj[...,0],dtype=np.float32)
    bm=nib.load(maskfile)
    from nibabel.processing import resample_from_to
    bm=resample_from_to(bm,(bg.shape[:3],bg.affine),order=0)
    mask=(bm.get_fdata()>.5) & (t1>45)
    labs,n=label(mask);counts=np.bincount(labs.ravel());counts[0]=0;mask=labs==counts.argmax()
    brain=mesh(gaussian_filter(mask.astype(np.float32),1.0),bg.affine,.5)
    if render_only:
        target=out/'cache/tract_renders';target.mkdir(parents=True,exist_ok=True)
        for fam in ['FA27','MD27']:
            rows=selected[selected.family.eq(fam)]
            plt.imsave(target/(fam+'.png'),render(brain,meshes,rows,(-1,-.35,.18)))
        print('Saved shared-scale tract maps:',out,flush=True)
        return
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42})
    # Main panel b: no colour is shown for nonsignificant MD.
    fig,axs=plt.subplots(1,2,figsize=(13,6.3))
    for ax,fam in zip(axs,['FA27','MD27']):
        rows=selected[selected.family.eq(fam)]
        rendered=render(brain,meshes,rows,(-1,-.35,.18))
        cache=out/'cache/tract_renders';cache.mkdir(parents=True,exist_ok=True)
        plt.imsave(cache/(fam+'.png'),rendered)
        ax.imshow(rendered);ax.axis('off')
        ax.set_title(f'{fam[:2]} | {len(rows)} of 27 tracts',fontsize=16,pad=4)
        if rows.empty:
            ax.text(.5,.52,'No association survives\nFDR correction',transform=ax.transAxes,ha='center',va='center',fontsize=16,
                    bbox=dict(facecolor='white',alpha=.9,edgecolor='none',pad=8))
    fig.text(.025,.97,'b',fontsize=25,fontweight='bold')
    cax=fig.add_axes([.15,.15,.27,.024]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=TRACT_NORM,cmap=TRACT_CMAP),cax=cax,orientation='horizontal')
    cb.set_label('MRI SD per FLT point',fontsize=11);cb.set_ticks([-BETA_LIMIT,0,BETA_LIMIT] if mixed else [0,BETA_LIMIT/2,BETA_LIMIT])
    fig.text(.5,.035,'YGB group-average tract surfaces | Within-modality FDR < 0.05',ha='center',fontsize=11)
    fig.subplots_adjust(left=.02,right=.98,top=.9,bottom=.22,wspace=.02);save(fig,out,'b_tract_anatomy')
    # Additional viewpoints expose overlap hidden by one camera.
    fa=selected[selected.family.eq('FA27')]
    fig,axs=plt.subplots(1,3,figsize=(17,6))
    for ax,direction,title in zip(axs,[(-1,-.35,.18),(1,-.35,.18),(0,0,1)],['Left oblique','Right oblique','Superior']):
        ax.imshow(render(brain,meshes,fa,direction));ax.axis('off');ax.set_title(title,fontsize=15)
    fig.suptitle(f'FLT and white-matter FA | {len(fa)} FDR-significant tracts',fontsize=18)
    fig.text(.5,.045,'Red: positive beta; blue: negative beta. Tract overlap varies with viewpoint.',ha='center')
    fig.subplots_adjust(left=0,right=1,bottom=.12,top=.87,wspace=.01);save(fig,out,'FA_three_views')
    # Separate groups make all six anatomical structures identifiable.
    groups=[(['fma'],'Forceps major'),(['ifo_l','ifo_r'],'Inferior fronto-occipital fasciculus'),
            (['ilf_l','ilf_r'],'Inferior longitudinal fasciculus'),(['ptr_l','ptr_r'],'Posterior thalamic radiation'),
            (['slf_l','slf_r'],'Superior longitudinal fasciculus'),(['unc_l','unc_r'],'Uncinate fasciculus')]
    fig,axs=plt.subplots(2,3,figsize=(15,10))
    for ax,(keys,title) in zip(axs.flat,groups):
        rows=fa[fa.atlas_key.isin(keys)];ax.imshow(render(brain,meshes,rows,(-1,-.45,.25)));ax.axis('off')
        ax.set_title(title,fontsize=12)
        text='; '.join(('L' if v.atlas_key.endswith('_l') else 'R' if v.atlas_key.endswith('_r') else '')+f' beta={v.beta:.4f}' for _,v in rows.iterrows())
        ax.text(.5,.06,text,transform=ax.transAxes,ha='center',fontsize=10)
    fig.suptitle('Significant FA associations by tract group',fontsize=18)
    fig.subplots_adjust(left=.01,right=.99,bottom=.03,top=.92,wspace=.01,hspace=.06);save(fig,out,'FA_tract_groups')
    info={'atlas_source':BASE+NAMES[0],'background_source':BASE+NAMES[1],'background_brain_mask':MASK_URL,
          'license':'CC BY 4.0; University of Oxford / UK Biobank group-average resources',
          'threshold_atlas_intensity':args.threshold,'smoothing_sigma_voxels':.7,
          'geometry':'Marching-cubes surfaces of real group-average tract maps, not individual streamlines',
          'FA_significant':int((selected.family=='FA27').sum()),'MD_significant':int((selected.family=='MD27').sum()),
          'sha256':{name:hashlib.sha256((atlas/name).read_bytes()).hexdigest() for name in NAMES}}
    (out/'tract_provenance.json').write_text(json.dumps(info,indent=2))
    print('Finished real atlas-based tract rendering.',flush=True)

