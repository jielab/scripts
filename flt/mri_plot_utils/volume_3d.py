#!/usr/bin/env python3
"""FLT regional volume associations; retain the existing 139-test FDR family."""
import argparse
import json
import hashlib
import os
import re
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage.measure import marching_cubes

ASSETS = {
 'Cortical': ('HarvardOxford-Cortical.xml', 'HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz'),
 'Subcortical': ('HarvardOxford-Subcortical.xml', 'HarvardOxford/HarvardOxford-sub-maxprob-thr25-2mm.nii.gz'),
 'Cerebellar': ('Cerebellum_MNIfnirt.xml', 'Cerebellum/Cerebellum-MNIfnirt-maxprob-thr25-2mm.nii.gz')}
MIRROR='https://raw.githubusercontent.com/fithisux/FSL/master/data/atlases/'

def norm(s):
 return re.sub(r'[^a-z0-9]', '', str(s).lower())

def target_name(row):
 s=re.sub(r'^Volume of (grey matter in )?', '', row.label)
 s=re.sub(r'\s*\((left|right|vermis)\)$', '', s, flags=re.I)
 if row.category=='Cortical': return s
 if row.category=='Brainstem': return 'Brain-Stem'
 if row.category=='Subcortical': return f'{row.side} {s}'
 if row.category=='Cerebellar': return f'{row.side} {s.replace(" Cerebellum", "")}'
 raise ValueError(f'Unknown category: {row.category}')

def read_results(path, exposure, model, qcol):
 d=path.copy() if isinstance(path,pd.DataFrame) else pd.read_csv(path)
 required={'exposure','model','family','field_id','category','label','side','beta','P',qcol}
 if not required.issubset(d): raise ValueError(f'Missing columns: {required-set(d)}')
 d=d[(d.exposure==exposure)&(d.model==model)&(d.family=='Volume139')].copy()
 if len(d)!=139 or d.field_id.nunique()!=139:
  raise ValueError(f'Expected 139 unique Volume139 rows, found {len(d)}; do not silently select a subset.')
 for c in ['beta','P',qcol]: d[c]=pd.to_numeric(d[c],errors='raise')
 if not np.isfinite(d[['beta','P',qcol]].to_numpy()).all(): raise ValueError('Non-finite results; resolve failed models before plotting.')
 if not d[['P',qcol]].apply(lambda x:x.between(0,1).all()).all(): raise ValueError('Invalid P/q values.')
 # Audit only, never replace supplied q values.
 p=d.P.to_numpy(); order=np.argsort(p); q=np.empty(len(p)); q[order]=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(0,1)
 if qcol=='FDR_family' and not np.allclose(q,d[qcol],atol=1e-8,rtol=1e-5):
  raise ValueError('FDR_family differs from BH across these 139 tests. Check the analysis family.')
 return d

def load_atlases(root, download):
 ans={}
 for group,(xml,nii) in ASSETS.items():
  for rel in (xml,nii):
   p=root/Path(rel).name
   if not p.exists():
    if not download: raise FileNotFoundError(f'{p}\nSet --atlas-dir to FSL data/atlases, or use --download-atlases (FSL mirror).')
    p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_name(p.name+'.partial')
    print('Download:',rel,flush=True)
    urllib.request.urlretrieve(MIRROR+rel,tmp); tmp.replace(p)
  tree=ET.parse(root/xml)
  offset=1 if tree.findtext('./header/type').strip()=='Probabilistic' else 0
  labels={norm(e.text):(int(e.attrib['index'])+offset,e.text.strip()) for e in tree.findall('./data/label')}
  img=nib.load(root/Path(nii).name); data=np.asarray(img.dataobj)
  if data.ndim!=3 or not np.allclose(data,np.round(data)): raise ValueError('Expected discrete 3D atlas')
  ans[group]=(img,data.astype(np.int16),labels)
 return ans

def make_mapping(d, atlases):
 rows=[]
 for r in d.itertuples():
  group='Subcortical' if r.category=='Brainstem' else r.category
  img,data,labels=atlases[group]; name=target_name(r)
  # Explicit spelling aliases only; no fuzzy matching.
  key=norm(name)
  aliases={'leftaccumbens':'leftaccumbens','rightaccumbens':'rightaccumbens'}
  key=aliases.get(key,key)
  if key not in labels: raise ValueError(f'No exact atlas match: {r.field_id}: {name}; candidates={list(labels)}')
  code,matched=labels[key]; mask=data==code
  if group=='Cortical':
   ijk=np.argwhere(mask); x=nib.affines.apply_affine(img.affine,ijk)[:,0]
   keep=x<0 if r.side=='Left' else x>0
   mask[:]=False; ijk=ijk[keep]; mask[tuple(ijk.T)]=True
  if not mask.any(): raise ValueError(f'Empty ROI: {r.field_id}')
  centroid=nib.affines.apply_affine(img.affine,np.argwhere(mask)).mean(0)
  if r.side=='Left' and centroid[0]>=0 or r.side=='Right' and centroid[0]<=0:raise ValueError(f'Hemisphere mismatch: {r.field_id}')
  rows.append(dict(centroid_x=float(centroid[0]),field_id=r.field_id,atlas_group=group,atlas_label=matched,atlas_code=code,n_voxels=int(mask.sum()),side=r.side))
 return d.merge(pd.DataFrame(rows),on=['field_id','side'],validate='one_to_one')

def roi_mesh(img, data, code):
 mask=data==code
 inds=np.argwhere(mask); lo=np.maximum(inds.min(0)-1,0); hi=np.minimum(inds.max(0)+2,data.shape)
 crop=mask[tuple(slice(a,b) for a,b in zip(lo,hi))]
 # Padding ensures closed meshes even if an ROI touches a crop boundary.
 v,f,_,_=marching_cubes(np.pad(crop,1).astype(float),.5)
 v=nib.affines.apply_affine(img.affine,v+lo-1)
 return v,f

def setup_axis(ax,xyz,elev,azim,title):
 lo=xyz.min(0); hi=xyz.max(0); centre=(hi+lo)/2; radius=max(hi-lo)/2*1.06
 ax.set(xlim=(centre[0]-radius,centre[0]+radius),ylim=(centre[1]-radius,centre[1]+radius),zlim=(centre[2]-radius,centre[2]+radius))
 ax.set_box_aspect((1,1,1),zoom=1.5);ax.set_proj_type('ortho');ax.view_init(elev=elev,azim=azim);ax.set_axis_off();ax.set_title(title,fontsize=11,pad=0)

def collection(v,f,colors):
 c=Poly3DCollection(v[f],facecolors=colors,edgecolors=(0,0,0,0),linewidth=0,rasterized=True,shade=True)
 return c

def render(d,atlases,args,out):
 from nilearn import datasets,surface
 scale=args.vmax or float(np.ceil(np.abs(d.beta).max()*1000)/1000)
 if scale<=0 or scale<np.abs(d.beta).max(): raise ValueError('--vmax must include all 139 beta values.')
 cm=plt.get_cmap('RdBu_r'); nm=Normalize(-scale,scale)
 neutral=(.83,.83,.83,1); background=(.95,.95,.95,1)
 def color(r): return cm(nm(r.beta)) if args.show_all or getattr(r,args.q_col)<args.alpha else neutral
 fs=datasets.fetch_surf_fsaverage('fsaverage5',data_dir=args.cache_dir)
 fig=plt.figure(figsize=(15,11));gs=fig.add_gridspec(3,4,left=.01,right=.99,top=.92,bottom=.10,hspace=.14,wspace=.01)
 cortex=atlases['Cortical'][0]; projected=[]
 for col,(side,view,azim) in enumerate([('Left','lateral',180),('Left','medial',0),('Right','medial',180),('Right','lateral',0)]):
  hemi=side.lower();pial=surface.load_surf_mesh(fs['pial_'+hemi]);mesh=surface.load_surf_mesh(fs['infl_'+hemi])
  labels=surface.vol_to_surf(cortex,pial,inner_mesh=fs['white_'+hemi],interpolation='nearest_most_frequent',depth=[.1,.3,.5,.7,.9])
  labels=np.nan_to_num(labels).astype(int)
  face_labels=labels[mesh.faces];a,b,c=face_labels.T
  # Majority vertex label, stable lowest label for a three-way boundary.
  codes=np.where((a==b)|(a==c),a,np.where(b==c,b,np.min(face_labels,axis=1)))
  colors=np.tile(background,(len(codes),1))
  rows=d[(d.category=='Cortical')&(d.side==side)]
  for r in rows.itertuples():
   colors[codes==r.atlas_code]=color(r)
   if view=='lateral': projected.append({'field_id':r.field_id,'surface_faces':int((codes==r.atlas_code).sum())})
  ax=fig.add_subplot(gs[0,col],projection='3d');ax.add_collection3d(collection(mesh.coordinates,mesh.faces,colors));setup_axis(ax,mesh.coordinates,0,azim,f'{side} {view}')
 pd.DataFrame(projected).to_csv(out/'cortical_projection_coverage.csv',index=False)
 absent=pd.DataFrame(projected).query('surface_faces == 0')
 if len(absent): print('WARNING: cortical ROIs absent from projection:',absent.field_id.tolist(),flush=True)
 for row,group in [(1,'Subcortical'),(2,'Cerebellar')]:
  sub=d[d.atlas_group==group];img,arr,_=atlases[group]
  meshes={r.field_id:roi_mesh(img,arr,r.atlas_code) for r in sub.itertuples()}
  views=[('Left structures',0,180,'Left'),('Right structures',0,0,'Right'),('Posterior',0,-90,None),('Inferior',-90,-90,None)] if row==1 else [('Left lateral',0,180,None),('Right lateral',0,0,None),('Posterior',0,-90,None),('Inferior',-90,-90,None)]
  for col,(title,elev,azim,side) in enumerate(views):
   ax=fig.add_subplot(gs[row,col],projection='3d')
   chosen=sub if side is None else sub[sub.side==side]
   # One collection for all meshes allows depth sorting across structures.
   polys=[];colors=[];points=[]
   for r in chosen.itertuples():
    v,f=meshes[r.field_id];polys.extend(v[f]);colors.extend([color(r)]*len(f));points.append(v)
   coll=Poly3DCollection(polys,facecolors=colors,edgecolors=(0,0,0,0),linewidth=0,rasterized=True,shade=True)
   ax.add_collection3d(coll);setup_axis(ax,np.concatenate(points),elev,azim,title)
 fig.text(.025,.965,'a  Cortical volume',weight='bold',fontsize=14)
 fig.text(.025,.663,'b  Subcortical and brainstem volume',weight='bold',fontsize=14)
 fig.text(.025,.38,'c  Cerebellar volume',weight='bold',fontsize=14)
 cax=fig.add_axes([.35,.072,.30,.016]);fig.colorbar(ScalarMappable(norm=nm,cmap=cm),cax=cax,orientation='horizontal',ticks=np.linspace(-scale,scale,5),label='MRI SD per FLT point (positive = red; negative = blue)')
 selection='All coefficients' if args.show_all else f'Coloured: {args.q_col} < {args.alpha:g}; grey: not significant'
 fig.text(.5,.008,selection+' | Existing FDR retained across all 139 volumes',ha='center',fontsize=10)
 for ext in ['png','pdf']: fig.savefig(out/f'volume_brain_maps.{ext}',dpi=300)
 plt.close(fig)
 # A complete significant-ROI forest plot preserves effects hidden in 3D views.
 sig=d[d[args.q_col]<args.alpha].sort_values(['category','field_id'])
 if len(sig) and {'LCI','HCI'}.issubset(sig):
  fig,ax=plt.subplots(figsize=(10,max(3,.34*len(sig)+1)))
  for i,r in enumerate(sig.itertuples()):
   ax.plot([r.LCI,r.HCI],[i,i],color=cm(nm(r.beta)));ax.plot(r.beta,i,'o',color=cm(nm(r.beta)))
  ax.set_yticks(range(len(sig)),[re.sub(r'^Volume of (grey matter in )?','',s) for s in sig.label]);ax.invert_yaxis();ax.axvline(0,color='.5',lw=.8);ax.set_xlabel('MRI SD per FLT point (95% CI)');ax.set_title(f'{args.exposure}: {len(sig)} / 139 volumes, FDR < {args.alpha:g}')
  fig.tight_layout()
  for ext in ['png','pdf']:fig.savefig(out/f'volume_significant_forest.{ext}',dpi=300)
  plt.close(fig)

def generate(results,out,atlas_dir):
 from types import SimpleNamespace
 out=Path(out);atlas_dir=Path(atlas_dir)
 args=SimpleNamespace(exposure=EXPOSURE,model='covs.all',q_col='FDR_family',alpha=.05,
    vmax=None,show_all=False,cache_dir=str(out/'cache/nilearn'))
 d=read_results(results,args.exposure,args.model,args.q_col)
 atlases=load_atlases(atlas_dir,False);d=make_mapping(d,atlases)
 d.to_csv(out/'volume_mapping_and_results.csv',index=False)
 counts=d.assign(significant=d[args.q_col]<args.alpha).groupby('category').agg(tested=('field_id','size'),significant=('significant','sum'))
 counts.to_csv(out/'volume_summary.csv');print(counts.to_string(),flush=True)
 render(d,atlases,args,out)
 files=[atlas_dir/Path(name).name for pair in ASSETS.values() for name in pair]
 info=dict(vars(args),atlas_dir=str(atlas_dir),
    source='FSL Harvard-Oxford cortical/subcortical and Cerebellum MNI-FNIRT atlases',
    documentation='https://fsl.fmrib.ox.ac.uk/fsl/docs/other/datasets.html',
    subcortical_interpretation='Same-named reference anatomy, NOT participant FIRST segmentations',
    sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
 (out/'volume_run_settings.json').write_text(json.dumps(info,indent=2),encoding='utf-8')
 (out/'volume_figure_caption.txt').write_text(
    'FLT and regional brain volumes. Cortical: Harvard-Oxford labels projected to fsaverage5 inflated surfaces; '
    'subcortical/brainstem and cerebellar: surfaces from standard atlas masks. '
    'FIRST-derived volumes use same-named Harvard-Oxford regions for anatomical reference, not individual segmentation boundaries. '
    'Colours are regression beta (MRI SD per FLT point); grey indicates FDR >= 0.05. '
    'BH correction retained jointly across 139 volumes. Colours represent regional associations, not shape changes or vertexwise tests. '
    'covs.all with ICV adjustment for volumes. '    'Some internal regions may be occluded in 3D; refer to the full significant-volume forest plot.',encoding='utf-8')
