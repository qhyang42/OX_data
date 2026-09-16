#!/usr/bin/env python3
"""Anatomical nine-panel schematic of six bilateral standard-space ROIs."""
from pathlib import Path
from itertools import combinations
import hashlib,json
import nibabel as nib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch

root=Path(__file__).resolve().parents[1];out=root/'ROIs/schematic_figures';out.mkdir(parents=True,exist_ok=True)
refpath=root/'ROIs/CZ_ROIs_MNI_1mm/MNI152_T1_1mm_brain.nii';ref=nib.as_closest_canonical(nib.load(refpath));brain=ref.get_fdata()
names=['TU','AON','PirF','PirT','olfAMG','olfOFC']
colors=['#e69f00','#00b8d9','#e64b35','#8b65d6','#4dbb66','#e66bb2']
inputs={}
for name,stem in [('TU','TU'),('AON','AON'),('PirF','pirF'),('PirT','pirT')]:inputs[name]=[root/f'ROIs/CZ_ROIs_MNI_1mm/olfROI/{h}_{stem}.nii.gz' for h in ['L','R']]
inputs['olfAMG']=[root/f'ROIs/CZ_ROIs_MNI_1mm/amgROI/{h}_{c}.nii.gz' for h in ['L','R'] for c in ['MeA','ACo','CeA','PAC']]
inputs['olfOFC']=[root/'ROIs/OFC_small_MNI152_1mm.nii.gz']
masks=[];records=[]
for name in names:
 mask=np.zeros(ref.shape,bool)
 for path in inputs[name]:
  im=nib.as_closest_canonical(nib.load(path));a=im.get_fdata()
  assert im.shape==ref.shape and np.allclose(im.affine,ref.affine,atol=1e-4) and np.isin(a,[0,1]).all()
  mask|=a>0
 assert mask.any();masks.append(mask)
 records.append(dict(roi=name,voxel_count=int(mask.sum()),sources=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in inputs[name]]))
masks=np.array(masks);union=np.any(masks,axis=0)
# Select triples maximizing representation of each ROI, with at least 4-mm spacing.
# Normalize by each ROI's peak section area, so large masks do not dominate.
selection={}
for axis in [2,1,0]:
 area=np.array([np.sum(m,axis=tuple(i for i in range(3) if i!=axis)) for m in masks])
 candidates=np.where(np.any(area>0,axis=0))[0]
 if axis==0: candidates=candidates[(ref.affine[0,0]*candidates+ref.affine[0,3])<0] # left hemisphere representative sagittal views
 normalized=area/np.max(area,axis=1)[:,None]
 best=None;score=-1
 for triple in combinations(candidates,3):
  if min(np.diff(triple))<4:continue
  sample=normalized[:,triple]
  value=10*np.sum(np.max(sample,axis=1)>.05)+np.sum(np.max(sample,axis=1))+.12*np.sum(sample)
  if value>score:score=value;best=triple
 assert best is not None
 selection[axis]=list(map(int,best))
 print(axis,[float(ref.affine[axis,axis]*i+ref.affine[axis,3]) for i in best], 'ROIs visible',np.any(area[:,best]>0,axis=1).tolist())
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'svg.fonttype':'none'})
fig,axes=plt.subplots(3,3,figsize=(12,10.7),facecolor='white')
fig.subplots_adjust(left=.07,right=.985,top=.89,bottom=.15,hspace=.22,wspace=.06)
planes=[(2,'Axial','z'),(1,'Coronal','y'),(0,'Sagittal · left','x')]
# Fixed anatomical crops per plane based on all six masks plus generous margin.
bounds=np.array(np.where(union));lo=bounds.min(axis=1);hi=bounds.max(axis=1)
for row,(axis,label,coord) in enumerate(planes):
 remaining=[i for i in range(3) if i!=axis]
 for col,k in enumerate(selection[axis]):
  ax=axes[row,col];a=np.take(brain,k,axis=axis).T
  ax.imshow(a,origin='lower',cmap='gray',vmin=0,vmax=np.percentile(brain[brain>0],99),interpolation='bilinear')
  sections=np.array([np.take(m,k,axis=axis).T for m in masks]);total=sections.sum(axis=0)
  rgb=np.einsum('nxy,nc->xyc',sections.astype(float),np.array([to_rgb(c) for c in colors]))/np.maximum(total,1)[...,None]
  rgba=np.concatenate([rgb,(.88*(total>0))[...,None]],axis=-1)
  ax.imshow(rgba,origin='lower',interpolation='nearest')
  for j,section in enumerate(sections):
   if section.any():ax.contour(section,levels=[.5],colors=[colors[j]],linewidths=.75)
  ax.set_xlim(lo[remaining[0]]-18,hi[remaining[0]]+18);ax.set_ylim(lo[remaining[1]]-16,hi[remaining[1]]+16)
  ax.set_axis_off();world=ref.affine[axis,axis]*k+ref.affine[axis,3]
  ax.set_title(f'{coord} = {world:+.0f} mm',fontsize=13,pad=7)
  ax.text(.02,.94,'P' if axis==0 else 'L',transform=ax.transAxes,color='white',fontsize=11,va='top',weight='bold')
  ax.text(.98,.94,'A' if axis==0 else 'R',transform=ax.transAxes,color='white',fontsize=11,va='top',ha='right',weight='bold')
  if col==0:ax.text(-.08,.5,label,transform=ax.transAxes,rotation=90,va='center',ha='right',fontsize=14,weight='bold')
fig.suptitle('Olfactory regions of interest',fontsize=23,weight='bold',y=.97)
fig.text(.5,.925,'Bilateral anatomical masks in MNI152 space',ha='center',fontsize=13,color='#444444')
fig.legend(handles=[Patch(facecolor=c,label=n) for n,c in zip(names,colors)],loc='lower center',bbox_to_anchor=(.5,.078),ncol=6,frameon=False,fontsize=13,handlelength=1.2,columnspacing=1.7)
fig.text(.5,.025,'Axial/coronal: neurological orientation (L = left). Sagittal: left hemisphere; P = posterior, A = anterior.\nolfAMG = MeA + ACo + CeA + PAC; olfOFC = small olfactory OFC mask. Anatomical masks, without functional restriction.\nColors blend at overlapping voxels; colored outlines retain ROI boundaries.',ha='center',fontsize=9,linespacing=1.5)
for ext in ['png','pdf','svg']:
 target=out/f'olfactory_ROIs_MNI_3planes_3slices.{ext}'
 if target.exists():raise FileExistsError(target)
 fig.savefig(target,dpi=300,facecolor='white')
plt.close(fig)
(out/'olfactory_ROIs_MNI_schematic_manifest.json').write_text(json.dumps(dict(reference=str(refpath),rois=records,colors=dict(zip(names,colors)),slices={label:[float(ref.affine[ax,ax]*k+ref.affine[ax,3]) for k in selection[ax]] for ax,label,_ in planes},selection='Three slices per plane maximize ROI-normalized cross-section coverage, minimum spacing 4 mm; sagittal restricted to left hemisphere',resampling='none: canonical reorientation only; all source grids verified identical',overlap_voxels=int(np.sum(masks.sum(axis=0)>1))),indent=2))
