#!/usr/bin/env python3
"""Paired 3-plane x 3-slice comparisons of two participant-count maps."""
from pathlib import Path
from itertools import combinations
import hashlib,json
import nibabel as nib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm

root=Path(__file__).resolve().parents[1];out=root/'results/map_comparison';out.mkdir(parents=True,exist_ok=True)
refpath=root/'ROIs/CZ_ROIs_MNI_1mm/MNI152_T1_1mm_brain.nii';ref=nib.as_closest_canonical(nib.load(refpath));brain=ref.get_fdata()
palettes={
'PERSON':[['#ffb4cd','#ff80ad','#e94e88','#c51b61','#8e0045'],['#ffb38a','#ff855a','#f65335','#d72c1f','#a50f15']],
'FOOD':[['#a5eff5','#6bd9e8','#30bed7','#0799bf','#007396'],['#adc8ff','#7da6f7','#4e7de2','#3154bc','#20317f']],
'LOCATION':[['#c3ee87','#9cda55','#70bd2b','#489b15','#287400'],['#a0e4cd','#68c9aa','#35ac88','#118465','#005940']]}
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'svg.fonttype':'none'})
norm=BoundaryNorm(np.arange(.5,6),5);vmax=np.percentile(brain[brain>0],99)
for context in ['PERSON','FOOD','LOCATION']:
 sources=[root/f'results/searchlight_context_split_half_similarity_results/{context}_split_half_n_significant_MNI.nii.gz',root/f'results/contrast_map/{context}_gt_OTHER_CONTEXTS_n_significant_MNI.nii.gz']
 maps=[]
 for path in sources:
  im=nib.as_closest_canonical(nib.load(path));a=im.get_fdata()
  assert im.shape==ref.shape and np.allclose(im.affine,ref.affine) and np.isin(a,np.arange(6)).all();maps.append(a)
 selected={}
 for axis in [2,1,0]:
  # Equal weighting of the two maps; counts weight more reproducible locations.
  area=np.array([np.sum(a,axis=tuple(i for i in range(3) if i!=axis)) for a in maps])
  scaled=area/np.maximum(area.max(axis=1),1)[:,None]
  candidates=np.where(np.any(area>0,axis=0))[0]
  score=-np.inf;best=None
  for triple in combinations(candidates,3):
   if min(np.diff(triple))<12:continue
   world=ref.affine[axis,axis]*np.array(triple)+ref.affine[axis,3]
   if axis==0 and not (world.min()<0<world.max()):continue
   # High total mapped signal, balanced between analyses, separated by >=12 mm.
   per_map=scaled[:,triple].sum(axis=1)
   value=np.sqrt(per_map).sum()
   if value>score:score=value;best=triple
  assert best is not None;selected[axis]=list(map(int,best))
 fig=plt.figure(figsize=(20,14.5),facecolor='#080808')
 outer=fig.add_gridspec(3,3,left=.055,right=.985,bottom=.15,top=.87,wspace=.14,hspace=.22)
 planes=[(2,'Axial','z'),(1,'Coronal','y'),(0,'Sagittal','x')]
 # Full brain extent with modest surrounding space; no differing crops between paired maps.
 points=np.array(np.where(brain>0));lo=points.min(axis=1);hi=points.max(axis=1)
 for row,(axis,label,coord) in enumerate(planes):
  remaining=[i for i in range(3) if i!=axis]
  for col,k in enumerate(selected[axis]):
   inner=outer[row,col].subgridspec(1,2,wspace=.025)
   for j in range(2):
    ax=fig.add_subplot(inner[0,j]);a=np.take(maps[j],k,axis=axis).T
    ax.imshow(np.take(brain,k,axis=axis).T,origin='lower',cmap='gray',vmin=0,vmax=vmax,interpolation='bilinear')
    ax.imshow(np.ma.masked_equal(a,0),origin='lower',cmap=ListedColormap(palettes[context][j]),norm=norm,interpolation='nearest')
    ax.set_xlim(lo[remaining[0]]-3,hi[remaining[0]]+3);ax.set_ylim(lo[remaining[1]]-3,hi[remaining[1]]+3);ax.set_axis_off()
    world=ref.affine[axis,axis]*k+ref.affine[axis,3]
    ax.set_title(f'{coord} = {world:+.0f} mm\n'+('Searchlight' if j==0 else 'Activation'),fontsize=12,color='white',pad=7)
    ax.text(.02,.5,'P' if axis==0 else 'L',transform=ax.transAxes,color='white',fontsize=9)
    ax.text(.98,.5,'A' if axis==0 else 'R',transform=ax.transAxes,color='white',ha='right',fontsize=9)
    if col==0 and j==0:ax.text(-.17,.5,label,rotation=90,transform=ax.transAxes,color='white',fontsize=17,ha='right',va='center',weight='bold')
 fig.suptitle(f'{context} | Context similarity and activation',color='white',fontsize=28,weight='bold',y=.97)
 fig.text(.5,.918,'Paired views at nine MNI locations · Each voxel shows the number of significant participants (1–5)',color='#dddddd',fontsize=16,ha='center')
 for j in range(2):
  cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=ListedColormap(palettes[context][j])),cax=fig.add_axes([.17+j*.43,.09,.26,.018]),orientation='horizontal',ticks=range(1,6))
  cb.outline.set_visible(False);cb.ax.tick_params(colors='white',length=0,labelsize=12)
  cb.set_label(('Searchlight: joint max-stat FWE p < .05' if j==0 else 'Activation: voxelwise FWE p < .001'),color='white',fontsize=13,labelpad=7)
 fig.text(.5,.02,'Searchlight: context-specific diagonal − off-diagonal similarity. Activation: context > mean of other three contexts (including CONTROL).\nMNI152 · Neurological orientation · Same slice within each pair · Different source tests and thresholds; descriptive participant overlap.',ha='center',color='#cccccc',fontsize=11,linespacing=1.6)
 stem=context+'_searchlight_vs_activation_3planes_3slices'
 for ext in ['png','pdf']:
  target=out/f'{stem}.{ext}'
  if target.exists():raise FileExistsError(target)
  fig.savefig(target,dpi=220,facecolor=fig.get_facecolor())
 plt.close(fig)
 (out/f'{stem}_manifest.json').write_text(json.dumps(dict(context=context,sources=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sources],underlay=str(refpath),palettes=dict(searchlight=palettes[context][0],activation=palettes[context][1]),slices={label:[float(ref.affine[axis,axis]*k+ref.affine[axis,3]) for k in selected[axis]] for axis,label,_ in planes},selection='Each plane: maximize sum of square-root normalized map slice counts across two maps; >=12-mm separation; sagittal views include both hemispheres',display='Counts 1–5, zero transparent; paired views prevent overlay occlusion',dimensions=[4400,3190]),indent=2))
 print(context,selected,flush=True)
