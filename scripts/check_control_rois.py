#!/usr/bin/env python3
"""Validate saved controls against native GM and finite sniff GLMsingle data; render QC."""
from pathlib import Path
import os,json,csv,sys
os.environ.setdefault('MPLCONFIGDIR','/tmp/ox_control_mpl')
import nibabel as nib
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'ROIs/control_rois'
wm_mode='--white-matter' in sys.argv
manifest=json.loads((out/('white_matter_manifest.json' if wm_mode else 'manifest.json')).read_text());rows=[]
names=['control_WM'] if wm_mode else ['control_GM_putamen','control_A1']
for s in ([] if wm_mode else range(2,7)):
 base=ROOT/f'MRI/subj_{s}/nifti';gmim=nib.load(base/'coreg/gm_mask_thr05_func.nii');gm=gmim.get_fdata()>0
 gmflat=np.flatnonzero(gm.ravel(order='F'));p=base/'sniff_single_trial_by_category_physio/TYPED_FITHRF_GLMDENOISE_RR.mat'
 with h5py.File(p) as f:
  d=f['modelmd'];assert d.shape==(800,1,1,len(gmflat))
  for r in [r for r in manifest['native'] if r['subject']==f'subj_{s}']:
   im=nib.load(r['output']);a=im.get_fdata();assert np.isin(a,[0,1]).all() and a.shape==gm.shape and np.allclose(im.affine,gmim.affine)
   features=np.flatnonzero(a.ravel(order='F')[gmflat]>0)
   values=np.asarray(d[:,:, :,features]).reshape(800,-1);finite=np.isfinite(values).all(0)
   row={'subject':f'subj_{s}','roi':r['roi'],'gm_voxels':len(features),'finite_sniff_features':int(finite.sum())}
   print(row,flush=True)
   row['status_anatomical_gm']='sufficient' if finite.sum()>=10 else 'insufficient'
   for suffix,name in [('p001','spmT_0001_uncorrected_p001.nii'),('FWE_p001','spmT_0001_FWE_p001.nii')]:
    restrict=nib.load(base/'first_level_model_sniff_physio'/name).get_fdata().ravel(order='F')[gmflat][features]>0
    n=int((finite&restrict).sum());row['finite_odor_'+suffix]=n;row['status_odor_'+suffix]='sufficient' if n>=10 else 'insufficient'
   rows.append(row)
if rows:
 with (out/'feature_qc.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
# Each row shows axial/coronal sections through each ROI's largest section.
fig,axes=plt.subplots(6,2*len(names),figsize=(6*len(names),16),facecolor='white')
for row,s in enumerate([None,2,3,4,5,6]):
 refpath=ROOT/'ROIs/CZ_ROIs_MNI_1mm/MNI152_T1_1mm_brain.nii' if s is None else Path(next(r['reference']['path'] for r in manifest['native'] if r['subject']==f'subj_{s}'))
 ref=nib.as_closest_canonical(nib.load(refpath));brain=ref.get_fdata().squeeze()
 for n,name in enumerate(names):
  p=out/f'{name}_bilateral_MNI152_1mm.nii.gz' if s is None else ROOT/f'MRI/subj_{s}/nifti/coreg/roi_decoding/control/{name}_bilateral_func_thr02.nii'
  im=nib.as_closest_canonical(nib.load(p));m=im.get_fdata()>0
  for j,axis in enumerate([2,1]):
   ax=axes[row,n*2+j];k=int(np.argmax(m.sum(axis=tuple(i for i in range(3) if i!=axis))))
   a=np.take(brain,k,axis=axis).T;z=np.take(m,k,axis=axis).T
   ax.imshow(a,origin='lower',cmap='gray',vmin=0,vmax=np.percentile(brain[brain>0],99))
   overlay=np.zeros((*z.shape,4));overlay[z]=[1,.3,.1,.85] if n==0 else [0,.8,1,.85];ax.imshow(overlay,origin='lower')
   ax.set_axis_off();ax.set_title(f"{'MNI' if s is None else 'subj_'+str(s)} · {name.replace('control_','')} · {'axial' if axis==2 else 'coronal'}",fontsize=10)
fig.suptitle('Control ROI anatomical QC · neurological orientation',fontsize=11 if wm_mode else 17)
fig.tight_layout(rect=(0,0,1,.98));fig.savefig(out/('white_matter_QC.png' if wm_mode else 'control_roi_QC.png'),dpi=150);plt.close(fig)
print(json.dumps(rows,indent=2))
