#!/usr/bin/env python3
"""Register saved binary global searchlight FWE masks and sum participants."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import nibabel as nib
import numpy as np

root=Path(__file__).resolve().parents[1]
out=root/'results/searchlight_context_split_half_similarity_results'
refpath=root/'ROIs/CZ_ROIs_MNI_1mm/MNI152_T1_1mm_brain.nii'
ref=nib.load(refpath)
contexts=['PERSON','FOOD','LOCATION']
def digest(path):return dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
records=[]
# Verify source geometry, discrete masks, and agreement with saved corrected p maps.
for s in range(2,7):
 base=root/f'MRI/subj_{s}/nifti'; folder=base/'olf_context_similarity_searchlight/global'
 prefix=f'olf_context_similarity_subj{s}'
 support=nib.load(folder/f'{prefix}_valid_center_mask.nii');valid=support.get_fdata()
 assert np.isin(valid,[0,1]).all()
 for name in contexts:
  path=folder/f'{prefix}_{name.lower()}_sig_fwe_p05.nii'
  im=nib.load(path);a=im.get_fdata()
  pimg=nib.load(folder/f'{prefix}_{name.lower()}_p_fwe_maxstat.nii')
  assert im.shape==support.shape==pimg.shape and np.allclose(im.affine,support.affine) and np.allclose(im.affine,pimg.affine)
  assert np.isin(a,[0,1]).all() and np.all(a<=valid)
  assert np.array_equal(a>0,(pimg.get_fdata()<.05)&(valid>0))
  records.append(dict(subject=s,context=name,source=digest(path),p_source=digest(folder/f'{prefix}_{name.lower()}_p_fwe_maxstat.nii'),native_significant_voxels=int(a.sum())))
 for name in ['func2T1mat','std2T1mat']:
  m=np.loadtxt(base/'coreg'/name);assert m.shape==(4,4) and np.isfinite(m).all() and abs(np.linalg.det(m))>1e-8
if out.exists():assert not any(out.iterdir()),'Preserving existing output directory'
out.mkdir(parents=True,exist_ok=True)
commands=[];transforms=[];support_sources=[]
env=dict(os.environ,FSLOUTPUTTYPE='NIFTI_GZ',OMP_NUM_THREADS='1')
def run(cmd):
 subprocess.run(cmd,check=True,env=env,capture_output=True,text=True);commands.append(cmd)
def warp(src,target,mat):
 run(['flirt','-in',str(src),'-ref',str(refpath),'-applyxfm','-init',str(mat),'-interp','nearestneighbour','-out',str(target)])
 im=nib.load(target);a=im.get_fdata()
 assert im.shape==ref.shape and np.allclose(im.affine,ref.affine,atol=1e-4) and np.isin(a,[0,1]).all()
 return a.astype(np.uint8)
def save(a,path):
 header=ref.header.copy();header.set_data_dtype(np.uint8)
 nib.save(nib.Nifti1Image(a,ref.affine,header),path)
counts={c:np.zeros(ref.shape,np.uint8) for c in contexts};coverage=np.zeros(ref.shape,np.uint8)
for s in range(2,7):
 print(f'Projecting subject {s}',flush=True)
 base=root/f'MRI/subj_{s}/nifti'; folder=base/'olf_context_similarity_searchlight/global';prefix=f'olf_context_similarity_subj{s}'
 dest=out/f'subj_{s}';dest.mkdir()
 inv=dest/'T1_to_MNI.mat';mat=dest/'func_to_MNI.mat'
 run(['convert_xfm','-omat',str(inv),'-inverse',str(base/'coreg/std2T1mat')])
 run(['convert_xfm','-omat',str(mat),'-concat',str(inv),str(base/'coreg/func2T1mat')])
 transforms.append(dict(subject=s,inputs=[digest(base/'coreg'/name) for name in ['func2T1mat','std2T1mat']],composed=digest(mat)))
 supportpath=folder/f'{prefix}_valid_center_mask.nii';support_sources.append(digest(supportpath))
 support=warp(supportpath,dest/'valid_centers_MNI.nii.gz',mat);coverage+=support
 for name in contexts:
  a=warp(folder/f'{prefix}_{name.lower()}_sig_fwe_p05.nii',dest/f'{name}_sig_fwe_p05_MNI.nii.gz',mat)
  assert np.all(a<=support);counts[name]+=a
  next(r for r in records if r['subject']==s and r['context']==name)['mni_significant_voxels']=int(a.sum())
save(coverage,out/'n_valid_subjects_MNI.nii.gz')
for name in contexts:
 assert np.all(counts[name]<=coverage)
 save(counts[name],out/f'{name}_split_half_n_significant_MNI.nii.gz')
 # Independent check against sum of saved individual binary masks.
 summed=sum(nib.load(out/f'subj_{s}/{name}_sig_fwe_p05_MNI.nii.gz').get_fdata() for s in range(2,7))
 assert np.array_equal(summed,counts[name])
(out/'manifest.json').write_text(json.dumps(dict(subjects=list(range(2,7)),source_analysis='global 6-mm searchlight, sniff-aligned physiology-regressed betas; semantic trials only; 200 splits; 5000 permutations',threshold='saved joint max-stat FWE p < .05 across centers x three contexts',reference=digest(refpath),script=digest(Path(__file__)),transforms=transforms,inputs=records,support_sources=support_sources,commands=commands,interpolation='nearestneighbour',voxel_counts={c:{str(n):int(np.sum(counts[c]==n)) for n in range(1,6)} for c in contexts}),indent=2))
(out/'README.md').write_text('''# Global context split-half similarity: MNI participant overlap

Sources: subjects 2–6, saved global searchlight binary significance masks. Each participant used 6-mm spheres, at least 10 usable features, 200 session-balanced splits and 5,000 within-run context permutations. The tested effect is context-specific diagonal minus bidirectional off-diagonal Fisher-z similarity. Masks use joint studentized max-stat FWE p < .05 across valid centers × three contexts, rather than the activation-map FWE p < .001 threshold.

Existing inverse(std2T1mat) × func2T1mat transforms project masks to repository MNI152 1-mm space using nearest-neighbour interpolation. No registration or statistical model is re-estimated. Each context's integer image is the sum of five binary masks, 0–5. n_valid_subjects_MNI records searchlight-center coverage. Global support remains GM ∩ positive Odor > Rest uncorrected p < .001 ∩ finite beta support within the acquired functional grid.

Plots: coronal slices every 10 mm, neurological orientation, MNI brain underlay, counts 1–5; PERSON red, FOOD blue, LOCATION green. Combined overlapping-context colors are averaged in RGB; consult individual maps for unambiguous counts. These are descriptive overlaps of participant-level inference, not group-level significance tests. Individual MNI masks, transformations and hashes are retained. Checks verify binary source/output values, corrected-p agreement, MNI geometry, support inclusion, all five subjects, and exact sums of saved masks.
''')
for name in [*contexts,'COMBINED']:
 subprocess.run([sys.executable,str(root/'scripts/plot_odor_overlap_committee.py'),'--plane','coronal','--context',name,'--analysis','searchlight'],check=True,env=env)
