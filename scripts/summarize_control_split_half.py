#!/usr/bin/env python3
"""Audit production control split-half CSVs and add prespecified BH-FDR families."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'MRI/group/context_split_half_similarity_physio_control'
def read(name):
 with (OUT/name).open() as f:return list(csv.DictReader(f))
def correct(filename,pcol,outname):
 rows=read(filename);p=np.array([float(r[pcol]) for r in rows]);assert np.isfinite(p).all()
 order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
 for r,v in zip(rows,q):r['q_BH_FDR']=float(v)
 with (OUT/outname).open('w') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 return rows
rows=read('context_split_half_similarity_summary.csv')
assert len(rows)==10 and {int(r['subject']) for r in rows}==set(range(2,7))
assert {r['roi'] for r in rows}=={'control_GM_putamen','control_A1'}
for r in rows:
 assert r['status']=='ok' and int(r['voxel_count'])>=10 and int(r['n_valid_splits'])==200
 assert np.isfinite(float(r['mean_delta_z']))
overall=correct('context_split_half_similarity_group_permutation_summary.csv','p_group_mean_delta_z','control_overall_group_FDR.csv')
context=correct('context_split_half_similarity_context_group_summary.csv','p_group_context_delta_z','control_context_group_FDR.csv')
assert len(overall)==2 and len(context)==6 and all(int(r['n_permutations'])==5000 for r in overall+context)
paths=[ROOT/'scripts/context_split_half_similarity_control_analysis.m',ROOT/'utils/OX_utilities/OX_roi_context_split_half_similarity.m',ROOT/'utils/OX_utilities/OX_context_split_half_similarity.m',ROOT/'utils/OX_utilities/OX_group_context_split_half_similarity.m',Path(__file__)]
for s in range(2,7):
 base=ROOT/f'MRI/subj_{s}/nifti';paths += [base/'coreg/gm_mask_thr05_func.nii',base/'sniff_single_trial_by_category_physio/TYPED_FITHRF_GLMDENOISE_RR.mat']
 paths += [base/f'coreg/roi_decoding/control/{n}_bilateral_func_thr02.nii' for n in ['control_GM_putamen','control_A1']]
provenance={'subjects':list(range(2,7)),'controls':['control_GM_putamen','control_A1'],'GM_intersection':True,'functional_restriction':False,'splits':200,'permutations':5000,'split_seed':1,'permutation_seed_base':1001,'FDR_families':{'overall':'2 control ROI mean delta-z tests','context':'6 control ROI x semantic-context delta-z tests'},'inference':'Equal-weight aggregate within-subject permutation null; pertains to measured participants, not population random effects. Upper-tail add-one p-values. Maximum-across-split statistic remains supplemental, not the primary test.','files':[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths]}
(OUT/'control_analysis_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
print(json.dumps({'overall':overall,'context':context,'participants':rows},indent=2))
