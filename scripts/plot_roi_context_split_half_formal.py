#!/usr/bin/env python3
"""Plot saved group Fisher-mean correlations without recomputing analysis."""
import argparse
import csv
import hashlib
import json
import shutil
from datetime import datetime
from scipy.stats import false_discovery_control
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

parser=argparse.ArgumentParser()
parser.add_argument('--annotation',choices=['none','diagonal'],default='diagonal')
args=parser.parse_args()
root=Path(__file__).resolve().parents[1]
out=root/'results/ROI_context_split_half_similarity_results'; out.mkdir(parents=True,exist_ok=True)
order=[('TU','olf_TU'),('AON','AON'),('PirF','PirF'),('PirT','PirT'),('olfAMG','olfAMG'),('OFC','olfOFC')]
contexts=['PERSON','FOOD','LOCATION']; results={}; sources=[]
for suffix in ['', '_TU']:
 folder=root/'MRI/group'/('context_split_half_similarity_physio'+suffix)
 source=folder/'context_split_half_similarity_group_results.mat'
 sources.append(dict(path=str(source),sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
 with h5py.File(source) as f:
  g=f['group_results']
  names=[''.join(chr(x) for x in f[r][:].flat) for r in g['roi_names'][:].flat]
  labels=[''.join(chr(x) for x in f[r][:].flat) for r in g['context_order'][:].flat]
  assert labels==contexts
  matrix=g['group_fisher_mean_r_matrix'][:].T
  subjects=g['subject_mean_fisher_z_matrix'][:].T
  ids=g['subject_ids'][:].ravel().astype(int)
  for j,name in enumerate(names):
   valid=np.isfinite(subjects[:,j]).all(axis=(1,2))
   assert np.allclose(np.tanh(np.mean(subjects[valid,j],axis=0)),matrix[j])
   results[name]=dict(matrix=matrix[j],subject_z=subjects[:,j],ids=ids[valid].tolist(),p=np.ones(3),delta_z=np.full(3,np.nan))
 with open(folder/'context_split_half_similarity_context_group_summary.csv') as f:
  for row in csv.DictReader(f):
   item=results[row['roi']]; c=contexts.index(row['context'])
   item['p'][c]=float(row['p_group_context_delta_z'])
   item['delta_z'][c]=float(row['observed_group_context_delta_z'])
# Prespecified display family: 6 ROIs x 3 context-selectivity tests.
p=np.concatenate([results[key]['p'] for _,key in order]); ranked=np.argsort(p)
q=np.empty_like(p);q[ranked]=np.minimum(1,np.minimum.accumulate((p[ranked]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
assert len(p)==18 and np.isfinite(p).all() and np.all((p>=0)&(p<=1))
assert np.allclose(q, false_discovery_control(p, method='bh'))
for i,(_,key) in enumerate(order):results[key]['q']=q[i*3:i*3+3]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'svg.fonttype':'none'})
norm=TwoSlopeNorm(vmin=-.30,vcenter=0,vmax=.30)
fig,axes=plt.subplots(2,3,figsize=(12,8.6))
fig.subplots_adjust(left=.10,right=.89,bottom=.22,top=.87,wspace=.48,hspace=.55)
rows=[]
for ax,(label,key) in zip(axes.flat,order):
 item=results[key];m=item['matrix']
 im=ax.imshow(m,cmap='RdBu_r',norm=norm)
 ax.set_xticks(range(3),contexts,rotation=25,ha='right',fontsize=10)
 ax.set_yticks(range(3),contexts,fontsize=10)
 ax.set_xlabel('Half B',fontsize=10);ax.set_ylabel('Half A',fontsize=10)
 ax.set_title(f'{label}  |  n = {len(item["ids"])}',weight='bold',pad=10)
 ax.set_xticks(np.arange(-.5,3,1),minor=True);ax.set_yticks(np.arange(-.5,3,1),minor=True)
 ax.grid(which='minor',color='white',linewidth=2);ax.tick_params(which='both',length=0)
 for spine in ax.spines.values():spine.set_visible(False)
 for a in range(3):
  for b in range(3):
   star=''
   if args.annotation=='diagonal' and a==b:
    value=item['q'][a];star='***' if value<.001 else '**' if value<.01 else '*' if value<.05 else 'ns'
   if a==b:
    ax.text(b,a,f'{item["delta_z"][a]:.2f}'+('\n'+star if star else ''),ha='center',va='center',fontsize=11,color='white' if abs(m[a,b])>.20 else '#202020')
   rows.append(dict(roi=label,source_roi=key,context_a=contexts[a],context_b=contexts[b],fisher_mean_r=m[a,b],n=len(item['ids']),subjects=' '.join(map(str,item['ids'])),context_delta_z=item['delta_z'][a] if a==b else '',context_selectivity_p=item['p'][a] if a==b else '',context_selectivity_q=item['q'][a] if a==b else '',annotation=star))
fig.suptitle('Semantic-context split-half similarity',fontsize=20,weight='bold',y=.96)
fig.text(.5,.912,'Equal-weight participant means in Fisher-z space, back-transformed to r',ha='center',fontsize=11,color='#444444')
cb=fig.colorbar(im,cax=fig.add_axes([.925,.29,.018,.45]));cb.set_label('Cross-half Pearson correlation (r)');cb.outline.set_visible(False)
foot='PERSON, FOOD and LOCATION only; 200 session-balanced splits. OFC = olfOFC.\nTU: subjects 2, 3, 4, 6; subject 5 has <10 usable voxels. Other ROIs: subjects 2–6.'
if args.annotation=='diagonal':foot+='\nDiagonal text: Δz = diagonal minus bidirectional off-diagonal mean in Fisher-z space; cell colors show r.\n5,000 within-run permutations; BH-FDR across 18 ROI × context tests. * q < .05, ** q < .01, *** q < .001; ns = not significant.'
else:foot+='\nDraft: significance annotations pending selection of the inferential test.'
fig.text(.10,.035,foot,fontsize=9,linespacing=1.5)
stem='ROI_context_split_half_group_heatmaps'+('_draft' if args.annotation=='none' else '')
existing=[out/f'{stem}{suffix}' for suffix in ['.png','.pdf','.svg','_data.csv','_provenance.json'] if (out/f'{stem}{suffix}').exists()]
if existing:
 archive=out/'old'/datetime.now().strftime('%Y%m%dT%H%M%S_%f');archive.mkdir(parents=True)
 for path in existing:shutil.copy2(path,archive/path.name)
for ext in ['png','pdf','svg']:fig.savefig(out/f'{stem}.{ext}',dpi=300,facecolor='white')
plt.close(fig)
with open(out/f'{stem}_data.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
(out/f'{stem}_provenance.json').write_text(json.dumps(dict(sources=sources,roi_order=[x[0] for x in order],annotation=args.annotation,diagonal_text='Saved group context delta z; off-diagonal text omitted; colors retain Fisher-mean r',multiple_comparison_family='Benjamini-Hochberg FDR across 18 ROI x context-selectivity tests; stars use adjusted p-values (q)',interpretation='Inference about measured participants, not population random effects'),indent=2))
print(out/f'{stem}.png')
