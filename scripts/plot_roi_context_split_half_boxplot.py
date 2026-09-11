#!/usr/bin/env python3
"""Overall split-half context selectivity: subject boxes and group permutation FDR."""
import csv
import hashlib
import json
from pathlib import Path
from datetime import datetime
import shutil
import numpy as np
from scipy.stats import false_discovery_control
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

root=Path(__file__).resolve().parents[1]
out=root/'results/ROI_context_split_half_similarity_results'
out.mkdir(parents=True,exist_ok=True)
order=[('TU','olf_TU'),('AON','AON'),('PirF','PirF'),('PirT','PirT'),('olfAMG','olfAMG'),('OFC','olfOFC')]
subject_rows=[];group={};sources=[]
for suffix in ['', '_TU']:
 folder=root/'MRI/group'/('context_split_half_similarity_physio'+suffix)
 for filename in ['context_split_half_similarity_summary.csv','context_split_half_similarity_group_permutation_summary.csv']:
  path=folder/filename
  sources.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
  with path.open() as f:rows=list(csv.DictReader(f))
  if filename.endswith('group_permutation_summary.csv'):
   group.update({r['roi']:r for r in rows})
  else:subject_rows.extend(rows)
p=np.array([float(group[key]['p_group_mean_delta_z']) for _,key in order])
q=false_discovery_control(p,method='bh')
assert len(p)==6 and np.isfinite(p).all()
values=[];records=[];subject_table=[]
for i,(label,key) in enumerate(order):
 rows=sorted([r for r in subject_rows if r['roi']==key],key=lambda r:int(r['subject']))
 assert [int(r['subject']) for r in rows]==[2,3,4,5,6]
 valid=[r for r in rows if r['status']=='ok' and np.isfinite(float(r['mean_delta_z']))]
 ids=[int(r['subject']) for r in valid]
 assert ids==([2,3,4,6] if key=='olf_TU' else [2,3,4,5,6])
 vals=np.array([float(r['mean_delta_z']) for r in valid]);values.append(vals)
 observed=float(group[key]['observed_group_mean_delta_z'])
 assert np.isclose(vals.mean(),observed)
 for r in rows:
  subject_table.append(dict(roi=label,source_roi=key,subject=r['subject'],delta_z=r['mean_delta_z'],voxel_count=r['voxel_count'],status=r['status']))
 records.append(dict(roi=label,source_roi=key,n=len(vals),group_mean_delta_z=observed,p_permutation=p[i],q_bh_fdr=q[i],significant_fdr=bool(q[i]<.05),subjects=' '.join(map(str,ids))))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(11,7.3));fig.subplots_adjust(left=.12,right=.97,top=.80,bottom=.27)
box=ax.boxplot(values,positions=np.arange(6),widths=.52,patch_artist=True,showfliers=False,
 boxprops=dict(facecolor='#e5eaf0',edgecolor='#4a5866',linewidth=1.4),
 medianprops=dict(color='#24384b',linewidth=2),whiskerprops=dict(color='#4a5866',linewidth=1.3),capprops=dict(color='#4a5866',linewidth=1.3))
colors=['#0072B2','#D55E00','#009E73','#CC79A7','#80621A']
markers=['o','s','^','D','P'];offsets=np.linspace(-.16,.16,5)
for i,(label,key) in enumerate(order):
 for row in subject_table:
  if row['source_roi']!=key or row['status']!='ok':continue
  sid=int(row['subject']);j=sid-2
  ax.scatter(i+offsets[j],float(row['delta_z']),s=68,marker=markers[j],color=colors[j],edgecolor='white',linewidth=.7,zorder=3)
 if q[i]<.05:
  ax.text(i,max(values[i])+.035,'*',ha='center',va='bottom',fontsize=24,color='#111111')
ax.axhline(0,color='#8a8a8a',linewidth=1,linestyle='--',zorder=0)
ax.set_xticks(range(6),[f'{label}\n(n = {len(values[i])})' for i,(label,_) in enumerate(order)])
ax.set_ylabel('Overall context selectivity (Δz)',labelpad=12)
ax.set_ylim(min(-.12,min(map(min,values))-.06),max(map(max,values))+.13)
ax.set_xlim(-.65,5.65);ax.spines[['top','right']].set_visible(False)
ax.spines[['bottom','left']].set_color('#888888');ax.tick_params(axis='x',length=0,pad=10)
ax.yaxis.grid(True,color='#eeeeee');ax.set_axisbelow(True)
fig.suptitle('Overall semantic-context split-half selectivity',fontsize=21,weight='bold',y=.96)
fig.text(.5,.90,'Mean diagonal − mean off-diagonal Fisher-z similarity',ha='center',fontsize=13,color='#444444')
handles=[Line2D([0],[0],marker=m,color='none',markerfacecolor=c,markeredgecolor='white',markersize=9,label=f'Subject {sid}') for sid,c,m in zip(range(2,7),colors,markers)]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.15),ncol=5,frameon=False,fontsize=11)
fig.text(.12,.025,'Boxes: median and interquartile range; whiskers: 1.5 × IQR. Dots: individual participants.\n* BH-FDR q < .05 across six overall ROI tests; 5,000 within-run permutations, equal-weight group means.\nTU excludes subject 5 (<10 usable voxels); OFC = olfOFC. Inference pertains to the measured participants.',fontsize=10,linespacing=1.6)
stem='ROI_context_split_half_overall_delta_z_boxplot'
existing=list(out.glob(stem+'*'))
if existing:
 archive=out/'old'/datetime.now().strftime('%Y%m%dT%H%M%S_%f');archive.mkdir(parents=True)
 for path in existing:shutil.copy2(path,archive/path.name)
for ext in ['png','pdf','svg']:fig.savefig(out/f'{stem}.{ext}',dpi=300,facecolor='white')
plt.close(fig)
for suffix,rows in [('_group_statistics.csv',records),('_subject_values.csv',subject_table)]:
 with (out/f'{stem}{suffix}').open('w') as f:
  w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
(out/f'{stem}_provenance.json').write_text(json.dumps(dict(sources=sources,script=str(Path(__file__)),correction='Benjamini-Hochberg across six overall mean delta-z ROI tests',stars='one star for q < .05',subject_values='mean over 200 splits per participant; splits are not independent observations',box_definition='25th/50th/75th percentiles; 1.5 IQR whiskers; all subjects overlaid'),indent=2))
print(json.dumps(records,indent=2))
