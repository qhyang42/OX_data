#!/usr/bin/env python3
"""Descriptive MNI averages and participant overlap of existing FWE contrasts.
Run with $FSLDIR/bin/python scripts/group_fwe_contrasts_mni.py.
Preserves existing inputs and refuses to overwrite the destination.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import h5py
import nibabel as nib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

ROOT = Path(__file__).resolve().parents[1]
CONTRASTS = [('ODOR_gt_REST', 'first_level_model_sniff_physio', 1, 'Odor > Rest'),
             *[(c+'_gt_OTHER_CONTEXTS', 'first_level_model_category_sniff_physio', i, c.title()+' > Other Contexts')
               for i,c in enumerate(('PERSON','FOOD','LOCATION'),1)]]

def digest(p):
    return dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/contrast_map')
    args=parser.parse_args(); out=args.output
    refpath=ROOT/'ROIs/CZ_ROIs_MNI_1mm/MNI152_T1_1mm_brain.nii'
    ref=nib.load(refpath)
    records=[]
    # Preflight all subjects before writing any output.
    for s in range(2,7):
        base=ROOT/f'MRI/subj_{s}/nifti'
        for t in ('func2T1mat','std2T1mat'):
            m=np.loadtxt(base/'coreg'/t)
            assert m.shape==(4,4) and np.isfinite(m).all() and abs(np.linalg.det(m))>1e-8
        for key,folder,idx,name in CONTRASTS:
            src=base/folder/f'spmT_{idx:04d}_FWE_p001.nii'
            a=nib.load(src).get_fdata()
            assert not np.isinf(a).any() and not np.any(a[np.isfinite(a)]<0)
            with h5py.File(base/folder/'SPM.mat') as f:
                ds=f['SPM/xCon/name']
                names=([ ''.join(chr(x) for x in f[r][:].flat) for r in ds[:].flat]
                       if h5py.check_dtype(ref=ds.dtype) else [''.join(chr(x) for x in ds[:].flat)])
                assert names[idx-1] in (name, name+' - All Sessions'), names
                assert f['SPM/Sess/row'].size==80
            records.append(dict(subject=s,contrast=key,source=digest(src),native_significant_voxels=int(np.sum(a>0)),native_nan_background=int(np.isnan(a).sum())))
    if out.exists():
        assert all(p.name in ('old', '.DS_Store') for p in out.iterdir()), 'Preserving existing outputs'
    out.mkdir(parents=True,exist_ok=True)
    env=dict(os.environ,FSLOUTPUTTYPE='NIFTI_GZ',OMP_NUM_THREADS='1')
    commands=[]
    def run(cmd):
        subprocess.run(cmd,check=True,env=env,capture_output=True,text=True); commands.append(cmd)
    def save(a,p,dtype=np.float32):
        h=ref.header.copy(); h.set_data_dtype(dtype)
        nib.save(nib.Nifti1Image(a.astype(dtype),ref.affine,h),p)
    def warp(src,dst,mat):
        run(['flirt','-in',str(src),'-ref',str(refpath),'-applyxfm','-init',str(mat),'-interp','nearestneighbour','-out',str(dst)])
        im=nib.load(dst); assert im.shape==ref.shape and np.allclose(im.affine,ref.affine,atol=1e-4)
        a=im.get_fdata(dtype=np.float32); assert np.isfinite(a).all(); return a
    transforms=[]
    sums={k:np.zeros(ref.shape,np.float32) for k,*_ in CONTRASTS}
    counts={k:np.zeros(ref.shape,np.uint8) for k,*_ in CONTRASTS}
    coverage={k:np.zeros(ref.shape,np.uint8) for k,*_ in CONTRASTS}
    for s in range(2,7):
        print(f'Projecting subj_{s}',flush=True)
        dest=out/f'subj_{s}'; dest.mkdir(); base=ROOT/f'MRI/subj_{s}/nifti'
        inv=dest/'T1_to_MNI.mat'; mat=dest/'func_to_MNI.mat'
        run(['convert_xfm','-omat',str(inv),'-inverse',str(base/'coreg/std2T1mat')])
        run(['convert_xfm','-omat',str(mat),'-concat',str(inv),str(base/'coreg/func2T1mat')])
        transforms.append(dict(subject=s,inputs=[digest(base/'coreg'/t) for t in ('func2T1mat','std2T1mat')],composed=digest(mat)))
        for key,folder,idx,name in CONTRASTS:
            src=base/folder/f'spmT_{idx:04d}_FWE_p001.nii'
            native=nib.load(src); a=native.get_fdata(dtype=np.float32); clean=np.nan_to_num(a,nan=0)
            temp=dest/'clean_native.nii.gz'
            nib.save(nib.Nifti1Image(clean,native.affine,native.header),temp)
            v=warp(temp,dest/f'{key}_FWE_p001_t_MNI.nii.gz',mat); temp.unlink()
            assert np.isin(np.unique(v),np.append(np.unique(clean),0)).all()
            sig=v>0; sums[key]+=v; counts[key]+=sig
            save(sig,dest/f'{key}_significant_MNI.nii.gz',np.uint8)
            mask=warp(base/folder/'mask.nii',dest/f'{key}_coverage_MNI.nii.gz',mat)>0
            assert not np.any(sig & ~mask)
            coverage[key]+=mask
            next(r for r in records if r['subject']==s and r['contrast']==key)['mni_significant_voxels']=int(sig.sum())
    colors=['#377eb8','#4daf4a','#ffff33','#ff7f00','#e41a1c']
    cmap=ListedColormap(colors); norm=BoundaryNorm(np.arange(.5,6),5)
    bg=nib.as_closest_canonical(ref); anatomy=bg.get_fdata()
    zs=list(range(-30,71,10))
    fig,axes=plt.subplots(4,len(zs),figsize=(22,9),facecolor='black')
    summary={}
    for row,(key,*_) in enumerate(CONTRASTS):
        save(sums[key]/5,out/f'{key}_mean_FWE_p001_t_MNI.nii.gz')
        save(counts[key],out/f'{key}_n_significant_MNI.nii.gz',np.uint8)
        save(counts[key]/5,out/f'{key}_fraction_significant_MNI.nii.gz')
        save(coverage[key],out/f'{key}_n_coverage_MNI.nii.gz',np.uint8)
        assert np.all(counts[key]<=coverage[key])
        summary[key]={str(n):int(np.sum(counts[key]==n)) for n in range(1,6)}
        data=nib.as_closest_canonical(nib.load(out/f'{key}_n_significant_MNI.nii.gz')).get_fdata()
        for col,z in enumerate(zs):
            ax=axes[row,col]; iz=int(round((z-bg.affine[2,3])/bg.affine[2,2]))
            ax.imshow(anatomy[:,:,iz].T,origin='lower',cmap='gray',vmin=0,vmax=np.percentile(anatomy[anatomy>0],99))
            ax.imshow(np.ma.masked_less(data[:,:,iz].T,1),origin='lower',cmap=cmap,norm=norm,interpolation='nearest')
            ax.set_axis_off(); ax.set_title(f'z = {z}',color='white',fontsize=9)
            if col==0: ax.text(-.08,.5,key.replace('_gt_',' >\n').replace('_',' '),transform=ax.transAxes,ha='right',va='center',color='white',fontsize=10)
        # Save every axial level containing signal as an additional complete montage.
        active=np.where(np.any(data>0,axis=(0,1)))[0][::5]
        f,axs=plt.subplots(max(1,int(np.ceil(len(active)/10))),10,figsize=(20,max(2,2*np.ceil(len(active)/10))),facecolor='black')
        for ax in axs.flat: ax.set_axis_off()
        for ax,iz in zip(axs.flat,active):
            ax.imshow(anatomy[:,:,iz].T,origin='lower',cmap='gray')
            ax.imshow(np.ma.masked_less(data[:,:,iz].T,1),origin='lower',cmap=cmap,norm=norm,interpolation='nearest')
            ax.set_title(f'z={bg.affine[2,2]*iz+bg.affine[2,3]:.0f}',color='white',fontsize=8)
        f.suptitle(key+' | significant participants: 1 blue, 2 green, 3 yellow, 4 orange, 5 red',color='white')
        f.savefig(out/f'{key}_overlap_axial.png',dpi=130,facecolor='black',bbox_inches='tight'); plt.close(f)
    fig.suptitle('Participant overlap | positive contrasts, subject-level FWE p < .001 | n = 5',color='white',fontsize=16)
    fig.subplots_adjust(left=.10,right=.98,bottom=.12,top=.9,wspace=.02,hspace=.18)
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=fig.add_axes([.35,.05,.35,.02]),orientation='horizontal',ticks=range(1,6))
    cb.ax.tick_params(colors='white'); cb.set_label('Number of significant participants',color='white')
    fig.savefig(out/'group_overlap.png',dpi=170,facecolor='black'); plt.close(fig)
    manifest=dict(subjects=list(range(2,7)),reference=digest(refpath),script=digest(Path(__file__)),interpolation='nearestneighbour',transforms=transforms,inputs=records,commands=commands,voxel_counts_by_n=summary)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    (out/'README.md').write_text('''# Group FWE contrast displays in MNI space

Subjects 2–6; sniff-aligned, physiology-regressed SPM first-level models. Existing positive voxelwise FWE p < .001 maps, no extent threshold. PERSON, FOOD and LOCATION each contrast against the equally weighted other three contexts, including CONTROL. FWE correction is inherited per subject and contrast; no new correction across these four displays is applied.

Existing affine chain: inverse(coreg/std2T1mat) × coreg/func2T1mat. Reference is the repository MNI152 T1 1-mm brain. Nearest-neighbour resampling preserves original thresholded t values and significance decisions. Native NaN background is explicitly set to zero; source nonfinite counts are recorded. No registration is re-estimated, smoothing applied, or subject excluded.

* mean_FWE_p001_t: sum of registered thresholded t values / 5, with nonsignificant and unsupported voxels zero. Descriptive mean of t statistics, not mean contrast beta or a group t test.
* n_significant: integer 0–5; color scale blue=1, green=2, yellow=3, orange=4, red=5; zero transparent.
* fraction_significant: n_significant / 5.
* n_coverage: number of transformed SPM estimation masks containing the voxel; distinguish absence of coverage from nonsignificance.

PNG axial views use neurological convention (image left = anatomical left), MNI z coordinates in mm. Individual registered maps and masks are retained in subj_N directories. These overlap/average displays are descriptive summaries of the five measured participants, not group-level FWE inference. manifest.json records input hashes, transforms, commands, and counts. Validation checks all 20 contrasts, 80 SPM sessions per subject/model, expected saved contrast names, geometry, finite output, source-value preservation and significance within coverage.
''')
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':
    main()
