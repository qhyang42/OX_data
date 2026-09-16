#!/usr/bin/env python3
"""Build bilateral, anatomically selected, median-volume OX control ROIs.
Run with FSL's Python. Existing outputs are never overwritten.
"""
from pathlib import Path
import os, json, hashlib, subprocess, tempfile, heapq, csv
import numpy as np
import nibabel as nib
from nibabel.processing import resample_from_to
from scipy.ndimage import label
ROOT = Path(__file__).resolve().parents[1]
FSL = Path(os.environ.get('FSLDIR', Path.home()/'fsl'))
OUT = ROOT/'ROIs/control_rois'
SUBJECTS = range(2,7)

def provenance(p):
    return {'path':str(p), 'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def run(*args):
    env = dict(os.environ, FSLDIR=str(FSL), FSLOUTPUTTYPE='NIFTI')
    subprocess.run([str(FSL/'bin'/args[0]), *map(str,args[1:])],check=True,env=env,stdout=subprocess.DEVNULL)

def save(mask, ref, path):
    if path.exists(): raise FileExistsError(path)
    hdr=ref.header.copy(); hdr.set_data_dtype(np.uint8)
    nib.save(nib.Nifti1Image(mask.astype('uint8'),ref.affine,hdr),path)

def compact(mask, n, ref, coverage=None):
    original_mask=mask.copy()
    labs,k=label(mask); counts=np.bincount(labs.ravel()); counts[0]=0
    mask=labs==counts.argmax(); assert mask.sum()>=n
    coords=np.argwhere(mask); world=nib.affines.apply_affine(ref.affine,coords)
    center=world.mean(0)
    if coverage is not None:
        covered=np.argwhere(original_mask & coverage)
        assert len(covered)>0
        center=nib.affines.apply_affine(ref.affine,covered).mean(0)
    seed=tuple(coords[np.argmin(((world-center)**2).sum(1))])
    queue=[(0.,seed)]; seen={seed}; out=np.zeros(mask.shape,bool)
    while out.sum()<n:
        _,p=heapq.heappop(queue); out[p]=True
        for ax in range(3):
            for delta in [-1,1]:
                q=list(p);q[ax]+=delta;q=tuple(q)
                if all(0<=q[i]<mask.shape[i] for i in range(3)) and mask[q] and q not in seen:
                    seen.add(q); d=np.linalg.norm(ref.affine[:3,:3]@(np.array(q)-seed))
                    heapq.heappush(queue,(float(d),q))
    return out, nib.affines.apply_affine(ref.affine,seed).tolist()

def main():
    assert not OUT.exists(), 'Output directory already exists; preserve provenance.'
    refpath=ROOT/'ROIs/CZ_ROIs_MNI_1mm/MNI152_T1_1mm_brain.nii'
    ref=nib.load(refpath)
    sources={name:[ROOT/f'ROIs/CZ_ROIs_MNI_1mm/olfROI/{h}_{stem}.nii.gz' for h in ['L','R']] for name,stem in [('TU','TU'),('AON','AON'),('PirF','pirF'),('PirT','pirT')]}
    sources['olfAMG']=[ROOT/f'ROIs/CZ_ROIs_MNI_1mm/amgROI/{h}_{c}.nii.gz' for h in ['L','R'] for c in ['MeA','ACo','CeA','PAC']]
    sources['olfOFC']=[ROOT/'ROIs/OFC_small_MNI152_1mm.nii.gz']
    def aligned(p):
        im=nib.load(p)
        return resample_from_to(im,ref,order=0).get_fdata() if im.shape!=ref.shape or not np.allclose(im.affine,ref.affine) else im.get_fdata()
    union=np.zeros(ref.shape,bool);sizes={}
    for name,paths in sources.items():
        m=np.zeros(ref.shape,bool)
        for p in paths:
            a=aligned(p);assert np.isin(a,[0,1]).all();m|=a>0
        sizes[name]=int(m.sum());union|=m
    median=float(np.median(list(sizes.values())));target=int(2*np.floor(median/2+.5)); half=target//2
    atlas_paths=[FSL/'data/atlases'/p for p in ['Juelich/Juelich-maxprob-thr50-1mm.nii.gz','HarvardOxford/HarvardOxford-cort-maxprob-thr50-1mm.nii.gz','HarvardOxford/HarvardOxford-sub-maxprob-thr50-1mm.nii.gz']]
    j,h,g=map(aligned,atlas_paths)
    x=ref.affine[0,0]*np.arange(ref.shape[0])[:,None,None]+ref.affine[0,3]
    gm=np.isin(g,[2,13]); masks={}; records={}
    common=np.ones(ref.shape,bool)
    with tempfile.TemporaryDirectory(prefix='ox_control_coverage_') as td:
        w=Path(td)
        for s in SUBJECTS:
            c=ROOT/f'MRI/subj_{s}/nifti/coreg'
            run('convert_xfm','-omat',w/'sw.mat','-concat',c/'T12wbmat',c/'std2T1mat')
            run('convert_xfm','-omat',w/'sf.mat','-concat',c/'wb2funcmat',w/'sw.mat')
            run('convert_xfm','-omat',w/'fs.mat','-inverse',w/'sf.mat')
            run('flirt','-in',c/'gm_mask_thr05_func.nii','-ref',refpath,'-applyxfm','-init',w/'fs.mat','-interp','nearestneighbour','-out',w/'gm.nii')
            common &= nib.load(w/'gm.nii').get_fdata()>0
    # XML probability indices + 1 give maximum-probability labels.
    # TE1.0 is the A1 core; use it before broadening to TE1.1/1.2 or V1.
    for name,candidates in [('control_GM_putamen',[np.isin(g,[6,17])&~union]),('control_A1',[np.isin(j,[41,42,43,44,45,46])&gm&~union])]:
        eligible=candidates[0];m=np.zeros(ref.shape,bool);seeds=[]
        for hemi in [x<0,x>0]:
            patch,seed=compact(eligible&hemi,half,ref,common if name=='control_A1' else None);m|=patch;seeds.append(seed)
        assert m.sum()==target and not (m&union).any()
        masks[name]=m;records[name]={'voxel_count':target,'volume_mm3':float(target*np.prod(ref.header.get_zooms()[:3])),'seed_mni_mm_L_R':seeds,'components':int(label(m)[1])}
    print('Source sizes:',sizes,'median:',median,'target:',target,flush=True)
    native=[]
    with tempfile.TemporaryDirectory(prefix='ox_controls_') as td:
        work=Path(td)
        for name,m in masks.items():save(m,ref,work/f'{name}.nii')
        for s in SUBJECTS:
            base=ROOT/f'MRI/subj_{s}/nifti';coreg=base/'coreg'
            refs=list((base/'func').glob('mean2*.nii'))+list((base/'func').glob('mean2*.nii.gz'));assert len(refs)==1
            func=nib.load(refs[0]);gm_native=nib.load(coreg/'gm_mask_thr05_func.nii');assert np.allclose(func.affine,gm_native.affine)
            dest=coreg/'roi_decoding/control';assert not dest.exists()
            run('convert_xfm','-omat',work/'std2wb.mat','-concat',coreg/'T12wbmat',coreg/'std2T1mat')
            run('convert_xfm','-omat',work/'std2func.mat','-concat',coreg/'wb2funcmat',work/'std2wb.mat')
            for name in masks:
                temp=work/f'subj_{s}_{name}.nii'
                run('flirt','-in',work/f'{name}.nii','-ref',refs[0],'-applyxfm','-init',work/'std2func.mat','-interp','trilinear','-out',temp)
                im=nib.load(temp);a=im.get_fdata();assert np.isfinite(a).all() and im.shape==func.shape[:3] and np.allclose(im.affine,func.affine)
                m=a>=.2;gm_m=m&(gm_native.get_fdata()>0)
                assert gm_m.sum()>=10, f'Insufficient native GM features: {s} {name}'
                row={'subject':f'subj_{s}','roi':name,'anatomical_voxels':int(m.sum()),'gm_voxels':int(gm_m.sum()),'volume_mm3':float(m.sum()*np.prod(im.header.get_zooms()[:3]))}
                for suffix,f in [('odor_p001','spmT_0001_uncorrected_p001.nii'),('odor_FWE_p001','spmT_0001_FWE_p001.nii')]:
                    restriction=nib.load(base/'first_level_model_sniff_physio'/f);assert np.allclose(restriction.affine,func.affine)
                    row[suffix+'_gm_voxels']=int((gm_m&(restriction.get_fdata()>0)).sum())
                row['output']=str(dest/f'{name}_bilateral_func_thr02.nii');row['transforms']=[provenance(coreg/t) for t in ['std2T1mat','T12wbmat','wb2funcmat']];row['reference']=provenance(refs[0])
                native.append((row,m,im)); print({k:v for k,v in row.items() if k not in ['transforms','reference']},flush=True)
        # Publish only once every participant has passed feasibility checks.
        OUT.mkdir()
        for name,m in masks.items():save(m,ref,OUT/f'{name}_bilateral_MNI152_1mm.nii.gz')
        for row,m,im in native:
            dest=Path(row['output']);dest.parent.mkdir(exist_ok=True);save(m,im,dest)
        for s in SUBJECTS:
            rows=[r for r,_,_ in native if r['subject']==f'subj_{s}']
            dest=Path(rows[0]['output']).parent
            with (dest/'roi_manifest.tsv').open('w') as f:
                f.write('subject\tgroup\troi_id\tsource\tcontributing_labels\themisphere\toutput_file\tvoxel_count\n')
                for r in rows:f.write(f"{r['subject']}\tcontrol\t{r['roi']}\tFSL anatomical atlas; see ROIs/control_rois/manifest.json\tmedian-volume compact patch\tbilateral\t{Path(r['output']).name}\t{r['anatomical_voxels']}\n")
        manifest={'source_sizes_1mm_voxels':sizes,'median_mm3':median,'target_voxels':target,'method':'407 voxels per hemisphere; largest 6-connected component; seed nearest component centroid; deterministic 6-connected growth by Euclidean distance from seed. Atlas GM maxprob >50%; exclude six olfactory masks. GM control: HO putamen labels 6/17. SPL initially considered but failed subject 2 native GM coverage; no SPL outputs published. A1: Juelich TE1.0/1.1/1.2 labels 41-46; seed nearest centroid of eligible atlas support intersected with all five inverse-mapped native GM masks. Unrestricted centroid A1 and V1 candidates failed native GM feasibility. No task-based seed selection. No task-based selection or functional restriction applied to saved masks. Native: existing affine chain, trilinear interpolation >=0.2.','standard':records,'native':[r for r,_,_ in native],'inputs':[provenance(p) for paths in sources.values() for p in paths]+[provenance(p) for p in atlas_paths]+[provenance(ROOT/f'MRI/subj_{s}/nifti/coreg/gm_mask_thr05_func.nii') for s in SUBJECTS]+[provenance(refpath),provenance(Path(__file__))],'atlas_docs':['https://fsl.fmrib.ox.ac.uk/fsl/docs/other/datasets.html','https://www.fz-juelich.de/en/inm/inm-7/resources/tools/copy_of_jubrain-anatomy-toolbox']}
        (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        rows=[{k:v for k,v in r.items() if k not in ['transforms','reference']} for r,_,_ in native]
        with (OUT/'native_voxel_counts.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)

def add_white_matter():
    """Append a deep cerebral WM control without rebuilding established masks."""
    from scipy.ndimage import distance_transform_edt
    name='control_WM'
    refpath=ROOT/'ROIs/CZ_ROIs_MNI_1mm/MNI152_T1_1mm_brain.nii'
    ref=nib.load(refpath); target=int(json.loads((OUT/'manifest.json').read_text())['target_voxels'])
    atlaspath=FSL/'data/atlases/HarvardOxford/HarvardOxford-sub-prob-1mm.nii.gz'
    atlas=nib.load(atlaspath)
    assert atlas.shape[:3]==ref.shape and np.allclose(atlas.affine,ref.affine)
    output=OUT/f'{name}_bilateral_MNI152_1mm.nii.gz'
    assert not output.exists()
    wm=np.zeros(ref.shape,bool);seeds=[];probabilities=[]
    for volume,anchor in [(0,[-30,-25,0]),(11,[30,-25,0])]:
        prob=np.asarray(atlas.dataobj[...,volume]); support=prob>=90
        interior=distance_transform_edt(support,sampling=ref.header.get_zooms()[:3])>=3
        coords=np.argwhere(interior);world=nib.affines.apply_affine(ref.affine,coords)
        nearby=np.linalg.norm(world-np.array(anchor),axis=1)<=20
        candidate=np.zeros(ref.shape,bool);candidate[tuple(coords[nearby].T)]=True
        # Anatomical target only; do not use GM or task activation to place WM.
        p=coords[nearby][np.argmin(np.linalg.norm(world[nearby]-anchor,axis=1))]
        anchor_mask=np.zeros(ref.shape,bool);anchor_mask[tuple(p)]=True
        patch,seed=compact(candidate,target//2,ref,anchor_mask)
        wm|=patch;seeds.append(seed);probabilities.append([float(prob[patch].min()),float(prob[patch].mean())])
    assert wm.sum()==target and label(wm)[1]==2
    for p in OUT.glob('*_bilateral_MNI152_1mm.nii.gz'):
        assert not (wm&(nib.load(p).get_fdata()>0)).any(), f'Overlaps existing control: {p}'
    prior=json.loads((OUT/'manifest.json').read_text())
    for p in prior['inputs']:
        if '/ROIs/' in p['path'] and ('/olfROI/' in p['path'] or '/amgROI/' in p['path'] or 'OFC_small' in p['path']):
            im=nib.load(p['path']);assert np.allclose(im.affine,ref.affine)
            assert not (wm&(im.get_fdata()>0)).any()
    native=[]
    with tempfile.TemporaryDirectory(prefix='ox_wm_') as td:
        w=Path(td);save(wm,ref,w/'standard.nii')
        for s in SUBJECTS:
            base=ROOT/f'MRI/subj_{s}/nifti';c=base/'coreg'
            refs=list((base/'func').glob('mean2*.nii'))+list((base/'func').glob('mean2*.nii.gz'));assert len(refs)==1
            func=nib.load(refs[0]);dest=c/f'roi_decoding/control/{name}_bilateral_func_thr02.nii';assert not dest.exists()
            run('convert_xfm','-omat',w/'sw.mat','-concat',c/'T12wbmat',c/'std2T1mat')
            run('convert_xfm','-omat',w/'sf.mat','-concat',c/'wb2funcmat',w/'sw.mat')
            run('flirt','-in',w/'standard.nii','-ref',refs[0],'-applyxfm','-init',w/'sf.mat','-interp','trilinear','-out',w/'native.nii')
            im=nib.load(w/'native.nii');a=im.get_fdata();m=a>=.2
            assert np.isfinite(a).all() and m.sum()>=10 and im.shape==func.shape[:3] and np.allclose(im.affine,func.affine)
            # WM tissue check is descriptive; no GM/functional intersection.
            run('convert_xfm','-omat',w/'tf.mat','-concat',c/'wb2funcmat',c/'T12wbmat')
            tissue=c/'T1_fast_pve_2.nii.gz'
            if not tissue.exists(): tissue=base/'anat/betT1brain_pve_2.nii.gz'
            assert tissue.exists()
            run('flirt','-in',tissue,'-ref',refs[0],'-applyxfm','-init',w/'tf.mat','-interp','trilinear','-out',w/f'wm_pve_{s}.nii')
            pve=nib.load(w/f'wm_pve_{s}.nii').get_fdata()
            row={'subject':f'subj_{s}','roi':name,'anatomical_voxels':int(m.sum()),'volume_mm3':float(m.sum()*np.prod(im.header.get_zooms()[:3])),'mean_native_WM_PVE':float(pve[m].mean()),'native_WM_PVE_ge05_voxels':int((m&(pve>=.5)).sum()),'output':str(dest),'reference':provenance(refs[0]),'wm_tissue_input':provenance(tissue),'transforms':[provenance(c/t) for t in ['std2T1mat','T12wbmat','wb2funcmat']]}
            native.append((row,m,im));print({k:v for k,v in row.items() if k not in ['reference','wm_tissue_input','transforms']},flush=True)
        save(wm,ref,output)
        for row,m,im in native:
            dest=Path(row['output']);save(m,im,dest);row['output_sha256']=provenance(dest)['sha256']
            with (dest.parent/'roi_manifest.tsv').open('a') as f:
                f.write(f"{row['subject']}\tcontrol\t{name}\tHarvardOxford cerebral WM probability >=90%; 3mm interior\tL:volume0;R:volume11; 407 voxels per hemisphere\tbilateral\t{dest.name}\t{row['anatomical_voxels']}\n")
    record={'roi':name,'target_volume_mm3':target,'hemisphere_voxels':target//2,'seed_mni_mm_L_R':seeds,'atlas_WM_probability_min_mean_L_R':probabilities,'method':'HarvardOxford cerebral WM probability >=90%, interior distance >=3 mm; within 20mm of fixed anatomical targets (+/-30,-25,0) mm. Deterministic compact 6-connected 407-voxel patch per hemisphere, seed nearest target. Standard-to-native existing affine chain, trilinear >=0.2. No native GM or functional restriction.','standard_output':provenance(output),'native':[r for r,_,_ in native],'inputs':[provenance(atlaspath),provenance(refpath),provenance(Path(__file__))],'qc_note':'WM PVE overlap is descriptive tissue QA, not an intersection. Existing canonical GLMsingle outputs are GM-only; they cannot supply full WM control betas.'}
    (OUT/'white_matter_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    rows=[{k:v for k,v in r.items() if k not in ['reference','wm_tissue_input','transforms','output_sha256']} for r,_,_ in native]
    with (OUT/'white_matter_native_qc.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--add-white-matter',action='store_true')
    args=parser.parse_args()
    add_white_matter() if args.add_white_matter else main()

