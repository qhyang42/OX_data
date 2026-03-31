#!/bin/bash
### produce gm mask and register to func space. require the results of coreg_gen.sh

subjname=$1
datapath=/Users/qhyang/Desktop/OX_DATA/MRI/${subjname}/nifti
funcdir=func
anatdir=anat
coregdir=coreg
wbdir=wb

t1brain=${datapath}/${anatdir}/betT1brain
fast -t 1 -o "${datapath}/${coregdir}/T1_fast" -n 3 -B "$t1brain"

fslmaths "${datapath}/${coregdir}/T1_fast_seg.nii.gz" -thr 2 -uthr 2 -bin "${datapath}/${coregdir}/T1_gm_mask"

flirt -in "${datapath}/${coregdir}/T1_gm_mask.nii.gz" -ref "$datapath/$coregdir/wb2func.nii.gz" -applyxfm -init "$datapath/$coregdir/T12funcmat" -out "${datapath}/${coregdir}/gm_mask_func"

fslmaths  "${datapath}/${coregdir}/gm_mask_func" -thr 0.5 -bin  "${datapath}/${coregdir}/gm_mask_thr05_func"

gunzip "${datapath}/${coregdir}/gm_mask_thr05_func.nii.gz" 
