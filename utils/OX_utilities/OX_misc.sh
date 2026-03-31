#!/bin/bash

# misc command for OX project

# segmentation
fast -g betT1brain.nii.gz

# register wb to T1 with epi_reg
epi_reg --epi=../wb/mean240814_fMRI_OX_NWU_JN_1_WB_MB3_ax_2iso_1750TR_amyg-Brainstem_20240814140715_5.nii --t1=../anat/mean240814_fMRI_OX_NWU_JN_1_Accelerated_Sagittal_MPRAGE_20240814140715_25.nii --t1brain=../anat/betT1brain.nii.gz --out=epiregtest
# in anat -- merge transformation matrix

convert_xfm -omat T12wbmat -inverse wb2T1mat

convert_xfm -omat T12funcmat -concat wb2funcmat T12wbmat
convert_xfm -omat func2T1mat -concat wb2T1mat func2wbmat

#
flirt -in ../anat/betT1brain_seg_1.nii.gz -applyxfm -init T12funcmat -ref ../func/mean240814_fMRI_OX_NWU_JN_1_Run1_MB3_ax_2iso_760TR_amyg-Brainstem_20240814140715_15.nii -out gm_mask_func

fslmaths gm_mask_func.nii.gz -thr 0.2 -bin gm_mask_thr02_func.nii.gz

## coreg pipeline:
# wb to T1 using epi_reg, inverse wb2T1 for T12wb, wb to func, func to wb with flirt. combine to make func2T1 and T12func

