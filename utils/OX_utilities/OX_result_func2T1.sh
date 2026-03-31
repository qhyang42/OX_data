#!/bin/bash
### register any result map from func to T1

subjname=$1
fname=$2
outfname=$3

datapath=/Users/qhyang/Desktop/OX_DATA/MRI/${subjname}/nifti
anatdir=anat
funcdir=func
coregdir=coreg

flirt -in ${fname} -ref ${datapath}/${anatdir}/mean*MPRAGE*.nii -applyxfm -init ${datapath}/${coregdir}func2T1mat -out ${outfname}.nii

