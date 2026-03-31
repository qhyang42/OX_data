#!/bin/bash
# dicom to nifti converter for OX project
# convert all files in $FILEPATH at a time and put them into correct folders

SUBJNAME=$1
DIRNAME=$2

#SUBJDIR=/Volumes/ExtremeSSD/OX_DATA/MRI/"$SUBJNAME"/
#FILEPATH=/Volumes/ExtremeSSD/OX_DATA/MRI/"$SUBJNAME"/"$DIRNAME"

SUBJDIR=/Users/qhyang/Desktop/OX_DATA/MRI/"$SUBJNAME"/
FILEPATH=/Users/qhyang/Desktop/OX_DATA/MRI/"$SUBJNAME"/"$DIRNAME"
# convert files and put
dcm2niix -o "$SUBJDIR"/nifti $FILEPATH

# make folders if not exist in "$SUBJDIR"/nifti. The folders needed are info, localizer, func, anat, wb
mkdir -p "$SUBJDIR"/nifti/info
mkdir -p "$SUBJDIR"/nifti/localizer
mkdir -p "$SUBJDIR"/nifti/func
mkdir -p "$SUBJDIR"/nifti/anat
mkdir -p "$SUBJDIR"/nifti/wb

# move all .json files to info
mv "$SUBJDIR"/nifti/*.json "$SUBJDIR"/nifti/info/

# move *AA_Head*.nii to localizer
mv "$SUBJDIR"/nifti/*AAHead*.nii "$SUBJDIR"/nifti/localizer/

# move *Run*.nii to func
mv "$SUBJDIR"/nifti/*Run*.nii "$SUBJDIR"/nifti/func/

# move *MPRAGE*.nii to anat
mv "$SUBJDIR"/nifti/*MPRAGE*.nii "$SUBJDIR"/nifti/anat/

# move *WB*.nii to wb
mv "$SUBJDIR"/nifti/*WB*.nii "$SUBJDIR"/nifti/wb/




