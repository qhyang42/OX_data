#!/bin/bash

### master script for registering amgROIs for OX project
### register std to T1 to wb func to functional

### required input : subjname

subjname=$1

#datapath=/Volumes/ExtremeSSD/OX_DATA/MRI/${subjname}/nifti
datapath=/Users/qhyang/Desktop/OX_DATA/MRI/${subjname}/nifti
funcdir=func
anatdir=anat
coregdir=coreg
wbdir=wb
#roidir=/Volumes/ExtremeSSD/OX_DATA/CZ_ROIs_MNI_1mm/amgROI
roidir=/Users/qhyang/Desktop/OX_DATA/CZ_ROIs_MNI_1mm/amgROI
logfile="$datapath/$coregdir/amg_log.txt"

# Check if datapath/coregdir exist. mkdir if not
mkdir -p "$datapath/$coregdir"

# Clear or create log file
> "$logfile"

# Register all standard space AMG ROIs to func space
combined_roi_folder="$datapath/$coregdir/combined_amg_rois"
mkdir -p "$combined_roi_folder"

# Combine left and right ROIs for each region
for roi in ACo BMA BLA CeA LA MeA PAC PCo; do
  echo "Combining left and right ROIs for $roi..." >> "$logfile"
  fslmaths "$roidir/L_${roi}.nii.gz" -add "$roidir/R_${roi}.nii.gz" "$combined_roi_folder/${roi}_combined.nii.gz"
  echo "Combined ROI created for $roi." >> "$logfile"

  echo "Registering ${roi}_combined from std to T1 space..." >> "$logfile"
  # Register combined ROI from std to T1
  flirt -in "$combined_roi_folder/${roi}_combined.nii.gz" -ref "$datapath/$anatdir/betT1brain" -out "$combined_roi_folder/${roi}_T1" -applyxfm -init "$datapath/$coregdir/std2T1mat"
  echo "Registered ${roi}_combined from std to T1 space." >> "$logfile"

  echo "Registering ${roi}_combined from T1 to wb space..." >> "$logfile"
  # Register combined ROI from T1 to wb
  wbimg=$(ls $datapath/$wbdir/mean24*)
  flirt -in "$combined_roi_folder/${roi}_T1.nii.gz" -ref "$wbimg" -out "$combined_roi_folder/${roi}_wb" -applyxfm -init "$datapath/$coregdir/T12wbmat"
  echo "Registered ${roi}_combined from T1 to wb space." >> "$logfile"

  echo "Registering ${roi}_combined from wb to func space..." >> "$logfile"
  # Register combined ROI from wb to func and store in coregdir
  funcimg=$(ls $datapath/$funcdir/mean24*)
  flirt -in "$combined_roi_folder/${roi}_wb.nii.gz" -ref "$funcimg" -out "$datapath/$coregdir/${roi}_func" -applyxfm -init "$datapath/$coregdir/wb2funcmat"
  echo "Registered ${roi}_combined from wb to func space." >> "$logfile"

  # Threshold, binarize, and add suffix
  echo "Thresholding and binarizing ${roi}_func..." >> "$logfile"
  fslmaths "$datapath/$coregdir/${roi}_func" -thr 0.2 -bin "$datapath/$coregdir/${roi}_func_thr02"
  echo "Thresholded and binarized ${roi}_func to ${roi}_func_thr02." >> "$logfile"
done

echo "Coregistration complete. Results are in the $datapath/$coregdir directory." >> "$logfile"
