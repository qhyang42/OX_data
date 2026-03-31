#!/bin/bash

### coregistration for OX project. olf cortex ROIs and hand drawn T1 space ROIs. 
### register std to T1 to wb func to functional

### required input : subjname

subjname=$1

#datapath=/Volumes/ExtremeSSD/OX_DATA/MRI/${subjname}/nifti
datapath=/Users/qhyang/Desktop/OX_DATA/MRI/${subjname}/nifti
funcdir=func
anatdir=anat
coregdir=coreg
wbdir=wb
#roidir=/Volumes/ExtremeSSD/OX_DATA/CZ_ROIs_MNI_1mm
roidir=/Users/qhyang/Desktop/OX_DATA/CZ_ROIs_MNI_1mm
logfile="$datapath/$coregdir/coreg_olf_log.txt"
t1brain=${datapath}/${anatdir}/betT1brain
wbimg=$(ls $datapath/$wbdir/mean24*)
funcimg=$(ls $datapath/$funcdir/mean24*)


# Check if datapath/coregdir exist. mkdir if not
mkdir -p "$datapath/$coregdir"

# Clear or create log file if not a rerun
  > "$logfile"

 

 
  # Register all standard space ROI to func space
  roi_folder="$roidir/olfROI"
  combined_roi_folder="$datapath/$coregdir/combined_rois"
  mkdir -p "$combined_roi_folder"

  # Combine left and right ROIs for each region
  for roi in AON pirF pirT TU; do
    echo "Combining left and right ROIs for $roi..." >> "$logfile"
    fslmaths "$roi_folder/L_${roi}.nii.gz" -add "$roi_folder/R_${roi}.nii.gz" "$combined_roi_folder/${roi}_combined.nii.gz"
    echo "Combined ROI created for $roi." >> "$logfile"

    echo "Registering ${roi}_combined from std to T1 space..." >> "$logfile"
    # Register combined ROI from std to T1
    flirt -in "$combined_roi_folder/${roi}_combined.nii.gz" -ref "$t1brain" -out "$combined_roi_folder/${roi}_T1" -applyxfm -init "$datapath/$coregdir/std2T1mat"
    echo "Registered ${roi}_combined from std to T1 space." >> "$logfile"

    echo "Registering ${roi}_combined from T1 to wb space..." >> "$logfile"
    # Register combined ROI from T1 to wb
    flirt -in "$combined_roi_folder/${roi}_T1.nii.gz" -ref "$wbimg" -out "$combined_roi_folder/${roi}_wb" -applyxfm -init "$datapath/$coregdir/T12wbmat"
    echo "Registered ${roi}_combined from T1 to wb space." >> "$logfile"

    echo "Registering ${roi}_combined from wb to func space..." >> "$logfile"
    # Register combined ROI from wb to func and store in coregdir
    flirt -in "$combined_roi_folder/${roi}_wb.nii.gz" -ref "$funcimg" -out "$datapath/$coregdir/${roi}_func" -applyxfm -init "$datapath/$coregdir/wb2funcmat"
    echo "Registered ${roi}_combined from wb to func space." >> "$logfile"

    fslmaths "$datapath/$coregdir/${roi}_func" -thr 0.2 -bin "$datapath/$coregdir/${roi}_func_thr02"
  done


# Register T1 space ROIs to func space
for roi_mask in meanT1_pir_mask meanT1_mea_mask; do
  if [ -f "$datapath/$anatdir/${roi_mask}.nii.gz" ]; then
    echo "Registering $roi_mask from T1 to wb space..." >> "$logfile"
    # Register ROI from T1 to wb
    flirt -in "$datapath/$anatdir/${roi_mask}.nii.gz" -ref "$wbimg" -out "$datapath/$coregdir/${roi_mask}_wb" -applyxfm -init "$datapath/$coregdir/T12wbmat"
    echo "Registered ${roi_mask} from T1 to wb space." >> "$logfile"

    echo "Registering $roi_mask from wb to func space..." >> "$logfile"
    # Register ROI from wb to func
    flirt -in "$datapath/$coregdir/${roi_mask}_wb.nii.gz" -ref "$funcimg" -out "$datapath/$coregdir/${roi_mask}_func" -applyxfm -init "$datapath/$coregdir/wb2funcmat"
    echo "Registered ${roi_mask} from wb to func space." >> "$logfile"
  else
    echo "$roi_mask not found in $datapath/$anatdir, skipping..." >> "$logfile"
  fi
done

echo "Coregistration complete. Results are in the $datapath/$coregdir directory." >> "$logfile"

