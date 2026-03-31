#!/bin/bash

### generate all necessary coreg matrices
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
logfile="$datapath/$coregdir/coreg_gen_log.txt"

# Check if datapath/coregdir exist. mkdir if not
mkdir -p "$datapath/$coregdir"

# Clear or create log file

  > "$logfile"


  echo "Starting brain extraction for mean T1 image..." >> "$logfile"
  # Do brain extraction for mean T1 image
  t1img=$(ls $datapath/$anatdir/mean2*)
  echo "T1 image: $t1img" >> "$logfile"
  t1brain=${datapath}/${anatdir}/betT1brain
  bet "$t1img" "$t1brain" -f 0.2 -B 
  echo "Brain extraction completed for $t1img." >> "$logfile"

  echo "Registering std brain to T1 brain..." >> "$logfile"
  # Register std brain to T1 brain
  flirt -in "$roidir/MNI152_T1_1mm_brain.nii" -ref "$t1brain" -out "$datapath/$coregdir/std2T1_brain" -omat "$datapath/$coregdir/std2T1mat"
  echo "Transformation matrix std2T1mat created." >> "$logfile"

  #echo "Registering T1 to wb..." >> "$logfile"
  # Register T1 to wb
  #wbimg=$(ls $datapath/$wbdir/mean24*)
  #echo "Whole brain image: $wbimg" >> "$logfile"
  #flirt -in "$t1img" -ref "$wbimg" -out "$datapath/$coregdir/T12wb" -omat "$datapath/$coregdir/T12wbmat"
  #echo "Transformation matrix T12wbmat created." >> "$logfile"


  echo "Registering wb to T1..." >> "$logfile"
  # Register wb to T1
  wbimg=$(ls $datapath/$wbdir/mean2*)
  echo "Whole brain image: $wbimg" >> "$logfile"
  epi_reg --epi="$wbimg" --t1="$t1img" --t1brain="$t1brain" --out="$datapath/$coregdir/wb2T1"
  mv "$datapath/$coregdir/wb2T1.mat" "$datapath/$coregdir/wb2T1mat"
  echo "Transformation matrix wb2T1mat created." >> "$logfile"

  
  echo "Registering wb to mean func..." >> "$logfile"
  # Register wb to mean func
  funcimg=$(ls $datapath/$funcdir/mean2*)
  echo "Functional image: $funcimg" >> "$logfile"
  flirt -in "$wbimg" -ref "$funcimg" -out "$datapath/$coregdir/wb2func" -omat "$datapath/$coregdir/wb2funcmat"
  echo "Transformation matrix wb2funcmat created." >> "$logfile"

  echo "Registering mean func to wb..." >> "$logfile"
  # Register mean func to wb
  flirt -ref "$wbimg" -in "$funcimg" -out "$datapath/$coregdir/func2wb" -omat "$datapath/$coregdir/func2wbmat"
  echo "Transformation matrix func2wbmat created." >> "$logfile"

  # flip wb to T1
  convert_xfm -omat "$datapath/$coregdir/T12wbmat" -inverse "$datapath/$coregdir/wb2T1mat"  

  # combine transformation matrices
  convert_xfm -omat "$datapath/$coregdir/T12funcmat" -concat "$datapath/$coregdir/wb2funcmat" "$datapath/$coregdir/T12wbmat"
  convert_xfm -omat "$datapath/$coregdir/func2T1mat" -concat "$datapath/$coregdir/wb2T1mat" "$datapath/$coregdir/func2wbmat" 
