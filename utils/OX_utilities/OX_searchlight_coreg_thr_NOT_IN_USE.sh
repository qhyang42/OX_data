#!/bin/bash
# coreg and thresholding of searchlight result maps -- use inside searchlight result dir 
basename=$1
accthr=$2
itp=$3

flirt -in "$basename".nii -ref ../anat/betT1brain.nii.gz -applyxfm -init ../coreg/func2T1mat -interp $itp -out "$basename"_T1
fslmaths "$basename"_T1 -nan "$basename"_T1
fslmaths "$basename"_T1 -thr $accthr "$basename"_fdr05_T1




# for each odor under each category 
#flirt -in person_odor_searchlight_accuracy.nii -ref ../anat/betT1brain.nii.gz -applyxfm -init ../coreg/func2T1mat -interp nearestneighbour -out person_odor_searchlight_accuracy_T1.nii

#flirt -in food_odor_searchlight_accuracy.nii -ref ../anat/betT1brain.nii.gz -applyxfm -init ../coreg/func2T1mat -interp nearestneighbour -out food_odor_searchlight_accuracy_T1.nii

#flirt -in loc_odor_searchlight_accuracy.nii -ref ../anat/betT1brain.nii.gz -applyxfm -init ../coreg/func2T1mat -interp nearestneighbour -out loc_odor_searchlight_accuracy_T1.nii

#fslmaths person_odor_searchlight_accuracy_T1.nii.gz -nan person_odor_searchlight_accuracy_T1.nii.gz 
#fslmaths food_odor_searchlight_accuracy_T1.nii.gz -nan food_odor_searchlight_accuracy_T1.nii.gz
#fslmaths loc_odor_searchlight_accuracy_T1.nii.gz -nan loc_odor_searchlight_accuracy_T1.nii.gz
#fslmaths person_odor_searchlight_accuracy_T1.nii.gz -thr 0.1 person_odor_searchlight_accuracy_fdr05_T1.nii.gz
#fslmaths food_odor_searchlight_accuracy_T1.nii.gz -thr 0.1 food_odor_searchlight_accuracy_fdr05_T1.nii.gz
#fslmaths loc_odor_searchlight_accuracy_T1.nii.gz -thr 0.1 loc_odor_searchlight_accuracy_fdr05_T1.nii.gz


# for each odor under each category -- trilinear interp 
#flirt -in person_odor_searchlight_accuracy.nii -ref ../anat/betT1brain.nii.gz -applyxfm -init ../coreg/func2T1mat -out person_odor_searchlight_accuracy_T1.nii

#flirt -in food_odor_searchlight_accuracy.nii -ref ../anat/betT1brain.nii.gz -applyxfm -init ../coreg/func2T1mat -out food_odor_searchlight_accuracy_T1.nii

#flirt -in loc_odor_searchlight_accuracy.nii -ref ../anat/betT1brain.nii.gz -applyxfm -init ../coreg/func2T1mat -out loc_odor_searchlight_accuracy_T1.nii

#fslmaths person_odor_searchlight_accuracy_T1.nii.gz -nan person_odor_searchlight_accuracy_T1.nii.gz 
#fslmaths food_odor_searchlight_accuracy_T1.nii.gz -nan food_odor_searchlight_accuracy_T1.nii.gz
#fslmaths loc_odor_searchlight_accuracy_T1.nii.gz -nan loc_odor_searchlight_accuracy_T1.nii.gz
#fslmaths person_odor_searchlight_accuracy_T1.nii.gz -thr 0.1 person_odor_searchlight_accuracy_fdr05_T1.nii.gz
#fslmaths food_odor_searchlight_accuracy_T1.nii.gz -thr 0.1 food_odor_searchlight_accuracy_fdr05_T1.nii.gz
#fslmaths loc_odor_searchlight_accuracy_T1.nii.gz -thr 0.1 loc_odor_searchlight_accuracy_fdr05_T1.nii.gz


# for each category 
#flirt -in person_searchlight_accuracy.nii -ref ../anat/betT1brain.nii.gz -applyxfm -init ../coreg/func2T1mat -interp nearestneighbour -out person_searchlight_accuracy_T1.nii

#flirt -in food_searchlight_accuracy.nii -ref ../anat/betT1brain.nii.gz -applyxfm -init ../coreg/func2T1mat -interp nearestneighbour -out food_searchlight_accuracy_T1.nii

#flirt -in loc_searchlight_accuracy.nii -ref ../anat/betT1brain.nii.gz -applyxfm -init ../coreg/func2T1mat -interp nearestneighbour -out loc_searchlight_accuracy_T1.nii

#fslmaths person_searchlight_accuracy_T1.nii.gz -nan person_searchlight_accuracy_T1.nii.gz 
#fslmaths food_searchlight_accuracy_T1.nii.gz -nan food_searchlight_accuracy_T1.nii.gz
#fslmaths loc_searchlight_accuracy_T1.nii.gz -nan loc_searchlight_accuracy_T1.nii.gz
#fslmaths person_searchlight_accuracy_T1.nii.gz -thr 0.52 person_searchlight_accuracy_fdr05_T1.nii.gz
#fslmaths food_searchlight_accuracy_T1.nii.gz -thr 0.56 food_searchlight_accuracy_fdr05_T1.nii.gz
#fslmaths loc_searchlight_accuracy_T1.nii.gz -thr 0.5 loc_searchlight_accuracy_fdr05_T1.nii.gz
