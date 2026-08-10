#!/usr/bin/env bash
# Coregister and significance-mask revised searchlight outputs in T1 space.
#
# Run this script from the default searchlight output directory:
#   <MRIRoot>/subj_<ID>/nifti/single_trial_by_category/searchlight_<target>_loro
#
# Usage:
#   bash OX_searchlight_coreg_thr_updated.sh <target> <subject_id> [correction] [accuracy_interp]
#
# Arguments:
#   target          context or odor
#   subject_id      numeric subject index, e.g. 5
#   correction      fwe (default), fdr, or both
#   accuracy_interp FLIRT interpolation for the accuracy map (default: trilinear)
#
# Examples:
#   bash OX_searchlight_coreg_thr_updated.sh context 5
#   bash OX_searchlight_coreg_thr_updated.sh odor 5 fdr nearestneighbour
#   bash OX_searchlight_coreg_thr_updated.sh context 5 both trilinear

set -euo pipefail

if [[ $# -lt 2 || $# -gt 4 ]]; then
    echo "Usage: $0 <context|odor> <subject_id> [fwe|fdr|both] [accuracy_interp]" >&2
    exit 2
fi

target=$(printf '%s' "$1" | tr '[:upper:]' '[:lower:]')
subject_id=$2
correction=${3:-fwe}
accuracy_interp=${4:-trilinear}

case "$target" in
    context|odor) ;;
    *)
        echo "Error: target must be 'context' or 'odor'; received '$target'." >&2
        exit 2
        ;;
esac

if [[ ! "$subject_id" =~ ^[0-9]+$ ]]; then
    echo "Error: subject_id must be a positive integer; received '$subject_id'." >&2
    exit 2
fi

case "$correction" in
    fwe|fdr|both) ;;
    *)
        echo "Error: correction must be 'fwe', 'fdr', or 'both'; received '$correction'." >&2
        exit 2
        ;;
esac

prefix="${target}_subj${subject_id}"
accuracy_file="${prefix}_accuracy.nii"

# The revised MATLAB functions save two directories below nifti/.
ref_file="../../anat/betT1brain.nii.gz"
xfm_file="../../coreg/func2T1mat"

for required_file in "$accuracy_file" "$ref_file" "$xfm_file"; do
    if [[ ! -e "$required_file" ]]; then
        echo "Error: required file not found: $required_file" >&2
        echo "Run this script from the default searchlight_${target}_loro output directory." >&2
        exit 1
    fi
done

accuracy_t1="${prefix}_accuracy_T1.nii.gz"

# Accuracy is continuous, so use the requested interpolation.
flirt \
    -in "$accuracy_file" \
    -ref "$ref_file" \
    -applyxfm \
    -init "$xfm_file" \
    -interp "$accuracy_interp" \
    -out "$accuracy_t1"

# Replace NaNs after resampling while retaining a separate significance mask.
fslmaths "$accuracy_t1" -nan "$accuracy_t1"

apply_significance_mask() {
    local correction_name=$1
    local mask_suffix

    case "$correction_name" in
        fwe) mask_suffix="sig_fwe_p05" ;;
        fdr) mask_suffix="sig_fdr_q05" ;;
        *)
            echo "Internal error: unsupported correction '$correction_name'." >&2
            exit 2
            ;;
    esac

    local mask_file="${prefix}_${mask_suffix}.nii"
    local mask_t1="${prefix}_${mask_suffix}_T1.nii.gz"
    local masked_accuracy_t1="${prefix}_accuracy_${mask_suffix}_T1.nii.gz"

    if [[ ! -e "$mask_file" ]]; then
        echo "Error: significance mask not found: $mask_file" >&2
        echo "The MATLAB analysis must be run with NumPermutations > 0." >&2
        exit 1
    fi

    # Significance maps are binary labels, so always use nearest-neighbour interpolation.
    flirt \
        -in "$mask_file" \
        -ref "$ref_file" \
        -applyxfm \
        -init "$xfm_file" \
        -interp nearestneighbour \
        -out "$mask_t1"

    fslmaths "$mask_t1" -nan -bin "$mask_t1"
    fslmaths "$accuracy_t1" -mas "$mask_t1" "$masked_accuracy_t1"

    echo "Wrote: $mask_t1"
    echo "Wrote: $masked_accuracy_t1"
}

case "$correction" in
    fwe)
        apply_significance_mask fwe
        ;;
    fdr)
        apply_significance_mask fdr
        ;;
    both)
        apply_significance_mask fwe
        apply_significance_mask fdr
        ;;
esac

echo "Wrote: $accuracy_t1"
echo "Done."
