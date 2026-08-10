#!/usr/bin/env bash
# Generate MNI-space group maps for context and odor-score searchlight analyses.
#
# Spatial transforms and image arithmetic use FSL.  The final Stouffer and
# Benjamini-Hochberg calculations are delegated to group_searchlight_stouffer.m.
#
# Usage (from anywhere):
#   bash scripts/group_searchlight_mni.sh

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
MRI_DIR="${ROOT_DIR}/MRI"
RESULT_DIR="${ROOT_DIR}/results/group_searchlight_mni"
REFERENCE="${ROOT_DIR}/ROIs/CZ_ROIs_MNI_1mm/MNI152_T1_1mm_brain.nii"
MATLAB_BIN="${MATLAB_BIN:-/Applications/MATLAB_R2025b.app/bin/matlab}"
SUBJECTS=(2 3 4 5 6)

for command_name in flirt convert_xfm fslmaths fslval; do
    command -v "${command_name}" >/dev/null 2>&1 || {
        echo "Required FSL command not found on PATH: ${command_name}" >&2
        exit 1
    }
done
test -x "${MATLAB_BIN}" || {
    echo "MATLAB executable not found: ${MATLAB_BIN}" >&2
    exit 1
}
test -f "${REFERENCE}" || {
    echo "MNI reference not found: ${REFERENCE}" >&2
    exit 1
}

# The installed SPM reader cannot ingest .nii.gz, so use uncompressed NIfTI
# for all generated maps that MATLAB subsequently reads.
export FSLOUTPUTTYPE=NIFTI
mkdir -p "${RESULT_DIR}"/{context,odor_score,transforms,logs}

TIMESTAMP=$(date '+%Y%m%dT%H%M%S')
LOG_FILE="${RESULT_DIR}/logs/group_searchlight_mni_${TIMESTAMP}.log"
MANIFEST_FILE="${RESULT_DIR}/logs/group_searchlight_mni_${TIMESTAMP}_manifest.tsv"
PARAMETERS_FILE="${RESULT_DIR}/logs/group_searchlight_parameters.tsv"
exec > >(tee "${LOG_FILE}") 2>&1

printf 'subject\tanalysis\teffect_input\tp_input\tfunc2T1\tstd2T1\tT1_to_MNI\tfunc_to_MNI\n' > "${MANIFEST_FILE}"
{
    printf 'parameter\tvalue\n'
    printf 'subjects\t%s\n' "${SUBJECTS[*]}"
    printf 'mni_reference\t%s\n' "${REFERENCE}"
    printf 'mni_resolution\t1 mm\n'
    printf 'effect_and_p_interpolation\tnearestneighbour\n'
    printf 'support_interpolation\tnearestneighbour\n'
    printf 'group_validity\tfinite warped effect and p maps from all 5 subjects\n'
    printf 'stouffer\tone-sided unweighted; z_i=norminv(1-p_i), Z=sum(z_i)/sqrt(5)\n'
    printf 'fdr\tBenjamini-Hochberg, separate context and odor_score families, q<0.05\n'
} > "${PARAMETERS_FILE}"

echo "Started: $(date -Iseconds)"
echo "Repository root: ${ROOT_DIR}"
echo "MNI reference: ${REFERENCE}"
echo "Subjects: ${SUBJECTS[*]}"

require_matrix() {
    local matrix_file=$1
    test -f "${matrix_file}" || {
        echo "Missing transform matrix: ${matrix_file}" >&2
        exit 1
    }
    test "$(awk 'END { print NR }' "${matrix_file}")" -eq 4 || {
        echo "Transform does not have four rows: ${matrix_file}" >&2
        exit 1
    }
    awk 'NF != 4 { exit 1 }' "${matrix_file}" || {
        echo "Transform is not a 4x4 whitespace-delimited matrix: ${matrix_file}" >&2
        exit 1
    }
}

verify_mni_geometry() {
    local image_file=$1
    local axis
    for axis in 1 2 3; do
        test "$(fslval "${image_file}" "dim${axis}")" = "$(fslval "${REFERENCE}" "dim${axis}")" || {
            echo "MNI geometry check failed for ${image_file} (dim${axis})" >&2
            exit 1
        }
    done
}

make_native_support() {
    local effect_file=$1
    local support_file=$2
    # Decoding outputs are finite only at valid searchlight centers.  The
    # arithmetic preserves finite voxels as 1 and converts NaNs to 0.
    fslmaths "${effect_file}" -mul 0 -add 1 -nan -bin "${support_file}"
}

process_analysis() {
    local analysis_name=$1
    local effect_suffix=$2
    local source_directory
    local subject effect_input p_input coreg_dir func_to_t1 std_to_t1
    local t1_to_mni func_to_mni effect_mni p_mni support_native support_nearest_mni support_mni support_p_mni

    case "${analysis_name}" in
        context) source_directory='searchlight_context_loo' ;;
        odor_score) source_directory='searchlight_score_odor_loo' ;;
        *) echo "Unsupported analysis: ${analysis_name}" >&2; exit 1 ;;
    esac

    echo "Processing ${analysis_name}"
    for subject in "${SUBJECTS[@]}"; do
        case "${analysis_name}" in
            context)
                effect_input="${MRI_DIR}/subj_${subject}/nifti/single_trial_by_category/${source_directory}/context_subj${subject}_${effect_suffix}.nii"
                p_input="${MRI_DIR}/subj_${subject}/nifti/single_trial_by_category/${source_directory}/context_subj${subject}_p_uncorrected.nii"
                ;;
            odor_score)
                effect_input="${MRI_DIR}/subj_${subject}/nifti/single_trial_by_category/${source_directory}/odor_score_subj${subject}_${effect_suffix}.nii"
                p_input="${MRI_DIR}/subj_${subject}/nifti/single_trial_by_category/${source_directory}/odor_score_subj${subject}_p_uncorrected.nii"
                ;;
        esac

        for required_file in "${effect_input}" "${p_input}"; do
            test -f "${required_file}" || {
                echo "Missing input map: ${required_file}" >&2
                exit 1
            }
        done

        coreg_dir="${MRI_DIR}/subj_${subject}/nifti/coreg"
        func_to_t1="${coreg_dir}/func2T1mat"
        std_to_t1="${coreg_dir}/std2T1mat"
        require_matrix "${func_to_t1}"
        require_matrix "${std_to_t1}"

        t1_to_mni="${RESULT_DIR}/transforms/subj${subject}_T1_to_MNI.mat"
        func_to_mni="${RESULT_DIR}/transforms/subj${subject}_func_to_MNI.mat"
        convert_xfm -omat "${t1_to_mni}" -inverse "${std_to_t1}"
        # FSL concatenation applies the second transform first:
        # functional -> T1 -> MNI.
        convert_xfm -omat "${func_to_mni}" -concat "${t1_to_mni}" "${func_to_t1}"

        effect_mni="${RESULT_DIR}/${analysis_name}/subj${subject}_${effect_suffix}_mni.nii"
        p_mni="${RESULT_DIR}/${analysis_name}/subj${subject}_p_uncorrected_mni.nii"
        support_native="${RESULT_DIR}/${analysis_name}/subj${subject}_support_native.nii"
        support_nearest_mni="${RESULT_DIR}/${analysis_name}/subj${subject}_support_nearest_mni.nii"
        support_mni="${RESULT_DIR}/${analysis_name}/subj${subject}_support_mni.nii"
        support_p_mni="${RESULT_DIR}/${analysis_name}/subj${subject}_support_p_mni.nii"

        # Searchlight values are defined at discrete native-space centers;
        # nearest-neighbour resampling preserves empirical values and avoids
        # creating interpolated p-values.
        flirt -in "${effect_input}" -ref "${REFERENCE}" -applyxfm -init "${func_to_mni}" \
            -interp nearestneighbour -out "${effect_mni}"
        flirt -in "${p_input}" -ref "${REFERENCE}" -applyxfm -init "${func_to_mni}" \
            -interp nearestneighbour -out "${p_mni}"
        make_native_support "${effect_input}" "${support_native}"
        flirt -in "${support_native}" -ref "${REFERENCE}" -applyxfm -init "${func_to_mni}" \
            -interp nearestneighbour -out "${support_nearest_mni}"
        # Define inference support from the actual finite warped effect and
        # p values, not merely from the nearest native support image.
        fslmaths "${effect_mni}" -mul 0 -add 1 -nan -bin "${support_mni}"
        fslmaths "${p_mni}" -mul 0 -add 1 -nan -bin "${support_p_mni}"
        fslmaths "${support_mni}" -mul "${support_p_mni}" "${support_mni}"

        verify_mni_geometry "${effect_mni}"
        verify_mni_geometry "${p_mni}"
        verify_mni_geometry "${support_mni}"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${subject}" "${analysis_name}" "${effect_input}" "${p_input}" \
            "${func_to_t1}" "${std_to_t1}" "${t1_to_mni}" "${func_to_mni}" >> "${MANIFEST_FILE}"
    done

    local coverage="${RESULT_DIR}/${analysis_name}/group_coverage.nii"
    local valid_mask="${RESULT_DIR}/${analysis_name}/group_valid_all_subjects.nii"
    local mean_effect="${RESULT_DIR}/${analysis_name}/group_mean_${effect_suffix}.nii"
    local subject

    fslmaths "${REFERENCE}" -mul 0 "${coverage}"
    for subject in "${SUBJECTS[@]}"; do
        fslmaths "${coverage}" -add "${RESULT_DIR}/${analysis_name}/subj${subject}_support_mni.nii" "${coverage}"
    done
    fslmaths "${coverage}" -thr "${#SUBJECTS[@]}" -bin "${valid_mask}"

    fslmaths "${RESULT_DIR}/${analysis_name}/subj${SUBJECTS[0]}_${effect_suffix}_mni.nii" \
        -mul 0 "${mean_effect}"
    for subject in "${SUBJECTS[@]}"; do
        fslmaths "${mean_effect}" -add "${RESULT_DIR}/${analysis_name}/subj${subject}_${effect_suffix}_mni.nii" "${mean_effect}"
    done
    fslmaths "${mean_effect}" -div "${#SUBJECTS[@]}" -mas "${valid_mask}" "${mean_effect}"

    if test "${analysis_name}" = 'context'; then
        fslmaths "${mean_effect}" -add 0.25 -mas "${valid_mask}" \
            "${RESULT_DIR}/context/group_mean_raw_accuracy.nii"
    fi

    "${MATLAB_BIN}" -batch "addpath('${SCRIPT_DIR}'); group_searchlight_stouffer('${RESULT_DIR}', '${analysis_name}', [${SUBJECTS[*]}]);"
    echo "Completed ${analysis_name}"
}

process_analysis context accuracy_minus_chance
process_analysis odor_score meanEvidence

echo "Completed: $(date -Iseconds)"
echo "Log: ${LOG_FILE}"
echo "Manifest: ${MANIFEST_FILE}"
