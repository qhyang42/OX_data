#!/usr/bin/env bash

# Register pooled Brainnetome and in-house ROI masks from MNI space to each
# subject's functional space. Usage:
#   bash scripts/register_decoding_rois_to_func.sh [subject_id ...]
# Examples:
#   bash scripts/register_decoding_rois_to_func.sh        # subjects 2–6
#   bash scripts/register_decoding_rois_to_func.sh 2 5
#   bash scripts/register_decoding_rois_to_func.sh subj_2 subj_5

set -euo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

ATLAS_FILE="${ROOT_DIR}/ROIs/BN_Atlas_246_1mm.nii.gz"
AMYGDALA_DIR="${ROOT_DIR}/ROIs/CZ_ROIs_MNI_1mm/amgROI"
OLFACTORY_DIR="${ROOT_DIR}/ROIs/CZ_ROIs_MNI_1mm/olfROI"
DEFAULT_SUBJECTS=(2 3 4 5 6)

BN_ROI_NAMES=(
    medial_anterior_pfc
    dlpfc
    lateral_frontopolar_pfc
    ofc
    temporal_pole
    parahippocampal_mtc
    insula
    hippocampus
    amygdala
    nacc
)
BN_LEFT_LABELS=(
    "1 11 13"
    "5 15 19 21"
    "27"
    "41 43 45 47 49 51"
    "69 77"
    "109 111 113 115 117 119"
    "163 165 167 169 171 173"
    "215 217"
    "211 213"
    "223"
)
BN_RIGHT_LABELS=(
    "2 12 14"
    "6 16 20 22"
    "28"
    "42 44 46 48 50 52"
    "70 78"
    "110 112 114 116 118 120"
    "164 166 168 170 172 174"
    "216 218"
    "212 214"
    "224"
)
AMYGDALA_ROIS=(ACo BLA BMA CeA LA MeA PAC PCo)
OLFACTORY_ROIS=(AON pirF pirT TU)

for command_name in flirt convert_xfm fslmaths fslstats fslval; do
    command -v "${command_name}" >/dev/null 2>&1 || {
        echo "Required FSL command not found: ${command_name}" >&2
        exit 1
    }
done

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

join_with_commas() {
    local labels="$1"
    echo "${labels// /,}"
}

normalise_subject_id() {
    local input_id="$1"
    if [[ "${input_id}" == subj_* ]]; then
        echo "${input_id}"
    elif [[ "${input_id}" =~ ^[0-9]+$ ]]; then
        echo "subj_${input_id}"
    else
        fail "Subject IDs must be numeric (for example, 2) or subj_<number>: ${input_id}"
    fi
}

assert_nonempty_mask() {
    local image_file="$1"
    fslstats "${image_file}" -V | awk '$1 > 0 { found = 1 } END { exit !found }' ||
        fail "Mask is empty: ${image_file}"
}

assert_functional_geometry() {
    local image_file="$1"
    local reference_file="$2"
    local key
    for key in dim1 dim2 dim3 pixdim1 pixdim2 pixdim3; do
        [[ "$(fslval "${image_file}" "${key}")" == "$(fslval "${reference_file}" "${key}")" ]] ||
            fail "Geometry mismatch (${key}) between ${image_file} and ${reference_file}"
    done
}

assert_binary_mask() {
    local image_file="$1"
    fslstats "${image_file}" -R | awk '$1 != 0 || ($2 != 0 && $2 != 1) { exit 1 }' ||
        fail "Mask is not binary: ${image_file}"
}

mask_voxel_count() {
    local image_file="$1"
    fslstats "${image_file}" -V | awk '{ print $1 }'
}

build_bn_mask() {
    local labels="$1"
    local output_file="$2"
    local first_label="${labels%% *}"
    local label sum_file

    fslmaths "${ATLAS_FILE}" -thr "${first_label}" -uthr "${first_label}" -bin "${output_file}"
    for label in ${labels}; do
        [[ "${label}" == "${first_label}" ]] && continue
        sum_file="${output_file%.nii.gz}_sum.nii.gz"
        fslmaths "${ATLAS_FILE}" -thr "${label}" -uthr "${label}" -bin "${sum_file}"
        fslmaths "${output_file}" -add "${sum_file}" "${output_file}"
    done
    fslmaths "${output_file}" -bin "${output_file}"
    assert_nonempty_mask "${output_file}"
}

write_manifest_row() {
    local roi_id="$1"
    local source="$2"
    local labels="$3"
    local hemisphere="$4"
    local output_file="$5"
    local voxels
    voxels="$(mask_voxel_count "${output_file}")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${CURRENT_SUBJECT}" "${roi_id}" "${source}" "${labels}" "${hemisphere}" \
        "${output_file}" "${voxels}" >> "${CURRENT_MANIFEST}"
    if [[ "${voxels}" == "0" ]]; then
        echo "WARNING: ${CURRENT_SUBJECT} ${roi_id} ${hemisphere} has zero functional voxels." >&2
    fi
}

register_unilateral_mask() {
    local source_file="$1"
    local roi_id="$2"
    local hemisphere="$3"
    local source_name="$4"
    local labels="$5"
    local output_file="${CURRENT_OUTPUT_DIR}/${roi_id}_${hemisphere}_func_thr02.nii.gz"
    local resampled_file="${CURRENT_WORK_DIR}/${roi_id}_${hemisphere}_resampled.nii.gz"

    flirt -in "${source_file}" -ref "${CURRENT_FUNC_REF}" -applyxfm \
        -init "${CURRENT_STD2FUNC_MAT}" -interp trilinear -out "${resampled_file}"
    fslmaths "${resampled_file}" -thr 0.2 -bin "${output_file}"
    assert_functional_geometry "${output_file}" "${CURRENT_FUNC_REF}"
    assert_binary_mask "${output_file}"
    write_manifest_row "${roi_id}" "${source_name}" "${labels}" "${hemisphere}" "${output_file}"
}

make_bilateral_mask() {
    local roi_id="$1"
    local source_name="$2"
    local labels="$3"
    local left_file="${CURRENT_OUTPUT_DIR}/${roi_id}_L_func_thr02.nii.gz"
    local right_file="${CURRENT_OUTPUT_DIR}/${roi_id}_R_func_thr02.nii.gz"
    local output_file="${CURRENT_OUTPUT_DIR}/${roi_id}_bilateral_func_thr02.nii.gz"

    fslmaths "${left_file}" -add "${right_file}" -bin "${output_file}"
    assert_functional_geometry "${output_file}" "${CURRENT_FUNC_REF}"
    assert_binary_mask "${output_file}"
    write_manifest_row "${roi_id}" "${source_name}" "${labels}" "bilateral" "${output_file}"
}

subjects=()
if (( $# == 0 )); then
    subjects=("${DEFAULT_SUBJECTS[@]}")
else
    for subject in "$@"; do
        subjects+=("${subject}")
    done
fi

[[ -f "${ATLAS_FILE}" ]] || fail "Brainnetome atlas not found: ${ATLAS_FILE}"
for roi in "${AMYGDALA_ROIS[@]}"; do
    [[ -f "${AMYGDALA_DIR}/L_${roi}.nii.gz" ]] || fail "Missing amygdala ROI: L_${roi}.nii.gz"
    [[ -f "${AMYGDALA_DIR}/R_${roi}.nii.gz" ]] || fail "Missing amygdala ROI: R_${roi}.nii.gz"
done
for roi in "${OLFACTORY_ROIS[@]}"; do
    [[ -f "${OLFACTORY_DIR}/L_${roi}.nii.gz" ]] || fail "Missing olfactory ROI: L_${roi}.nii.gz"
    [[ -f "${OLFACTORY_DIR}/R_${roi}.nii.gz" ]] || fail "Missing olfactory ROI: R_${roi}.nii.gz"
done

# Complete input validation before creating any subject output directories.
normalised_subjects=()
for subject in "${subjects[@]}"; do
    subject_id="$(normalise_subject_id "${subject}")"
    subject_dir="${ROOT_DIR}/MRI/${subject_id}/nifti"
    coreg_dir="${subject_dir}/coreg"
    [[ -d "${subject_dir}" ]] || fail "Subject directory not found: ${subject_dir}"
    for matrix in std2T1mat T12wbmat wb2funcmat; do
        [[ -f "${coreg_dir}/${matrix}" ]] || fail "Missing ${matrix} for ${subject_id}: ${coreg_dir}/${matrix}"
    done
    functional_references=("${subject_dir}"/func/mean2*.nii "${subject_dir}"/func/mean2*.nii.gz)
    (( ${#functional_references[@]} == 1 )) ||
        fail "Expected exactly one mean functional image for ${subject_id}; found ${#functional_references[@]}"
    normalised_subjects+=("${subject_id}")
done

WORK_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/ox_roi_registration.XXXXXX")"
trap 'rm -rf "${WORK_ROOT}"' EXIT
STANDARD_MASK_DIR="${WORK_ROOT}/standard_masks"
mkdir -p "${STANDARD_MASK_DIR}"

echo "Building pooled Brainnetome masks in standard space..."
for index in "${!BN_ROI_NAMES[@]}"; do
    roi="${BN_ROI_NAMES[${index}]}"
    build_bn_mask "${BN_LEFT_LABELS[${index}]}" "${STANDARD_MASK_DIR}/bn_${roi}_L.nii.gz"
    build_bn_mask "${BN_RIGHT_LABELS[${index}]}" "${STANDARD_MASK_DIR}/bn_${roi}_R.nii.gz"
done

for CURRENT_SUBJECT in "${normalised_subjects[@]}"; do
    subject_dir="${ROOT_DIR}/MRI/${CURRENT_SUBJECT}/nifti"
    coreg_dir="${subject_dir}/coreg"
    functional_references=("${subject_dir}"/func/mean2*.nii "${subject_dir}"/func/mean2*.nii.gz)
    CURRENT_FUNC_REF="${functional_references[0]}"
    CURRENT_OUTPUT_DIR="${coreg_dir}/roi_decoding"
    CURRENT_WORK_DIR="${WORK_ROOT}/${CURRENT_SUBJECT}"
    FINAL_MANIFEST="${CURRENT_OUTPUT_DIR}/roi_manifest.tsv"
    CURRENT_MANIFEST="${CURRENT_WORK_DIR}/roi_manifest.tsv"
    CURRENT_STD2FUNC_MAT="${CURRENT_WORK_DIR}/std2func.mat"
    mkdir -p "${CURRENT_OUTPUT_DIR}" "${CURRENT_WORK_DIR}"

    std2wb_mat="${CURRENT_WORK_DIR}/std2wb.mat"
    convert_xfm -omat "${std2wb_mat}" -concat "${coreg_dir}/T12wbmat" "${coreg_dir}/std2T1mat"
    convert_xfm -omat "${CURRENT_STD2FUNC_MAT}" -concat "${coreg_dir}/wb2funcmat" "${std2wb_mat}"

    printf 'subject\troi_id\tsource\tcontributing_labels\themisphere\toutput_file\tvoxel_count\n' > "${CURRENT_MANIFEST}"
    echo "Registering decoding ROIs for ${CURRENT_SUBJECT}..."

    for index in "${!BN_ROI_NAMES[@]}"; do
        roi="${BN_ROI_NAMES[${index}]}"
        left_labels="$(join_with_commas "${BN_LEFT_LABELS[${index}]}")"
        right_labels="$(join_with_commas "${BN_RIGHT_LABELS[${index}]}")"
        register_unilateral_mask "${STANDARD_MASK_DIR}/bn_${roi}_L.nii.gz" "bn_${roi}" L \
            Brainnetome "${left_labels}"
        register_unilateral_mask "${STANDARD_MASK_DIR}/bn_${roi}_R.nii.gz" "bn_${roi}" R \
            Brainnetome "${right_labels}"
        make_bilateral_mask "bn_${roi}" Brainnetome "L:${left_labels};R:${right_labels}"
    done

    for roi in "${AMYGDALA_ROIS[@]}"; do
        register_unilateral_mask "${AMYGDALA_DIR}/L_${roi}.nii.gz" "amygsub_${roi}" L \
            in_house_amygdala_subregion "${roi}"
        register_unilateral_mask "${AMYGDALA_DIR}/R_${roi}.nii.gz" "amygsub_${roi}" R \
            in_house_amygdala_subregion "${roi}"
        make_bilateral_mask "amygsub_${roi}" in_house_amygdala_subregion "${roi}"
    done

    for roi in "${OLFACTORY_ROIS[@]}"; do
        register_unilateral_mask "${OLFACTORY_DIR}/L_${roi}.nii.gz" "olf_${roi}" L \
            in_house_olfactory_roi "${roi}"
        register_unilateral_mask "${OLFACTORY_DIR}/R_${roi}.nii.gz" "olf_${roi}" R \
            in_house_olfactory_roi "${roi}"
        make_bilateral_mask "olf_${roi}" in_house_olfactory_roi "${roi}"
    done

    output_masks=("${CURRENT_OUTPUT_DIR}"/*_func_thr02.nii.gz)
    output_count="${#output_masks[@]}"
    [[ "${output_count}" == "66" ]] || fail "Expected 66 masks for ${CURRENT_SUBJECT}; found ${output_count}"
    mv "${CURRENT_MANIFEST}" "${FINAL_MANIFEST}"
    echo "Completed ${CURRENT_SUBJECT}: ${output_count} masks written to ${CURRENT_OUTPUT_DIR}"
done

echo "ROI registration completed for: ${normalised_subjects[*]}"
