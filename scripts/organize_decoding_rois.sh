#!/usr/bin/env bash

# Organize the OX decoding masks into bilateral primary and secondary sets.
#
# On the first run, every file at the top level of each subject's existing
# coreg/roi_decoding directory is moved unchanged into old_rois. The active
# masks are then rebuilt from that archive. The script is safe to rerun: once
# old_rois exists, it remains the immutable source archive.

set -euo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
OFC_SOURCE="${ROOT_DIR}/ROIs/OFC_small_MNI152_1mm.nii.gz"
DEFAULT_SUBJECTS=(2 3 4 5 6)

# Active masks are consumed by SPM, which requires uncompressed NIfTI files.
export FSLOUTPUTTYPE=NIFTI

for command_name in convert_xfm flirt fslmaths fslstats fslval; do
    command -v "${command_name}" >/dev/null 2>&1 || {
        echo "Required FSL command not found: ${command_name}" >&2
        exit 1
    }
done

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

normalise_subject_id() {
    local input_id="$1"
    if [[ "${input_id}" == subj_* ]]; then
        echo "${input_id}"
    elif [[ "${input_id}" =~ ^[0-9]+$ ]]; then
        echo "subj_${input_id}"
    else
        fail "Subject IDs must be numeric or subj_<number>: ${input_id}"
    fi
}

resolve_old_mask() {
    local archive_dir="$1"
    local stem="$2"
    local candidates=("${archive_dir}/${stem}.nii" "${archive_dir}/${stem}.nii.gz")
    local candidate
    for candidate in "${candidates[@]}"; do
        if [[ -f "${candidate}" ]]; then
            echo "${candidate}"
            return 0
        fi
    done
    fail "Missing archived source mask ${stem} in ${archive_dir}"
}

assert_mask() {
    local image_file="$1"
    local reference_file="$2"
    local key
    # The sform is the authoritative voxel-to-world transform used by SPM for
    # these files. (The NIfTI quaternion can round at the sixth decimal place.)
    for key in dim1 dim2 dim3 pixdim1 pixdim2 pixdim3 \
        'sto_xyz:1' 'sto_xyz:2' 'sto_xyz:3' 'sto_xyz:4'; do
        [[ "$(fslval "${image_file}" "${key}")" == "$(fslval "${reference_file}" "${key}")" ]] ||
            fail "Geometry mismatch (${key}) between ${image_file} and ${reference_file}"
    done
    fslstats "${image_file}" -R | awk '$1 == 0 && $2 == 1 { ok = 1 } END { exit !ok }' ||
        fail "Mask is empty or nonbinary: ${image_file}"
}

write_manifest_header() {
    local manifest="$1"
    printf 'subject\tgroup\troi_id\tsource\tcontributing_labels\themisphere\toutput_file\tvoxel_count\n' > "${manifest}"
}

write_manifest_row() {
    local manifest="$1"
    local subject="$2"
    local group="$3"
    local roi_id="$4"
    local source="$5"
    local labels="$6"
    local output_file="$7"
    local voxels
    voxels="$(fslstats "${output_file}" -V | awk '{print $1}')"
    printf '%s\t%s\t%s\t%s\t%s\tbilateral\t%s\t%s\n' \
        "${subject}" "${group}" "${roi_id}" "${source}" "${labels}" \
        "$(basename "${output_file}")" "${voxels}" >> "${manifest}"
}

copy_binary_mask() {
    local source_file="$1"
    local output_file="$2"
    local reference_file="$3"
    fslmaths "${source_file}" -bin "${output_file}"
    assert_mask "${output_file}" "${reference_file}"
}

merge_binary_masks() {
    local output_file="$1"
    local reference_file="$2"
    shift 2
    (( $# >= 2 )) || fail "merge_binary_masks needs at least two input masks"
    local first_file="$1"
    shift
    local command=(fslmaths "${first_file}")
    local input_file
    for input_file in "$@"; do
        command+=(-add "${input_file}")
    done
    command+=(-bin "${output_file}")
    "${command[@]}"
    assert_mask "${output_file}" "${reference_file}"
}

subjects=()
if (( $# == 0 )); then
    subjects=("${DEFAULT_SUBJECTS[@]}")
else
    for subject in "$@"; do
        subjects+=("$(normalise_subject_id "${subject}")")
    done
fi

[[ -f "${OFC_SOURCE}" ]] || fail "Missing small OFC source mask: ${OFC_SOURCE}"

# Validate all subject inputs before changing the directory layout.
for subject in "${subjects[@]}"; do
    [[ "${subject}" == subj_* ]] || subject="$(normalise_subject_id "${subject}")"
    nifti_dir="${ROOT_DIR}/MRI/${subject}/nifti"
    coreg_dir="${nifti_dir}/coreg"
    roi_dir="${coreg_dir}/roi_decoding"
    [[ -d "${roi_dir}" ]] || fail "Missing ROI directory: ${roi_dir}"
    for matrix in std2T1mat T12wbmat wb2funcmat; do
        [[ -s "${coreg_dir}/${matrix}" ]] || fail "Missing transform: ${coreg_dir}/${matrix}"
    done
    functional_references=("${nifti_dir}"/func/mean2*.nii "${nifti_dir}"/func/mean2*.nii.gz)
    (( ${#functional_references[@]} == 1 )) ||
        fail "Expected one mean functional image for ${subject}; found ${#functional_references[@]}"

    source_dir="${roi_dir}/old_rois"
    if [[ ! -d "${source_dir}" ]]; then
        source_dir="${roi_dir}"
    fi
    required_stems=(
        olf_AON_bilateral_func_thr02
        olf_pirF_bilateral_func_thr02
        olf_pirT_bilateral_func_thr02
        amygsub_MeA_bilateral_func_thr02
        amygsub_ACo_bilateral_func_thr02
        amygsub_CeA_bilateral_func_thr02
        amygsub_PAC_bilateral_func_thr02
        bn_medial_anterior_pfc_bilateral_func_thr02
        bn_hippocampus_bilateral_func_thr02
        bn_parahippocampal_mtc_bilateral_func_thr02
        bn_temporal_pole_bilateral_func_thr02
        bn_insula_bilateral_func_thr02
        bn_nacc_bilateral_func_thr02
        bn_dlpfc_bilateral_func_thr02
        bn_lateral_frontopolar_pfc_bilateral_func_thr02
        bn_ofc_bilateral_func_thr02
        bn_amygdala_bilateral_func_thr02
        olf_TU_bilateral_func_thr02
        amygsub_BMA_bilateral_func_thr02
        amygsub_BLA_bilateral_func_thr02
        amygsub_LA_bilateral_func_thr02
        amygsub_PCo_bilateral_func_thr02
    )
    for stem in "${required_stems[@]}"; do
        resolve_old_mask "${source_dir}" "${stem}" >/dev/null
    done
done

WORK_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/ox_roi_groups.XXXXXX")"
trap 'rm -rf "${WORK_ROOT}"' EXIT

for subject in "${subjects[@]}"; do
    [[ "${subject}" == subj_* ]] || subject="$(normalise_subject_id "${subject}")"
    nifti_dir="${ROOT_DIR}/MRI/${subject}/nifti"
    coreg_dir="${nifti_dir}/coreg"
    roi_dir="${coreg_dir}/roi_decoding"
    archive_dir="${roi_dir}/old_rois"
    primary_dir="${roi_dir}/primary"
    secondary_dir="${roi_dir}/secondary"
    functional_references=("${nifti_dir}"/func/mean2*.nii "${nifti_dir}"/func/mean2*.nii.gz)
    func_ref="${functional_references[0]}"

    if [[ ! -d "${archive_dir}" ]]; then
        mkdir -p "${archive_dir}"
        original_files=("${roi_dir}"/*)
        for original_file in "${original_files[@]}"; do
            [[ -f "${original_file}" ]] && mv "${original_file}" "${archive_dir}/"
        done
    fi
    mkdir -p "${primary_dir}" "${secondary_dir}"

    primary_manifest="${primary_dir}/roi_manifest.tsv"
    secondary_manifest="${secondary_dir}/roi_manifest.tsv"
    write_manifest_header "${primary_manifest}"
    write_manifest_header "${secondary_manifest}"

    # Primary: requested aliases, all bilateral and in native functional space.
    while IFS='|' read -r output_name source_stem source_name labels; do
        output_file="${primary_dir}/${output_name}_bilateral_func_thr02.nii"
        copy_binary_mask "$(resolve_old_mask "${archive_dir}" "${source_stem}")" "${output_file}" "${func_ref}"
        write_manifest_row "${primary_manifest}" "${subject}" primary "${output_name}" "${source_name}" "${labels}" "${output_file}"
    done <<'PRIMARY_MASKS'
AON|olf_AON_bilateral_func_thr02|in_house_olfactory_roi|AON
PirF|olf_pirF_bilateral_func_thr02|in_house_olfactory_roi|pirF
PirT|olf_pirT_bilateral_func_thr02|in_house_olfactory_roi|pirT
maPFC|bn_medial_anterior_pfc_bilateral_func_thr02|Brainnetome|L:1,11,13;R:2,12,14
HIPP|bn_hippocampus_bilateral_func_thr02|Brainnetome|L:215,217;R:216,218
paraHIPP|bn_parahippocampal_mtc_bilateral_func_thr02|Brainnetome|L:109,111,113,115,117,119;R:110,112,114,116,118,120
TP|bn_temporal_pole_bilateral_func_thr02|Brainnetome|L:69,77;R:70,78
insula|bn_insula_bilateral_func_thr02|Brainnetome|L:163,165,167,169,171,173;R:164,166,168,170,172,174
nACC|bn_nacc_bilateral_func_thr02|Brainnetome|L:223;R:224
PRIMARY_MASKS

    olf_amg_file="${primary_dir}/olfAMG_bilateral_func_thr02.nii"
    merge_binary_masks "${olf_amg_file}" "${func_ref}" \
        "$(resolve_old_mask "${archive_dir}" amygsub_MeA_bilateral_func_thr02)" \
        "$(resolve_old_mask "${archive_dir}" amygsub_ACo_bilateral_func_thr02)" \
        "$(resolve_old_mask "${archive_dir}" amygsub_CeA_bilateral_func_thr02)" \
        "$(resolve_old_mask "${archive_dir}" amygsub_PAC_bilateral_func_thr02)"
    write_manifest_row "${primary_manifest}" "${subject}" primary olfAMG \
        in_house_amygdala_subregion 'MeA+ACo+CeA+PAC' "${olf_amg_file}"

    subject_work="${WORK_ROOT}/${subject}"
    mkdir -p "${subject_work}"
    std2wb_mat="${subject_work}/std2wb.mat"
    std2func_mat="${subject_work}/std2func.mat"
    ofc_resampled="${subject_work}/olfOFC_resampled.nii"
    convert_xfm -omat "${std2wb_mat}" -concat "${coreg_dir}/T12wbmat" "${coreg_dir}/std2T1mat"
    convert_xfm -omat "${std2func_mat}" -concat "${coreg_dir}/wb2funcmat" "${std2wb_mat}"
    flirt -in "${OFC_SOURCE}" -ref "${func_ref}" -applyxfm -init "${std2func_mat}" \
        -interp trilinear -out "${ofc_resampled}"
    olf_ofc_file="${primary_dir}/olfOFC_bilateral_func_thr02.nii"
    fslmaths "${ofc_resampled}" -thr 0.2 -bin "${olf_ofc_file}"
    assert_mask "${olf_ofc_file}" "${func_ref}"
    write_manifest_row "${primary_manifest}" "${subject}" primary olfOFC \
        OFC_small_MNI152_1mm 'bilateral small OFC; trilinear resampling; threshold=0.2' "${olf_ofc_file}"

    # Secondary: the remaining conceptual bilateral decoding masks. Amygdala
    # subregion components are represented by the requested nonolfAMG union.
    while IFS='|' read -r output_name source_stem source_name labels; do
        output_file="${secondary_dir}/${output_name}_bilateral_func_thr02.nii"
        copy_binary_mask "$(resolve_old_mask "${archive_dir}" "${source_stem}")" "${output_file}" "${func_ref}"
        write_manifest_row "${secondary_manifest}" "${subject}" secondary "${output_name}" "${source_name}" "${labels}" "${output_file}"
    done <<'SECONDARY_MASKS'
bn_dlpfc|bn_dlpfc_bilateral_func_thr02|Brainnetome|L:5,15,19,21;R:6,16,20,22
bn_lateral_frontopolar_pfc|bn_lateral_frontopolar_pfc_bilateral_func_thr02|Brainnetome|L:27;R:28
bn_ofc|bn_ofc_bilateral_func_thr02|Brainnetome|L:41,43,45,47,49,51;R:42,44,46,48,50,52
bn_amygdala|bn_amygdala_bilateral_func_thr02|Brainnetome|L:211,213;R:212,214
olf_TU|olf_TU_bilateral_func_thr02|in_house_olfactory_roi|TU
SECONDARY_MASKS

    nonolf_amg_file="${secondary_dir}/nonolfAMG_bilateral_func_thr02.nii"
    merge_binary_masks "${nonolf_amg_file}" "${func_ref}" \
        "$(resolve_old_mask "${archive_dir}" amygsub_BMA_bilateral_func_thr02)" \
        "$(resolve_old_mask "${archive_dir}" amygsub_BLA_bilateral_func_thr02)" \
        "$(resolve_old_mask "${archive_dir}" amygsub_LA_bilateral_func_thr02)" \
        "$(resolve_old_mask "${archive_dir}" amygsub_PCo_bilateral_func_thr02)"
    write_manifest_row "${secondary_manifest}" "${subject}" secondary nonolfAMG \
        in_house_amygdala_subregion 'BMA+BLA+LA+PCo' "${nonolf_amg_file}"

    primary_masks=("${primary_dir}"/*_bilateral_func_thr02.nii)
    secondary_masks=("${secondary_dir}"/*_bilateral_func_thr02.nii)
    (( ${#primary_masks[@]} == 11 )) || fail "Expected 11 primary masks for ${subject}; found ${#primary_masks[@]}"
    (( ${#secondary_masks[@]} == 6 )) || fail "Expected 6 secondary masks for ${subject}; found ${#secondary_masks[@]}"
    echo "${subject}: archived old ROIs; wrote 11 primary and 6 secondary bilateral masks"
done

echo "ROI organization completed."
