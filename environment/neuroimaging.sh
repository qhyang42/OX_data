#!/bin/zsh
# Shared MATLAB, FSL, and FreeSurfer environment for the OX_DATA project.
#
# Source this file rather than executing it:
#   source /Users/qhyang/Desktop/OX_DATA/environment/neuroimaging.sh
#
# Set any of the variables below before sourcing to override the defaults.

if [[ "${OX_NEURO_ENV_LOADED:-}" == "1" ]]; then
    return 0 2>/dev/null || exit 0
fi

_ox_env_fail() {
    print -u2 -- "OX_DATA environment error: $1"
    return 1
}

export MATLAB_ROOT="${MATLAB_ROOT:-/Applications/MATLAB_R2025b.app}"
export MATLAB_BIN="${MATLAB_BIN:-${MATLAB_ROOT}/bin/matlab}"
export FSLDIR="${FSLDIR:-/Users/qhyang/fsl}"
export FREESURFER_HOME="${FREESURFER_HOME:-/Applications/freesurfer}"
export SUBJECTS_DIR="${SUBJECTS_DIR:-${FREESURFER_HOME}/subjects}"
export FS_LICENSE="${FS_LICENSE:-${FREESURFER_HOME}/license.txt}"

[[ -x "${MATLAB_BIN}" ]] || {
    _ox_env_fail "MATLAB executable not found: ${MATLAB_BIN}"
    return 1 2>/dev/null || exit 1
}
[[ -f "${FSLDIR}/etc/fslconf/fsl.sh" ]] || {
    _ox_env_fail "FSL setup script not found under: ${FSLDIR}"
    return 1 2>/dev/null || exit 1
}
[[ -f "${FREESURFER_HOME}/SetUpFreeSurfer.sh" ]] || {
    _ox_env_fail "FreeSurfer setup script not found under: ${FREESURFER_HOME}"
    return 1 2>/dev/null || exit 1
}
[[ -f "${FS_LICENSE}" ]] || {
    _ox_env_fail "FreeSurfer license not found: ${FS_LICENSE}"
    return 1 2>/dev/null || exit 1
}

# MATLAB should take precedence when more than one release is installed.
export PATH="${MATLAB_ROOT}/bin:${PATH}"

# These vendor scripts populate the remaining command paths and variables.
source "${FSLDIR}/etc/fslconf/fsl.sh" || {
    _ox_env_fail "FSL initialization failed"
    return 1 2>/dev/null || exit 1
}
source "${FREESURFER_HOME}/SetUpFreeSurfer.sh" >/dev/null || {
    _ox_env_fail "FreeSurfer initialization failed"
    return 1 2>/dev/null || exit 1
}

# Normalize command precedence after the vendor scripts modify PATH. zsh's
# unique path array also removes duplicates introduced by repeated setup.
typeset -U path PATH
path=(
    "${MATLAB_ROOT}/bin"
    "${FREESURFER_HOME}/bin"
    "${FREESURFER_HOME}/fsfast/bin"
    "${FREESURFER_HOME}/mni/bin"
    "${FSLDIR}/share/fsl/bin"
    $path
)

export OX_NEURO_ENV_LOADED=1
unset -f _ox_env_fail
