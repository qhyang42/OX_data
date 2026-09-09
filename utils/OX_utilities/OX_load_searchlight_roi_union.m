function [roi_union, metadata, selection] = OX_load_searchlight_roi_union( ...
        nifti_dir, reference_header, requested)
%OX_LOAD_SEARCHLIGHT_ROI_UNION Resolve searchlight ROI names or mask paths.
% Empty input preserves the original seven-mask olfactory network. GLOBAL
% returns the full reference grid; the caller still applies GM, functional,
% and finite-data restrictions. Names resolve bilateral thr02 masks in the
% primary/secondary decoding directories. Explicit paths must match the grid.

if nargin < 3 || isempty(requested)
    requested = "";
end
names = strip(string(requested));
names = names(:);
assert(~any(ismissing(names)), 'ROI must not contain missing strings.');
is_default = isscalar(names) && strlength(names) == 0;
if is_default
    names = ["TU"; "AON"; "PirF"; "PirT"; "olfAMG"; "EC"; "HIPP"];
end
assert(all(strlength(names) > 0), 'ROI names/paths must be nonempty.');
assert(~any(strcmpi(names, 'global')) || isscalar(names), ...
    'ROI global cannot be combined with other ROIs.');
roi_union = false(reference_header.dim);
if isscalar(names) && strcmpi(names, 'global')
    roi_union(:) = true;
    metadata = table("global", "global", "", "", nnz(roi_union), ...
        'VariableNames', {'roi', 'selection_group', 'basename', 'file', ...
        'n_mask_voxels'});
    selection = "global";
    return;
end

roi_root = fullfile(nifti_dir, 'coreg', 'roi_decoding');
groups = strings(numel(names), 1);
basenames = strings(numel(names), 1);
files = strings(numel(names), 1);
n_mask_voxels = zeros(numel(names), 1);
for index = 1:numel(names)
    if is_default
        groups(index) = "primary";
        stem = names(index);
        if names(index) == "TU"
            groups(index) = "secondary";
            stem = "olf_TU";
        end
        files(index) = fullfile(roi_root, groups(index), ...
            stem + '_bilateral_func_thr02.nii');
        assert(isfile(files(index)), 'Missing ROI mask: %s', files(index));
    elseif isfile(names(index))
        files(index) = names(index);
        groups(index) = "custom";
    else
        candidates = strings(0, 1);
        candidate_groups = strings(0, 1);
        for group = ["primary", "secondary"]
            entries = dir(fullfile(roi_root, group, '*_bilateral_func_thr02.nii'));
            for entry = entries'
                stem = erase(string(entry.name), '_bilateral_func_thr02.nii');
                if strcmpi(stem, names(index)) || ...
                        strcmpi(regexprep(stem, '^olf_', ''), names(index))
                    candidates(end + 1, 1) = fullfile(entry.folder, entry.name); %#ok<AGROW>
                    candidate_groups(end + 1, 1) = group; %#ok<AGROW>
                end
            end
        end
        assert(isscalar(candidates), ...
            'ROI %s resolved to %d masks; pass an existing explicit NIfTI path.', ...
            names(index), numel(candidates));
        files(index) = candidates(1);
        groups(index) = candidate_groups(1);
    end
    [~, base, extension] = fileparts(files(index));
    basenames(index) = base + extension;
    header = spm_vol(char(files(index)));
    assert(isscalar(header), 'ROI must be a single 3-D mask: %s', files(index));
    assert(isequal(header.dim, reference_header.dim), ...
        'Dimension mismatch for ROI %s.', files(index));
    assert(max(abs(header.mat(:) - reference_header.mat(:))) < 1e-4, ...
        'Affine mismatch for ROI %s.', files(index));
    values = spm_read_vols(header);
    mask = isfinite(values) & values > 0;
    n_mask_voxels(index) = nnz(mask);
    roi_union = roi_union | mask;
end
metadata = table(names, groups, basenames, files, n_mask_voxels, ...
    'VariableNames', {'roi', 'selection_group', 'basename', 'file', ...
    'n_mask_voxels'});
if is_default
    selection = "olfactory";
else
    % Include a path-derived digest to distinguish masks with equal basenames.
    digest = java.security.MessageDigest.getInstance('SHA-256');
    digest.update(unicode2native(char(strjoin(files, newline)), 'UTF-8'));
    key = lower(reshape(dec2hex(typecast(digest.digest(), 'uint8'), 2)', 1, []));
    selection = "roi_" + string(key(1:16));
end
end
