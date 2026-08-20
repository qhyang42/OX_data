function [selection, roi_dirs] = OX_resolve_decoding_roi_selection( ...
        nifti_dir, requested_selection, roi_dir_override)
%OX_RESOLVE_DECODING_ROI_SELECTION Resolve selected decoding-mask folders.
%
% OLD uses old_rois (or a legacy flat ROIDir when old_rois is absent).
% PRIMARY and SECONDARY use their named subfolders. ALL combines both.

selection = OX_normalize_decoding_roi_selection(requested_selection);
if nargin < 3 || strlength(string(roi_dir_override)) == 0
    roi_root = fullfile(nifti_dir, 'coreg', 'roi_decoding');
else
    roi_root = char(string(roi_dir_override));
end
assert(isfolder(roi_root), 'Missing ROI directory: %s', roi_root);

switch selection
    case 'old'
        archived_dir = fullfile(roi_root, 'old_rois');
        if isfolder(archived_dir)
            roi_dirs = {archived_dir};
        else
            % Backward compatibility for an explicit legacy flat ROIDir.
            roi_dirs = {roi_root};
        end
    case 'primary'
        roi_dirs = {resolve_group_dir(roi_root, 'primary')};
    case 'secondary'
        roi_dirs = {resolve_group_dir(roi_root, 'secondary')};
    case 'all'
        roi_dirs = {resolve_group_dir(roi_root, 'primary'), ...
                    resolve_group_dir(roi_root, 'secondary')};
end

for index = 1:numel(roi_dirs)
    assert(isfolder(roi_dirs{index}), ...
        'ROISelection %s requires directory: %s', selection, roi_dirs{index});
end
end

function group_dir = resolve_group_dir(roi_root, group_name)
[~, root_name] = fileparts(roi_root);
if strcmpi(root_name, group_name)
    % Allow ROIDir to point directly at the selected group.
    group_dir = roi_root;
else
    group_dir = fullfile(roi_root, group_name);
end
end
