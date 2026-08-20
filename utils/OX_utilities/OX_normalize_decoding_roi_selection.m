function selection = OX_normalize_decoding_roi_selection(requested)
%OX_NORMALIZE_DECODING_ROI_SELECTION Validate an ROI-set selector.

selection = lower(strtrim(char(string(requested))));
allowed = {'primary', 'secondary', 'all', 'old'};
assert(isscalar(string(requested)) && ismember(selection, allowed), ...
    'ROISelection must be primary, secondary, all, or old.');
end
