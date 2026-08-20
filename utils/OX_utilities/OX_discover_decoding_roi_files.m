function [files, stems, manifest] = OX_discover_decoding_roi_files( ...
        roi_dirs, bilateral_only, max_rois)
%OX_DISCOVER_DECODING_ROI_FILES Index one or more selected ROI directories.

if ischar(roi_dirs) || (isstring(roi_dirs) && isscalar(roi_dirs))
    roi_dirs = {char(string(roi_dirs))};
end
assert(iscell(roi_dirs) && ~isempty(roi_dirs), ...
    'roi_dirs must contain at least one directory.');

entries = struct([]);
manifest = empty_manifest();
if bilateral_only
    patterns = {'*_bilateral_func_thr02.nii', '*_bilateral_func_thr02.nii.gz'};
else
    patterns = {'*_func_thr02.nii', '*_func_thr02.nii.gz'};
end

for directory_index = 1:numel(roi_dirs)
    roi_dir = char(string(roi_dirs{directory_index}));
    assert(isfolder(roi_dir), 'Missing selected ROI directory: %s', roi_dir);
    for pattern_index = 1:numel(patterns)
        found = dir(fullfile(roi_dir, patterns{pattern_index}));
        if isempty(entries)
            entries = found;
        else
            entries = [entries; found]; %#ok<AGROW>
        end
    end
    manifest = [manifest; read_normalized_manifest(roi_dir)]; %#ok<AGROW>
end

assert(~isempty(entries), 'No selected ROI masks found in: %s', ...
    strjoin(string(roi_dirs), ', '));
names = string({entries.name})';
stems = regexprep(names, '\.nii(\.gz)?$', '', 'ignorecase');
[unique_stems, ~, group] = unique(lower(stems));
duplicates = unique_stems(accumarray(group, 1) > 1);
assert(isempty(duplicates), 'Duplicate selected ROI basenames: %s', ...
    strjoin(duplicates, ', '));

[~, order] = sort(lower(names));
entries = entries(order);
names = names(order);
stems = stems(order);
if isfinite(max_rois) && max_rois < numel(entries)
    warning('MaxROIs=%d selects a deterministic smoke-test subset.', max_rois);
    keep = 1:min(max_rois, numel(entries));
    entries = entries(keep);
    names = names(keep);
    stems = stems(keep);
end
files = cellstr(fullfile(string({entries.folder})', names));
stems = cellstr(stems);
end

function manifest = read_normalized_manifest(roi_dir)
manifest_file = fullfile(roi_dir, 'roi_manifest.tsv');
if ~isfile(manifest_file)
    warning('ROI manifest is missing: %s', manifest_file);
    manifest = empty_manifest();
    return;
end
loaded = readtable(manifest_file, 'FileType', 'text', 'Delimiter', '\t', ...
    'TextType', 'string');
required = ["source", "contributing_labels", "hemisphere", "output_file"];
assert(all(ismember(required, string(loaded.Properties.VariableNames))), ...
    'ROI manifest lacks required columns: %s', manifest_file);
output_basename = regexprep(loaded.output_file, '^.*[/\\]', '');
output_stem = regexprep(output_basename, '\.nii(\.gz)?$', '', 'ignorecase');
manifest = table(string(loaded.source), string(loaded.contributing_labels), ...
    string(loaded.hemisphere), string(loaded.output_file), string(output_stem), ...
    'VariableNames', {'source', 'contributing_labels', 'hemisphere', ...
    'output_file', 'output_stem'});
end

function manifest = empty_manifest()
manifest = table(strings(0, 1), strings(0, 1), strings(0, 1), ...
    strings(0, 1), strings(0, 1), ...
    'VariableNames', {'source', 'contributing_labels', 'hemisphere', ...
    'output_file', 'output_stem'});
end
