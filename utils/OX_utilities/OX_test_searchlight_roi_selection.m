function OX_test_searchlight_roi_selection()
%OX_TEST_SEARCHLIGHT_ROI_SELECTION Synthetic NIfTI ROI selection regression.
root = tempname;
mkdir(root);
cleanup = onCleanup(@() rmdir(root, 's'));
header = struct('fname', '', 'dim', [7 7 7], 'mat', diag([2 2 2 1]), ...
    'dt', [2 0], 'pinfo', [1; 0; 0]);
names = ["olf_TU", "AON", "PirF", "PirT", "olfAMG", "EC", "HIPP"];
expected = false(header.dim);
for index = 1:numel(names)
    group = 'primary';
    if index == 1
        group = 'secondary';
    end
    folder = fullfile(root, 'coreg', 'roi_decoding', group);
    if ~isfolder(folder), mkdir(folder); end
    header.fname = char(fullfile(folder, names(index) + '_bilateral_func_thr02.nii'));
    mask = false(header.dim);
    mask(index, 3:5, 3:5) = true;
    expected = expected | mask;
    spm_write_vol(header, double(mask));
end
[mask, metadata, label] = OX_load_searchlight_roi_union(root, header);
assert(isequal(mask, expected) && height(metadata) == 7 && label == "olfactory");
empty_mask = OX_load_searchlight_roi_union(root, header, []);
assert(isequal(empty_mask, expected));
[selected, metadata] = OX_load_searchlight_roi_union(root, header, ["TU", "PirF"]);
assert(nnz(selected) == 18 && height(metadata) == 2);
[explicit, ~, label] = OX_load_searchlight_roi_union(root, header, metadata.file);
assert(isequal(selected, explicit) && startsWith(label, 'roi_'));
% GLOBAL must work even when no ROI directory exists.
[global_mask, ~, label] = OX_load_searchlight_roi_union('nonexistent', header, 'GLOBAL');
assert(all(global_mask(:)) && label == "global");
gm = true(header.dim); gm(1, :, :) = false;
functional = true(header.dim); functional(:, 1, :) = false;
finite = true(header.dim); finite(:, :, 1) = false;
final = global_mask & gm & functional & finite;
assert(nnz(final) == 216 && any(final(:) & ~expected(:)));
spheres = OX_build_masked_searchlights(final, header.mat, 4, 2);
assert(all(final(spheres.center_linear_indices)));
for index = 1:numel(spheres.neighborhoods)
    assert(all(final(spheres.feature_linear_indices(spheres.neighborhoods{index}))));
end
must_fail(@() OX_load_searchlight_roi_union(root, header, ["global", "TU"]));
must_fail(@() OX_load_searchlight_roi_union(root, header, 'unknown'));
wrong = header; wrong.mat(1, 4) = 10;
must_fail(@() OX_load_searchlight_roi_union(root, wrong, 'TU'));
wrong = header; wrong.dim = [8 7 7];
must_fail(@() OX_load_searchlight_roi_union(root, wrong, 'TU'));
fprintf('Searchlight ROI selection tests passed.\n');
end

function must_fail(action)
failed = false;
try
    action();
catch
    failed = true;
end
assert(failed, 'Invalid ROI selection or geometry was accepted.');
end
