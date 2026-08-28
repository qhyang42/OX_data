function results = OX_roi_context_split_half_similarity(subjidx, varargin)
%OX_ROI_CONTEXT_SPLIT_HALF_SIMILARITY Semantic context-pattern stability.
%
%   results = OX_roi_context_split_half_similarity(subjidx, Name, Value, ...)
%
% Loads the physio-regressed odor-aligned GLMsingle estimates, restricts
% features to the requested bilateral ROI, gray matter, and uncorrected
% Odor > Rest functional mask, and runs repeated independent-run split-half
% similarity for PERSON, FOOD, and LOCATION. CONTROL trials are excluded.
%
% ROISelection is compatible with the existing decoding ROI resolver.
% ROINames defaults to PirF, PirT, AON, olfOFC, and olfAMG and can be expanded.

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'subjidx', @(x) isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x == round(x));
addParameter(p, 'MRIRoot', '', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'CueRoot', '', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'ROISelection', 'primary', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'ROIDir', '', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'ROINames', ["PirF", "PirT", "AON", "olfOFC", "olfAMG"], ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(p, 'MinVoxels', 10, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 2 && x == round(x));
addParameter(p, 'NumSplits', 200, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == round(x));
addParameter(p, 'MaxSplitAttempts', 100000, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == round(x));
addParameter(p, 'RandomSeed', 1, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'NumPermutations', 5000, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0 && x == round(x));
addParameter(p, 'PermutationSeed', 1001, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'UseParallel', true, ...
    @(x) islogical(x) && isscalar(x));
addParameter(p, 'PermutationBatchSize', 250, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == round(x));
addParameter(p, 'MaxPermutationLabelAttempts', 1000, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == round(x));
addParameter(p, 'StratifyBySession', true, ...
    @(x) islogical(x) && isscalar(x));
addParameter(p, 'OutputDir', '', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'SaveOutputs', true, ...
    @(x) islogical(x) && isscalar(x));
parse(p, subjidx, varargin{:});
opts = p.Results;
opts.ROISelection = OX_normalize_decoding_roi_selection(opts.ROISelection);
opts.ROINames = string(opts.ROINames);
opts.ROINames = opts.ROINames(:);
assert(~isempty(opts.ROINames) && all(strlength(opts.ROINames) > 0), ...
    'ROINames must contain at least one nonempty name.');
assert(numel(unique(lower(opts.ROINames))) == numel(opts.ROINames), ...
    'ROINames contains duplicates.');

subject_real_names = {'240711_fMRI_OX_NWU_AS', ...
                      '240723_fMRI_OX_NWU_LS', ...
                      '240814_fMRI_OX_NWU_JN', ...
                      '240816_fMRI_OX_NWU_RR', ...
                      '241018_fMRI_OX_NWU_BN', ...
                      '250117_fMRI_OX_NWU_VS'};
assert(subjidx >= 1 && subjidx <= numel(subject_real_names), ...
    'subjidx must be between 1 and %d.', numel(subject_real_names));

analysis_tic = tic;
subject_name = sprintf('subj_%d', subjidx);
mri_root = resolve_mri_root(opts.MRIRoot);
nifti_dir = fullfile(mri_root, subject_name, 'nifti');
beta_dir = fullfile(nifti_dir, 'sniff_single_trial_by_category_physio');
fit_file = fullfile(beta_dir, 'TYPED_FITHRF_GLMDENOISE_RR.mat');
gm_mask_file = fullfile(nifti_dir, 'coreg', 'gm_mask_thr05_func.nii');
functional_mask_file = fullfile(nifti_dir, ...
    'first_level_model_sniff_physio', ...
    'spmT_0001_uncorrected_p001.nii');
[roi_selection, roi_dirs] = OX_resolve_decoding_roi_selection( ...
    nifti_dir, opts.ROISelection, opts.ROIDir);

assert(isfile(fit_file), 'Missing GLMsingle file: %s', fit_file);
assert(isfile(gm_mask_file), 'Missing gray-matter mask: %s', gm_mask_file);
assert(isfile(functional_mask_file), ...
    'Missing functional restriction mask: %s', functional_mask_file);
assert(exist('spm_vol', 'file') == 2 && exist('spm_read_vols', 'file') == 2, ...
    'SPM must be on the MATLAB path. Run setup_ox first.');

fprintf('\n[CONTEXT SPLIT-HALF] Loading %s\n', fit_file);
loaded = load(fit_file, 'modelmd');
assert(isfield(loaded, 'modelmd'), 'GLMsingle file lacks modelmd: %s', fit_file);
modelmd = squeeze(loaded.modelmd);
clear loaded
assert(ismatrix(modelmd), ...
    'squeeze(modelmd) must produce a two-dimensional matrix.');

metadata_args = {};
if strlength(string(opts.CueRoot)) > 0
    metadata_args = {'CueRoot', char(string(opts.CueRoot))};
end
all_trial_metadata = OX_load_trial_metadata(subjidx, metadata_args{:});
assert(size(modelmd, 2) == height(all_trial_metadata), ...
    'modelmd has %d trials but metadata has %d.', ...
    size(modelmd, 2), height(all_trial_metadata));
if strlength(string(opts.CueRoot)) == 0 && exist('OX_get_odor', 'file') == 2
    [legacy_odors, legacy_contexts] = OX_get_odor(subject_name);
    assert(isequal(double(legacy_odors(:)), all_trial_metadata.odor) && ...
        isequal(upper(strtrim(string(legacy_contexts(:)))), ...
        all_trial_metadata.context), ...
        'Exact trial metadata order differs from OX_get_odor.');
end

context_order = ["PERSON", "FOOD", "LOCATION"];
semantic_mask = ismember(all_trial_metadata.context, context_order);
assert(nnz(semantic_mask) == 600, ...
    'Expected 600 semantic trials; found %d.', nnz(semantic_mask));
trial_metadata = all_trial_metadata(semantic_mask, :);
modelmd = modelmd(:, semantic_mask);
% To honor the explicit CONTROL exclusion, even the run nuisance mean is
% estimated only from the three semantic contexts.
modelmd = center_voxels_within_runs(modelmd, trial_metadata.run_id);

fprintf('Subject %s | %d semantic trials | %d exact metadata run IDs\n', ...
    subject_name, height(trial_metadata), ...
    numel(unique(trial_metadata.run_id)));
fprintf(['Preprocessing: semantic-only within-run voxel centering; ' ...
    'Pearson spatial demeaning; no voxel scaling or whitening.\n']);
fprintf('Functional restriction: %s\n', functional_mask_file);
fprintf('ROI selection: %s | requested: %s\n', ...
    roi_selection, strjoin(opts.ROINames, ', '));

gm_header = spm_vol(gm_mask_file);
gm_mask = spm_read_vols(gm_header) > 0;
functional_header = spm_vol(functional_mask_file);
assert_same_geometry(functional_header, gm_header, ...
    functional_mask_file, gm_mask_file);
functional_mask = spm_read_vols(functional_header) > 0;
restriction_mask = gm_mask & functional_mask;
gm_indices = find(gm_mask);
assert(size(modelmd, 1) == numel(gm_indices), ...
    'modelmd rows (%d) do not match gray-matter voxels (%d).', ...
    size(modelmd, 1), numel(gm_indices));
model_index_volume = zeros(size(gm_mask), 'uint32');
model_index_volume(gm_indices) = uint32(1:numel(gm_indices));

[all_roi_files, all_roi_stems, manifest] = ...
    OX_discover_decoding_roi_files(roi_dirs, true, Inf);
[roi_files, roi_stems, roi_names] = select_requested_rois( ...
    all_roi_files, all_roi_stems, opts.ROINames);
n_rois = numel(roi_files);

roi_features = cell(n_rois, 1);
roi_sources = strings(n_rois, 1);
roi_labels = strings(n_rois, 1);
n_mask_voxels = zeros(n_rois, 1);
n_gm_overlap = zeros(n_rois, 1);
n_functional_overlap = zeros(n_rois, 1);
n_features = zeros(n_rois, 1);
status = repmat("ok", n_rois, 1);

for roi_idx = 1:n_rois
    roi_header = spm_vol(roi_files{roi_idx});
    assert_same_geometry(roi_header, gm_header, roi_files{roi_idx}, gm_mask_file);
    roi_mask = spm_read_vols(roi_header) > 0;
    n_mask_voxels(roi_idx) = nnz(roi_mask);
    n_gm_overlap(roi_idx) = nnz(roi_mask & gm_mask);
    n_functional_overlap(roi_idx) = nnz(roi_mask & restriction_mask);
    feature_indices = double(model_index_volume(roi_mask & restriction_mask));
    feature_indices = feature_indices(feature_indices > 0);
    if ~isempty(feature_indices)
        feature_indices = feature_indices( ...
            all(isfinite(modelmd(feature_indices, :)), 2));
    end
    feature_indices = unique(feature_indices(:)', 'stable');
    roi_features{roi_idx} = feature_indices;
    n_features(roi_idx) = numel(feature_indices);
    [roi_sources(roi_idx), roi_labels(roi_idx)] = ...
        manifest_metadata(manifest, roi_stems{roi_idx});
    if n_features(roi_idx) < opts.MinVoxels
        status(roi_idx) = "insufficient_features";
    end
end
assert(any(status == "ok"), ...
    'No requested ROI contains at least %d usable voxels.', opts.MinVoxels);

roi_results = repmat(struct(), n_rois, 1);
roi_X = cell(n_rois, 1);
reference_splits = [];
for roi_idx = 1:n_rois
    roi_results(roi_idx).roi_name = char(roi_names(roi_idx));
    roi_results(roi_idx).mask_file = roi_files{roi_idx};
    roi_results(roi_idx).n_voxels = n_features(roi_idx);
    roi_results(roi_idx).status = char(status(roi_idx));
    if status(roi_idx) ~= "ok"
        roi_results(roi_idx).split_half = [];
        continue;
    end
    fprintf('  ROI %s: %d usable voxels | %d splits\n', ...
        roi_names(roi_idx), n_features(roi_idx), opts.NumSplits);
    X = double(modelmd(roi_features{roi_idx}, :)');
    roi_X{roi_idx} = X;
    split_args = {'NumSplits', opts.NumSplits, ...
        'RandomSeed', opts.RandomSeed, ...
        'MaxSplitAttempts', opts.MaxSplitAttempts, ...
        'ExpectedOdors', 20};
    if opts.StratifyBySession
        split_args = [split_args, ...
            {'StratifyGroups', trial_metadata.session_id}]; %#ok<AGROW>
    end
    roi_results(roi_idx).split_half = ...
        OX_context_split_half_similarity(X, trial_metadata.odor, ...
        trial_metadata.context, trial_metadata.run_id, split_args{:});
    if isempty(reference_splits)
        reference_splits = ...
            roi_results(roi_idx).split_half.split_run_membership;
    else
        assert(isequal(reference_splits, ...
            roi_results(roi_idx).split_half.split_run_membership), ...
            'ROIs did not use identical run splits despite the fixed seed.');
    end
end

permutation = run_context_permutations(roi_X, status == "ok", ...
    trial_metadata, opts, subjidx);

subject_column = repmat(subjidx, n_rois, 1);
voxel_count = n_features;
n_valid_splits = zeros(n_rois, 1);
mean_diagonal_z = nan(n_rois, 1);
mean_off_diagonal_z = nan(n_rois, 1);
mean_delta_z = nan(n_rois, 1);
max_delta_z = nan(n_rois, 1);
sd_delta_z = nan(n_rois, 1);
p_permutation_mean_delta_z = nan(n_rois, 1);
p_permutation_max_delta_z = nan(n_rois, 1);
for roi_idx = 1:n_rois
    if isempty(roi_results(roi_idx).split_half)
        continue;
    end
    item = roi_results(roi_idx).split_half;
    n_valid_splits(roi_idx) = item.n_valid_splits;
    mean_diagonal_z(roi_idx) = item.aggregate.mean_diagonal_z;
    mean_off_diagonal_z(roi_idx) = item.aggregate.mean_off_diagonal_z;
    mean_delta_z(roi_idx) = item.aggregate.mean_delta_z;
    max_delta_z(roi_idx) = item.aggregate.max_delta_z;
    sd_delta_z(roi_idx) = item.aggregate.sd_delta_z;
    if opts.NumPermutations > 0
        p_permutation_mean_delta_z(roi_idx) = empirical_upper_p( ...
            mean_delta_z(roi_idx), ...
            permutation.roi_mean_delta_z(roi_idx, :));
        p_permutation_max_delta_z(roi_idx) = empirical_upper_p( ...
            max_delta_z(roi_idx), ...
            permutation.roi_max_split_delta_z(roi_idx, :));
    end
end
summary = table(subject_column, roi_names, voxel_count, n_valid_splits, ...
    mean_diagonal_z, mean_off_diagonal_z, mean_delta_z, max_delta_z, ...
    sd_delta_z, p_permutation_mean_delta_z, p_permutation_max_delta_z, ...
    status, ...
    'VariableNames', {'subject', 'roi', 'voxel_count', 'n_valid_splits', ...
    'mean_diagonal_z', 'mean_off_diagonal_z', 'mean_delta_z', ...
    'max_delta_z', 'sd_delta_z', 'p_permutation_mean_delta_z', ...
    'p_permutation_max_delta_z', 'status'});

n_contexts = numel(context_order);
n_context_rows = n_rois * n_contexts;
context_subject = repmat(subjidx, n_context_rows, 1);
context_roi = strings(n_context_rows, 1);
context_name = strings(n_context_rows, 1);
context_delta_z = nan(n_context_rows, 1);
context_p_permutation = nan(n_context_rows, 1);
context_voxel_count = nan(n_context_rows, 1);
context_status = strings(n_context_rows, 1);
context_row = 0;
for roi_idx = 1:n_rois
    for context_idx = 1:n_contexts
        context_row = context_row + 1;
        context_roi(context_row) = roi_names(roi_idx);
        context_name(context_row) = context_order(context_idx);
        context_voxel_count(context_row) = n_features(roi_idx);
        context_status(context_row) = status(roi_idx);
        if isempty(roi_results(roi_idx).split_half)
            continue;
        end
        context_delta_z(context_row) = ...
            roi_results(roi_idx).split_half.aggregate.mean_context_delta_z( ...
            context_idx);
        if opts.NumPermutations > 0
            context_p_permutation(context_row) = empirical_upper_p( ...
                context_delta_z(context_row), squeeze( ...
                permutation.roi_context_mean_delta_z( ...
                roi_idx, context_idx, :)));
        end
    end
end
context_summary = table(context_subject, context_roi, context_name, ...
    context_voxel_count, context_delta_z, context_p_permutation, ...
    repmat(opts.NumPermutations, n_context_rows, 1), context_status, ...
    'VariableNames', {'subject', 'roi', 'context', 'voxel_count', ...
    'context_delta_z', 'p_context_permutation', 'n_permutations', 'status'});

roi_file = string(roi_files(:));
roi_metadata = table(roi_names, roi_sources, roi_labels, roi_file, status, ...
    n_mask_voxels, n_gm_overlap, n_functional_overlap, n_features, ...
    'VariableNames', {'roi_name', 'source', 'contributing_labels', ...
    'mask_file', 'status', 'n_mask_voxels', 'n_gm_overlap', ...
    'n_functional_gm_overlap', 'n_features_used'});

if strlength(string(opts.OutputDir)) == 0
    output_dir = fullfile(mri_root, 'group', ...
        'context_split_half_similarity_physio');
else
    output_dir = char(string(opts.OutputDir));
end

results = struct();
results.subject = struct('index', subjidx, 'name', subject_name, ...
    'name_real', subject_real_names{subjidx});
results.analysis = struct( ...
    'name', 'semantic_context_cross_half_pattern_similarity', ...
    'contexts', {{'PERSON', 'FOOD', 'LOCATION'}}, ...
    'context_effect', ['mean Fisher-z of 3 diagonal cells minus mean ' ...
        'Fisher-z of 6 off-diagonal cells'], ...
    'split_unit', 'acquisition run', ...
    'runs_per_half', 40, ...
    'session_stratified', opts.StratifyBySession, ...
    'template_weighting', ['mean trials within context-by-odor, then ' ...
        'equal mean across 20 odors'], ...
    'similarity', 'Pearson correlation across voxels', ...
    'inference', ['semantic context labels permuted within exact run; ' ...
        'the complete repeated split-half statistic is rerun']);
results.preprocessing = struct( ...
    'within_run_voxel_centering', true, ...
    'within_run_centering_trials', 'PERSON, FOOD, and LOCATION only', ...
    'demean_spatial_patterns', true, ...
    'scale_voxels', false, ...
    'whiten_voxels', false, ...
    'control_trials_used', false, ...
    'gray_matter_mask', gm_mask_file, ...
    'functional_mask', functional_mask_file);
results.inputs = struct('fit_file', fit_file, ...
    'roi_selection', roi_selection, 'roi_directories', {roi_dirs});
results.task_structure = summarize_task_structure(all_trial_metadata);
results.trial_metadata = trial_metadata;
results.roi_metadata = roi_metadata;
results.roi_results = roi_results;
results.null = permutation;
results.summary = summary;
results.context_summary = context_summary;
results.dimension_order = struct( ...
    'split_matrix', {{'half_a_context', 'half_b_context', 'split'}}, ...
    'cell_trial_counts', {{'split', 'half', 'context', 'odor'}});
results.options = opts;
results.output_dir = output_dir;
results.runtime_minutes = toc(analysis_tic) / 60;
results.matlab_version = version;

if opts.SaveOutputs
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end
    output_file = fullfile(output_dir, sprintf( ...
        'context_split_half_similarity_subj%d_results.mat', subjidx));
    save(output_file, 'results', '-v7.3');
    fprintf('Saved subject result to %s\n', output_file);
end
end

function permutation = run_context_permutations(roi_X, valid_roi, ...
        trial_metadata, opts, subjidx)
n_rois = numel(roi_X);
n_permutations = opts.NumPermutations;
permutation = struct( ...
    'n_permutations', n_permutations, ...
    'roi_mean_diagonal_z', nan(n_rois, n_permutations), ...
    'roi_mean_off_diagonal_z', nan(n_rois, n_permutations), ...
    'roi_mean_delta_z', nan(n_rois, n_permutations), ...
    'roi_context_mean_delta_z', nan(n_rois, 3, n_permutations), ...
    'roi_max_split_delta_z', nan(n_rois, n_permutations), ...
    'n_attempted_splits', zeros(1, n_permutations), ...
    'n_rejected_coverage', zeros(1, n_permutations), ...
    'n_rejected_duplicate', zeros(1, n_permutations), ...
    'n_label_shuffle_attempts', zeros(1, n_permutations), ...
    'context_permutation_seed_used', nan(1, n_permutations));
permutation.context_permutation_seed_initial = opts.PermutationSeed + ...
    subjidx * 100000 + (1:n_permutations);
permutation.split_sampling_seed = opts.RandomSeed + ...
    1000000 + (1:n_permutations);
context_permutation_seeds = permutation.context_permutation_seed_initial;
split_sampling_seeds = permutation.split_sampling_seed;
permutation.exchangeability = ['PERSON/FOOD/LOCATION labels shuffled ' ...
    'within exact acquisition run; odor, run, and imaging data fixed; ' ...
    'shuffles without global 3-by-20 coverage are resampled'];
permutation.max_delta_definition = ['for each subject and ROI separately, ' ...
    'maximum delta-z across that permutation''s valid repeated splits'];
permutation.roi_pooling = false;
if n_permutations == 0
    permutation.runtime_minutes = 0;
    return;
end

use_parallel = opts.UseParallel && license('test', 'Distrib_Computing_Toolbox');
if opts.UseParallel && ~use_parallel
    warning('Parallel Computing Toolbox unavailable; permutations are serial.');
end
if use_parallel && isempty(gcp('nocreate'))
    parpool;
end
if use_parallel
    roi_X_constant = parallel.pool.Constant(roi_X);
    metadata_constant = parallel.pool.Constant(trial_metadata);
    constant_cleanup = onCleanup(@() delete_parallel_constants( ...
        roi_X_constant, metadata_constant));
end

fprintf('Running %d within-run semantic-context permutations...\n', ...
    n_permutations);
permutation_tic = tic;
batch_size = min(opts.PermutationBatchSize, n_permutations);
for batch_start = 1:batch_size:n_permutations
    batch_indices = batch_start:min(batch_start + batch_size - 1, n_permutations);
    batch_context_seeds = context_permutation_seeds(batch_indices);
    batch_split_seeds = split_sampling_seeds(batch_indices);
    batch_results = cell(numel(batch_indices), 1);
    if use_parallel
        parfor batch_idx = 1:numel(batch_indices)
            batch_results{batch_idx} = evaluate_context_permutation( ...
                roi_X_constant.Value, valid_roi, metadata_constant.Value, opts, ...
                batch_context_seeds(batch_idx), ...
                batch_split_seeds(batch_idx));
        end
    else
        for batch_idx = 1:numel(batch_indices)
            batch_results{batch_idx} = evaluate_context_permutation( ...
                roi_X, valid_roi, trial_metadata, opts, ...
                batch_context_seeds(batch_idx), ...
                batch_split_seeds(batch_idx));
        end
    end
    for batch_idx = 1:numel(batch_indices)
        permutation_idx = batch_indices(batch_idx);
        item = batch_results{batch_idx};
        permutation.roi_mean_diagonal_z(:, permutation_idx) = ...
            item.mean_diagonal_z;
        permutation.roi_mean_off_diagonal_z(:, permutation_idx) = ...
            item.mean_off_diagonal_z;
        permutation.roi_mean_delta_z(:, permutation_idx) = item.mean_delta_z;
        permutation.roi_context_mean_delta_z(:, :, permutation_idx) = ...
            item.context_mean_delta_z;
        permutation.roi_max_split_delta_z(:, permutation_idx) = ...
            item.max_split_delta_z;
        permutation.n_attempted_splits(permutation_idx) = ...
            item.n_attempted_splits;
        permutation.n_rejected_coverage(permutation_idx) = ...
            item.n_rejected_coverage;
        permutation.n_rejected_duplicate(permutation_idx) = ...
            item.n_rejected_duplicate;
        permutation.n_label_shuffle_attempts(permutation_idx) = ...
            item.n_label_shuffle_attempts;
        permutation.context_permutation_seed_used(permutation_idx) = ...
            item.context_permutation_seed_used;
    end
    fprintf('  permutations %d/%d | %.1f min\n', batch_indices(end), ...
        n_permutations, toc(permutation_tic) / 60);
end
permutation.runtime_minutes = toc(permutation_tic) / 60;
if use_parallel
    clear constant_cleanup
end
end

function delete_parallel_constants(varargin)
for constant_idx = 1:nargin
    delete(varargin{constant_idx});
end
end

function item = evaluate_context_permutation(roi_X, valid_roi, ...
        trial_metadata, opts, context_seed, split_seed)
permuted_contexts = strings(height(trial_metadata), 1);
label_attempt = 0;
coverage_ok = false;
valid_indices = find(valid_roi(:))';
assert(~isempty(valid_indices), 'No valid ROI available for permutation.');
first_result = [];
while ~coverage_ok && label_attempt < opts.MaxPermutationLabelAttempts
    label_attempt = label_attempt + 1;
    context_seed_used = context_seed + (label_attempt - 1) * 10000000;
    permuted_contexts = OX_permute_semantic_context_within_runs( ...
        trial_metadata.context, trial_metadata.run_id, context_seed_used);
    coverage_ok = has_split_eligible_context_odors(permuted_contexts, ...
        trial_metadata.odor, trial_metadata.run_id);
    if coverage_ok
        try
            first_result = score_permuted_roi(roi_X{valid_indices(1)}, ...
                trial_metadata, permuted_contexts, opts, split_seed);
        catch exception
            if strcmp(exception.identifier, 'OX_DATA:InsufficientValidSplits')
                coverage_ok = false;
            else
                rethrow(exception);
            end
        end
    end
end
assert(coverage_ok, ...
    ['Unable to generate a split-eligible context permutation in %d ' ...
     'attempts.'], ...
    opts.MaxPermutationLabelAttempts);
n_rois = numel(roi_X);
item = struct( ...
    'mean_diagonal_z', nan(n_rois, 1), ...
    'mean_off_diagonal_z', nan(n_rois, 1), ...
    'mean_delta_z', nan(n_rois, 1), ...
    'context_mean_delta_z', nan(n_rois, 3), ...
    'max_split_delta_z', nan(n_rois, 1), ...
    'n_attempted_splits', 0, ...
    'n_rejected_coverage', 0, ...
    'n_rejected_duplicate', 0, ...
    'n_label_shuffle_attempts', label_attempt, ...
    'context_permutation_seed_used', context_seed_used);
for valid_idx = 1:numel(valid_indices)
    roi_idx = valid_indices(valid_idx);
    if valid_idx == 1
        null_result = first_result;
    else
        null_result = score_permuted_roi(roi_X{roi_idx}, ...
            trial_metadata, permuted_contexts, opts, split_seed);
    end
    item.mean_diagonal_z(roi_idx) = ...
        null_result.aggregate.mean_diagonal_z;
    item.mean_off_diagonal_z(roi_idx) = ...
        null_result.aggregate.mean_off_diagonal_z;
    item.mean_delta_z(roi_idx) = null_result.aggregate.mean_delta_z;
    item.context_mean_delta_z(roi_idx, :) = ...
        null_result.aggregate.mean_context_delta_z;
    item.max_split_delta_z(roi_idx) = null_result.aggregate.max_delta_z;
    if item.n_attempted_splits == 0
        item.n_attempted_splits = null_result.n_attempted_splits;
        item.n_rejected_coverage = null_result.n_rejected_coverage;
        item.n_rejected_duplicate = null_result.n_rejected_duplicate;
    end
end
end

function result = score_permuted_roi(X, trial_metadata, ...
        permuted_contexts, opts, split_seed)
split_args = {'NumSplits', opts.NumSplits, ...
    'RandomSeed', split_seed, ...
    'MaxSplitAttempts', opts.MaxSplitAttempts, ...
    'ExpectedOdors', 20, ...
    'CollectDetails', false};
if opts.StratifyBySession
    split_args = [split_args, ...
        {'StratifyGroups', trial_metadata.session_id}];
end
result = OX_context_split_half_similarity(X, trial_metadata.odor, ...
    permuted_contexts, trial_metadata.run_id, split_args{:});
end

function complete = has_split_eligible_context_odors(contexts, odors, run_ids)
context_order = ["PERSON", "FOOD", "LOCATION"];
odor_values = unique(odors, 'sorted');
complete = true;
for context_idx = 1:3
    for odor_idx = 1:numel(odor_values)
        cell_mask = contexts == context_order(context_idx) & ...
            odors == odor_values(odor_idx);
        if numel(unique(run_ids(cell_mask))) < 2
            complete = false;
            return;
        end
    end
end
end

function p_value = empirical_upper_p(observed, null_values)
null_values = null_values(isfinite(null_values));
if isempty(null_values) || ~isfinite(observed)
    p_value = NaN;
else
    p_value = (1 + sum(null_values >= observed)) / ...
        (numel(null_values) + 1);
end
end

function mri_root = resolve_mri_root(requested)
if strlength(string(requested)) > 0
    mri_root = char(string(requested));
    assert(isfolder(mri_root), 'MRIRoot does not exist: %s', mri_root);
    return;
end
candidates = {'/Users/qhyang/Desktop/OX_DATA/MRI', ...
              '/Volumes/ExtremeSSD/OX_DATA/MRI'};
mri_root = '';
for candidate_idx = 1:numel(candidates)
    if isfolder(candidates{candidate_idx})
        mri_root = candidates{candidate_idx};
        break;
    end
end
assert(~isempty(mri_root), 'Could not auto-detect MRIRoot.');
end

function modelmd = center_voxels_within_runs(modelmd, run_ids)
runs = unique(run_ids, 'stable');
for run_idx = 1:numel(runs)
    mask = run_ids == runs(run_idx);
    assert(any(mask), 'Internal error: empty run encountered.');
    modelmd(:, mask) = modelmd(:, mask) - mean(modelmd(:, mask), 2);
end
end

function [files, stems, names] = select_requested_rois( ...
        all_files, all_stems, requested_names)
all_stems = string(all_stems(:));
canonical = regexprep(all_stems, '_bilateral_func_thr02$', '', 'ignorecase');
n_requested = numel(requested_names);
files = cell(n_requested, 1);
stems = cell(n_requested, 1);
names = strings(n_requested, 1);
for requested_idx = 1:n_requested
    match = find(strcmpi(canonical, requested_names(requested_idx)));
    assert(isscalar(match), ...
        'Requested ROI %s matched %d bilateral files. Available: %s', ...
        requested_names(requested_idx), numel(match), strjoin(canonical, ', '));
    files{requested_idx} = all_files{match};
    stems{requested_idx} = char(all_stems(match));
    names(requested_idx) = canonical(match);
end
end

function [source, labels] = manifest_metadata(manifest, roi_stem)
source = "unknown";
labels = "";
if isempty(manifest)
    return;
end
match = find(strcmpi(manifest.output_stem, string(roi_stem)));
assert(numel(match) <= 1, 'ROI manifest has duplicate rows for %s.', roi_stem);
if isscalar(match)
    source = manifest.source(match);
    labels = manifest.contributing_labels(match);
end
end

function assert_same_geometry(image_header, reference_header, ...
        image_file, reference_file)
assert(isequal(image_header.dim, reference_header.dim), ...
    'Dimension mismatch between %s and %s.', image_file, reference_file);
assert(max(abs(image_header.mat(:) - reference_header.mat(:))) < 1e-4, ...
    'Affine mismatch between %s and %s.', image_file, reference_file);
end

function structure = summarize_task_structure(metadata)
contexts = ["PERSON", "FOOD", "LOCATION", "CONTROL"];
runs = unique(metadata.run_id, 'stable');
run_id = runs;
session_id = zeros(numel(runs), 1);
counts = zeros(numel(runs), numel(contexts));
for run_idx = 1:numel(runs)
    run_mask = metadata.run_id == runs(run_idx);
    sessions = unique(metadata.session_id(run_mask));
    assert(isscalar(sessions), 'A run spans multiple acquisition sessions.');
    session_id(run_idx) = sessions;
    for context_idx = 1:numel(contexts)
        counts(run_idx, context_idx) = nnz( ...
            run_mask & metadata.context == contexts(context_idx));
    end
end
run_context_counts = table(run_id, session_id, counts(:, 1), counts(:, 2), ...
    counts(:, 3), counts(:, 4), 'VariableNames', ...
    {'run_id', 'session_id', 'PERSON', 'FOOD', 'LOCATION', 'CONTROL'});
structure = struct();
structure.run_context_counts = run_context_counts;
structure.contexts_are_run_blocked = all(sum(counts > 0, 2) == 1);
structure.runs_with_context = sum(counts > 0, 1);
structure.total_trials_by_context = sum(counts, 1);
structure.context_order = cellstr(contexts(:));
structure.n_sessions = numel(unique(session_id));
structure.runs_per_session = groupcounts(session_id);
structure.note = ['Runs contain uneven mixtures of contexts; splits are ' ...
    'balanced by session and explicitly checked for all 60 semantic ' ...
    'context-by-odor cells in both halves.'];
end
