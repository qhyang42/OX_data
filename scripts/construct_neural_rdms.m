%% Construct simple and crossnobis neural RDMs in five olfactory ROIs
% Subjects: subj_2--subj_6
% Conditions: PERSON, FOOD, LOCATION, CONTROL; odors 1:20 within context
% Outputs: RDMs/neural/subj_N_neural_RDMs.mat

script_path = mfilename('fullpath');
assert(~isempty(script_path), 'Run this file as a script, not as pasted code.');
project_root = fileparts(fileparts(script_path));
mri_root = fullfile(project_root, 'MRI');
output_root = fullfile(project_root, 'RDMs', 'neural');
utility_root = fullfile(project_root, 'utils', 'OX_utilities');
addpath(utility_root);

if exist('spm_vol', 'file') ~= 2 || exist('spm_read_vols', 'file') ~= 2
    spm_root = getenv('SPM_DIR');
    if isempty(spm_root)
        spm_root = '/Users/qhyang/Desktop/Utilities/spm';
    end
    assert(isfolder(spm_root), 'SPM directory not found: %s', spm_root);
    addpath(spm_root);
end
assert(exist('spm_vol', 'file') == 2 && exist('spm_read_vols', 'file') == 2, ...
    'SPM must be available to read the NIfTI masks.');
assert(exist('OX_load_trial_metadata', 'file') == 2, ...
    'OX_load_trial_metadata.m was not found under %s.', utility_root);
if ~isfolder(output_root)
    mkdir(output_root);
end

subject_ids = 2:6;
roi_names = ["AON", "PirF", "PirT", "olfAMG", "olfOFC"];
context_order = ["PERSON", "FOOD", "LOCATION", "CONTROL"];
odor_ids = (1:20)';
n_conditions = numel(context_order) * numel(odor_ids);
n_splits = 200;
n_candidate_splits = 100000;
split_seed = 1;
min_voxels = 10;

condition_index = (1:n_conditions)';
condition_context = repelem(context_order(:), numel(odor_ids));
condition_odor = repmat(odor_ids, numel(context_order), 1);
condition_label = condition_context + "_odor" + string(condition_odor);
condition_metadata = table(condition_index, condition_context, ...
    condition_odor, condition_label);

run_crossnobis_synthetic_checks();

for subject_id = subject_ids
    subject_tic = tic;
    subject_name = sprintf('subj_%d', subject_id);
    nifti_dir = fullfile(mri_root, subject_name, 'nifti');
    fit_file = fullfile(nifti_dir, 'sniff_single_trial_by_category_physio', ...
        'TYPED_FITHRF_GLMDENOISE_RR.mat');
    gm_mask_file = fullfile(nifti_dir, 'coreg', 'gm_mask_thr05_func.nii');
    functional_mask_file = fullfile(nifti_dir, ...
        'first_level_model_sniff_physio', 'spmT_0001_uncorrected_p001.nii');
    roi_dir = fullfile(nifti_dir, 'coreg', 'roi_decoding', 'primary');

    assert(isfile(fit_file), 'Missing GLMsingle file: %s', fit_file);
    assert(isfile(gm_mask_file), 'Missing gray-matter mask: %s', gm_mask_file);
    assert(isfile(functional_mask_file), ...
        'Missing functional restriction mask: %s', functional_mask_file);
    assert(isfolder(roi_dir), 'Missing primary ROI directory: %s', roi_dir);

    fprintf('\n[NEURAL RDM] Loading %s\n', subject_name);
    loaded = load(fit_file, 'modelmd');
    assert(isfield(loaded, 'modelmd'), 'GLMsingle file lacks modelmd: %s', fit_file);
    modelmd = squeeze(loaded.modelmd);
    clear loaded
    assert(ismatrix(modelmd), 'squeeze(modelmd) must be two-dimensional.');

    trial_metadata = OX_load_trial_metadata(subject_id);
    n_trials = height(trial_metadata);
    assert(n_trials == 800, 'Expected 800 trials for %s; found %d.', ...
        subject_name, n_trials);
    assert(size(modelmd, 2) == n_trials, ...
        'modelmd has %d trials but metadata has %d for %s.', ...
        size(modelmd, 2), n_trials, subject_name);

    [known_context, trial_context_id] = ismember( ...
        trial_metadata.context, context_order);
    [known_odor, trial_odor_id] = ismember(trial_metadata.odor, odor_ids);
    assert(all(known_context) && all(known_odor), ...
        'Unknown context or odor label found for %s.', subject_name);
    trial_condition_id = (trial_context_id - 1) * numel(odor_ids) + trial_odor_id;
    assert(all(trial_condition_id >= 1 & trial_condition_id <= n_conditions), ...
        'Condition encoding failed for %s.', subject_name);
    total_condition_counts = accumarray(trial_condition_id, 1, ...
        [n_conditions, 1]);
    assert(all(total_condition_counts > 0), ...
        'At least one condition has no trials for %s.', subject_name);

    split_seed_subject = split_seed + subject_id * 1000;
    split_metadata = select_balanced_run_splits( ...
        trial_metadata.run_id, trial_metadata.session_id, ...
        trial_condition_id, n_conditions, n_splits, ...
        n_candidate_splits, split_seed_subject);
    fprintf(['  Selected %d unique session-balanced run splits | ' ...
        'condition counts per half %d--%d\n'], n_splits, ...
        min(split_metadata.condition_counts, [], 'all'), ...
        max(split_metadata.condition_counts, [], 'all'));

    gm_header = spm_vol(gm_mask_file);
    gm_mask = spm_read_vols(gm_header) > 0;
    functional_header = spm_vol(functional_mask_file);
    assert_same_geometry(functional_header, gm_header, ...
        functional_mask_file, gm_mask_file);
    functional_mask = spm_read_vols(functional_header) > 0;
    restriction_mask = gm_mask & functional_mask;
    gm_indices = find(gm_mask);
    assert(size(modelmd, 1) == numel(gm_indices), ...
        ['GLMsingle/gray-matter mapping failed: modelmd has %d rows but ' ...
         'find(gm_mask) has %d voxels.'], size(modelmd, 1), numel(gm_indices));
    model_index_volume = zeros(size(gm_mask), 'uint32');
    model_index_volume(gm_indices) = uint32(1:numel(gm_indices));

    % This basis removes both run intercepts and condition means from the
    % trial patterns. Its rows provide orthonormal residual contrasts for
    % estimating voxel-by-voxel noise covariance.
    noise_design = make_noise_design(trial_metadata.run_id, ...
        trial_condition_id, n_conditions);
    noise_basis = null(noise_design', 1e-10);
    noise_degrees_of_freedom = size(noise_basis, 2);
    expected_df = n_trials - rank(noise_design);
    assert(noise_degrees_of_freedom == expected_df && expected_df > 1, ...
        'Noise residual-basis rank check failed for %s.', subject_name);

    n_rois = numel(roi_names);
    roi_results = repmat(struct(), n_rois, 1);
    roi_mask_files = strings(n_rois, 1);
    roi_mask_voxels = zeros(n_rois, 1);
    roi_gm_voxels = zeros(n_rois, 1);
    roi_restricted_voxels = zeros(n_rois, 1);
    roi_usable_voxels = zeros(n_rois, 1);

    for roi_index = 1:n_rois
        roi_name = roi_names(roi_index);
        roi_file = char(fullfile(roi_dir, ...
            roi_name + "_bilateral_func_thr02.nii"));
        assert(isfile(roi_file), 'Missing ROI mask: %s', roi_file);
        roi_header = spm_vol(roi_file);
        assert_same_geometry(roi_header, gm_header, roi_file, gm_mask_file);
        roi_mask = spm_read_vols(roi_header) > 0;
        restricted_roi_mask = roi_mask & restriction_mask;
        feature_indices = double(model_index_volume(restricted_roi_mask));
        feature_indices = feature_indices(feature_indices > 0);
        feature_indices = unique(feature_indices(:), 'stable');
        if ~isempty(feature_indices)
            feature_indices = feature_indices( ...
                all(isfinite(modelmd(feature_indices, :)), 2));
        end
        n_voxels = numel(feature_indices);
        assert(n_voxels >= min_voxels, ...
            '%s %s has %d usable voxels; at least %d are required.', ...
            subject_name, roi_name, n_voxels, min_voxels);

        roi_mask_files(roi_index) = string(roi_file);
        roi_mask_voxels(roi_index) = nnz(roi_mask);
        roi_gm_voxels(roi_index) = nnz(roi_mask & gm_mask);
        roi_restricted_voxels(roi_index) = nnz(restricted_roi_mask);
        roi_usable_voxels(roi_index) = n_voxels;
        fprintf('  ROI %s: %d usable voxels\n', roi_name, n_voxels);

        X_raw = double(modelmd(feature_indices, :)');
        X = center_voxels_within_runs(X_raw, trial_metadata.run_id);
        condition_patterns = average_condition_patterns( ...
            X, trial_condition_id, n_conditions);

        pattern_correlations = corr(condition_patterns');
        assert(all(isfinite(pattern_correlations), 'all'), ...
            '%s %s has an undefined condition-pattern correlation.', ...
            subject_name, roi_name);
        D_neural_simple = 1 - pattern_correlations;
        D_neural_simple = (D_neural_simple + D_neural_simple') / 2;
        D_neural_simple(1:n_conditions+1:end) = 0;
        assert(min(D_neural_simple, [], 'all') >= -1e-10 && ...
            max(D_neural_simple, [], 'all') <= 2 + 1e-10, ...
            'Simple correlation distance is outside [0,2] for %s %s.', ...
            subject_name, roi_name);
        D_neural_simple = min(max(D_neural_simple, 0), 2);

        noise_samples = noise_basis' * X_raw;
        [noise_covariance, shrinkage, sample_covariance] = ...
            schafer_strimmer_covariance(noise_samples);
        [whitening_matrix, whitening_error, covariance_eigenvalues] = ...
            symmetric_whitener(noise_covariance);
        assert(whitening_error < 1e-8, ...
            'Whitening check failed for %s %s (max error %.3g).', ...
            subject_name, roi_name, whitening_error);

        X_whitened = X * whitening_matrix;
        [split_patterns_a, split_patterns_b] = make_split_patterns( ...
            X_whitened, trial_metadata.run_id, trial_condition_id, ...
            split_metadata.run_membership, split_metadata.run_ids, ...
            split_metadata.condition_counts);

        D_neural_crossnobis_splits = zeros( ...
            n_conditions, n_conditions, n_splits, 'single');
        crossnobis_sum = zeros(n_conditions, n_conditions);
        for split_index = 1:n_splits
            rows = (split_index - 1) * n_conditions + (1:n_conditions);
            split_D = crossnobis_from_patterns( ...
                split_patterns_a(rows, :), split_patterns_b(rows, :));
            validate_rdm(split_D, n_conditions, false, ...
                sprintf('%s %s crossnobis split %d', ...
                subject_name, roi_name, split_index));
            D_neural_crossnobis_splits(:, :, split_index) = single(split_D);
            crossnobis_sum = crossnobis_sum + split_D;
        end
        D_neural_crossnobis = crossnobis_sum / n_splits;
        D_neural_crossnobis = (D_neural_crossnobis + ...
            D_neural_crossnobis') / 2;
        D_neural_crossnobis(1:n_conditions+1:end) = 0;
        D_neural_crossnobis_sd = std( ...
            double(D_neural_crossnobis_splits), 0, 3);

        validate_rdm(D_neural_simple, n_conditions, true, ...
            sprintf('%s %s simple', subject_name, roi_name));
        validate_rdm(D_neural_crossnobis, n_conditions, false, ...
            sprintf('%s %s mean crossnobis', subject_name, roi_name));
        validate_rdm(D_neural_crossnobis_sd, n_conditions, false, ...
            sprintf('%s %s crossnobis SD', subject_name, roi_name));

        roi_results(roi_index).roi_name = char(roi_name);
        roi_results(roi_index).mask_file = char(roi_file);
        roi_results(roi_index).n_voxels = n_voxels;
        roi_results(roi_index).model_feature_indices = feature_indices;
        roi_results(roi_index).condition_patterns = condition_patterns;
        roi_results(roi_index).D_neural_simple = D_neural_simple;
        roi_results(roi_index).D_neural_crossnobis = D_neural_crossnobis;
        roi_results(roi_index).D_neural_crossnobis_splits = ...
            D_neural_crossnobis_splits;
        roi_results(roi_index).D_neural_crossnobis_sd = ...
            D_neural_crossnobis_sd;
        roi_results(roi_index).noise = struct( ...
            'method', 'Schaefer-Strimmer off-diagonal shrinkage to diagonal', ...
            'degrees_of_freedom', noise_degrees_of_freedom, ...
            'shrinkage', shrinkage, ...
            'sample_covariance', sample_covariance, ...
            'regularized_covariance', noise_covariance, ...
            'whitening_matrix', whitening_matrix, ...
            'regularized_covariance_eigenvalues', covariance_eigenvalues, ...
            'max_whitening_identity_error', whitening_error);
    end

    roi_metadata = table(roi_names(:), roi_mask_files, roi_mask_voxels, ...
        roi_gm_voxels, roi_restricted_voxels, roi_usable_voxels, ...
        'VariableNames', {'roi_name', 'mask_file', 'n_anatomical_voxels', ...
        'n_gray_matter_voxels', 'n_functionally_restricted_voxels', ...
        'n_usable_voxels'});

    results = struct();
    results.subject_id = subject_id;
    results.subject_name = subject_name;
    results.condition_metadata = condition_metadata;
    results.trial_metadata = trial_metadata;
    results.total_condition_counts = total_condition_counts;
    results.split_metadata = split_metadata;
    results.roi_metadata = roi_metadata;
    results.roi_results = roi_results;
    results.preprocessing = struct( ...
        'beta_source', fit_file, ...
        'within_run_voxel_centering', true, ...
        'condition_order', 'PERSON, FOOD, LOCATION, CONTROL; odors 1:20', ...
        'roi_selection', 'primary bilateral', ...
        'gray_matter_mask', gm_mask_file, ...
        'functional_mask', functional_mask_file, ...
        'functional_threshold', 'Odor > Rest, uncorrected p < .001', ...
        'simple_distance', '1 - Pearson correlation across voxels', ...
        'crossnobis_scaling', 'crossvalidated whitened inner product / n_voxels', ...
        'negative_crossnobis_retained', true, ...
        'noise_fixed_effects_removed', 'exact acquisition run and 80-condition means', ...
        'noise_covariance', 'Schaefer-Strimmer off-diagonal shrinkage to diagonal');
    results.options = struct( ...
        'subject_ids', subject_ids, 'roi_names', roi_names, ...
        'context_order', context_order, 'odor_ids', odor_ids, ...
        'n_splits', n_splits, 'n_candidate_splits', n_candidate_splits, ...
        'split_seed', split_seed_subject, 'min_voxels', min_voxels);
    results.dimension_order = struct( ...
        'condition_patterns', {{'condition', 'voxel'}}, ...
        'crossnobis_splits', {{'condition_i', 'condition_j', 'split'}}, ...
        'split_condition_counts', {{'split', 'half', 'condition'}}, ...
        'split_run_membership', {{'split', 'run'}});
    results.runtime_minutes = toc(subject_tic) / 60;
    results.matlab_version = version;

    output_file = fullfile(output_root, sprintf( ...
        '%s_neural_RDMs.mat', subject_name));
    save(output_file, 'results', '-v7.3');
    fprintf('Saved %s (%.1f minutes)\n', output_file, results.runtime_minutes);
    clear modelmd noise_basis roi_results results
end

fprintf('\nFinished neural RDM construction for %d subjects in %s\n', ...
    numel(subject_ids), output_root);

function X_centered = center_voxels_within_runs(X, run_ids)
X_centered = X;
runs = unique(run_ids, 'stable');
for run_index = 1:numel(runs)
    in_run = run_ids == runs(run_index);
    X_centered(in_run, :) = X_centered(in_run, :) - ...
        mean(X_centered(in_run, :), 1);
end
end

function patterns = average_condition_patterns(X, condition_ids, n_conditions)
patterns = zeros(n_conditions, size(X, 2));
for condition_index = 1:n_conditions
    in_condition = condition_ids == condition_index;
    assert(any(in_condition), 'Condition %d has no observations.', condition_index);
    patterns(condition_index, :) = mean(X(in_condition, :), 1);
end
end

function design = make_noise_design(run_ids, condition_ids, n_conditions)
runs = unique(run_ids, 'stable');
[known_runs, run_index] = ismember(run_ids, runs);
assert(all(known_runs), 'Unable to encode one or more run IDs.');
n_trials = numel(run_ids);
run_design = sparse((1:n_trials)', run_index, 1, n_trials, numel(runs));
condition_design = sparse((1:n_trials)', condition_ids, 1, ...
    n_trials, n_conditions);
% Keep the complete column space. The run/condition incidence structure can
% create more than one exact dependency, so downstream rank and null-space
% calculations must not assume that dropping a single column is sufficient.
design = full([run_design, condition_design]);
end

function split = select_balanced_run_splits(run_ids, session_ids, ...
        condition_ids, n_conditions, n_splits, n_candidates, random_seed)
runs = unique(run_ids, 'stable');
n_runs = numel(runs);
assert(mod(n_runs, 2) == 0, 'An even number of runs is required.');
runs_per_half = n_runs / 2;
[known_runs, trial_run_index] = ismember(run_ids, runs);
assert(all(known_runs), 'Unable to encode one or more run IDs.');

run_sessions = zeros(n_runs, 1);
for run_index = 1:n_runs
    values = unique(session_ids(trial_run_index == run_index));
    assert(isscalar(values), 'Run %g spans multiple sessions.', runs(run_index));
    run_sessions(run_index) = values;
end
session_values = unique(run_sessions, 'stable');
session_sizes = arrayfun(@(x) nnz(run_sessions == x), session_values);
base_take = floor(session_sizes / 2);
n_extra = runs_per_half - sum(base_take);
odd_sessions = find(mod(session_sizes, 2) == 1);
assert(n_extra >= 0 && n_extra <= numel(odd_sessions), ...
    'Cannot construct a session-balanced half of %d runs.', runs_per_half);

run_condition_counts = zeros(n_runs, n_conditions);
for trial_index = 1:numel(run_ids)
    run_condition_counts(trial_run_index(trial_index), ...
        condition_ids(trial_index)) = ...
        run_condition_counts(trial_run_index(trial_index), ...
        condition_ids(trial_index)) + 1;
end
total_counts = sum(run_condition_counts, 1);

previous_rng = rng;
rng_cleanup = onCleanup(@() rng(previous_rng));
rng(random_seed, 'twister');
candidate_membership = false(n_candidates, n_runs);
for candidate_index = 1:n_candidates
    take = base_take;
    if n_extra > 0
        chosen_odd = odd_sessions(randperm(numel(odd_sessions), n_extra));
        take(chosen_odd) = take(chosen_odd) + 1;
    end
    membership = false(1, n_runs);
    for session_index = 1:numel(session_values)
        indices = find(run_sessions == session_values(session_index));
        selected = indices(randperm(numel(indices), take(session_index)));
        membership(selected) = true;
    end
    assert(nnz(membership) == runs_per_half, ...
        'Internal error constructing a run half.');
    if ~membership(1)
        membership = ~membership;
    end
    candidate_membership(candidate_index, :) = membership;
end

candidate_membership = unique(candidate_membership, 'rows', 'stable');
counts_a = double(candidate_membership) * run_condition_counts;
counts_b = total_counts - counts_a;
has_coverage = all(counts_a > 0 & counts_b > 0, 2);
candidate_membership = candidate_membership(has_coverage, :);
counts_a = counts_a(has_coverage, :);
counts_b = counts_b(has_coverage, :);
assert(size(candidate_membership, 1) >= n_splits, ...
    ['Only %d unique candidate partitions represented all conditions in ' ...
     'both halves; %d are required.'], size(candidate_membership, 1), n_splits);

deviation = abs(counts_a - total_counts / 2);
max_deviation = max(deviation, [], 2);
sum_squared_deviation = sum(deviation .^ 2, 2);
candidate_order = (1:size(candidate_membership, 1))';
[~, ranking] = sortrows( ...
    [max_deviation, sum_squared_deviation, candidate_order], [1, 2, 3]);
selected = ranking(1:n_splits);
membership = candidate_membership(selected, :);
selected_counts_a = counts_a(selected, :);
selected_counts_b = counts_b(selected, :);

session_counts_a = zeros(n_splits, numel(session_values));
session_counts_b = zeros(n_splits, numel(session_values));
for session_index = 1:numel(session_values)
    in_session = run_sessions == session_values(session_index);
    session_counts_a(:, session_index) = sum(membership(:, in_session), 2);
    session_counts_b(:, session_index) = sum(~membership(:, in_session), 2);
end
assert(all(sum(membership, 2) == runs_per_half), ...
    'A retained split does not have %d A runs.', runs_per_half);
assert(all(sum(~membership, 2) == runs_per_half), ...
    'A retained split does not have %d B runs.', runs_per_half);
assert(all(abs(session_counts_a - session_counts_b) <= 1, 'all'), ...
    'A retained split is not balanced within session.');
assert(all(selected_counts_a > 0 & selected_counts_b > 0, 'all'), ...
    'A retained split lacks condition coverage.');
assert(size(unique(membership, 'rows'), 1) == n_splits, ...
    'Retained run partitions are not unique.');

half_a_run_ids = zeros(n_splits, runs_per_half);
half_b_run_ids = zeros(n_splits, runs_per_half);
for split_index = 1:n_splits
    half_a_run_ids(split_index, :) = runs(membership(split_index, :));
    half_b_run_ids(split_index, :) = runs(~membership(split_index, :));
end
condition_counts = zeros(n_splits, 2, n_conditions, 'uint16');
condition_counts(:, 1, :) = uint16(reshape(selected_counts_a, ...
    n_splits, 1, n_conditions));
condition_counts(:, 2, :) = uint16(reshape(selected_counts_b, ...
    n_splits, 1, n_conditions));

split = struct();
split.n_requested_candidates = n_candidates;
split.n_unique_candidates = size(candidate_membership, 1);
split.random_seed = random_seed;
split.run_ids = runs;
split.run_session_ids = run_sessions;
split.session_ids = session_values;
split.runs_per_half = runs_per_half;
split.run_membership = membership;
split.half_a_run_ids = half_a_run_ids;
split.half_b_run_ids = half_b_run_ids;
split.condition_counts = condition_counts;
split.session_run_counts_a = session_counts_a;
split.session_run_counts_b = session_counts_b;
split.max_condition_deviation_from_half = max_deviation(selected);
split.sum_squared_condition_deviation = sum_squared_deviation(selected);
split.total_condition_counts = total_counts(:);
split.selection_rule = ['rank unique, session-balanced, full-coverage ' ...
    'partitions by maximum then summed-squared condition-count deviation'];
end

function [covariance_shrunk, lambda, covariance_sample] = ...
        schafer_strimmer_covariance(noise_samples)
% Analytic shrinkage of off-diagonal correlations toward zero, retaining
% each voxel's sample variance (Schaefer and Strimmer, 2005).
noise_samples = double(noise_samples);
[n_samples, n_voxels] = size(noise_samples);
assert(n_samples >= 2 && n_voxels >= 2, ...
    'Noise covariance needs at least two samples and two voxels.');
assert(all(isfinite(noise_samples), 'all'), ...
    'Noise samples contain non-finite values.');

covariance_sample = (noise_samples' * noise_samples) / n_samples;
variances = diag(covariance_sample);
assert(all(isfinite(variances) & variances > 0), ...
    'At least one voxel has zero or invalid residual variance.');
standardized = noise_samples ./ sqrt(variances)';
sample_correlation = (standardized' * standardized) / n_samples;
sample_correlation = (sample_correlation + sample_correlation') / 2;
sample_correlation(1:n_voxels+1:end) = 1;

% The variance of each sample-correlation entry is estimated from the
% variance of the voxel-product observations. Matrix multiplication avoids
% constructing one voxel-by-voxel product matrix per residual sample.
sum_squared_products = (standardized .^ 2)' * (standardized .^ 2);
correlation_variance = max(sum_squared_products - ...
    n_samples * sample_correlation .^ 2, 0) / ...
    (n_samples * max(n_samples - 1, 1));
upper_triangle = triu(true(n_voxels), 1);
denominator = sum(sample_correlation(upper_triangle) .^ 2);
if denominator <= eps
    lambda = 1;
else
    lambda = min(1, max(0, ...
        sum(correlation_variance(upper_triangle)) / denominator));
end

correlation_shrunk = (1 - lambda) * sample_correlation;
correlation_shrunk(1:n_voxels+1:end) = 1;
scale = sqrt(variances);
covariance_shrunk = correlation_shrunk .* (scale * scale');
covariance_shrunk = (covariance_shrunk + covariance_shrunk') / 2;
end

function [whitener, max_error, eigenvalues] = symmetric_whitener(covariance)
[vectors, values] = eig((covariance + covariance') / 2, 'vector');
eigenvalues = real(values);
assert(all(isfinite(eigenvalues)) && min(eigenvalues) > 0, ...
    'Regularized covariance is not positive definite.');
whitener = vectors * diag(1 ./ sqrt(eigenvalues)) * vectors';
whitener = real((whitener + whitener') / 2);
identity_check = whitener' * covariance * whitener;
max_error = max(abs(identity_check - eye(size(covariance))), [], 'all');
end

function [patterns_a, patterns_b] = make_split_patterns(X, run_ids, ...
        condition_ids, run_membership, runs, condition_counts)
n_splits = size(run_membership, 1);
n_conditions = size(condition_counts, 3);
n_trials = size(X, 1);
[known_runs, trial_run_index] = ismember(run_ids, runs);
assert(all(known_runs), 'Unable to map trials to split run IDs.');

max_entries = n_splits * n_trials;
row_a = zeros(max_entries, 1);
column_a = zeros(max_entries, 1);
value_a = zeros(max_entries, 1);
row_b = zeros(max_entries, 1);
column_b = zeros(max_entries, 1);
value_b = zeros(max_entries, 1);
cursor_a = 0;
cursor_b = 0;

for split_index = 1:n_splits
    in_a = run_membership(split_index, trial_run_index)';
    counts_a = double(reshape(condition_counts( ...
        split_index, 1, :), n_conditions, 1));
    counts_b = double(reshape(condition_counts( ...
        split_index, 2, :), n_conditions, 1));
    assert(all(counts_a > 0) && all(counts_b > 0), ...
        'Split %d has an empty condition.', split_index);

    indices_a = find(in_a);
    rows_a = (split_index - 1) * n_conditions + ...
        condition_ids(indices_a);
    destination_a = cursor_a + (1:numel(indices_a));
    row_a(destination_a) = rows_a;
    column_a(destination_a) = indices_a;
    value_a(destination_a) = 1 ./ counts_a(condition_ids(indices_a));
    cursor_a = cursor_a + numel(indices_a);

    indices_b = find(~in_a);
    rows_b = (split_index - 1) * n_conditions + ...
        condition_ids(indices_b);
    destination_b = cursor_b + (1:numel(indices_b));
    row_b(destination_b) = rows_b;
    column_b(destination_b) = indices_b;
    value_b(destination_b) = 1 ./ counts_b(condition_ids(indices_b));
    cursor_b = cursor_b + numel(indices_b);
end

weights_a = sparse(row_a(1:cursor_a), column_a(1:cursor_a), ...
    value_a(1:cursor_a), n_splits * n_conditions, n_trials);
weights_b = sparse(row_b(1:cursor_b), column_b(1:cursor_b), ...
    value_b(1:cursor_b), n_splits * n_conditions, n_trials);
patterns_a = weights_a * X;
patterns_b = weights_b * X;
assert(all(isfinite(patterns_a), 'all') && all(isfinite(patterns_b), 'all'), ...
    'A split pattern contains non-finite values.');
end

function distance = crossnobis_from_patterns(patterns_a, patterns_b)
assert(isequal(size(patterns_a), size(patterns_b)), ...
    'Crossnobis A and B pattern matrices must have equal size.');
n_voxels = size(patterns_a, 2);
cross_gram = (patterns_a * patterns_b') / n_voxels;
self_cross = diag(cross_gram);
distance = self_cross + self_cross' - cross_gram - cross_gram';
distance = (distance + distance') / 2;
distance(1:size(distance, 1)+1:end) = 0;
end

function validate_rdm(distance, n_conditions, bounded, label)
assert(isequal(size(distance), [n_conditions, n_conditions]), ...
    '%s RDM is not %d-by-%d.', label, n_conditions, n_conditions);
assert(all(isfinite(distance), 'all'), '%s RDM contains non-finite values.', label);
assert(max(abs(distance - distance'), [], 'all') < 1e-10, ...
    '%s RDM is not symmetric.', label);
assert(all(diag(distance) == 0), '%s RDM diagonal is not exactly zero.', label);
if bounded
    assert(min(distance, [], 'all') >= 0 && ...
        max(distance, [], 'all') <= 2, ...
        '%s correlation-distance RDM is outside [0,2].', label);
end
end

function assert_same_geometry(image_header, reference_header, ...
        image_file, reference_file)
assert(isequal(image_header.dim, reference_header.dim), ...
    'Dimension mismatch between %s and %s.', image_file, reference_file);
assert(max(abs(image_header.mat(:) - reference_header.mat(:))) < 1e-4, ...
    'Affine mismatch between %s and %s.', image_file, reference_file);
end

function run_crossnobis_synthetic_checks()
patterns_a = [0, 0; 1, 0; 0, 2];
patterns_b = patterns_a;
distance = crossnobis_from_patterns(patterns_a, patterns_b);
expected = [0, 0.5, 2; 0.5, 0, 2.5; 2, 2.5, 0];
assert(max(abs(distance - expected), [], 'all') < 1e-12, ...
    'Synthetic crossnobis scaling check failed.');

negative_a = [0, 0; 1, 0];
negative_b = [0, 0; -1, 0];
negative_distance = crossnobis_from_patterns(negative_a, negative_b);
assert(abs(negative_distance(1, 2) + 0.5) < 1e-12, ...
    'Synthetic crossnobis negative-distance check failed.');
assert(isequal(distance, distance') && all(diag(distance) == 0), ...
    'Synthetic crossnobis symmetry or diagonal check failed.');
fprintf('Crossnobis synthetic checks passed.\n');
end
