function results = OX_roi_odor_context_transfer(subjidx, varargin)
%OX_ROI_ODOR_CONTEXT_TRANSFER Decode odor identity across semantic contexts.
%
%   results = OX_roi_odor_context_transfer(subjidx, Name, Value, ...)
%
% The confirmatory analysis uses the bilateral ROI masks and a 3-by-3
% PERSON/FOOD/LOCATION source-to-target matrix. CONTROL-derived odor
% templates provide a separate transfer measurement into each semantic
% context. All decoding uses leave-one-run-out folds and Pearson template
% similarity across voxels.
% ROISelection may be 'old' (default), 'primary', 'secondary', or 'all'.

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'subjidx', @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x == round(x));
addParameter(p, 'MRIRoot', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'CueRoot', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'ROISelection', 'old', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'ROIDir', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'MinVoxels', 10, @(x) isnumeric(x) && isscalar(x) && x >= 2 && x == round(x));
addParameter(p, 'MaxROIs', Inf, @(x) isnumeric(x) && isscalar(x) && ...
    (isinf(x) || (x >= 1 && x == round(x))));
addParameter(p, 'NumPermutations', 1000, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x == round(x));
addParameter(p, 'RandomSeed', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'UseParallel', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'OutputDir', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'SaveOutputs', true, @(x) islogical(x) && isscalar(x));
parse(p, subjidx, varargin{:});
opts = p.Results;
opts.ROISelection = OX_normalize_decoding_roi_selection(opts.ROISelection);

subject_names = {'240711_fMRI_OX_NWU_AS', ...
                 '240723_fMRI_OX_NWU_LS', ...
                 '240814_fMRI_OX_NWU_JN', ...
                 '240816_fMRI_OX_NWU_RR', ...
                 '241018_fMRI_OX_NWU_BN', ...
                 '250117_fMRI_OX_NWU_VS'};
assert(subjidx >= 1 && subjidx <= numel(subject_names), ...
    'subjidx must be between 1 and %d.', numel(subject_names));

analysis_tic = tic;
subject_name = sprintf('subj_%d', subjidx);
mri_root = resolve_mri_root(opts.MRIRoot);
nifti_dir = fullfile(mri_root, subject_name, 'nifti');
beta_dir = fullfile(nifti_dir, 'single_trial_by_category');
fit_file = fullfile(beta_dir, 'TYPED_FITHRF_GLMDENOISE_RR.mat');
gm_mask_file = fullfile(nifti_dir, 'coreg', 'gm_mask_thr05_func.nii');
[roi_selection, roi_dirs] = OX_resolve_decoding_roi_selection( ...
    nifti_dir, opts.ROISelection, opts.ROIDir);
assert(isfile(fit_file), 'Missing odor-aligned GLMsingle file: %s', fit_file);
assert(isfile(gm_mask_file), 'Missing gray-matter mask: %s', gm_mask_file);
assert(exist('spm_vol', 'file') == 2 && exist('spm_read_vols', 'file') == 2, ...
    'SPM must be on the MATLAB path. Run setup_ox first.');

fprintf('\n[ODOR CONTEXT TRANSFER] Loading %s\n', fit_file);
loaded = load(fit_file, 'modelmd');
assert(isfield(loaded, 'modelmd'), 'GLMsingle file lacks modelmd: %s', fit_file);
modelmd = squeeze(loaded.modelmd);
clear loaded
assert(ismatrix(modelmd), 'squeeze(modelmd) must produce a two-dimensional matrix.');

metadata_args = {};
if strlength(string(opts.CueRoot)) > 0
    metadata_args = {'CueRoot', char(string(opts.CueRoot))};
end
trial_metadata = OX_load_trial_metadata(subjidx, metadata_args{:});
n_trials = height(trial_metadata);
assert(size(modelmd, 2) == n_trials, ...
    'modelmd has %d trials but metadata has %d.', size(modelmd, 2), n_trials);
if strlength(string(opts.CueRoot)) == 0 && exist('OX_get_odor', 'file') == 2
    [legacy_odors, legacy_contexts] = OX_get_odor(subject_name);
    assert(isequal(double(legacy_odors(:)), trial_metadata.odor) && ...
        isequal(upper(strtrim(string(legacy_contexts(:)))), trial_metadata.context), ...
        'Exact trial metadata order differs from the labels used by OX_get_odor.');
end
odor_values = unique(trial_metadata.odor, 'sorted');
n_odors = numel(odor_values);

modelmd = center_voxels_within_runs(modelmd, trial_metadata.run_id);
fprintf('Subject %s | %d trials | %d runs | %d odors\n', ...
    subject_name, n_trials, numel(unique(trial_metadata.run_id)), n_odors);
fprintf('Preprocessing: within-run voxel centering; spatial pattern demeaning; no voxel scaling.\n');
fprintf('ROI selection: %s | Directories: %s\n', ...
    roi_selection, strjoin(string(roi_dirs), ', '));

gm_header = spm_vol(gm_mask_file);
gm_mask = spm_read_vols(gm_header) > 0;
gm_indices = find(gm_mask);
assert(size(modelmd, 1) == numel(gm_indices), ...
    'modelmd rows (%d) do not match gray-matter voxels (%d).', ...
    size(modelmd, 1), numel(gm_indices));

[roi_files, roi_stems, manifest] = OX_discover_decoding_roi_files( ...
    roi_dirs, true, opts.MaxROIs);
n_rois = numel(roi_files);
model_index_volume = zeros(size(gm_mask), 'uint32');
model_index_volume(gm_indices) = uint32(1:numel(gm_indices));

roi_names = strings(n_rois, 1);
roi_sources = strings(n_rois, 1);
roi_labels = strings(n_rois, 1);
roi_features = cell(n_rois, 1);
n_mask_voxels = zeros(n_rois, 1);
n_gm_overlap = zeros(n_rois, 1);
n_features = zeros(n_rois, 1);
status = repmat("ok", n_rois, 1);

for roi_idx = 1:n_rois
    roi_header = spm_vol(roi_files{roi_idx});
    assert_same_geometry(roi_header, gm_header, roi_files{roi_idx}, gm_mask_file);
    roi_mask = spm_read_vols(roi_header) > 0;
    n_mask_voxels(roi_idx) = nnz(roi_mask);
    overlap = find(roi_mask & gm_mask);
    n_gm_overlap(roi_idx) = numel(overlap);
    feature_indices = double(model_index_volume(overlap));
    feature_indices = feature_indices(feature_indices > 0);
    if ~isempty(feature_indices)
        feature_indices = feature_indices(all(isfinite(modelmd(feature_indices, :)), 2));
    end
    roi_features{roi_idx} = feature_indices(:)';
    n_features(roi_idx) = numel(feature_indices);
    roi_names(roi_idx) = erase(string(roi_stems{roi_idx}), "_func_thr02");
    [roi_sources(roi_idx), roi_labels(roi_idx)] = ...
        manifest_metadata(manifest, roi_stems{roi_idx});
    if n_gm_overlap(roi_idx) < opts.MinVoxels
        status(roi_idx) = "insufficient_gm_overlap";
    elseif n_features(roi_idx) < opts.MinVoxels
        status(roi_idx) = "insufficient_features";
    end
end
valid_roi = status == "ok";
assert(any(valid_roi), 'No bilateral ROI has at least %d usable voxels.', opts.MinVoxels);
fprintf('Bilateral ROIs: %d total, %d valid.\n', n_rois, nnz(valid_roi));

use_parallel = opts.UseParallel && license('test', 'Distrib_Computing_Toolbox');
if opts.UseParallel && ~use_parallel
    warning('Parallel Computing Toolbox unavailable; running serially.');
end
if use_parallel && isempty(gcp('nocreate'))
    parpool;
end

observed_cells = cell(n_rois, 1);
fprintf('Running observed leave-one-run-out transfer analysis...\n');
observed_tic = tic;
odor_labels = trial_metadata.odor;
context_labels = trial_metadata.context;
run_labels = trial_metadata.run_id;
minimum_voxels = opts.MinVoxels;
if use_parallel
    model_constant = parallel.pool.Constant(modelmd);
    model_cleanup = onCleanup(@() delete(model_constant));
    parfor roi_idx = 1:n_rois
        if valid_roi(roi_idx)
            X = double(model_constant.Value(roi_features{roi_idx}, :)');
            observed_cells{roi_idx} = OX_template_context_transfer(X, ...
                odor_labels, context_labels, run_labels, ...
                'MinVoxels', minimum_voxels, 'DemeanPatterns', true, ...
                'CollectDetails', true);
        end
    end
else
    for roi_idx = 1:n_rois
        if valid_roi(roi_idx)
            X = double(modelmd(roi_features{roi_idx}, :)');
            observed_cells{roi_idx} = OX_template_context_transfer(X, ...
                odor_labels, context_labels, run_labels, ...
                'MinVoxels', minimum_voxels, 'DemeanPatterns', true, ...
                'CollectDetails', true);
        end
    end
end
fprintf('Observed analysis completed in %.2f minutes.\n', toc(observed_tic) / 60);

observed = package_observed(observed_cells, n_rois, n_trials, ...
    numel(unique(trial_metadata.run_id)));
failed = valid_roi & ~isfinite(observed.primary_evidence_effect);
status(failed) = "decoder_failed";
valid_results = isfinite(observed.primary_evidence_effect);
assert(any(valid_results), 'No ROI produced a finite context-transfer result.');

null = initialize_null(n_rois, opts.NumPermutations);
if opts.NumPermutations > 0
    fprintf('Running %d context-by-run odor-label permutations...\n', opts.NumPermutations);
    rng(opts.RandomSeed, 'twister');
    permutation_tic = tic;
    for permutation_idx = 1:opts.NumPermutations
        permuted_odors = permute_odors_within_context_run( ...
            trial_metadata.odor, trial_metadata.context, trial_metadata.run_id);
        permutation_cells = cell(n_rois, 1);
        if use_parallel
            parfor roi_idx = 1:n_rois
                if valid_results(roi_idx)
                    X = double(model_constant.Value(roi_features{roi_idx}, :)');
                    permutation_cells{roi_idx} = OX_template_context_transfer(X, ...
                        permuted_odors, context_labels, run_labels, ...
                        'MinVoxels', minimum_voxels, 'DemeanPatterns', true, ...
                        'CollectDetails', false);
                end
            end
        else
            for roi_idx = 1:n_rois
                if valid_results(roi_idx)
                    X = double(modelmd(roi_features{roi_idx}, :)');
                    permutation_cells{roi_idx} = OX_template_context_transfer(X, ...
                        permuted_odors, context_labels, run_labels, ...
                        'MinVoxels', minimum_voxels, 'DemeanPatterns', true, ...
                        'CollectDetails', false);
                end
            end
        end
        null = store_permutation(null, permutation_cells, permutation_idx, n_rois);
        if permutation_idx == 1 || mod(permutation_idx, 10) == 0 || ...
                permutation_idx == opts.NumPermutations
            fprintf('  permutation %d/%d | %.1f min\n', permutation_idx, ...
                opts.NumPermutations, toc(permutation_tic) / 60);
        end
    end
end
null = finish_null_inference(null, observed, valid_results);

hemisphere = repmat("bilateral", n_rois, 1);
roi_file = string(roi_files(:));
roi_metadata = table(roi_names, roi_sources, hemisphere, roi_labels, roi_file, ...
    status, n_mask_voxels, n_gm_overlap, n_features, ...
    'VariableNames', {'roi_name', 'source', 'hemisphere', 'contributing_labels', ...
    'mask_file', 'status', 'n_mask_voxels', 'n_gm_overlap', 'n_features_used'});
summary = [roi_metadata, table( ...
    observed.within_evidence, observed.cross_evidence, observed.primary_evidence_effect, ...
    observed.within_accuracy, observed.cross_accuracy, observed.primary_accuracy_effect, ...
    observed.control_transfer_mean_evidence, observed.control_transfer_mean_accuracy, ...
    observed.control_within_evidence, observed.control_within_accuracy, ...
    null.p_primary_evidence, null.p_primary_accuracy, ...
    null.p_control_mean_evidence, null.p_control_mean_accuracy, ...
    'VariableNames', {'within_semantic_evidence', 'cross_semantic_evidence', ...
    'semantic_evidence_effect', 'within_semantic_accuracy', 'cross_semantic_accuracy', ...
    'semantic_accuracy_effect', 'control_transfer_mean_evidence', ...
    'control_transfer_mean_accuracy', 'control_within_evidence', ...
    'control_within_accuracy', 'p_primary_evidence', 'p_primary_accuracy', ...
    'p_control_mean_evidence', 'p_control_mean_accuracy'})];

if strlength(string(opts.OutputDir)) == 0
    if strcmp(roi_selection, 'old')
        output_name = 'roi_odor_context_template_loro';
    else
        output_name = sprintf('roi_%s_odor_context_template_loro', roi_selection);
    end
    output_dir = fullfile(beta_dir, output_name);
else
    output_dir = char(string(opts.OutputDir));
end

results = struct();
results.subject = struct('index', subjidx, 'name', subject_name, ...
    'name_real', subject_names{subjidx});
results.analysis = struct( ...
    'name', 'semantic_context_modulation_of_odor_identity', ...
    'roi_selection', roi_selection, ...
    'beta_alignment', 'odor', ...
    'cross_validation', 'leave-one-run-out', ...
    'classifier', 'nearest odor template by Pearson correlation', ...
    'evidence', 'correct-template similarity minus mean incorrect-template similarity', ...
    'semantic_context_order', {{'PERSON', 'FOOD', 'LOCATION'}}, ...
    'control_role', 'secondary training context only');
results.preprocessing = struct( ...
    'within_run_voxel_centering', true, ...
    'demean_patterns', true, ...
    'scale_voxels', false, ...
    'restrict_features_to_r2', false, ...
    'minimum_voxels', opts.MinVoxels, ...
    'held_out_run_excluded_from_all_templates', true);
results.odor_values = odor_values;
results.chance_accuracy = 1 / n_odors;
results.trial_metadata = trial_metadata;
results.roi_metadata = roi_metadata;
results.summary = summary;
results.semantic = struct('evidence', observed.semantic_evidence, ...
    'accuracy', observed.semantic_accuracy);
results.primary = struct( ...
    'within_context_evidence', observed.within_evidence, ...
    'cross_context_evidence', observed.cross_evidence, ...
    'evidence_effect', observed.primary_evidence_effect, ...
    'within_context_accuracy', observed.within_accuracy, ...
    'cross_context_accuracy', observed.cross_accuracy, ...
    'accuracy_effect', observed.primary_accuracy_effect);
results.control_transfer = struct( ...
    'evidence', observed.control_transfer_evidence, ...
    'accuracy', observed.control_transfer_accuracy, ...
    'mean_evidence', observed.control_transfer_mean_evidence, ...
    'mean_accuracy', observed.control_transfer_mean_accuracy);
results.control_within = struct('evidence', observed.control_within_evidence, ...
    'accuracy', observed.control_within_accuracy);
results.trial_predictions = observed.trial_predictions;
results.trial_evidence = observed.trial_evidence;
results.fold_accuracy = observed.fold_accuracy;
results.fold_evidence = observed.fold_evidence;
results.dimension_order = struct( ...
    'semantic', {{'roi', 'source_semantic_context', 'target_semantic_context'}}, ...
    'control_transfer', {{'roi', 'target_semantic_context'}}, ...
    'trial', {{'roi', 'source_context', 'target_context', 'trial'}}, ...
    'fold', {{'roi', 'source_context', 'target_context', 'fold'}});
results.null = null;
results.permutation = ['odor labels shuffled independently within context-by-run blocks; ' ...
    'the complete leave-one-run-out template decoder rerun'];
results.options = opts;
results.output_dir = output_dir;
results.runtime_minutes = toc(analysis_tic) / 60;
results.matlab_version = version;

if opts.SaveOutputs
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end
    output_file = fullfile(output_dir, ...
        sprintf('odor_context_template_subj%d_loro_results.mat', subjidx));
    save(output_file, 'results', '-v7.3');
    writetable(summary, fullfile(output_dir, ...
        sprintf('odor_context_template_subj%d_roi_summary.tsv', subjidx)), ...
        'FileType', 'text', 'Delimiter', '\t');
    fprintf('Saved subject results to %s\n', output_file);
end
end

function mri_root = resolve_mri_root(requested)
if strlength(string(requested)) > 0
    mri_root = char(string(requested));
    assert(isfolder(mri_root), 'MRIRoot does not exist: %s', mri_root);
    return;
end
candidates = {'/Users/qhyang/Desktop/OX_DATA/MRI', '/Volumes/ExtremeSSD/OX_DATA/MRI'};
mri_root = '';
for idx = 1:numel(candidates)
    if isfolder(candidates{idx})
        mri_root = candidates{idx};
        break;
    end
end
assert(~isempty(mri_root), 'Could not auto-detect MRIRoot.');
end

function modelmd = center_voxels_within_runs(modelmd, run_ids)
runs = unique(run_ids, 'stable');
for run_idx = 1:numel(runs)
    mask = run_ids == runs(run_idx);
    modelmd(:, mask) = modelmd(:, mask) - mean(modelmd(:, mask), 2);
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

function assert_same_geometry(image_header, reference_header, image_file, reference_file)
assert(isequal(image_header.dim, reference_header.dim), ...
    'Dimension mismatch between %s and %s.', image_file, reference_file);
assert(max(abs(image_header.mat(:) - reference_header.mat(:))) < 1e-4, ...
    'Affine mismatch between %s and %s.', image_file, reference_file);
end

function observed = package_observed(cells, n_rois, n_trials, n_folds)
observed.semantic_evidence = nan(n_rois, 3, 3);
observed.semantic_accuracy = nan(n_rois, 3, 3);
observed.within_evidence = nan(n_rois, 1);
observed.cross_evidence = nan(n_rois, 1);
observed.primary_evidence_effect = nan(n_rois, 1);
observed.within_accuracy = nan(n_rois, 1);
observed.cross_accuracy = nan(n_rois, 1);
observed.primary_accuracy_effect = nan(n_rois, 1);
observed.control_transfer_evidence = nan(n_rois, 3);
observed.control_transfer_accuracy = nan(n_rois, 3);
observed.control_transfer_mean_evidence = nan(n_rois, 1);
observed.control_transfer_mean_accuracy = nan(n_rois, 1);
observed.control_within_evidence = nan(n_rois, 1);
observed.control_within_accuracy = nan(n_rois, 1);
observed.trial_predictions = nan(n_rois, 4, 4, n_trials);
observed.trial_evidence = nan(n_rois, 4, 4, n_trials);
observed.fold_accuracy = nan(n_rois, 4, 4, n_folds);
observed.fold_evidence = nan(n_rois, 4, 4, n_folds);
for roi_idx = 1:n_rois
    if isempty(cells{roi_idx})
        continue;
    end
    item = cells{roi_idx};
    observed.semantic_evidence(roi_idx, :, :) = item.semantic.evidence;
    observed.semantic_accuracy(roi_idx, :, :) = item.semantic.accuracy;
    observed.within_evidence(roi_idx) = item.primary.within_context_evidence;
    observed.cross_evidence(roi_idx) = item.primary.cross_context_evidence;
    observed.primary_evidence_effect(roi_idx) = item.primary.evidence_effect;
    observed.within_accuracy(roi_idx) = item.primary.within_context_accuracy;
    observed.cross_accuracy(roi_idx) = item.primary.cross_context_accuracy;
    observed.primary_accuracy_effect(roi_idx) = item.primary.accuracy_effect;
    observed.control_transfer_evidence(roi_idx, :) = item.control_transfer.evidence;
    observed.control_transfer_accuracy(roi_idx, :) = item.control_transfer.accuracy;
    observed.control_transfer_mean_evidence(roi_idx) = item.control_transfer.mean_evidence;
    observed.control_transfer_mean_accuracy(roi_idx) = item.control_transfer.mean_accuracy;
    observed.control_within_evidence(roi_idx) = item.control_within.evidence;
    observed.control_within_accuracy(roi_idx) = item.control_within.accuracy;
    observed.trial_predictions(roi_idx, :, :, :) = item.trial_predictions;
    observed.trial_evidence(roi_idx, :, :, :) = item.trial_evidence;
    observed.fold_accuracy(roi_idx, :, :, :) = item.fold_accuracy;
    observed.fold_evidence(roi_idx, :, :, :) = item.fold_evidence;
end
end

function null = initialize_null(n_rois, n_permutations)
null.primary_evidence_effect = nan(n_rois, n_permutations);
null.primary_accuracy_effect = nan(n_rois, n_permutations);
null.control_transfer_mean_evidence = nan(n_rois, n_permutations);
null.control_transfer_mean_accuracy = nan(n_rois, n_permutations);
null.control_transfer_evidence = nan(n_rois, 3, n_permutations);
null.control_transfer_accuracy = nan(n_rois, 3, n_permutations);
end

function permuted = permute_odors_within_context_run(odors, contexts, run_ids)
permuted = odors(:);
contexts = string(contexts(:));
runs = unique(run_ids, 'stable');
context_values = ["PERSON", "FOOD", "LOCATION", "CONTROL"];
for run_idx = 1:numel(runs)
    for context_idx = 1:numel(context_values)
        indices = find(run_ids == runs(run_idx) & contexts == context_values(context_idx));
        if numel(indices) > 1
            permuted(indices) = permuted(indices(randperm(numel(indices))));
        end
    end
end
end

function null = store_permutation(null, cells, permutation_idx, n_rois)
for roi_idx = 1:n_rois
    if isempty(cells{roi_idx})
        continue;
    end
    item = cells{roi_idx};
    null.primary_evidence_effect(roi_idx, permutation_idx) = item.primary.evidence_effect;
    null.primary_accuracy_effect(roi_idx, permutation_idx) = item.primary.accuracy_effect;
    null.control_transfer_mean_evidence(roi_idx, permutation_idx) = ...
        item.control_transfer.mean_evidence;
    null.control_transfer_mean_accuracy(roi_idx, permutation_idx) = ...
        item.control_transfer.mean_accuracy;
    null.control_transfer_evidence(roi_idx, :, permutation_idx) = item.control_transfer.evidence;
    null.control_transfer_accuracy(roi_idx, :, permutation_idx) = item.control_transfer.accuracy;
end
end

function null = finish_null_inference(null, observed, valid_results)
n_rois = numel(valid_results);
null.p_primary_evidence = nan(n_rois, 1);
null.p_primary_accuracy = nan(n_rois, 1);
null.p_control_mean_evidence = nan(n_rois, 1);
null.p_control_mean_accuracy = nan(n_rois, 1);
null.p_control_target_evidence = nan(n_rois, 3);
null.p_control_target_accuracy = nan(n_rois, 3);
if isempty(null.primary_evidence_effect)
    return;
end
for roi_idx = find(valid_results)'
    null.p_primary_evidence(roi_idx) = empirical_upper_p( ...
        observed.primary_evidence_effect(roi_idx), ...
        null.primary_evidence_effect(roi_idx, :));
    null.p_primary_accuracy(roi_idx) = empirical_upper_p( ...
        observed.primary_accuracy_effect(roi_idx), ...
        null.primary_accuracy_effect(roi_idx, :));
    null.p_control_mean_evidence(roi_idx) = empirical_upper_p( ...
        observed.control_transfer_mean_evidence(roi_idx), ...
        null.control_transfer_mean_evidence(roi_idx, :));
    null.p_control_mean_accuracy(roi_idx) = empirical_upper_p( ...
        observed.control_transfer_mean_accuracy(roi_idx), ...
        null.control_transfer_mean_accuracy(roi_idx, :));
    for target_idx = 1:3
        null.p_control_target_evidence(roi_idx, target_idx) = empirical_upper_p( ...
            observed.control_transfer_evidence(roi_idx, target_idx), ...
            squeeze(null.control_transfer_evidence(roi_idx, target_idx, :)));
        null.p_control_target_accuracy(roi_idx, target_idx) = empirical_upper_p( ...
            observed.control_transfer_accuracy(roi_idx, target_idx), ...
            squeeze(null.control_transfer_accuracy(roi_idx, target_idx, :)));
    end
end
null.mean_primary_evidence_effect = mean(null.primary_evidence_effect, 2, 'omitnan');
null.mean_primary_accuracy_effect = mean(null.primary_accuracy_effect, 2, 'omitnan');
null.mean_control_transfer_evidence = mean(null.control_transfer_mean_evidence, 2, 'omitnan');
null.mean_control_transfer_accuracy = mean(null.control_transfer_mean_accuracy, 2, 'omitnan');
null.mean_control_target_evidence = mean(null.control_transfer_evidence, 3, 'omitnan');
null.mean_control_target_accuracy = mean(null.control_transfer_accuracy, 3, 'omitnan');
end

function p_value = empirical_upper_p(observed, null_values)
null_values = null_values(isfinite(null_values));
if isempty(null_values) || ~isfinite(observed)
    p_value = NaN;
else
    p_value = (1 + sum(null_values >= observed)) / (numel(null_values) + 1);
end
end
