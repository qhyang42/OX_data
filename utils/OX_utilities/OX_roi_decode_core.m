function results = OX_roi_decode_core(subjidx, target, varargin)
%OX_ROI_DECODE_CORE ROI-based context or odor decoding for the OX project.
%
%   results = OX_roi_decode_core(subjidx, target, Name, Value, ...)
%
% The function discovers all *_func_thr02 ROI masks in
% coreg/roi_decoding, maps them to GLMsingle rows through the functional
% gray-matter mask, and evaluates one decoder with one cross-validation
% scheme. Data preparation and permutation inference follow the current
% searchlight implementation.
%
% Required inputs:
%   subjidx                 Numeric subject index.
%   target                  'context' or 'odor'.
%
% Important name-value options:
%   'Decoder'               'template' (default), 'score', or 'svm'.
%   'CrossValidation'       'leave-one-run-out' (default),
%                           'leave-one-out', or '10-fold'.
%   'MRIRoot'               MRI root; auto-detected when empty.
%   'ROIDir'                ROI directory; defaults to coreg/roi_decoding.
%   'NumRuns'               Number of runs used to infer labels (80).
%   'RunLabels'             Exact trial-wise run labels; preferred.
%   'R2Threshold'           GLMsingle R2 cutoff (0.5).
%   'RestrictFeaturesToR2' Apply R2 cutoff to ROI features (true).
%   'MinVoxels'             Minimum usable features per ROI (10).
%   'MaxROIs'               Deterministic smoke-test subset (Inf).
%   'NumPermutations'       Within-run label permutations (1000).
%   'RandomSeed'            CV/permutation RNG seed (1).
%   'UseParallel'           Use Parallel Computing Toolbox (true).
%   'DemeanPatterns'        Remove each trial's spatial mean (true).
%   'ScaleVoxels'           Training-fold-only voxel z-scoring (true).
%   'SVMBoxConstraint'      Linear SVM box constraint (1).
%   'OutputDir'             Output directory; generated when empty.
%   'SaveOutputs'           Save results MAT file (true).
%
% Decoders:
%   template  Nearest training-class template by Pearson correlation.
%   score     Correct-template Pearson similarity minus the mean
%             similarity to all incorrect templates.
%   svm       Linear one-vs-one ECOC SVM with uniform class priors.
%
% Statistical inference is one-sided. Labels are shuffled within run and
% the complete cross-validation pipeline is rerun. Benjamini-Hochberg FDR
% and max-statistic FWE corrections are applied across all valid ROIs in
% this subject/target/decoder/CV analysis.

%% Parse inputs
p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'subjidx', @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x == round(x));
addRequired(p, 'target', @(x) any(strcmpi(string(x), ["context", "odor"])));
addParameter(p, 'Decoder', 'template', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'CrossValidation', 'leave-one-run-out', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'MRIRoot', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'ROIDir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'NumRuns', 80, @(x) isnumeric(x) && isscalar(x) && x >= 2 && x == round(x));
addParameter(p, 'RunLabels', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
addParameter(p, 'R2Threshold', 0.5, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'RestrictFeaturesToR2', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'MinVoxels', 10, @(x) isnumeric(x) && isscalar(x) && x >= 2 && x == round(x));
addParameter(p, 'MaxROIs', Inf, @(x) isnumeric(x) && isscalar(x) && (isinf(x) || (x >= 1 && x == round(x))));
addParameter(p, 'NumPermutations', 1000, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x == round(x));
addParameter(p, 'RandomSeed', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'UseParallel', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'DemeanPatterns', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ScaleVoxels', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'SVMBoxConstraint', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'OutputDir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'SaveOutputs', true, @(x) islogical(x) && isscalar(x));
parse(p, subjidx, target, varargin{:});
opts = p.Results;

target = lower(char(string(target)));
[decoder, metric_name, metric_null] = normalize_decoder(opts.Decoder);
[cv_method, cv_short_name] = normalize_cv_method(opts.CrossValidation);
opts.Decoder = decoder;
opts.CrossValidation = cv_method;

SUBJNAMES = {'240711_fMRI_OX_NWU_AS', ...
             '240723_fMRI_OX_NWU_LS', ...
             '240814_fMRI_OX_NWU_JN', ...
             '240816_fMRI_OX_NWU_RR', ...
             '241018_fMRI_OX_NWU_BN', ...
             '250117_fMRI_OX_NWU_VS'};
assert(subjidx >= 1 && subjidx <= numel(SUBJNAMES), ...
    'subjidx must be between 1 and %d.', numel(SUBJNAMES));
subjname = sprintf('subj_%d', subjidx);
subjname_real = SUBJNAMES{subjidx};

%% Resolve paths and load GLMsingle data
analysis_tic = tic;
mriroot = resolve_mri_root(opts.MRIRoot);
mridatapath = fullfile(mriroot, subjname, 'nifti');
base_outdir = fullfile(mridatapath, 'single_trial_by_category');
fit_file = fullfile(base_outdir, 'TYPED_FITHRF_GLMDENOISE_RR.mat');
gm_mask_file = fullfile(mridatapath, 'coreg', 'gm_mask_thr05_func.nii');
if strlength(string(opts.ROIDir)) == 0
    roi_dir = fullfile(mridatapath, 'coreg', 'roi_decoding');
else
    roi_dir = char(string(opts.ROIDir));
end

assert(isfile(fit_file), 'Missing GLMsingle file: %s', fit_file);
assert(isfile(gm_mask_file), 'Missing gray-matter mask: %s', gm_mask_file);
assert(isfolder(roi_dir), 'Missing ROI directory: %s', roi_dir);
if strcmp(decoder, 'svm')
    assert(exist('fitcecoc', 'file') == 2 && exist('templateSVM', 'file') == 2, ...
        'SVM decoding requires Statistics and Machine Learning Toolbox.');
end

fprintf('\n[%s | %s | %s] Loading %s\n', upper(target), upper(decoder), upper(cv_short_name), fit_file);
S = load(fit_file, 'modelmd', 'R2');
assert(isfield(S, 'modelmd') && isfield(S, 'R2'), ...
    'GLMsingle file must contain modelmd and R2.');
modelmd = squeeze(S.modelmd);
R2 = S.R2(:);
clear S

[odor_raw, category_raw] = OX_get_odor(subjname);
[y, class_values] = encode_target(target, odor_raw, category_raw);
nTrials = numel(y);
nClasses = numel(unique(y));
chance = 1 / nClasses;
if isnan(metric_null)
    metric_null = chance;
end

assert(ismatrix(modelmd), 'squeeze(modelmd) must produce a 2-D matrix.');
assert(size(modelmd, 2) == nTrials, ...
    'modelmd has %d columns but labels contain %d trials.', size(modelmd, 2), nTrials);
assert(numel(R2) == size(modelmd, 1), ...
    'R2 has %d elements but modelmd has %d voxel rows.', numel(R2), size(modelmd, 1));

run_ids = make_run_ids(nTrials, opts.NumRuns, opts.RunLabels);
unique_runs = unique(run_ids, 'stable');
nRuns = numel(unique_runs);
modelmd = center_model_within_runs(modelmd, run_ids);

rng(opts.RandomSeed, 'twister');
[cv_ids, run_fold_table] = make_cv_ids(cv_method, run_ids, 10);
validate_training_folds(y, cv_ids, nClasses);
nFolds = numel(unique(cv_ids, 'stable'));

fprintf('Subject: %s (%s) | Trials: %d | Runs: %d | Classes: %d\n', ...
    subjname, subjname_real, nTrials, nRuns, nClasses);
fprintf('Decoder: %s | Metric: %s | Null: %.4f | CV folds: %d\n', ...
    decoder, metric_name, metric_null, nFolds);
fprintf('Subtracted each voxel''s within-run mean before ROI extraction.\n');
print_class_counts(y, class_values);

if opts.NumPermutations > 0
    movable_runs = false(nRuns, 1);
    for ri = 1:nRuns
        movable_runs(ri) = numel(unique(y(run_ids == unique_runs(ri)))) > 1;
    end
    assert(any(movable_runs), ...
        'Within-run permutation cannot change labels in any run.');
    fprintf('Runs contributing to within-run label shuffling: %d/%d\n', sum(movable_runs), nRuns);
end

if strcmp(decoder, 'svm') && strcmp(cv_method, 'leave-one-out')
    projected_fits = nTrials * max(1, opts.NumPermutations + 1);
    warning(['SVM leave-one-out requires %d model fits per valid ROI for this run. ' ...
             'Use leave-one-run-out or 10-fold for substantially lower runtime.'], projected_fits);
end

%% Load gray-matter mask and ROI feature sets
gm_header = spm_vol(gm_mask_file);
gm_mask = spm_read_vols(gm_header) > 0;
gm_inds = find(gm_mask);
assert(size(modelmd, 1) == numel(gm_inds), ...
    ['Mapping check failed: modelmd has %d rows but find(gm_mask) has %d voxels. ' ...
     'modelmd row i must correspond to gm_inds(i).'], size(modelmd, 1), numel(gm_inds));

[roi_files, roi_stems] = discover_roi_files(roi_dir, opts.MaxROIs);
nROIs = numel(roi_files);
manifest = load_roi_manifest(fullfile(roi_dir, 'roi_manifest.tsv'));

model_index_volume = zeros(size(gm_mask), 'uint32');
model_index_volume(gm_inds) = uint32(1:numel(gm_inds));
roi_features = cell(nROIs, 1);
roi_names = strings(nROIs, 1);
roi_sources = strings(nROIs, 1);
roi_hemispheres = strings(nROIs, 1);
roi_labels = strings(nROIs, 1);
n_mask_voxels = zeros(nROIs, 1);
n_gm_overlap = zeros(nROIs, 1);
n_features_used = zeros(nROIs, 1);
status = repmat("ok", nROIs, 1);

fprintf('Loading %d ROI masks from %s\n', nROIs, roi_dir);
for roi_idx = 1:nROIs
    roi_file = roi_files{roi_idx};
    roi_header = spm_vol(roi_file);
    assert_same_geometry(roi_header, gm_header, roi_file, gm_mask_file);
    roi_mask = spm_read_vols(roi_header) > 0;
    n_mask_voxels(roi_idx) = nnz(roi_mask);

    overlap_linear = find(roi_mask & gm_mask);
    n_gm_overlap(roi_idx) = numel(overlap_linear);
    feature_idx = double(model_index_volume(overlap_linear));
    feature_idx = feature_idx(feature_idx > 0);
    if opts.RestrictFeaturesToR2
        feature_idx = feature_idx(isfinite(R2(feature_idx)) & R2(feature_idx) > opts.R2Threshold);
    end
    if ~isempty(feature_idx)
        feature_idx = feature_idx(all(isfinite(modelmd(feature_idx, :)), 2));
    end
    roi_features{roi_idx} = feature_idx(:)';
    n_features_used(roi_idx) = numel(feature_idx);

    roi_names(roi_idx) = erase(string(roi_stems{roi_idx}), "_func_thr02");
    [roi_sources(roi_idx), roi_hemispheres(roi_idx), roi_labels(roi_idx)] = ...
        lookup_manifest_metadata(manifest, roi_stems{roi_idx}, roi_names(roi_idx));

    if n_gm_overlap(roi_idx) < opts.MinVoxels
        status(roi_idx) = "insufficient_gm_overlap";
    elseif n_features_used(roi_idx) < opts.MinVoxels
        status(roi_idx) = "insufficient_features";
    end
end

valid_roi = status == "ok";
fprintf('Valid ROIs: %d/%d | MinVoxels: %d\n', sum(valid_roi), nROIs, opts.MinVoxels);
assert(any(valid_roi), 'No ROI has at least %d usable features.', opts.MinVoxels);

%% Parallel setup and observed analysis
use_parallel = opts.UseParallel && license('test', 'Distrib_Computing_Toolbox');
if opts.UseParallel && ~use_parallel
    warning('Parallel Computing Toolbox unavailable; running serially.');
end
if use_parallel && isempty(gcp('nocreate'))
    parpool;
end

observed_metric = nan(nROIs, 1);
accuracy = nan(nROIs, 1);
balanced_accuracy = nan(nROIs, 1);
predictions_cell = cell(nROIs, 1);
trial_scores_cell = cell(nROIs, 1);
fold_metrics_cell = cell(nROIs, 1);
confusion_cell = cell(nROIs, 1);

fprintf('Running observed ROI analysis...\n');
observed_tic = tic;
min_voxels = opts.MinVoxels;
demean_patterns = opts.DemeanPatterns;
scale_voxels = opts.ScaleVoxels;
svm_box_constraint = opts.SVMBoxConstraint;
if use_parallel
    model_const = parallel.pool.Constant(modelmd);
    model_cleanup = onCleanup(@() delete(model_const));
    parfor roi_idx = 1:nROIs
        if valid_roi(roi_idx)
            [observed_metric(roi_idx), accuracy(roi_idx), balanced_accuracy(roi_idx), ...
             predictions_cell{roi_idx}, trial_scores_cell{roi_idx}, ...
             fold_metrics_cell{roi_idx}, confusion_cell{roi_idx}] = evaluate_one_roi( ...
                model_const.Value, roi_features{roi_idx}, y, cv_ids, nClasses, ...
                decoder, min_voxels, demean_patterns, scale_voxels, ...
                svm_box_constraint, true);
        end
    end
else
    for roi_idx = 1:nROIs
        if valid_roi(roi_idx)
            [observed_metric(roi_idx), accuracy(roi_idx), balanced_accuracy(roi_idx), ...
             predictions_cell{roi_idx}, trial_scores_cell{roi_idx}, ...
             fold_metrics_cell{roi_idx}, confusion_cell{roi_idx}] = evaluate_one_roi( ...
                modelmd, roi_features{roi_idx}, y, cv_ids, nClasses, decoder, ...
                min_voxels, demean_patterns, scale_voxels, ...
                svm_box_constraint, true);
        end
    end
end
fprintf('Observed ROI analysis completed in %.2f minutes.\n', toc(observed_tic) / 60);

decoder_failed = valid_roi & ~isfinite(observed_metric);
status(decoder_failed) = "decoder_failed";
valid_results = isfinite(observed_metric);
assert(any(valid_results), 'No ROI produced a finite observed metric.');

predictions = cells_to_matrix(predictions_cell, nROIs, nTrials);
trial_scores = cells_to_matrix(trial_scores_cell, nROIs, nTrials);
fold_metrics = cells_to_matrix(fold_metrics_cell, nROIs, nFolds);
confusion_matrices = cells_to_confusions(confusion_cell, nROIs, nClasses);
confusion_matrices_normalized = normalize_confusions(confusion_matrices);

%% Permutation inference
null_distribution = nan(nROIs, opts.NumPermutations);
p_empirical = nan(nROIs, 1);
q_fdr = nan(nROIs, 1);
p_fwe_maxstat = nan(nROIs, 1);
sig_fdr05 = false(nROIs, 1);
sig_fwe05 = false(nROIs, 1);
max_null = nan(opts.NumPermutations, 1);
fwe_metric_threshold05 = NaN;

if opts.NumPermutations > 0
    fprintf('Running %d within-run label permutations...\n', opts.NumPermutations);
    permutation_tic = tic;
    for perm_idx = 1:opts.NumPermutations
        y_perm = permute_within_blocks(y, run_ids);
        perm_metric = nan(nROIs, 1);
        if use_parallel
            parfor roi_idx = 1:nROIs
                if valid_results(roi_idx)
                    perm_metric(roi_idx) = evaluate_one_roi( ...
                        model_const.Value, roi_features{roi_idx}, y_perm, cv_ids, ...
                        nClasses, decoder, min_voxels, demean_patterns, ...
                        scale_voxels, svm_box_constraint, false);
                end
            end
        else
            for roi_idx = 1:nROIs
                if valid_results(roi_idx)
                    perm_metric(roi_idx) = evaluate_one_roi( ...
                        modelmd, roi_features{roi_idx}, y_perm, cv_ids, nClasses, ...
                        decoder, min_voxels, demean_patterns, ...
                        scale_voxels, svm_box_constraint, false);
                end
            end
        end
        null_distribution(:, perm_idx) = perm_metric;
        comparable = valid_results & isfinite(perm_metric);
        if any(comparable)
            max_null(perm_idx) = max(perm_metric(comparable) - metric_null);
        end
        if perm_idx == 1 || mod(perm_idx, 10) == 0 || perm_idx == opts.NumPermutations
            fprintf('  permutation %d/%d | elapsed %.1f min\n', ...
                perm_idx, opts.NumPermutations, toc(permutation_tic) / 60);
        end
    end

    for roi_idx = find(valid_results)'
        roi_null = null_distribution(roi_idx, :);
        roi_null = roi_null(isfinite(roi_null));
        p_empirical(roi_idx) = (1 + sum(roi_null >= observed_metric(roi_idx))) / ...
            (numel(roi_null) + 1);
        valid_max = max_null(isfinite(max_null));
        p_fwe_maxstat(roi_idx) = (1 + sum(valid_max >= ...
            (observed_metric(roi_idx) - metric_null))) / (numel(valid_max) + 1);
    end
    q_fdr(valid_results) = bh_fdr_qvalues(p_empirical(valid_results));
    sig_fdr05 = q_fdr < 0.05;
    sig_fwe05 = p_fwe_maxstat < 0.05;
    valid_max = sort(max_null(isfinite(max_null)));
    if ~isempty(valid_max)
        threshold_index = min(numel(valid_max), max(1, ceil(0.95 * (numel(valid_max) + 1))));
        fwe_metric_threshold05 = metric_null + valid_max(threshold_index);
    end
    fprintf('Permutation inference: FDR q<.05 in %d ROIs; max-stat FWE p<.05 in %d ROIs.\n', ...
        sum(sig_fdr05), sum(sig_fwe05));
end

null_mean = mean(null_distribution, 2, 'omitnan');
null_sd = std(null_distribution, 0, 2, 'omitnan');
null_percentile95 = row_percentile(null_distribution, 95);
metric_minus_null = observed_metric - metric_null;
accuracy_minus_chance = accuracy - chance;
mean_evidence = nan(nROIs, 1);
if strcmp(decoder, 'score')
    mean_evidence = observed_metric;
end

%% Package and save results
roi_file_strings = string(roi_files(:));
summary = table(roi_names, roi_sources, roi_hemispheres, roi_labels, roi_file_strings, ...
    status, n_mask_voxels, n_gm_overlap, n_features_used, observed_metric, ...
    metric_minus_null, accuracy, accuracy_minus_chance, balanced_accuracy, ...
    mean_evidence, p_empirical, q_fdr, p_fwe_maxstat, sig_fdr05, sig_fwe05, ...
    null_mean, null_sd, null_percentile95, ...
    'VariableNames', {'roi_name', 'source', 'hemisphere', 'contributing_labels', ...
    'mask_file', 'status', 'n_mask_voxels', 'n_gm_overlap', 'n_features_used', ...
    'observed_metric', 'metric_minus_null', 'accuracy', 'accuracy_minus_chance', ...
    'balanced_accuracy', 'mean_evidence', 'p_empirical', 'q_fdr', ...
    'p_fwe_maxstat', 'sig_fdr05', 'sig_fwe05', 'null_mean', 'null_sd', ...
    'null_percentile95'});

if strlength(string(opts.OutputDir)) == 0
    output_dir = fullfile(base_outdir, ...
        sprintf('roi_decoding_%s_%s_%s', target, decoder, cv_short_name));
else
    output_dir = char(string(opts.OutputDir));
end

results = struct();
results.subject = struct('index', subjidx, 'name', subjname, 'name_real', subjname_real);
results.analysis = struct('target', target, 'decoder', decoder, ...
    'cross_validation', cv_method, 'metric', metric_name, 'metric_null', metric_null);
results.classes = struct('values', {class_values}, ...
    'counts', accumarray(y, 1, [nClasses, 1]), 'chance', chance);
results.preprocessing = struct( ...
    'within_run_voxel_centering', true, ...
    'demean_patterns', opts.DemeanPatterns, ...
    'scale_voxels_with_training_fold_only', opts.ScaleVoxels, ...
    'restrict_features_to_r2', opts.RestrictFeaturesToR2, ...
    'r2_threshold', opts.R2Threshold, ...
    'mapping_assertion', 'modelmd row i == find(gm_mask)(i)');
results.roi_metadata = summary(:, 1:9);
results.summary = summary;
results.true_labels = y;
results.run_ids = run_ids;
results.cv_fold_ids = cv_ids;
results.run_fold_table = run_fold_table;
results.predictions = predictions;
results.trial_scores = trial_scores;
results.fold_metrics = fold_metrics;
results.confusion_matrices = confusion_matrices;
results.confusion_matrices_normalized = confusion_matrices_normalized;
results.confusion_matrix_dimension_order = {'roi', 'true_class', 'predicted_class'};
results.null_distribution = null_distribution;
results.max_null = max_null;
results.fwe_metric_threshold05 = fwe_metric_threshold05;
results.options = opts;
results.output_dir = output_dir;
results.permutation = ['labels shuffled independently within run; complete ROI feature ' ...
    'preprocessing, decoder fitting, and cross-validation pipeline rerun'];
results.fdr_family = 'all valid ROI masks in this subject/target/decoder/CV result';
results.svm = struct('kernel', 'linear', 'coding', 'one-vs-one ECOC', ...
    'box_constraint', opts.SVMBoxConstraint, 'class_prior', 'uniform', ...
    'internal_standardization', false);
results.runtime_minutes = toc(analysis_tic) / 60;
results.matlab_version = version;

if opts.SaveOutputs
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end
    output_file = fullfile(output_dir, ...
        sprintf('%s_%s_subj%d_%s_results.mat', target, decoder, subjidx, cv_short_name));
    save(output_file, 'results', '-v7.3');
    fprintf('Saved ROI results to %s\n', output_file);
end
end

%% ------------------------------------------------------------------------
function [decoder, metric_name, metric_null] = normalize_decoder(requested)
decoder = lower(strtrim(char(string(requested))));
switch decoder
    case {'template', 'template-comparison', 'correlation'}
        decoder = 'template';
        metric_name = 'accuracy';
        metric_null = NaN; % set to chance after labels are known
    case {'score', 'evidence'}
        decoder = 'score';
        metric_name = 'mean_evidence';
        metric_null = 0;
    case {'svm', 'linear-svm'}
        decoder = 'svm';
        metric_name = 'accuracy';
        metric_null = NaN; % set to chance after labels are known
    otherwise
        error('Decoder must be ''template'', ''score'', or ''svm''.');
end
end

%% ------------------------------------------------------------------------
function [cv_method, cv_short_name] = normalize_cv_method(requested)
requested = lower(strtrim(char(string(requested))));
switch requested
    case {'leave-one-run-out', 'loro'}
        cv_method = 'leave-one-run-out';
        cv_short_name = 'loro';
    case {'leave-one-out', 'loo'}
        cv_method = 'leave-one-out';
        cv_short_name = 'loo';
    case {'10-fold', '10fold', 'kfold10', 'ten-fold'}
        cv_method = '10-fold';
        cv_short_name = 'kfold10';
    otherwise
        error(['CrossValidation must be ''leave-one-run-out'' (''loro''), ' ...
               '''leave-one-out'' (''loo''), or ''10-fold'' (''kfold10'').']);
end
end

%% ------------------------------------------------------------------------
function mriroot = resolve_mri_root(requested_root)
if strlength(string(requested_root)) > 0
    mriroot = char(string(requested_root));
    assert(isfolder(mriroot), 'MRIRoot does not exist: %s', mriroot);
    return;
end
candidates = {'/Users/qhyang/Desktop/OX_DATA/MRI', '/Volumes/ExtremeSSD/OX_DATA/MRI'};
mriroot = '';
for i = 1:numel(candidates)
    if isfolder(candidates{i})
        mriroot = candidates{i};
        break;
    end
end
assert(~isempty(mriroot), 'Could not auto-detect MRIRoot. Pass ''MRIRoot'' explicitly.');
end

%% ------------------------------------------------------------------------
function [y, class_values] = encode_target(target, odor_raw, category_raw)
if strcmp(target, 'context')
    class_order = ["PERSON", "FOOD", "LOCATION", "CONTROL"];
    labels = upper(strtrim(string(category_raw(:))));
    [known, y] = ismember(labels, class_order);
    assert(all(known), 'Unknown context labels: %s', strjoin(unique(labels(~known)), ', '));
    class_values = cellstr(class_order(:));
else
    odor_raw = odor_raw(:);
    if isnumeric(odor_raw) || islogical(odor_raw)
        class_values = unique(odor_raw, 'sorted');
        [known, y] = ismember(odor_raw, class_values);
    else
        labels = string(odor_raw);
        values = unique(labels, 'sorted');
        [known, y] = ismember(labels, values);
        class_values = cellstr(values);
    end
    assert(all(known), 'Unable to encode one or more odor labels.');
end
y = double(y(:));
end

%% ------------------------------------------------------------------------
function run_ids = make_run_ids(nTrials, nRuns, supplied)
if ~isempty(supplied)
    run_ids = supplied(:);
    assert(numel(run_ids) == nTrials, 'RunLabels has %d values; expected %d.', numel(run_ids), nTrials);
    assert(all(isfinite(run_ids)), 'RunLabels must be finite.');
    return;
end
assert(mod(nTrials, nRuns) == 0, ...
    'Cannot infer equal contiguous runs: %d trials is not divisible by %d runs.', nTrials, nRuns);
trials_per_run = nTrials / nRuns;
run_ids = repelem((1:nRuns)', trials_per_run);
warning(['RunLabels inferred as %d contiguous equal-sized runs (%d trials/run). ' ...
         'Pass exact RunLabels for the final analysis when available.'], nRuns, trials_per_run);
end

%% ------------------------------------------------------------------------
function modelmd = center_model_within_runs(modelmd, run_ids)
runs = unique(run_ids, 'stable');
for ri = 1:numel(runs)
    trial_mask = run_ids == runs(ri);
    modelmd(:, trial_mask) = modelmd(:, trial_mask) - mean(modelmd(:, trial_mask), 2);
end
end

%% ------------------------------------------------------------------------
function [cv_ids, run_fold_table] = make_cv_ids(method, run_ids, nKfold)
runs = unique(run_ids, 'stable');
switch method
    case 'leave-one-run-out'
        cv_ids = run_ids;
        run_fold_table = table(runs, runs, 'VariableNames', {'run_id', 'fold_id'});
    case 'leave-one-out'
        cv_ids = (1:numel(run_ids))';
        run_fold_table = table(runs, nan(size(runs)), 'VariableNames', {'run_id', 'fold_id'});
    case '10-fold'
        assert(numel(runs) >= nKfold, '10-fold CV requires at least 10 runs.');
        order = randperm(numel(runs));
        folds_by_position = zeros(numel(runs), 1);
        folds_by_position(order) = mod((0:numel(runs)-1)', nKfold) + 1;
        [known, run_position] = ismember(run_ids, runs);
        assert(all(known), 'Failed to map one or more trials to a run.');
        cv_ids = folds_by_position(run_position);
        run_fold_table = table(runs, folds_by_position, 'VariableNames', {'run_id', 'fold_id'});
end
cv_ids = double(cv_ids(:));
end

%% ------------------------------------------------------------------------
function validate_training_folds(y, cv_ids, nClasses)
folds = unique(cv_ids, 'stable');
for fi = 1:numel(folds)
    training_labels = unique(y(cv_ids ~= folds(fi)));
    missing = setdiff((1:nClasses)', training_labels);
    assert(isempty(missing), ...
        'Fold %g training data lack class(es) %s.', folds(fi), mat2str(missing'));
end
end

%% ------------------------------------------------------------------------
function print_class_counts(y, class_values)
counts = accumarray(y, 1);
fprintf('Class counts:\n');
for ci = 1:numel(counts)
    if isnumeric(class_values)
        label = num2str(class_values(ci));
    else
        label = char(string(class_values{ci}));
    end
    fprintf('  %2d  %-20s  %d\n', ci, label, counts(ci));
end
end

%% ------------------------------------------------------------------------
function [roi_files, roi_stems] = discover_roi_files(roi_dir, max_rois)
entries = [dir(fullfile(roi_dir, '*_func_thr02.nii')); ...
           dir(fullfile(roi_dir, '*_func_thr02.nii.gz'))];
assert(~isempty(entries), 'No *_func_thr02 NIfTI masks found in %s.', roi_dir);
names = string({entries.name})';
stems = arrayfun(@nifti_stem, names);
[unique_stems, ~, group] = unique(lower(stems));
duplicates = unique_stems(accumarray(group, 1) > 1);
assert(isempty(duplicates), 'Duplicate ROI basenames found: %s', strjoin(duplicates, ', '));
[~, order] = sort(lower(names));
entries = entries(order);
names = names(order);
stems = stems(order);
if isfinite(max_rois) && max_rois < numel(entries)
    entries = entries(1:max_rois);
    names = names(1:max_rois);
    stems = stems(1:max_rois);
    warning('MaxROIs=%d: running a smoke-test subset, not the full ROI analysis.', max_rois);
end
roi_files = cellstr(fullfile(string({entries.folder})', names));
roi_stems = cellstr(stems);
end

%% ------------------------------------------------------------------------
function stem = nifti_stem(filename)
stem = regexprep(string(filename), '\.nii(\.gz)?$', '', 'ignorecase');
end

%% ------------------------------------------------------------------------
function manifest = load_roi_manifest(manifest_file)
if isfile(manifest_file)
    manifest = readtable(manifest_file, 'FileType', 'text', 'Delimiter', '\t', 'TextType', 'string');
    required = ["source", "contributing_labels", "hemisphere", "output_file"];
    assert(all(ismember(required, string(manifest.Properties.VariableNames))), ...
        'ROI manifest lacks one or more required columns: %s', manifest_file);
    manifest.output_stem = arrayfun(@nifti_stem, manifest.output_file);
else
    warning('ROI manifest not found; source and label metadata will be inferred where possible.');
    manifest = table();
end
end

%% ------------------------------------------------------------------------
function [source, hemisphere, labels] = lookup_manifest_metadata(manifest, roi_stem, roi_name)
source = "unknown";
labels = "";
hemisphere = infer_hemisphere(roi_name);
if isempty(manifest)
    return;
end
match = find(strcmpi(manifest.output_stem, string(roi_stem)));
assert(numel(match) <= 1, 'ROI manifest has duplicate rows for %s.', roi_stem);
if isscalar(match)
    source = string(manifest.source(match));
    hemisphere = string(manifest.hemisphere(match));
    labels = string(manifest.contributing_labels(match));
end
end

%% ------------------------------------------------------------------------
function hemisphere = infer_hemisphere(roi_name)
if endsWith(roi_name, "_L")
    hemisphere = "L";
elseif endsWith(roi_name, "_R")
    hemisphere = "R";
elseif endsWith(roi_name, "_bilateral")
    hemisphere = "bilateral";
else
    hemisphere = "unknown";
end
end

%% ------------------------------------------------------------------------
function assert_same_geometry(image_header, reference_header, image_file, reference_file)
assert(isequal(image_header.dim, reference_header.dim), ...
    'Dimension mismatch between %s and %s.', image_file, reference_file);
assert(max(abs(image_header.mat(:) - reference_header.mat(:))) < 1e-4, ...
    'Affine mismatch between %s and %s.', image_file, reference_file);
end

%% ------------------------------------------------------------------------
function [metric, accuracy, balanced_accuracy, predictions, trial_scores, fold_metrics, confusion] = ...
        evaluate_one_roi(modelmd, feature_idx, y, cv_ids, nClasses, decoder, ...
        min_voxels, demean_patterns, scale_voxels, svm_box_constraint, collect_details)
metric = NaN;
accuracy = NaN;
balanced_accuracy = NaN;
predictions = [];
trial_scores = [];
fold_metrics = [];
confusion = [];
if numel(feature_idx) < min_voxels
    return;
end

X = double(modelmd(feature_idx, :)'); % trials x voxels
finite_features = all(isfinite(X), 1);
X = X(:, finite_features);
if size(X, 2) < min_voxels
    return;
end
if demean_patterns
    X = X - mean(X, 2);
end

is_loo = numel(unique(cv_ids)) == numel(y);
switch decoder
    case 'template'
        if is_loo
            predictions = template_leave_one_out(X, y, nClasses, min_voxels, scale_voxels);
        else
            predictions = template_grouped_folds(X, y, cv_ids, min_voxels, scale_voxels);
        end
        if all(isfinite(predictions))
            accuracy = mean(predictions == y);
            balanced_accuracy = calculate_balanced_accuracy(y, predictions, nClasses);
            metric = accuracy;
        end
    case 'score'
        if is_loo
            trial_scores = score_leave_one_out(X, y, nClasses, min_voxels, scale_voxels);
        else
            trial_scores = score_grouped_folds(X, y, cv_ids, min_voxels, scale_voxels);
        end
        if all(isfinite(trial_scores))
            metric = mean(trial_scores);
        end
    case 'svm'
        predictions = svm_grouped_folds(X, y, cv_ids, nClasses, min_voxels, ...
            scale_voxels, svm_box_constraint);
        if all(isfinite(predictions))
            accuracy = mean(predictions == y);
            balanced_accuracy = calculate_balanced_accuracy(y, predictions, nClasses);
            metric = accuracy;
        end
end

if ~collect_details
    predictions = [];
    trial_scores = [];
    return;
end

folds = unique(cv_ids, 'stable');
fold_metrics = nan(1, numel(folds));
for fi = 1:numel(folds)
    test_mask = cv_ids == folds(fi);
    if strcmp(decoder, 'score')
        fold_metrics(fi) = mean(trial_scores(test_mask), 'omitnan');
    else
        fold_metrics(fi) = mean(predictions(test_mask) == y(test_mask), 'omitnan');
    end
end
if ~strcmp(decoder, 'score') && all(isfinite(predictions))
    confusion = accumarray([y, predictions], 1, [nClasses, nClasses]);
end
end

%% ------------------------------------------------------------------------
function predictions = template_leave_one_out(X, y, nClasses, min_voxels, scale_voxels)
nTrials = size(X, 1);
predictions = nan(nTrials, 1);
class_counts = accumarray(y, 1, [nClasses, 1]);
assert(all(class_counts >= 2), 'Leave-one-out template decoding requires at least two trials per class.');
class_sums = zeros(nClasses, size(X, 2));
for ci = 1:nClasses
    class_sums(ci, :) = sum(X(y == ci, :), 1);
end
class_means = class_sums ./ class_counts;
[global_mu, total_centered_ss] = global_scaling_terms(X, scale_voxels);

for ti = 1:nTrials
    templates = class_means;
    test_class = y(ti);
    templates(test_class, :) = (class_sums(test_class, :) - X(ti, :)) ./ ...
        (class_counts(test_class) - 1);
    [templates, Xtest, usable] = apply_loo_scaling( ...
        templates, X(ti, :), global_mu, total_centered_ss, ti, X, scale_voxels);
    if sum(usable) < min_voxels
        continue;
    end
    predictions(ti) = predict_from_templates(Xtest, templates, (1:nClasses)');
end
end

%% ------------------------------------------------------------------------
function predictions = template_grouped_folds(X, y, cv_ids, min_voxels, scale_voxels)
predictions = nan(size(y));
folds = unique(cv_ids, 'stable');
for fi = 1:numel(folds)
    test_mask = cv_ids == folds(fi);
    train_mask = ~test_mask;
    [Xtrain, Xtest, usable] = prepare_fold_data(X(train_mask, :), X(test_mask, :), scale_voxels);
    if sum(usable) < min_voxels
        continue;
    end
    ytrain = y(train_mask);
    class_labels = unique(ytrain, 'sorted');
    templates = class_templates(Xtrain, ytrain, class_labels);
    predictions(test_mask) = predict_from_templates(Xtest, templates, class_labels);
end
end

%% ------------------------------------------------------------------------
function scores = score_leave_one_out(X, y, nClasses, min_voxels, scale_voxels)
nTrials = size(X, 1);
scores = nan(nTrials, 1);
class_counts = accumarray(y, 1, [nClasses, 1]);
assert(all(class_counts >= 2), 'Leave-one-out score decoding requires at least two trials per class.');
class_sums = zeros(nClasses, size(X, 2));
for ci = 1:nClasses
    class_sums(ci, :) = sum(X(y == ci, :), 1);
end
class_means = class_sums ./ class_counts;
[global_mu, total_centered_ss] = global_scaling_terms(X, scale_voxels);

for ti = 1:nTrials
    templates = class_means;
    test_class = y(ti);
    templates(test_class, :) = (class_sums(test_class, :) - X(ti, :)) ./ ...
        (class_counts(test_class) - 1);
    [templates, Xtest, usable] = apply_loo_scaling( ...
        templates, X(ti, :), global_mu, total_centered_ss, ti, X, scale_voxels);
    if sum(usable) < min_voxels
        continue;
    end
    scores(ti) = evidence_from_templates(Xtest, templates, (1:nClasses)', test_class);
end
end

%% ------------------------------------------------------------------------
function scores = score_grouped_folds(X, y, cv_ids, min_voxels, scale_voxels)
scores = nan(size(y));
folds = unique(cv_ids, 'stable');
for fi = 1:numel(folds)
    test_mask = cv_ids == folds(fi);
    train_mask = ~test_mask;
    [Xtrain, Xtest, usable] = prepare_fold_data(X(train_mask, :), X(test_mask, :), scale_voxels);
    if sum(usable) < min_voxels
        continue;
    end
    ytrain = y(train_mask);
    class_labels = unique(ytrain, 'sorted');
    templates = class_templates(Xtrain, ytrain, class_labels);
    scores(test_mask) = evidence_from_templates(Xtest, templates, class_labels, y(test_mask));
end
end

%% ------------------------------------------------------------------------
function predictions = svm_grouped_folds(X, y, cv_ids, nClasses, min_voxels, scale_voxels, box_constraint)
predictions = nan(size(y));
folds = unique(cv_ids, 'stable');
for fi = 1:numel(folds)
    test_mask = cv_ids == folds(fi);
    train_mask = ~test_mask;
    [Xtrain, Xtest, usable] = prepare_fold_data(X(train_mask, :), X(test_mask, :), scale_voxels);
    if sum(usable) < min_voxels
        continue;
    end
    ytrain = y(train_mask);
    learner = templateSVM('KernelFunction', 'linear', ...
        'BoxConstraint', box_constraint, 'Standardize', false);
    model = fitcecoc(Xtrain, ytrain, 'Learners', learner, ...
        'Coding', 'onevsone', 'ClassNames', (1:nClasses)', 'Prior', 'uniform');
    predictions(test_mask) = predict(model, Xtest);
end
end

%% ------------------------------------------------------------------------
function [Xtrain, Xtest, usable] = prepare_fold_data(Xtrain, Xtest, scale_voxels)
if scale_voxels
    train_mu = mean(Xtrain, 1);
    train_sd = std(Xtrain, 0, 1);
    usable = isfinite(train_mu) & isfinite(train_sd) & train_sd > eps;
    Xtrain = (Xtrain(:, usable) - train_mu(usable)) ./ train_sd(usable);
    Xtest = (Xtest(:, usable) - train_mu(usable)) ./ train_sd(usable);
else
    usable = all(isfinite(Xtrain), 1) & all(isfinite(Xtest), 1);
    Xtrain = Xtrain(:, usable);
    Xtest = Xtest(:, usable);
end
end

%% ------------------------------------------------------------------------
function [global_mu, total_centered_ss] = global_scaling_terms(X, scale_voxels)
if scale_voxels
    global_mu = mean(X, 1);
    total_centered_ss = sum((X - global_mu).^2, 1);
else
    global_mu = [];
    total_centered_ss = [];
end
end

%% ------------------------------------------------------------------------
function [templates, Xtest, usable] = apply_loo_scaling( ...
        templates, Xtest, global_mu, total_centered_ss, test_idx, X, scale_voxels)
if scale_voxels
    nTrials = size(X, 1);
    centered_test = X(test_idx, :) - global_mu;
    train_mu = global_mu - centered_test ./ (nTrials - 1);
    train_ss = total_centered_ss - (nTrials / (nTrials - 1)) .* centered_test.^2;
    train_sd = sqrt(max(train_ss, 0) ./ (nTrials - 2));
    usable = isfinite(train_mu) & isfinite(train_sd) & train_sd > eps;
    templates = (templates(:, usable) - train_mu(usable)) ./ train_sd(usable);
    Xtest = (Xtest(usable) - train_mu(usable)) ./ train_sd(usable);
else
    usable = all(isfinite(templates), 1) & isfinite(Xtest);
    templates = templates(:, usable);
    Xtest = Xtest(usable);
end
end

%% ------------------------------------------------------------------------
function templates = class_templates(Xtrain, ytrain, class_labels)
templates = nan(numel(class_labels), size(Xtrain, 2));
for ci = 1:numel(class_labels)
    templates(ci, :) = mean(Xtrain(ytrain == class_labels(ci), :), 1);
end
end

%% ------------------------------------------------------------------------
function predictions = predict_from_templates(Xtest, templates, class_labels)
[similarity, valid_tests] = template_similarity(Xtest, templates);
predictions = nan(size(Xtest, 1), 1);
if isempty(similarity)
    return;
end
similarity(~isfinite(similarity)) = -Inf;
[~, best] = max(similarity, [], 2);
predictions(valid_tests) = class_labels(best);
end

%% ------------------------------------------------------------------------
function evidence = evidence_from_templates(Xtest, templates, class_labels, ytest)
[similarity, valid_tests] = template_similarity(Xtest, templates);
evidence = nan(size(Xtest, 1), 1);
if isempty(similarity)
    return;
end
[known, correct_idx] = ismember(ytest(valid_tests), class_labels);
finite_scores = all(isfinite(similarity), 2);
usable = known & finite_scores;
if ~any(usable)
    return;
end
similarity = similarity(usable, :);
correct_idx = correct_idx(usable);
linear_idx = sub2ind(size(similarity), (1:size(similarity, 1))', correct_idx);
correct = similarity(linear_idx);
incorrect_mean = (sum(similarity, 2) - correct) ./ (numel(class_labels) - 1);
valid_indices = find(valid_tests);
evidence(valid_indices(usable)) = correct - incorrect_mean;
end

%% ------------------------------------------------------------------------
function [similarity, valid_tests] = template_similarity(Xtest, templates)
templates = templates - mean(templates, 2);
Xtest = Xtest - mean(Xtest, 2);
template_norm = sqrt(sum(templates.^2, 2));
test_norm = sqrt(sum(Xtest.^2, 2));
valid_templates = isfinite(template_norm) & template_norm > eps;
valid_tests = isfinite(test_norm) & test_norm > eps;
if ~all(valid_templates)
    similarity = [];
    return;
end
templates = templates ./ template_norm;
Xtest(valid_tests, :) = Xtest(valid_tests, :) ./ test_norm(valid_tests);
similarity = Xtest(valid_tests, :) * templates';
end

%% ------------------------------------------------------------------------
function value = calculate_balanced_accuracy(y, predictions, nClasses)
recall = nan(nClasses, 1);
for ci = 1:nClasses
    class_mask = y == ci;
    if any(class_mask)
        recall(ci) = mean(predictions(class_mask) == ci);
    end
end
value = mean(recall, 'omitnan');
end

%% ------------------------------------------------------------------------
function y_perm = permute_within_blocks(y, block_ids)
y_perm = y;
blocks = unique(block_ids, 'stable');
for bi = 1:numel(blocks)
    idx = find(block_ids == blocks(bi));
    if numel(idx) > 1
        y_perm(idx) = y(idx(randperm(numel(idx))));
    end
end
end

%% ------------------------------------------------------------------------
function q = bh_fdr_qvalues(pvals)
pvals = pvals(:);
q = nan(size(pvals));
valid = isfinite(pvals);
p = pvals(valid);
if isempty(p)
    return;
end
[p_sorted, order] = sort(p, 'ascend');
m = numel(p_sorted);
q_sorted = p_sorted .* m ./ (1:m)';
q_sorted = flipud(cummin(flipud(q_sorted)));
q_sorted = min(q_sorted, 1);
q_valid = nan(m, 1);
q_valid(order) = q_sorted;
q(valid) = q_valid;
end

%% ------------------------------------------------------------------------
function output = cells_to_matrix(values, nRows, nColumns)
output = nan(nRows, nColumns);
for row = 1:nRows
    if ~isempty(values{row})
        assert(numel(values{row}) == nColumns, 'Detailed result has an unexpected length.');
        output(row, :) = values{row}(:)';
    end
end
end

%% ------------------------------------------------------------------------
function output = cells_to_confusions(values, nROIs, nClasses)
output = nan(nROIs, nClasses, nClasses);
for roi_idx = 1:nROIs
    if ~isempty(values{roi_idx})
        output(roi_idx, :, :) = reshape(values{roi_idx}, [1, nClasses, nClasses]);
    end
end
end

%% ------------------------------------------------------------------------
function normalized = normalize_confusions(confusions)
normalized = nan(size(confusions));
for roi_idx = 1:size(confusions, 1)
    matrix = squeeze(confusions(roi_idx, :, :));
    if all(isfinite(matrix(:)))
        row_totals = sum(matrix, 2);
        valid_rows = row_totals > 0;
        matrix(valid_rows, :) = matrix(valid_rows, :) ./ row_totals(valid_rows);
        normalized(roi_idx, :, :) = reshape(matrix, [1, size(matrix, 1), size(matrix, 2)]);
    end
end
end

%% ------------------------------------------------------------------------
function percentile = row_percentile(values, requested_percentile)
percentile = nan(size(values, 1), 1);
for row = 1:size(values, 1)
    valid = values(row, isfinite(values(row, :)));
    if ~isempty(valid)
        percentile(row) = prctile(valid, requested_percentile);
    end
end
end
