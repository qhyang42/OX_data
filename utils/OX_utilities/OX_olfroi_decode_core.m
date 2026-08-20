function results = OX_olfroi_decode_core(subjidx, target, varargin)
%OX_OLFROI_DECODE_CORE Focused olfactory ROI decoding for the OX project.
%
%   results = OX_olfroi_decode_core(subjidx, target, Name, Value, ...)
%
% This function is intentionally independent of OX_roi_decode_core. It
% decodes either odor identity or context in bilateral primary olfactory
% cortex and olfactory-amygdala ROIs. The primary decoder uses nested,
% run-blocked cross-validation to choose between a regularized linear LDA
% and a linear one-vs-one ECOC SVM. A correlation-template decoder is saved
% as a leave-one-run-out benchmark.
%
% Trials are assumed to be ordered as 80 consecutive runs with an equal
% number of trials per run. The primary analysis uses every finite voxel in
% the intersection of each ROI and the functional gray-matter mask. Low
% gray-matter overlap is flagged but is not an exclusion criterion.
%
% Required inputs:
%   subjidx                 Numeric subject index (2--6 for final data).
%   target                  'context' or 'odor'.
%
% Important name-value options:
%   'MRIRoot'               MRI root; auto-detected when empty.
%   'ROISelection'          'old' (default), 'primary', 'secondary', or
%                           'all' (primary + secondary).
%   'NumRuns'               Independent consecutive runs (default 80).
%   'NumOuterFolds'         Run-blocked outer folds (default 10).
%   'NumInnerFolds'         Run-blocked inner folds (default 5).
%   'SVMBoxConstraints'     SVM C candidates (default 10.^(-3:2)).
%   'LDAGammas'             LDA covariance shrinkage candidates.
%   'MinVoxels'             Minimum finite GM-overlap features (10).
%   'LowGMOverlapThreshold' Flag overlap fractions below this value (0.5).
%   'ROINames'              Optional subset of canonical output ROI names.
%   'NumPermutations'       Within-run label permutations (1000).
%   'RunTemplateBackup'     Save LORO template results (true).
%   'UseParallel'           Parallelize across ROIs (true).
%   'RandomSeed'            Fold and permutation seed (1).
%   'OutputDir'             Output directory; generated when empty.
%   'SaveOutputs'           Save results MAT file (true).
%
% Legacy hierarchy used when ROISelection='old':
%   olf_primary_bilateral   union(AON, TU, pirF, pirT)
%   olf_amygdala_bilateral union(ACo, MeA, PAC, PCo)
% followed by the eight bilateral component masks.

%% Parse inputs
p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'subjidx', @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x == round(x));
addRequired(p, 'target', @(x) any(strcmpi(string(x), ["context", "odor"])));
addParameter(p, 'MRIRoot', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'ROISelection', 'old', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'NumRuns', 80, @(x) isnumeric(x) && isscalar(x) && x >= 2 && x == round(x));
addParameter(p, 'NumOuterFolds', 10, @(x) isnumeric(x) && isscalar(x) && x >= 2 && x == round(x));
addParameter(p, 'NumInnerFolds', 5, @(x) isnumeric(x) && isscalar(x) && x >= 2 && x == round(x));
addParameter(p, 'SVMBoxConstraints', 10.^(-3:2), @is_positive_vector);
addParameter(p, 'LDAGammas', [0, 0.25, 0.5, 0.75, 1], @is_lda_gamma_vector);
addParameter(p, 'DemeanPatterns', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ScaleVoxels', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'MinVoxels', 10, @(x) isnumeric(x) && isscalar(x) && x >= 2 && x == round(x));
addParameter(p, 'LowGMOverlapThreshold', 0.5, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 1);
addParameter(p, 'ROINames', strings(0, 1), @is_text_vector);
addParameter(p, 'NumPermutations', 1000, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0 && x == round(x));
addParameter(p, 'RunTemplateBackup', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'UseParallel', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'RandomSeed', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'OutputDir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'SaveOutputs', true, @(x) islogical(x) && isscalar(x));
parse(p, subjidx, target, varargin{:});
opts = p.Results;
opts.SVMBoxConstraints = unique(double(opts.SVMBoxConstraints(:)'), 'stable');
opts.LDAGammas = unique(double(opts.LDAGammas(:)'), 'stable');
opts.ROINames = string(opts.ROINames);
opts.ROINames = opts.ROINames(:);
opts.ROISelection = OX_normalize_decoding_roi_selection(opts.ROISelection);
target = lower(char(string(target)));

assert(~isempty(opts.SVMBoxConstraints) || ~isempty(opts.LDAGammas), ...
    'At least one SVMBoxConstraints or LDAGammas candidate is required.');
assert(opts.NumOuterFolds <= opts.NumRuns, ...
    'NumOuterFolds cannot exceed NumRuns.');
assert(opts.NumInnerFolds < opts.NumRuns, ...
    'NumInnerFolds must be smaller than NumRuns.');

%% Resolve subject paths and labels
analysis_tic = tic;
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

mriroot = resolve_mri_root(opts.MRIRoot);
mridatapath = fullfile(mriroot, subjname, 'nifti');
base_outdir = fullfile(mridatapath, 'single_trial_by_category');
fit_file = fullfile(base_outdir, 'TYPED_FITHRF_GLMDENOISE_RR.mat');
gm_mask_file = fullfile(mridatapath, 'coreg', 'gm_mask_thr05_func.nii');
[roi_selection, roi_dirs] = OX_resolve_decoding_roi_selection( ...
    mridatapath, opts.ROISelection, '');
assert(isfile(fit_file), 'Missing GLMsingle file: %s', fit_file);
assert(isfile(gm_mask_file), 'Missing gray-matter mask: %s', gm_mask_file);
assert(exist('fitcecoc', 'file') == 2 && exist('templateSVM', 'file') == 2, ...
    'Nested SVM decoding requires Statistics and Machine Learning Toolbox.');
assert(exist('fitcdiscr', 'file') == 2, ...
    'Nested LDA decoding requires Statistics and Machine Learning Toolbox.');

fprintf('\n[%s | OPTIMIZED OLF ROI] Loading %s\n', upper(target), fit_file);
S = load(fit_file, 'modelmd');
assert(isfield(S, 'modelmd'), 'GLMsingle file must contain modelmd.');
modelmd = squeeze(S.modelmd);
clear S

[odor_raw, category_raw] = OX_get_odor(subjname);
[y, class_values] = encode_target(target, odor_raw, category_raw);
nTrials = numel(y);
nClasses = numel(unique(y));
chance = 1 / nClasses;
assert(ismatrix(modelmd) && size(modelmd, 2) == nTrials, ...
    'modelmd must have one column per trial (%d expected).', nTrials);
assert(mod(nTrials, opts.NumRuns) == 0, ...
    '%d trials cannot be divided into %d equal consecutive runs.', nTrials, opts.NumRuns);
trials_per_run = nTrials / opts.NumRuns;
run_ids = repelem((1:opts.NumRuns)', trials_per_run);
modelmd = center_model_within_runs(modelmd, run_ids);

rng(opts.RandomSeed, 'twister');
outer_fold_ids = make_valid_group_folds(run_ids, y, opts.NumOuterFolds, ...
    nClasses, opts.RandomSeed);
inner_fold_ids = make_inner_fold_matrix(run_ids, y, outer_fold_ids, ...
    opts.NumInnerFolds, nClasses, opts.RandomSeed + 1000);

fprintf('Subject: %s (%s) | Trials: %d | Runs: %d x %d trials\n', ...
    subjname, subjname_real, nTrials, opts.NumRuns, trials_per_run);
fprintf('Target: %s | Classes: %d | Chance: %.4f | Outer/inner folds: %d/%d\n', ...
    target, nClasses, chance, opts.NumOuterFolds, opts.NumInnerFolds);
fprintf('ROI selection: %s | Directories: %s\n', ...
    roi_selection, strjoin(string(roi_dirs), ', '));
print_class_counts(y, class_values);

%% Load focused ROI feature sets
gm_header = spm_vol(gm_mask_file);
gm_mask = spm_read_vols(gm_header) > 0;
gm_inds = find(gm_mask);
assert(size(modelmd, 1) == numel(gm_inds), ...
    ['Mapping check failed: modelmd has %d rows but find(gm_mask) has %d voxels. ' ...
     'modelmd row i must correspond to gm_inds(i).'], size(modelmd, 1), numel(gm_inds));
model_index_volume = zeros(size(gm_mask), 'uint32');
model_index_volume(gm_inds) = uint32(1:numel(gm_inds));

if strcmp(roi_selection, 'old')
    [roi_defs, component_masks] = define_and_load_rois( ...
        roi_dirs{1}, gm_header, opts.ROINames);
else
    [roi_defs, component_masks] = define_and_load_selected_rois( ...
        roi_dirs, gm_header, opts.ROINames, roi_selection);
end
nROIs = height(roi_defs);
roi_features = cell(nROIs, 1);
n_mask_voxels = zeros(nROIs, 1);
n_gm_overlap = zeros(nROIs, 1);
n_features_used = zeros(nROIs, 1);
gm_overlap_fraction = nan(nROIs, 1);
low_gm_overlap = false(nROIs, 1);
status = repmat("ok", nROIs, 1);

for ri = 1:nROIs
    roi_mask = build_roi_mask(roi_defs.component_keys{ri}, component_masks);
    n_mask_voxels(ri) = nnz(roi_mask);
    overlap_linear = find(roi_mask & gm_mask);
    n_gm_overlap(ri) = numel(overlap_linear);
    if n_mask_voxels(ri) > 0
        gm_overlap_fraction(ri) = n_gm_overlap(ri) / n_mask_voxels(ri);
    end
    low_gm_overlap(ri) = gm_overlap_fraction(ri) < opts.LowGMOverlapThreshold;

    feature_idx = double(model_index_volume(overlap_linear));
    feature_idx = feature_idx(feature_idx > 0);
    if ~isempty(feature_idx)
        feature_idx = feature_idx(all(isfinite(modelmd(feature_idx, :)), 2));
    end
    feature_idx = unique(feature_idx(:)', 'stable');
    roi_features{ri} = feature_idx;
    n_features_used(ri) = numel(feature_idx);
    if n_gm_overlap(ri) < opts.MinVoxels
        status(ri) = "insufficient_gm_overlap";
    elseif n_features_used(ri) < opts.MinVoxels
        status(ri) = "insufficient_finite_features";
    end
end

roi_metadata = table(roi_defs.roi_name, roi_defs.roi_level, roi_defs.roi_family, ...
    roi_defs.source_masks, status, n_mask_voxels, n_gm_overlap, ...
    gm_overlap_fraction, low_gm_overlap, n_features_used, ...
    'VariableNames', {'roi_name', 'roi_level', 'roi_family', 'source_masks', ...
    'status', 'n_mask_voxels', 'n_gm_overlap', 'gm_overlap_fraction', ...
    'low_gm_overlap', 'n_features_used'});
valid_roi = status == "ok";
assert(any(valid_roi), 'No focused ROI contains at least %d finite GM voxels.', opts.MinVoxels);
fprintf('Focused ROIs: %d total, %d valid, %d flagged for low GM overlap (<%.2f).\n', ...
    nROIs, sum(valid_roi), sum(low_gm_overlap), opts.LowGMOverlapThreshold);

%% Set up candidate models and parallel execution
candidates = make_candidate_table(opts.LDAGammas, opts.SVMBoxConstraints);
use_parallel = opts.UseParallel && license('test', 'Distrib_Computing_Toolbox');
if opts.UseParallel && ~use_parallel
    warning('Parallel Computing Toolbox unavailable; running serially.');
end
if use_parallel && isempty(gcp('nocreate'))
    parpool;
end

%% Observed optimized decoding
fprintf('Running nested optimized decoding across %d candidate models...\n', height(candidates));
observed_tic = tic;
optimized_cells = cell(nROIs, 1);
if use_parallel
    model_const = parallel.pool.Constant(modelmd);
    model_cleanup = onCleanup(@() delete(model_const));
    parfor ri = 1:nROIs
        if valid_roi(ri)
            X = double(model_const.Value(roi_features{ri}, :)');
            optimized_cells{ri} = evaluate_optimized_roi(X, y, outer_fold_ids, ...
                inner_fold_ids, candidates, nClasses, opts, true);
        end
    end
else
    for ri = 1:nROIs
        if valid_roi(ri)
            X = double(modelmd(roi_features{ri}, :)');
            optimized_cells{ri} = evaluate_optimized_roi(X, y, outer_fold_ids, ...
                inner_fold_ids, candidates, nClasses, opts, true);
        end
    end
end
fprintf('Observed optimized decoding completed in %.2f minutes.\n', toc(observed_tic) / 60);
optimized = package_observed_results(optimized_cells, nROIs, nTrials, ...
    opts.NumOuterFolds, nClasses, chance);
decoder_failed = valid_roi & ~isfinite(optimized.balanced_accuracy);
status(decoder_failed) = "decoder_failed";
valid_results = isfinite(optimized.balanced_accuracy);
assert(any(valid_results), 'No focused ROI produced a finite optimized result.');
roi_metadata.status = status;

%% Within-run permutation inference for optimized decoder
null_balanced_accuracy = nan(nROIs, opts.NumPermutations);
p_empirical = nan(nROIs, 1);
q_fdr = nan(nROIs, 1);
p_fwe_maxstat = nan(nROIs, 1);
max_null = nan(opts.NumPermutations, 1);

if opts.NumPermutations > 0
    fprintf('Running %d within-run permutations of the complete nested decoder...\n', ...
        opts.NumPermutations);
    permutation_tic = tic;
    rng(opts.RandomSeed + 2000, 'twister');
    for pi = 1:opts.NumPermutations
        y_perm = permute_within_runs(y, run_ids);
        perm_metric = nan(nROIs, 1);
        if use_parallel
            parfor ri = 1:nROIs
                if valid_results(ri)
                    X = double(model_const.Value(roi_features{ri}, :)');
                    perm_result = evaluate_optimized_roi(X, y_perm, ...
                        outer_fold_ids, inner_fold_ids, candidates, nClasses, opts, false);
                    perm_metric(ri) = perm_result.balanced_accuracy;
                end
            end
        else
            for ri = 1:nROIs
                if valid_results(ri)
                    X = double(modelmd(roi_features{ri}, :)');
                    perm_result = evaluate_optimized_roi(X, y_perm, ...
                        outer_fold_ids, inner_fold_ids, candidates, nClasses, opts, false);
                    perm_metric(ri) = perm_result.balanced_accuracy;
                end
            end
        end
        null_balanced_accuracy(:, pi) = perm_metric;
        comparable = valid_results & isfinite(perm_metric);
        if any(comparable)
            max_null(pi) = max(perm_metric(comparable) - chance);
        end
        if pi == 1 || mod(pi, 10) == 0 || pi == opts.NumPermutations
            fprintf('  permutation %d/%d | elapsed %.1f min\n', ...
                pi, opts.NumPermutations, toc(permutation_tic) / 60);
        end
    end

    valid_max = max_null(isfinite(max_null));
    for ri = find(valid_results)'
        roi_null = null_balanced_accuracy(ri, :);
        roi_null = roi_null(isfinite(roi_null));
        p_empirical(ri) = (1 + sum(roi_null >= optimized.balanced_accuracy(ri))) / ...
            (numel(roi_null) + 1);
        p_fwe_maxstat(ri) = (1 + sum(valid_max >= ...
            (optimized.balanced_accuracy(ri) - chance))) / (numel(valid_max) + 1);
    end
    q_fdr(valid_results) = bh_fdr_qvalues(p_empirical(valid_results));
end

null_mean = mean(null_balanced_accuracy, 2, 'omitnan');
null_sd = std(null_balanced_accuracy, 0, 2, 'omitnan');
optimized.null_balanced_accuracy = null_balanced_accuracy;
optimized.null_mean = null_mean;
optimized.null_sd = null_sd;
optimized.effect_vs_null_mean = optimized.balanced_accuracy - null_mean;
optimized.p_empirical = p_empirical;
optimized.q_fdr = q_fdr;
optimized.p_fwe_maxstat = p_fwe_maxstat;
optimized.max_null = max_null;

%% Correlation-template backup
template_backup = struct('enabled', opts.RunTemplateBackup);
if opts.RunTemplateBackup
    fprintf('Running leave-one-run-out correlation-template backup...\n');
    template_cells = cell(nROIs, 1);
    demean_patterns = opts.DemeanPatterns;
    if use_parallel
        parfor ri = 1:nROIs
            if valid_roi(ri)
                X = double(model_const.Value(roi_features{ri}, :)');
                template_cells{ri} = evaluate_template_roi(X, y, run_ids, ...
                    nClasses, demean_patterns);
            end
        end
    else
        for ri = 1:nROIs
            if valid_roi(ri)
                X = double(modelmd(roi_features{ri}, :)');
                template_cells{ri} = evaluate_template_roi(X, y, run_ids, ...
                    nClasses, demean_patterns);
            end
        end
    end
    template_backup = package_template_results(template_cells, nROIs, ...
        nTrials, opts.NumRuns, nClasses, chance);
    template_backup.enabled = true;
    template_backup.cross_validation = 'leave-one-run-out';
    template_backup.scale_voxels = false;
end
if use_parallel
    clear model_cleanup
end

%% Package and save
summary = roi_metadata;
summary.optimized_accuracy = optimized.accuracy;
summary.optimized_balanced_accuracy = optimized.balanced_accuracy;
summary.optimized_balanced_accuracy_minus_chance = ...
    optimized.balanced_accuracy_minus_chance;
summary.optimized_null_mean = optimized.null_mean;
summary.optimized_effect_vs_null_mean = optimized.effect_vs_null_mean;
summary.p_empirical = optimized.p_empirical;
summary.q_fdr = optimized.q_fdr;
summary.p_fwe_maxstat = optimized.p_fwe_maxstat;
if opts.RunTemplateBackup
    summary.template_accuracy = template_backup.accuracy;
    summary.template_balanced_accuracy = template_backup.balanced_accuracy;
    summary.template_balanced_accuracy_minus_chance = ...
        template_backup.balanced_accuracy_minus_chance;
end

if strlength(string(opts.OutputDir)) == 0
    if strcmp(roi_selection, 'old')
        output_name = sprintf('olfroi_decoding_%s_nested_kfold%d', ...
            target, opts.NumOuterFolds);
    else
        output_name = sprintf('olfroi_decoding_%s_%s_nested_kfold%d', ...
            roi_selection, target, opts.NumOuterFolds);
    end
    output_dir = fullfile(base_outdir, output_name);
else
    output_dir = char(string(opts.OutputDir));
end

results = struct();
results.subject = struct('index', subjidx, 'name', subjname, 'name_real', subjname_real);
results.analysis = struct('target', target, 'roi_selection', roi_selection, ...
    'primary_decoder', 'nested_model_selection', ...
    'metric', 'balanced_accuracy', 'chance', chance, ...
    'run_structure', sprintf('%d consecutive independent runs x %d trials', ...
    opts.NumRuns, trials_per_run));
results.classes = struct('values', {class_values}, ...
    'counts', accumarray(y, 1, [nClasses, 1]), 'chance', chance);
results.roi_metadata = roi_metadata;
results.summary = summary;
results.true_labels = y;
results.run_ids = run_ids;
results.outer_fold_ids = outer_fold_ids;
results.inner_fold_ids = inner_fold_ids;
results.candidates = candidates;
results.optimized = optimized;
results.template_backup = template_backup;
results.options = opts;
results.permutation = ['labels shuffled independently within each of 80 runs; ' ...
    'nested model and hyperparameter selection rerun'];
results.individual_inference_family = ...
    'all valid focused bilateral composite and component ROIs in this result';
results.output_dir = output_dir;
results.runtime_minutes = toc(analysis_tic) / 60;
results.matlab_version = version;

if opts.SaveOutputs
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end
    output_file = fullfile(output_dir, ...
        sprintf('%s_olfroi_nested_subj%d_results.mat', target, subjidx));
    results.output_file = output_file;
    save(output_file, 'results', '-v7.3');
    fprintf('Saved individual olfactory ROI results to %s\n', output_file);
else
    results.output_file = '';
end
end

%% ------------------------------------------------------------------------
function tf = is_positive_vector(x)
tf = isnumeric(x) && isvector(x) && all(isfinite(x)) && all(x > 0);
end

%% ------------------------------------------------------------------------
function tf = is_lda_gamma_vector(x)
tf = isnumeric(x) && isvector(x) && all(isfinite(x)) && all(x >= 0) && all(x <= 1);
end

%% ------------------------------------------------------------------------
function tf = is_text_vector(x)
tf = isempty(x) || ischar(x) || isstring(x) || iscellstr(x);
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
assert(~isempty(mriroot), 'Could not auto-detect MRIRoot. Pass MRIRoot explicitly.');
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
function modelmd = center_model_within_runs(modelmd, run_ids)
runs = unique(run_ids, 'stable');
for ri = 1:numel(runs)
    mask = run_ids == runs(ri);
    modelmd(:, mask) = modelmd(:, mask) - mean(modelmd(:, mask), 2);
end
end

%% ------------------------------------------------------------------------
function fold_ids = make_valid_group_folds(run_ids, y, nFolds, nClasses, seed)
runs = unique(run_ids, 'stable');
assert(nFolds <= numel(runs), 'Cannot make %d folds from %d runs.', nFolds, numel(runs));
stream = RandStream('mt19937ar', 'Seed', seed);
for attempt = 1:10000
    order = runs(randperm(stream, numel(runs)));
    fold_by_position = mod((0:numel(runs)-1)', nFolds) + 1;
    run_fold = zeros(numel(runs), 1);
    [~, order_position] = ismember(order, runs);
    run_fold(order_position) = fold_by_position;
    [known, run_position] = ismember(run_ids, runs);
    assert(all(known), 'Failed to map trials to runs.');
    candidate = run_fold(run_position);
    if folds_have_all_classes(candidate, y, nFolds, nClasses)
        fold_ids = double(candidate(:));
        return;
    end
end
error('Could not create %d run-blocked folds with all %d classes represented.', ...
    nFolds, nClasses);
end

%% ------------------------------------------------------------------------
function tf = folds_have_all_classes(fold_ids, y, nFolds, nClasses)
tf = true;
for fi = 1:nFolds
    if numel(unique(y(fold_ids == fi))) < nClasses || ...
            numel(unique(y(fold_ids ~= fi))) < nClasses
        tf = false;
        return;
    end
end
end

%% ------------------------------------------------------------------------
function inner_matrix = make_inner_fold_matrix(run_ids, y, outer_ids, ...
        nInnerFolds, nClasses, seed)
nTrials = numel(y);
nOuterFolds = numel(unique(outer_ids));
inner_matrix = nan(nTrials, nOuterFolds);
for oi = 1:nOuterFolds
    train_mask = outer_ids ~= oi;
    train_runs = run_ids(train_mask);
    train_y = y(train_mask);
    inner_matrix(train_mask, oi) = make_valid_group_folds(train_runs, train_y, ...
        nInnerFolds, nClasses, seed + oi);
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
function [defs, component_masks] = define_and_load_rois(roi_dir, ref_header, requested)
keys = ["AON", "TU", "pirF", "pirT", "ACo", "MeA", "PAC", "PCo"];
stems = ["olf_AON_bilateral", "olf_TU_bilateral", ...
    "olf_pirF_bilateral", "olf_pirT_bilateral", ...
    "amygsub_ACo_bilateral", "amygsub_MeA_bilateral", ...
    "amygsub_PAC_bilateral", "amygsub_PCo_bilateral"];
component_masks = struct();
for i = 1:numel(keys)
    mask_file = find_nifti_by_stem(roi_dir, stems(i) + "_func_thr02");
    header = spm_vol(mask_file);
    assert_same_geometry(header, ref_header, mask_file);
    component_masks.(char(keys(i))) = spm_read_vols(header) > 0;
end

roi_name = ["olf_primary_bilateral"; "olf_amygdala_bilateral"; ...
    "olf_AON_bilateral"; "olf_TU_bilateral"; "olf_pirF_bilateral"; ...
    "olf_pirT_bilateral"; "amygsub_ACo_bilateral"; ...
    "amygsub_MeA_bilateral"; "amygsub_PAC_bilateral"; ...
    "amygsub_PCo_bilateral"];
roi_level = ["composite"; "composite"; repmat("subregion", 8, 1)];
roi_family = ["primary_olfactory"; "olfactory_amygdala"; ...
    repmat("primary_olfactory", 4, 1); repmat("olfactory_amygdala", 4, 1)];
source_masks = ["AON,TU,pirF,pirT"; "ACo,MeA,PAC,PCo"; ...
    "AON"; "TU"; "pirF"; "pirT"; "ACo"; "MeA"; "PAC"; "PCo"];
component_keys = {{'AON', 'TU', 'pirF', 'pirT'}; {'ACo', 'MeA', 'PAC', 'PCo'}; ...
    {'AON'}; {'TU'}; {'pirF'}; {'pirT'}; {'ACo'}; {'MeA'}; {'PAC'}; {'PCo'}};
defs = table(roi_name, roi_level, roi_family, source_masks, component_keys);

if ~isempty(requested)
    [known, position] = ismember(lower(requested), lower(defs.roi_name));
    assert(all(known), 'Unknown ROINames: %s', strjoin(requested(~known), ', '));
    assert(numel(unique(position)) == numel(position), 'ROINames contains duplicates.');
    defs = defs(position, :);
end
end

%% ------------------------------------------------------------------------
function [defs, component_masks] = define_and_load_selected_rois( ...
        roi_dirs, ref_header, requested, selection)
[mask_files, stems] = OX_discover_decoding_roi_files(roi_dirs, true, Inf);
n_masks = numel(mask_files);
component_masks = struct();
roi_name = strings(n_masks, 1);
roi_level = repmat("selected_set", n_masks, 1);
roi_family = repmat(string(selection), n_masks, 1);
source_masks = strings(n_masks, 1);
component_keys = cell(n_masks, 1);

for index = 1:n_masks
    key = sprintf('mask_%03d', index);
    header = spm_vol(mask_files{index});
    assert_same_geometry(header, ref_header, mask_files{index});
    component_masks.(key) = spm_read_vols(header) > 0;
    roi_name(index) = erase(string(stems{index}), "_func_thr02");
    source_masks(index) = string(mask_files{index});
    component_keys{index} = {key};
end
defs = table(roi_name, roi_level, roi_family, source_masks, component_keys);

if ~isempty(requested)
    [known, position] = ismember(lower(requested), lower(defs.roi_name));
    assert(all(known), 'Unknown ROINames: %s', strjoin(requested(~known), ', '));
    assert(numel(unique(position)) == numel(position), 'ROINames contains duplicates.');
    defs = defs(position, :);
end
end

%% ------------------------------------------------------------------------
function mask_file = find_nifti_by_stem(folder, stem)
candidates = {fullfile(folder, char(stem + ".nii")), ...
              fullfile(folder, char(stem + ".nii.gz"))};
exists = cellfun(@isfile, candidates);
assert(sum(exists) == 1, 'Expected exactly one NIfTI for %s; found %d.', stem, sum(exists));
mask_file = candidates{find(exists, 1)};
end

%% ------------------------------------------------------------------------
function assert_same_geometry(header, reference, image_file)
assert(isequal(header.dim, reference.dim), ...
    'Dimension mismatch between %s and gray-matter mask.', image_file);
assert(max(abs(header.mat(:) - reference.mat(:))) < 1e-4, ...
    'Affine mismatch between %s and gray-matter mask.', image_file);
end

%% ------------------------------------------------------------------------
function roi_mask = build_roi_mask(component_keys, component_masks)
first = component_masks.(component_keys{1});
roi_mask = false(size(first));
for i = 1:numel(component_keys)
    roi_mask = roi_mask | component_masks.(component_keys{i});
end
end

%% ------------------------------------------------------------------------
function candidates = make_candidate_table(lda_gammas, svm_cs)
model_type = strings(0, 1);
parameter_name = strings(0, 1);
parameter_value = zeros(0, 1);
% Put LDA first and strongest regularization first so exact ties favor it.
for gamma = sort(lda_gammas, 'descend')
    model_type(end+1, 1) = "lda"; %#ok<AGROW>
    parameter_name(end+1, 1) = "Gamma"; %#ok<AGROW>
    parameter_value(end+1, 1) = gamma; %#ok<AGROW>
end
for c = sort(svm_cs, 'ascend')
    model_type(end+1, 1) = "svm"; %#ok<AGROW>
    parameter_name(end+1, 1) = "BoxConstraint"; %#ok<AGROW>
    parameter_value(end+1, 1) = c; %#ok<AGROW>
end
candidates = table(model_type, parameter_name, parameter_value);
end

%% ------------------------------------------------------------------------
function result = evaluate_optimized_roi(X, y, outer_ids, inner_matrix, ...
        candidates, nClasses, opts, collect_details)
if opts.DemeanPatterns
    X = X - mean(X, 2);
end
nTrials = numel(y);
nOuter = numel(unique(outer_ids));
predictions = nan(nTrials, 1);
fold_balanced_accuracy = nan(nOuter, 1);
selected_model = strings(nOuter, 1);
selected_parameter = nan(nOuter, 1);
selected_inner_balanced_accuracy = nan(nOuter, 1);

for oi = 1:nOuter
    outer_train = outer_ids ~= oi;
    outer_test = outer_ids == oi;
    inner_ids = inner_matrix(:, oi);
    candidate_scores = nan(height(candidates), 1);

    for ci = 1:height(candidates)
        inner_predictions = nan(sum(outer_train), 1);
        train_positions = find(outer_train);
        failed = false;
        for ii = 1:opts.NumInnerFolds
            inner_train_global = outer_train & inner_ids ~= ii;
            inner_valid_global = outer_train & inner_ids == ii;
            [Xtrain, Xvalid, usable] = preprocess_fold( ...
                X(inner_train_global, :), X(inner_valid_global, :), ...
                opts.ScaleVoxels, opts.MinVoxels);
            if ~usable
                failed = true;
                break;
            end
            pred = fit_predict_candidate(Xtrain, y(inner_train_global), Xvalid, ...
                candidates(ci, :));
            if any(~isfinite(pred))
                failed = true;
                break;
            end
            [~, local_positions] = ismember(find(inner_valid_global), train_positions);
            inner_predictions(local_positions) = pred;
        end
        if ~failed && all(isfinite(inner_predictions))
            candidate_scores(ci) = balanced_accuracy_score( ...
                y(outer_train), inner_predictions, nClasses);
        end
    end

    finite_candidates = find(isfinite(candidate_scores));
    if isempty(finite_candidates)
        continue;
    end
    best_value = max(candidate_scores(finite_candidates));
    best_idx = find(candidate_scores >= best_value - 1e-12, 1, 'first');
    [Xtrain, Xtest, usable] = preprocess_fold(X(outer_train, :), X(outer_test, :), ...
        opts.ScaleVoxels, opts.MinVoxels);
    if ~usable
        continue;
    end
    pred = fit_predict_candidate(Xtrain, y(outer_train), Xtest, candidates(best_idx, :));
    if any(~isfinite(pred))
        continue;
    end
    predictions(outer_test) = pred;
    fold_balanced_accuracy(oi) = balanced_accuracy_score(y(outer_test), pred, nClasses);
    selected_model(oi) = candidates.model_type(best_idx);
    selected_parameter(oi) = candidates.parameter_value(best_idx);
    selected_inner_balanced_accuracy(oi) = candidate_scores(best_idx);
end

result = empty_observed_result();
if all(isfinite(predictions))
    result.accuracy = mean(predictions == y);
    result.balanced_accuracy = balanced_accuracy_score(y, predictions, nClasses);
end
if collect_details
    result.predictions = predictions;
    result.fold_balanced_accuracy = fold_balanced_accuracy;
    result.selected_model = selected_model;
    result.selected_parameter = selected_parameter;
    result.selected_inner_balanced_accuracy = selected_inner_balanced_accuracy;
    if all(isfinite(predictions))
        result.confusion_matrix = accumarray([y, predictions], 1, [nClasses, nClasses]);
    end
end
end

%% ------------------------------------------------------------------------
function result = empty_observed_result()
result = struct('accuracy', NaN, 'balanced_accuracy', NaN, 'predictions', [], ...
    'fold_balanced_accuracy', [], 'selected_model', strings(0, 1), ...
    'selected_parameter', [], 'selected_inner_balanced_accuracy', [], ...
    'confusion_matrix', []);
end

%% ------------------------------------------------------------------------
function [Xtrain, Xtest, usable] = preprocess_fold(Xtrain, Xtest, scale_voxels, min_voxels)
finite_features = all(isfinite(Xtrain), 1) & all(isfinite(Xtest), 1);
Xtrain = Xtrain(:, finite_features);
Xtest = Xtest(:, finite_features);
if scale_voxels
    mu = mean(Xtrain, 1);
    sigma = std(Xtrain, 0, 1);
    variable = isfinite(sigma) & sigma > sqrt(eps);
    Xtrain = Xtrain(:, variable);
    Xtest = Xtest(:, variable);
    mu = mu(variable);
    sigma = sigma(variable);
    Xtrain = (Xtrain - mu) ./ sigma;
    Xtest = (Xtest - mu) ./ sigma;
end
usable = size(Xtrain, 2) >= min_voxels && all(isfinite(Xtrain), 'all') && ...
    all(isfinite(Xtest), 'all');
end

%% ------------------------------------------------------------------------
function predictions = fit_predict_candidate(Xtrain, ytrain, Xtest, candidate)
persistent warned_messages
predictions = nan(size(Xtest, 1), 1);
try
    if candidate.model_type == "lda"
        model = fitcdiscr(Xtrain, ytrain, 'DiscrimType', 'linear', ...
            'Gamma', candidate.parameter_value, 'Delta', 0, ...
            'Prior', 'uniform', 'FillCoeffs', 'off');
    else
        learner = templateSVM('KernelFunction', 'linear', ...
            'BoxConstraint', candidate.parameter_value, 'Standardize', false);
        model = fitcecoc(Xtrain, ytrain, 'Learners', learner, ...
            'Coding', 'onevsone', 'Prior', 'uniform', 'FitPosterior', false);
    end
    predictions = double(predict(model, Xtest));
catch err
    if isempty(warned_messages)
        warned_messages = strings(0, 1);
    end
    signature = string(err.identifier) + ":" + string(err.message);
    if ~any(warned_messages == signature) && numel(warned_messages) < 10
        warning('OX_olfroi_decode_core:CandidateFailed', ...
            'Candidate model failed and was skipped: %s', err.message);
        warned_messages(end+1, 1) = signature;
    end
end
end

%% ------------------------------------------------------------------------
function score = balanced_accuracy_score(ytrue, ypred, nClasses)
recall = nan(nClasses, 1);
for ci = 1:nClasses
    class_mask = ytrue == ci;
    if any(class_mask)
        recall(ci) = mean(ypred(class_mask) == ci);
    end
end
score = mean(recall, 'omitnan');
end

%% ------------------------------------------------------------------------
function optimized = package_observed_results(cells, nROIs, nTrials, nFolds, nClasses, chance)
optimized = struct();
optimized.accuracy = nan(nROIs, 1);
optimized.balanced_accuracy = nan(nROIs, 1);
optimized.predictions = nan(nROIs, nTrials);
optimized.fold_balanced_accuracy = nan(nROIs, nFolds);
optimized.selected_model = strings(nROIs, nFolds);
optimized.selected_parameter = nan(nROIs, nFolds);
optimized.selected_inner_balanced_accuracy = nan(nROIs, nFolds);
optimized.confusion_matrices = nan(nROIs, nClasses, nClasses);
for ri = 1:nROIs
    if isempty(cells{ri})
        continue;
    end
    value = cells{ri};
    optimized.accuracy(ri) = value.accuracy;
    optimized.balanced_accuracy(ri) = value.balanced_accuracy;
    optimized.predictions(ri, :) = value.predictions(:)';
    optimized.fold_balanced_accuracy(ri, :) = value.fold_balanced_accuracy(:)';
    optimized.selected_model(ri, :) = value.selected_model(:)';
    optimized.selected_parameter(ri, :) = value.selected_parameter(:)';
    optimized.selected_inner_balanced_accuracy(ri, :) = ...
        value.selected_inner_balanced_accuracy(:)';
    if ~isempty(value.confusion_matrix)
        optimized.confusion_matrices(ri, :, :) = value.confusion_matrix;
    end
end
optimized.accuracy_minus_chance = optimized.accuracy - chance;
optimized.balanced_accuracy_minus_chance = optimized.balanced_accuracy - chance;
end

%% ------------------------------------------------------------------------
function result = evaluate_template_roi(X, y, run_ids, nClasses, demean_patterns)
if demean_patterns
    X = X - mean(X, 2);
end
runs = unique(run_ids, 'stable');
predictions = nan(numel(y), 1);
fold_accuracy = nan(numel(runs), 1);
for ri = 1:numel(runs)
    train = run_ids ~= runs(ri);
    test = ~train;
    templates = nan(nClasses, size(X, 2));
    for ci = 1:nClasses
        templates(ci, :) = mean(X(train & y == ci, :), 1);
    end
    similarity = row_correlation(X(test, :), templates);
    [~, predictions(test)] = max(similarity, [], 2);
    fold_accuracy(ri) = mean(predictions(test) == y(test));
end
result = struct('accuracy', NaN, 'balanced_accuracy', NaN, ...
    'predictions', predictions, 'fold_accuracy', fold_accuracy, 'confusion_matrix', []);
if all(isfinite(predictions))
    result.accuracy = mean(predictions == y);
    result.balanced_accuracy = balanced_accuracy_score(y, predictions, nClasses);
    result.confusion_matrix = accumarray([y, predictions], 1, [nClasses, nClasses]);
end
end

%% ------------------------------------------------------------------------
function similarity = row_correlation(test_patterns, templates)
test_patterns = test_patterns - mean(test_patterns, 2);
templates = templates - mean(templates, 2);
test_norm = sqrt(sum(test_patterns.^2, 2));
template_norm = sqrt(sum(templates.^2, 2));
denominator = test_norm * template_norm';
similarity = (test_patterns * templates') ./ denominator;
similarity(~isfinite(similarity)) = -Inf;
end

%% ------------------------------------------------------------------------
function output = package_template_results(cells, nROIs, nTrials, nFolds, nClasses, chance)
output = struct();
output.accuracy = nan(nROIs, 1);
output.balanced_accuracy = nan(nROIs, 1);
output.predictions = nan(nROIs, nTrials);
output.fold_accuracy = nan(nROIs, nFolds);
output.confusion_matrices = nan(nROIs, nClasses, nClasses);
for ri = 1:nROIs
    if isempty(cells{ri})
        continue;
    end
    value = cells{ri};
    output.accuracy(ri) = value.accuracy;
    output.balanced_accuracy(ri) = value.balanced_accuracy;
    output.predictions(ri, :) = value.predictions(:)';
    output.fold_accuracy(ri, :) = value.fold_accuracy(:)';
    if ~isempty(value.confusion_matrix)
        output.confusion_matrices(ri, :, :) = value.confusion_matrix;
    end
end
output.accuracy_minus_chance = output.accuracy - chance;
output.balanced_accuracy_minus_chance = output.balanced_accuracy - chance;
end

%% ------------------------------------------------------------------------
function y_perm = permute_within_runs(y, run_ids)
y_perm = y;
runs = unique(run_ids, 'stable');
for ri = 1:numel(runs)
    positions = find(run_ids == runs(ri));
    y_perm(positions) = y(positions(randperm(numel(positions))));
end
end

%% ------------------------------------------------------------------------
function q = bh_fdr_qvalues(p)
p = p(:);
q = nan(size(p));
valid = isfinite(p);
pv = p(valid);
if isempty(pv)
    return;
end
[sorted_p, order] = sort(pv);
m = numel(sorted_p);
adjusted = sorted_p .* m ./ (1:m)';
adjusted = flipud(cummin(flipud(adjusted)));
adjusted = min(adjusted, 1);
restored = nan(m, 1);
restored(order) = adjusted;
q(valid) = restored;
end
