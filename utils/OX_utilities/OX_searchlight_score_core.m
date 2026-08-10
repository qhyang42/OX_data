function results = OX_searchlight_score_core(subjidx, target, varargin)
%OX_SEARCHLIGHT_SCORE_CORE Shared evidence-scoring implementation for context/odor searchlights.
%
% This function uses the same inputs, preprocessing, and cross-validation
% conventions as OX_searchlight_decode_core.
%
% Important default assumptions:
%   1. modelmd is [gray-matter voxels x trials].
%   2. modelmd rows follow find(gm_mask) exactly.
%   3. If RunLabels is omitted, trials are ordered as contiguous, equally
%      sized runs. For 800 trials and NumRuns=80, trials 1:10 are run 1,
%      trials 11:20 are run 2, and so on.
%   4. Each voxel is mean-centered across trials within each run before
%      decoding.
%   5. R2Threshold is in the same units as GLMsingle's R2 variable.
%
% Name-value options:
%   'MRIRoot'             MRI root directory. Auto-detected if empty.
%   'NumRuns'             Number of runs used to infer run labels (80).
%   'RunLabels'           Exact [nTrials x 1] run labels; preferred when
%                         available.
%   'CrossValidation'     'leave-one-run-out' (default) or 'leave-one-out'.
%   'R2Threshold'         Include voxels with R2 > threshold (0.5).
%   'SearchlightRadius'   Radius in voxels (2).
%   'RestrictFeaturesToR2' Restrict both centers and neighborhood features
%                         to R2-selected voxels (true). Set false to match
%                         the old script, which restricted centers only.
%   'MinVoxels'           Minimum usable features per searchlight (10).
%   'MaxCenters'          Optional deterministic subset for smoke tests
%                         (Inf = all centers). Do not use for final maps.
%   'NumPermutations'     Number of full-map permutations (1000).
%   'RandomSeed'          RNG seed (1).
%   'UseParallel'         Use parfor when available (true).
%   'DemeanPatterns'      Remove each trial's spatial mean (true). This is
%                         mathematically redundant with Pearson correlation,
%                         but makes the intended preprocessing explicit.
%   'ScaleVoxels'         Z-score each voxel using training-fold statistics
%                         only (true).
%   'OutputDir'           Output directory; auto-generated if empty.
%   'SaveOutputs'         Write NIfTI and MAT outputs (true).
%
% For each test trial, evidence is its Pearson similarity to the correct
% training-class template minus its mean similarity to all incorrect
% templates. Output NIfTIs include mean evidence, uncorrected p, FDR q,
% max-stat FWE p, and q<.05 / FWE<.05 masks.

%% Parse inputs
p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'subjidx', @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x == round(x));
addRequired(p, 'target', @(x) any(strcmpi(string(x), ["context", "odor"])));
addParameter(p, 'MRIRoot', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'NumRuns', 80, @(x) isnumeric(x) && isscalar(x) && x >= 2 && x == round(x));
addParameter(p, 'RunLabels', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
addParameter(p, 'CrossValidation', 'leave-one-run-out', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'R2Threshold', 0.5, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'SearchlightRadius', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'RestrictFeaturesToR2', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'MinVoxels', 10, @(x) isnumeric(x) && isscalar(x) && x >= 2 && x == round(x));
addParameter(p, 'MaxCenters', Inf, @(x) isnumeric(x) && isscalar(x) && (isinf(x) || (x >= 1 && x == round(x))));
addParameter(p, 'NumPermutations', 1000, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x == round(x));
addParameter(p, 'RandomSeed', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'UseParallel', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'DemeanPatterns', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ScaleVoxels', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'OutputDir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'SaveOutputs', true, @(x) islogical(x) && isscalar(x));
parse(p, subjidx, target, varargin{:});
opts = p.Results;
target = lower(char(string(target)));
[cv_method, cv_short_name] = normalize_cv_method(opts.CrossValidation);
opts.CrossValidation = cv_method;

SUBJNAMES = {'240711_fMRI_OX_NWU_AS', ...
             '240723_fMRI_OX_NWU_LS', ...
             '240814_fMRI_OX_NWU_JN', ...
             '240816_fMRI_OX_NWU_RR', ...
             '241018_fMRI_OX_NWU_BN', ...
             '250117_fMRI_OX_NWU_VS'};

if subjidx < 1 || subjidx > numel(SUBJNAMES)
    error('subjidx must be between 1 and %d.', numel(SUBJNAMES));
end
subjname = sprintf('subj_%d', subjidx);
subjname_real = SUBJNAMES{subjidx};

%% Resolve paths and load data
mriroot = resolve_mri_root(opts.MRIRoot);
mridatapath = fullfile(mriroot, subjname, 'nifti');
base_outdir = fullfile(mridatapath, 'single_trial_by_category');
fit_file = fullfile(base_outdir, 'TYPED_FITHRF_GLMDENOISE_RR.mat');
mask_file = fullfile(mridatapath, 'coreg', 'gm_mask_thr05_func.nii');

assert(isfile(fit_file), 'Missing GLMsingle file: %s', fit_file);
assert(isfile(mask_file), 'Missing gray-matter mask: %s', mask_file);

fprintf('\n[%s] Loading %s\n', upper(target), fit_file);
S = load(fit_file, 'modelmd', 'R2');
assert(isfield(S, 'modelmd') && isfield(S, 'R2'), ...
    'GLMsingle file must contain modelmd and R2.');
modelmd = squeeze(S.modelmd);
R2 = S.R2(:);
clear S

[odor_raw, category_raw] = OX_get_odor(subjname);
category_raw = string(category_raw(:));
nTrials = numel(category_raw);

assert(ndims(modelmd) == 2, 'squeeze(modelmd) must produce a 2-D matrix.');
assert(size(modelmd, 2) == nTrials, ...
    ['Trial mismatch: modelmd has %d columns but OX_get_odor returned %d trials. ' ...
     'The code assumes columns of modelmd align exactly with odor/category order.'], ...
    size(modelmd, 2), nTrials);
assert(numel(odor_raw) == nTrials, 'odor and category must have the same number of trials.');

if strcmp(target, 'context')
    class_order = ["PERSON", "FOOD", "LOCATION", "CONTROL"];
    category_upper = upper(strtrim(category_raw));
    [is_known, y] = ismember(category_upper, class_order);
    assert(all(is_known), 'Unknown context labels found: %s', ...
        strjoin(unique(category_upper(~is_known)), ', '));
    class_values = cellstr(class_order(:));
else
    odor_raw = odor_raw(:);
    if isnumeric(odor_raw) || islogical(odor_raw)
        class_values_numeric = unique(odor_raw, 'sorted');
        [is_known, y] = ismember(odor_raw, class_values_numeric);
        assert(all(is_known), 'Unable to encode one or more odor labels.');
        class_values = class_values_numeric;
    else
        odor_string = string(odor_raw);
        class_values_string = unique(odor_string, 'sorted');
        [is_known, y] = ismember(odor_string, class_values_string);
        assert(all(is_known), 'Unable to encode one or more odor labels.');
        class_values = cellstr(class_values_string);
    end
end
y = double(y(:));
nClasses = numel(unique(y));
chance = 1 / nClasses;
evidence_null = 0;
assert(nClasses >= 2, 'Evidence scoring requires at least two classes.');

run_ids = make_run_ids(nTrials, opts.NumRuns, opts.RunLabels);
unique_runs = unique(run_ids, 'stable');
nRuns = numel(unique_runs);
modelmd = center_model_within_runs(modelmd, run_ids);
if strcmp(cv_method, 'leave-one-run-out')
    cv_ids = run_ids;
else
    cv_ids = (1:nTrials)';
end
is_loo = strcmp(cv_method, 'leave-one-out');

fprintf('Subject: %s (%s)\n', subjname, subjname_real);
fprintf('Trials: %d | Runs: %d | Classes: %d | Evidence null: %.1f | CV: %s\n', ...
    nTrials, nRuns, nClasses, evidence_null, cv_method);
fprintf('Subtracted each voxel''s mean across trials within each run.\n');
print_class_counts(y, class_values);

% Every fold must retain every class in training.
unique_cv_ids = unique(cv_ids, 'stable');
for fi = 1:numel(unique_cv_ids)
    train_mask = cv_ids ~= unique_cv_ids(fi);
    missing = setdiff((1:nClasses)', unique(y(train_mask)));
    if ~isempty(missing)
        error(['CV fold %g cannot be held out because the remaining training data lack ' ...
               'class(es): %s. Use a different cross-validation method or grouping.'], ...
               unique_cv_ids(fi), mat2str(missing'));
    end
end

if opts.NumPermutations > 0
    movable_runs = false(nRuns, 1);
    for ri = 1:nRuns
        movable_runs(ri) = numel(unique(y(run_ids == unique_runs(ri)))) > 1;
    end
    if ~any(movable_runs)
        error(['Within-run permutation cannot change any labels because every run contains ' ...
               'only one class. The exchangeability scheme must be changed after checking ' ...
               'the task structure.']);
    end
    fprintf('Runs containing >1 label and contributing to within-run shuffling: %d/%d\n', ...
        sum(movable_runs), nRuns);
end

%% Load mask and verify row-to-volume mapping
gmmask = spm_vol(mask_file);
gm_mask = spm_read_vols(gmmask) > 0;
gm_inds = find(gm_mask);

assert(size(modelmd, 1) == numel(gm_inds), ...
    ['Mapping check failed: modelmd has %d rows but find(gm_mask) has %d voxels. ' ...
     'This implementation requires modelmd row i to correspond to gm_inds(i).'], ...
    size(modelmd, 1), numel(gm_inds));
assert(numel(R2) == size(modelmd, 1), ...
    'R2 has %d elements but modelmd has %d voxel rows.', numel(R2), size(modelmd, 1));

selected_model_idx = find(isfinite(R2) & R2 > opts.R2Threshold);
assert(~isempty(selected_model_idx), 'No voxels survive R2 > %.4g.', opts.R2Threshold);
if isfinite(opts.MaxCenters) && opts.MaxCenters < numel(selected_model_idx)
    subset_position = unique(round(linspace(1, numel(selected_model_idx), opts.MaxCenters)));
    selected_model_idx = selected_model_idx(subset_position);
    warning('MaxCenters=%d: running a sparse smoke-test map, not a final whole-brain analysis.', ...
        numel(selected_model_idx));
end
center_linear_idx = gm_inds(selected_model_idx);

fprintf('Gray-matter voxels: %d | R2-selected searchlight centers/features: %d\n', ...
    numel(gm_inds), numel(selected_model_idx));

%% Build spherical neighborhoods
fprintf('Building radius-%g voxel searchlights...\n', opts.SearchlightRadius);
[neighborhoods, neighborhood_sizes] = build_neighborhoods( ...
    size(gm_mask), gm_inds, center_linear_idx, selected_model_idx, ...
    opts.SearchlightRadius, opts.RestrictFeaturesToR2);

fprintf('Searchlights: %d | median size: %.0f | range: %d-%d voxels\n', ...
    numel(neighborhoods), median(neighborhood_sizes), ...
    min(neighborhood_sizes), max(neighborhood_sizes));

%% Parallel setup
use_parallel = opts.UseParallel && license('test', 'Distrib_Computing_Toolbox');
if opts.UseParallel && ~use_parallel
    warning('Parallel Computing Toolbox unavailable; running serially.');
end
if use_parallel && isempty(gcp('nocreate'))
    parpool;
end

%% Observed evidence scoring
fprintf('Running observed %s %s evidence scoring...\n', cv_method, target);
tic;
if use_parallel
    model_const = parallel.pool.Constant(modelmd);
    model_cleanup = onCleanup(@() delete(model_const)); %#ok<NASGU>
    voxel_mean_evidence = score_searchlight_map_parallel( ...
        model_const, neighborhoods, y, cv_ids, is_loo, nClasses, ...
        opts.MinVoxels, opts.DemeanPatterns, opts.ScaleVoxels);
else
    voxel_mean_evidence = score_searchlight_map_serial( ...
        modelmd, neighborhoods, y, cv_ids, is_loo, nClasses, ...
        opts.MinVoxels, opts.DemeanPatterns, opts.ScaleVoxels);
end
fprintf('Observed map complete in %.1f minutes.\n', toc / 60);

valid_centers = isfinite(voxel_mean_evidence);
assert(any(valid_centers), 'No searchlights produced a finite mean evidence value.');

%% Full-map permutation inference
p_uncorrected = nan(size(voxel_mean_evidence));
q_fdr = nan(size(voxel_mean_evidence));
p_fwe = nan(size(voxel_mean_evidence));
sig_fdr05 = false(size(voxel_mean_evidence));
sig_fwe05 = false(size(voxel_mean_evidence));
max_null = [];
fwe_meanEvidence_threshold05 = NaN;

if opts.NumPermutations > 0
    fprintf(['Running %d full-map label permutations. Labels are shuffled within run; ' ...
             'the complete CV pipeline is rerun each time.\n'], opts.NumPermutations);
    rng(opts.RandomSeed, 'twister');
    exceed_count = zeros(size(voxel_mean_evidence), 'uint32');
    max_null = nan(opts.NumPermutations, 1);

    perm_tic = tic;
    for perm_idx = 1:opts.NumPermutations
        y_perm = permute_within_blocks(y, run_ids);

        if use_parallel
            perm_evidence = score_searchlight_map_parallel( ...
                model_const, neighborhoods, y_perm, cv_ids, is_loo, nClasses, ...
                opts.MinVoxels, opts.DemeanPatterns, opts.ScaleVoxels);
        else
            perm_evidence = score_searchlight_map_serial( ...
                modelmd, neighborhoods, y_perm, cv_ids, is_loo, nClasses, ...
                opts.MinVoxels, opts.DemeanPatterns, opts.ScaleVoxels);
        end

        comparable = valid_centers & isfinite(perm_evidence);
        exceed_count(comparable) = exceed_count(comparable) + ...
            uint32(perm_evidence(comparable) >= voxel_mean_evidence(comparable));
        max_null(perm_idx) = max(perm_evidence(comparable) - evidence_null);

        if perm_idx == 1 || mod(perm_idx, 10) == 0 || perm_idx == opts.NumPermutations
            elapsed_min = toc(perm_tic) / 60;
            fprintf('  permutation %d/%d | elapsed %.1f min\n', ...
                perm_idx, opts.NumPermutations, elapsed_min);
        end
    end

    p_uncorrected(valid_centers) = ...
        (double(exceed_count(valid_centers)) + 1) ./ (opts.NumPermutations + 1);
    q_fdr(valid_centers) = bh_fdr_qvalues(p_uncorrected(valid_centers));

    max_null_sorted = sort(max_null);
    for vi = find(valid_centers)'
        p_fwe(vi) = (1 + sum(max_null >= ...
            (voxel_mean_evidence(vi) - evidence_null))) / ...
            (opts.NumPermutations + 1);
    end

    sig_fdr05 = q_fdr < 0.05;
    sig_fwe05 = p_fwe < 0.05;

    % A global max-stat threshold is valid because every map uses the same
    % mean-evidence metric and null value. The p_FWE map remains the
    % preferred output for reporting.
    threshold_index = min(opts.NumPermutations, ...
        max(1, ceil(0.95 * (opts.NumPermutations + 1))));
    fwe_meanEvidence_threshold05 = evidence_null + max_null_sorted(threshold_index);

    fprintf('Permutation inference complete. FDR q<.05: %d centers | max-stat FWE p<.05: %d centers\n', ...
        sum(sig_fdr05), sum(sig_fwe05));
    fprintf('Max-stat FWE .05 mean evidence threshold: %.4f\n', ...
        fwe_meanEvidence_threshold05);
end

%% Package results
if strlength(string(opts.OutputDir)) == 0
    output_dir = fullfile(base_outdir, ...
        sprintf('searchlight_score_%s_%s', target, cv_short_name));
else
    output_dir = char(string(opts.OutputDir));
end

results = struct();
results.subject_index = subjidx;
results.subject_name = subjname;
results.subject_name_real = subjname_real;
results.target = target;
results.class_values = class_values;
results.class_counts = accumarray(y, 1, [nClasses, 1]);
results.chance = chance;
results.evidence_null = evidence_null;
results.run_ids = run_ids;
results.center_model_indices = selected_model_idx;
results.center_linear_indices = center_linear_idx;
results.neighborhood_sizes = neighborhood_sizes;
results.meanEvidence = voxel_mean_evidence;
results.meanEvidence_minus_null = voxel_mean_evidence - evidence_null;
results.p_uncorrected = p_uncorrected;
results.q_fdr = q_fdr;
results.p_fwe = p_fwe;
results.sig_fdr05 = sig_fdr05;
results.sig_fwe05 = sig_fwe05;
results.max_null = max_null;
results.fwe_meanEvidence_threshold05 = fwe_meanEvidence_threshold05;
results.options = opts;
results.output_dir = output_dir;
results.mapping_assertion = 'modelmd row i == find(gm_mask)(i)';
results.cv = cv_method;
results.permutation = 'labels shuffled within run; complete CV pipeline rerun';
results.metric = 'meanEvidence';
results.scoring = ['Pearson similarity to correct template minus mean Pearson ' ...
    'similarity to all incorrect templates, averaged across test trials'];

%% Save outputs
if opts.SaveOutputs
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end
    prefix = sprintf('%s_score_subj%d', target, subjidx);

    write_float_map(gmmask, size(gm_mask), center_linear_idx, voxel_mean_evidence, ...
        fullfile(output_dir, [prefix '_meanEvidence.nii']));

    if opts.NumPermutations > 0
        write_float_map(gmmask, size(gm_mask), center_linear_idx, p_uncorrected, ...
            fullfile(output_dir, [prefix '_p_uncorrected.nii']));
        write_float_map(gmmask, size(gm_mask), center_linear_idx, q_fdr, ...
            fullfile(output_dir, [prefix '_q_fdr.nii']));
        write_float_map(gmmask, size(gm_mask), center_linear_idx, p_fwe, ...
            fullfile(output_dir, [prefix '_p_fwe_maxstat.nii']));
        write_binary_map(gmmask, size(gm_mask), center_linear_idx, sig_fdr05, ...
            fullfile(output_dir, [prefix '_sig_fdr_q05.nii']));
        write_binary_map(gmmask, size(gm_mask), center_linear_idx, sig_fwe05, ...
            fullfile(output_dir, [prefix '_sig_fwe_p05.nii']));
    end

    save(fullfile(output_dir, [prefix '_results.mat']), 'results', '-v7.3');
    fprintf('Saved outputs to %s\n', output_dir);
end
end

%% ------------------------------------------------------------------------
function [cv_method, cv_short_name] = normalize_cv_method(requested_method)
requested_method = lower(strtrim(char(string(requested_method))));
switch requested_method
    case {'leave-one-run-out', 'loro'}
        cv_method = 'leave-one-run-out';
        cv_short_name = 'loro';
    case {'leave-one-out', 'loo'}
        cv_method = 'leave-one-out';
        cv_short_name = 'loo';
    otherwise
        error(['CrossValidation must be ''leave-one-run-out'' (or ''loro'') ' ...
               'or ''leave-one-out'' (or ''loo'').']);
end
end

%% ------------------------------------------------------------------------
function mriroot = resolve_mri_root(requested_root)
if strlength(string(requested_root)) > 0
    mriroot = char(string(requested_root));
    assert(isfolder(mriroot), 'MRIRoot does not exist: %s', mriroot);
    return;
end

candidates = { ...
    '/Users/qhyang/Desktop/OX_DATA/MRI', ...
    '/Volumes/ExtremeSSD/OX_DATA/MRI'};

mriroot = '';
for i = 1:numel(candidates)
    if isfolder(candidates{i})
        mriroot = candidates{i};
        break;
    end
end
assert(~isempty(mriroot), ...
    'Could not auto-detect MRIRoot. Pass it explicitly with ''MRIRoot'', path.');
end

%% ------------------------------------------------------------------------
function run_ids = make_run_ids(nTrials, nRuns, supplied_run_ids)
if ~isempty(supplied_run_ids)
    run_ids = supplied_run_ids(:);
    assert(numel(run_ids) == nTrials, ...
        'RunLabels has %d elements; expected %d.', numel(run_ids), nTrials);
    assert(all(isfinite(run_ids)), 'RunLabels must be finite.');
    return;
end

assert(mod(nTrials, nRuns) == 0, ...
    ['Cannot infer equal contiguous runs: %d trials is not divisible by %d runs. ' ...
     'Supply exact RunLabels.'], nTrials, nRuns);
trials_per_run = nTrials / nRuns;
run_ids = repelem((1:nRuns)', trials_per_run);
warning(['RunLabels were inferred as %d contiguous equal-sized runs with %d trials/run. ' ...
         'Verify that GLMsingle trial order is run 1 trials, then run 2 trials, etc. ' ...
         'For the final analysis, passing exact RunLabels is safer.'], ...
         nRuns, trials_per_run);
end

%% ------------------------------------------------------------------------
function modelmd = center_model_within_runs(modelmd, run_ids)
runs = unique(run_ids, 'stable');
for ri = 1:numel(runs)
    run_mask = run_ids == runs(ri);
    modelmd(:, run_mask) = modelmd(:, run_mask) - mean(modelmd(:, run_mask), 2);
end
end

%% ------------------------------------------------------------------------
function print_class_counts(y, class_values)
counts = accumarray(y, 1);
fprintf('Class counts:\n');
for c = 1:numel(counts)
    if isnumeric(class_values)
        label_text = num2str(class_values(c));
    else
        label_text = char(string(class_values{c}));
    end
    fprintf('  %2d  %-20s  %d\n', c, label_text, counts(c));
end
end

%% ------------------------------------------------------------------------
function [neighborhoods, neighborhood_sizes] = build_neighborhoods(vol_size, gm_inds, center_linear_idx, selected_model_idx, radius, restrict_features_to_r2)
idx_vol = nan(vol_size);
if restrict_features_to_r2
    % Both centers and neighborhood features obey the same GM + R2 rule.
    idx_vol(center_linear_idx) = selected_model_idx;
else
    % Reproduce the old script: centers pass R2, but features may be any GM voxel.
    idx_vol(gm_inds) = 1:numel(gm_inds);
end

[xc, yc, zc] = ind2sub(vol_size, center_linear_idx);
nCenters = numel(center_linear_idx);
neighborhoods = cell(nCenters, 1);
neighborhood_sizes = zeros(nCenters, 1);

for vi = 1:nCenters
    x_range = floor(xc(vi) - radius):ceil(xc(vi) + radius);
    y_range = floor(yc(vi) - radius):ceil(yc(vi) + radius);
    z_range = floor(zc(vi) - radius):ceil(zc(vi) + radius);

    [xg, yg, zg] = ndgrid(x_range, y_range, z_range);
    in_bounds = xg >= 1 & xg <= vol_size(1) & ...
                yg >= 1 & yg <= vol_size(2) & ...
                zg >= 1 & zg <= vol_size(3);
    xg = xg(in_bounds);
    yg = yg(in_bounds);
    zg = zg(in_bounds);

    in_sphere = ((xg - xc(vi)).^2 + (yg - yc(vi)).^2 + (zg - zc(vi)).^2) <= radius^2;
    xg = xg(in_sphere);
    yg = yg(in_sphere);
    zg = zg(in_sphere);

    lin_idx = sub2ind(vol_size, xg, yg, zg);
    model_idx = idx_vol(lin_idx);
    model_idx = model_idx(isfinite(model_idx));

    neighborhoods{vi} = model_idx(:)';
    neighborhood_sizes(vi) = numel(model_idx);
end
end

%% ------------------------------------------------------------------------
function meanEvidence = score_searchlight_map_parallel(model_const, neighborhoods, y, cv_ids, is_loo, nClasses, min_voxels, demean_patterns, scale_voxels)
nCenters = numel(neighborhoods);
meanEvidence = nan(nCenters, 1);
parfor vi = 1:nCenters
    meanEvidence(vi) = score_one_searchlight( ...
        model_const.Value, neighborhoods{vi}, y, cv_ids, is_loo, nClasses, ...
        min_voxels, demean_patterns, scale_voxels);
end
end

%% ------------------------------------------------------------------------
function meanEvidence = score_searchlight_map_serial(modelmd, neighborhoods, y, cv_ids, is_loo, nClasses, min_voxels, demean_patterns, scale_voxels)
nCenters = numel(neighborhoods);
meanEvidence = nan(nCenters, 1);
for vi = 1:nCenters
    meanEvidence(vi) = score_one_searchlight( ...
        modelmd, neighborhoods{vi}, y, cv_ids, is_loo, nClasses, ...
        min_voxels, demean_patterns, scale_voxels);
end
end

%% ------------------------------------------------------------------------
function meanEvidence = score_one_searchlight(modelmd, voxel_inds, y, cv_ids, is_loo, nClasses, min_voxels, demean_patterns, scale_voxels)
if numel(voxel_inds) < min_voxels
    meanEvidence = NaN;
    return;
end

X = double(modelmd(voxel_inds, :)'); % trials x voxels
finite_features = all(isfinite(X), 1);
X = X(:, finite_features);
if size(X, 2) < min_voxels
    meanEvidence = NaN;
    return;
end

if demean_patterns
    X = X - mean(X, 2);
end

if is_loo
    trial_evidence = score_leave_one_out(X, y, nClasses, min_voxels, scale_voxels);
else
    trial_evidence = score_grouped_folds(X, y, cv_ids, min_voxels, scale_voxels);
end

if ~all(isfinite(trial_evidence))
    meanEvidence = NaN;
else
    meanEvidence = mean(trial_evidence);
end
end

%% ------------------------------------------------------------------------
function trial_evidence = score_leave_one_out(X, y, nClasses, min_voxels, scale_voxels)
nTrials = size(X, 1);
trial_evidence = nan(nTrials, 1);
class_counts = accumarray(y, 1, [nClasses, 1]);
class_sums = zeros(nClasses, size(X, 2));
for ci = 1:nClasses
    class_sums(ci, :) = sum(X(y == ci, :), 1);
end
class_means = class_sums ./ class_counts;

if scale_voxels
    global_mu = mean(X, 1);
    X_centered = X - global_mu;
    total_centered_ss = sum(X_centered.^2, 1);
end

for ti = 1:nTrials
    templates = class_means;
    test_class = y(ti);
    templates(test_class, :) = ...
        (class_sums(test_class, :) - X(ti, :)) ./ (class_counts(test_class) - 1);

    if scale_voxels
        train_mu = global_mu - X_centered(ti, :) ./ (nTrials - 1);
        train_ss = total_centered_ss - ...
            (nTrials / (nTrials - 1)) .* X_centered(ti, :).^2;
        train_sd = sqrt(max(train_ss, 0) ./ (nTrials - 2));
        usable = isfinite(train_mu) & isfinite(train_sd) & train_sd > eps;
        if sum(usable) < min_voxels
            continue;
        end
        templates = (templates(:, usable) - train_mu(usable)) ./ train_sd(usable);
        Xtest = (X(ti, usable) - train_mu(usable)) ./ train_sd(usable);
    else
        Xtest = X(ti, :);
    end

    trial_evidence(ti) = evidence_from_templates( ...
        Xtest, templates, (1:nClasses)', test_class);
end
end

%% ------------------------------------------------------------------------
function trial_evidence = score_grouped_folds(X, y, cv_ids, min_voxels, scale_voxels)
folds = unique(cv_ids, 'stable');
trial_evidence = nan(size(y));

for fi = 1:numel(folds)
    test_mask = cv_ids == folds(fi);
    train_mask = ~test_mask;

    Xtrain = X(train_mask, :);
    Xtest = X(test_mask, :);
    ytrain = y(train_mask);

    if scale_voxels
        train_mu = mean(Xtrain, 1);
        train_sd = std(Xtrain, 0, 1);
        usable = isfinite(train_mu) & isfinite(train_sd) & train_sd > eps;
        if sum(usable) < min_voxels
            continue;
        end
        Xtrain = (Xtrain(:, usable) - train_mu(usable)) ./ train_sd(usable);
        Xtest = (Xtest(:, usable) - train_mu(usable)) ./ train_sd(usable);
    end

    class_labels = unique(ytrain, 'sorted');
    templates = nan(numel(class_labels), size(Xtrain, 2));
    for ci = 1:numel(class_labels)
        templates(ci, :) = mean(Xtrain(ytrain == class_labels(ci), :), 1);
    end

    trial_evidence(test_mask) = evidence_from_templates( ...
        Xtest, templates, class_labels, y(test_mask));
end
end

%% ------------------------------------------------------------------------
function evidence = evidence_from_templates(Xtest, templates, class_labels, ytest)
% Pearson correlation across voxels, implemented as centered cosine similarity.
templates = templates - mean(templates, 2);
Xtest_centered = Xtest - mean(Xtest, 2);

template_norm = sqrt(sum(templates.^2, 2));
test_norm = sqrt(sum(Xtest_centered.^2, 2));
valid_templates = isfinite(template_norm) & template_norm > eps;
valid_tests = isfinite(test_norm) & test_norm > eps;
evidence = nan(size(Xtest, 1), 1);

if ~all(valid_templates)
    return;
end

templates = templates ./ template_norm;
Xtest_centered(valid_tests, :) = Xtest_centered(valid_tests, :) ./ test_norm(valid_tests);
scores = Xtest_centered(valid_tests, :) * templates';
[is_known, correct_idx] = ismember(ytest(valid_tests), class_labels);
finite_scores = all(isfinite(scores), 2);
usable_tests = is_known & finite_scores;
if ~any(usable_tests)
    return;
end

scores = scores(usable_tests, :);
correct_idx = correct_idx(usable_tests);
correct_linear_idx = sub2ind(size(scores), (1:size(scores, 1))', correct_idx);
correct_similarity = scores(correct_linear_idx);
incorrect_mean = (sum(scores, 2) - correct_similarity) ./ (numel(class_labels) - 1);
valid_test_idx = find(valid_tests);
valid_test_idx = valid_test_idx(usable_tests);
evidence(valid_test_idx) = correct_similarity - incorrect_mean;
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
function write_float_map(header, vol_size, center_linear_idx, values, filename)
vol = nan(vol_size, 'single');
vol(center_linear_idx) = single(values);
out_header = header;
out_header.fname = filename;
out_header.dt = [16, 0]; % float32
out_header.pinfo = [1; 0; 0];
spm_write_vol(out_header, vol);
end

%% ------------------------------------------------------------------------
function write_binary_map(header, vol_size, center_linear_idx, values, filename)
vol = zeros(vol_size, 'single');
vol(center_linear_idx) = single(values);
out_header = header;
out_header.fname = filename;
out_header.dt = [16, 0]; % float32 for broad SPM compatibility
out_header.pinfo = [1; 0; 0];
spm_write_vol(out_header, vol);
end
