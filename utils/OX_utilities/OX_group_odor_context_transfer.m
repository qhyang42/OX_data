function group_results = OX_group_odor_context_transfer(subject_ids, varargin)
%OX_GROUP_ODOR_CONTEXT_TRANSFER Aggregate ROI context-transfer results.
%
%   group_results = OX_group_odor_context_transfer(subject_ids, Name, Value, ...)
%
% Exact one-sided sign-flip inference is performed across subjects. The
% primary max-statistic family contains the bilateral ROIs. CONTROL target
% follow-ups form a separate family across all ROI-by-target combinations.

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'subject_ids', @(x) isnumeric(x) && isvector(x) && ...
    all(isfinite(x)) && all(x == round(x)));
addParameter(p, 'MRIRoot', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'ROISelection', 'old', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'OutputDir', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'MakePlots', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'SaveOutputs', true, @(x) islogical(x) && isscalar(x));
parse(p, subject_ids, varargin{:});
opts = p.Results;
roi_selection = OX_normalize_decoding_roi_selection(opts.ROISelection);
opts.ROISelection = roi_selection;

run_inference_self_tests();

subject_ids = subject_ids(:)';
n_subjects = numel(subject_ids);
assert(n_subjects >= 2, 'At least two subjects are required.');
mri_root = resolve_mri_root(opts.MRIRoot);
if strlength(string(opts.OutputDir)) == 0
    if strcmp(roi_selection, 'old')
        output_name = 'roi_odor_context_template_loro';
    else
        output_name = sprintf('roi_%s_odor_context_template_loro', roi_selection);
    end
    output_dir = fullfile(mri_root, 'group', output_name);
else
    output_dir = char(string(opts.OutputDir));
end

subject_results = cell(n_subjects, 1);
for subject_idx = 1:n_subjects
    subject_id = subject_ids(subject_idx);
    if strcmp(roi_selection, 'old')
        subject_output_name = 'roi_odor_context_template_loro';
    else
        subject_output_name = sprintf('roi_%s_odor_context_template_loro', roi_selection);
    end
    filename = fullfile(mri_root, sprintf('subj_%d', subject_id), 'nifti', ...
        'single_trial_by_category', subject_output_name, ...
        sprintf('odor_context_template_subj%d_loro_results.mat', subject_id));
    assert(isfile(filename), 'Missing subject result: %s', filename);
    loaded = load(filename, 'results');
    subject_results{subject_idx} = loaded.results;
end

reference = subject_results{1};
assert_result_selection(reference, roi_selection, subject_ids(1));
roi_metadata = reference.roi_metadata;
roi_names = string(roi_metadata.roi_name);
n_rois = height(roi_metadata);
semantic_contexts = string(reference.analysis.semantic_context_order);
assert(isequal(semantic_contexts(:), ["PERSON"; "FOOD"; "LOCATION"]), ...
    'Unexpected semantic context ordering.');

semantic_evidence = nan(n_subjects, n_rois, 3, 3);
semantic_accuracy = nan(n_subjects, n_rois, 3, 3);
within_evidence = nan(n_subjects, n_rois);
cross_evidence = nan(n_subjects, n_rois);
primary_evidence_raw = nan(n_subjects, n_rois);
within_accuracy = nan(n_subjects, n_rois);
cross_accuracy = nan(n_subjects, n_rois);
primary_accuracy_raw = nan(n_subjects, n_rois);
control_evidence = nan(n_subjects, n_rois, 3);
control_accuracy = nan(n_subjects, n_rois, 3);
control_mean_evidence_raw = nan(n_subjects, n_rois);
control_mean_accuracy_raw = nan(n_subjects, n_rois);
control_within_evidence = nan(n_subjects, n_rois);
control_within_accuracy = nan(n_subjects, n_rois);
primary_evidence_null_mean = nan(n_subjects, n_rois);
primary_accuracy_null_mean = nan(n_subjects, n_rois);
control_mean_evidence_null_mean = nan(n_subjects, n_rois);
control_mean_accuracy_null_mean = nan(n_subjects, n_rois);
control_evidence_null_mean = nan(n_subjects, n_rois, 3);
control_accuracy_null_mean = nan(n_subjects, n_rois, 3);

for subject_idx = 1:n_subjects
    item = subject_results{subject_idx};
    assert_result_selection(item, roi_selection, subject_ids(subject_idx));
    assert(isequal(string(item.roi_metadata.roi_name), roi_names), ...
        'ROI ordering differs for subject %d.', subject_ids(subject_idx));
    assert(isequal(string(item.analysis.semantic_context_order(:)), semantic_contexts(:)), ...
        'Context ordering differs for subject %d.', subject_ids(subject_idx));
    assert(~isempty(item.null.primary_evidence_effect), ...
        'Subject %d has no permutation null distribution.', subject_ids(subject_idx));

    semantic_evidence(subject_idx, :, :, :) = item.semantic.evidence;
    semantic_accuracy(subject_idx, :, :, :) = item.semantic.accuracy;
    within_evidence(subject_idx, :) = item.primary.within_context_evidence;
    cross_evidence(subject_idx, :) = item.primary.cross_context_evidence;
    primary_evidence_raw(subject_idx, :) = item.primary.evidence_effect;
    within_accuracy(subject_idx, :) = item.primary.within_context_accuracy;
    cross_accuracy(subject_idx, :) = item.primary.cross_context_accuracy;
    primary_accuracy_raw(subject_idx, :) = item.primary.accuracy_effect;
    control_evidence(subject_idx, :, :) = item.control_transfer.evidence;
    control_accuracy(subject_idx, :, :) = item.control_transfer.accuracy;
    control_mean_evidence_raw(subject_idx, :) = item.control_transfer.mean_evidence;
    control_mean_accuracy_raw(subject_idx, :) = item.control_transfer.mean_accuracy;
    control_within_evidence(subject_idx, :) = item.control_within.evidence;
    control_within_accuracy(subject_idx, :) = item.control_within.accuracy;
    primary_evidence_null_mean(subject_idx, :) = item.null.mean_primary_evidence_effect;
    primary_accuracy_null_mean(subject_idx, :) = item.null.mean_primary_accuracy_effect;
    control_mean_evidence_null_mean(subject_idx, :) = item.null.mean_control_transfer_evidence;
    control_mean_accuracy_null_mean(subject_idx, :) = item.null.mean_control_transfer_accuracy;
    control_evidence_null_mean(subject_idx, :, :) = item.null.mean_control_target_evidence;
    control_accuracy_null_mean(subject_idx, :, :) = item.null.mean_control_target_accuracy;
end

all_group_values = [primary_evidence_raw, primary_accuracy_raw, ...
    control_mean_evidence_raw, control_mean_accuracy_raw];
assert(all(isfinite(all_group_values), 'all'), ...
    'All selected bilateral ROIs must have finite results in every subject.');

primary_evidence_centered = primary_evidence_raw - primary_evidence_null_mean;
primary_accuracy_centered = primary_accuracy_raw - primary_accuracy_null_mean;
control_mean_evidence_centered = control_mean_evidence_raw - control_mean_evidence_null_mean;
control_mean_accuracy_centered = control_mean_accuracy_raw - control_mean_accuracy_null_mean;
control_evidence_centered = control_evidence - control_evidence_null_mean;
control_accuracy_centered = control_accuracy - control_accuracy_null_mean;

primary_evidence_inference = exact_sign_flip(primary_evidence_centered);
primary_accuracy_inference = exact_sign_flip(primary_accuracy_centered);
control_mean_evidence_inference = exact_sign_flip(control_mean_evidence_centered);
control_mean_accuracy_inference = exact_sign_flip(control_mean_accuracy_centered);
control_target_evidence_flat = reshape(control_evidence_centered, n_subjects, []);
control_target_accuracy_flat = reshape(control_accuracy_centered, n_subjects, []);
control_target_evidence_inference = exact_sign_flip(control_target_evidence_flat);
control_target_accuracy_inference = exact_sign_flip(control_target_accuracy_flat);

primary_table = make_primary_table(roi_metadata, n_subjects, ...
    within_evidence, cross_evidence, primary_evidence_raw, primary_evidence_null_mean, ...
    primary_evidence_inference, within_accuracy, cross_accuracy, ...
    primary_accuracy_raw, primary_accuracy_null_mean, primary_accuracy_inference);
control_mean_table = make_control_mean_table(roi_metadata, n_subjects, ...
    control_mean_evidence_raw, control_mean_evidence_null_mean, ...
    control_mean_evidence_inference, control_mean_accuracy_raw, ...
    control_mean_accuracy_null_mean, control_mean_accuracy_inference, ...
    control_within_evidence, control_within_accuracy);
control_target_table = make_control_target_table(roi_metadata, semantic_contexts, ...
    n_subjects, control_evidence, control_evidence_null_mean, ...
    control_target_evidence_inference, control_accuracy, ...
    control_accuracy_null_mean, control_target_accuracy_inference);
subject_table = make_subject_table(subject_ids, roi_names, within_evidence, ...
    cross_evidence, primary_evidence_raw, within_accuracy, cross_accuracy, ...
    primary_accuracy_raw, control_evidence, control_accuracy, ...
    control_within_evidence, control_within_accuracy, semantic_contexts);

group_results = struct();
group_results.analysis = struct( ...
    'name', 'semantic_context_modulation_of_odor_identity', ...
    'roi_selection', roi_selection, ...
    'subjects', subject_ids, ...
    'n_subjects', n_subjects, ...
    'semantic_context_order', {cellstr(semantic_contexts)}, ...
    'primary_alternative', 'semantic diagonal evidence > semantic off-diagonal evidence', ...
    'secondary_alternative', 'CONTROL-trained odor identity transfers above empirical null', ...
    'inference', 'exact one-sided subject sign flips', ...
    'n_exact_sign_flips', 2^n_subjects, ...
    'minimum_attainable_p', 1 / (2^n_subjects), ...
    'multiple_comparison_control', 'max-statistic FWE; BH-FDR also reported');
group_results.roi_metadata = roi_metadata;
group_results.subject_level = struct( ...
    'semantic_evidence', semantic_evidence, ...
    'semantic_accuracy', semantic_accuracy, ...
    'within_evidence', within_evidence, ...
    'cross_evidence', cross_evidence, ...
    'primary_evidence_raw', primary_evidence_raw, ...
    'primary_evidence_null_centered', primary_evidence_centered, ...
    'within_accuracy', within_accuracy, ...
    'cross_accuracy', cross_accuracy, ...
    'primary_accuracy_raw', primary_accuracy_raw, ...
    'primary_accuracy_null_centered', primary_accuracy_centered, ...
    'control_transfer_evidence', control_evidence, ...
    'control_transfer_accuracy', control_accuracy, ...
    'control_within_evidence', control_within_evidence, ...
    'control_within_accuracy', control_within_accuracy);
group_results.inference = struct( ...
    'primary_evidence', primary_evidence_inference, ...
    'primary_accuracy', primary_accuracy_inference, ...
    'control_mean_evidence', control_mean_evidence_inference, ...
    'control_mean_accuracy', control_mean_accuracy_inference, ...
    'control_target_evidence', control_target_evidence_inference, ...
    'control_target_accuracy', control_target_accuracy_inference);
group_results.primary_table = primary_table;
group_results.control_mean_table = control_mean_table;
group_results.control_target_table = control_target_table;
group_results.subject_table = subject_table;
group_results.output_dir = output_dir;
group_results.matlab_version = version;

if opts.SaveOutputs || opts.MakePlots
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end
end
if opts.SaveOutputs
    save(fullfile(output_dir, 'group_odor_context_template_results.mat'), ...
        'group_results', '-v7.3');
    writetable(primary_table, fullfile(output_dir, 'group_primary_semantic_effect.tsv'), ...
        'FileType', 'text', 'Delimiter', '\t');
    writetable(control_mean_table, fullfile(output_dir, 'group_control_transfer_average.tsv'), ...
        'FileType', 'text', 'Delimiter', '\t');
    writetable(control_target_table, fullfile(output_dir, 'group_control_transfer_by_target.tsv'), ...
        'FileType', 'text', 'Delimiter', '\t');
    writetable(subject_table, fullfile(output_dir, 'subject_level_long_summary.tsv'), ...
        'FileType', 'text', 'Delimiter', '\t');
end
if opts.MakePlots
    create_group_plots(output_dir, roi_names, semantic_contexts, ...
        semantic_evidence, semantic_accuracy, within_evidence, cross_evidence, ...
        within_accuracy, cross_accuracy, control_evidence, control_accuracy, ...
        reference.chance_accuracy);
end
fprintf(['Group inference complete: %d subjects, %d exact sign flips; ' ...
    'minimum attainable one-sided p = %.5f.\n'], ...
    n_subjects, 2^n_subjects, 1 / (2^n_subjects));
end

function assert_result_selection(result, expected, subject_id)
if isfield(result, 'analysis') && isfield(result.analysis, 'roi_selection')
    actual = OX_normalize_decoding_roi_selection(result.analysis.roi_selection);
else
    actual = 'old';
end
assert(strcmp(actual, expected), ...
    'Subject %d result uses ROISelection=%s; expected %s.', ...
    subject_id, actual, expected);
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

function inference = exact_sign_flip(effects)
[n_subjects, n_tests] = size(effects);
assert(all(isfinite(effects), 'all'), 'Sign-flip effects must be finite.');
n_signs = 2^n_subjects;
sign_matrix = ones(n_signs, n_subjects);
for subject_idx = 1:n_subjects
    sign_matrix(:, subject_idx) = 2 * bitget((0:n_signs-1)', subject_idx) - 1;
end
null_statistics = (sign_matrix * effects) / n_subjects;
observed = mean(effects, 1);
p_uncorrected = sum(null_statistics >= observed, 1) / n_signs;
max_null = max(null_statistics, [], 2);
p_fwe = sum(max_null >= observed, 1) / n_signs;
q_fdr = bh_fdr(p_uncorrected(:))';
inference = struct( ...
    'observed_null_centered_mean', observed, ...
    'p_uncorrected', p_uncorrected, ...
    'q_fdr', q_fdr, ...
    'p_fwe_maxstat', p_fwe, ...
    'significant_fwe05', p_fwe <= 0.05, ...
    'null_statistics', null_statistics, ...
    'max_null', max_null, ...
    'n_exact_sign_flips', n_signs, ...
    'alternative', 'greater');
assert(numel(observed) == n_tests, 'Sign-flip output size mismatch.');
end

function run_inference_self_tests()
% Five uniformly positive subject effects attain the exact one-sided floor.
test_inference = exact_sign_flip(ones(5, 2));
assert(all(test_inference.p_uncorrected == 1/32), ...
    'Exact sign-flip p-value resolution self-test failed.');
assert(all(test_inference.p_fwe_maxstat >= test_inference.p_uncorrected), ...
    'Max-statistic FWE self-test failed.');
test_q = bh_fdr([0.01; 0.04; 0.03]);
assert(max(abs(test_q - [0.03; 0.04; 0.04])) < 1e-12, ...
    'BH-FDR self-test failed.');
end

function q = bh_fdr(p_values)
p_values = p_values(:);
[sorted_p, order] = sort(p_values, 'ascend');
n = numel(sorted_p);
sorted_q = sorted_p .* n ./ (1:n)';
sorted_q = flipud(cummin(flipud(sorted_q)));
sorted_q = min(sorted_q, 1);
q = nan(n, 1);
q(order) = sorted_q;
end

function output = make_primary_table(metadata, n_subjects, within_evidence, ...
        cross_evidence, effect_evidence, null_evidence, inference_evidence, ...
        within_accuracy, cross_accuracy, effect_accuracy, null_accuracy, inference_accuracy)
output = metadata;
output.n_subjects = repmat(n_subjects, height(metadata), 1);
output.mean_within_evidence = mean(within_evidence, 1)';
output.mean_cross_evidence = mean(cross_evidence, 1)';
output.mean_evidence_effect_raw = mean(effect_evidence, 1)';
output.mean_evidence_null = mean(null_evidence, 1)';
output.mean_evidence_effect_null_centered = inference_evidence.observed_null_centered_mean';
output.p_evidence = inference_evidence.p_uncorrected';
output.q_evidence_fdr = inference_evidence.q_fdr';
output.p_evidence_fwe = inference_evidence.p_fwe_maxstat';
output.sig_evidence_fwe05 = inference_evidence.significant_fwe05';
output.mean_within_accuracy = mean(within_accuracy, 1)';
output.mean_cross_accuracy = mean(cross_accuracy, 1)';
output.mean_accuracy_effect_raw = mean(effect_accuracy, 1)';
output.mean_accuracy_null = mean(null_accuracy, 1)';
output.mean_accuracy_effect_null_centered = inference_accuracy.observed_null_centered_mean';
output.p_accuracy = inference_accuracy.p_uncorrected';
output.q_accuracy_fdr = inference_accuracy.q_fdr';
output.p_accuracy_fwe = inference_accuracy.p_fwe_maxstat';
output.sig_accuracy_fwe05 = inference_accuracy.significant_fwe05';
end

function output = make_control_mean_table(metadata, n_subjects, evidence, null_evidence, ...
        inference_evidence, accuracy, null_accuracy, inference_accuracy, ...
        within_evidence, within_accuracy)
output = metadata;
output.n_subjects = repmat(n_subjects, height(metadata), 1);
output.mean_control_transfer_evidence = mean(evidence, 1)';
output.mean_control_transfer_evidence_null = mean(null_evidence, 1)';
output.mean_control_transfer_evidence_null_centered = ...
    inference_evidence.observed_null_centered_mean';
output.p_evidence = inference_evidence.p_uncorrected';
output.q_evidence_fdr = inference_evidence.q_fdr';
output.p_evidence_fwe = inference_evidence.p_fwe_maxstat';
output.sig_evidence_fwe05 = inference_evidence.significant_fwe05';
output.mean_control_transfer_accuracy = mean(accuracy, 1)';
output.mean_control_transfer_accuracy_null = mean(null_accuracy, 1)';
output.mean_control_transfer_accuracy_null_centered = ...
    inference_accuracy.observed_null_centered_mean';
output.p_accuracy = inference_accuracy.p_uncorrected';
output.q_accuracy_fdr = inference_accuracy.q_fdr';
output.p_accuracy_fwe = inference_accuracy.p_fwe_maxstat';
output.sig_accuracy_fwe05 = inference_accuracy.significant_fwe05';
output.mean_control_within_evidence = mean(within_evidence, 1)';
output.mean_control_within_accuracy = mean(within_accuracy, 1)';
end

function output = make_control_target_table(metadata, contexts, n_subjects, ...
        evidence, null_evidence, inference_evidence, accuracy, null_accuracy, inference_accuracy)
n_rois = height(metadata);
output = metadata(repmat((1:n_rois)', 3, 1), :);
output.target_context = repelem(contexts(:), n_rois);
output.n_subjects = repmat(n_subjects, height(output), 1);
output.mean_evidence = reshape(squeeze(mean(evidence, 1)), [], 1);
output.mean_evidence_null = reshape(squeeze(mean(null_evidence, 1)), [], 1);
output.mean_evidence_null_centered = inference_evidence.observed_null_centered_mean';
output.p_evidence = inference_evidence.p_uncorrected';
output.q_evidence_fdr = inference_evidence.q_fdr';
output.p_evidence_fwe = inference_evidence.p_fwe_maxstat';
output.sig_evidence_fwe05 = inference_evidence.significant_fwe05';
output.mean_accuracy = reshape(squeeze(mean(accuracy, 1)), [], 1);
output.mean_accuracy_null = reshape(squeeze(mean(null_accuracy, 1)), [], 1);
output.mean_accuracy_null_centered = inference_accuracy.observed_null_centered_mean';
output.p_accuracy = inference_accuracy.p_uncorrected';
output.q_accuracy_fdr = inference_accuracy.q_fdr';
output.p_accuracy_fwe = inference_accuracy.p_fwe_maxstat';
output.sig_accuracy_fwe05 = inference_accuracy.significant_fwe05';
end

function output = make_subject_table(subject_ids, roi_names, within_evidence, ...
        cross_evidence, evidence_effect, within_accuracy, cross_accuracy, ...
        accuracy_effect, control_evidence, control_accuracy, ...
        control_within_evidence, control_within_accuracy, contexts)
n_subjects = numel(subject_ids);
n_rois = numel(roi_names);
subject_id = repelem(subject_ids(:), n_rois);
roi_name = repmat(roi_names(:), n_subjects, 1);
output = table(subject_id, roi_name, ...
    reshape(within_evidence', [], 1), reshape(cross_evidence', [], 1), ...
    reshape(evidence_effect', [], 1), reshape(within_accuracy', [], 1), ...
    reshape(cross_accuracy', [], 1), reshape(accuracy_effect', [], 1), ...
    reshape(control_within_evidence', [], 1), reshape(control_within_accuracy', [], 1), ...
    'VariableNames', {'subject_id', 'roi_name', 'within_semantic_evidence', ...
    'cross_semantic_evidence', 'semantic_evidence_effect', ...
    'within_semantic_accuracy', 'cross_semantic_accuracy', ...
    'semantic_accuracy_effect', 'control_within_evidence', 'control_within_accuracy'});
for context_idx = 1:3
    label = lower(char(contexts(context_idx)));
    output.(['control_to_' label '_evidence']) = ...
        reshape(squeeze(control_evidence(:, :, context_idx))', [], 1);
    output.(['control_to_' label '_accuracy']) = ...
        reshape(squeeze(control_accuracy(:, :, context_idx))', [], 1);
end
end

function create_group_plots(output_dir, roi_names, contexts, semantic_evidence, ...
        semantic_accuracy, within_evidence, cross_evidence, within_accuracy, ...
        cross_accuracy, control_evidence, control_accuracy, chance)
mean_semantic_evidence = squeeze(mean(semantic_evidence, 1));
mean_semantic_accuracy = squeeze(mean(semantic_accuracy, 1));
plot_heatmaps(mean_semantic_evidence, roi_names, contexts, ...
    'Semantic odor-template evidence', fullfile(output_dir, 'semantic_evidence_3x3'));
plot_heatmaps(mean_semantic_accuracy, roi_names, contexts, ...
    'Semantic odor-decoding accuracy', fullfile(output_dir, 'semantic_accuracy_3x3'));
plot_paired(within_evidence, cross_evidence, roi_names, ...
    'Template evidence', 'Same versus cross semantic context: evidence', ...
    fullfile(output_dir, 'semantic_same_vs_cross_evidence'), NaN);
plot_paired(within_accuracy, cross_accuracy, roi_names, ...
    'Accuracy', 'Same versus cross semantic context: accuracy', ...
    fullfile(output_dir, 'semantic_same_vs_cross_accuracy'), chance);
plot_control(control_evidence, roi_names, contexts, 'Template evidence', ...
    'CONTROL-trained transfer evidence', ...
    fullfile(output_dir, 'control_transfer_evidence'), NaN);
plot_control(control_accuracy, roi_names, contexts, 'Accuracy', ...
    'CONTROL-trained transfer accuracy', ...
    fullfile(output_dir, 'control_transfer_accuracy'), chance);
end

function plot_heatmaps(values, roi_names, contexts, figure_title, filename_stem)
n_rois = numel(roi_names);
figure_handle = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1500 1100]);
layout = tiledlayout(5, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
value_limits = [min(values(:)), max(values(:))];
if value_limits(1) == value_limits(2)
    value_limits = value_limits + [-1 1] * eps(max(1, abs(value_limits(1))));
end
for roi_idx = 1:n_rois
    nexttile;
    imagesc(squeeze(values(roi_idx, :, :)), value_limits);
    axis image;
    xticks(1:3); yticks(1:3);
    xticklabels(contexts); yticklabels(contexts);
    xtickangle(35);
    title(roi_names(roi_idx), 'Interpreter', 'none', 'FontSize', 8);
end
xlabel(layout, 'Test context');
ylabel(layout, 'Template context');
title(layout, figure_title);
colorbar;
export_plot(figure_handle, filename_stem);
end

function plot_paired(within_values, cross_values, roi_names, y_label, figure_title, filename_stem, reference)
n_rois = numel(roi_names);
figure_handle = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1500 1100]);
layout = tiledlayout(5, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
for roi_idx = 1:n_rois
    nexttile; hold on;
    plot([1 2], [within_values(:, roi_idx), cross_values(:, roi_idx)]', ...
        '-', 'Color', [0.75 0.75 0.75]);
    scatter(ones(size(within_values, 1), 1), within_values(:, roi_idx), 18, 'filled');
    scatter(2 * ones(size(cross_values, 1), 1), cross_values(:, roi_idx), 18, 'filled');
    plot([1 2], [mean(within_values(:, roi_idx)), mean(cross_values(:, roi_idx))], ...
        '-ok', 'LineWidth', 1.3, 'MarkerFaceColor', 'k');
    if isfinite(reference)
        yline(reference, '--k');
    end
    xlim([0.7 2.3]); xticks([1 2]); xticklabels({'Same', 'Cross'});
    title(roi_names(roi_idx), 'Interpreter', 'none', 'FontSize', 8);
end
ylabel(layout, y_label);
title(layout, figure_title);
export_plot(figure_handle, filename_stem);
end

function plot_control(values, roi_names, contexts, y_label, figure_title, filename_stem, reference)
n_rois = numel(roi_names);
figure_handle = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1500 1100]);
layout = tiledlayout(5, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
for roi_idx = 1:n_rois
    nexttile; hold on;
    roi_values = squeeze(values(:, roi_idx, :));
    plot(1:3, roi_values', '-', 'Color', [0.75 0.75 0.75]);
    scatter(repmat(1:3, size(roi_values, 1), 1), roi_values, 16, 'filled');
    plot(1:3, mean(roi_values, 1), '-ok', 'LineWidth', 1.3, 'MarkerFaceColor', 'k');
    if isfinite(reference)
        yline(reference, '--k');
    end
    xlim([0.7 3.3]); xticks(1:3); xticklabels(contexts); xtickangle(35);
    title(roi_names(roi_idx), 'Interpreter', 'none', 'FontSize', 8);
end
xlabel(layout, 'Test context');
ylabel(layout, y_label);
title(layout, figure_title);
export_plot(figure_handle, filename_stem);
end

function export_plot(figure_handle, filename_stem)
exportgraphics(figure_handle, [filename_stem '.png'], 'Resolution', 200);
exportgraphics(figure_handle, [filename_stem '.pdf'], 'ContentType', 'vector');
close(figure_handle);
end
