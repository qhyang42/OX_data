function group_results = OX_group_context_split_half_similarity(subject_ids, varargin)
%OX_GROUP_CONTEXT_SPLIT_HALF_SIMILARITY Collect descriptive subject results.
%
% Loads subject-level cross-half context-similarity MAT files, writes one
% compact CSV, saves a full group MAT, and optionally creates descriptive
% QC panels. No test treats the repeated splits as independent samples.

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'subject_ids', @(x) isnumeric(x) && isvector(x) && ...
    all(isfinite(x)) && all(x == round(x)));
addParameter(p, 'MRIRoot', '', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'OutputDir', '', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'MakePlots', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'SaveOutputs', true, @(x) islogical(x) && isscalar(x));
parse(p, subject_ids, varargin{:});
opts = p.Results;
subject_ids = double(subject_ids(:));
assert(numel(unique(subject_ids)) == numel(subject_ids), ...
    'subject_ids contains duplicates.');

mri_root = resolve_mri_root(opts.MRIRoot);
if strlength(string(opts.OutputDir)) == 0
    output_dir = fullfile(mri_root, 'group', ...
        'context_split_half_similarity_physio');
else
    output_dir = char(string(opts.OutputDir));
end

n_subjects = numel(subject_ids);
subject_results = cell(n_subjects, 1);
summary_cells = cell(n_subjects, 1);
context_summary_cells = cell(n_subjects, 1);
for subject_idx = 1:n_subjects
    subject_file = fullfile(output_dir, sprintf( ...
        'context_split_half_similarity_subj%d_results.mat', ...
        subject_ids(subject_idx)));
    assert(isfile(subject_file), 'Missing subject result: %s', subject_file);
    loaded = load(subject_file, 'results');
    assert(isfield(loaded, 'results'), ...
        'Subject file lacks results: %s', subject_file);
    subject_results{subject_idx} = loaded.results;
    summary_cells{subject_idx} = loaded.results.summary;
    assert(isfield(loaded.results, 'context_summary'), ...
        'Subject file lacks context_summary: %s', subject_file);
    context_summary_cells{subject_idx} = loaded.results.context_summary;
end

reference_rois = string(subject_results{1}.summary.roi);
context_order = string(subject_results{1}.analysis.contexts);
n_rois = numel(reference_rois);
assert(numel(context_order) == 3, 'Expected three semantic contexts.');
mean_fisher_z_matrix = nan(n_subjects, n_rois, 3, 3);
fisher_mean_r_matrix = nan(n_subjects, n_rois, 3, 3);
mean_delta_z = nan(n_subjects, n_rois);
max_delta_z = nan(n_subjects, n_rois);
sd_delta_z = nan(n_subjects, n_rois);
n_permutations = subject_results{1}.null.n_permutations;
subject_null_mean_delta_z = nan(n_subjects, n_rois, n_permutations);
subject_null_max_delta_z = nan(n_subjects, n_rois, n_permutations);
subject_context_delta_z = nan(n_subjects, n_rois, 3);
subject_null_context_delta_z = nan(n_subjects, n_rois, 3, n_permutations);

for subject_idx = 1:n_subjects
    item = subject_results{subject_idx};
    assert(isequal(string(item.summary.roi), reference_rois), ...
        'ROI order differs for subject %d.', subject_ids(subject_idx));
    assert(isequal(string(item.analysis.contexts), context_order), ...
        'Context order differs for subject %d.', subject_ids(subject_idx));
    assert(item.null.n_permutations == n_permutations, ...
        'Permutation count differs for subject %d.', subject_ids(subject_idx));
    if n_permutations > 0
        subject_null_mean_delta_z(subject_idx, :, :) = ...
            item.null.roi_mean_delta_z;
        subject_null_max_delta_z(subject_idx, :, :) = ...
            item.null.roi_max_split_delta_z;
        subject_null_context_delta_z(subject_idx, :, :, :) = reshape( ...
            item.null.roi_context_mean_delta_z, ...
            [1, n_rois, 3, n_permutations]);
    end
    for roi_idx = 1:n_rois
        split = item.roi_results(roi_idx).split_half;
        if isempty(split)
            continue;
        end
        mean_fisher_z_matrix(subject_idx, roi_idx, :, :) = ...
            split.aggregate.mean_fisher_z_matrix;
        fisher_mean_r_matrix(subject_idx, roi_idx, :, :) = ...
            split.aggregate.fisher_mean_r_matrix;
        mean_delta_z(subject_idx, roi_idx) = split.aggregate.mean_delta_z;
        max_delta_z(subject_idx, roi_idx) = split.aggregate.max_delta_z;
        subject_context_delta_z(subject_idx, roi_idx, :) = ...
            split.aggregate.mean_context_delta_z;
        sd_delta_z(subject_idx, roi_idx) = split.aggregate.sd_delta_z;
    end
end

summary = vertcat(summary_cells{:});
context_summary = vertcat(context_summary_cells{:});
group_mean_fisher_z_matrix = reshape(mean( ...
    mean_fisher_z_matrix, 1, 'omitnan'), [n_rois, 3, 3]);
group_fisher_mean_r_matrix = tanh(group_mean_fisher_z_matrix);
group_observed_mean_delta_z = mean(mean_delta_z, 1, 'omitnan')';
group_observed_mean_max_delta_z = mean(max_delta_z, 1, 'omitnan')';
if n_permutations > 0
    group_null_mean_delta_z = reshape(mean( ...
        subject_null_mean_delta_z, 1, 'omitnan'), [n_rois, n_permutations]);
    group_null_mean_max_delta_z = reshape(mean( ...
        subject_null_max_delta_z, 1, 'omitnan'), [n_rois, n_permutations]);
    group_null_context_delta_z = reshape(mean( ...
        subject_null_context_delta_z, 1, 'omitnan'), [n_rois, 3, n_permutations]);
else
    group_null_mean_delta_z = nan(n_rois, 0);
    group_null_mean_max_delta_z = nan(n_rois, 0);
    group_null_context_delta_z = nan(n_rois, 3, 0);
end
group_p_mean_delta_z = nan(n_rois, 1);
group_p_max_delta_z = nan(n_rois, 1);
for roi_idx = 1:n_rois
    group_p_mean_delta_z(roi_idx) = empirical_upper_p( ...
        group_observed_mean_delta_z(roi_idx), ...
        group_null_mean_delta_z(roi_idx, :));
    group_p_max_delta_z(roi_idx) = empirical_upper_p( ...
        group_observed_mean_max_delta_z(roi_idx), ...
        group_null_mean_max_delta_z(roi_idx, :));
end
group_permutation_summary = table(reference_rois, ...
    group_observed_mean_delta_z, group_p_mean_delta_z, ...
    group_observed_mean_max_delta_z, group_p_max_delta_z, ...
    repmat(n_permutations, n_rois, 1), ...
    'VariableNames', {'roi', 'observed_group_mean_delta_z', ...
    'p_group_mean_delta_z', 'observed_group_mean_max_delta_z', ...
    'p_group_max_delta_z', 'n_permutations'});

group_observed_context_delta_z = reshape(mean( ...
    subject_context_delta_z, 1, 'omitnan'), [n_rois, 3]);
n_group_context_rows = n_rois * 3;
group_context_roi = strings(n_group_context_rows, 1);
group_context_name = strings(n_group_context_rows, 1);
group_context_delta_z = nan(n_group_context_rows, 1);
group_context_p = nan(n_group_context_rows, 1);
group_context_row = 0;
for roi_idx = 1:n_rois
    for context_idx = 1:3
        group_context_row = group_context_row + 1;
        group_context_roi(group_context_row) = reference_rois(roi_idx);
        group_context_name(group_context_row) = context_order(context_idx);
        group_context_delta_z(group_context_row) = ...
            group_observed_context_delta_z(roi_idx, context_idx);
        group_context_p(group_context_row) = empirical_upper_p( ...
            group_context_delta_z(group_context_row), squeeze( ...
            group_null_context_delta_z(roi_idx, context_idx, :)));
    end
end
group_context_permutation_summary = table(group_context_roi, ...
    group_context_name, group_context_delta_z, group_context_p, ...
    repmat(n_permutations, n_group_context_rows, 1), ...
    'VariableNames', {'roi', 'context', 'observed_group_context_delta_z', ...
    'p_group_context_delta_z', 'n_permutations'});

group_results = struct();
group_results.analysis = struct( ...
    'name', 'semantic_context_cross_half_pattern_similarity', ...
    'level', 'equal-weight subject mean within ROI', ...
    'inference', ['within-run semantic-context permutation; permutation k ' ...
        'averages one subject-null statistic from each subject separately ' ...
        'for each ROI; ROIs are never pooled']);
group_results.subject_ids = subject_ids;
group_results.roi_names = cellstr(reference_rois);
group_results.context_order = cellstr(context_order(:));
group_results.summary = summary;
group_results.context_summary = context_summary;
group_results.subject_mean_fisher_z_matrix = mean_fisher_z_matrix;
group_results.subject_fisher_mean_r_matrix = fisher_mean_r_matrix;
group_results.subject_mean_delta_z = mean_delta_z;
group_results.subject_max_delta_z = max_delta_z;
group_results.subject_context_delta_z = subject_context_delta_z;
group_results.subject_sd_delta_z_across_splits = sd_delta_z;
group_results.group_mean_fisher_z_matrix = group_mean_fisher_z_matrix;
group_results.group_fisher_mean_r_matrix = group_fisher_mean_r_matrix;
group_results.group_permutation_summary = group_permutation_summary;
group_results.group_context_permutation_summary = ...
    group_context_permutation_summary;
group_results.null = struct( ...
    'n_permutations', n_permutations, ...
    'subject_roi_mean_delta_z', subject_null_mean_delta_z, ...
    'subject_roi_max_delta_z', subject_null_max_delta_z, ...
    'subject_roi_context_delta_z', subject_null_context_delta_z, ...
    'group_roi_mean_delta_z', group_null_mean_delta_z, ...
    'group_roi_mean_max_delta_z', group_null_mean_max_delta_z, ...
    'group_roi_context_delta_z', group_null_context_delta_z, ...
    'roi_pooling', false, ...
    'group_construction', ['for permutation k and ROI r, average the five ' ...
        'subject-specific null statistics at k for r']);
group_results.subject_results = subject_results;
group_results.dimension_order = struct( ...
    'subject_matrix', {{'subject', 'roi', 'context_a', 'context_b'}}, ...
    'group_matrix', {{'roi', 'context_a', 'context_b'}}, ...
    'delta', {{'subject', 'roi'}}, ...
    'subject_null', {{'subject', 'roi', 'permutation'}}, ...
    'group_null', {{'roi', 'permutation'}}, ...
    'subject_context_null', ...
        {{'subject', 'roi', 'context', 'permutation'}}, ...
    'group_context_null', {{'roi', 'context', 'permutation'}});
group_results.options = opts;
group_results.output_dir = output_dir;
group_results.matlab_version = version;

if opts.SaveOutputs && ~isfolder(output_dir)
    mkdir(output_dir);
end

if opts.MakePlots
    qc_file = fullfile(output_dir, ...
        'context_split_half_similarity_qc.png');
    fig_file = fullfile(output_dir, ...
        'context_split_half_similarity_qc.fig');
    make_qc_plot(group_fisher_mean_r_matrix, mean_delta_z, ...
        subject_ids, reference_rois, context_order, qc_file, fig_file, ...
        opts.SaveOutputs);
    if opts.SaveOutputs
        group_results.qc_plot_file = qc_file;
        group_results.qc_figure_file = fig_file;
    end
end

if opts.SaveOutputs
    csv_file = fullfile(output_dir, ...
        'context_split_half_similarity_summary.csv');
    group_csv_file = fullfile(output_dir, ...
        'context_split_half_similarity_group_permutation_summary.csv');
    subject_context_csv_file = fullfile(output_dir, ...
        'context_split_half_similarity_context_subject_summary.csv');
    group_context_csv_file = fullfile(output_dir, ...
        'context_split_half_similarity_context_group_summary.csv');
    mat_file = fullfile(output_dir, ...
        'context_split_half_similarity_group_results.mat');
    writetable(summary, csv_file);
    writetable(group_permutation_summary, group_csv_file);
    writetable(context_summary, subject_context_csv_file);
    writetable(group_context_permutation_summary, group_context_csv_file);
    save(mat_file, 'group_results', '-v7.3');
    fprintf('Saved group result to %s\n', mat_file);
    fprintf('Saved compact summary to %s\n', csv_file);
    fprintf('Saved group permutation summary to %s\n', group_csv_file);
    fprintf('Saved context subject summary to %s\n', ...
        subject_context_csv_file);
    fprintf('Saved context group summary to %s\n', group_context_csv_file);
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

function make_qc_plot(mean_r, mean_delta, subject_ids, roi_names, ...
        context_order, png_file, fig_file, save_outputs)
n_rois = numel(roi_names);
figure_width = max(1200, 320 * n_rois);
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, figure_width, 650]);
cleanup = onCleanup(@() close(fig));
tiledlayout(2, n_rois, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
finite_r = mean_r(isfinite(mean_r));
if isempty(finite_r)
    color_limit = 1;
else
    color_limit = max(0.1, max(abs(finite_r)));
end

for roi_idx = 1:n_rois
    nexttile(roi_idx);
    imagesc(squeeze(mean_r(roi_idx, :, :)));
    axis image
    clim([-color_limit, color_limit]);
    colormap(gca, parula);
    colorbar;
    xticks(1:3);
    yticks(1:3);
    xticklabels(context_order);
    yticklabels(context_order);
    xtickangle(35);
    title(sprintf('%s: Fisher mean r', roi_names(roi_idx)), ...
        'Interpreter', 'none');
    xlabel('Half B');
    ylabel('Half A');

    nexttile(n_rois + roi_idx);
    values = mean_delta(:, roi_idx);
    plot(subject_ids, values, 'o-', 'LineWidth', 1.2, ...
        'MarkerFaceColor', [0.2, 0.45, 0.8], ...
        'Color', [0.2, 0.45, 0.8]);
    hold on
    yline(0, 'k:');
    xlim([min(subject_ids) - 0.5, max(subject_ids) + 0.5]);
    xticks(subject_ids);
    xlabel('Subject');
    ylabel('\Delta z');
    title(sprintf('%s: diagonal - off-diagonal', roi_names(roi_idx)), ...
        'Interpreter', 'none');
    grid on
end
if save_outputs
    exportgraphics(fig, png_file, 'Resolution', 180);
    savefig(fig, fig_file);
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
