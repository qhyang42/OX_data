function similarity = OX_context_split_half_similarity( ...
        X, odor_labels, context_labels, run_ids, varargin)
%OX_CONTEXT_SPLIT_HALF_SIMILARITY Cross-half semantic-context similarity.
%
%   similarity = OX_context_split_half_similarity(X, odor_labels, ...
%       context_labels, run_ids, Name, Value, ...)
%
% X is trials-by-voxels. Repeated disjoint run splits are used to build
% equally odor-weighted PERSON, FOOD, and LOCATION templates in each half.
% The returned 3-by-3 matrices contain Pearson correlations across voxels
% between the half-A and half-B templates. CONTROL rows, when supplied, are
% excluded before split construction and template estimation.
%
% Relevant options:
%   NumSplits          Number of valid, unique splits (default 200).
%   RandomSeed         Fixed split-sampling seed (default 1).
%   StratifyGroups     Optional trial-level acquisition/session labels.
%                      Each run must belong to one group; halves differ by
%                      at most one run within every group.
%   MaxSplitAttempts   Maximum proposals before failure (default 100000).
%   CorrelationClamp   Clamp r to +/- (1-CorrelationClamp) (default 1e-7).
%   CollectDetails     Retain all split matrices/run details (default true).

% No inferential test is performed across repeated splits because the
% splits reuse the same observations.

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'X', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'odor_labels', @(x) isnumeric(x) && isvector(x));
addRequired(p, 'context_labels', ...
    @(x) iscellstr(x) || ischar(x) || isstring(x));
addRequired(p, 'run_ids', @(x) isnumeric(x) && isvector(x));
addParameter(p, 'NumSplits', 200, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == round(x));
addParameter(p, 'RandomSeed', 1, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'StratifyGroups', [], ...
    @(x) isempty(x) || (isvector(x) && (isnumeric(x) || islogical(x) || ...
    isstring(x) || iscellstr(x) || iscategorical(x))));
addParameter(p, 'MaxSplitAttempts', 100000, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == round(x));
addParameter(p, 'ExpectedOdors', 20, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 2 && x == round(x));
addParameter(p, 'CorrelationClamp', 1e-7, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x < 1);
addParameter(p, 'CollectDetails', true, ...
    @(x) islogical(x) && isscalar(x));
parse(p, X, odor_labels, context_labels, run_ids, varargin{:});
opts = p.Results;

X = double(X);
odor_labels = double(odor_labels(:));
context_labels = upper(strtrim(string(context_labels(:))));
run_ids = double(run_ids(:));
n_input_trials = size(X, 1);
assert(numel(odor_labels) == n_input_trials && ...
    numel(context_labels) == n_input_trials && numel(run_ids) == n_input_trials, ...
    'Labels and run IDs must match the rows of X.');
assert(size(X, 2) >= 2, 'At least two voxels are required.');
assert(all(isfinite(X), 'all'), 'X contains non-finite values.');
assert(all(isfinite(odor_labels)) && all(isfinite(run_ids)), ...
    'Odor labels and run IDs must be finite.');

semantic_contexts = ["PERSON", "FOOD", "LOCATION"];
allowed_contexts = [semantic_contexts, "CONTROL"];
assert(all(ismember(context_labels, allowed_contexts)), ...
    'Context labels must be PERSON, FOOD, LOCATION, or CONTROL.');

if isempty(opts.StratifyGroups)
    stratify_groups = strings(n_input_trials, 1);
    use_stratification = false;
else
    assert(numel(opts.StratifyGroups) == n_input_trials, ...
        'StratifyGroups must match the rows of X.');
    stratify_groups = string(opts.StratifyGroups(:));
    assert(all(~ismissing(stratify_groups)), ...
        'StratifyGroups cannot contain missing values.');
    use_stratification = true;
end

% CONTROL contributes neither observations nor preprocessing here.
semantic_mask = ismember(context_labels, semantic_contexts);
X = X(semantic_mask, :);
odor_labels = odor_labels(semantic_mask);
context_labels = context_labels(semantic_mask);
run_ids = run_ids(semantic_mask);
stratify_groups = stratify_groups(semantic_mask);

odor_values = unique(odor_labels, 'sorted');
assert(numel(odor_values) == opts.ExpectedOdors, ...
    'Found %d semantic odor identities; expected %d.', ...
    numel(odor_values), opts.ExpectedOdors);
[known_odors, odor_id] = ismember(odor_labels, odor_values);
[known_contexts, context_id] = ismember(context_labels, semantic_contexts);
assert(all(known_odors) && all(known_contexts), ...
    'Unable to encode one or more semantic trials.');

for context_idx = 1:3
    assert(isequal(unique(odor_id(context_id == context_idx), 'sorted'), ...
        (1:opts.ExpectedOdors)'), ...
        'Context %s does not contain all %d odors.', ...
        semantic_contexts(context_idx), opts.ExpectedOdors);
end

runs = unique(run_ids, 'stable');
n_runs = numel(runs);
assert(mod(n_runs, 2) == 0, ...
    'An even number of runs is required; found %d.', n_runs);
n_half_runs = n_runs / 2;
assert(opts.MaxSplitAttempts >= opts.NumSplits, ...
    'MaxSplitAttempts must be at least NumSplits.');

run_groups = strings(n_runs, 1);
if use_stratification
    for run_idx = 1:n_runs
        values = unique(stratify_groups(run_ids == runs(run_idx)), 'stable');
        assert(isscalar(values), ...
            'Run %g belongs to more than one StratifyGroups value.', runs(run_idx));
        run_groups(run_idx) = values;
    end
    group_values = unique(run_groups, 'stable');
else
    group_values = strings(0, 1);
end

previous_rng = rng;
rng_cleanup = onCleanup(@() rng(previous_rng));
rng(opts.RandomSeed, 'twister');

raw_r = nan(3, 3, opts.NumSplits);
fisher_z = nan(3, 3, opts.NumSplits);
mean_diagonal_z = nan(opts.NumSplits, 1);
mean_off_diagonal_z = nan(opts.NumSplits, 1);
delta_z = nan(opts.NumSplits, 1);
context_delta_z = nan(opts.NumSplits, 3);
half_a_run_ids = nan(opts.NumSplits, n_half_runs);
half_b_run_ids = nan(opts.NumSplits, n_half_runs);
split_membership = false(opts.NumSplits, n_runs);
cell_trial_counts = zeros(opts.NumSplits, 2, 3, opts.ExpectedOdors, 'uint16');

n_valid = 0;
n_attempted = 0;
n_rejected_coverage = 0;
n_rejected_duplicate = 0;
diagonal_mask = eye(3) == 1;
off_diagonal_mask = ~diagonal_mask;

while n_valid < opts.NumSplits && n_attempted < opts.MaxSplitAttempts
    n_attempted = n_attempted + 1;
    half_a_membership = propose_half(n_runs, run_groups, group_values, ...
        use_stratification);
    assert(nnz(half_a_membership) == n_half_runs, ...
        'Internal error: proposed half does not contain %d runs.', n_half_runs);

    % A/B swaps carry the same delta. Canonicalizing avoids duplicate
    % complementary partitions and makes uniqueness checks deterministic.
    if ~half_a_membership(1)
        half_a_membership = ~half_a_membership;
    end
    if n_valid > 0 && any(all(split_membership(1:n_valid, :) == ...
            half_a_membership, 2))
        n_rejected_duplicate = n_rejected_duplicate + 1;
        continue;
    end

    [counts, coverage_ok] = count_split_cells(run_ids, runs, ...
        half_a_membership, context_id, odor_id, opts.ExpectedOdors);
    if ~coverage_ok
        n_rejected_coverage = n_rejected_coverage + 1;
        continue;
    end

    n_valid = n_valid + 1;
    split_membership(n_valid, :) = half_a_membership;
    half_a_run_ids(n_valid, :) = runs(half_a_membership);
    half_b_run_ids(n_valid, :) = runs(~half_a_membership);
    cell_trial_counts(n_valid, :, :, :) = counts;

end

if n_valid ~= opts.NumSplits
    error('OX_DATA:InsufficientValidSplits', ...
        ['Only %d/%d valid unique splits were found after %d attempts ' ...
         '(%d coverage rejections, %d duplicate rejections).'], ...
        n_valid, opts.NumSplits, n_attempted, n_rejected_coverage, ...
        n_rejected_duplicate);
end

% Estimate all templates with two sparse weight-matrix multiplications.
% Each trial weight is 1/(20 * repetitions for its odor in that cell), so
% this is exactly the specified within-odor mean followed by equal odor mean.
[half_a_weights, half_b_weights] = make_template_weights(run_ids, runs, ...
    split_membership, context_id, odor_id, cell_trial_counts, ...
    opts.ExpectedOdors);
half_a_templates = half_a_weights * X;
half_b_templates = half_b_weights * X;
for split_idx = 1:opts.NumSplits
    template_rows = (split_idx - 1) * 3 + (1:3);
    template_a = half_a_templates(template_rows, :);
    template_b = half_b_templates(template_rows, :);
    split_r = correlate_template_rows(template_a, template_b);
    split_r = min(max(split_r, -1 + opts.CorrelationClamp), ...
        1 - opts.CorrelationClamp);
    split_z = atanh(split_r);

    raw_r(:, :, split_idx) = split_r;
    fisher_z(:, :, split_idx) = split_z;
    mean_diagonal_z(split_idx) = mean(split_z(diagonal_mask));
    mean_off_diagonal_z(split_idx) = mean(split_z(off_diagonal_mask));
    delta_z(split_idx) = mean_diagonal_z(split_idx) - ...
        mean_off_diagonal_z(split_idx);
    for context_idx = 1:3
        other_contexts = setdiff(1:3, context_idx);
        context_off_diagonal = [split_z(context_idx, other_contexts), ...
            split_z(other_contexts, context_idx)'];
        assert(numel(context_off_diagonal) == 4, ...
            'A context-specific contrast must contain four off-diagonal cells.');
        context_delta_z(split_idx, context_idx) = ...
            split_z(context_idx, context_idx) - mean(context_off_diagonal);
    end
end

mean_fisher_z_matrix = mean(fisher_z, 3);
assert(abs(mean(mean(context_delta_z, 1)) - mean(delta_z)) < 1e-10, ...
    'Context-specific effects do not reconstruct the overall delta-z.');
similarity = struct();
similarity.context_order = cellstr(semantic_contexts(:));
similarity.odor_values = odor_values;
similarity.n_input_trials = n_input_trials;
similarity.n_semantic_trials = size(X, 1);
similarity.n_runs = n_runs;
similarity.runs_per_half = n_half_runs;
similarity.n_valid_splits = n_valid;
similarity.n_attempted_splits = n_attempted;
similarity.n_rejected_coverage = n_rejected_coverage;
similarity.n_rejected_duplicate = n_rejected_duplicate;
similarity.raw_r = raw_r;
similarity.fisher_z = fisher_z;
similarity.mean_diagonal_z = mean_diagonal_z;
similarity.mean_off_diagonal_z = mean_off_diagonal_z;
similarity.delta_z = delta_z;
similarity.context_delta_z = context_delta_z;
similarity.half_a_run_ids = half_a_run_ids;
similarity.half_b_run_ids = half_b_run_ids;
similarity.split_run_membership = split_membership;
similarity.cell_trial_counts = cell_trial_counts;
similarity.cell_trial_count_dimension_order = ...
    {'split', 'half', 'context', 'odor'};
similarity.stratified = use_stratification;
similarity.run_ids = runs;
similarity.run_stratification_groups = run_groups;
similarity.aggregate = struct( ...
    'mean_diagonal_z', mean(mean_diagonal_z), ...
    'mean_off_diagonal_z', mean(mean_off_diagonal_z), ...
    'mean_delta_z', mean(delta_z), ...
    'max_delta_z', max(delta_z), ...
    'mean_context_delta_z', mean(context_delta_z, 1), ...
    'max_context_delta_z', max(context_delta_z, [], 1), ...
    'sd_delta_z', std(delta_z, 0), ...
    'mean_raw_r_matrix', mean(raw_r, 3), ...
    'mean_fisher_z_matrix', mean_fisher_z_matrix, ...
    'fisher_mean_r_matrix', tanh(mean_fisher_z_matrix));
similarity.options = opts;
similarity.inference = ['descriptive repeated splits only; no t-test or ' ...
    'other inferential test across dependent splits'];
if ~opts.CollectDetails
    similarity.raw_r = [];
    similarity.fisher_z = [];
    similarity.mean_diagonal_z = [];
    similarity.mean_off_diagonal_z = [];
    similarity.delta_z = [];
    similarity.context_delta_z = [];
    similarity.half_a_run_ids = [];
    similarity.half_b_run_ids = [];
    similarity.split_run_membership = [];
    similarity.cell_trial_counts = [];
end
end

function half_a = propose_half(n_runs, run_groups, group_values, stratified)
half_a = false(1, n_runs);
if ~stratified
    half_a(randperm(n_runs, n_runs / 2)) = true;
    return;
end

n_groups = numel(group_values);
group_sizes = zeros(n_groups, 1);
take_counts = zeros(n_groups, 1);
for group_idx = 1:n_groups
    group_sizes(group_idx) = nnz(run_groups == group_values(group_idx));
    take_counts(group_idx) = floor(group_sizes(group_idx) / 2);
end
n_extra = n_runs / 2 - sum(take_counts);
odd_groups = find(mod(group_sizes, 2) == 1);
assert(n_extra >= 0 && n_extra <= numel(odd_groups), ...
    'Internal error while balancing stratification groups.');
if n_extra > 0
    selected = odd_groups(randperm(numel(odd_groups), n_extra));
    take_counts(selected) = take_counts(selected) + 1;
end
for group_idx = 1:n_groups
    indices = find(run_groups == group_values(group_idx));
    chosen = indices(randperm(numel(indices), take_counts(group_idx)));
    half_a(chosen) = true;
end
end

function [counts, valid] = count_split_cells(run_ids, runs, half_a, ...
        context_id, odor_id, n_odors)
counts = zeros(2, 3, n_odors, 'uint16');
for half_idx = 1:2
    if half_idx == 1
        half_runs = runs(half_a);
    else
        half_runs = runs(~half_a);
    end
    in_half = ismember(run_ids, half_runs);
    for context_idx = 1:3
        for odor_idx = 1:n_odors
            counts(half_idx, context_idx, odor_idx) = uint16(nnz( ...
                in_half & context_id == context_idx & odor_id == odor_idx));
        end
    end
end
valid = all(counts(:) > 0);
end

function [weights_a, weights_b] = make_template_weights(run_ids, runs, ...
        split_membership, context_id, odor_id, cell_counts, n_odors)
n_splits = size(split_membership, 1);
n_trials = numel(run_ids);
n_entries_per_half = n_splits * n_trials;
row_a = zeros(n_entries_per_half, 1);
column_a = zeros(n_entries_per_half, 1);
value_a = zeros(n_entries_per_half, 1);
row_b = zeros(n_entries_per_half, 1);
column_b = zeros(n_entries_per_half, 1);
value_b = zeros(n_entries_per_half, 1);
cursor_a = 0;
cursor_b = 0;
for split_idx = 1:n_splits
    in_a = ismember(run_ids, runs(split_membership(split_idx, :)));
    for context_idx = 1:3
        template_row = (split_idx - 1) * 3 + context_idx;
        indices_a = find(in_a & context_id == context_idx);
        indices_b = find(~in_a & context_id == context_idx);
        counts_a = double(reshape(cell_counts( ...
            split_idx, 1, context_idx, :), n_odors, 1));
        counts_b = double(reshape(cell_counts( ...
            split_idx, 2, context_idx, :), n_odors, 1));
        assert(all(counts_a > 0) && all(counts_b > 0), ...
            'Internal error: template weights found an empty cell.');

        destination_a = cursor_a + (1:numel(indices_a));
        row_a(destination_a) = template_row;
        column_a(destination_a) = indices_a;
        value_a(destination_a) = 1 ./ ...
            (n_odors * counts_a(odor_id(indices_a)));
        cursor_a = cursor_a + numel(indices_a);

        destination_b = cursor_b + (1:numel(indices_b));
        row_b(destination_b) = template_row;
        column_b(destination_b) = indices_b;
        value_b(destination_b) = 1 ./ ...
            (n_odors * counts_b(odor_id(indices_b)));
        cursor_b = cursor_b + numel(indices_b);
    end
end
row_a = row_a(1:cursor_a);
column_a = column_a(1:cursor_a);
value_a = value_a(1:cursor_a);
row_b = row_b(1:cursor_b);
column_b = column_b(1:cursor_b);
value_b = value_b(1:cursor_b);
weights_a = sparse(row_a, column_a, value_a, n_splits * 3, n_trials);
weights_b = sparse(row_b, column_b, value_b, n_splits * 3, n_trials);
end

function correlations = correlate_template_rows(template_a, template_b)
template_a = template_a - mean(template_a, 2);
template_b = template_b - mean(template_b, 2);
norm_a = sqrt(sum(template_a .^ 2, 2));
norm_b = sqrt(sum(template_b .^ 2, 2));
assert(all(norm_a > eps & isfinite(norm_a)), ...
    'A half-A template has zero or invalid spatial variance.');
assert(all(norm_b > eps & isfinite(norm_b)), ...
    'A half-B template has zero or invalid spatial variance.');
correlations = (template_a ./ norm_a) * (template_b ./ norm_b)';
end
