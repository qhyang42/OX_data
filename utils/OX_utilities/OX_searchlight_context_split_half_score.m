function score = OX_searchlight_context_split_half_score( ...
        X, odor_labels, context_labels, run_ids, session_ids, ...
        neighborhoods, varargin)
%OX_SEARCHLIGHT_CONTEXT_SPLIT_HALF_SCORE Score many local feature sets.
%
%   SCORE = OX_SEARCHLIGHT_CONTEXT_SPLIT_HALF_SCORE(X, ODOR_LABELS,
%   CONTEXT_LABELS, RUN_IDS, SESSION_IDS, NEIGHBORHOODS, Name, Value, ...)
%
% X is trials-by-features. NEIGHBORHOODS is a cell vector whose entries are
% feature-column indices. The function constructs one common collection of
% session-balanced run splits, estimates all split-half context templates
% over the complete feature matrix, and then scores every neighborhood.
%
% The three returned context effects are, for each context, the Fisher-z
% diagonal similarity minus the mean of its four bidirectional
% off-diagonal similarities, averaged across repeated run splits.

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'X', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'odor_labels', @(x) isnumeric(x) && isvector(x));
addRequired(p, 'context_labels', ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addRequired(p, 'run_ids', @(x) isnumeric(x) && isvector(x));
addRequired(p, 'session_ids', ...
    @(x) isempty(x) || (isnumeric(x) && isvector(x)));
addRequired(p, 'neighborhoods', @(x) iscell(x) && isvector(x));
addParameter(p, 'NumSplits', 200, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == round(x));
addParameter(p, 'RandomSeed', 1, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'MaxSplitAttempts', 100000, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == round(x));
addParameter(p, 'ExpectedOdors', 20, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 2 && x == round(x));
addParameter(p, 'CorrelationClamp', 1e-7, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'CollectSplitDetails', true, ...
    @(x) islogical(x) && isscalar(x));
parse(p, X, odor_labels, context_labels, run_ids, session_ids, ...
    neighborhoods, varargin{:});
opts = p.Results;

X = double(X);
odor_labels = double(odor_labels(:));
context_labels = upper(strtrim(string(context_labels(:))));
run_ids = double(run_ids(:));
session_ids = double(session_ids(:));
n_trials = size(X, 1);
n_features = size(X, 2);
assert(n_features >= 2, 'X must contain at least two features.');
assert(numel(odor_labels) == n_trials && ...
    numel(context_labels) == n_trials && numel(run_ids) == n_trials, ...
    'Trial labels must match the rows of X.');
assert(isempty(session_ids) || numel(session_ids) == n_trials, ...
    'session_ids must be empty or match the rows of X.');
assert(all(isfinite(X), 'all'), 'X contains non-finite values.');
assert(all(isfinite(odor_labels)) && all(isfinite(run_ids)) && ...
    (isempty(session_ids) || all(isfinite(session_ids))), ...
    'Odor, run, and session labels must be finite.');

context_order = ["PERSON", "FOOD", "LOCATION"];
assert(all(ismember(context_labels, context_order)), ...
    'Only PERSON, FOOD, and LOCATION trials may be scored.');
[known_contexts, context_id] = ismember(context_labels, context_order);
odor_values = unique(odor_labels, 'sorted');
[known_odors, odor_id] = ismember(odor_labels, odor_values);
assert(all(known_contexts) && all(known_odors), ...
    'Unable to encode one or more trial labels.');
assert(numel(odor_values) == opts.ExpectedOdors, ...
    'Found %d odors; expected %d.', numel(odor_values), opts.ExpectedOdors);
for context_idx = 1:3
    if ~isequal(unique(odor_id(context_id == context_idx), 'sorted'), ...
            (1:opts.ExpectedOdors)')
        error('OX_DATA:IncompleteContextOdorCoverage', ...
            'Context %s does not contain all %d odors.', ...
            context_order(context_idx), opts.ExpectedOdors);
    end
end

runs = unique(run_ids, 'stable');
n_runs = numel(runs);
assert(mod(n_runs, 2) == 0, ...
    'An even number of runs is required; found %d.', n_runs);
assert(opts.MaxSplitAttempts >= opts.NumSplits, ...
    'MaxSplitAttempts must be at least NumSplits.');
[run_groups, group_values, stratified] = encode_run_groups( ...
    runs, run_ids, session_ids);

n_centers = numel(neighborhoods);
for center_idx = 1:n_centers
    indices = neighborhoods{center_idx};
    assert(isnumeric(indices) && isvector(indices) && ...
        all(isfinite(indices)) && all(indices == round(indices)) && ...
        all(indices >= 1 & indices <= n_features), ...
        'Neighborhood %d contains an invalid feature index.', center_idx);
    neighborhoods{center_idx} = unique(double(indices(:)'), 'stable');
end

[split_membership, cell_counts, split_diagnostics] = build_splits( ...
    run_ids, runs, run_groups, group_values, stratified, context_id, ...
    odor_id, opts);
[weights_a, weights_b] = make_template_weights(run_ids, runs, ...
    split_membership, context_id, odor_id, cell_counts, ...
    opts.ExpectedOdors);
templates_a = weights_a * X;
templates_b = weights_b * X;
context_mean_delta_z = score_neighborhoods(templates_a, templates_b, ...
    neighborhoods, opts.NumSplits, opts.CorrelationClamp);

score = struct();
score.context_order = cellstr(context_order(:));
score.context_mean_delta_z = context_mean_delta_z;
score.n_centers = n_centers;
score.n_features = n_features;
score.n_splits = opts.NumSplits;
score.n_runs = n_runs;
score.runs_per_half = n_runs / 2;
score.n_attempted_splits = split_diagnostics.n_attempted;
score.n_rejected_coverage = split_diagnostics.n_rejected_coverage;
score.n_rejected_duplicate = split_diagnostics.n_rejected_duplicate;
score.options = opts;
if opts.CollectSplitDetails
    score.split_run_membership = split_membership;
    score.half_a_run_ids = split_membership_to_runs( ...
        split_membership, runs, true);
    score.half_b_run_ids = split_membership_to_runs( ...
        split_membership, runs, false);
    score.cell_trial_counts = cell_counts;
else
    score.split_run_membership = [];
    score.half_a_run_ids = [];
    score.half_b_run_ids = [];
    score.cell_trial_counts = [];
end
end

function [run_groups, group_values, stratified] = encode_run_groups( ...
        runs, run_ids, session_ids)
n_runs = numel(runs);
if isempty(session_ids)
    run_groups = zeros(n_runs, 1);
    group_values = zeros(0, 1);
    stratified = false;
    return;
end
run_groups = zeros(n_runs, 1);
for run_idx = 1:n_runs
    values = unique(session_ids(run_ids == runs(run_idx)), 'stable');
    assert(isscalar(values), 'Run %g spans multiple sessions.', runs(run_idx));
    run_groups(run_idx) = values;
end
group_values = unique(run_groups, 'stable');
stratified = true;
end

function [split_membership, cell_counts, diagnostics] = build_splits( ...
        run_ids, runs, run_groups, group_values, stratified, context_id, ...
        odor_id, opts)
n_runs = numel(runs);
n_half_runs = n_runs / 2;
split_membership = false(opts.NumSplits, n_runs);
cell_counts = zeros(opts.NumSplits, 2, 3, opts.ExpectedOdors, 'uint16');
stream = RandStream('mt19937ar', 'Seed', ...
    mod(floor(opts.RandomSeed), 2^32));
n_valid = 0;
n_attempted = 0;
n_rejected_coverage = 0;
n_rejected_duplicate = 0;
while n_valid < opts.NumSplits && n_attempted < opts.MaxSplitAttempts
    n_attempted = n_attempted + 1;
    half_a = propose_half(stream, n_runs, run_groups, group_values, ...
        stratified);
    assert(nnz(half_a) == n_half_runs, ...
        'Internal error: proposed split has the wrong number of runs.');
    if ~half_a(1)
        half_a = ~half_a;
    end
    if n_valid > 0 && any(all( ...
            split_membership(1:n_valid, :) == half_a, 2))
        n_rejected_duplicate = n_rejected_duplicate + 1;
        continue;
    end
    [counts, coverage_ok] = count_split_cells(run_ids, runs, half_a, ...
        context_id, odor_id, opts.ExpectedOdors);
    if ~coverage_ok
        n_rejected_coverage = n_rejected_coverage + 1;
        continue;
    end
    n_valid = n_valid + 1;
    split_membership(n_valid, :) = half_a;
    cell_counts(n_valid, :, :, :) = counts;
end
if n_valid ~= opts.NumSplits
    error('OX_DATA:InsufficientValidSplits', ...
        ['Only %d/%d valid unique splits were found after %d attempts ' ...
         '(%d coverage rejections, %d duplicate rejections).'], ...
        n_valid, opts.NumSplits, n_attempted, n_rejected_coverage, ...
        n_rejected_duplicate);
end
diagnostics = struct('n_attempted', n_attempted, ...
    'n_rejected_coverage', n_rejected_coverage, ...
    'n_rejected_duplicate', n_rejected_duplicate);
end

function half_a = propose_half(stream, n_runs, run_groups, ...
        group_values, stratified)
half_a = false(1, n_runs);
if ~stratified
    half_a(randperm(stream, n_runs, n_runs / 2)) = true;
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
    'Internal error while balancing session groups.');
if n_extra > 0
    selected = odd_groups(randperm(stream, numel(odd_groups), n_extra));
    take_counts(selected) = take_counts(selected) + 1;
end
for group_idx = 1:n_groups
    indices = find(run_groups == group_values(group_idx));
    chosen = indices(randperm(stream, numel(indices), take_counts(group_idx)));
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
n_entries = n_splits * n_trials;
row_a = zeros(n_entries, 1);
column_a = zeros(n_entries, 1);
value_a = zeros(n_entries, 1);
row_b = zeros(n_entries, 1);
column_b = zeros(n_entries, 1);
value_b = zeros(n_entries, 1);
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
            'Internal error: an accepted split contains an empty cell.');

        destination = cursor_a + (1:numel(indices_a));
        row_a(destination) = template_row;
        column_a(destination) = indices_a;
        value_a(destination) = 1 ./ ...
            (n_odors * counts_a(odor_id(indices_a)));
        cursor_a = cursor_a + numel(indices_a);

        destination = cursor_b + (1:numel(indices_b));
        row_b(destination) = template_row;
        column_b(destination) = indices_b;
        value_b(destination) = 1 ./ ...
            (n_odors * counts_b(odor_id(indices_b)));
        cursor_b = cursor_b + numel(indices_b);
    end
end
weights_a = sparse(row_a(1:cursor_a), column_a(1:cursor_a), ...
    value_a(1:cursor_a), n_splits * 3, n_trials);
weights_b = sparse(row_b(1:cursor_b), column_b(1:cursor_b), ...
    value_b(1:cursor_b), n_splits * 3, n_trials);
end

function context_mean_delta_z = score_neighborhoods( ...
        templates_a, templates_b, neighborhoods, n_splits, clamp)
n_centers = numel(neighborhoods);
context_mean_delta_z = nan(n_centers, 3);
for center_idx = 1:n_centers
    feature_indices = neighborhoods{center_idx};
    pattern_a = templates_a(:, feature_indices);
    pattern_b = templates_b(:, feature_indices);
    pattern_a = pattern_a - mean(pattern_a, 2);
    pattern_b = pattern_b - mean(pattern_b, 2);
    norm_a = sqrt(sum(pattern_a .^ 2, 2));
    norm_b = sqrt(sum(pattern_b .^ 2, 2));
    if any(norm_a <= eps | ~isfinite(norm_a)) || ...
            any(norm_b <= eps | ~isfinite(norm_b))
        error('OX_DATA:InvalidTemplateVariance', ...
            'Searchlight %d has a zero-variance split-half template.', ...
            center_idx);
    end
    fisher_z = nan(n_splits, 3, 3);
    for context_a = 1:3
        rows_a = context_a:3:(3 * n_splits);
        for context_b = 1:3
            rows_b = context_b:3:(3 * n_splits);
            correlation = sum(pattern_a(rows_a, :) .* ...
                pattern_b(rows_b, :), 2) ./ ...
                (norm_a(rows_a) .* norm_b(rows_b));
            correlation = min(max(correlation, -1 + clamp), 1 - clamp);
            fisher_z(:, context_a, context_b) = atanh(correlation);
        end
    end
    for context_idx = 1:3
        other = setdiff(1:3, context_idx);
        diagonal = reshape(fisher_z(:, context_idx, context_idx), ...
            n_splits, 1);
        outgoing = reshape(fisher_z(:, context_idx, other), ...
            n_splits, 2);
        incoming = reshape(fisher_z(:, other, context_idx), ...
            n_splits, 2);
        context_mean_delta_z(center_idx, context_idx) = mean( ...
            diagonal - mean([outgoing, incoming], 2));
    end
end
end

function half_run_ids = split_membership_to_runs(membership, runs, use_a)
n_splits = size(membership, 1);
n_half_runs = size(membership, 2) / 2;
half_run_ids = nan(n_splits, n_half_runs);
for split_idx = 1:n_splits
    selected = membership(split_idx, :);
    if ~use_a
        selected = ~selected;
    end
    half_run_ids(split_idx, :) = runs(selected);
end
end
