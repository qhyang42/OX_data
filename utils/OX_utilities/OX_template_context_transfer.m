function decoded = OX_template_context_transfer(X, odor_labels, context_labels, run_ids, varargin)
%OX_TEMPLATE_CONTEXT_TRANSFER Run semantic and control odor-template transfer.
%
%   decoded = OX_template_context_transfer(X, odor_labels, context_labels,
%       run_ids, Name, Value, ...)
%
% X is trials-by-features. For each held-out run, odor templates are built
% from one source context using every other run. Semantic source contexts
% are tested on PERSON, FOOD, and LOCATION. CONTROL templates are tested on
% those three contexts and on held-out CONTROL trials.

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'X', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'odor_labels', @(x) isnumeric(x) && isvector(x));
addRequired(p, 'context_labels', @(x) iscellstr(x) || ischar(x) || isstring(x));
addRequired(p, 'run_ids', @(x) isnumeric(x) && isvector(x));
addParameter(p, 'MinVoxels', 10, @(x) isnumeric(x) && isscalar(x) && x >= 2 && x == round(x));
addParameter(p, 'DemeanPatterns', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'CollectDetails', true, @(x) islogical(x) && isscalar(x));
parse(p, X, odor_labels, context_labels, run_ids, varargin{:});
opts = p.Results;

X = double(X);
odor_labels = double(odor_labels(:));
context_labels = upper(strtrim(string(context_labels(:))));
run_ids = double(run_ids(:));
n_trials = size(X, 1);
assert(numel(odor_labels) == n_trials && numel(context_labels) == n_trials && ...
    numel(run_ids) == n_trials, 'Labels and run IDs must match the rows of X.');
assert(all(isfinite(X), 'all'), 'X contains non-finite values.');
assert(all(isfinite(odor_labels)) && all(isfinite(run_ids)), ...
    'Odor labels and run IDs must be finite.');
assert(size(X, 2) >= opts.MinVoxels, ...
    'X has %d features; MinVoxels is %d.', size(X, 2), opts.MinVoxels);

context_order = ["PERSON", "FOOD", "LOCATION", "CONTROL"];
assert(all(ismember(context_labels, context_order)), 'Unknown context label found.');
[known_context, context_id] = ismember(context_labels, context_order);
assert(all(known_context), 'Unable to encode one or more context labels.');
odor_values = unique(odor_labels, 'sorted');
[known_odor, odor_id] = ismember(odor_labels, odor_values);
assert(all(known_odor), 'Unable to encode one or more odor labels.');
n_odors = numel(odor_values);
assert(n_odors >= 2, 'At least two odor identities are required.');

if opts.DemeanPatterns
    X = X - mean(X, 2);
end

runs = unique(run_ids, 'stable');
n_runs = numel(runs);
evidence_matrix = nan(4, 4);
accuracy_matrix = nan(4, 4);
trial_evidence = nan(4, 4, n_trials);
trial_predictions = nan(4, 4, n_trials);
fold_evidence = nan(4, 4, n_runs);
fold_accuracy = nan(4, 4, n_runs);

for source_context = 1:4
    if source_context == 4
        target_contexts = 1:4;
    else
        target_contexts = 1:3;
    end

    for fold_idx = 1:n_runs
        held_out_run = runs(fold_idx);
        train_mask = run_ids ~= held_out_run & context_id == source_context;
        assert(~any(train_mask & run_ids == held_out_run), ...
            'Held-out run %g leaked into a training template.', held_out_run);
        training_classes = unique(odor_id(train_mask), 'sorted');
        assert(isequal(training_classes, (1:n_odors)'), ...
            ['Held-out run %g leaves source context %s without all %d odors. ' ...
             'No test-run trial may be added to repair this fold.'], ...
            held_out_run, context_order(source_context), n_odors);
        templates = make_templates(X(train_mask, :), odor_id(train_mask), n_odors);

        for target_context = target_contexts
            test_mask = run_ids == held_out_run & context_id == target_context;
            if ~any(test_mask)
                continue;
            end
            test_indices = find(test_mask);
            [predictions, evidence] = score_templates( ...
                X(test_mask, :), templates, odor_id(test_mask));
            trial_predictions(source_context, target_context, test_indices) = ...
                reshape(predictions, 1, 1, []);
            trial_evidence(source_context, target_context, test_indices) = ...
                reshape(evidence, 1, 1, []);
            fold_accuracy(source_context, target_context, fold_idx) = ...
                mean(predictions == odor_id(test_mask));
            fold_evidence(source_context, target_context, fold_idx) = mean(evidence);
        end
    end

    for target_context = target_contexts
        eligible = context_id == target_context;
        cell_predictions = squeeze(trial_predictions(source_context, target_context, eligible));
        cell_evidence = squeeze(trial_evidence(source_context, target_context, eligible));
        assert(all(isfinite(cell_predictions)) && all(isfinite(cell_evidence)), ...
            'One or more eligible test trials were not decoded for %s -> %s.', ...
            context_order(source_context), context_order(target_context));
        accuracy_matrix(source_context, target_context) = ...
            mean(cell_predictions(:) == odor_id(eligible));
        evidence_matrix(source_context, target_context) = mean(cell_evidence(:));
    end
end

semantic_evidence = evidence_matrix(1:3, 1:3);
semantic_accuracy = accuracy_matrix(1:3, 1:3);
diagonal_mask = eye(3) == 1;
off_diagonal_mask = ~diagonal_mask;

decoded = struct();
decoded.context_order = cellstr(context_order(:));
decoded.semantic_contexts = cellstr(context_order(1:3)');
decoded.odor_values = odor_values;
decoded.chance_accuracy = 1 / n_odors;
decoded.semantic = struct('evidence', semantic_evidence, 'accuracy', semantic_accuracy);
decoded.primary = struct( ...
    'within_context_evidence', mean(semantic_evidence(diagonal_mask)), ...
    'cross_context_evidence', mean(semantic_evidence(off_diagonal_mask)), ...
    'evidence_effect', mean(semantic_evidence(diagonal_mask)) - ...
        mean(semantic_evidence(off_diagonal_mask)), ...
    'within_context_accuracy', mean(semantic_accuracy(diagonal_mask)), ...
    'cross_context_accuracy', mean(semantic_accuracy(off_diagonal_mask)), ...
    'accuracy_effect', mean(semantic_accuracy(diagonal_mask)) - ...
        mean(semantic_accuracy(off_diagonal_mask)));
decoded.control_transfer = struct( ...
    'evidence', evidence_matrix(4, 1:3), ...
    'accuracy', accuracy_matrix(4, 1:3), ...
    'mean_evidence', mean(evidence_matrix(4, 1:3)), ...
    'mean_accuracy', mean(accuracy_matrix(4, 1:3)));
decoded.control_within = struct( ...
    'evidence', evidence_matrix(4, 4), ...
    'accuracy', accuracy_matrix(4, 4));
decoded.all_evidence = evidence_matrix;
decoded.all_accuracy = accuracy_matrix;

if opts.CollectDetails
    decoded.trial_predictions = trial_predictions;
    decoded.trial_evidence = trial_evidence;
    decoded.fold_accuracy = fold_accuracy;
    decoded.fold_evidence = fold_evidence;
    decoded.detail_dimension_order = struct( ...
        'trial', {{'source_context', 'target_context', 'trial'}}, ...
        'fold', {{'source_context', 'target_context', 'fold'}});
else
    decoded.trial_predictions = [];
    decoded.trial_evidence = [];
    decoded.fold_accuracy = [];
    decoded.fold_evidence = [];
end
end

function templates = make_templates(Xtrain, ytrain, n_odors)
templates = nan(n_odors, size(Xtrain, 2));
for odor_idx = 1:n_odors
    templates(odor_idx, :) = mean(Xtrain(ytrain == odor_idx, :), 1);
end
assert(all(isfinite(templates), 'all'), 'A template contains non-finite values.');
end

function [predictions, evidence] = score_templates(Xtest, templates, ytest)
Xtest = Xtest - mean(Xtest, 2);
templates = templates - mean(templates, 2);
test_norm = sqrt(sum(Xtest.^2, 2));
template_norm = sqrt(sum(templates.^2, 2));
assert(all(test_norm > eps & isfinite(test_norm)), ...
    'A test pattern has zero or invalid spatial variance.');
assert(all(template_norm > eps & isfinite(template_norm)), ...
    'A template has zero or invalid spatial variance.');
similarity = (Xtest ./ test_norm) * (templates ./ template_norm)';
[~, predictions] = max(similarity, [], 2);

n_tests = size(similarity, 1);
n_odors = size(similarity, 2);
correct_linear = sub2ind(size(similarity), (1:n_tests)', ytest(:));
correct_similarity = similarity(correct_linear);
incorrect_mean = (sum(similarity, 2) - correct_similarity) / (n_odors - 1);
evidence = correct_similarity - incorrect_mean;
end
