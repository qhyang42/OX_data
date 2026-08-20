function decoded = OX_template_cross_odor_context_transfer(X, odor_labels, context_labels, run_ids, varargin)
%OX_TEMPLATE_CROSS_ODOR_CONTEXT_TRANSFER Decode context in held-out odors.
%
%   decoded = OX_template_cross_odor_context_transfer(X, odor_labels,
%       context_labels, run_ids, Name, Value, ...)
%
% PERSON, FOOD, and LOCATION trials are decoded. For every target odor and
% held-out run, three context templates are formed from all other odors and
% all other runs. CONTROL trials never enter training, testing, or summary
% statistics.

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

semantic_contexts = ["PERSON", "FOOD", "LOCATION"];
known_contexts = [semantic_contexts, "CONTROL"];
assert(all(ismember(context_labels, known_contexts)), 'Unknown context label found.');
semantic_mask = ismember(context_labels, semantic_contexts);
assert(any(semantic_mask) && any(context_labels == "CONTROL"), ...
    'Both semantic and CONTROL trials are expected in the input metadata.');
[known_semantic, context_id] = ismember(context_labels, semantic_contexts);
context_id(~known_semantic) = 0;

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
trial_predictions = nan(n_trials, 1);
trial_evidence = nan(n_trials, 1);
fold_accuracy = nan(n_odors, n_runs);
fold_evidence = nan(n_odors, n_runs);

for odor_idx = 1:n_odors
    for fold_idx = 1:n_runs
        held_out_run = runs(fold_idx);
        train_mask = semantic_mask & odor_id ~= odor_idx & run_ids ~= held_out_run;
        assert(~any(train_mask & odor_id == odor_idx), ...
            'Target odor %g leaked into a training template.', odor_values(odor_idx));
        assert(~any(train_mask & run_ids == held_out_run), ...
            'Held-out run %g leaked into a training template.', held_out_run);
        training_classes = unique(context_id(train_mask), 'sorted');
        assert(isequal(training_classes, (1:3)'), ...
            'Training fold for target odor %g and run %g lacks a semantic context.', ...
            odor_values(odor_idx), held_out_run);

        test_mask = semantic_mask & odor_id == odor_idx & run_ids == held_out_run;
        if ~any(test_mask)
            continue;
        end
        templates = make_context_templates(X(train_mask, :), context_id(train_mask));
        [predictions, evidence] = score_context_templates( ...
            X(test_mask, :), templates, context_id(test_mask));
        trial_predictions(test_mask) = predictions;
        trial_evidence(test_mask) = evidence;
        fold_accuracy(odor_idx, fold_idx) = mean(predictions == context_id(test_mask));
        fold_evidence(odor_idx, fold_idx) = mean(evidence);
    end
end

assert(all(isfinite(trial_predictions(semantic_mask))) && ...
    all(isfinite(trial_evidence(semantic_mask))), ...
    'One or more semantic trials were not decoded exactly once.');
assert(all(isnan(trial_predictions(~semantic_mask))) && ...
    all(isnan(trial_evidence(~semantic_mask))), ...
    'CONTROL trials must not receive predictions or evidence.');

odor_context_evidence = nan(n_odors, 3);
odor_context_accuracy = nan(n_odors, 3);
for odor_idx = 1:n_odors
    for context_idx = 1:3
        cell_mask = semantic_mask & odor_id == odor_idx & context_id == context_idx;
        if ~any(cell_mask)
            continue;
        end
        odor_context_evidence(odor_idx, context_idx) = mean(trial_evidence(cell_mask));
        odor_context_accuracy(odor_idx, context_idx) = ...
            mean(trial_predictions(cell_mask) == context_idx);
    end
end
by_odor_evidence = mean(odor_context_evidence, 2, 'omitmissing');
by_odor_accuracy = mean(odor_context_accuracy, 2, 'omitmissing');
by_context_evidence = mean(odor_context_evidence, 1, 'omitmissing');
by_context_recall = mean(odor_context_accuracy, 1, 'omitmissing');

true_context = context_id(semantic_mask);
predicted_context = trial_predictions(semantic_mask);
confusion = accumarray([true_context, predicted_context], 1, [3, 3]);
confusion_normalized = confusion ./ sum(confusion, 2);

decoded = struct();
decoded.semantic_contexts = cellstr(semantic_contexts(:));
decoded.odor_values = odor_values;
decoded.chance_accuracy = 1/3;
decoded.overall = struct( ...
    'evidence', mean(by_odor_evidence), ...
    'accuracy', mean(by_odor_accuracy));
decoded.by_odor = struct('evidence', by_odor_evidence, 'accuracy', by_odor_accuracy);
decoded.by_context = struct('evidence', by_context_evidence, 'recall', by_context_recall);
decoded.odor_context = struct( ...
    'evidence', odor_context_evidence, 'accuracy', odor_context_accuracy);
decoded.confusion_matrix = confusion;
decoded.confusion_matrix_normalized = confusion_normalized;
decoded.aggregation = ['odor-context cells averaged within target odor; ' ...
    'target odors then weighted equally'];

if opts.CollectDetails
    decoded.trial_predictions = trial_predictions;
    decoded.trial_evidence = trial_evidence;
    decoded.fold_accuracy = fold_accuracy;
    decoded.fold_evidence = fold_evidence;
    decoded.dimension_order = struct( ...
        'odor_context', {{'target_odor', 'true_context'}}, ...
        'fold', {{'target_odor', 'held_out_run'}});
else
    decoded.trial_predictions = [];
    decoded.trial_evidence = [];
    decoded.fold_accuracy = [];
    decoded.fold_evidence = [];
end
end

function templates = make_context_templates(Xtrain, ytrain)
templates = nan(3, size(Xtrain, 2));
for context_idx = 1:3
    templates(context_idx, :) = mean(Xtrain(ytrain == context_idx, :), 1);
end
assert(all(isfinite(templates), 'all'), 'A context template contains non-finite values.');
end

function [predictions, evidence] = score_context_templates(Xtest, templates, ytest)
Xtest = Xtest - mean(Xtest, 2);
templates = templates - mean(templates, 2);
test_norm = sqrt(sum(Xtest.^2, 2));
template_norm = sqrt(sum(templates.^2, 2));
assert(all(test_norm > eps & isfinite(test_norm)), ...
    'A test pattern has zero or invalid spatial variance.');
assert(all(template_norm > eps & isfinite(template_norm)), ...
    'A context template has zero or invalid spatial variance.');
similarity = (Xtest ./ test_norm) * (templates ./ template_norm)';
[~, predictions] = max(similarity, [], 2);

n_tests = size(similarity, 1);
correct_linear = sub2ind(size(similarity), (1:n_tests)', ytest(:));
correct_similarity = similarity(correct_linear);
incorrect_mean = (sum(similarity, 2) - correct_similarity) / 2;
evidence = correct_similarity - incorrect_mean;
end
