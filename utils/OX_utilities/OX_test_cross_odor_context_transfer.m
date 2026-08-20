function test_results = OX_test_cross_odor_context_transfer()
%OX_TEST_CROSS_ODOR_CONTEXT_TRANSFER Synthetic decoder verification.

rng(19, 'twister');
n_odors = 20;
n_contexts = 4;
n_runs = 4;
n_voxels = 180;
context_order = ["PERSON", "FOOD", "LOCATION", "CONTROL"];
odors = repmat((1:n_odors)', n_contexts * n_runs, 1);
contexts = repmat(repelem(context_order(:), n_odors), n_runs, 1);
run_ids = repelem((1:n_runs)', n_odors * n_contexts);
n_trials = numel(odors);

% Shared semantic context code must generalize despite strong odor nuisance.
odor_patterns = 3 * randn(n_odors, n_voxels);
shared_context_patterns = 2 * randn(3, n_voxels);
shared_X = zeros(n_trials, n_voxels);
for trial_idx = 1:n_trials
    context_idx = find(context_order == contexts(trial_idx));
    signal = odor_patterns(odors(trial_idx), :);
    if context_idx <= 3
        signal = signal + shared_context_patterns(context_idx, :);
    end
    shared_X(trial_idx, :) = signal + 0.20 * randn(1, n_voxels);
end
shared = OX_template_cross_odor_context_transfer( ...
    shared_X, odors, contexts, run_ids);
assert(shared.overall.accuracy > 0.90 && shared.overall.evidence > 0.25, ...
    'Shared context patterns did not generalize across held-out odors.');

% Independent context patterns for every odor must not generalize.
odor_specific_contexts = randn(n_odors, 3, n_voxels);
specific_X = zeros(n_trials, n_voxels);
for trial_idx = 1:n_trials
    context_idx = find(context_order == contexts(trial_idx));
    signal = odor_patterns(odors(trial_idx), :);
    if context_idx <= 3
        signal = signal + 2 * reshape( ...
            odor_specific_contexts(odors(trial_idx), context_idx, :), 1, []);
    end
    specific_X(trial_idx, :) = signal + 0.20 * randn(1, n_voxels);
end
specific = OX_template_cross_odor_context_transfer( ...
    specific_X, odors, contexts, run_ids);
assert(abs(specific.overall.evidence) < 0.10 && specific.overall.accuracy < 0.50, ...
    'Odor-specific context patterns transferred unexpectedly well.');

% With no context component, decoding must remain near chance.
null_X = zeros(n_trials, n_voxels);
for trial_idx = 1:n_trials
    null_X(trial_idx, :) = odor_patterns(odors(trial_idx), :) + ...
        0.20 * randn(1, n_voxels);
end
null_result = OX_template_cross_odor_context_transfer( ...
    null_X, odors, contexts, run_ids);
assert(abs(null_result.overall.evidence) < 0.05 && ...
    abs(null_result.overall.accuracy - 1/3) < 0.10, ...
    'No-context synthetic data did not approach the expected null.');

% Within-run semantic-label shuffling must reproduce the empirical null.
n_shuffles = 20;
shuffled_evidence = nan(n_shuffles, 1);
shuffled_accuracy = nan(n_shuffles, 1);
for shuffle_idx = 1:n_shuffles
    shuffled_contexts = contexts;
    for run_idx = 1:n_runs
        indices = find(run_ids == run_idx & contexts ~= "CONTROL");
        shuffled_contexts(indices) = shuffled_contexts(indices(randperm(numel(indices))));
    end
    shuffled = OX_template_cross_odor_context_transfer( ...
        shared_X, odors, shuffled_contexts, run_ids);
    shuffled_evidence(shuffle_idx) = shuffled.overall.evidence;
    shuffled_accuracy(shuffle_idx) = shuffled.overall.accuracy;
end
mean_shuffled_evidence = mean(shuffled_evidence);
mean_shuffled_accuracy = mean(shuffled_accuracy);
assert(abs(mean_shuffled_evidence) < 0.05 && ...
    abs(mean_shuffled_accuracy - 1/3) < 0.05, ...
    'Within-run shuffled labels retained unexpected context information.');

assert(all(isnan(shared.trial_predictions(contexts == "CONTROL"))), ...
    'CONTROL trials received predictions.');
assert(isequal(string(shared.semantic_contexts), ...
    ["PERSON"; "FOOD"; "LOCATION"]), 'Semantic context ordering changed.');

test_results = struct( ...
    'passed', true, ...
    'shared_evidence', shared.overall.evidence, ...
    'shared_accuracy', shared.overall.accuracy, ...
    'odor_specific_evidence', specific.overall.evidence, ...
    'odor_specific_accuracy', specific.overall.accuracy, ...
    'null_evidence', null_result.overall.evidence, ...
    'null_accuracy', null_result.overall.accuracy, ...
    'shuffled_evidence', mean_shuffled_evidence, ...
    'shuffled_accuracy', mean_shuffled_accuracy);
fprintf(['Cross-odor context tests passed | shared acc %.3f | ' ...
    'odor-specific acc %.3f | null acc %.3f | shuffled acc %.3f\n'], ...
    test_results.shared_accuracy, test_results.odor_specific_accuracy, ...
    test_results.null_accuracy, test_results.shuffled_accuracy);
end
