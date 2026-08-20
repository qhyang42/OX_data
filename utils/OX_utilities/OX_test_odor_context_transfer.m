function test_results = OX_test_odor_context_transfer()
%OX_TEST_ODOR_CONTEXT_TRANSFER Synthetic verification for template transfer.

rng(11, 'twister');
n_odors = 20;
n_contexts = 4;
n_runs = 3;
n_voxels = 120;
context_order = ["PERSON", "FOOD", "LOCATION", "CONTROL"];

odors = repmat((1:n_odors)', n_contexts * n_runs, 1);
contexts_one_run = repelem(context_order(:), n_odors);
contexts = repmat(contexts_one_run, n_runs, 1);
run_ids = repelem((1:n_runs)', n_odors * n_contexts);
n_trials = numel(odors);

% Context-invariant identity should transfer equally across semantic cells.
odor_patterns = randn(n_odors, n_voxels);
invariant_X = zeros(n_trials, n_voxels);
for trial_idx = 1:n_trials
    invariant_X(trial_idx, :) = odor_patterns(odors(trial_idx), :) + ...
        0.15 * randn(1, n_voxels);
end
invariant = OX_template_context_transfer(invariant_X, odors, contexts, run_ids);
assert(invariant.primary.within_context_accuracy > 0.95, ...
    'Context-invariant within-context decoding is unexpectedly weak.');
assert(invariant.primary.cross_context_accuracy > 0.95, ...
    'Context-invariant cross-context decoding is unexpectedly weak.');
assert(abs(invariant.primary.evidence_effect) < 0.02, ...
    'Context-invariant patterns produced a spurious specificity effect.');
assert(invariant.control_transfer.mean_accuracy > 0.95, ...
    'Shared CONTROL patterns failed to transfer into semantic contexts.');

% Independent semantic identity patterns should favor matrix diagonals.
semantic_patterns = randn(3, n_odors, n_voxels);
specific_X = zeros(n_trials, n_voxels);
for trial_idx = 1:n_trials
    context_idx = find(context_order == contexts(trial_idx));
    if context_idx <= 3
        signal = reshape(semantic_patterns(context_idx, odors(trial_idx), :), 1, []);
    else
        signal = odor_patterns(odors(trial_idx), :);
    end
    specific_X(trial_idx, :) = signal + 0.15 * randn(1, n_voxels);
end
specific = OX_template_context_transfer(specific_X, odors, contexts, run_ids);
assert(specific.primary.within_context_accuracy > 0.95, ...
    'Context-specific within-context decoding is unexpectedly weak.');
assert(specific.primary.cross_context_accuracy < 0.20, ...
    'Independent context-specific patterns transferred unexpectedly well.');
assert(specific.primary.evidence_effect > 0.70, ...
    'Context-specific patterns did not produce the expected positive effect.');

% Labels shuffled separately within context and run should approach null.
shuffled_odors = odors;
for run_idx = 1:n_runs
    for context_idx = 1:n_contexts
        indices = find(run_ids == run_idx & contexts == context_order(context_idx));
        shuffled_odors(indices) = shuffled_odors(indices(randperm(numel(indices))));
    end
end
shuffled = OX_template_context_transfer(invariant_X, shuffled_odors, contexts, run_ids);
assert(abs(mean(shuffled.semantic.evidence, 'all')) < 0.10, ...
    'Shuffled-label evidence is too far from zero.');
assert(mean(shuffled.semantic.accuracy, 'all') < 0.15, ...
    'Shuffled-label accuracy is too far above 5%% chance.');

% Every eligible trial must receive exactly one held-out-run prediction.
for source_idx = 1:3
    for target_idx = 1:3
        predictions = squeeze(invariant.trial_predictions(source_idx, target_idx, :));
        expected = contexts == context_order(target_idx);
        assert(all(isfinite(predictions(expected))) && all(isnan(predictions(~expected))), ...
            'Trial prediction coverage is incorrect for semantic cell %d -> %d.', ...
            source_idx, target_idx);
    end
end

test_results = struct( ...
    'passed', true, ...
    'invariant_evidence_effect', invariant.primary.evidence_effect, ...
    'specific_evidence_effect', specific.primary.evidence_effect, ...
    'shuffled_mean_evidence', mean(shuffled.semantic.evidence, 'all'), ...
    'shuffled_mean_accuracy', mean(shuffled.semantic.accuracy, 'all'));
fprintf(['OX odor-context transfer tests passed | invariant effect %.4f | ' ...
    'specific effect %.4f | shuffled accuracy %.4f\n'], ...
    test_results.invariant_evidence_effect, test_results.specific_evidence_effect, ...
    test_results.shuffled_mean_accuracy);
end
