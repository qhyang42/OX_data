function permuted_contexts = OX_permute_semantic_context_within_runs( ...
        context_labels, run_ids, random_seed)
%OX_PERMUTE_SEMANTIC_CONTEXT_WITHIN_RUNS Shuffle semantic labels by run.
%
% The PERSON/FOOD/LOCATION label vector is independently permuted within
% every exact acquisition run. Thus each run retains its observed context
% counts while the association between trials and semantic labels is broken.

assert(isnumeric(run_ids) && isvector(run_ids), ...
    'run_ids must be a numeric vector.');
assert(isnumeric(random_seed) && isscalar(random_seed) && ...
    isfinite(random_seed) && random_seed >= 0, ...
    'random_seed must be a nonnegative finite scalar.');
contexts = upper(strtrim(string(context_labels(:))));
run_ids = double(run_ids(:));
assert(numel(contexts) == numel(run_ids), ...
    'Context labels and run IDs must have equal lengths.');
semantic_contexts = ["PERSON", "FOOD", "LOCATION"];
assert(all(ismember(contexts, semantic_contexts)), ...
    'Only PERSON, FOOD, and LOCATION labels may be permuted.');
assert(all(isfinite(run_ids)), 'run_ids must be finite.');

stream = RandStream('mt19937ar', 'Seed', mod(floor(random_seed), 2^32));
permuted_contexts = contexts;
runs = unique(run_ids, 'stable');
for run_idx = 1:numel(runs)
    indices = find(run_ids == runs(run_idx));
    if numel(indices) > 1
        order = randperm(stream, numel(indices));
        permuted_contexts(indices) = contexts(indices(order));
    end
end

% Assert the intended exchangeability block was preserved exactly.
for run_idx = 1:numel(runs)
    run_mask = run_ids == runs(run_idx);
    assert(isequal(sort(permuted_contexts(run_mask)), ...
        sort(contexts(run_mask))), ...
        'Context counts changed in run %g.', runs(run_idx));
end
end
