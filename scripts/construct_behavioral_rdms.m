%% Construct categorical and rating-based 80-by-80 RDMs
% Conditions are ordered by context, then odor:
% PERSON odors 1:20, FOOD odors 1:20, LOCATION odors 1:20, and
% CONTROL odors 1:20. The canonical analysis sample is subj_2--subj_6.

script_path = mfilename('fullpath');
assert(~isempty(script_path), 'Run this file as a script, not pasted code.');
project_root = fileparts(fileparts(script_path));
behavior_root = fullfile(project_root, 'behavior');
output_root = fullfile(project_root, 'RDMs');
if ~isfolder(output_root)
    mkdir(output_root);
end

subject_ids = 2:6;
context_order = ["PERSON", "FOOD", "LOCATION", "CONTROL"];
odor_ids = (1:20)';
n_contexts = numel(context_order);
n_odors = numel(odor_ids);
n_conditions = n_contexts * n_odors;
expected_repetitions = 10;

% Fixed condition order: all odors within context 1, then context 2, etc.
condition_context = repelem(context_order(:), n_odors);
condition_odor = repmat(odor_ids, n_contexts, 1);
condition_index = (1:n_conditions)';

% These two matrices intentionally follow the requested coding: 1 means
% "same" and 0 means "different" (so they are categorical similarity
% models even though the output variables are named RDMs).
D_odor = double(condition_odor == condition_odor');
D_context = double(condition_context == condition_context');

% Use subject 2's labels to attach human-readable odor names to the common
% condition metadata. Odor identity itself is always represented by ID.
label_data = load(fullfile(behavior_root, 'subj_2', 'behavior.mat'), 'odorlabels');
assert(isfield(label_data, 'odorlabels') && numel(label_data.odorlabels) == n_odors, ...
    'Expected 20 odor labels in subj_2/behavior.mat.');
odor_labels = string(label_data.odorlabels(:));
condition_odor_label = odor_labels(condition_odor);
condition_label = condition_context + "_odor" + string(condition_odor);
condition_metadata = table(condition_index, condition_context, condition_odor, ...
    condition_odor_label, condition_label);

save(fullfile(output_root, 'categorical_RDMs.mat'), ...
    'D_odor', 'D_context', 'condition_metadata', 'context_order', 'odor_ids');

for subject_id = subject_ids
    subject_name = sprintf('subj_%d', subject_id);
    input_file = fullfile(behavior_root, subject_name, 'behavior.mat');
    assert(isfile(input_file), 'Missing behavioral input: %s', input_file);
    data = load(input_file, 'valence_all', 'intensity_all', 'odor', 'category');

    required_fields = {'valence_all', 'intensity_all', 'odor', 'category'};
    assert(all(isfield(data, required_fields)), ...
        '%s is missing one or more required variables.', input_file);

    pleasantness = double(data.valence_all(:));
    intensity = double(data.intensity_all(:));
    trial_odor = double(data.odor(:));
    trial_context = upper(strtrim(string(data.category(:))));
    n_trials = numel(pleasantness);

    assert(numel(intensity) == n_trials && numel(trial_odor) == n_trials && ...
        numel(trial_context) == n_trials, 'Trial vectors differ in length for %s.', subject_name);
    assert(n_trials == n_conditions * expected_repetitions, ...
        'Expected %d trials for %s, found %d.', ...
        n_conditions * expected_repetitions, subject_name, n_trials);
    assert(all(ismember(trial_odor, odor_ids)), 'Unexpected odor ID for %s.', subject_name);
    assert(all(ismember(trial_context, context_order)), 'Unexpected context for %s.', subject_name);

    % Standardize each rating across all non-NaN trials within subject.
    pleasantness_z = standardize_omitnan(pleasantness, subject_name, 'pleasantness');
    intensity_z = standardize_omitnan(intensity, subject_name, 'intensity');

    mean_pleasantness_z = nan(n_conditions, 1);
    mean_intensity_z = nan(n_conditions, 1);
    total_trials_by_condition = zeros(n_conditions, 1);
    valid_pleasantness_trials = zeros(n_conditions, 1);
    valid_intensity_trials = zeros(n_conditions, 1);

    for condition_id = 1:n_conditions
        condition_mask = trial_context == condition_context(condition_id) & ...
            trial_odor == condition_odor(condition_id);
        total_trials_by_condition(condition_id) = sum(condition_mask);
        assert(total_trials_by_condition(condition_id) > 0, ...
            '%s condition %s has no trials.', subject_name, condition_label(condition_id));

        pleasantness_values = pleasantness_z(condition_mask);
        intensity_values = intensity_z(condition_mask);
        valid_pleasantness_trials(condition_id) = sum(~isnan(pleasantness_values));
        valid_intensity_trials(condition_id) = sum(~isnan(intensity_values));
        mean_pleasantness_z(condition_id) = mean(pleasantness_values, 'omitnan');
        mean_intensity_z(condition_id) = mean(intensity_values, 'omitnan');
    end

    if any(total_trials_by_condition ~= expected_repetitions)
        warning('%s has %d conditions with counts other than %d (range %d--%d); averaging all available trials.', ...
            subject_name, sum(total_trials_by_condition ~= expected_repetitions), ...
            expected_repetitions, min(total_trials_by_condition), max(total_trials_by_condition));
    end

    assert(all(isfinite(mean_pleasantness_z)), ...
        'At least one %s condition has no valid pleasantness ratings.', subject_name);
    assert(all(isfinite(mean_intensity_z)), ...
        'At least one %s condition has no valid intensity ratings.', subject_name);

    % Absolute differences are required for symmetric RDMs. Signed
    % (i-minus-j) differences would be antisymmetric instead.
    D_pleasantness = abs(mean_pleasantness_z - mean_pleasantness_z');
    D_intensity = abs(mean_intensity_z - mean_intensity_z');
    D_pleasantness(1:n_conditions+1:end) = 0;
    D_intensity(1:n_conditions+1:end) = 0;

    assert(isequal(size(D_pleasantness), [80, 80]) && ...
        isequal(size(D_intensity), [80, 80]), 'RDM size check failed for %s.', subject_name);
    assert(isequaln(D_pleasantness, D_pleasantness') && ...
        isequaln(D_intensity, D_intensity'), 'RDM symmetry check failed for %s.', subject_name);

    output_file = fullfile(output_root, sprintf('%s_behavioral_RDMs.mat', subject_name));
    save(output_file, 'D_pleasantness', 'D_intensity', ...
        'mean_pleasantness_z', 'mean_intensity_z', ...
        'total_trials_by_condition', ...
        'valid_pleasantness_trials', 'valid_intensity_trials', ...
        'condition_metadata', 'context_order', 'odor_ids', 'subject_id');
    fprintf('Saved %s\n', output_file);
end

fprintf('Saved shared categorical RDMs and %d subject files in %s\n', ...
    numel(subject_ids), output_root);

function standardized = standardize_omitnan(values, subject_name, measure_name)
%STANDARDIZE_OMITNAN Z-score finite observations and retain missing values.
valid = ~isnan(values);
assert(sum(valid) >= 2, 'Too few valid %s ratings for %s.', measure_name, subject_name);
mu = mean(values(valid));
sigma = std(values(valid), 0);
assert(isfinite(sigma) && sigma > 0, ...
    '%s ratings have zero or invalid variance for %s.', measure_name, subject_name);
standardized = nan(size(values));
standardized(valid) = (values(valid) - mu) ./ sigma;
end
