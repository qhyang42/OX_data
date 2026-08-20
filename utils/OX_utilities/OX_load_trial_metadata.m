function trial_metadata = OX_load_trial_metadata(subjidx, varargin)
%OX_LOAD_TRIAL_METADATA Load trial labels and exact run IDs in GLMsingle order.
%
%   trial_metadata = OX_load_trial_metadata(subjidx, Name, Value, ...)
%
% Cue files are traversed by numeric session and run number, matching the
% ordering used by OX_get_odor and the single-trial design construction.
% The returned table has one row per trial and the variables trial_index,
% session_id, run_in_session, run_id, trial_in_run, odor, context, and
% cue_file.

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'subjidx', @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x == round(x));
addParameter(p, 'CueRoot', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'ExpectedTrials', 800, @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == round(x));
addParameter(p, 'ExpectedRuns', 80, @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == round(x));
addParameter(p, 'ExpectedOdors', 20, @(x) isnumeric(x) && isscalar(x) && x >= 2 && x == round(x));
parse(p, subjidx, varargin{:});
opts = p.Results;

assert(subjidx >= 1 && subjidx <= 6, 'subjidx must be between 1 and 6.');
subject_name = sprintf('subj_%d', subjidx);
if strlength(string(opts.CueRoot)) == 0
    utility_dir = fileparts(mfilename('fullpath'));
    project_root = fileparts(fileparts(utility_dir));
    cue_root = fullfile(project_root, 'cuelist');
else
    cue_root = char(string(opts.CueRoot));
end
subject_cue_dir = fullfile(cue_root, subject_name);
assert(isfolder(subject_cue_dir), 'Missing cue directory: %s', subject_cue_dir);

trial_index = zeros(0, 1);
session_id = zeros(0, 1);
run_in_session = zeros(0, 1);
run_id = zeros(0, 1);
trial_in_run = zeros(0, 1);
odor = zeros(0, 1);
context = strings(0, 1);
cue_file = strings(0, 1);

global_run_id = 0;
global_trial_id = 0;
for session = 1:16
    for run = 1:10
        filename = fullfile(subject_cue_dir, sprintf('session%d', session), ...
            sprintf('cuelist_sess%d_run%d.mat', session, run));
        if ~isfile(filename)
            continue;
        end

        loaded = load(filename, 'cuelist');
        assert(isfield(loaded, 'cuelist') && isstruct(loaded.cuelist), ...
            'Cue file lacks a cuelist structure: %s', filename);
        assert(isfield(loaded.cuelist, 'odor') && isfield(loaded.cuelist, 'category'), ...
            'Cue file lacks odor or category: %s', filename);

        run_odors = loaded.cuelist.odor(:);
        run_contexts = upper(strtrim(string(loaded.cuelist.category(:))));
        assert(isnumeric(run_odors) && all(isfinite(run_odors)), ...
            'Odor labels must be finite numeric values: %s', filename);
        assert(numel(run_odors) == numel(run_contexts), ...
            'Odor/context length mismatch in %s.', filename);
        assert(~isempty(run_odors), 'Cue file contains no trials: %s', filename);

        n_run_trials = numel(run_odors);
        global_run_id = global_run_id + 1;
        new_trials = global_trial_id + (1:n_run_trials)';
        trial_index = [trial_index; new_trials]; %#ok<AGROW>
        session_id = [session_id; repmat(session, n_run_trials, 1)]; %#ok<AGROW>
        run_in_session = [run_in_session; repmat(run, n_run_trials, 1)]; %#ok<AGROW>
        run_id = [run_id; repmat(global_run_id, n_run_trials, 1)]; %#ok<AGROW>
        trial_in_run = [trial_in_run; (1:n_run_trials)']; %#ok<AGROW>
        odor = [odor; double(run_odors)]; %#ok<AGROW>
        context = [context; run_contexts]; %#ok<AGROW>
        cue_file = [cue_file; repmat(string(filename), n_run_trials, 1)]; %#ok<AGROW>
        global_trial_id = global_trial_id + n_run_trials;
    end
end

trial_metadata = table(trial_index, session_id, run_in_session, run_id, ...
    trial_in_run, odor, context, cue_file);

context_order = ["PERSON", "FOOD", "LOCATION", "CONTROL"];
assert(height(trial_metadata) == opts.ExpectedTrials, ...
    'Loaded %d trials for %s; expected %d.', ...
    height(trial_metadata), subject_name, opts.ExpectedTrials);
assert(numel(unique(run_id)) == opts.ExpectedRuns, ...
    'Loaded %d runs for %s; expected %d.', ...
    numel(unique(run_id)), subject_name, opts.ExpectedRuns);
assert(numel(unique(odor)) == opts.ExpectedOdors, ...
    'Loaded %d odor identities for %s; expected %d.', ...
    numel(unique(odor)), subject_name, opts.ExpectedOdors);
assert(isequal(sort(unique(context)), sort(context_order(:))), ...
    'Contexts for %s must be PERSON, FOOD, LOCATION, and CONTROL.', subject_name);

odor_values = unique(odor, 'sorted');
for context_idx = 1:numel(context_order)
    context_mask = context == context_order(context_idx);
    assert(isequal(unique(odor(context_mask), 'sorted'), odor_values), ...
        'Context %s does not contain all odors for %s.', context_order(context_idx), subject_name);
end

% Every leave-one-run-out source-context fold must retain all odor classes.
runs = unique(run_id, 'stable');
for context_idx = 1:numel(context_order)
    for run_idx = 1:numel(runs)
        train_mask = run_id ~= runs(run_idx) & context == context_order(context_idx);
        assert(isequal(unique(odor(train_mask), 'sorted'), odor_values), ...
            'Training fold run %d lacks an odor in context %s for %s.', ...
            runs(run_idx), context_order(context_idx), subject_name);
    end
end
end
