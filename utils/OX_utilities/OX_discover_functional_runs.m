function runs = OX_discover_functional_runs(subjidx, varargin)
%OX_DISCOVER_FUNCTIONAL_RUNS Find GLMsingle-input EPIs and motion files.
%   RUNS = OX_DISCOVER_FUNCTIONAL_RUNS(SUBJIDX) discovers the smoothed,
%   realigned 4-D functional images (sr*.nii) used by the current
%   GLMsingle input code and pairs each with its SPM rp_*.txt file. Runs
%   are ordered numerically by session and then within-session run ID.
%
%   This function rejects duplicate session/run IDs, missing pairs,
%   non-contiguous run IDs within a session, and unexpected run counts.

p = inputParser;
p.addRequired('subjidx', @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x == fix(x));
p.addParameter('ProjectRoot', localProjectRoot(), @(x) ischar(x) || isstring(x));
p.addParameter('ExpectedRuns', 80, @(x) isnumeric(x) && isscalar(x) && x > 0 && x == fix(x));
p.parse(subjidx, varargin{:});

projectRoot = char(p.Results.ProjectRoot);
subject = sprintf('subj_%d', subjidx);
funcDir = fullfile(projectRoot, 'MRI', subject, 'nifti', 'func');
if ~isfolder(funcDir)
    error('OX:MotionQC:MissingFunctionalDirectory', ...
        'Functional directory does not exist: %s', funcDir);
end

funcListing = dir(fullfile(funcDir, 'sr*.nii'));
motionListing = dir(fullfile(funcDir, 'rp_*.txt'));
if isempty(funcListing)
    error('OX:MotionQC:MissingFunctionalFiles', ...
        'No GLMsingle-input sr*.nii files found in %s.', funcDir);
end
if isempty(motionListing)
    error('OX:MotionQC:MissingMotionFiles', ...
        'No rp_*.txt files found in %s.', funcDir);
end

funcEntries = parseEntries(funcListing, 'functional', '^sr(.+)_([0-9]+)_Run([0-9]+)_.+\.nii$');
motionEntries = parseEntries(motionListing, 'motion', '^rp_(.+)_([0-9]+)_Run([0-9]+)_.+\.txt$');
validateUniqueKeys(funcEntries, 'functional');
validateUniqueKeys(motionEntries, 'motion');
if ~strcmpi(funcEntries(1).subjectPrefix, motionEntries(1).subjectPrefix)
    error('OX:MotionQC:SubjectPrefixMismatch', ...
        'Functional and motion subject prefixes do not match (%s vs %s).', ...
        funcEntries(1).subjectPrefix, motionEntries(1).subjectPrefix);
end

funcKeys = [funcEntries.key];
motionKeys = [motionEntries.key];
missingMotion = setdiff(funcKeys, motionKeys);
missingFunc = setdiff(motionKeys, funcKeys);
if ~isempty(missingMotion)
    error('OX:MotionQC:MissingMotionPair', ...
        'Missing motion file(s) for session/run key(s): %s', numberList(missingMotion));
end
if ~isempty(missingFunc)
    error('OX:MotionQC:MissingFunctionalPair', ...
        'Missing functional file(s) for session/run key(s): %s', numberList(missingFunc));
end

[~, order] = sort(funcKeys);
funcEntries = funcEntries(order);
sessions = [funcEntries.session];
runNumbers = [funcEntries.run];
uniqueSessions = unique(sessions, 'stable');
if ~isequal(uniqueSessions, min(uniqueSessions):max(uniqueSessions))
    error('OX:MotionQC:AmbiguousSessionOrder', ...
        'Session IDs are not contiguous: %s', numberList(uniqueSessions));
end
for session = uniqueSessions
    observed = runNumbers(sessions == session);
    expected = 1:max(observed);
    if ~isequal(observed, expected)
        error('OX:MotionQC:AmbiguousRunOrder', ...
            'Session %d run IDs are not exactly 1:%d (observed: %s).', ...
            session, max(observed), numberList(observed));
    end
end
if numel(funcEntries) ~= p.Results.ExpectedRuns
    error('OX:MotionQC:UnexpectedRunCount', ...
        'Expected %d runs for %s, found %d.', ...
        p.Results.ExpectedRuns, subject, numel(funcEntries));
end

emptyRun = struct('ordinal', [], 'id', '', 'session', [], 'run', [], ...
    'functional_file', '', 'motion_file', '', 'source_stem_match', false);
runs = repmat(emptyRun, numel(funcEntries), 1);
for runIndex = 1:numel(funcEntries)
    funcEntry = funcEntries(runIndex);
    motionIndex = find(motionKeys == funcEntry.key, 1, 'first');
    motionEntry = motionEntries(motionIndex);
    runs(runIndex).ordinal = runIndex;
    runs(runIndex).id = sprintf('session%02d_run%02d', funcEntry.session, funcEntry.run);
    runs(runIndex).session = funcEntry.session;
    runs(runIndex).run = funcEntry.run;
    runs(runIndex).functional_file = fullfile(funcEntry.folder, funcEntry.name);
    runs(runIndex).motion_file = fullfile(motionEntry.folder, motionEntry.name);
    funcStem = funcEntry.name(3:end-4);
    motionStem = motionEntry.name(4:end-4);
    runs(runIndex).source_stem_match = strcmp(funcStem, motionStem);
end
end

function entries = parseEntries(listing, kind, expression)
emptyEntry = struct('name', '', 'folder', '', 'subjectPrefix', '', ...
    'session', [], 'run', [], 'key', []);
entries = repmat(emptyEntry, numel(listing), 1);
for index = 1:numel(listing)
    tokens = regexp(listing(index).name, expression, 'tokens', 'once');
    if isempty(tokens)
        error('OX:MotionQC:UnparseableFilename', ...
            'Cannot unambiguously parse %s filename: %s', kind, listing(index).name);
    end
    entries(index).name = listing(index).name;
    entries(index).folder = listing(index).folder;
    entries(index).subjectPrefix = tokens{1};
    entries(index).session = str2double(tokens{2});
    entries(index).run = str2double(tokens{3});
    entries(index).key = entries(index).session * 1000 + entries(index).run;
end
prefixes = unique(lower(string({entries.subjectPrefix})));
if numel(prefixes) ~= 1
    error('OX:MotionQC:AmbiguousSubjectPrefix', ...
        'Found multiple subject prefixes among %s files: %s', ...
        kind, strjoin(prefixes, ', '));
end
end

function validateUniqueKeys(entries, kind)
keys = [entries.key];
if numel(unique(keys)) ~= numel(keys)
    [uniqueKeys, ~, group] = unique(keys);
    counts = accumarray(group(:), 1);
    duplicateKeys = uniqueKeys(counts > 1);
    error('OX:MotionQC:DuplicateRunID', ...
        'Duplicate session/run IDs among %s files: %s', kind, numberList(duplicateKeys));
end
end

function root = localProjectRoot()
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

function value = numberList(numbers)
value = strjoin(arrayfun(@num2str, numbers, 'UniformOutput', false), ', ');
end
