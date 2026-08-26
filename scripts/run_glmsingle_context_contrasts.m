function manifest = run_glmsingle_context_contrasts(varargin)
%RUN_GLMSINGLE_CONTEXT_CONTRASTS Build all planned GLMsingle context maps.
%   MANIFEST = RUN_GLMSINGLE_CONTEXT_CONTRASTS() writes three context
%   contrast images for subjects 2:6 and for sniff and countdown GLMsingle
%   estimates. Each image is saved directly in its source estimates folder.
%
%   Optional name-value inputs:
%     SubjectIds  - integer subject IDs (default: 2:6)
%     EventTypes  - sniff and/or countdown (default: both)
%     ProjectRoot - OX_DATA project root (default: inferred from this file)
%     Overwrite   - replace existing outputs (default: false)


driverFile = mfilename('fullpath');
defaultProjectRoot = fileparts(fileparts(driverFile));
addpath(fullfile(defaultProjectRoot, 'environment'));
setup_ox();

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'SubjectIds', 2:6, @(x) isnumeric(x) && isvector(x) && ...
    all(isfinite(x)) && all(x == round(x)) && all(x >= 1));
addParameter(p, 'EventTypes', ["sniff", "countdown"], ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(p, 'ProjectRoot', defaultProjectRoot, ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'Overwrite', false, ...
    @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
parse(p, varargin{:});
opts = p.Results;

subjectIds = double(opts.SubjectIds(:)');
eventTypes = lower(strtrim(string(opts.EventTypes(:)')));
assert(~isempty(eventTypes) && all(ismember(eventTypes, ["sniff", "countdown"])), ...
    'OX:GLMsingleContrasts:InvalidEventTypes', ...
    'EventTypes may contain only sniff and countdown.');
assert(numel(unique(eventTypes)) == numel(eventTypes), ...
    'OX:GLMsingleContrasts:DuplicateEventTypes', ...
    'EventTypes must not contain duplicates.');
projectRoot = char(string(opts.ProjectRoot));

nAnalyses = numel(subjectIds) * numel(eventTypes);
manifests = cell(nAnalyses, 1);
analysisIndex = 0;
for subjectId = subjectIds
    for eventType = eventTypes
        analysisIndex = analysisIndex + 1;
        manifests{analysisIndex} = OX_glmsingle_context_contrasts( ...
            subjectId, eventType, 'ProjectRoot', projectRoot, ...
            'Overwrite', logical(opts.Overwrite));
    end
end

manifest = vertcat(manifests{:});
fprintf('Completed %d subject/event analyses and wrote %d contrast maps.\n', ...
    nAnalyses, height(manifest));
end
