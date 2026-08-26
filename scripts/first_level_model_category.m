%% Run the sniff-category first-level SPM model for OX subjects 2-6.
% Sniff events use BreathMetrics inhale onsets. Each functional run is an
% independent SPM session with its own motion, respiratory, and bad-volume
% nuisance regressors. The semantic contrasts compare Person, Food, and
% Location with the equally weighted mean of the other three categories.

info = setup_ox;

subjects = 2:6;
TR = 0.76;
highPassSeconds = 128;
fweP = 0.001;

spm('Defaults', 'fMRI');
spm_jobman('initcfg');

for subjidx = subjects
    subject = sprintf('subj_%d', subjidx);
    niftiDir = fullfile(info.projectRoot, 'MRI', subject, 'nifti');
    eventFile = fullfile(info.projectRoot, 'labchart', 'extracted_events', ...
        sprintf('subj%d_events_bm.mat', subjidx));
    outputDir = fullfile(niftiDir, ...
        'first_level_model_category_sniff_physio');

    events = loadBreathMetricsEvents(eventFile, subject);
    runs = OX_discover_functional_runs(subjidx, ...
        'ProjectRoot', info.projectRoot, ...
        'ExpectedRuns', numel(events.nframes));
    categories = loadRunCategories(subject, size(events.event_onsets));

    matlabbatch = buildCategoryBatch(subject, runs, events, categories, ...
        niftiDir, outputDir, TR, highPassSeconds, fweP);

    if ~isfolder(outputDir)
        mkdir(outputDir);
    end

    fprintf('Running sniff-category first-level model for %s (%d runs).\n', ...
        subject, numel(runs));
    spm_jobman('run', matlabbatch);
end

%%
function events = loadBreathMetricsEvents(eventFile, subject)
requiredFields = {'event_onsets', 'nframes', 'bm_processing'};
if ~isfile(eventFile)
    error('OX:FirstLevelCategory:MissingBreathMetricsEvents', ...
        'BreathMetrics event file not found for %s: %s', subject, eventFile);
end

events = load(eventFile, requiredFields{:});
missingFields = requiredFields(~isfield(events, requiredFields));
if ~isempty(missingFields)
    error('OX:FirstLevelCategory:InvalidBreathMetricsEvents', ...
        '%s is missing required field(s): %s', eventFile, ...
        strjoin(missingFields, ', '));
end
if ~isnumeric(events.event_onsets) || ~isnumeric(events.nframes) || ...
        isempty(events.event_onsets) || isempty(events.nframes) || ...
        any(~isfinite(events.event_onsets(:))) || ...
        any(~isfinite(events.nframes(:)))
    error('OX:FirstLevelCategory:InvalidBreathMetricsEvents', ...
        '%s contains invalid event_onsets or nframes.', eventFile);
end

events.nframes = events.nframes(:);
if size(events.event_onsets, 2) ~= numel(events.nframes)
    error('OX:FirstLevelCategory:EventRunCountMismatch', ...
        ['%s contains %d event-onset run(s), but nframes contains ', ...
         '%d run(s).'], eventFile, size(events.event_onsets, 2), ...
        numel(events.nframes));
end
if any(events.nframes <= 0) || any(events.nframes ~= fix(events.nframes))
    error('OX:FirstLevelCategory:InvalidFrameCounts', ...
        '%s contains invalid run frame counts.', eventFile);
end
end

function categories = loadRunCategories(subject, expectedSize)
[~, categoryVector] = OX_get_odor(subject);
if numel(categoryVector) ~= prod(expectedSize)
    error('OX:FirstLevelCategory:CategoryCountMismatch', ...
        '%s has %d category labels; expected %d.', subject, ...
        numel(categoryVector), prod(expectedSize));
end

categories = reshape(upper(string(categoryVector)), expectedSize);
categoryNames = ["PERSON", "FOOD", "LOCATION", "CONTROL"];
unexpected = setdiff(unique(categories), categoryNames);
if ~isempty(unexpected)
    error('OX:FirstLevelCategory:UnexpectedCategory', ...
        '%s has unexpected category label(s): %s', subject, ...
        strjoin(unexpected, ', '));
end
if any(ismissing(categories), 'all')
    error('OX:FirstLevelCategory:MissingCategory', ...
        '%s has missing category labels.', subject);
end
end

function matlabbatch = buildCategoryBatch(subject, runs, events, ...
        categories, niftiDir, outputDir, TR, highPassSeconds, fweP)
categoryLabels = ["PERSON", "FOOD", "LOCATION", "CONTROL"];
conditionNames = {'Person', 'Food', 'Location', 'Control'};
semanticContrastNames = {'Person > Other Contexts', ...
    'Food > Other Contexts', 'Location > Other Contexts'};
contrastTemplate = eye(4) - (ones(4) - eye(4)) / 3;

numberOfRuns = numel(runs);
if numberOfRuns ~= numel(events.nframes) || ...
        size(categories, 2) ~= numberOfRuns
    error('OX:FirstLevelCategory:RunCountMismatch', ...
        '%s run, event, and category counts do not agree.', subject);
end

matlabbatch = cell(1, 4);
matlabbatch{1}.spm.stats.fmri_spec.dir = {outputDir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 30;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;

categoryColumns = cell(1, numel(categoryLabels));
designColumn = 0;
for runIndex = 1:numberOfRuns
    runInfo = runs(runIndex);
    scanList = cellstr(spm_select('expand', runInfo.functional_file));
    expectedFrames = events.nframes(runIndex);
    if numel(scanList) ~= expectedFrames
        error('OX:FirstLevelCategory:FunctionalFrameMismatch', ...
            '%s %s has %d functional frames; expected %d.', subject, ...
            runInfo.id, numel(scanList), expectedFrames);
    end

    onsets = events.event_onsets(:, runIndex);
    runDuration = expectedFrames * TR;
    if any(onsets < 0) || any(onsets >= runDuration)
        error('OX:FirstLevelCategory:OnsetOutsideRun', ...
            '%s %s has a sniff onset outside [0, %.3f) seconds.', ...
            subject, runInfo.id, runDuration);
    end

    confoundFile = fullfile(niftiDir, 'glmsingle_confounds', sprintf( ...
        '%s_session%02d_run%02d_confounds.mat', subject, ...
        runInfo.session, runInfo.run));
    [confounds, confoundNames] = loadRunConfounds( ...
        confoundFile, subject, runInfo, runIndex, expectedFrames);

    runCategories = categories(:, runIndex);
    conditions = struct('name', {}, 'onset', {}, 'duration', {}, ...
        'tmod', {}, 'pmod', {}, 'orth', {});
    for categoryIndex = 1:numel(categoryLabels)
        inCategory = runCategories == categoryLabels(categoryIndex);
        if ~any(inCategory)
            continue;
        end

        conditionIndex = numel(conditions) + 1;
        conditions(conditionIndex).name = conditionNames{categoryIndex};
        conditions(conditionIndex).onset = onsets(inCategory);
        conditions(conditionIndex).duration = 0;
        conditions(conditionIndex).tmod = 0;
        conditions(conditionIndex).pmod = ...
            struct('name', {}, 'param', {}, 'poly', {});
        conditions(conditionIndex).orth = 1;

        categoryColumns{categoryIndex}(end + 1) = ...
            designColumn + conditionIndex;
    end

    session = struct();
    session.scans = scanList;
    session.cond = conditions;
    session.multi = {''};
    session.regress = struct('name', confoundNames, ...
        'val', num2cell(confounds, 1));
    session.multi_reg = {''};
    session.hpf = highPassSeconds;
    matlabbatch{1}.spm.stats.fmri_spec.sess(runIndex) = session;

    designColumn = designColumn + numel(conditions) + size(confounds, 2);
end

matlabbatch{1}.spm.stats.fmri_spec.fact = ...
    struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

spmMat = fullfile(outputDir, 'SPM.mat');
matlabbatch{2}.spm.stats.fmri_est.spmmat = {spmMat};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

numberOfDesignColumns = designColumn + numberOfRuns;
matlabbatch{3}.spm.stats.con.spmmat = {spmMat};
for contrastIndex = 1:numel(semanticContrastNames)
    weights = zeros(1, numberOfDesignColumns);
    for categoryIndex = 1:numel(categoryLabels)
        columns = categoryColumns{categoryIndex};
        if isempty(columns)
            error('OX:FirstLevelCategory:CategoryAbsent', ...
                '%s has no %s events.', subject, categoryLabels(categoryIndex));
        end
        weights(columns) = contrastTemplate(contrastIndex, categoryIndex) / ...
            numel(columns);
    end

    matlabbatch{3}.spm.stats.con.consess{contrastIndex}.tcon.name = ...
        semanticContrastNames{contrastIndex};
    matlabbatch{3}.spm.stats.con.consess{contrastIndex}.tcon.weights = ...
        weights;
    matlabbatch{3}.spm.stats.con.consess{contrastIndex}.tcon.sessrep = ...
        'none';
end
matlabbatch{3}.spm.stats.con.delete = 1;

matlabbatch{4}.spm.stats.results.spmmat = {spmMat};
for contrastIndex = 1:numel(semanticContrastNames)
    matlabbatch{4}.spm.stats.results.conspec(contrastIndex).titlestr = ...
        semanticContrastNames{contrastIndex};
    matlabbatch{4}.spm.stats.results.conspec(contrastIndex).contrasts = ...
        contrastIndex;
    matlabbatch{4}.spm.stats.results.conspec(contrastIndex).threshdesc = ...
        'FWE';
    matlabbatch{4}.spm.stats.results.conspec(contrastIndex).thresh = fweP;
    matlabbatch{4}.spm.stats.results.conspec(contrastIndex).extent = 0;
    matlabbatch{4}.spm.stats.results.conspec(contrastIndex).conjunction = 1;
    matlabbatch{4}.spm.stats.results.conspec(contrastIndex).mask.none = 1;
end
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.tspm.basename = 'FWE_p001';
end

function [confounds, names] = loadRunConfounds(confoundFile, subject, ...
        runInfo, runIndex, expectedFrames)
requiredFields = {'confounds', 'confound_names', 'nframes', ...
    'run_ordinal', 'session_id', 'run_id'};
if ~isfile(confoundFile)
    error('OX:FirstLevelCategory:MissingConfounds', ...
        'Confound file not found for %s %s: %s', ...
        subject, runInfo.id, confoundFile);
end

data = load(confoundFile, requiredFields{:});
missingFields = requiredFields(~isfield(data, requiredFields));
if ~isempty(missingFields)
    error('OX:FirstLevelCategory:InvalidConfounds', ...
        '%s is missing required field(s): %s', confoundFile, ...
        strjoin(missingFields, ', '));
end

confounds = data.confounds;
names = data.confound_names;
if ~isnumeric(confounds) || size(confounds, 1) ~= expectedFrames || ...
        any(~isfinite(confounds(:)))
    error('OX:FirstLevelCategory:ConfoundFrameMismatch', ...
        '%s must contain a finite %d-by-N confound matrix.', ...
        confoundFile, expectedFrames);
end
if ~iscell(names) || numel(names) ~= size(confounds, 2) || ...
        ~all(cellfun(@(name) ischar(name) || ...
        (isstring(name) && isscalar(name)), names))
    error('OX:FirstLevelCategory:InvalidConfoundNames', ...
        '%s has invalid or mismatched confound_names.', confoundFile);
end

names = cellfun(@char, names(:)', 'UniformOutput', false);
if data.nframes ~= expectedFrames || data.run_ordinal ~= runIndex || ...
        data.session_id ~= runInfo.session || data.run_id ~= runInfo.run
    error('OX:FirstLevelCategory:ConfoundRunMismatch', ...
        '%s metadata does not match %s %s.', ...
        confoundFile, subject, runInfo.id);
end
end
