%% Run the sniff-only first-level SPM model for OX subjects 2-6.
% Sniff events use BreathMetrics inhale onsets. Each functional run is an
% independent SPM session with its own motion, respiratory, and bad-volume
% nuisance regressors. Rest is represented by SPM's implicit baseline.

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
    outputDir = fullfile(niftiDir, 'first_level_model_sniff_physio');

    events = loadBreathMetricsEvents(eventFile, subject);
    runs = OX_discover_functional_runs(subjidx, ...
        'ProjectRoot', info.projectRoot, ...
        'ExpectedRuns', numel(events.nframes));

    matlabbatch = buildFirstLevelBatch(subject, runs, events, niftiDir, ...
        outputDir, TR, highPassSeconds, fweP);

    if ~isfolder(outputDir)
        mkdir(outputDir);
    end

    fprintf('Running sniff-only first-level model for %s (%d runs).\n', ...
        subject, numel(runs));
    spm_jobman('run', matlabbatch);
end


%% 
function events = loadBreathMetricsEvents(eventFile, subject)
requiredFields = {'event_onsets', 'nframes', 'bm_processing'};
if ~isfile(eventFile)
    error('OX:FirstLevel:MissingBreathMetricsEvents', ...
        'BreathMetrics event file not found for %s: %s', subject, eventFile);
end

events = load(eventFile, requiredFields{:});
missingFields = requiredFields(~isfield(events, requiredFields));
if ~isempty(missingFields)
    error('OX:FirstLevel:InvalidBreathMetricsEvents', ...
        '%s is missing required field(s): %s', eventFile, ...
        strjoin(missingFields, ', '));
end

if ~isnumeric(events.event_onsets) || ~isnumeric(events.nframes) || ...
        isempty(events.event_onsets) || isempty(events.nframes) || ...
        any(~isfinite(events.event_onsets(:))) || ...
        any(~isfinite(events.nframes(:)))
    error('OX:FirstLevel:InvalidBreathMetricsEvents', ...
        '%s contains invalid event_onsets or nframes.', eventFile);
end

events.nframes = events.nframes(:);
if size(events.event_onsets, 2) ~= numel(events.nframes)
    error('OX:FirstLevel:EventRunCountMismatch', ...
        ['%s contains %d event-onset run(s), but nframes contains ', ...
         '%d run(s).'], eventFile, size(events.event_onsets, 2), ...
        numel(events.nframes));
end
if any(events.nframes <= 0) || any(events.nframes ~= fix(events.nframes))
    error('OX:FirstLevel:InvalidFrameCounts', ...
        '%s contains invalid run frame counts.', eventFile);
end
end

function matlabbatch = buildFirstLevelBatch(subject, runs, events, ...
        niftiDir, outputDir, TR, highPassSeconds, fweP)
numberOfRuns = numel(runs);
if numberOfRuns ~= numel(events.nframes)
    error('OX:FirstLevel:RunCountMismatch', ...
        '%s has %d discovered runs but %d event runs.', subject, ...
        numberOfRuns, numel(events.nframes));
end

matlabbatch = cell(1, 4);
matlabbatch{1}.spm.stats.fmri_spec.dir = {outputDir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 30;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;

for runIndex = 1:numberOfRuns
    runInfo = runs(runIndex);
    scanList = cellstr(spm_select('expand', runInfo.functional_file));
    expectedFrames = events.nframes(runIndex);
    if numel(scanList) ~= expectedFrames
        error('OX:FirstLevel:FunctionalFrameMismatch', ...
            '%s %s has %d functional frames; expected %d.', subject, ...
            runInfo.id, numel(scanList), expectedFrames);
    end

    onsets = events.event_onsets(:, runIndex);
    runDuration = expectedFrames * TR;
    if any(onsets < 0) || any(onsets >= runDuration)
        error('OX:FirstLevel:OnsetOutsideRun', ...
            '%s %s has a sniff onset outside [0, %.3f) seconds.', ...
            subject, runInfo.id, runDuration);
    end

    confoundFile = fullfile(niftiDir, 'glmsingle_confounds', sprintf( ...
        '%s_session%02d_run%02d_confounds.mat', subject, ...
        runInfo.session, runInfo.run));
    [confounds, confoundNames] = loadRunConfounds( ...
        confoundFile, subject, runInfo, runIndex, expectedFrames);

    session = struct();
    session.scans = scanList;
    session.cond.name = 'Odor';
    session.cond.onset = onsets(:);
    session.cond.duration = 0;
    session.cond.tmod = 0;
    session.cond.pmod = struct('name', {}, 'param', {}, 'poly', {});
    session.cond.orth = 1;
    session.multi = {''};
    session.regress = struct('name', confoundNames, ...
        'val', num2cell(confounds, 1));
    session.multi_reg = {''};
    session.hpf = highPassSeconds;
    matlabbatch{1}.spm.stats.fmri_spec.sess(runIndex) = session;
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

matlabbatch{3}.spm.stats.con.spmmat = {spmMat};
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Odor > Rest';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';
matlabbatch{3}.spm.stats.con.delete = 1;

matlabbatch{4}.spm.stats.results.spmmat = {spmMat};
matlabbatch{4}.spm.stats.results.conspec.titlestr = 'Odor > Rest';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = fweP;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.tspm.basename = 'FWE_p001';
end

function [confounds, names] = loadRunConfounds(confoundFile, subject, ...
        runInfo, runIndex, expectedFrames)
requiredFields = {'confounds', 'confound_names', 'nframes', ...
    'run_ordinal', 'session_id', 'run_id'};
if ~isfile(confoundFile)
    error('OX:FirstLevel:MissingConfounds', ...
        'Confound file not found for %s %s: %s', ...
        subject, runInfo.id, confoundFile);
end

data = load(confoundFile, requiredFields{:});
missingFields = requiredFields(~isfield(data, requiredFields));
if ~isempty(missingFields)
    error('OX:FirstLevel:InvalidConfounds', ...
        '%s is missing required field(s): %s', confoundFile, ...
        strjoin(missingFields, ', '));
end

confounds = data.confounds;
names = data.confound_names;
if ~isnumeric(confounds) || size(confounds, 1) ~= expectedFrames || ...
        any(~isfinite(confounds(:)))
    error('OX:FirstLevel:ConfoundFrameMismatch', ...
        '%s must contain a finite %d-by-N confound matrix.', ...
        confoundFile, expectedFrames);
end
if ~iscell(names) || numel(names) ~= size(confounds, 2) || ...
        ~all(cellfun(@(name) ischar(name) || ...
        (isstring(name) && isscalar(name)), names))
    error('OX:FirstLevel:InvalidConfoundNames', ...
        '%s has invalid or mismatched confound_names.', confoundFile);
end

names = cellfun(@char, names(:)', 'UniformOutput', false);
if data.nframes ~= expectedFrames || data.run_ordinal ~= runIndex || ...
        data.session_id ~= runInfo.session || data.run_id ~= runInfo.run
    error('OX:FirstLevel:ConfoundRunMismatch', ...
        '%s metadata does not match %s %s.', ...
        confoundFile, subject, runInfo.id);
end
end
