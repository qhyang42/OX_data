function runConfounds = OX_build_glmsingle_confounds(subjidx, varargin)
%OX_BUILD_GLMSINGLE_CONFOUNDS Build strict run-wise OX nuisance regressors.
%   RUNCONFOUNDS = OX_BUILD_GLMSINGLE_CONFOUNDS(SUBJIDX) validates and
%   builds all 80 runs for one subject, then saves one standalone MAT-file
%   per run under MRI/subj_N/nifti/glmsingle_confounds. Existing files are
%   never overwritten. No fMRI, event, preprocessing, or GLMsingle output
%   is modified.
%
%   Each run contains [motion24, airflow_z, airflow_sq_z, resp_volume_z,
%   one-hot spike columns]. Spikes use FD > 0.4 mm OR run-wise robust DVARS
%   z > 5. Airflow is polarity-corrected, low-pass filtered at 10 Hz, and
%   linearly detrended over the acquisition segment. Squaring occurs at
%   1000 Hz before frame averaging. Respiratory volume is the continuous
%   cumtrapz integral of cleaned airflow divided by 1000 Hz, followed by a
%   linear detrend of the integrated trace before frame averaging.
%
%   Name-value options:
%     ProjectRoot  OX_DATA root (inferred from this file by default)
%     OutputDir    Standalone derivative directory (subject default)
%     ExpectedRuns Required run count (80)
%     SaveOutput   Save per-run files after every run validates (true)
%
%   MRI frame windows use individual TTL rising edges when exactly N or
%   N+1 regularly spaced edges are present. With N edges, the final edge is
%   extrapolated by one TR. Otherwise, the first recorded MRI onset plus
%   the nominal TR is used only when the complete acquisition is covered.
%   If and only if missing coverage is confined to volume N, its three
%   respiratory values copy volume N-1 and volume N receives a forced spike.
%   Missing coverage that reaches any earlier volume is a fatal error.

p = inputParser;
p.addRequired('subjidx', @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x == fix(x));
p.addParameter('ProjectRoot', localProjectRoot(), @(x) ischar(x) || isstring(x));
p.addParameter('OutputDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('ExpectedRuns', 80, @(x) isnumeric(x) && isscalar(x) && x > 0 && x == fix(x));
p.addParameter('SaveOutput', true, @(x) islogical(x) && isscalar(x));
p.parse(subjidx, varargin{:});

projectRoot = char(p.Results.ProjectRoot);
subject = sprintf('subj_%d', subjidx);
outputDir = char(p.Results.OutputDir);
if isempty(outputDir)
    outputDir = fullfile(projectRoot, 'MRI', subject, 'nifti', ...
        'glmsingle_confounds');
end

constants = processingConstants();
runs = OX_discover_functional_runs(subjidx, 'ProjectRoot', projectRoot, ...
    'ExpectedRuns', p.Results.ExpectedRuns);
nRuns = numel(runs);

labchartDir = fullfile(projectRoot, 'labchart', 'extracted_events');
eventFile = fullfile(labchartDir, sprintf('subj%d_events_bm.mat', subjidx));
qcFile = fullfile(projectRoot, 'motion_qc', sprintf('%s_motion_qc.mat', subject));
assertFile(eventFile, 'BreathMetrics event metadata');
assertFile(qcFile, 'validated FD/DVARS QC');

events = load(eventFile);
if ~isfield(events, 'bm_processing') || ...
        ~isfield(events.bm_processing, 'source_raw_file')
    error('OX:GLMsingleConfounds:InvalidBreathMetricsSource', ...
        '%s lacks bm_processing.source_raw_file.', eventFile);
end
rawFile = fullfile(labchartDir, events.bm_processing.source_raw_file);
assertFile(rawFile, 'raw waveform referenced by BreathMetrics metadata');
raw = load(rawFile, 'eventdata');
qc = load(qcFile, 'results');
validateTopLevelSources(raw, events, qc, runs, nRuns, subject);

labchartRuns = raw.eventdata(:);
labchartProvenance = breathMetricsProvenance(events, rawFile, nRuns);
[polarityMultipliers, polarityFlippedRuns, polarityReason] = ...
    OX_get_respiration_polarity(subjidx, nRuns);
if isfield(events.bm_processing, 'respiration_polarity_multiplier_by_run')
    if ~isequal(events.bm_processing.respiration_polarity_multiplier_by_run(:), ...
            polarityMultipliers)
        error('OX:GLMsingleConfounds:PolarityMetadataMismatch', ...
            '%s BreathMetrics polarity metadata disagrees with the shared OX definition.', ...
            subject);
    end
end
[gitCommit, gitDirty] = gitState(projectRoot);

template = emptyRunOutput();
runConfounds = repmat(template, nRuns, 1);
fprintf('Building GLMsingle confounds for %s (%d runs; validation before save)\n', ...
    subject, nRuns);

for runIndex = 1:nRuns
    runInfo = runs(runIndex);
    nframes = functionalNFrames(runInfo.functional_file);
    if nframes ~= events.nframes(runIndex)
        error('OX:GLMsingleConfounds:NFramesOrderMismatch', ...
            ['%s %s has %d functional frames but legacy nframes(%d)=%d. ', ...
             'Run ordering is not safe.'], subject, runInfo.id, nframes, ...
            runIndex, events.nframes(runIndex));
    end

    mp = readmatrix(runInfo.motion_file, 'FileType', 'text');
    validateFiniteMatrix(mp, nframes, 6, 'motion', runInfo.id);
    dmp = [zeros(1, 6); diff(mp, 1, 1)];
    motion24 = [mp, dmp, mp .^ 2, dmp .^ 2];
    validateFiniteMatrix(motion24, nframes, 24, 'motion24', runInfo.id);

    FD = columnVector(qc.results.fd_mm{runIndex});
    robust_DVARS_z = columnVector(qc.results.robust_dvars_z{runIndex});
    if numel(FD) ~= nframes || numel(robust_DVARS_z) ~= nframes
        error('OX:GLMsingleConfounds:QCLengthMismatch', ...
            '%s FD/DVARS lengths (%d/%d) do not equal nframes=%d.', ...
            runInfo.id, numel(FD), numel(robust_DVARS_z), nframes);
    end
    if any(~isfinite(FD)) || any(~isfinite(robust_DVARS_z))
        error('OX:GLMsingleConfounds:NonfiniteQC', ...
            '%s contains nonfinite FD or robust DVARS z.', runInfo.id);
    end

    runData = labchartRuns{runIndex};
    bmMRIOnset = breathMetricsMRIOnset(events, runIndex, constants, runInfo.id);
    validateLabchartRun(runData, events, runIndex, bmMRIOnset, ...
        runInfo.id, constants);
    airflow = double(runData(1, :)) .* polarityMultipliers(runIndex);
    mriTTL = double(runData(3, :));
    if any(~isfinite(airflow))
        error('OX:GLMsingleConfounds:NonfiniteAirflow', ...
            ['%s airflow contains nonfinite samples. No interpolation rule ', ...
             'was authorized, so the run cannot be built.'], runInfo.id);
    end

    [frameEdges, alignment] = frameEdgesFromMRI( ...
        mriTTL, nframes, bmMRIOnset, constants, runInfo.id);
    filteredAirflow = lowpass(airflow, constants.lowpass_cutoff_hz, ...
        constants.labchart_sampling_rate_hz);
    if any(~isfinite(filteredAirflow))
        error('OX:GLMsingleConfounds:NonfiniteFilteredAirflow', ...
            '%s low-pass filtering produced nonfinite values.', runInfo.id);
    end

    acquisitionEnd = min(frameEdges(end) - 1, numel(mriTTL));
    acquisitionSamples = frameEdges(1):acquisitionEnd;
    cleanedAirflow = detrend(filteredAirflow(acquisitionSamples(:)), 'linear');
    airflowSquared = cleanedAirflow .^ 2;
    respVolume = cumtrapz(cleanedAirflow) ./ constants.labchart_sampling_rate_hz;
    respVolume = detrend(respVolume, 'linear');

    [airflow_TR_raw, airflow_sq_TR_raw, resp_volume_TR_raw] = ...
        averageInFrameWindows(cleanedAirflow, airflowSquared, respVolume, ...
        frameEdges, nframes, alignment.final_volume_breath_samples_incomplete);
    physiologyRaw = [airflow_TR_raw, airflow_sq_TR_raw, resp_volume_TR_raw];
    validatePhysiology(physiologyRaw, nframes, runInfo.id, 'pre-z-scored');
    physiologyZ = zscoreColumns(physiologyRaw, runInfo.id);
    validatePhysiology(physiologyZ, nframes, runInfo.id, 'standardized');

    airflow_TR_z = physiologyZ(:, 1);
    airflow_sq_TR_z = physiologyZ(:, 2);
    resp_volume_TR_z = physiologyZ(:, 3);

    fdFlag = FD > constants.fd_threshold_mm;
    dvarsFlag = robust_DVARS_z > constants.robust_dvars_z_threshold;
    respirationFlag = false(nframes, 1);
    if alignment.final_volume_breath_samples_incomplete
        respirationFlag(end) = true;
    end
    bad_volume = fdFlag | dvarsFlag | respirationFlag;
    badIndices = find(bad_volume);
    spikeRegressors = zeros(nframes, numel(badIndices));
    if ~isempty(badIndices)
        spikeRegressors(sub2ind(size(spikeRegressors), ...
            badIndices, (1:numel(badIndices))')) = 1;
    end
    bad_volume_reason = repmat({''}, nframes, 1);
    bad_volume_reason(fdFlag & ~dvarsFlag) = {'FD only'};
    bad_volume_reason(~fdFlag & dvarsFlag) = {'DVARS only'};
    bad_volume_reason(fdFlag & dvarsFlag) = {'both'};
    if respirationFlag(end)
        if fdFlag(end) && dvarsFlag(end)
            bad_volume_reason{end} = 'FD + DVARS + incomplete respiration';
        elseif fdFlag(end)
            bad_volume_reason{end} = 'FD + incomplete respiration';
        elseif dvarsFlag(end)
            bad_volume_reason{end} = 'DVARS + incomplete respiration';
        else
            bad_volume_reason{end} = 'incomplete respiration';
        end
    end
    validateSpikes(spikeRegressors, bad_volume, badIndices, runInfo.id);

    confounds = [motion24, physiologyZ, spikeRegressors];
    confound_names = [motionNames(), ...
        {'airflow_z', 'airflow_sq_z', 'resp_volume_z'}, ...
        arrayfun(@(x) sprintf('spike_vol_%04d', x), badIndices', ...
            'UniformOutput', false)];
    validateFinalConfounds(confounds, confound_names, motion24, ...
        physiologyRaw, FD, robust_DVARS_z, nframes, badIndices, runInfo.id);

    [designRank, conditionNumber, singularValues] = designDiagnostics(confounds);
    processing_metadata = makeMetadata(constants, runInfo, rawFile, ...
        eventFile, qcFile, labchartProvenance(runIndex), ...
        polarityMultipliers(runIndex), polarityFlippedRuns, polarityReason, ...
        alignment, gitCommit, gitDirty, designRank, conditionNumber, ...
        singularValues);

    out = template;
    out.confounds = confounds;
    out.confound_names = confound_names;
    out.motion24 = motion24;
    out.airflow_TR_raw = airflow_TR_raw;
    out.airflow_sq_TR_raw = airflow_sq_TR_raw;
    out.resp_volume_TR_raw = resp_volume_TR_raw;
    out.airflow_TR_z = airflow_TR_z;
    out.airflow_sq_TR_z = airflow_sq_TR_z;
    out.resp_volume_TR_z = resp_volume_TR_z;
    out.FD = FD;
    out.robust_DVARS_z = robust_DVARS_z;
    out.bad_volume = bad_volume;
    out.bad_volume_reason = bad_volume_reason;
    out.bad_volume_indices = badIndices;
    out.respiratory_incomplete_volume = respirationFlag;
    out.run_id = runInfo.run;
    out.session_id = runInfo.session;
    out.run_ordinal = runIndex;
    out.run_label = runInfo.id;
    out.nframes = nframes;
    out.processing_metadata = processing_metadata;
    runConfounds(runIndex) = out;

    if runIndex == 1 || mod(runIndex, 10) == 0 || runIndex == nRuns
        fprintf('  %s: %d/%d runs validated\n', subject, runIndex, nRuns);
    end
end

if p.Results.SaveOutput
    saveOutputsAtomically(runConfounds, outputDir, subject);
    fprintf('Saved %d standalone confound files in %s\n', nRuns, outputDir);
end
end

function constants = processingConstants()
constants = struct();
constants.TR_seconds = 0.76;
constants.labchart_sampling_rate_hz = 1000;
constants.samples_per_TR = 760;
constants.lowpass_cutoff_hz = 10;
constants.lowpass_function = 'MATLAB lowpass';
constants.airflow_detrending = [ ...
    'MATLAB detrend(x,''linear'') over the complete acquisition-aligned ', ...
    '1-kHz segment after 10-Hz low-pass filtering'];
constants.integration = [ ...
    'cumtrapz(cleaned_airflow)/1000; then MATLAB ', ...
    'detrend(resp_volume,''linear'') over the continuous run trace; ', ...
    'the integral is never reset at breaths'];
constants.fd_threshold_mm = 0.4;
constants.robust_dvars_z_threshold = 5;
constants.mri_ttl_threshold = 1;
constants.event_ttl_threshold = 0.1;
constants.expected_event_ttl_rises = 23;
constants.trigger_interval_tolerance_samples = 2;
end

function out = emptyRunOutput()
out = struct('confounds', [], 'confound_names', {{}}, 'motion24', [], ...
    'airflow_TR_raw', [], 'airflow_sq_TR_raw', [], ...
    'resp_volume_TR_raw', [], 'airflow_TR_z', [], ...
    'airflow_sq_TR_z', [], 'resp_volume_TR_z', [], 'FD', [], ...
    'robust_DVARS_z', [], 'bad_volume', [], ...
    'bad_volume_reason', {{}}, 'bad_volume_indices', [], ...
    'respiratory_incomplete_volume', [], ...
    'run_id', [], 'session_id', [], 'run_ordinal', [], ...
    'run_label', '', 'nframes', [], 'processing_metadata', struct());
end

function validateTopLevelSources(raw, events, qc, runs, nRuns, subject)
if ~isfield(raw, 'eventdata') || ~iscell(raw.eventdata) || ...
        numel(raw.eventdata) ~= nRuns
    error('OX:GLMsingleConfounds:InvalidRawLabChart', ...
        '%s raw LabChart file must contain exactly %d eventdata cells.', ...
        subject, nRuns);
end
requiredEventFields = {'nframes', 'sessi', 'sessf', 'event_onsets', ...
    'cue_onsets', 'bm_processing', 'bm_sniff_ttl_samples', ...
    'bm_matched_inhale_samples'};
if ~all(isfield(events, requiredEventFields)) || numel(events.nframes) ~= nRuns
    error('OX:GLMsingleConfounds:InvalidEventMetadata', ...
        '%s event metadata is missing required run-wise fields.', subject);
end
if any(~isfinite(events.nframes)) || any(events.nframes ~= fix(events.nframes))
    error('OX:GLMsingleConfounds:InvalidNFrames', ...
        '%s event nframes must be finite positive integers.', subject);
end
if events.sessi ~= min([runs.session]) || events.sessf ~= max([runs.session])
    error('OX:GLMsingleConfounds:SessionRangeMismatch', ...
        '%s event session range %d-%d differs from functional range %d-%d.', ...
        subject, events.sessi, events.sessf, min([runs.session]), max([runs.session]));
end
if ~isfield(qc, 'results') || numel(qc.results.run_ids) ~= nRuns
    error('OX:GLMsingleConfounds:InvalidMotionQC', ...
        '%s motion QC file does not contain %d run results.', subject, nRuns);
end
if ~isequal(qc.results.run_ids(:), {runs.id}')
    error('OX:GLMsingleConfounds:QCRunOrderMismatch', ...
        '%s QC run IDs do not exactly match functional discovery order.', subject);
end
if ~isequal(qc.results.n_volumes_per_run(:), events.nframes(:))
    error('OX:GLMsingleConfounds:QCNFramesMismatch', ...
        '%s QC and event nframes vectors do not agree exactly.', subject);
end
if ~isfield(events.bm_processing, 'source_raw_file')
    error('OX:GLMsingleConfounds:IncompleteBreathMetricsMetadata', ...
        '%s BreathMetrics processing metadata is incomplete.', subject);
end
if ~isequal(size(events.bm_sniff_ttl_samples), [10, nRuns]) || ...
        ~isequal(size(events.bm_matched_inhale_samples), [10, nRuns]) || ...
        ~isequal(size(events.event_onsets), [10, nRuns])
    error('OX:GLMsingleConfounds:InvalidBreathMetricsDimensions', ...
        '%s BreathMetrics run/trial arrays must be 10-by-%d.', subject, nRuns);
end
end

function provenance = breathMetricsProvenance(events, rawFile, nRuns)
template = struct('raw_eventdata_ordinal', [], ...
    'source_breathmetrics_file', events.bm_processing.source_raw_file, ...
    'ordering_method', ['subjN_events_bm.mat run ordinal is authoritative; ', ...
    'raw waveform cell identity is validated against saved BreathMetrics TTL samples']);
provenance = repmat(template, nRuns, 1);
for runIndex = 1:nRuns
    provenance(runIndex).raw_eventdata_ordinal = runIndex;
    provenance(runIndex).source_breathmetrics_file = rawFile;
end
end

function onset = breathMetricsMRIOnset(events, runIndex, constants, runID)
candidate = events.bm_matched_inhale_samples(:, runIndex) - ...
    events.event_onsets(:, runIndex) .* constants.labchart_sampling_rate_hz;
if any(~isfinite(candidate)) || max(abs(candidate - candidate(1))) > 1e-6 || ...
        abs(candidate(1) - round(candidate(1))) > 1e-6
    error('OX:GLMsingleConfounds:InvalidBreathMetricsMRIOnset', ...
        '%s cannot recover one integer MRI onset from BreathMetrics metadata.', runID);
end
onset = round(candidate(1));
end

function validateLabchartRun(runData, events, runIndex, bmMRIOnset, runID, constants)
if ~isnumeric(runData) || size(runData, 1) ~= 3 || isempty(runData)
    error('OX:GLMsingleConfounds:InvalidLabChartRun', ...
        '%s LabChart data must be a numeric 3-by-N matrix.', runID);
end
if any(~isfinite(runData(:)))
    error('OX:GLMsingleConfounds:NonfiniteLabChart', ...
        '%s LabChart data contains NaN or Inf.', runID);
end
eventRises = find(diff([false, runData(2, :) > ...
    constants.event_ttl_threshold]) == 1);
if numel(eventRises) ~= constants.expected_event_ttl_rises
    error('OX:GLMsingleConfounds:UnexpectedEventTTLCount', ...
        '%s has %d event TTL rises; expected %d.', runID, ...
        numel(eventRises), constants.expected_event_ttl_rises);
end
if ~isequal(eventRises(5:2:end)', ...
        events.bm_sniff_ttl_samples(:, runIndex))
    error('OX:GLMsingleConfounds:BreathMetricsRawIdentityMismatch', ...
        ['%s raw waveform cell does not match bm_sniff_ttl_samples for ', ...
         'the same authoritative BreathMetrics run ordinal.'], runID);
end
mriOnset = find(runData(3, :) > constants.mri_ttl_threshold, 1, 'first');
if isempty(mriOnset) || mriOnset ~= bmMRIOnset
    error('OX:GLMsingleConfounds:BreathMetricsMRIOnsetMismatch', ...
        '%s raw MRI onset does not match the onset encoded by BreathMetrics.', runID);
end
end

function [edges, info] = frameEdgesFromMRI(mriTTL, nframes, ...
        authoritativeOnset, constants, runID)
rises = find(diff([false, mriTTL > constants.mri_ttl_threshold]) == 1);
if isempty(rises)
    error('OX:GLMsingleConfounds:MissingMRIOnset', ...
        '%s has no MRI TTL rise above %.3g.', runID, constants.mri_ttl_threshold);
end
regular = numel(rises) > 1 && all(abs(diff(rises) - ...
    constants.samples_per_TR) <= constants.trigger_interval_tolerance_samples);
if numel(rises) == nframes + 1 && regular
    edges = rises(1:(nframes + 1));
    method = 'individual MRI TTL boundaries (N+1 validated rising edges)';
    usedIndividual = true;
elseif numel(rises) == nframes && regular
    edges = [rises(:); rises(end) + constants.samples_per_TR]';
    method = ['individual MRI TTL onsets (N validated rising edges); ', ...
        'final boundary extrapolated by one nominal TR'];
    usedIndividual = true;
else
    edges = authoritativeOnset + (0:nframes) .* constants.samples_per_TR;
    method = ['BreathMetrics-authoritative first MRI TTL onset plus nominal TR; individual ', ...
        'trigger count did not validate against nframes'];
    usedIndividual = false;
end
edges = round(edges(:)');
if any(diff(edges) <= 0)
    error('OX:GLMsingleConfounds:InvalidFrameEdges', ...
        '%s produced non-increasing MRI frame edges.', runID);
end
if edges(1) < 1
    error('OX:GLMsingleConfounds:InvalidFirstFrameEdge', ...
        '%s produced a frame edge before the first LabChart sample.', runID);
end
missing = max(0, edges(end) - 1 - numel(mriTTL));
incompleteFinal = missing > 0;
if incompleteFinal && (nframes < 2 || edges(end - 1) > numel(mriTTL) + 1)
    error('OX:GLMsingleConfounds:IncompleteRespiratoryCoverage', ...
        ['%s has nframes=%d and %d MRI TTL rises, and the LabChart ', ...
         'shortage reaches a volume before the final volume. Refusing ', ...
         'to replace more than the final respiratory value.'], ...
        runID, nframes, numel(rises));
end
info = struct();
info.method = method;
info.used_individual_triggers = usedIndividual;
info.mri_ttl_threshold = constants.mri_ttl_threshold;
info.n_detected_rising_edges = numel(rises);
info.first_mri_onset_sample = edges(1);
info.frame_edge_samples = edges;
info.final_volume_breath_samples_incomplete = incompleteFinal;
info.final_volume_missing_samples = missing;
info.final_volume_available_samples = max(0, ...
    numel(mriTTL) - edges(end - 1) + 1);
info.final_volume_expected_samples = edges(end) - edges(end - 1);
if incompleteFinal
    info.final_volume_replacement_source_volume = nframes - 1;
    info.final_volume_replacement_method = [ ...
        'Copy all three pre-z-scored respiratory values from volume N-1 ', ...
        'into volume N; force a one-hot bad-volume regressor for volume N'];
else
    info.final_volume_replacement_source_volume = [];
    info.final_volume_replacement_method = '';
end
info.trigger_interval_tolerance_samples = ...
    constants.trigger_interval_tolerance_samples;
if numel(rises) > 1
    info.detected_interval_range_samples = [min(diff(rises)), max(diff(rises))];
else
    info.detected_interval_range_samples = [NaN, NaN];
end
end

function [a, a2, v] = averageInFrameWindows(airflow, airflowSquared, ...
        respVolume, frameEdges, nframes, incompleteFinal)
a = nan(nframes, 1);
a2 = nan(nframes, 1);
v = nan(nframes, 1);
origin = frameEdges(1);
lastMeasuredFrame = nframes - double(incompleteFinal);
for frame = 1:lastMeasuredFrame
    firstSample = frameEdges(frame) - origin + 1;
    lastSample = frameEdges(frame + 1) - origin;
    index = firstSample:lastSample;
    if isempty(index)
        error('OX:GLMsingleConfounds:EmptyFrameWindow', ...
            'MRI frame %d has an empty respiratory window.', frame);
    end
    a(frame) = mean(airflow(index));
    a2(frame) = mean(airflowSquared(index));
    v(frame) = mean(respVolume(index));
end
if incompleteFinal
    a(end) = a(end - 1);
    a2(end) = a2(end - 1);
    v(end) = v(end - 1);
end
end

function validatePhysiology(physiology, nframes, runID, stage)
if ~isequal(size(physiology), [nframes, 3]) || any(~isfinite(physiology(:)))
    error('OX:GLMsingleConfounds:InvalidPhysiology', ...
        '%s %s physiology must be a finite nframes-by-3 matrix.', runID, stage);
end
scales = std(physiology, 0, 1);
if any(scales <= 1e-12)
    error('OX:GLMsingleConfounds:ConstantPhysiology', ...
        '%s has a constant/degenerate %s physiology column.', runID, stage);
end
c = corr(physiology);
offDiagonal = abs(c(triu(true(3), 1)));
if any(~isfinite(c(:))) || any(offDiagonal > 1 - 1e-10) || rank(physiology) < 3
    error('OX:GLMsingleConfounds:DegeneratePhysiology', ...
        '%s respiratory regressors are identical or numerically degenerate.', runID);
end
end

function z = zscoreColumns(x, runID)
mu = mean(x, 1);
sigma = std(x, 0, 1);
if any(~isfinite(sigma)) || any(sigma <= 1e-12)
    error('OX:GLMsingleConfounds:CannotStandardize', ...
        '%s has a nonfinite or zero physiology scale.', runID);
end
z = (x - mu) ./ sigma;
end

function validateSpikes(spikes, badVolume, badIndices, runID)
if size(spikes, 2) ~= numel(unique(badIndices)) || ...
        size(spikes, 2) ~= nnz(badVolume)
    error('OX:GLMsingleConfounds:SpikeCountMismatch', ...
        '%s spike count does not equal the unique bad-volume count.', runID);
end
if ~isempty(spikes)
    if any(sum(spikes, 1) ~= 1) || any(sum(spikes, 2) > 1) || ...
            size(unique(spikes', 'rows'), 1) ~= size(spikes, 2)
        error('OX:GLMsingleConfounds:InvalidSpikes', ...
            '%s spike columns are not unique one-hot regressors.', runID);
    end
end
end

function validateFinalConfounds(confounds, names, motion24, physiologyRaw, ...
        FD, robustDvarsZ, nframes, badIndices, runID)
if size(confounds, 1) ~= nframes || size(motion24, 1) ~= nframes || ...
        numel(FD) ~= nframes || numel(robustDvarsZ) ~= nframes || ...
        size(physiologyRaw, 1) ~= nframes
    error('OX:GLMsingleConfounds:FinalLengthMismatch', ...
        '%s failed one or more required nframes assertions.', runID);
end
if size(confounds, 2) ~= 27 + numel(badIndices) || ...
        numel(names) ~= size(confounds, 2)
    error('OX:GLMsingleConfounds:FinalColumnMismatch', ...
        '%s must have 27 plus one column per bad volume.', runID);
end
if any(~isfinite(confounds(:)))
    error('OX:GLMsingleConfounds:NonfiniteFinalConfounds', ...
        '%s final confounds contain NaN or Inf.', runID);
end
continuousScale = std(confounds(:, 1:27), 0, 1);
if any(continuousScale <= 1e-12)
    badColumns = find(continuousScale <= 1e-12);
    error('OX:GLMsingleConfounds:ConstantContinuousColumn', ...
        '%s has unexpected constant continuous column(s): %s.', ...
        runID, mat2str(badColumns));
end
end

function [r, kappa, singularValues] = designDiagnostics(x)
singularValues = svd(x, 'econ');
if isempty(singularValues)
    r = 0;
    kappa = NaN;
    return
end
tolerance = max(size(x)) * eps(max(singularValues));
r = nnz(singularValues > tolerance);
if r < size(x, 2) || singularValues(end) == 0
    kappa = Inf;
else
    kappa = singularValues(1) / singularValues(end);
end
end

function metadata = makeMetadata(constants, runInfo, rawFile, eventFile, ...
        qcFile, provenance, polarityMultiplier, flippedRuns, polarityReason, ...
        alignment, gitCommit, gitDirty, designRank, conditionNumber, singularValues)
metadata = struct();
metadata.created_at = char(datetime('now', 'TimeZone', 'local', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
metadata.TR_seconds = constants.TR_seconds;
metadata.labchart_sampling_rate_hz = constants.labchart_sampling_rate_hz;
metadata.airflow_lowpass_cutoff_hz = constants.lowpass_cutoff_hz;
metadata.airflow_lowpass_function = constants.lowpass_function;
metadata.airflow_detrending = constants.airflow_detrending;
metadata.integration_and_volume_detrending = constants.integration;
metadata.squared_airflow_order = ...
    'cleaned 1-kHz airflow squared before averaging within MRI windows';
metadata.physiology_standardization = ...
    'independent within-run arithmetic z-score after MRI-frame averaging';
metadata.HRF_convolution_applied = false;
metadata.fd_threshold_mm = constants.fd_threshold_mm;
metadata.fd_operator = '>';
metadata.robust_dvars_z_threshold = constants.robust_dvars_z_threshold;
metadata.robust_dvars_operator = '>';
metadata.bad_volume_rule = [ ...
    'FD > 0.4 OR robust_DVARS_z > 5 OR final-volume respiratory ', ...
    'coverage replacement'];
metadata.respiration_polarity_multiplier = polarityMultiplier;
metadata.polarity_flipped_runs = flippedRuns;
metadata.polarity_flip_reason = polarityReason;
metadata.mri_frame_alignment = alignment;
metadata.labchart_run_provenance = provenance;
metadata.source_files = struct('functional', runInfo.functional_file, ...
    'motion', runInfo.motion_file, 'raw_labchart', rawFile, ...
    'event_metadata', eventFile, 'motion_qc', qcFile);
metadata.matlab_version = version;
metadata.git_commit = gitCommit;
metadata.git_dirty = gitDirty;
metadata.design_rank = designRank;
metadata.design_n_columns = numel(singularValues);
metadata.design_condition_number = conditionNumber;
metadata.design_singular_values = singularValues;
metadata.processing_code = mfilename('fullpath');
end

function saveOutputsAtomically(outputs, outputDir, subject)
targets = cell(numel(outputs), 1);
for runIndex = 1:numel(outputs)
    targets{runIndex} = fullfile(outputDir, sprintf( ...
        '%s_session%02d_run%02d_confounds.mat', subject, ...
        outputs(runIndex).session_id, outputs(runIndex).run_id));
end
existing = targets(cellfun(@isfile, targets));
if ~isempty(existing)
    error('OX:GLMsingleConfounds:WouldOverwrite', ...
        'Refusing to overwrite existing derivative: %s', existing{1});
end
if ~isfolder(outputDir)
    [ok, message] = mkdir(outputDir);
    if ~ok
        error('OX:GLMsingleConfounds:CannotCreateOutput', ...
            'Cannot create %s: %s', outputDir, message);
    end
end
for runIndex = 1:numel(outputs)
    out = outputs(runIndex);
    temporary = [tempname(outputDir), '.mat'];
    try
        save(temporary, '-struct', 'out', '-v7.3');
        [ok, message] = movefile(temporary, targets{runIndex});
        if ~ok
            error('OX:GLMsingleConfounds:MoveFailed', ...
                'Cannot move output to %s: %s', targets{runIndex}, message);
        end
    catch saveError
        if isfile(temporary)
            delete(temporary);
        end
        rethrow(saveError);
    end
end
end

function names = motionNames()
base = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
names = [base, strcat('d_', base), strcat(base, '_sq'), ...
    strcat('d_', base, '_sq')];
end

function n = functionalNFrames(filename)
info = niftiinfo(filename);
if numel(info.ImageSize) < 4
    error('OX:GLMsingleConfounds:NotFourDimensional', ...
        'Functional image is not 4-D: %s', filename);
end
n = info.ImageSize(4);
end

function validateFiniteMatrix(x, nRows, nColumns, label, runID)
if ~isnumeric(x) || ~isequal(size(x), [nRows, nColumns]) || any(~isfinite(x(:)))
    error('OX:GLMsingleConfounds:InvalidMatrix', ...
        '%s %s must be a finite %d-by-%d matrix.', runID, label, nRows, nColumns);
end
end

function x = columnVector(x)
x = x(:);
end

function [commit, dirty] = gitState(projectRoot)
[commitStatus, commitText] = system(sprintf( ...
    'git -C "%s" rev-parse HEAD 2>/dev/null', projectRoot));
if commitStatus == 0
    commit = strtrim(commitText);
else
    commit = '';
end
[dirtyStatus, dirtyText] = system(sprintf( ...
    'git -C "%s" status --porcelain 2>/dev/null', projectRoot));
dirty = dirtyStatus ~= 0 || ~isempty(strtrim(dirtyText));
end

function assertFile(filename, label)
if ~isfile(filename)
    error('OX:GLMsingleConfounds:MissingSource', ...
        'Missing %s file: %s', label, filename);
end
end

function root = localProjectRoot()
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
