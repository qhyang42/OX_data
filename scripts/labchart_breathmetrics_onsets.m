%% Re-estimate odor onsets from respiration using BreathMetrics
% This script preserves each subject's legacy event MAT-file structure,
% replaces the odor onset arrays with the nearest detected inhale onsets,
% and adds BreathMetrics matching QC fields.

scriptFile = mfilename('fullpath');
projectRoot = fileparts(fileparts(scriptFile));
labchartDir = fullfile(projectRoot, 'labchart');
breathmetricsDir = '~/Documents/Projects/breathmetrics/';

% Optionally define subjectNumbersToProcess before running this script to
% process only selected subjects, for example:
% subjectNumbersToProcess = 3;
if ~exist('subjectNumbersToProcess', 'var')
    subjectNumbersToProcess = [];
end

addpath(genpath(breathmetricsDir));
if isempty(which('breathmetrics'))
    error('OX_DATA:MissingBreathMetrics', ...
        'BreathMetrics was not found after adding: %s', breathmetricsDir);
end

samplingRate = 1000;
lowpassCutoff = 10;
dataType = 'humanAirflow';
eventThreshold = 0.1;
mriThreshold = 1;
expectedTTLRises = 23;
expectedTrials = 10;

% Subject 3 session 1 (concatenated runs 1-10) was acquired with reversed
% airflow polarity. Flip it so positive flow consistently represents inhale.
polarityFlipRunsBySubject = struct();
polarityFlipRunsBySubject.subj3 = 1:10;

rawFiles = dir(fullfile(labchartDir, 'subj*_events_raw.mat'));
rawFiles = rawFiles(~[rawFiles.isdir]);
if isempty(rawFiles)
    error('OX_DATA:NoRawLabChartFiles', ...
        'No subj*_events_raw.mat files were found in %s.', labchartDir);
end

if ~isempty(subjectNumbersToProcess)
    requestedNumbers = unique(subjectNumbersToProcess(:)');
    keepFile = false(size(rawFiles));
    for rawFileIndex = 1:numel(rawFiles)
        token = regexp(rawFiles(rawFileIndex).name, ...
            '^subj(\d+)_events_raw\.mat$', 'tokens', 'once');
        keepFile(rawFileIndex) = ~isempty(token) && ...
            ismember(str2double(token{1}), requestedNumbers);
    end
    rawFiles = rawFiles(keepFile);
    foundNumbers = arrayfun(@(file) str2double(regexp(file.name, ...
        '\d+', 'match', 'once')), rawFiles);
    missingNumbers = setdiff(requestedNumbers, foundNumbers);
    if ~isempty(missingNumbers)
        error('OX_DATA:RequestedRawFilesMissing', ...
            'Raw event files were not found for subject(s): %s', ...
            mat2str(missingNumbers));
    end
end

fprintf('Found %d raw LabChart subject files.\n', numel(rawFiles));

for subjectIndex = 1:numel(rawFiles)
    rawFile = rawFiles(subjectIndex);
    subjectToken = regexp(rawFile.name, ...
        '^subj(\d+)_events_raw\.mat$', 'tokens', 'once');
    if isempty(subjectToken)
        error('OX_DATA:UnexpectedRawFilename', ...
            'Unexpected raw filename: %s', rawFile.name);
    end

    subjectNumber = subjectToken{1};
    subjectName = ['subj', subjectNumber];
    legacyName = ['subj_', subjectNumber, '_events.mat'];
    outputName = [subjectName, '_events_bm.mat'];
    rawPath = fullfile(rawFile.folder, rawFile.name);
    legacyPath = fullfile(labchartDir, legacyName);
    outputPath = fullfile(labchartDir, outputName);

    if ~isfile(legacyPath)
        error('OX_DATA:MissingLegacyEvents', ...
            'Legacy event template not found: %s', legacyPath);
    end

    rawData = load(rawPath, 'eventdata');
    if ~isfield(rawData, 'eventdata') || ~iscell(rawData.eventdata)
        error('OX_DATA:InvalidRawEvents', ...
            '%s must contain a cell array named eventdata.', rawFile.name);
    end

    legacyData = load(legacyPath);
    validateLegacyTemplate(legacyData, numel(rawData.eventdata), legacyName, ...
        expectedTrials);

    nruns = numel(rawData.eventdata);
    eventOnsets = nan(expectedTrials, nruns);
    allInhaleOnsets = cell(nruns, 1);
    sniffTTLSamples = nan(expectedTrials, nruns);
    matchedInhaleSamples = nan(expectedTrials, nruns);
    matchDeltaSeconds = nan(expectedTrials, nruns);
    polarityMultipliers = ones(nruns, 1);
    if isfield(polarityFlipRunsBySubject, subjectName)
        polarityFlipRuns = polarityFlipRunsBySubject.(subjectName);
    else
        polarityFlipRuns = [];
    end
    if any(polarityFlipRuns < 1 | polarityFlipRuns > nruns | ...
            polarityFlipRuns ~= fix(polarityFlipRuns))
        error('OX_DATA:InvalidPolarityFlipRuns', ...
            'Invalid polarity-flip run indices configured for %s.', subjectName);
    end
    polarityMultipliers(polarityFlipRuns) = -1;

    fprintf('\n%s: processing %d runs\n', subjectName, nruns);
    for runIndex = 1:nruns
        runData = rawData.eventdata{runIndex};
        if ~isnumeric(runData) || size(runData, 1) ~= 3
            error('OX_DATA:InvalidRunData', ...
                '%s run %d must be a numeric 3-by-N array.', ...
                rawFile.name, runIndex);
        end

        respiration = runData(1, :);
        eventTTL = runData(2, :);
        mriTTL = runData(3, :);
        respiration = respiration .* polarityMultipliers(runIndex);

        ttlRises = find(diff([false, eventTTL > eventThreshold]) == 1);
        if numel(ttlRises) ~= expectedTTLRises
            error('OX_DATA:UnexpectedTTLCount', ...
                '%s run %d has %d TTL rises; expected %d.', ...
                rawFile.name, runIndex, numel(ttlRises), expectedTTLRises);
        end
        sniffSamples = ttlRises(5:2:end);
        if numel(sniffSamples) ~= expectedTrials
            error('OX_DATA:UnexpectedSniffCount', ...
                '%s run %d has %d sniff TTLs; expected %d.', ...
                rawFile.name, runIndex, numel(sniffSamples), expectedTrials);
        end

        mriOnset = find(mriTTL > mriThreshold, 1, 'first');
        if isempty(mriOnset)
            error('OX_DATA:MissingMRIOnset', ...
                '%s run %d has no MRI pulse above %.3g.', ...
                rawFile.name, runIndex, mriThreshold);
        end

        filteredRespiration = lowpass(respiration, lowpassCutoff, samplingRate);
        bmObject = breathmetrics(filteredRespiration, samplingRate, dataType);
        bmObject.estimateAllFeatures(0, 'sliding', 1, 0);
        inhaleSamples = bmObject.inhaleOnsets(:)';
        if isempty(inhaleSamples) || any(~isfinite(inhaleSamples))
            error('OX_DATA:InvalidInhaleOnsets', ...
                'BreathMetrics returned invalid inhale onsets for %s run %d.', ...
                rawFile.name, runIndex);
        end

        nearestSamples = nan(expectedTrials, 1);
        for trialIndex = 1:expectedTrials
            [~, nearestIndex] = min(abs(inhaleSamples - sniffSamples(trialIndex)));
            nearestSamples(trialIndex) = inhaleSamples(nearestIndex);

            % Independently verify the stored result is an absolute-nearest
            % BreathMetrics inhale onset (ties resolve to the first onset).
            minimumDistance = min(abs(inhaleSamples - sniffSamples(trialIndex)));
            if abs(nearestSamples(trialIndex) - sniffSamples(trialIndex)) ~= ...
                    minimumDistance
                error('OX_DATA:NearestMatchValidationFailed', ...
                    'Nearest-onset validation failed for %s run %d trial %d.', ...
                    rawFile.name, runIndex, trialIndex);
            end
        end

        allInhaleOnsets{runIndex} = inhaleSamples;
        sniffTTLSamples(:, runIndex) = sniffSamples(:);
        matchedInhaleSamples(:, runIndex) = nearestSamples;
        matchDeltaSeconds(:, runIndex) = ...
            (nearestSamples - sniffSamples(:)) ./ samplingRate;
        eventOnsets(:, runIndex) = ...
            (nearestSamples - mriOnset) ./ samplingRate;

        if runIndex == 1 || mod(runIndex, 10) == 0 || runIndex == nruns
            fprintf('  completed run %d/%d\n', runIndex, nruns);
        end
    end

    outputData = legacyData;
    outputData.event_onsets = eventOnsets;
    outputData.event_onsets_all = bsxfun(@plus, eventOnsets, legacyData.offsets);
    outputData.event_onsets_vec = outputData.event_onsets_all(:);
    outputData.bm_all_inhale_onsets_samples = allInhaleOnsets;
    outputData.bm_sniff_ttl_samples = sniffTTLSamples;
    outputData.bm_matched_inhale_samples = matchedInhaleSamples;
    outputData.bm_match_delta_seconds = matchDeltaSeconds;
    outputData.bm_processing = makeProcessingMetadata( ...
        samplingRate, lowpassCutoff, dataType, eventThreshold, mriThreshold, ...
        rawFile.name, legacyName, which('breathmetrics'), scriptFile, ...
        polarityMultipliers, polarityFlipRuns);

    validateOutput(outputData, legacyData, nruns, expectedTrials, ...
        allInhaleOnsets, subjectName);

    temporaryPath = [tempname(labchartDir), '.mat'];
    try
        save(temporaryPath, '-struct', 'outputData');
        savedData = load(temporaryPath);
        validateOutput(savedData, legacyData, nruns, expectedTrials, ...
            allInhaleOnsets, subjectName);
        [moveSucceeded, moveMessage] = movefile(temporaryPath, outputPath, 'f');
        if ~moveSucceeded
            error('OX_DATA:OutputMoveFailed', ...
                'Could not move output to %s: %s', outputPath, moveMessage);
        end
    catch processingError
        if isfile(temporaryPath)
            delete(temporaryPath);
        end
        rethrow(processingError);
    end

    absoluteDelta = abs(matchDeltaSeconds(:));
    fprintf(['%s saved: %d matches, median |delta| %.3f s, ', ...
        'maximum |delta| %.3f s\n'], outputName, numel(absoluteDelta), ...
        median(absoluteDelta), max(absoluteDelta));
end

fprintf('\nBreathMetrics onset reprocessing complete for %d subjects.\n', ...
    numel(rawFiles));

function validateLegacyTemplate(legacyData, nruns, legacyName, expectedTrials)
requiredFields = {'event_onsets', 'event_onsets_all', 'event_onsets_vec', ...
    'cue_onsets', 'cue_onsets_all', 'cue_onsets_vec', 'nframes', ...
    'offsets', 'sessi', 'sessf'};
missingFields = requiredFields(~isfield(legacyData, requiredFields));
if ~isempty(missingFields)
    error('OX_DATA:InvalidLegacyTemplate', ...
        '%s is missing fields: %s', legacyName, strjoin(missingFields, ', '));
end
if ~isequal(size(legacyData.event_onsets), [expectedTrials, nruns]) || ...
        ~isequal(size(legacyData.event_onsets_all), [expectedTrials, nruns]) || ...
        numel(legacyData.event_onsets_vec) ~= expectedTrials * nruns || ...
        numel(legacyData.offsets) ~= nruns
    error('OX_DATA:LegacyDimensionMismatch', ...
        '%s onset or offset dimensions do not match %d trials by %d runs.', ...
        legacyName, expectedTrials, nruns);
end
if ~isequal(size(legacyData.cue_onsets), [expectedTrials, nruns]) || ...
        ~isequal(size(legacyData.cue_onsets_all), [expectedTrials, nruns]) || ...
        numel(legacyData.cue_onsets_vec) ~= expectedTrials * nruns
    error('OX_DATA:LegacyCueDimensionMismatch', ...
        '%s cue dimensions do not match %d trials by %d runs.', ...
        legacyName, expectedTrials, nruns);
end
end

function metadata = makeProcessingMetadata(samplingRate, lowpassCutoff, ...
        dataType, eventThreshold, mriThreshold, rawName, legacyName, ...
        breathmetricsClassFile, scriptFile, polarityMultipliers, ...
        polarityFlipRuns)
metadata = struct();
metadata.created_at = char(datetime('now', 'TimeZone', 'local', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
metadata.sampling_rate_hz = samplingRate;
metadata.lowpass_cutoff_hz = lowpassCutoff;
metadata.lowpass_function = 'lowpass';
metadata.data_type = dataType;
metadata.event_ttl_threshold = eventThreshold;
metadata.mri_ttl_threshold = mriThreshold;
metadata.estimate_all_features = struct( ...
    'z_score', 0, ...
    'baseline_correction_method', 'sliding', ...
    'simplify', 1, ...
    'verbose', 0);
metadata.source_raw_file = rawName;
metadata.source_legacy_file = legacyName;
metadata.breathmetrics_class_file = breathmetricsClassFile;
metadata.processing_script = scriptFile;
metadata.matlab_version = version;
metadata.respiration_polarity_multiplier_by_run = polarityMultipliers;
metadata.polarity_flipped_runs = polarityFlipRuns;
if isempty(polarityFlipRuns)
    metadata.polarity_flip_reason = '';
else
    metadata.polarity_flip_reason = ...
        'Reversed acquisition polarity; confirmed from sniff-locked waveform';
end
end

function validateOutput(outputData, legacyData, nruns, expectedTrials, ...
        allInhaleOnsets, subjectName)
replacedFields = {'event_onsets', 'event_onsets_all', 'event_onsets_vec'};
legacyFields = fieldnames(legacyData);
preservedFields = setdiff(legacyFields, replacedFields);
for fieldIndex = 1:numel(preservedFields)
    fieldName = preservedFields{fieldIndex};
    if ~isfield(outputData, fieldName) || ...
            ~isequaln(outputData.(fieldName), legacyData.(fieldName))
        error('OX_DATA:LegacyFieldChanged', ...
            'Legacy field %s changed unexpectedly for %s.', ...
            fieldName, subjectName);
    end
end

if ~isequal(size(outputData.event_onsets), [expectedTrials, nruns]) || ...
        any(~isfinite(outputData.event_onsets(:)))
    error('OX_DATA:InvalidOutputOnsets', ...
        '%s does not contain %d finite event onsets per run.', ...
        subjectName, expectedTrials);
end
expectedAll = bsxfun(@plus, outputData.event_onsets, outputData.offsets);
if ~isequaln(outputData.event_onsets_all, expectedAll) || ...
        ~isequaln(outputData.event_onsets_vec, expectedAll(:))
    error('OX_DATA:InvalidGlobalOnsets', ...
        'Global or vectorized event onsets are inconsistent for %s.', ...
        subjectName);
end

qcFields = {'bm_all_inhale_onsets_samples', 'bm_sniff_ttl_samples', ...
    'bm_matched_inhale_samples', 'bm_match_delta_seconds', 'bm_processing'};
missingQCFields = qcFields(~isfield(outputData, qcFields));
if ~isempty(missingQCFields)
    error('OX_DATA:MissingQCFields', ...
        '%s is missing QC fields: %s', ...
        subjectName, strjoin(missingQCFields, ', '));
end
if ~isequal(outputData.bm_all_inhale_onsets_samples, allInhaleOnsets) || ...
        ~isequal(size(outputData.bm_sniff_ttl_samples), ...
            [expectedTrials, nruns]) || ...
        ~isequal(size(outputData.bm_matched_inhale_samples), ...
            [expectedTrials, nruns]) || ...
        ~isequal(size(outputData.bm_match_delta_seconds), ...
            [expectedTrials, nruns]) || ...
        any(~isfinite(outputData.bm_matched_inhale_samples(:))) || ...
        any(~isfinite(outputData.bm_match_delta_seconds(:)))
    error('OX_DATA:InvalidQCFields', ...
        'BreathMetrics QC fields are invalid for %s.', subjectName);
end
end
