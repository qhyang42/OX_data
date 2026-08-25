function summary = OX_summarize_glmsingle_confounds(varargin)
%OX_SUMMARIZE_GLMSINGLE_CONFOUNDS Create compact QC for completed derivatives.
%   SUMMARY = OX_SUMMARIZE_GLMSINGLE_CONFOUNDS('Subjects',2:6) requires
%   all 80 standalone files for every requested subject and refuses partial
%   datasets. Existing QC outputs are never overwritten.

p = inputParser;
p.addParameter('Subjects', 2:6, @(x) isnumeric(x) && isvector(x));
p.addParameter('ProjectRoot', localProjectRoot(), @(x) ischar(x) || isstring(x));
p.addParameter('OutputDir', '', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
subjects = p.Results.Subjects(:)';
projectRoot = char(p.Results.ProjectRoot);
labchartDir = fullfile(projectRoot, 'labchart', 'extracted_events');
outputDir = char(p.Results.OutputDir);
if isempty(outputDir)
    outputDir = fullfile(projectRoot, 'glmsingle_confounds_qc');
end
if isfolder(outputDir)
    listing = dir(outputDir);
    listing = listing(~ismember({listing.name}, {'.', '..'}));
    if ~isempty(listing)
        error('OX:GLMsingleConfoundsQC:WouldOverwrite', ...
            'QC output directory is not empty: %s', outputDir);
    end
end

nSubjects = numel(subjects);
physiologyNames = {'airflow_z', 'airflow_sq_z', 'resp_volume_z'};
allPhys = cell(nSubjects, 1);
meanCorr27 = cell(nSubjects, 1);
runRows = cell(0, 16);
taskRows = cell(0, 17);
exampleRuns = [1, 40, 80];
exampleData = cell(nSubjects, numel(exampleRuns));

for subjectIndex = 1:nSubjects
    subjidx = subjects(subjectIndex);
    subject = sprintf('subj_%d', subjidx);
    runs = OX_discover_functional_runs(subjidx, 'ProjectRoot', projectRoot, ...
        'ExpectedRuns', 80);
    bm = load(fullfile(labchartDir, ...
        sprintf('subj%d_events_bm.mat', subjidx)), 'event_onsets', ...
        'bm_processing');
    [~, categories] = OX_get_odor(subject);
    if numel(categories) ~= 800
        error('OX:GLMsingleConfoundsQC:TaskLabelCount', ...
            '%s must have exactly 800 context labels.', subject);
    end
    subjectPhys = cell(80, 1);
    subjectTask = cell(80, 1);
    corrSum = zeros(27);
    for runIndex = 1:80
        filename = fullfile(projectRoot, 'MRI', subject, 'nifti', ...
            'glmsingle_confounds', sprintf( ...
            '%s_session%02d_run%02d_confounds.mat', subject, ...
            runs(runIndex).session, runs(runIndex).run));
        if ~isfile(filename)
            error('OX:GLMsingleConfoundsQC:MissingRun', ...
                'Missing completed confound derivative: %s', filename);
        end
        c = load(filename);
        validateSavedRun(c, runs(runIndex), runIndex);
        physiology = [c.airflow_TR_z, c.airflow_sq_TR_z, c.resp_volume_TR_z];
        subjectPhys{runIndex} = physiology;
        corrSum = corrSum + corr(c.confounds(:, 1:27));
        nSpikes = nnz(c.bad_volume);
        driftR = corr(c.resp_volume_TR_raw, (1:c.nframes)');
        saturationRatio = range(c.resp_volume_TR_raw) / ...
            max(iqr(c.resp_volume_TR_raw), eps);
        alignmentMethod = string(c.processing_metadata.mri_frame_alignment.method);
        if isfield(c, 'respiratory_incomplete_volume')
            finalRespiratoryReplacement = any(c.respiratory_incomplete_volume);
        else
            finalRespiratoryReplacement = false;
        end
        runRows(end + 1, :) = {subjidx, runIndex, c.session_id, c.run_id, ...
            string(c.run_label), c.nframes, nSpikes, 100*nSpikes/c.nframes, ...
            size(c.confounds, 2), c.processing_metadata.design_rank, ...
            c.processing_metadata.design_condition_number, driftR, ...
            saturationRatio, alignmentMethod, abs(driftR) > 0.8, ...
            finalRespiratoryReplacement}; %#ok<AGROW>

        trialIndex = (runIndex - 1) * 10 + (1:10);
        onsetFrames = ceil(bm.event_onsets(:, runIndex) ./ 0.76);
        if any(onsetFrames < 1 | onsetFrames > c.nframes)
            error('OX:GLMsingleConfoundsQC:OnsetOutsideRun', ...
                '%s has an odor onset outside run %d.', subject, runIndex);
        end
        task = zeros(c.nframes, 5);
        task(onsetFrames, 1) = 1;
        contextNames = {'PERSON', 'FOOD', 'LOCATION', 'CONTROL'};
        runCategories = categories(trialIndex);
        for contextIndex = 1:4
            thisFrame = onsetFrames(strcmp(runCategories, contextNames{contextIndex}));
            task(thisFrame, contextIndex + 1) = 1;
        end
        subjectTask{runIndex} = task;

        exampleIndex = find(exampleRuns == runIndex, 1);
        if ~isempty(exampleIndex)
            exampleData{subjectIndex, exampleIndex} = struct( ...
                'run', c, 'raw_file', fullfile(labchartDir, ...
                bm.bm_processing.source_raw_file));
        end
    end
    allPhys{subjectIndex} = vertcat(subjectPhys{:});
    meanCorr27{subjectIndex} = corrSum ./ 80;
    pcat = vertcat(subjectPhys{:});
    tcat = vertcat(subjectTask{:});
    taskCorrelation = corr(pcat, tcat);
    for physiologyIndex = 1:3
        taskRows(end + 1, :) = {subjidx, physiologyIndex, ...
            string(physiologyNames{physiologyIndex}), ...
            taskCorrelation(physiologyIndex, 1), ...
            taskCorrelation(physiologyIndex, 2), ...
            taskCorrelation(physiologyIndex, 3), ...
            taskCorrelation(physiologyIndex, 4), ...
            taskCorrelation(physiologyIndex, 5), ...
            max(abs(taskCorrelation(physiologyIndex, :))), ...
            mean(pcat(:, physiologyIndex)), std(pcat(:, physiologyIndex)), ...
            min(pcat(:, physiologyIndex)), max(pcat(:, physiologyIndex)), ...
            prctile(pcat(:, physiologyIndex), 1), ...
            prctile(pcat(:, physiologyIndex), 50), ...
            prctile(pcat(:, physiologyIndex), 99), size(pcat, 1)}; %#ok<AGROW>
    end
end

runSummary = cell2table(runRows, 'VariableNames', {'subject', 'run_ordinal', ...
    'session_id', 'run_id', 'run_label', 'nframes', 'n_spikes', ...
    'percent_spikes', 'n_confound_columns', 'design_rank', ...
    'design_condition_number', 'resp_volume_time_correlation', ...
    'resp_volume_range_over_iqr', 'alignment_method', ...
    'integrated_airflow_drift_flag', 'final_volume_respiration_replaced'});
taskLocking = cell2table(taskRows, 'VariableNames', {'subject', ...
    'physiology_index', 'physiology_name', 'r_odor_all', 'r_person', ...
    'r_food', 'r_location', 'r_control', 'max_abs_task_r', 'mean', 'std', ...
    'minimum', 'maximum', 'p01', 'median', 'p99', 'n_frames'});
subjectSummary = makeSubjectSummary(runSummary, subjects);

mkdir(outputDir);
writetable(runSummary, fullfile(outputDir, 'run_summary.csv'));
writetable(subjectSummary, fullfile(outputDir, 'subject_summary.csv'));
writetable(taskLocking, fullfile(outputDir, 'respiratory_task_locking.csv'));
plotDistributions(allPhys, subjects, fullfile(outputDir, ...
    'physiology_distributions.png'));
plotSpikeCounts(runSummary, subjects, fullfile(outputDir, ...
    'runwise_spike_counts.png'));
plotContinuousCorrelations(meanCorr27, subjects, fullfile(outputDir, ...
    'mean_runwise_continuous_correlation.png'));
plotAlignmentExamples(exampleData, subjects, exampleRuns, fullfile(outputDir, ...
    'airflow_alignment_examples.png'));

summary = struct('run_summary', runSummary, 'subject_summary', subjectSummary, ...
    'task_locking', taskLocking, 'output_directory', outputDir);
save(fullfile(outputDir, 'qc_summary.mat'), 'summary', '-v7.3');
end

function validateSavedRun(c, runInfo, ordinal)
required = {'confounds','confound_names','motion24','airflow_TR_raw', ...
    'airflow_sq_TR_raw','resp_volume_TR_raw','airflow_TR_z', ...
    'airflow_sq_TR_z','resp_volume_TR_z','FD','robust_DVARS_z', ...
    'bad_volume','bad_volume_reason','run_id','session_id','nframes', ...
    'processing_metadata'};
if ~all(isfield(c, required)) || c.run_id ~= runInfo.run || ...
        c.session_id ~= runInfo.session || c.run_ordinal ~= ordinal || ...
        size(c.confounds, 1) ~= c.nframes || any(~isfinite(c.confounds(:)))
    error('OX:GLMsingleConfoundsQC:InvalidSavedRun', ...
        'Saved confound file failed validation for %s.', runInfo.id);
end
end

function T = makeSubjectSummary(runSummary, subjects)
rows = cell(numel(subjects), 9);
for i = 1:numel(subjects)
    z = runSummary(runSummary.subject == subjects(i), :);
    rows(i, :) = {subjects(i), sum(z.nframes), sum(z.n_spikes), ...
        100*sum(z.n_spikes)/sum(z.nframes), min(z.n_confound_columns), ...
        max(z.n_confound_columns), nnz(z.integrated_airflow_drift_flag), ...
        max(z.resp_volume_range_over_iqr), ...
        nnz(z.final_volume_respiration_replaced)};
end
T = cell2table(rows, 'VariableNames', {'subject','n_frames','n_bad_volumes', ...
    'percent_bad_volumes','min_nuisance_columns','max_nuisance_columns', ...
    'runs_with_integrated_airflow_drift_flag','max_resp_volume_range_over_iqr', ...
    'runs_with_final_volume_respiration_replaced'});
end

function plotDistributions(allPhys, subjects, filename)
f = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1200 900]);
tiledlayout(numel(subjects), 3, 'TileSpacing', 'compact');
names = {'airflow z', 'airflow squared z', 'respiratory volume z'};
for s = 1:numel(subjects)
    for j = 1:3
        nexttile; histogram(allPhys{s}(:,j), 60, 'Normalization', 'pdf');
        title(sprintf('S%d %s', subjects(s), names{j})); xlim([-6 6]);
    end
end
exportgraphics(f, filename, 'Resolution', 160); close(f);
end

function plotSpikeCounts(T, subjects, filename)
f = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1200 850]);
tiledlayout(numel(subjects), 1, 'TileSpacing', 'compact');
for s = 1:numel(subjects)
    z = T(T.subject == subjects(s), :); nexttile;
    bar(z.run_ordinal, z.percent_spikes, 1); ylabel(sprintf('S%d %%', subjects(s)));
    xlim([0 81]);
end
xlabel('run ordinal');
exportgraphics(f, filename, 'Resolution', 160); close(f);
end

function plotContinuousCorrelations(correlations, subjects, filename)
f = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1200 850]);
tiledlayout(2, ceil(numel(subjects)/2), 'TileSpacing', 'compact');
for s = 1:numel(subjects)
    nexttile; imagesc(correlations{s}, [-1 1]); axis image;
    title(sprintf('S%d mean run-wise r (27 continuous)', subjects(s)));
end
colorbar;
exportgraphics(f, filename, 'Resolution', 160); close(f);
end

function plotAlignmentExamples(examples, subjects, runOrdinals, filename)
f = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1400 1000]);
tiledlayout(numel(subjects), numel(runOrdinals), 'TileSpacing', 'compact');
for s = 1:numel(subjects)
    raw = load(examples{s,1}.raw_file, 'eventdata');
    for j = 1:numel(runOrdinals)
        e = examples{s,j}; c = e.run; runIndex = runOrdinals(j);
        x = double(raw.eventdata{runIndex}(1,:)) .* ...
            c.processing_metadata.respiration_polarity_multiplier;
        edges = c.processing_metadata.mri_frame_alignment.frame_edge_samples;
        last = min(edges(1) + 60*1000 - 1, edges(end)-1);
        idx = edges(1):10:last; rawShown = zscore(x(idx));
        frameCenters = ((edges(1:end-1)+edges(2:end)-1)/2 - edges(1))/1000;
        keepFrame = frameCenters <= 60;
        nexttile; plot((idx-edges(1))/1000, rawShown, 'Color', [.65 .65 .65]);
        hold on; plot(frameCenters(keepFrame), c.airflow_TR_z(keepFrame), ...
            '.-', 'LineWidth', 1); yline(0, ':');
        title(sprintf('S%d run %d', subjects(s), runIndex));
        if s == numel(subjects); xlabel('seconds from BM MRI onset'); end
    end
    clear raw
end
exportgraphics(f, filename, 'Resolution', 160); close(f);
end

function root = localProjectRoot()
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
