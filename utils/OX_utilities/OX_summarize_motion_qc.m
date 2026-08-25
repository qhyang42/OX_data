function summary = OX_summarize_motion_qc(varargin)
%OX_SUMMARIZE_MOTION_QC Create dataset-level tables and compact QC figures.
%   SUMMARY = OX_SUMMARIZE_MOTION_QC() loads saved motion QC results for
%   subjects 2:6. Existing summary files are never overwritten.

p = inputParser;
p.addParameter('Subjects', 2:6, @(x) isnumeric(x) && isvector(x));
p.addParameter('ProjectRoot', localProjectRoot(), @(x) ischar(x) || isstring(x));
p.addParameter('OutputDir', '', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

projectRoot = char(p.Results.ProjectRoot);
outputDir = char(p.Results.OutputDir);
if isempty(outputDir)
    outputDir = fullfile(projectRoot, 'motion_qc');
end
subjects = p.Results.Subjects(:)';
if ~isfolder(outputDir)
    error('OX:MotionQC:MissingOutputDirectory', ...
        'Motion QC output directory does not exist: %s', outputDir);
end

participantCsv = fullfile(outputDir, 'motion_qc_participant_summary.csv');
runCsv = fullfile(outputDir, 'motion_qc_run_summary.csv');
candidateCsv = fullfile(outputDir, 'motion_qc_candidate_summary.csv');
wideCsv = fullfile(outputDir, 'motion_qc_candidate_percent_wide.csv');
comparisonPng = fullfile(outputDir, 'motion_qc_candidate_comparison.png');
extremesPng = fullfile(outputDir, 'motion_qc_run_extremes.png');
planned = {participantCsv, runCsv, candidateCsv, wideCsv, comparisonPng, extremesPng};
for index = 1:numel(subjects)
    planned{end + 1} = fullfile(outputDir, sprintf('subj_%d_distributions.png', subjects(index))); %#ok<AGROW>
    planned{end + 1} = fullfile(outputDir, sprintf('subj_%d_runwise.png', subjects(index))); %#ok<AGROW>
end
existing = planned(cellfun(@isfile, planned));
if ~isempty(existing)
    error('OX:MotionQC:WouldOverwrite', ...
        'Refusing to overwrite existing summary output: %s', strjoin(existing, ', '));
end

loaded = cell(numel(subjects), 1);
for subjectIndex = 1:numel(subjects)
    file = fullfile(outputDir, sprintf('subj_%d_motion_qc.mat', subjects(subjectIndex)));
    if ~isfile(file)
        error('OX:MotionQC:MissingSubjectResult', 'Missing subject result: %s', file);
    end
    value = load(file, 'results');
    loaded{subjectIndex} = value.results;
end
validateCompatibleRules(loaded);
rules = loaded{1}.candidate_rules;
nRules = numel(rules);

participantRows = cell(numel(subjects), 1);
runTables = cell(numel(subjects), 1);
candidateTables = cell(numel(subjects), 1);
overallPercent = zeros(numel(subjects), nRules);
maxRunPercent = zeros(numel(subjects), nRules);

for subjectIndex = 1:numel(subjects)
    result = loaded{subjectIndex};
    fd = vertcat(result.fd_mm{:});
    raw = vertcat(result.raw_dvars{:});
    z = vertcat(result.robust_dvars_z{:});
    nonInitialCells = cellfun(@(x) [false; true(numel(x)-1, 1)], ...
        result.fd_mm, 'UniformOutput', false);
    nonInitial = vertcat(nonInitialCells{:});
    correlation = correlationValue(fd(nonInitial), raw(nonInitial));
    correlationZ = correlationValue(fd(nonInitial), z(nonInitial));
    fdStats = distributionStats(fd);
    rawStats = distributionStats(raw);
    zStats = distributionStats(z);
    participantRows{subjectIndex} = table(string(result.subject), numel(result.run_ids), ...
        sum(result.n_volumes_per_run), correlation, correlationZ, ...
        fdStats(1), fdStats(2), fdStats(3), fdStats(4), fdStats(5), ...
        rawStats(1), rawStats(2), rawStats(3), rawStats(4), rawStats(5), ...
        zStats(1), zStats(2), zStats(3), zStats(4), zStats(5), ...
        'VariableNames', {'subject', 'n_runs', 'n_volumes', ...
        'fd_raw_dvars_correlation', 'fd_robust_dvars_z_correlation', ...
        'fd_median', 'fd_p75', 'fd_p95', 'fd_p99', 'fd_max', ...
        'raw_dvars_median', 'raw_dvars_p75', 'raw_dvars_p95', ...
        'raw_dvars_p99', 'raw_dvars_max', 'robust_dvars_z_median', ...
        'robust_dvars_z_p75', 'robust_dvars_z_p95', ...
        'robust_dvars_z_p99', 'robust_dvars_z_max'});

    runTables{subjectIndex} = makeRunTable(result);
    thisCandidate = result.subject_rule_summary;
    thisCandidate.subject = repmat(string(result.subject), height(thisCandidate), 1);
    thisCandidate = movevars(thisCandidate, 'subject', 'Before', 1);
    candidateTables{subjectIndex} = thisCandidate;
    overallPercent(subjectIndex, :) = thisCandidate.percent_censored';
    maxRunPercent(subjectIndex, :) = thisCandidate.max_run_percent';

    makeDistributionFigure(result, fullfile(outputDir, ...
        sprintf('%s_distributions.png', result.subject)));
    makeRunwiseFigure(result, fullfile(outputDir, ...
        sprintf('%s_runwise.png', result.subject)));
end

participantSummary = vertcat(participantRows{:});
runSummary = vertcat(runTables{:});
candidateSummary = vertcat(candidateTables{:});
wideNames = matlab.lang.makeValidName({rules.id});
candidatePercentWide = array2table(overallPercent, 'VariableNames', wideNames);
subjectNames = string(cellfun(@(x) x.subject, loaded, 'UniformOutput', false));
candidatePercentWide = addvars(candidatePercentWide, subjectNames, ...
    'Before', 1, 'NewVariableNames', 'subject');

writetable(participantSummary, participantCsv);
writetable(runSummary, runCsv);
writetable(candidateSummary, candidateCsv);
writetable(candidatePercentWide, wideCsv);
makeComparisonFigure(subjects, rules, overallPercent, maxRunPercent, comparisonPng);
makeExtremesFigure(loaded, extremesPng);

summary = struct('participant', participantSummary, 'run', runSummary, ...
    'candidate', candidateSummary, 'candidate_percent_wide', candidatePercentWide);
summary.files = planned;
fprintf('Saved dataset motion QC summary in %s\n', outputDir);
end

function tableOut = makeRunTable(result)
nRuns = numel(result.run_ids);
nRules = numel(result.candidate_rules);
subject = repmat(string(result.subject), nRuns, 1);
ordinal = result.run_order(:, 1);
session = result.run_order(:, 2);
run = result.run_order(:, 3);
runId = string(result.run_ids);
nVolumes = result.n_volumes_per_run;
fdMedian = zeros(nRuns, 1); fdP95 = zeros(nRuns, 1); fdMax = zeros(nRuns, 1);
rawMedian = zeros(nRuns, 1); rawP95 = zeros(nRuns, 1); rawMax = zeros(nRuns, 1);
zP95 = zeros(nRuns, 1); zMax = zeros(nRuns, 1); correlation = zeros(nRuns, 1);
for index = 1:nRuns
    fdStats = distributionStats(result.fd_mm{index});
    rawStats = distributionStats(result.raw_dvars{index});
    zStats = distributionStats(result.robust_dvars_z{index});
    fdMedian(index) = fdStats(1); fdP95(index) = fdStats(3); fdMax(index) = fdStats(5);
    rawMedian(index) = rawStats(1); rawP95(index) = rawStats(3); rawMax(index) = rawStats(5);
    zP95(index) = zStats(3); zMax(index) = zStats(5);
    correlation(index) = correlationValue(result.fd_mm{index}(2:end), result.raw_dvars{index}(2:end));
end
tableOut = table(subject, ordinal, session, run, runId, nVolumes, ...
    fdMedian, fdP95, fdMax, rawMedian, rawP95, rawMax, zP95, zMax, correlation, ...
    'VariableNames', {'subject', 'ordinal', 'session', 'run', 'run_id', ...
    'n_volumes', 'fd_median', 'fd_p95', 'fd_max', 'raw_dvars_median', ...
    'raw_dvars_p95', 'raw_dvars_max', 'robust_dvars_z_p95', ...
    'robust_dvars_z_max', 'fd_raw_dvars_correlation'});
for ruleIndex = 1:nRules
    variable = matlab.lang.makeValidName(['pct_' result.candidate_rules(ruleIndex).id]);
    tableOut.(variable) = result.censor_percent_by_run(:, ruleIndex);
end
end

function stats = distributionStats(values)
values = sort(values(:));
stats = [percentile(values, 50), percentile(values, 75), ...
    percentile(values, 95), percentile(values, 99), values(end)];
end

function value = percentile(sortedValues, percent)
n = numel(sortedValues);
position = 1 + (n - 1) * percent / 100;
lower = floor(position);
upper = ceil(position);
if lower == upper
    value = sortedValues(lower);
else
    value = sortedValues(lower) + (position - lower) * ...
        (sortedValues(upper) - sortedValues(lower));
end
end

function value = correlationValue(x, y)
c = corrcoef(x(:), y(:));
if numel(c) < 4 || ~isfinite(c(1, 2))
    value = NaN;
else
    value = c(1, 2);
end
end

function makeDistributionFigure(result, outputFile)
fd = vertcat(result.fd_mm{:});
raw = vertcat(result.raw_dvars{:});
z = vertcat(result.robust_dvars_z{:});
figureHandle = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1200 760]);
cleanup = onCleanup(@() close(figureHandle));
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; histogram(fd, 60); xlabel('FD (mm)'); ylabel('Volumes'); title('Power FD');
set(gca, 'YScale', 'log'); xlim([0 max(percentile(sort(fd), 99.5) * 1.1, eps)]);
nexttile; histogram(raw, 60); xlabel('Raw DVARS'); ylabel('Volumes'); title('Raw DVARS');
set(gca, 'YScale', 'log'); xlim([0 max(percentile(sort(raw), 99.5) * 1.1, eps)]);
nexttile; histogram(z, 60); xlabel('Robust DVARS z'); ylabel('Volumes'); title('Run-wise robust standardized DVARS');
set(gca, 'YScale', 'log'); xlim([min(-2, min(z)) max(8, percentile(sort(z), 99.5) * 1.1)]);
nexttile; scatter(fd, z, 7, '.', 'MarkerEdgeAlpha', 0.15); xlabel('FD (mm)'); ylabel('Robust DVARS z');
title(sprintf('FD vs DVARS z (r = %.2f)', correlationValue(fd, z))); grid on;
sgtitle(sprintf('%s motion/artifact distributions', strrep(result.subject, '_', '\_')));
exportgraphics(figureHandle, outputFile, 'Resolution', 180);
end

function makeRunwiseFigure(result, outputFile)
nRuns = numel(result.run_ids);
maxVolumes = max(result.n_volumes_per_run);
fdMatrix = nan(nRuns, maxVolumes);
zMatrix = nan(nRuns, maxVolumes);
exceedanceCount = nan(nRuns, maxVolumes);
for runIndex = 1:nRuns
    n = result.n_volumes_per_run(runIndex);
    fdMatrix(runIndex, 1:n) = result.fd_mm{runIndex};
    zMatrix(runIndex, 1:n) = result.robust_dvars_z{runIndex};
    count = zeros(n, 1);
    for ruleIndex = 1:numel(result.candidate_rules)
        count = count + result.censor_vectors.(result.candidate_rules(ruleIndex).id){runIndex};
    end
    exceedanceCount(runIndex, 1:n) = count;
end
figureHandle = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1450 900]);
cleanup = onCleanup(@() close(figureHandle));
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; imagesc(fdMatrix); colorbar; xlabel('Volume within run'); ylabel('Run ordinal'); title('FD (mm)');
clim([0 max(percentile(sort(fdMatrix(isfinite(fdMatrix))), 99), eps)]);
nexttile; imagesc(zMatrix); colorbar; xlabel('Volume within run'); ylabel('Run ordinal'); title('Robust DVARS z');
clim([0 max(percentile(sort(zMatrix(isfinite(zMatrix))), 99), eps)]);
nexttile; imagesc(exceedanceCount); colorbar; xlabel('Volume within run'); ylabel('Run ordinal');
title('Number of candidate rules exceeded'); clim([0 numel(result.candidate_rules)]);
nexttile; imagesc(result.censor_percent_by_run); colorbar; xlabel('Candidate rule'); ylabel('Run ordinal');
title('Censored volumes per run (%)'); xticks(1:numel(result.candidate_rules));
xticklabels({result.candidate_rules.id}); xtickangle(40);
sgtitle(sprintf('%s run-wise motion/artifact overview', strrep(result.subject, '_', '\_')));
exportgraphics(figureHandle, outputFile, 'Resolution', 180);
end

function makeComparisonFigure(subjects, rules, overallPercent, maxRunPercent, outputFile)
figureHandle = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1500 650]);
cleanup = onCleanup(@() close(figureHandle));
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; imagesc(overallPercent); colorbar; title('Overall censored volumes (%)');
ylabel('Participant'); xlabel('Candidate rule'); yticks(1:numel(subjects));
yticklabels(compose('subj_%d', subjects)); xticks(1:numel(rules));
xticklabels({rules.id}); xtickangle(40);
nexttile; imagesc(maxRunPercent); colorbar; title('Maximum censored in any run (%)');
ylabel('Participant'); xlabel('Candidate rule'); yticks(1:numel(subjects));
yticklabels(compose('subj_%d', subjects)); xticks(1:numel(rules));
xticklabels({rules.id}); xtickangle(40);
sgtitle('Diagnostic candidate censor-rule comparison (no rule selected)');
exportgraphics(figureHandle, outputFile, 'Resolution', 180);
end

function makeExtremesFigure(results, outputFile)
figureHandle = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1100 750]);
cleanup = onCleanup(@() close(figureHandle));
hold on;
colors = lines(numel(results));
legendHandles = gobjects(numel(results), 1);
allRows = [];
for subjectIndex = 1:numel(results)
    result = results{subjectIndex};
    maxFd = cellfun(@max, result.fd_mm);
    maxZ = cellfun(@max, result.robust_dvars_z);
    legendHandles(subjectIndex) = scatter(maxFd, maxZ, 36, colors(subjectIndex, :), 'filled', ...
        'MarkerFaceAlpha', 0.65);
    score = tiedRank(maxFd) + tiedRank(maxZ);
    allRows = [allRows; [repmat(subjectIndex, numel(maxFd), 1), ...
        (1:numel(maxFd))', maxFd, maxZ, score]]; %#ok<AGROW>
end
[~, order] = sort(allRows(:, 5), 'descend');
for index = order(1:min(12, size(allRows, 1)))'
    subjectIndex = allRows(index, 1);
    runIndex = allRows(index, 2);
    label = sprintf(' %s:%s', results{subjectIndex}.subject, results{subjectIndex}.run_ids{runIndex});
    text(allRows(index, 3), allRows(index, 4), label, 'FontSize', 7);
end
xlabel('Maximum FD in run (mm)'); ylabel('Maximum robust DVARS z in run');
title('Run extremes (labels show joint FD/DVARS outliers)'); grid on;
legend(legendHandles, cellfun(@(x) x.subject, results, 'UniformOutput', false), ...
    'Location', 'best');
exportgraphics(figureHandle, outputFile, 'Resolution', 180);
end

function ranks = tiedRank(values)
[sorted, order] = sort(values(:));
ranksSorted = zeros(size(sorted));
first = 1;
while first <= numel(sorted)
    offset = find(sorted(first:end) ~= sorted(first), 1, 'first');
    if isempty(offset)
        last = numel(sorted) + 1;
    else
        last = first + offset - 1;
    end
    last = last - 1;
    ranksSorted(first:last) = mean(first:last);
    first = last + 1;
end
ranks = zeros(size(values(:)));
ranks(order) = ranksSorted;
end

function validateCompatibleRules(results)
reference = {results{1}.candidate_rules.id};
for index = 2:numel(results)
    if ~isequal(reference, {results{index}.candidate_rules.id})
        error('OX:MotionQC:IncompatibleRules', ...
            'Candidate rule definitions differ across subject result files.');
    end
end
end

function root = localProjectRoot()
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
