function summary = OX_supplement_motion_qc_fd(varargin)
%OX_SUPPLEMENT_MOTION_QC_FD Evaluate extra FD-only censor thresholds.
%   SUMMARY = OX_SUPPLEMENT_MOTION_QC_FD() evaluates strict FD > 0.4 mm
%   and FD > 0.45 mm rules from existing run-wise motion QC results. It
%   does not recompute FD/DVARS or overwrite any existing output.
%
%   Name-value options:
%     Subjects     Subject indices (default 2:6)
%     Thresholds   Additional strict FD thresholds in mm (default [.4 .45])
%     ProjectRoot  OX_DATA root (inferred from this file by default)
%     OutputDir    Existing motion_qc directory (default under project root)

p = inputParser;
p.addParameter('Subjects', 2:6, @(x) isnumeric(x) && isvector(x));
p.addParameter('Thresholds', [0.4 0.45], ...
    @(x) isnumeric(x) && isvector(x) && all(isfinite(x)) && all(x > 0));
p.addParameter('ProjectRoot', localProjectRoot(), @(x) ischar(x) || isstring(x));
p.addParameter('OutputDir', '', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

subjects = p.Results.Subjects(:)';
thresholds = p.Results.Thresholds(:)';
if numel(unique(thresholds)) ~= numel(thresholds)
    error('OX:MotionQC:DuplicateThreshold', 'FD thresholds must be unique.');
end
thresholds = sort(thresholds);
projectRoot = char(p.Results.ProjectRoot);
outputDir = char(p.Results.OutputDir);
if isempty(outputDir)
    outputDir = fullfile(projectRoot, 'motion_qc');
end
if ~isfolder(outputDir)
    error('OX:MotionQC:MissingOutputDirectory', ...
        'Motion QC output directory does not exist: %s', outputDir);
end

thresholdTokens = arrayfun(@numberToken, thresholds, 'UniformOutput', false);
tag = ['fd_' strjoin(thresholdTokens, '_')];
longCsv = fullfile(outputDir, sprintf('motion_qc_%s_summary.csv', tag));
wideCsv = fullfile(outputDir, sprintf('motion_qc_%s_percent_wide.csv', tag));
figureFile = fullfile(outputDir, sprintf('motion_qc_%s_comparison.png', tag));
subjectFiles = arrayfun(@(s) fullfile(outputDir, ...
    sprintf('subj_%d_motion_qc_%s_supplement.mat', s, tag)), ...
    subjects, 'UniformOutput', false);
planned = [{longCsv, wideCsv, figureFile}, subjectFiles];
existing = planned(cellfun(@isfile, planned));
if ~isempty(existing)
    error('OX:MotionQC:WouldOverwrite', ...
        'Refusing to overwrite existing supplemental output: %s', ...
        strjoin(existing, ', '));
end

nSubjects = numel(subjects);
nThresholds = numel(thresholds);
overallPercent = zeros(nSubjects, nThresholds);
maxRunPercent = zeros(nSubjects, nThresholds);
longTables = cell(nSubjects, 1);

for subjectIndex = 1:nSubjects
    subject = sprintf('subj_%d', subjects(subjectIndex));
    sourceFile = fullfile(outputDir, sprintf('%s_motion_qc.mat', subject));
    if ~isfile(sourceFile)
        error('OX:MotionQC:MissingSubjectResult', ...
            'Missing source motion QC result: %s', sourceFile);
    end
    loaded = load(sourceFile, 'results');
    result = loaded.results;
    if numel(result.run_ids) ~= 80
        error('OX:MotionQC:UnexpectedRunCount', ...
            '%s source result contains %d rather than 80 runs.', ...
            subject, numel(result.run_ids));
    end

    supplement = struct();
    supplement.subject = subject;
    supplement.subject_index = subjects(subjectIndex);
    supplement.source_result_file = sourceFile;
    supplement.run_ids = result.run_ids;
    supplement.n_volumes_per_run = result.n_volumes_per_run;
    supplement.rule_definitions = repmat(struct( ...
        'id', '', 'label', '', 'comparison', 'strict_greater_than', ...
        'fd_threshold_mm', NaN), nThresholds, 1);
    supplement.censor_vectors = struct();
    supplement.censor_count_by_run = zeros(80, nThresholds);
    supplement.censor_percent_by_run = zeros(80, nThresholds);

    candidateId = strings(nThresholds, 1);
    candidateLabel = strings(nThresholds, 1);
    nCensored = zeros(nThresholds, 1);
    percentCensored = zeros(nThresholds, 1);
    maximumRunPercent = zeros(nThresholds, 1);
    runsOver5 = zeros(nThresholds, 1);
    runsOver10 = zeros(nThresholds, 1);
    runsOver20 = zeros(nThresholds, 1);

    for thresholdIndex = 1:nThresholds
        threshold = thresholds(thresholdIndex);
        id = ['fd_gt_' thresholdTokens{thresholdIndex}];
        label = sprintf('FD > %.6g mm', threshold);
        candidateId(thresholdIndex) = id;
        candidateLabel(thresholdIndex) = label;
        supplement.rule_definitions(thresholdIndex).id = id;
        supplement.rule_definitions(thresholdIndex).label = label;
        supplement.rule_definitions(thresholdIndex).fd_threshold_mm = threshold;
        supplement.censor_vectors.(id) = cell(80, 1);
        for runIndex = 1:80
            fd = result.fd_mm{runIndex};
            if numel(fd) ~= result.n_volumes_per_run(runIndex) || any(~isfinite(fd))
                error('OX:MotionQC:InvalidSourceFD', ...
                    '%s %s has invalid saved FD.', subject, result.run_ids{runIndex});
            end
            censor = logical(fd > threshold);
            supplement.censor_vectors.(id){runIndex} = censor;
            supplement.censor_count_by_run(runIndex, thresholdIndex) = nnz(censor);
            supplement.censor_percent_by_run(runIndex, thresholdIndex) = 100 * mean(censor);
        end
        nCensored(thresholdIndex) = sum(supplement.censor_count_by_run(:, thresholdIndex));
        percentCensored(thresholdIndex) = 100 * nCensored(thresholdIndex) / ...
            sum(result.n_volumes_per_run);
        maximumRunPercent(thresholdIndex) = max( ...
            supplement.censor_percent_by_run(:, thresholdIndex));
        runsOver5(thresholdIndex) = nnz( ...
            supplement.censor_percent_by_run(:, thresholdIndex) > 5);
        runsOver10(thresholdIndex) = nnz( ...
            supplement.censor_percent_by_run(:, thresholdIndex) > 10);
        runsOver20(thresholdIndex) = nnz( ...
            supplement.censor_percent_by_run(:, thresholdIndex) > 20);
    end

    supplement.subject_rule_summary = table(candidateId, candidateLabel, ...
        nCensored, percentCensored, maximumRunPercent, runsOver5, ...
        runsOver10, runsOver20, 'VariableNames', {'candidate_id', ...
        'candidate_label', 'n_censored', 'percent_censored', ...
        'max_run_percent', 'runs_over_5pct', 'runs_over_10pct', ...
        'runs_over_20pct'});
    supplement.generated_at = char(datetime('now', 'TimeZone', 'local', ...
        'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
    save(subjectFiles{subjectIndex}, 'supplement', '-v7.3');

    thisTable = supplement.subject_rule_summary;
    thisTable.subject = repmat(string(subject), nThresholds, 1);
    longTables{subjectIndex} = movevars(thisTable, 'subject', 'Before', 1);
    overallPercent(subjectIndex, :) = percentCensored';
    maxRunPercent(subjectIndex, :) = maximumRunPercent';
end

candidateSummary = vertcat(longTables{:});
wideNames = matlab.lang.makeValidName( ...
    arrayfun(@(x) sprintf('fd_gt_%s', numberToken(x)), thresholds, ...
    'UniformOutput', false));
percentWide = array2table(overallPercent, 'VariableNames', wideNames);
percentWide = addvars(percentWide, string(compose('subj_%d', subjects))', ...
    'Before', 1, 'NewVariableNames', 'subject');
writetable(candidateSummary, longCsv);
writetable(percentWide, wideCsv);
makeFigure(subjects, thresholds, overallPercent, maxRunPercent, figureFile);

summary = struct('candidate', candidateSummary, ...
    'percent_wide', percentWide, 'subject_files', {subjectFiles}, ...
    'summary_csv', longCsv, 'wide_csv', wideCsv, 'figure', figureFile);
fprintf('Saved supplemental FD threshold analysis in %s\n', outputDir);
end

function makeFigure(subjects, thresholds, overallPercent, maxRunPercent, outputFile)
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [50 50 1100 600]);
cleanup = onCleanup(@() close(figureHandle));
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
labels = compose('FD > %.6g', thresholds);
nexttile;
bar(overallPercent);
xlabel('Participant'); ylabel('Censored volumes (%)');
title('Overall censor percentage');
xticks(1:numel(subjects)); xticklabels(compose('subj_%d', subjects));
legend(labels, 'Location', 'best'); grid on;
nexttile;
bar(maxRunPercent);
xlabel('Participant'); ylabel('Maximum run censor percentage');
title('Maximum in any individual run');
xticks(1:numel(subjects)); xticklabels(compose('subj_%d', subjects));
legend(labels, 'Location', 'best'); grid on;
sgtitle('Supplemental diagnostic FD thresholds (no rule selected)');
exportgraphics(figureHandle, outputFile, 'Resolution', 180);
end

function token = numberToken(value)
token = strrep(sprintf('%.6g', value), '.', 'p');
end

function root = localProjectRoot()
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
