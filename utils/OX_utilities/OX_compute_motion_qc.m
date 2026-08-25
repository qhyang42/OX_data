function results = OX_compute_motion_qc(subjidx, varargin)
%OX_COMPUTE_MOTION_QC Calculate run-wise FD, raw DVARS, and censor diagnostics.
%   RESULTS = OX_COMPUTE_MOTION_QC(SUBJIDX) processes the 80 functional
%   runs for one OX participant and saves motion_qc/subj_N_motion_qc.mat.
%   Existing result files are never overwritten.
%
%   Robust standardized DVARS is calculated independently for every run.
%   Let x be raw DVARS for frames 2:N (frame 1 has no temporal difference),
%   m = median(x), and s = 1.4826 * median(abs(x-m)). Then robust DVARS z is
%   (DVARS-m)/s for frames 2:N, with frame 1 explicitly set to zero. Raw
%   DVARS is retained unchanged. A zero/nonfinite robust scale is an error.
%
%   Name-value options:
%     ProjectRoot  OX_DATA root (inferred from this file by default)
%     OutputDir    Output directory (PROJECTROOT/motion_qc by default)
%     SaveResult   Save the .mat result (true by default)
%     ChunkVoxels  Maximum masked voxels read at once (20000 by default)

p = inputParser;
p.addRequired('subjidx', @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x == fix(x));
p.addParameter('ProjectRoot', localProjectRoot(), @(x) ischar(x) || isstring(x));
p.addParameter('OutputDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('SaveResult', true, @(x) islogical(x) && isscalar(x));
p.addParameter('ChunkVoxels', 20000, @(x) isnumeric(x) && isscalar(x) && x > 0 && x == fix(x));
p.parse(subjidx, varargin{:});

projectRoot = char(p.Results.ProjectRoot);
outputDir = char(p.Results.OutputDir);
if isempty(outputDir)
    outputDir = fullfile(projectRoot, 'motion_qc');
end
subject = sprintf('subj_%d', subjidx);
resultFile = fullfile(outputDir, sprintf('%s_motion_qc.mat', subject));
if p.Results.SaveResult && isfile(resultFile)
    error('OX:MotionQC:WouldOverwrite', ...
        'Refusing to overwrite existing result: %s', resultFile);
end

assertSpmAvailable();
runs = OX_discover_functional_runs(subjidx, 'ProjectRoot', projectRoot, 'ExpectedRuns', 80);
maskFile = fullfile(projectRoot, 'MRI', subject, 'nifti', 'coreg', 'gm_mask_thr05_func.nii');
if ~isfile(maskFile)
    error('OX:MotionQC:MissingMask', 'GM mask is missing: %s', maskFile);
end
maskHeader = spm_vol(maskFile);
if numel(maskHeader) ~= 1
    error('OX:MotionQC:InvalidMask', 'GM mask must be a single 3-D volume: %s', maskFile);
end
maskVolume = spm_read_vols(maskHeader);
maskIndex = find(isfinite(maskVolume) & maskVolume > 0);
if isempty(maskIndex)
    error('OX:MotionQC:NoMaskVoxels', 'No finite positive voxels in GM mask: %s', maskFile);
end
[maskI, maskJ, maskK] = ind2sub(maskHeader.dim, maskIndex);
maskXYZ = [maskI'; maskJ'; maskK'];

rules = candidateRules();
nRuns = numel(runs);
nRules = numel(rules);
fd = cell(nRuns, 1);
rawDvars = cell(nRuns, 1);
robustDvarsZ = cell(nRuns, 1);
censorVectors = struct();
for ruleIndex = 1:nRules
    censorVectors.(rules(ruleIndex).id) = cell(nRuns, 1);
end
nVolumes = zeros(nRuns, 1);
nValidVoxels = zeros(nRuns, 1);
dvarsMedian = zeros(nRuns, 1);
dvarsMad = zeros(nRuns, 1);
dvarsRobustScale = zeros(nRuns, 1);
censorCountByRun = zeros(nRuns, nRules);
censorPercentByRun = zeros(nRuns, nRules);

fprintf('Computing motion QC for %s (%d runs)\n', subject, nRuns);
for runIndex = 1:nRuns
    runInfo = runs(runIndex);
    if ~isfile(runInfo.functional_file)
        error('OX:MotionQC:MissingFunctionalFile', ...
            'Run %s functional file is missing: %s', runInfo.id, runInfo.functional_file);
    end
    if ~isfile(runInfo.motion_file)
        error('OX:MotionQC:MissingMotionFile', ...
            'Run %s motion file is missing: %s', runInfo.id, runInfo.motion_file);
    end

    functionalHeaders = spm_vol(runInfo.functional_file);
    if isempty(functionalHeaders)
        error('OX:MotionQC:EmptyFunctionalFile', ...
            'No volumes found in %s.', runInfo.functional_file);
    end
    validateGeometry(maskHeader, functionalHeaders(1), runInfo.id);
    nVolumes(runIndex) = numel(functionalHeaders);

    motion = readmatrix(runInfo.motion_file, 'FileType', 'text');
    if ~isnumeric(motion) || size(motion, 2) ~= 6 || any(~isfinite(motion(:)))
        error('OX:MotionQC:InvalidMotionFile', ...
            'Motion file must contain only a finite N-by-6 matrix: %s', runInfo.motion_file);
    end
    if size(motion, 1) ~= nVolumes(runIndex)
        error('OX:MotionQC:MotionVolumeMismatch', ...
            'Run %s has %d motion rows but %d functional volumes.', ...
            runInfo.id, size(motion, 1), nVolumes(runIndex));
    end

    deltaMotion = [zeros(1, 6); diff(motion, 1, 1)];
    deltaMotion(:, 4:6) = deltaMotion(:, 4:6) * 50;
    thisFd = sum(abs(deltaMotion), 2);
    thisFd(1) = 0;

    [thisDvars, thisNValid] = calculateDvars(functionalHeaders, maskXYZ, p.Results.ChunkVoxels);
    if numel(thisDvars) ~= nVolumes(runIndex) || numel(thisDvars) ~= size(motion, 1)
        error('OX:MotionQC:DvarsLengthMismatch', ...
            'Run %s FD/DVARS/motion lengths do not agree.', runInfo.id);
    end
    if thisNValid == 0
        error('OX:MotionQC:NoValidVoxels', ...
            'Run %s has no voxels finite at every frame inside the GM mask.', runInfo.id);
    end
    if any(~isfinite(thisFd)) || any(~isfinite(thisDvars))
        error('OX:MotionQC:NonfiniteMetric', ...
            'Run %s produced NaN or Inf in FD/DVARS.', runInfo.id);
    end

    baseline = thisDvars(2:end);
    thisMedian = median(baseline);
    thisMad = median(abs(baseline - thisMedian));
    thisScale = 1.4826 * thisMad;
    if ~isfinite(thisScale) || thisScale <= 0
        error('OX:MotionQC:InvalidDvarsScale', ...
            'Run %s has a zero or nonfinite robust DVARS scale.', runInfo.id);
    end
    thisZ = (thisDvars - thisMedian) ./ thisScale;
    thisZ(1) = 0;
    if any(~isfinite(thisZ))
        error('OX:MotionQC:NonfiniteStandardizedDvars', ...
            'Run %s produced NaN or Inf in robust standardized DVARS.', runInfo.id);
    end

    fd{runIndex} = thisFd;
    rawDvars{runIndex} = thisDvars;
    robustDvarsZ{runIndex} = thisZ;
    nValidVoxels(runIndex) = thisNValid;
    dvarsMedian(runIndex) = thisMedian;
    dvarsMad(runIndex) = thisMad;
    dvarsRobustScale(runIndex) = thisScale;
    for ruleIndex = 1:nRules
        censor = applyRule(rules(ruleIndex), thisFd, thisZ);
        censorVectors.(rules(ruleIndex).id){runIndex} = censor;
        censorCountByRun(runIndex, ruleIndex) = nnz(censor);
        censorPercentByRun(runIndex, ruleIndex) = 100 * mean(censor);
    end
    if mod(runIndex, 10) == 0 || runIndex == nRuns
        fprintf('  %s: %d/%d runs complete\n', subject, runIndex, nRuns);
    end
end

allFd = vertcat(fd{:});
allRawDvars = vertcat(rawDvars{:});
allRobustDvarsZ = vertcat(robustDvarsZ{:});
if any(~isfinite([allFd; allRawDvars; allRobustDvarsZ]))
    error('OX:MotionQC:NonfiniteSubjectMetric', ...
        '%s contains nonfinite values after concatenation.', subject);
end

subjectCensorCount = sum(censorCountByRun, 1);
subjectCensorPercent = 100 * subjectCensorCount / sum(nVolumes);
subjectRuleSummary = table(string({rules.id})', string({rules.label})', ...
    subjectCensorCount', subjectCensorPercent', max(censorPercentByRun, [], 1)', ...
    sum(censorPercentByRun > 5, 1)', sum(censorPercentByRun > 10, 1)', ...
    sum(censorPercentByRun > 20, 1)', ...
    'VariableNames', {'candidate_id', 'candidate_label', 'n_censored', ...
    'percent_censored', 'max_run_percent', 'runs_over_5pct', ...
    'runs_over_10pct', 'runs_over_20pct'});

results = struct();
results.subject = subject;
results.subject_index = subjidx;
results.TR_seconds = 0.76;
results.run_ids = {runs.id}';
results.run_order = [(1:nRuns)', [runs.session]', [runs.run]'];
results.run_order_columns = {'ordinal', 'session', 'run'};
results.n_volumes_per_run = nVolumes;
results.fd_mm = fd;
results.raw_dvars = rawDvars;
results.robust_dvars_z = robustDvarsZ;
results.candidate_rules = rules;
results.censor_vectors = censorVectors;
results.censor_count_by_run = censorCountByRun;
results.censor_percent_by_run = censorPercentByRun;
results.subject_rule_summary = subjectRuleSummary;
results.validation = struct();
results.validation.functional_files = {runs.functional_file}';
results.validation.motion_files = {runs.motion_file}';
results.validation.source_stem_match = [runs.source_stem_match]';
results.validation.mask_file = maskFile;
results.validation.mask_dimensions = maskHeader.dim;
results.validation.mask_affine = maskHeader.mat;
results.validation.n_positive_finite_mask_voxels = numel(maskIndex);
results.validation.n_valid_voxels_per_run = nValidVoxels;
results.validation.dvars_run_median = dvarsMedian;
results.validation.dvars_run_unscaled_mad = dvarsMad;
results.validation.dvars_run_robust_scale = dvarsRobustScale;
results.validation.fd_head_radius_mm = 50;
results.validation.frame_differences_reset_at_each_run = true;
results.validation.robust_dvars_definition = [ ...
    'For each run independently, using raw DVARS frames 2:N: ', ...
    'm=median(x), MAD=median(abs(x-m)), scale=1.4826*MAD, ', ...
    'robust_z=(DVARS-m)/scale; frame 1 robust_z=0.'];
results.validation.dvars_voxel_rule = [ ...
    'Mask values must be finite and >0; a voxel contributes only if its ', ...
    'functional values are finite at every frame of that run.'];
results.validation.geometry_tolerance = 1e-4;
results.validation.matlab_version = version;
results.validation.spm_version = spm('Ver');
results.validation.generated_at = char(datetime('now', 'TimeZone', 'local', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
results.validation.code_file = mfilename('fullpath');

if p.Results.SaveResult
    if ~isfolder(outputDir)
        [ok, message] = mkdir(outputDir);
        if ~ok
            error('OX:MotionQC:CannotCreateOutput', ...
                'Cannot create %s: %s', outputDir, message);
        end
    end
    if isfile(resultFile)
        error('OX:MotionQC:WouldOverwrite', ...
            'Refusing to overwrite existing result: %s', resultFile);
    end
    save(resultFile, 'results', '-v7.3');
    fprintf('Saved %s\n', resultFile);
end
end

function [dvars, nValid] = calculateDvars(headers, maskXYZ, chunkSize)
nVolumes = numel(headers);
sumSquaredDifference = zeros(nVolumes, 1);
nValid = 0;
for firstVoxel = 1:chunkSize:size(maskXYZ, 2)
    lastVoxel = min(firstVoxel + chunkSize - 1, size(maskXYZ, 2));
    data = spm_get_data(headers, maskXYZ(:, firstVoxel:lastVoxel));
    valid = all(isfinite(data), 1);
    data = data(:, valid);
    if isempty(data)
        continue
    end
    difference = [zeros(1, size(data, 2)); diff(data, 1, 1)];
    sumSquaredDifference = sumSquaredDifference + sum(difference .^ 2, 2);
    nValid = nValid + size(data, 2);
end
if nValid == 0
    dvars = nan(nVolumes, 1);
else
    dvars = sqrt(sumSquaredDifference ./ nValid);
    dvars(1) = 0;
end
end

function rules = candidateRules()
template = struct('id', '', 'label', '', 'type', '', ...
    'fd_threshold_mm', NaN, 'robust_dvars_z_threshold', NaN);
rules = repmat(template, 10, 1);
rules(1) = makeRule('fd_gt_0p2', 'FD > 0.2 mm', 'fd', 0.2, NaN);
rules(2) = makeRule('fd_gt_0p3', 'FD > 0.3 mm', 'fd', 0.3, NaN);
rules(3) = makeRule('fd_gt_0p5', 'FD > 0.5 mm', 'fd', 0.5, NaN);
rules(4) = makeRule('dvarsz_gt_3', 'robust DVARS z > 3', 'dvars', NaN, 3);
rules(5) = makeRule('dvarsz_gt_4', 'robust DVARS z > 4', 'dvars', NaN, 4);
rules(6) = makeRule('dvarsz_gt_5', 'robust DVARS z > 5', 'dvars', NaN, 5);
rules(7) = makeRule('fd0p2_or_dvarsz3', 'FD > 0.2 OR robust DVARS z > 3', 'or', 0.2, 3);
rules(8) = makeRule('fd0p3_or_dvarsz3', 'FD > 0.3 OR robust DVARS z > 3', 'or', 0.3, 3);
rules(9) = makeRule('fd0p3_or_dvarsz4', 'FD > 0.3 OR robust DVARS z > 4', 'or', 0.3, 4);
rules(10) = makeRule('fd0p5_or_dvarsz5', 'FD > 0.5 OR robust DVARS z > 5', 'or', 0.5, 5);
end

function rule = makeRule(id, label, type, fdThreshold, dvarsThreshold)
rule = struct('id', id, 'label', label, 'type', type, ...
    'fd_threshold_mm', fdThreshold, ...
    'robust_dvars_z_threshold', dvarsThreshold);
end

function censor = applyRule(rule, fd, dvarsZ)
switch rule.type
    case 'fd'
        censor = fd > rule.fd_threshold_mm;
    case 'dvars'
        censor = dvarsZ > rule.robust_dvars_z_threshold;
    case 'or'
        censor = fd > rule.fd_threshold_mm | dvarsZ > rule.robust_dvars_z_threshold;
    otherwise
        error('OX:MotionQC:UnknownRuleType', 'Unknown censor rule type: %s', rule.type);
end
censor = logical(censor(:));
end

function validateGeometry(maskHeader, functionalHeader, runId)
if ~isequal(maskHeader.dim, functionalHeader.dim)
    error('OX:MotionQC:GeometryMismatch', ...
        'Run %s dimensions [%s] do not match GM mask [%s].', runId, ...
        num2str(functionalHeader.dim), num2str(maskHeader.dim));
end
if max(abs(maskHeader.mat(:) - functionalHeader.mat(:))) > 1e-4
    error('OX:MotionQC:GeometryMismatch', ...
        'Run %s affine geometry does not match the GM mask (tolerance 1e-4).', runId);
end
end

function assertSpmAvailable()
if exist('spm_vol', 'file') ~= 2 || exist('spm_get_data', 'file') ~= 2
    error('OX:MotionQC:MissingSPM', ...
        'SPM is not on the MATLAB path. Run environment/setup_ox.m first.');
end
end

function root = localProjectRoot()
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
