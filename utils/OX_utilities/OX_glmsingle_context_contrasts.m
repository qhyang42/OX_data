function manifest = OX_glmsingle_context_contrasts(subjectId, eventType, varargin)
%OX_GLMSINGLE_CONTEXT_CONTRASTS Write context contrasts from GLMsingle betas.
%   MANIFEST = OX_GLMSINGLE_CONTEXT_CONTRASTS(SUBJECTID, EVENTTYPE) pools
%   the final GLMsingle single-trial betas within PERSON, FOOD, LOCATION,
%   and CONTROL and writes three target-versus-mean-other-context contrast
%   images in the corresponding GLMsingle estimates directory.
%
%   EVENTTYPE must be 'sniff' or 'countdown'. Optional name-value inputs:
%     ProjectRoot - OX_DATA project directory (inferred from this file)
%     Overwrite   - replace existing contrast outputs (default: false)

%   These are contrast-estimate images, analogous to SPM con images. They
%   are not SPM t-statistic images because GLMsingle's saved beta output
%   does not contain the residual covariance needed for an SPM t contrast.

%   Invalid GLMsingle voxel rows and voxels outside the gray-matter mask are
%   written as zero. contrast_valid_voxels.nii identifies valid estimates.

%   See also OX_LOAD_TRIAL_METADATA.


p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'subjectId', @(x) isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x == round(x) && x >= 1);
addRequired(p, 'eventType', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'ProjectRoot', '', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'Overwrite', false, ...
    @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
parse(p, subjectId, eventType, varargin{:});
opts = p.Results;

eventType = lower(strtrim(string(eventType)));
assert(ismember(eventType, ["sniff", "countdown"]), ...
    'OX:GLMsingleContrasts:InvalidEventType', ...
    'eventType must be sniff or countdown.');

if strlength(string(opts.ProjectRoot)) == 0
    utilityDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(fileparts(utilityDir));
else
    projectRoot = char(string(opts.ProjectRoot));
end
assert(isfolder(projectRoot), 'OX:GLMsingleContrasts:MissingProjectRoot', ...
    'Project root does not exist: %s', projectRoot);
assert(exist('spm_vol', 'file') == 2 && exist('spm_write_vol', 'file') == 2, ...
    'OX:GLMsingleContrasts:MissingSPM', ...
    'SPM must be on the MATLAB path before running %s.', mfilename);

run_contrast_self_test();

subjectName = sprintf('subj_%d', subjectId);
subjectNiftiDir = fullfile(projectRoot, 'MRI', subjectName, 'nifti');
estimateDir = fullfile(subjectNiftiDir, sprintf( ...
    '%s_single_trial_by_category_physio', char(eventType)));
betaFile = fullfile(estimateDir, 'TYPED_FITHRF_GLMDENOISE_RR.mat');
designFile = fullfile(estimateDir, 'DESIGNINFO.mat');
maskFile = fullfile(subjectNiftiDir, 'coreg', 'gm_mask_thr05_func.nii');

requiredFiles = {betaFile, designFile, maskFile};
for fileIndex = 1:numel(requiredFiles)
    assert(isfile(requiredFiles{fileIndex}), ...
        'OX:GLMsingleContrasts:MissingInput', ...
        'Required input does not exist: %s', requiredFiles{fileIndex});
end

contrastFiles = {
    fullfile(estimateDir, 'con_0001_person_gt_mean_others.nii')
    fullfile(estimateDir, 'con_0002_food_gt_mean_others.nii')
    fullfile(estimateDir, 'con_0003_location_gt_mean_others.nii')
};
validMaskFile = fullfile(estimateDir, 'contrast_valid_voxels.nii');
manifestFile = fullfile(estimateDir, 'contrast_manifest.tsv');
outputFiles = [contrastFiles; {validMaskFile}; {manifestFile}];
existingOutputs = outputFiles(cellfun(@isfile, outputFiles));
if ~logical(opts.Overwrite) && ~isempty(existingOutputs)
    error('OX:GLMsingleContrasts:WouldOverwrite', ...
        'Refusing to overwrite existing output(s):\n%s', ...
        strjoin(existingOutputs, newline));
end

designInfo = load(designFile, 'stimorder', 'condcounts');
assert(isfield(designInfo, 'stimorder') && isfield(designInfo, 'condcounts'), ...
    'OX:GLMsingleContrasts:InvalidDesignInfo', ...
    'DESIGNINFO.mat must contain stimorder and condcounts: %s', designFile);
stimOrder = double(designInfo.stimorder(:));
conditionCounts = double(designInfo.condcounts(:)');
assert(numel(stimOrder) == 800, 'OX:GLMsingleContrasts:TrialCount', ...
    'Expected 800 trials in stimorder; found %d.', numel(stimOrder));
assert(isequal(conditionCounts, [200, 200, 200, 200]), ...
    'OX:GLMsingleContrasts:ConditionCounts', ...
    'Expected DESIGNINFO condcounts [200 200 200 200]; found %s.', ...
    mat2str(conditionCounts));
assert(all(ismember(stimOrder, 1:4)), ...
    'OX:GLMsingleContrasts:ConditionLabels', ...
    'stimorder contains labels outside 1:4.');
observedCounts = arrayfun(@(condition) nnz(stimOrder == condition), 1:4);
assert(isequal(observedCounts, [200, 200, 200, 200]), ...
    'OX:GLMsingleContrasts:ConditionCounts', ...
    'Expected 200 stimorder entries per context; found %s.', ...
    mat2str(observedCounts));

trialMetadata = OX_load_trial_metadata(subjectId, ...
    'CueRoot', fullfile(projectRoot, 'cuelist'), ...
    'ExpectedTrials', 800, 'ExpectedRuns', 80, 'ExpectedOdors', 20);
contextOrder = ["PERSON", "FOOD", "LOCATION", "CONTROL"];
metadataOrder = zeros(height(trialMetadata), 1);
for conditionIndex = 1:numel(contextOrder)
    metadataOrder(trialMetadata.context == contextOrder(conditionIndex)) = conditionIndex;
end
assert(isequal(stimOrder, metadataOrder), ...
    'OX:GLMsingleContrasts:MetadataMismatch', ...
    ['GLMsingle stimorder does not match OX_load_trial_metadata. ' ...
     'Contrast labels cannot be assigned safely.']);

loadedBetas = load(betaFile, 'modelmd');
assert(isfield(loadedBetas, 'modelmd'), ...
    'OX:GLMsingleContrasts:MissingModel', ...
    'GLMsingle file does not contain modelmd: %s', betaFile);
modelmd = reshape(loadedBetas.modelmd, size(loadedBetas.modelmd, 1), []);
clear loadedBetas
assert(size(modelmd, 2) == 800, 'OX:GLMsingleContrasts:BetaTrialCount', ...
    'Expected 800 modelmd trial columns; found %d.', size(modelmd, 2));

maskHeader = spm_vol(maskFile);
assert(isscalar(maskHeader), 'OX:GLMsingleContrasts:MaskVolumes', ...
    'Gray-matter mask must contain exactly one volume: %s', maskFile);
grayMatterMask = spm_read_vols(maskHeader) > 0;
grayMatterIndices = find(grayMatterMask);
assert(size(modelmd, 1) == numel(grayMatterIndices), ...
    'OX:GLMsingleContrasts:MaskMapping', ...
    ['modelmd has %d voxel rows but the gray-matter mask contains %d voxels. ' ...
     'Row-to-voxel mapping is unsafe.'], size(modelmd, 1), numel(grayMatterIndices));

finiteByTrial = isfinite(modelmd);
validRows = all(finiteByTrial, 2);
partiallyInvalidRows = any(finiteByTrial, 2) & ~validRows;
assert(~any(partiallyInvalidRows), ...
    'OX:GLMsingleContrasts:PartialNonfiniteRows', ...
    ['Found %d modelmd rows containing a mixture of finite and non-finite ' ...
     'trials; refusing to average incomplete context data.'], ...
    nnz(partiallyInvalidRows));
assert(any(validRows), 'OX:GLMsingleContrasts:NoValidVoxels', ...
    'No modelmd rows contain finite estimates for every trial.');
clear finiteByTrial partiallyInvalidRows

categoryMeans = zeros(size(modelmd, 1), 4, 'single');
for conditionIndex = 1:4
    categoryMeans(validRows, conditionIndex) = mean( ...
        modelmd(validRows, stimOrder == conditionIndex), 2);
end

[contrastNames, contrastWeights] = contrast_definitions();
contrastValues = zeros(size(modelmd, 1), 3, 'single');
contrastValues(validRows, :) = single( ...
    double(categoryMeans(validRows, :)) * contrastWeights');

validVolume = zeros(maskHeader.dim, 'uint8');
validVolume(grayMatterIndices(validRows)) = 1;
write_volume(maskHeader, validMaskFile, validVolume, 2, ...
    sprintf('Valid GLMsingle context contrast voxels: %s %s', ...
    subjectName, char(eventType)));

for contrastIndex = 1:3
    contrastVolume = zeros(maskHeader.dim, 'single');
    contrastVolume(grayMatterIndices(validRows)) = ...
        contrastValues(validRows, contrastIndex);
    write_volume(maskHeader, contrastFiles{contrastIndex}, contrastVolume, 16, ...
        sprintf('GLMsingle %% signal contrast: %s', contrastNames(contrastIndex)));
end

manifest = build_manifest(subjectId, eventType, contrastFiles, contrastNames, ...
    contrastWeights, observedCounts, nnz(validRows), nnz(~validRows), ...
    betaFile, designFile, maskFile, validMaskFile);
writetable(manifest, manifestFile, 'FileType', 'text', 'Delimiter', '\t');

validate_written_outputs(maskHeader, grayMatterMask, grayMatterIndices, ...
    validRows, modelmd, stimOrder, contrastWeights, contrastFiles, validMaskFile);

fprintf('[%s | %s] Wrote 3 contrasts (%d valid, %d invalid mask voxels).\n', ...
    subjectName, upper(char(eventType)), nnz(validRows), nnz(~validRows));
fprintf('  %s\n', estimateDir);
end

function [names, weights] = contrast_definitions()
names = [
    "person_gt_mean_others"
    "food_gt_mean_others"
    "location_gt_mean_others"
];
weights = [
     1,    -1/3, -1/3, -1/3
    -1/3,   1,   -1/3, -1/3
    -1/3,  -1/3,  1,   -1/3
];
end

function run_contrast_self_test()
[~, weights] = contrast_definitions();
syntheticMeans = [4, 1, 1, 1; 1, 4, 1, 1; 1, 1, 4, 1];
observed = diag(syntheticMeans * weights');
assert(max(abs(observed - 3)) < 1e-12, ...
    'OX:GLMsingleContrasts:SelfTest', ...
    'Target-versus-mean-other contrast arithmetic self-test failed.');
assert(max(abs(sum(weights, 2))) < 1e-12, ...
    'OX:GLMsingleContrasts:SelfTest', ...
    'Context contrast weights must sum to zero.');
end

function write_volume(referenceHeader, filename, values, datatype, description)
header = referenceHeader;
header.fname = filename;
header.dt = [datatype, 0];
header.pinfo = [1; 0; 0];
header.descrip = description;
spm_write_vol(header, values);
end

function manifest = build_manifest(subjectId, eventType, contrastFiles, ...
        contrastNames, weights, counts, nValid, nInvalid, betaFile, ...
        designFile, maskFile, validMaskFile)
contrast_index = (1:3)';
subject_id = repmat(subjectId, 3, 1);
event_type = repmat(eventType, 3, 1);
contrast_name = contrastNames;
filename = string(contrastFiles);
weight_person = weights(:, 1);
weight_food = weights(:, 2);
weight_location = weights(:, 3);
weight_control = weights(:, 4);
n_person = repmat(counts(1), 3, 1);
n_food = repmat(counts(2), 3, 1);
n_location = repmat(counts(3), 3, 1);
n_control = repmat(counts(4), 3, 1);
valid_voxels = repmat(nValid, 3, 1);
invalid_mask_voxels = repmat(nInvalid, 3, 1);
units = repmat("percent_signal_change", 3, 1);
source_beta_file = repmat(string(betaFile), 3, 1);
source_design_file = repmat(string(designFile), 3, 1);
source_gray_matter_mask = repmat(string(maskFile), 3, 1);
valid_voxel_mask = repmat(string(validMaskFile), 3, 1);

manifest = table(contrast_index, subject_id, event_type, contrast_name, ...
    filename, weight_person, weight_food, weight_location, weight_control, ...
    n_person, n_food, n_location, n_control, valid_voxels, ...
    invalid_mask_voxels, units, source_beta_file, source_design_file, ...
    source_gray_matter_mask, valid_voxel_mask);
end

function validate_written_outputs(referenceHeader, grayMatterMask, ...
        grayMatterIndices, validRows, modelmd, stimOrder, weights, ...
        contrastFiles, validMaskFile)
validHeader = spm_vol(validMaskFile);
assert_geometry(referenceHeader, validHeader, validMaskFile);
assert(validHeader.dt(1) == 2, 'OX:GLMsingleContrasts:OutputDatatype', ...
    'Expected uint8 valid mask: %s', validMaskFile);
writtenValidMask = spm_read_vols(validHeader) > 0;
expectedValidMask = false(referenceHeader.dim);
expectedValidMask(grayMatterIndices(validRows)) = true;
assert(isequal(writtenValidMask, expectedValidMask), ...
    'OX:GLMsingleContrasts:ValidMaskMismatch', ...
    'Written valid-voxel mask does not match expected GLMsingle validity.');

validRowIndices = find(validRows);
nSamples = min(25, numel(validRowIndices));
samplePositions = unique(round(linspace(1, numel(validRowIndices), nSamples)));
sampleRows = validRowIndices(samplePositions);
sampleVoxelIndices = grayMatterIndices(sampleRows);

for contrastIndex = 1:3
    outputHeader = spm_vol(contrastFiles{contrastIndex});
    assert_geometry(referenceHeader, outputHeader, contrastFiles{contrastIndex});
    assert(outputHeader.dt(1) == 16, ...
        'OX:GLMsingleContrasts:OutputDatatype', ...
        'Expected float32 contrast image: %s', contrastFiles{contrastIndex});
    outputVolume = spm_read_vols(outputHeader);
    assert(all(isfinite(outputVolume(:))), ...
        'OX:GLMsingleContrasts:NonfiniteOutput', ...
        'Contrast output contains non-finite values: %s', ...
        contrastFiles{contrastIndex});
    assert(all(outputVolume(~expectedValidMask) == 0), ...
        'OX:GLMsingleContrasts:InvalidVoxelValues', ...
        'Contrast output is nonzero outside the valid mask: %s', ...
        contrastFiles{contrastIndex});
    assert(all(outputVolume(~grayMatterMask) == 0), ...
        'OX:GLMsingleContrasts:OutsideMaskValues', ...
        'Contrast output is nonzero outside the gray-matter mask: %s', ...
        contrastFiles{contrastIndex});

    directValues = zeros(numel(sampleRows), 1);
    for conditionIndex = 1:4
        directValues = directValues + weights(contrastIndex, conditionIndex) .* ...
            mean(double(modelmd(sampleRows, stimOrder == conditionIndex)), 2);
    end
    writtenValues = outputVolume(sampleVoxelIndices);
    tolerance = 1e-5 * max(1, max(abs(directValues)));
    assert(max(abs(writtenValues - directValues)) <= tolerance, ...
        'OX:GLMsingleContrasts:ValueMismatch', ...
        'Written values do not match direct beta contrasts: %s', ...
        contrastFiles{contrastIndex});
end
end

function assert_geometry(referenceHeader, outputHeader, filename)
assert(isequal(outputHeader.dim, referenceHeader.dim), ...
    'OX:GLMsingleContrasts:OutputDimensions', ...
    'Output dimensions do not match the reference mask: %s', filename);
assert(max(abs(outputHeader.mat(:) - referenceHeader.mat(:))) < 1e-6, ...
    'OX:GLMsingleContrasts:OutputAffine', ...
    'Output affine does not match the reference mask: %s', filename);
end
