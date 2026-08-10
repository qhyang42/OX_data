function group_searchlight_stouffer(outputRoot, analysisName, subjectIds)
%GROUP_SEARCHLIGHT_STOUFFER Combine warped empirical p maps at group level.
%   Uses one-sided, unweighted Stouffer inference and Benjamini-Hochberg
%   correction across voxels valid for every supplied subject.

arguments
    outputRoot (1, :) char
    analysisName (1, :) char {mustBeMember(analysisName, {'context', 'odor_score'})}
    subjectIds (1, :) double {mustBeInteger, mustBePositive}
end

run_self_tests();

analysisDir = fullfile(outputRoot, analysisName);
validMaskFile = fullfile(analysisDir, 'group_valid_all_subjects.nii');
assert(isfile(validMaskFile), 'Missing valid-group mask: %s', validMaskFile);

Vmask = spm_vol(validMaskFile);
validMask = spm_read_vols(Vmask) > 0.5;
assert(any(validMask(:)), 'No all-subject-valid MNI voxels for %s.', analysisName);

nSubjects = numel(subjectIds);
pValues = nan([size(validMask), nSubjects]);
for subjectIndex = 1:nSubjects
    subjectId = subjectIds(subjectIndex);
    pFile = fullfile(analysisDir, sprintf('subj%d_p_uncorrected_mni.nii', subjectId));
    assert(isfile(pFile), 'Missing warped p map: %s', pFile);

    Vp = spm_vol(pFile);
    assert(isequal(Vp.dim, Vmask.dim) && max(abs(Vp.mat(:) - Vmask.mat(:))) < 1e-5, ...
        'MNI geometry mismatch for %s.', pFile);
    pValues(:, :, :, subjectIndex) = spm_read_vols(Vp);
end

pAtValid = reshape(pValues(repmat(validMask, 1, 1, 1, nSubjects)), [], nSubjects);
assert(all(isfinite(pAtValid), 'all'), 'Non-finite p values within the valid group mask.');
assert(all(pAtValid >= 0 & pAtValid <= 1, 'all'), ...
    'Warped p values must fall within [0, 1] before clipping.');

pAtValid = min(max(pAtValid, eps('double')), 1 - eps('double'));
zAtValid = sum(norminv(1 - pAtValid), 2) / sqrt(nSubjects);
groupPAtValid = 1 - normcdf(zAtValid);
qAtValid = bh_fdr_qvalues(groupPAtValid);

zMap = nan(size(validMask));
pMap = nan(size(validMask));
qMap = nan(size(validMask));
sigMap = nan(size(validMask));
zMap(validMask) = zAtValid;
pMap(validMask) = groupPAtValid;
qMap(validMask) = qAtValid;
sigMap(validMask) = qAtValid < 0.05;

write_float_map(Vmask, zMap, fullfile(analysisDir, 'group_stouffer_z.nii'));
write_float_map(Vmask, pMap, fullfile(analysisDir, 'group_stouffer_p_uncorrected.nii'));
write_float_map(Vmask, qMap, fullfile(analysisDir, 'group_stouffer_q_fdr.nii'));
write_float_map(Vmask, sigMap, fullfile(analysisDir, 'group_stouffer_sig_fdr_q05.nii'));
write_summary(analysisDir, analysisName, nSubjects, validMask, zAtValid, groupPAtValid, qAtValid);
end

function q = bh_fdr_qvalues(p)
p = p(:);
[sortedP, order] = sort(p, 'ascend');
m = numel(sortedP);
sortedQ = sortedP .* m ./ (1:m)';
sortedQ = flipud(cummin(flipud(sortedQ)));
sortedQ = min(sortedQ, 1);
q = nan(m, 1);
q(order) = sortedQ;
end

function write_float_map(referenceHeader, values, filename)
header = referenceHeader;
header.fname = filename;
header.dt = [16, 0];
header.pinfo = [1; 0; 0];
spm_write_vol(header, single(values));
end

function write_summary(analysisDir, analysisName, nSubjects, validMask, z, p, q)
summaryFile = fullfile(analysisDir, 'group_statistics_summary.tsv');
fid = fopen(summaryFile, 'w');
assert(fid ~= -1, 'Could not write summary: %s', summaryFile);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'analysis\t%s\n', analysisName);
fprintf(fid, 'subjects\t%d\n', nSubjects);
fprintf(fid, 'valid_voxels\t%d\n', nnz(validMask));
fprintf(fid, 'fdr_alpha\t0.05\n');
fprintf(fid, 'significant_voxels_q_lt_0.05\t%d\n', nnz(q < 0.05));
fprintf(fid, 'stouffer_z_min\t%.9g\n', min(z));
fprintf(fid, 'stouffer_z_max\t%.9g\n', max(z));
fprintf(fid, 'group_p_min\t%.9g\n', min(p));
fprintf(fid, 'group_p_max\t%.9g\n', max(p));
fprintf(fid, 'group_q_min\t%.9g\n', min(q));
fprintf(fid, 'group_q_max\t%.9g\n', max(q));
end

function run_self_tests()
% Directional Stouffer: two p=.01 values must yield stronger evidence.
testP = [0.01; 0.01];
testGroupP = 1 - normcdf(sum(norminv(1 - testP)) / sqrt(numel(testP)));
assert(testGroupP < 0.01, 'Stouffer directionality self-test failed.');

% Clipping keeps boundary p-values finite before conversion to z scores.
edgeP = min(max([0; 1], eps('double')), 1 - eps('double'));
assert(all(isfinite(norminv(1 - edgeP))), 'P-value clipping self-test failed.');

% BH calculation has known adjusted values for this small example.
testQ = bh_fdr_qvalues([0.01; 0.04; 0.03]);
assert(max(abs(testQ - [0.03; 0.04; 0.04])) < 1e-12, 'BH-FDR self-test failed.');
end
