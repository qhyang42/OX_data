function results = run_omnibus_rsa_control()
% Putamen/A1 omnibus RSA: GM intersection, no functional restriction.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'environment')); setup_ox;
addpath(fullfile(root,'scripts'));
neural_dir=fullfile(root,'RDMs','neural_control');
output_dir=fullfile(root,'RDMs','omnibus_RSA_control');
assert(~isfolder(neural_dir) && ~isfolder(output_dir), ...
    'Control output directory already exists; preserve existing outputs.');
subjects=2:6; rois=["control_GM_putamen","control_A1"];
OX_test_formal_rsa();
build_control_rdms(root,neural_dir,subjects,rois);
expected=[59 78;37 95;46 66;107 42;103 13];
for s=subjects
    N=load(fullfile(neural_dir,sprintf('subj_%d_neural_RDMs.mat',s)),'results');
    assert(isequal(N.results.roi_metadata.n_usable_voxels',expected(s-1,:)));
    assert(isempty(N.results.preprocessing.functional_mask));
    assert(~N.results.preprocessing.use_functional_restriction);
    assert(height(N.results.trial_metadata)==800);
    assert(numel(unique(N.results.trial_metadata.run_id))==80);
end
results=run_omnibus_rsa('SubjectIDs',subjects,'ROINames',rois, ...
    'NeuralDir',neural_dir,'RequireExploratoryReference',false, ...
    'NumPermutations',5000,'PermutationSeed',1001,'OutputDir',output_dir);
results.fdr_family='Eight control ROI x scientific predictor tests; separately per subject and group, independent of the prior olfactory family.';
results.control_restriction='Anatomical ROI AND native gray matter AND finite betas; no functional restriction.';
% Verify identical permutation mappings to the established omnibus analysis.
old=load(fullfile(root,'RDMs','omnibus_RSA','results.mat'),'results');
assert(isequal(results.permutation_mappings,old.results.permutation_mappings));
assert(all(results.models.fit_statistics.n_valid==3160));
assert(all(results.models.fit_statistics.n_dropped==0));
assert(all(isfinite(results.models.null_beta),'all'));
assert(isequal(size(results.models.null_beta),[5 2 4 5000]));
save(fullfile(output_dir,'results.mat'),'results','-v7.3');
fid=fopen(fullfile(output_dir,'README.md'),'a');
fprintf(fid,'\n## Control-specific run\n\nPutamen and A1; subjects 2–6. New RDMs in `RDMs/neural_control/`. Native GM intersection with no functional restriction. All 800 trials and all four contexts enter within-run centering and the 80-condition RDMs. This differs from semantic-only split-half normalization. Simple-distance RDMs enter the omnibus model; crossnobis estimation is disabled because it is not used by this model.\n\nBH-FDR covers eight control ROI × predictor tests, separately within each subject and at group level, separate from the original olfactory family. All 5000 permutation mappings match the existing omnibus analysis. No subject or condition is excluded. Input QC records pre-existing missing behavioral ratings and condition-count imbalance; existing behavioral RDMs are reused. Reproduce with `scripts/run_omnibus_rsa_control.m`.\n');
fclose(fid);
disp(results.models.group_statistics);
end
function build_control_rdms(root,neural_dir,subjects,rois)
neural_rdm_options=struct('SubjectIDs',subjects,'ROINames',rois, ...
    'ROISelection','control','OutputDir',neural_dir,'UseFunctionalRestriction',false,'ComputeCrossnobis',false);
run(fullfile(root,'scripts','construct_neural_rdms.m'));
end
