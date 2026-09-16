function results = run_rating_displacement_rsa_TU()
% TU-only overall displacement RSA using established TU RDMs.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'scripts'),fullfile(root,'utils','OX_utilities'));
neural_dir=fullfile(root,'RDMs','neural_TU');
qc=readtable(fullfile(root,'MRI','group','context_split_half_similarity_physio_TU', ...
    'context_split_half_similarity_summary.csv'),'TextType','string');
assert(isequal(qc.subject(:)',2:6) && all(qc.roi=="olf_TU"));
subjects=qc.subject(qc.voxel_count>=10)';
assert(isequal(subjects,[2 3 4 6]));
assert(all(arrayfun(@(s) isfile(fullfile(neural_dir,sprintf('subj_%d_neural_RDMs.mat',s))),subjects)));
for sid=subjects
    N=load(fullfile(neural_dir,sprintf('subj_%d_neural_RDMs.mat',sid)),'results'); N=N.results;
    assert(N.subject_id==sid && numel(N.trial_metadata.run_id)==800);
    assert(numel(unique(N.trial_metadata.run_id))==80 && all(groupcounts(N.trial_metadata.run_id)==10));
    R=N.roi_results(strcmp(string({N.roi_results.roi_name}),"olf_TU"));
    assert(isscalar(R) && R.n_voxels==qc.voxel_count(qc.subject==sid) && R.n_voxels>=10);
end
out=fullfile(root,'RDMs','rating_displacement_RSA_TU');
assert(~isfolder(out),'Output exists; preserve established results.');
results=run_rating_displacement_rsa('SubjectIDs',subjects,'ROINames',"olf_TU", ...
    'NeuralDir',neural_dir,'RequireExploratoryReference',false, ...
    'NumPermutations',5000,'PermutationSeed',1001,'RunMixedModels',false,'OutputDir',out);
assert(~isfield(results,'mixed_models'));
results.excluded_subjects=qc(qc.voxel_count<10,{'subject','roi','voxel_count','status'});
results.usable_features=qc(:,{'subject','roi','voxel_count','status'});
results.fdr_family='Two TU rating predictors per overall adjustment model, separate from original five-ROI families.';
% Verify identical condition permutations to the existing production TU run.
O=load(fullfile(root,'RDMs','omnibus_RSA_TU','results.mat'),'results');
assert(isequal(results.subject_ids,O.results.subject_ids));
assert(isequal(results.permutation_mappings,O.results.permutation_mappings));
assert(height(results.observations)==480 && all(groupcounts(results.observations.subject_id)==120));
for m=1:numel(results.models)
    M=results.models(m);
    assert(isequal(size(M.null_beta),[4 1 2 5000]) && all(isfinite(M.null_beta),'all'));
    assert(all(M.fit_statistics.n_valid==120) && all(M.fit_statistics.n_dropped==0));
end
save(fullfile(out,'results.mat'),'results','-v7.3');
writetable(results.excluded_subjects,fullfile(out,'tables','excluded_subjects.csv'));
writetable(results.usable_features,fullfile(out,'tables','usable_features.csv'));
fid=fopen(fullfile(out,'README.md'),'a'); assert(fid>=0);
fprintf(fid,['\n## TU run\n\nSubjects 2, 3, 4, 6: 23, 12, 11, 39 usable voxels respectively. ' ...
    'Subject 5 excluded with 6 usable voxels (unchanged minimum 10). TU uses the established gray-matter and positive Odor > Rest p < .001 restriction. ' ...
    'BH-FDR covers two TU rating predictors separately per overall adjustment model and separately from the original five-ROI families. ' ...
    'All 5,000 condition permutations exactly match the previous TU omnibus run. ' ...
    'Overall models retain both context-pair-plus-odor and regular context-pair adjustment. No context-specific fits ran. ' ...
    'Original results are preserved. Reproduce with scripts/run_rating_displacement_rsa_TU.m.\n']);
fclose(fid);
end
