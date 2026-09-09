function results = run_omnibus_rsa_TU()
% TU-only omnibus RSA; preserve all original RDM and RSA outputs.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'environment')); setup_ox;
addpath(fullfile(root,'scripts'));
neural_dir = fullfile(root,'RDMs','neural_TU');
% Subject 5 has six voxels after the same GM/functional-mask restriction.
% Keep the existing minimum of ten and record this exclusion explicitly.
qc = readtable(fullfile(root,'MRI','group', ...
    'context_split_half_similarity_physio_TU', ...
    'context_split_half_similarity_summary.csv'),'TextType','string');
assert(isequal(qc.subject(:)',2:6) && all(qc.roi=="olf_TU"));
subjects = qc.subject(qc.voxel_count>=10)';
assert(isequal(subjects,[2 3 4 6]));
missing = subjects(~arrayfun(@(s) isfile(fullfile(neural_dir, ...
    sprintf('subj_%d_neural_RDMs.mat',s))),subjects));
if ~isempty(missing)
    build_tu_rdms(root,neural_dir,missing);
end
results = run_omnibus_rsa('SubjectIDs',subjects,'ROINames',"olf_TU", ...
    'NeuralDir',neural_dir,'RequireExploratoryReference',false, ...
    'NumPermutations',5000,'PermutationSeed',1001, ...
    'OutputDir',fullfile(root,'RDMs','omnibus_RSA_TU'));
results.excluded_subjects = qc(qc.voxel_count<10,{'subject','roi','voxel_count','status'});
results.fdr_family = 'Four TU scientific predictors; separate from the original five-ROI family.';
save(fullfile(results.output_dir,'results.mat'),'results','-v7.3');
writetable(results.excluded_subjects,fullfile(results.output_dir,'tables','excluded_subjects.csv'));
fid=fopen(fullfile(results.output_dir,'README.md'),'a');
fprintf(fid,'\nTU run: subject 5 excluded (6 usable voxels; minimum 10). Group N=4. BH-FDR covers the four TU predictors, separately from the original five-ROI results. Existing ROI analyses were not rerun. TU lacks exploratory reference coefficients; independent SVD verification was used.\n');
fclose(fid);
end
function build_tu_rdms(root,neural_dir,subjects)
neural_rdm_options = struct('SubjectIDs',subjects,'ROINames',"olf_TU", ...
    'ROISelection','secondary','OutputDir',neural_dir);
run(fullfile(root,'scripts','construct_neural_rdms.m'));
end
