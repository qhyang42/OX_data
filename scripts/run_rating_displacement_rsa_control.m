function results = run_rating_displacement_rsa_control()
% Existing control RDMs: native GM intersection, no functional restriction.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'environment')); setup_ox;
addpath(fullfile(root,'scripts'));
out=fullfile(root,'RDMs','rating_displacement_RSA_control');
assert(~isfile(fullfile(out,'control_run_complete.txt')),'Control run already complete.');
subjects=2:6; rois=["control_GM_putamen","control_A1"];
neural_dir=fullfile(root,'RDMs','neural_control');
for sid=subjects
 L=load(fullfile(neural_dir,sprintf('subj_%d_neural_RDMs.mat',sid)),'results');N=L.results;
 assert(~N.preprocessing.use_functional_restriction && isempty(N.preprocessing.functional_mask));
 assert(isequal(string({N.roi_results.roi_name}),rois));
 assert(all(N.roi_metadata.n_usable_voxels>=10));
end
OX_test_formal_rsa();
if isfile(fullfile(out,'results.mat'))
 L=load(fullfile(out,'results.mat'),'results');results=L.results;
 assert(isequal(results.roi_names,rois)&&isequal(results.subject_ids,subjects));
 assert(results.options.NumPermutations==5000&&string(results.options.NeuralDir)==string(neural_dir));
else
results=run_rating_displacement_rsa('SubjectIDs',subjects,'ROINames',rois, ...
 'NeuralDir',neural_dir,'RequireExploratoryReference',false, ...
 'NumPermutations',5000,'PermutationSeed',1001,'OutputDir',out, ...
 'MakePlots',true,'RunMixedModels',false);
end
results.mixed_models=OX_rsa_context_displacement_lme(results.observations,out,true,true);
results.options.RunMixedModels=true;
OX_rsa_write_readme(results);
L=load(fullfile(root,'RDMs','omnibus_RSA_control','results.mat'),'results');
assert(isequal(results.permutation_mappings,L.results.permutation_mappings));
for k=1:2
 M=results.models(k);
 assert(isequal(size(M.null_beta),[5 2 2 5000]));
 assert(all(isfinite(M.null_beta),'all'));
 assert(all(M.fit_statistics.n_valid==120) && all(M.fit_statistics.n_dropped==0));
end
assert(height(results.mixed_models.slopes)+height(results.mixed_models.failed_models)==12);
results.control_restriction='Native anatomical ROI AND GM AND finite betas; no functional restriction.';
results.fdr_family='Four control ROI x rating coefficients per model, separately by subject/group; 12 pooled mixed-model slopes as a separate family.';
save(fullfile(out,'results.mat'),'results','-v7.3');
fid=fopen(fullfile(out,'README.md'),'a');
fprintf(fid,'\n## Control-specific run\n\nPutamen and A1, subjects 2–6. Reuses `RDMs/neural_control/` with GM intersection and no functional restriction. No subjects or conditions dropped. FDR: four control ROI × rating tests per RSA model, separately within each subject and at group level; twelve mixed-model slopes in a separate family. All 5000 condition-correspondence mappings match control omnibus RSA. Original olfactory analyses are unchanged. Reproduce with `scripts/run_rating_displacement_rsa_control.m`.\n');fclose(fid);
for k=1:2, disp(results.models(k).name);disp(results.models(k).group_statistics);end
disp(results.mixed_models.slopes);
fid=fopen(fullfile(out,'control_run_complete.txt'),'w');fprintf(fid,'RSA complete; mixed-model successes=%d, failures=%d. See diagnostics.\n',height(results.mixed_models.slopes),height(results.mixed_models.failed_models));fclose(fid);
end
