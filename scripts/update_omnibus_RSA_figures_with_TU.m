function update_omnibus_RSA_figures_with_TU()
% Combine saved estimates for plotting only; retain original inferential families.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'utils','OX_utilities'));
folder=fullfile(root,'RDMs','omnibus_RSA');
A=load(fullfile(folder,'results.mat'),'results');
B=load(fullfile(root,'RDMs','omnibus_RSA_TU','results.mat'),'results');
A=A.results; B=B.results;
assert(isequal(B.roi_names,"olf_TU") && ~any(A.roi_names=="olf_TU"));
assert(isequal(A.models.predictors,B.models.predictors));
assert(isequal(A.subject_ids,2:6) && isequal(B.subject_ids,[2 3 4 6]));
assert(isequal(A.options.NumPermutations,B.options.NumPermutations));
plot_results=struct('analysis',"omnibus",'roi_names',[A.roi_names,B.roi_names], ...
    'subject_ids',A.subject_ids,'output_dir',folder);
plot_results.models=struct('name',A.models.name,'predictors',A.models.predictors, ...
    'beta',nan(5,6,4),'group_statistics',[A.models.group_statistics;B.models.group_statistics]);
plot_results.models.beta(:,1:5,:)=A.models.beta;
[found,rows]=ismember(B.subject_ids,A.subject_ids); assert(all(found));
plot_results.models.beta(rows,6,:)=B.models.beta;
assert(all(isnan(plot_results.models.beta(A.subject_ids==5,6,:)),'all'));
assert(isequaln(plot_results.models.beta(:,1:5,:),A.models.beta));
assert(isequaln(plot_results.models.beta(rows,6,:),B.models.beta));
plot_results.figure_note=["Signed-beta BH-FDR: original ROIs, 20 tests; TU, 4 separate tests"; ...
    "* q < 0.05, ** q < 0.01, *** q < 0.001 | TU excludes subject 5 (<10 voxels)"];
% Display TU first, followed by the original ROI order.
roi_order=[6,1:5];
plot_results.roi_names=plot_results.roi_names(roi_order);
plot_results.models.beta=plot_results.models.beta(:,roi_order,:);
OX_plot_formal_rsa(plot_results);
save(fullfile(folder,'figures','combined_with_TU_plot_data.mat'),'plot_results','-v7.3');
fid=fopen(fullfile(folder,'figures','README.md'),'w'); assert(fid>=0);
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,['# Omnibus RSA figures including TU\n\n' ...
    'Generated from saved omnibus_RSA/results.mat and omnibus_RSA_TU/results.mat. ' ...
    'Original five ROIs retain subjects 2–6; TU uses subjects 2, 3, 4, and 6. ' ...
    'Subject 5 TU is missing (6 usable voxels, minimum 10), never zero-filled.\n\n' ...
    'Plots display absolute standardized subject betas. Stars retain saved signed-beta BH-FDR q-values: ' ...
    '20 original ROI × predictor tests and a separate four-predictor TU family. ' ...
    'No RSA or permutations were rerun and original result MAT files/tables were not modified.\n\n' ...
    'Reproduce with scripts/update_omnibus_RSA_figures_with_TU.m. ' ...
    'combined_with_TU_plot_data.mat records the exact plotted values and q-values.\n']);
end
