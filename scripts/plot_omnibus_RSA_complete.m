function plot_omnibus_RSA_complete(signed_display)
% Pass true to save the signed-beta display in results/omnibus_RSA_complete/signed.
if nargin<1, signed_display=false; end
% Combine saved omnibus fits for display, retaining their inferential families.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'utils','OX_utilities'));
sources=["omnibus_RSA","omnibus_RSA_TU","omnibus_RSA_control"];
S=cell(1,3);
for k=1:3
 L=load(fullfile(root,'RDMs',sources(k),'results.mat'),'results');S{k}=L.results;
end
A=S{1};T=S{2};C=S{3};
assert(isequal(A.subject_ids,2:6)&&isequal(T.subject_ids,[2 3 4 6])&&isequal(C.subject_ids,2:6));
assert(isequal(A.models.predictors,T.models.predictors,C.models.predictors));
assert(isequal(T.roi_names,"olf_TU"));
assert(isequal(C.roi_names,["control_GM_putamen","control_A1"]));
out=fullfile(root,'results','omnibus_RSA_complete');
if signed_display, out=fullfile(root,'results','omnibus_RSA_complete','signed'); end
assert(~isfolder(out),'Output exists; preserve previous figure.');mkdir(out);
P=struct('analysis',"omnibus",'subject_ids',2:6, ...
 'roi_names',[T.roi_names,A.roi_names,C.roi_names], ...
 'roi_display_names',["TU",A.roi_names,"Putamen","A1"], ...
 'control_roi_names',C.roi_names,'output_dir',out,'figure_output_dir',out,'plot_by_roi',false);
P.signed_beta_display=signed_display;
P.models=struct('name',"omnibus",'predictors',A.models.predictors, ...
 'beta',nan(5,8,4),'group_statistics',[T.models.group_statistics;A.models.group_statistics;C.models.group_statistics]);
[found,rows]=ismember(T.subject_ids,P.subject_ids);assert(all(found));
P.models.beta(rows,1,:)=T.models.beta;
P.models.beta(:,2:6,:)=A.models.beta;
P.models.beta(:,7:8,:)=C.models.beta;
assert(all(isnan(P.models.beta(4,1,:)),'all'));
assert(isequaln(P.models.beta(:,2:6,:),A.models.beta));
assert(isequaln(P.models.beta(:,7:8,:),C.models.beta));
P.figure_note=["Gray boxes/dots: controls (GM only); colored ROIs: GM + Odor > Rest p < .001"; ...
 "Signed-beta FDR families: olfactory 20 tests; TU 4; controls 8 | * q < .05, ** q < .01, *** q < .001"];
OX_plot_formal_rsa(P);
save(fullfile(out,'plot_data.mat'),'P','sources','-v7.3');
fid=fopen(fullfile(out,'README.md'),'w');assert(fid>=0);
display_description='absolute standardized subject betas, as in the existing omnibus figures';
reproduce='scripts/plot_omnibus_RSA_complete.m';
if signed_display
 display_description=['signed standardized subject betas. Categorical predictors are similarity-coded (1 = same); ' ...
  'negative odor/context betas indicate smaller neural distances for same-identity pairs after adjustment'];
 reproduce='addpath scripts; plot_omnibus_RSA_complete(true)';
end
fprintf(fid,['# Complete omnibus RSA by predictor\n\n' ...
 'Saved results combined from RDMs/omnibus_RSA, omnibus_RSA_TU and omnibus_RSA_control; no analyses rerun. ' ...
 'Order: TU, AON, PirF, PirT, olfAMG, olfOFC, Putamen, A1. Controls use gray boxes and gray subject dots. ' ...
 'Dots retain fixed subject offsets (subjects 2–6, left to right). Colored dots identify subjects for other ROIs.\n\n' ...
 'Displays %s. Stars use saved signed-beta group tests ' ...
 'and retain separate BH-FDR families: 20 original olfactory tests, four TU tests, eight control tests. No combined 32-test correction was recomputed. ' ...
 'Original/TU ROIs use GM + Odor > Rest uncorrected p < .001; controls use GM only. ' ...
 'TU omits subject 5 (<10 usable voxels), represented by NaN rather than zero; all other ROIs use subjects 2–6.\n\n' ...
 'Files: betas_by_predictor.png, vector betas_by_predictor.pdf, plot_data.mat. ' ...
 'Reproduce with `%s`.\n'],display_description,reproduce);fclose(fid);
end
