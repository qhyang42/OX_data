function plot_rating_displacement_RSA_formal()
% Formal display of saved overall context-pair-adjusted coefficients.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'utils','OX_utilities'));
source=fullfile(root,'RDMs','rating_displacement_RSA','results.mat');
L=load(source,'results'); P=L.results;
keep=arrayfun(@(m) string(m.name)=="context_pair",P.models);
assert(nnz(keep)==1); P.models=P.models(keep);
assert(isequal(P.subject_ids,2:6));
assert(isequal(size(P.models.beta),[5 5 2]) && all(isfinite(P.models.beta),'all'));
G=P.models.group_statistics;
assert(height(G)==10 && all(G.n_subjects==5) && all(G.n_permutations==5000));
for r=1:5
    for j=1:2
        idx=G.roi==P.roi_names(r) & G.predictor==P.models.predictors(j);
        assert(nnz(idx)==1 && abs(mean(P.models.beta(:,r,j))-G.mean_beta(idx))<1e-12);
    end
end
out=fullfile(root,'results','rating_displacement_RSA');
assert(~isfolder(out),'Output exists; preserve previous figure.'); mkdir(out);
P.output_dir=out; P.figure_output_dir=out;
P.figure_note="Group permutation test · BH-FDR across 10 coefficients | * q < .05, ** q < .01, *** q < .001";
OX_plot_formal_rsa(P);
% Store only display inputs; complete nulls/provenance remain in source.
subject_ids=P.subject_ids; roi_names=P.roi_names; model=P.models;
save(fullfile(out,'plot_data.mat'),'subject_ids','roi_names','model','source','-v7.3');
writetable(G,fullfile(out,'group_statistics.csv'));
copyfile(fullfile(root,'RDMs','rating_displacement_RSA','tables','context_pair_subject_statistics.csv'),fullfile(out,'subject_statistics.csv'));
fid=fopen(fullfile(out,'README.md'),'w'); assert(fid>=0);
fprintf(fid,['# Rating displacement RSA: context-pair adjustment\n\n' ...
 'Overall same-odor cross-context displacement coefficients from the saved context_pair model. ' ...
 'Both absolute pleasantness and intensity changes enter simultaneously, with context-pair fixed effects and no odor fixed effects. ' ...
 'All six pairs among PERSON, FOOD, LOCATION, CONTROL contribute (20 odors; 120 observations per subject/ROI). ' ...
 'Outcome (1 - Pearson r) and rating predictors are sample-SD standardized. Signed betas are displayed.\n\n' ...
 'Boxes summarize five subjects (2–6): median and interquartile range, whiskers to the most extreme values within 1.5 IQR. ' ...
 'All five colored subject dots are shown, with fixed offsets by subject. ' ...
 'Stars use saved two-sided inclusive add-one permutation p-values and BH-FDR across five ROIs × two predictors: ' ...
 '* q < .05, ** q < .01, *** q < .001; unmarked coefficients do not pass q < .05. ' ...
 'Inference concerns the measured subjects: 5,000 global neural condition-correspondence permutations per subject, ' ...
 'with equal-weight subject null averaging. These are not coefficient-specific conditional nulls. ' ...
 'See RDMs/rating_displacement_RSA/README.md and results.mat for complete methods and input provenance.\n\n' ...
 'No fits or permutations were rerun. This requested context-pair-only model is the documented sensitivity model. ' ...
 'Reproduce with scripts/plot_rating_displacement_RSA_formal.m. PNG and vector PDF share the stem coefficients_context_pair.\n']);
fclose(fid);
end
