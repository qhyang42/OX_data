function results = run_rating_displacement_mixed_models()
% Update only mixed models using the saved formal RSA observations/results.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'utils','OX_utilities'));
out=fullfile(root,'RDMs','rating_displacement_RSA');
L=load(fullfile(out,'results.mat'),'results'); results=L.results;
if isfield(results,'mixed_models') && ~isfield(results.mixed_models,'grouping') && ...
        ~isfolder(fullfile(out,'archive','pair_specific_mixed_models'))
    archive=fullfile(out,'archive','pair_specific_mixed_models');
    mkdir(archive);
    for folder=["tables","figures"]
        mkdir(fullfile(archive,folder));
        if folder=="tables", patterns=["mixed_model_*.csv","mixed_model_*.md"];
        else, patterns=["displacement_vs_*_mixed_model.*","mixed_model_diagnostics_*.*"]; end
        for pattern=patterns
            files=dir(fullfile(out,folder,pattern));
            for k=1:numel(files), copyfile(fullfile(files(k).folder,files(k).name),fullfile(archive,folder,files(k).name)); end
        end
    end
    copyfile(fullfile(out,'mixed_models.mat'),fullfile(archive,'mixed_models.mat'));
end
previous_models=results.models;
results.mixed_models=OX_rsa_context_displacement_lme(results.observations,out,true);
assert(isequaln(previous_models,results.models),'Permutation coefficient results changed unexpectedly.');
save(fullfile(out,'results.mat'),'results','-v7.3'); OX_rsa_write_readme(results);
% Remove superseded active files only after the new results are complete and archived.
if isfolder(fullfile(out,'archive','pair_specific_mixed_models'))
    stale=["tables/mixed_model_pair_slopes.csv","tables/mixed_model_pair_slopes.md", ...
        "figures/displacement_vs_pleasantness_mixed_model.png","figures/displacement_vs_pleasantness_mixed_model.pdf", ...
        "figures/displacement_vs_intensity_mixed_model.png","figures/displacement_vs_intensity_mixed_model.pdf", ...
        "figures/mixed_model_diagnostics_pleasantness.png","figures/mixed_model_diagnostics_pleasantness.pdf", ...
        "figures/mixed_model_diagnostics_intensity.png","figures/mixed_model_diagnostics_intensity.pdf"];
    for file=stale, if isfile(fullfile(out,file)), delete(fullfile(out,file)); end, end
end
end
