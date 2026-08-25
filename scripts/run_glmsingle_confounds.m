%% Build final run-wise GLMsingle nuisance regressors for OX subjects 2-6.
% This does not launch GLMsingle and does not modify existing derivatives.
% A complete source preflight is required before any subject is built.

setup_ox;
subjects = 2:6;
preflight = OX_preflight_glmsingle_confounds('Subjects', subjects);
problems = preflight(preflight.status ~= "ok", :);
if ~isempty(problems)
    disp(problems);
    error('OX:GLMsingleConfounds:PreflightFailed', ...
        ['GLMsingle confound preflight found %d run(s) that cannot be ', ...
         'built without inventing respiratory samples. No subjects were run.'], ...
        height(problems));
end

for subjidx = subjects
    OX_build_glmsingle_confounds(subjidx);
end
OX_summarize_glmsingle_confounds('Subjects', subjects);
