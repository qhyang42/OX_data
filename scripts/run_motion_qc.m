%% Run read-only motion/artifact QC for the completed OX dataset.
% This stage calculates FD/DVARS and diagnostic censor-rule comparisons.
% It does not modify preprocessing, GLMsingle, nuisance regressors, or GLMs.

setup_ox;
for subjidx = 2:6
    OX_compute_motion_qc(subjidx);
end
OX_summarize_motion_qc('Subjects', 2:6);
