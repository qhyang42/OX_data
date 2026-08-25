function [multipliers, flippedRuns, reason] = OX_get_respiration_polarity(subjidx, nruns)
%OX_GET_RESPIRATION_POLARITY Return the validated OX airflow polarity map.
%   MULTIPLIERS = OX_GET_RESPIRATION_POLARITY(SUBJIDX, NRUNS) returns one
%   multiplier per numerically ordered run. Positive airflow represents
%   inhale after applying the multiplier.

%   This is the single project definition used by both BreathMetrics onset
%   processing and the GLMsingle confound builder.

%   Subject 3 session 1 (runs 1-10) was acquired with reversed airflow
%   polarity, confirmed from the sniff-locked waveform.

%   See also OX_BUILD_GLMSINGLE_CONFOUNDS.

validateattributes(subjidx, {'numeric'}, {'scalar', 'integer', 'finite'});
validateattributes(nruns, {'numeric'}, {'scalar', 'integer', 'positive', 'finite'});

flippedRuns = [];
reason = '';
if subjidx == 3
    flippedRuns = 1:10;
    reason = 'Reversed acquisition polarity; confirmed from sniff-locked waveform';
end
if any(flippedRuns > nruns)
    error('OX:Respiration:PolarityRunOutOfRange', ...
        'Configured polarity-flip run exceeds the %d runs for subject %d.', ...
        nruns, subjidx);
end

multipliers = ones(nruns, 1);
multipliers(flippedRuns) = -1;
end
