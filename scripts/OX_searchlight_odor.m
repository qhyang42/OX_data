function results = OX_searchlight_odor(subjidx, varargin)
%OX_SEARCHLIGHT_ODOR Run odor-identity searchlight decoding.
%
%   results = OX_searchlight_odor(subjidx)
%   results = OX_searchlight_odor(subjidx, Name, Value, ...)
%
% Odor identity is decoded across all context conditions. The decoder uses
% leave-one-run-out cross-validation and a nearest-template Pearson
% correlation classifier. Statistical inference uses label permutations
% within run, voxelwise empirical p-values, Benjamini-Hochberg FDR q-values,
% and max-statistic FWE-corrected p-values.
%
% Fast smoke test (sparse centers; do not use as a final map):
%   OX_searchlight_odor(5, 'MaxCenters', 100, 'NumPermutations', 2, 'UseParallel', false)
%
% Common test call:
%   OX_searchlight_odor(5, 'NumPermutations', 10, 'UseParallel', false)
%
% Final analysis example:
%   OX_searchlight_odor(5, 'NumPermutations', 1000)

results = OX_searchlight_decode_core(subjidx, 'odor', varargin{:});
end
