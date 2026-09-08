function results = run_rating_displacement_rsa(varargin)
%RUN_RATING_DISPLACEMENT_RSA Permutation coefficients and pooled mixed models.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'utils','OX_utilities'));
results = OX_formal_rsa('displacement',root,varargin{:});
end
