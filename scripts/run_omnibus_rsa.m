function results = run_omnibus_rsa(varargin)
%RUN_OMNIBUS_RSA Formal simple-distance, linear omnibus RSA.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'utils','OX_utilities'));
results = OX_formal_rsa('omnibus',root,varargin{:});
end
