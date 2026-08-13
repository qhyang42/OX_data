function info = setup_ox()
%SETUP_OX Configure MATLAB for the OX_DATA project.
%   SETUP_OX adds the project code, SPM, FSL, and FreeSurfer MATLAB paths.
%   It also configures the environment inherited by commands launched with
%   SYSTEM. Existing environment variables override the machine defaults.
%
%   INFO = SETUP_OX returns the resolved paths.

setupFile = mfilename('fullpath');
environmentDir = fileparts(setupFile);
projectRoot = fileparts(environmentDir);

matlabRoot = envOrDefault('MATLAB_ROOT', '/Applications/MATLAB_R2025b.app');
fslDir = envOrDefault('FSLDIR', '/Users/qhyang/fsl');
freesurferHome = envOrDefault('FREESURFER_HOME', '/Applications/freesurfer');
subjectsDir = envOrDefault('SUBJECTS_DIR', fullfile(freesurferHome, 'subjects'));
fsLicense = envOrDefault('FS_LICENSE', fullfile(freesurferHome, 'license.txt'));
spmDir = envOrDefault('SPM_DIR', '/Users/qhyang/Desktop/Utilities/spm');
glmSingleDir = envOrDefault('GLMSINGLE_DIR', ...
    '/Users/qhyang/Desktop/Utilities/GLMsingle');

requiredDirectories = {
    projectRoot, ...
    fullfile(projectRoot, 'scripts'), ...
    fullfile(projectRoot, 'utils', 'OX_utilities'), ...
    matlabRoot, ...
    fslDir, ...
    freesurferHome, ...
    spmDir, ...
    glmSingleDir
};
for directoryIndex = 1:numel(requiredDirectories)
    assertDirectory(requiredDirectories{directoryIndex});
end
if ~isfile(fsLicense)
    error('OX_DATA:MissingFreeSurferLicense', ...
        'FreeSurfer license not found: %s', fsLicense);
end
glmSingleSetup = fullfile(glmSingleDir, 'setup.m');
if ~isfile(glmSingleSetup)
    error('OX_DATA:MissingGLMsingleSetup', ...
        'GLMsingle setup script not found: %s', glmSingleSetup);
end

setenv('MATLAB_ROOT', matlabRoot);
setenv('MATLAB_BIN', fullfile(matlabRoot, 'bin', 'matlab'));
setenv('FSLDIR', fslDir);
setenv('FREESURFER_HOME', freesurferHome);
setenv('SUBJECTS_DIR', subjectsDir);
setenv('FS_LICENSE', fsLicense);
setenv('SPM_DIR', spmDir);
setenv('GLMSINGLE_DIR', glmSingleDir);
if isempty(getenv('FSLOUTPUTTYPE'))
    setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
end
if isempty(getenv('FSLMULTIFILEQUIT'))
    setenv('FSLMULTIFILEQUIT', 'TRUE');
end

% Make neuroimaging commands available to MATLAB child processes. The FSL
% wrapper directory is preferred over FSLDIR/bin to avoid exposing FSL's
% bundled copies of unrelated programs.
commandDirectories = {
    fullfile(matlabRoot, 'bin'), ...
    fullfile(freesurferHome, 'bin'), ...
    fullfile(freesurferHome, 'fsfast', 'bin'), ...
    fullfile(freesurferHome, 'mni', 'bin'), ...
    fullfile(fslDir, 'share', 'fsl', 'bin')
};
prependEnvironmentPath(commandDirectories);

% Add only the intended entry points for external toolboxes. In particular,
% do not use genpath on SPM because it contains compatibility copies and
% private implementation directories that should not all be on MATLAB path.
addpath(fullfile(projectRoot, 'scripts'));
addpath(fullfile(projectRoot, 'utils'));
addpath(fullfile(projectRoot, 'utils', 'OX_utilities'));
addpath(fullfile(fslDir, 'etc', 'matlab'));
addpath(fullfile(freesurferHome, 'matlab'));
addpath(spmDir);

% Let GLMsingle configure its own paths and validate its fracridge
% dependency. Running the project's entry point keeps this setup aligned
% with future changes to the installed GLMsingle checkout.
run(glmSingleSetup);

info = struct( ...
    'projectRoot', projectRoot, ...
    'matlabRoot', matlabRoot, ...
    'fslDir', fslDir, ...
    'freesurferHome', freesurferHome, ...
    'subjectsDir', subjectsDir, ...
    'fsLicense', fsLicense, ...
    'spmDir', spmDir, ...
    'glmSingleDir', glmSingleDir);

if nargout == 0
    fprintf('OX_DATA MATLAB environment configured.\n');
    fprintf('  Project:    %s\n', projectRoot);
    fprintf('  MATLAB:     %s\n', version);
    fprintf('  SPM:        %s\n', which('spm'));
    fprintf('  GLMsingle:  %s\n', which('GLMestimatesingletrial'));
    fprintf('  FSLDIR:     %s\n', fslDir);
    fprintf('  FreeSurfer: %s\n', freesurferHome);
    clear info
end
end

function value = envOrDefault(variableName, defaultValue)
value = getenv(variableName);
if isempty(value)
    value = defaultValue;
end
end

function assertDirectory(directoryPath)
if ~isfolder(directoryPath)
    error('OX_DATA:MissingDirectory', ...
        'Required directory not found: %s', directoryPath);
end
end

function prependEnvironmentPath(directories)
currentPath = strsplit(getenv('PATH'), pathsep);
for directoryIndex = numel(directories):-1:1
    directoryPath = directories{directoryIndex};
    if isfolder(directoryPath)
        currentPath(strcmp(currentPath, directoryPath)) = [];
        currentPath = [{directoryPath}, currentPath]; %#ok<AGROW>
    end
end
setenv('PATH', strjoin(currentPath, pathsep));
end
