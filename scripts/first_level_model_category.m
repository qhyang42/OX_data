%% run first level analysis in SPM 

%% 
SUBJNAMES = {'240711_fMRI_OX_NWU_AS', ...
        '240723_fMRI_OX_NWU_LS', ...
        '240814_fMRI_OX_NWU_JN', ...
        '240816_fMRI_OX_NWU_RR', ...
        '241018_fMRI_OX_NWU_BN', ...
        '250117_fMRI_OX_NWU_VS'}; 

session_count = [4, 13, 10, 15, 16, 10]; % number of sessions for each subj so far. EDIT as needed. 
TR= 0.76; 

subjidx = 2; % enter subjidx here 
subjname = ['subj_', num2str(subjidx)];  
subjname_real = SUBJNAMES{subjidx}; 

% mridir = '/Volumes/ExtremeSSD/OX_DATA/MRI'; 
mridir = '/Users/qhyang/Desktop/OX_DATA/MRI'; 

mridatapath = fullfile(mridir, ['subj_', num2str(subjidx)], 'nifti'); % sn:subject's name

% evdir = '/Volumes/ExtremeSSD/OX_DATA/labchart'; 
evdir = '/Users/qhyang/Desktop/OX_DATA/labchart'; 


%% make multiple regressors -- this only needs to be done for each subjct once 
R = make_motion_regressor (subjname_real, mridatapath, 1, session_count(subjidx)); 
%%
if ~exist(fullfile(mridatapath, 'motion_param'), 'file')
    save(fullfile(mridatapath, 'motion_param'), "R");
else 
    fprintf('motion parameter file already exists!')
end

%% specify
% Initialize matlabbatch
matlabbatch = [];

% Specify the directory to save the first-level model
matlabbatch{1}.spm.stats.fmri_spec.dir = {mridatapath};  % Change to your desired directory

% Specify the timing parameters
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';  % 'secs' or 'scans' for onset times
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;  % Repetition time (TR) in seconds
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 30;  % Number of slices per volume (for temporal resampling)
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;  % Reference slice (usually the first slice)

%% Specify the preprocessed functional data (e.g., smoothed and realigned images) -- 4D volume does not work 
% filename = func_list(subjname_real, mridatapath, 1, session_count(subjidx)); 
% %%% concatenate 4D files 
% if ~exist(fullfile(mridatapath, 'func', 'merged_runs.nii'), 'file')
%     spm_file_merge(filename, fullfile(mridatapath, 'func', 'merged_runs.nii'));
% end

% %% try merging with fslmerge 
% [~, filename_fsl] = func_list(subjname_real, mridatapath, 1, session_count(subjidx)); 
% fsloutput = fullfile(mridatapath, 'func' , 'merged_runs_fsl.nii');
% fslmerge_command = sprintf('fslmerge -t %s %s', fsloutput, filename_fsl);
% %% 
% system(fslmerge_command);

%% specify functional data as 3d images 
filename = func_lsit_3D(subjname_real, mridatapath, 1, session_count(subjidx)); 


%% load events 
load(fullfile(evdir, [subjname, '_events.mat'])); % all events 
[odor, category] = OX_get_odor(subjname); 

% goodtrials = find(goodtrials); % this is specific to cue. odor is always fine. 
personi = find(strcmp(category, 'PERSON'));
foodi = find(strcmp(category, 'FOOD'));
loci = find(strcmp(category, 'LOCATION')); 
controli = find(strcmp(category, 'CONTROL')); 

% categoryidx = zeros(length(category), 1); 
% categoryidx(personi) = 1; 
% categoryidx(foodi) = 2;
% categoryidx(loci) = 3; 
% categoryidx(controli) = 4; 
% 
% catvec = reshape(categoryidx, [numel(categoryidx), 1]); 
% odorvec = reshape(odor, [numel(categoryidx), 1]); 
%% 
% matlabbatch{1}.spm.stats.fmri_spec.sess.scans = {fullfile(mridatapath, 'func', 'merged_runs.nii')};  % List of preprocessed NIfTI files for the run/session
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = filename; 
% Specify the conditions (onsets and durations)

%% odor contrast 
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'Person';  % Name of the condition
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset = event_onsets_vec(personi);    % Onsets in seconds
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = 2;    % Durations in seconds (or 0 for events)
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;               % Temporal modulation (0 = none)

matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'Food';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset = event_onsets_vec(foodi);    % Onsets in seconds
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = 2; % use cue word onset as cue event    
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;  

matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name = 'Location';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset = event_onsets_vec(loci);    % Onsets in seconds
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).duration = 2; % use cue word onset as cue event    
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;  

matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name = 'Control';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset = event_onsets_vec(controli);    % Onsets in seconds
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration = 2; % use cue word onset as cue event    
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;  

%% cue contrast 
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'Person';  % Name of the condition
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset = cue_onsets_all(personi);    % Onsets in seconds
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = 5;    % Durations in seconds (or 0 for events)
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;               % Temporal modulation (0 = none)

matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'Food';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset = cue_onsets_all(foodi);    % Onsets in seconds
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = 5; % use cue word onset as cue event    
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;  

matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name = 'Location';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset = cue_onsets_all(loci);    % Onsets in seconds
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).duration = 5; % use cue word onset as cue event    
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;  

matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name = 'Control';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset = cue_onsets_all(controli);    % Onsets in seconds
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration = 5; % use cue word onset as cue event    
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;  
%%

matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128; 

% regressors 
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {fullfile(mridatapath, 'motion_param.mat')}; 

%%  
matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(mridatapath, 'SPM.mat')};

%% specify
spm_jobman('run', matlabbatch(1));
%% add run specific regressors 
spm_fmri_concatenate(fullfile(mridatapath, 'SPM.mat'), nframes');

%% estimate 
spm_jobman('run', matlabbatch(2));

%% contrast
% load(fullfile(mridatapath, 'SPM.mat'));

% % Specify a t-contrast: Odor > Rest
% matlabbatch{3}.spm.stats.con.spmmat = {fullfile(mridatapath, 'SPM.mat')};
% matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Odor > Rest';
% matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 0];
% matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
% 
% % Specify another t-contrast: Cue > Rest
% matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Cue > Rest';
% matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 1];
% matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';


% matlabbatch{3}.spm.stats.con.spmmat = {fullfile(mridatapath, 'SPM.mat')};
% 
% matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'person > rest';
% matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1, 0, 0, 0];
% matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
% 
% matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'food > rest';
% matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0, 1, 0, 0];
% matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
% 
% matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'location > rest';
% matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0, 0, 1, 0];
% matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
% 
% matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'control > rest';
% matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0, 0, 0, 1];
% matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
% 
% %% 
% spm_jobman('run', matlabbatch(3));
% 
% %% contrast with control. 
% matlabbatch{4}.spm.stats.con.spmmat = {fullfile(mridatapath, 'SPM.mat')};
% 
% matlabbatch{4}.spm.stats.con.consess{1}.tcon.name = 'person > control';
% matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights = [1, 0, 0, -1];
% matlabbatch{4}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
% 
% matlabbatch{4}.spm.stats.con.consess{2}.tcon.name = 'food > control';
% matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [0, 1, 0, -1];
% matlabbatch{4}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
% 
% matlabbatch{4}.spm.stats.con.consess{3}.tcon.name = 'location > control';
% matlabbatch{4}.spm.stats.con.consess{3}.tcon.weights = [0, 0, 1, -1];
% matlabbatch{4}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
% %% 
% spm_jobman('run', matlabbatch(4));

%% contrast with control. 
matlabbatch{5}.spm.stats.con.spmmat = {fullfile(mridatapath, 'SPM.mat')};

matlabbatch{5}.spm.stats.con.consess{1}.tcon.name = 'person > else';
matlabbatch{5}.spm.stats.con.consess{1}.tcon.weights = [3, -1, -1, -1];
matlabbatch{5}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

matlabbatch{5}.spm.stats.con.consess{2}.tcon.name = 'food > else';
matlabbatch{5}.spm.stats.con.consess{2}.tcon.weights = [-1, 3, -1, -1];
matlabbatch{5}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

matlabbatch{5}.spm.stats.con.consess{3}.tcon.name = 'location > else';
matlabbatch{5}.spm.stats.con.consess{3}.tcon.weights = [-1, -1, 3, -1];
matlabbatch{5}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
%% 
spm_jobman('run', matlabbatch(5));


%% calculate fwe corrected T threshold 
load(fullfile(mridatapath, 'SPM.mat'));

df = [SPM.xX.erdf SPM.xX.trRV];  % Error df and residual variance

% Compute FWE-corrected voxel-level threshold
p_fwe = 0.001;  % Desired FWE-corrected significance level
STAT = 'T';    % Statistic type ('T' for T-maps)
R = SPM.xVol.R;  % Resels (spatial smoothness of the data)
n_voxels = prod(SPM.xVol.DIM);  % Number of voxels
u = spm_uc(p_fwe, df, STAT, R, 1, n_voxels);  % FWE threshold


%% contrast by category 

%% test register Tmap to standard space 
% cimages = {fullfile(mridatapath, 'spmT_0001.nii'); ...
%     fullfile(mridatapath, 'spmT_0002.nii')}; 
% anatfile = dir(fullfile(mridatapath, 'anat', 'mean*.nii')); 
% anatfile = anatfile.name; 
% fname = fullfile(mridatapath, 'anat', sprintf('y_%s', anatfile));
% matlabbatch{4}.spm.spatial.normalise.write.subj.def = {fname};
% matlabbatch{4}.spm.spatial.normalise.write.subj.resample = cimages;
% matlabbatch{4}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
% matlabbatch{4}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
% matlabbatch{4}.spm.spatial.normalise.write.woptions.interp = 4;

%% steps after T contrast 
% brain extraction using bet2 
% func to wb to T1 brain coreg 
% register contrast to T1
% threshold contrast 



%% functions 
%%
function [filename, filename_fsl] = func_list(subjname, datapath, nsess_i, nsess_f)

%%%%% NOTE. This is now working for 4D nifti volumes. 

% rnctr is run counter
rcntr = 0;
filename = [];
filename_fsl = []; 

for sess = nsess_i:nsess_f

    path_ = fullfile(datapath, 'func/');
    n = dir(fullfile(path_, [subjname, '_', num2str(sess), '*.nii']));
    nfiles = length(n); 

    %     if rcntr==1
    %         funcpath = path_;
    %         funcfile = n(1).name;
    %     end

    for i=1:length(n)
        rcntr = rcntr+1; % sess_2_run_1 is counted as 5 if sess_1 had 4 runs
        thisfile = dir(fullfile(path_, ['sr', subjname, '_', num2str(sess), '_Run', num2str(i), '_*.nii']));
        fname = fullfile(path_, thisfile.name);
        % Different runs in different cells, add ",1" for spm
%         filename{rcntr}{i,1} = sprintf('%s,1', fname);
%         filename{rcntr, 1}{1,1} = sprintf('%s', fname);
        filename{rcntr, 1} = sprintf('%s', fname);

        filename_fsl = [filename_fsl, ' ', fname]; 
    end

end


end

%% 
function filename = func_lsit_3D(subjname, datapath, nsess_i, nsess_f)
%%% list all functional volumes in one cell array 

vcntr = 0; % run counter
filename = [];

for sess = nsess_i:nsess_f
    % sess_2_run_1 is counted as 5 if sess_1 had 4 runs
    path_ = fullfile(datapath, 'func/');
    n = dir(fullfile(path_, ['sr', subjname, '_', num2str(sess), '_*.nii']));
    nfiles = length(n);

    for i=1:nfiles

        thisfile = dir(fullfile(path_, ['sr', subjname, '_', num2str(sess), '_Run', num2str(i), '_*.nii']));
        fname = fullfile(path_, thisfile.name);

        nifti_header = niftiinfo(fname);
        nvols = nifti_header.ImageSize(4);

        for vol = 1:nvols
            vcntr = vcntr+1;
            % list all volumes in a cell array 
            filename{vcntr, 1} = sprintf('%s,%d', fname, vol);
        end
    end

end
end

%%
function R = make_motion_regressor (subjname, datapath, nsess_i, nsess_f)
%%% make motion regressors
motion_regressor = [];
path_ = fullfile(datapath, 'func/');

for sess = nsess_i:nsess_f

%     path_ = fullfile(datapath, 'func/');
    n = dir(fullfile(path_, [subjname, '_', num2str(sess), '_*.nii']));
    nfiles = length(n);

    for i = 1:nfiles
        mpfile = dir(fullfile(path_, ['rp_', subjname, '_', num2str(sess), '_Run', num2str(i), '_*.txt']));
        mp = load(fullfile(path_, mpfile.name));
        current_reg = [mp, [zeros(1,6); diff(mp)], mp.^2, [zeros(1,6); diff(mp).^2]];
        motion_regressor = [motion_regressor; current_reg];
    end
end
R = motion_regressor; 
end


