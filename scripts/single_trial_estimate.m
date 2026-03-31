%% estimate single trial beta values with glmsingle 

%% dont do this. this does not distinguish categories 
SUBJNAMES = {'240711_fMRI_OX_NWU_AS', ...
        '240723_fMRI_OX_NWU_LS', ...
        '240814_fMRI_OX_NWU_JN', ...
        '240816_fMRI_OX_NWU_RR', ...
        '241018_fMRI_OX_NWU_BN', ...
        '250117_fMRI_OX_NWU_VS'}; 

session_count = [4, 13, 10, 15, 16, 10]; % number of sessions for each subj so far. EDIT as needed. 
TR= 0.76; 

subjidx = 5; % enter subjidx here 
subjname = ['subj_', num2str(subjidx)];  
subjname_real = SUBJNAMES{subjidx}; 

% mridir = '/Volumes/ExtremeSSD/OX_DATA/MRI'; 
mridir = '/Users/qhyang/Desktop/OX_DATA/MRI'; 

mridatapath = fullfile(mridir, ['subj_', num2str(subjidx)], 'nifti'); % sn:subject's name

% evdir = '/Volumes/ExtremeSSD/OX_DATA/labchart'; 
evdir = '/Users/qhyang/Desktop/OX_DATA/labchart'; 

nruns = 80; 
% stimdur = 2;
stimdur = 0; 

outdir = fullfile(mridatapath, 'single_trials'); 

%% get design matrix 
%%% design matrix should be stim type * TR. 
%%% we want a cue beta and an odor beta for each voxel to look at patterns 
load(fullfile(evdir, [subjname, '_events.mat'])); % all events 

designmat = zeros(sum(nframes), 2); % TR by stimuli 

% get n of TR from event_onset and cue_onset 
eventTR = ceil(event_onsets_vec/TR); 
cueTR = ceil(cue_onsets_vec/TR); 

designmat(eventTR, 1) = 1;
designmat(cueTR, 2) = 1; 

%% glm single only takes one stim duration
designmat = designmat(:, 1); % odor
% designmat = designmat(:, 2); % cue

%% get data. use a functional mask/grey matter mask
%%% use a smaller mask for proof of concept 

mask = spm_vol(fullfile(mridatapath, 'coreg', 'gm_mask_thr05_func.nii')); 
maskvol = spm_read_vols(mask); 
maskvol = logical(maskvol);
%%% read func data
funcfiles = func_list(subjname_real, mridatapath, 1, session_count(subjidx));

% Get data from all voxels in the mask -- from vivek's sample script
% ARC_createsingletrials.m 
Res_V = spm_vol(funcfiles);
[~, XYZmm] = spm_read_vols(Res_V{1});
XYZvx = round(Res_V{1}(1).mat\[XYZmm; ones(1,size(XYZmm,2))]);
Res_Vol = cell2mat(Res_V);
mask1D = maskvol(:);
vxl = XYZvx(:,mask1D);
voxel_act = spm_get_data(Res_Vol,vxl)';
voxel_act = single(voxel_act);
[r,c] = find(isnan(voxel_act));
voxel_act(r,:) = []; % this is the input data for glmsingle
% t_axis = (0:1:sum(nframes))*TR;

%% get nuissance regressors 
load(fullfile(mridatapath, 'motion_param.mat')); 

%% put data into run-wise blocks 
% Change to runwise blocks
% can't seem to run without this. Let vivek know if i figure out why 
design_m2 = cell(1,nruns);
voxel_act2 = cell(1,nruns);
mp2 = cell(1,nruns);
for zz = 1:nruns
    if zz==1
        idx = 1:nframes(1);
    else 
        idx = sum(nframes(1:zz-1))+1:sum(nframes(1:zz)); 
    end 

    design_m2{zz} = designmat(idx,:);
    voxel_act2{zz} = voxel_act(:,idx);
    temp = R(idx,:);
    mp2{zz} = temp;%(:,any(temp)); % Include non-zero regressors
end
opt.extraregressors = mp2;



%% run glmsingle 

% try
    results = GLMestimatesingletrial(design_m2,voxel_act2,stimdur,TR,outdir,opt);
% catch
%     save(fullfile(outdir,'error_rep.mat'));
% end



%% try plot beta across voxels 

% clear; 
% SUBJNAMES = {'240711_fMRI_OX_NWU_AS', ...
%         '240723_fMRI_OX_NWU_LS', ...
%         '240814_fMRI_OX_NWU_JN', ...
%         '240816_fMRI_OX_NWU_RR', ...
%         '240816_fMRI_OX_NWU_BN'}; 
% 
% session_count = [4, 13, 10, 8, 6]; % number of sessions for each subj so far. EDIT as needed. 
% TR= 0.76; 
% 
% subjidx = 3; % enter subjidx here 
% subjname = ['subj_', num2str(subjidx)];  
% subjname_real = SUBJNAMES{subjidx}; 
% 
% % mridir = '/Volumes/ExtremeSSD/OX_DATA/MRI'; 
% mridir = '/Users/qhyang/Desktop/OX_DATA/MRI'; 
% 
% mridatapath = fullfile(mridir, ['subj_', num2str(subjidx)], 'nifti'); % sn:subject's name
% 
% % evdir = '/Volumes/ExtremeSSD/OX_DATA/labchart'; 
% evdir = '/Users/qhyang/Desktop/OX_DATA/labchart'; 
% 
% nruns = 80; 
% stimdur = 2;
% 
% outdir = fullfile(mridatapath, 'single_trials'); 
% 
% load(fullfile(outdir, 'TYPED_FITHRF_GLMDENOISE_RR.mat')); 
% [odor, category] = OX_get_odor(subjname); 
% 
% personi = find(strcmp(category, 'PERSON'));
% foodi = find(strcmp(category, 'FOOD'));
% loci = find(strcmp(category, 'LOCATION')); 
% controli = find(strcmp(category, 'CONTROL')); 
% 
% modelmd = squeeze(modelmd);
% bperson = modelmd(:, personi);
% bfood = modelmd(:, foodi);
% bloc = modelmd(:, loci);
% bcontrol = modelmd(:, controli);
% 
% %% 
% line_colors = [[0, 0.4470, 0.7410];...	          
%           	[0.8500, 0.3250, 0.0980];...	       
%           	[0.9290, 0.6940, 0.1250];...	          
%           	[0.4940, 0.1840, 0.5560];...	          
%           	[0.4660, 0.6740, 0.1880];...	          
%           	[0.3010, 0.7450, 0.9330];...	          
%           	[0.6350, 0.0780, 0.1840]]; 	
% 
% v = 1:1:size(modelmd,1); 
% 
% figure; 
% hold on
% [~, lobj(1)] = patchPlot(bperson', v, line_colors(1, :), '-', 1, 2);
% [~, lobj(2)] = patchPlot(bfood', v, line_colors(2, :), '-', 1, 2);
% [~, lobj(3)] = patchPlot(bloc', v, line_colors(4, :), '-', 1, 2);
% [~, lobj(4)] = patchPlot(bcontrol', v, line_colors(5, :), '-', 1, 2);
% 
% xlabel('voxel');
% ylabel('beta');
% 
% lgd = legend(lobj(1:4), {'PERSON', 'FOOD', 'LOCATION', 'CONTROL'});

%% alternatively -- treat odors from different category as 4 different types of events
%%%% this is the correct method (at least what vivek did). 
SUBJNAMES = {'240711_fMRI_OX_NWU_AS', ...
        '240723_fMRI_OX_NWU_LS', ...
        '240814_fMRI_OX_NWU_JN', ...
        '240816_fMRI_OX_NWU_RR', ...
        '241018_fMRI_OX_NWU_BN', ...
        '250117_fMRI_OX_NWU_VS'}; 

session_count = [4, 13, 10, 15, 16, 10]; % number of sessions for each subj so far. EDIT as needed. 
TR= 0.76; 

subjidx = 5; % enter subjidx here 
subjname = ['subj_', num2str(subjidx)];  
subjname_real = SUBJNAMES{subjidx}; 

% mridir = '/Volumes/ExtremeSSD/OX_DATA/MRI'; 
mridir = '/Users/qhyang/Desktop/OX_DATA/MRI'; 

mridatapath = fullfile(mridir, ['subj_', num2str(subjidx)], 'nifti'); % sn:subject's name

% evdir = '/Volumes/ExtremeSSD/OX_DATA/labchart'; 
evdir = '/Users/qhyang/Desktop/OX_DATA/labchart'; 

nruns = 80; 
stimdur = 2;
% stimdur = 5;

outdir = fullfile(mridatapath, 'single_trial_by_category'); 

%% get design matrix 
%%% design matrix should be stim type * TR. 
%%% we want a cue beta and an odor beta for each voxel to look at patterns 
load(fullfile(evdir, [subjname, '_events.mat'])); % all events 

designmat = zeros(sum(nframes), 4); % TR by stimuli 

% get n of TR from event_onset and cue_onset 
eventTR = ceil(event_onsets_vec/TR); 
% cueTR = ceil(cue_onsets_vec/TR); 


[odor, category] = OX_get_odor(subjname); 
personi = find(strcmp(category, 'PERSON'));
foodi = find(strcmp(category, 'FOOD'));
loci = find(strcmp(category, 'LOCATION')); 
controli = find(strcmp(category, 'CONTROL')); 

designmat(eventTR(personi), 1) = 1;
designmat(eventTR(foodi), 2) = 1;
designmat(eventTR(loci), 3) = 1;
designmat(eventTR(controli), 4) = 1;

% designmat(cueTR(personi), 1) = 1;
% designmat(cueTR(foodi), 2) = 1;
% designmat(cueTR(loci), 3) = 1;
% designmat(cueTR(controli), 4) = 1;

%% get data. use a functional mask/grey matter mask
%%% use a smaller mask for proof of concept 

mask = spm_vol(fullfile(mridatapath, 'coreg', 'gm_mask_thr05_func.nii')); 
maskvol = spm_read_vols(mask); 
maskvol = logical(maskvol);
%%% read func data
funcfiles = func_list(subjname_real, mridatapath, 1, session_count(subjidx));

% Get data from all voxels in the mask -- from vivek's sample script
% ARC_createsingletrials.m 
Res_V = spm_vol(funcfiles);
[~, XYZmm] = spm_read_vols(Res_V{1});
XYZvx = round(Res_V{1}(1).mat\[XYZmm; ones(1,size(XYZmm,2))]);
Res_Vol = cell2mat(Res_V);
mask1D = maskvol(:);
vxl = XYZvx(:,mask1D);
voxel_act = spm_get_data(Res_Vol,vxl)';
voxel_act = single(voxel_act);
[r,c] = find(isnan(voxel_act));
voxel_act(r,:) = []; % this is the input data for 
t_axis = (0:1:sum(nframes))*TR;

%% get nuissance regressors 
load(fullfile(mridatapath, 'motion_param.mat')); 

%% put data into run-wise blocks 
% Change to runwise blocks
% can't seem to run without this. Let vivek know if i figure out why 
design_m2 = cell(1,nruns);
voxel_act2 = cell(1,nruns);
mp2 = cell(1,nruns);
for zz = 1:nruns
    if zz==1
        idx = 1:nframes(1);
    else 
        idx = sum(nframes(1:zz-1))+1:sum(nframes(1:zz)); 
    end 

    design_m2{zz} = designmat(idx,:);
    voxel_act2{zz} = voxel_act(:,idx);
    temp = R(idx,:);
    mp2{zz} = temp;%(:,any(temp)); % Include non-zero regressors
end
opt.extraregressors = mp2;

% %% scale runwise data as percentage change -- this is not necessary.
% design_m2 = cell(1, nruns);
% voxel_act2 = cell(1, nruns);
% mp2 = cell(1, nruns);
% 
% for zz = 1:nruns
%     if zz == 1
%         idx = 1:nframes(1);
%     else 
%         idx = sum(nframes(1:zz-1)) + 1 : sum(nframes(1:zz)); 
%     end 
% 
%     design_m2{zz} = designmat(idx,:);
%     raw_run_data = voxel_act(:, idx);  % size: [nVoxels x nTimepoints]
% 
%     % Convert to percent signal change
%     mean_signal = mean(raw_run_data, 2); % mean across time, per voxel
%     % Avoid divide-by-zero
%     mean_signal(mean_signal == 0) = eps;
%     percent_change = 100 * (raw_run_data - mean_signal) ./ mean_signal;
% 
%     voxel_act2{zz} = percent_change;
% 
%     temp = R(idx,:);
%     mp2{zz} = temp; % extra regressors
% end
% 
% opt.extraregressors = mp2;

%% run glmsingle 

% try
    results = GLMestimatesingletrial(design_m2,voxel_act2,stimdur,TR,outdir,opt);
% catch
%     save(fullfile(outdir,'error_rep.mat'));
% end


%% quick check 
gmmask = spm_vol(fullfile(mridatapath, 'coreg', 'gm_mask_thr05_func.nii'));
gmvol = spm_read_vols(gmmask);  % size: [X Y Z]
gmvol = logical(gmvol);         % binarize
R2_map = zeros(size(gmvol));     % full 3D volume

dataloc = find(gmvol);           % linear indices of voxels in mask
R2_map(dataloc) = mean(squeeze(modelmd), 2);     % assign R² values into volume

%% 
% Reuse the header from your mask
outvol = gmmask;                 
outvol.fname = fullfile(pwd,'quick_check', 'avg_beta_map.nii');  % or a temp path
spm_write_vol(outvol, R2_map);

% Now view it
% spm_image('Display', outvol.fname);



%% 


%% functions 
function filename = func_list(subjname, datapath, nsess_i, nsess_f)

%%%%% read all realigned and smoothed func file 

% rnctr is run counter
rcntr = 0;
filename = [];

for sess = nsess_i:nsess_f

    path_ = fullfile(datapath, 'func/');
    n = dir(fullfile(path_, [subjname, '_', num2str(sess), '_*.nii']));
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
        filename{rcntr, 1} = sprintf('%s', fname);
    end

end


end


