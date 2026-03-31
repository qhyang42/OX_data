%% single trial analysis. perform single trial estimate before this step. 
%% load data 

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

outdir = fullfile(mridatapath, 'single_trial_by_category');

load(fullfile(outdir, 'TYPED_FITHRF_GLMDENOISE_RR.mat')); 
[odor, category] = OX_get_odor(subjname); 

personi = find(strcmp(category, 'PERSON'));
foodi = find(strcmp(category, 'FOOD'));
loci = find(strcmp(category, 'LOCATION')); 
controli = find(strcmp(category, 'CONTROL')); 

modelmd = squeeze(modelmd);

%% set nan value to 0
% modelmd(isnan(modelmd)) = 0; 

%%
bperson = modelmd(:, personi);
bfood = modelmd(:, foodi);
bloc = modelmd(:, loci); 
bcontrol = modelmd(:, controli);

%% R2 mask 
rmask = R2>0.5; 
bperson = bperson(rmask); 
bfood = bfood(rmask); 
bloc = bloc(rmask); 
bcontrol = bcontrol(rmask); 

%% 













%% map trial averaged beta back to mask 

% gmmask = spm_vol(fullfile(mridatapath, 'coreg', 'gm_mask_thr05_func.nii')); 
% gmvol = spm_read_vols(gmmask); 
% gmvol = logical(gmvol);
% 
% dataloc = find(gmvol); 
% 
% controlmap = zeros(size(gmvol)); 
% foodmap = zeros(size(gmvol)); 
% locmap = zeros(size(gmvol)); 
% personmap = zeros(size(gmvol)); 
% 
% for i = 1: length(dataloc)
%     currentloc = dataloc(i);
%     cc = mean(bcontrol(i, :));
%     cp = mean(bperson(i, :));
%     cf = mean(bfood(i, :));
%     cl = mean(bloc(i, :));
%     
%     controlmap(dataloc(i)) = cc; 
%     foodmap(dataloc(i)) = cf; 
%     locmap(dataloc(i)) = cl; 
%     personmap(dataloc(i)) = cp; 
% 
% end 

%% save map
% personvol = spm_create_vol(gmmask); 
% personvol.fname = fullfile(outdir, 'glmsingle_person_beta.nii');
% spm_write_vol(personvol,personmap); 
% 
% foodvol = spm_create_vol(gmmask); 
% foodvol.fname = fullfile(outdir, 'glmsingle_food_beta.nii');
% spm_write_vol(foodvol,foodmap); 
% 
% locvol = spm_create_vol(gmmask); 
% locvol.fname = fullfile(outdir, 'glmsingle_loc_beta.nii');
% spm_write_vol(locvol,locmap); 
% 
% controlvol = spm_create_vol(gmmask); 
% controlvol.fname = fullfile(outdir, 'glmsingle_control_beta.nii');
% spm_write_vol(controlvol,controlmap); 

%%% beta map mostly overlapping like single HRF models 

%% nan values in output? 
%%% checked map. a thin strip on occipital lobe. not our primary area of
%%% interest 
%%% the nan values are set to 0 for the map

%% 
betas = [bperson, bfood, bloc, bcontrol]; 

%%
for vidx = 1: size(betas, 1)
    

end 






%% plot cov matrix -- find voxels that have higher within context correlation than between context correlation 


nanmask = ~isnan(modelmd(:, 1)); % voxels with value 
nanidx = find(nanmask);  % voxel index with value. use this to map result back to grey matter mask 


% data = modelmd(nanmask, :);
data = modelmd(nanmask & rmask, :);

nrep = 100; 
ntrials = 100; % select trials

pcorr = zeros(size(data, 1), nrep); 
fcorr = zeros(size(data, 1), nrep); 
lcorr = zeros(size(data, 1), nrep); 
ccorr = zeros(size(data, 1), nrep); 


for vidx = 1: size(data, 1)
    if mod(vidx, 100) == 0
    fprintf('voxel %d   \n', vidx);
    end 
    for repidx = 1:nrep
        tselect1 = randperm(200, ntrials); 
        tselect2 = randperm(200, ntrials);
        pcorr(vidx, repidx )  = corr(data(vidx, personi(tselect1))', data(vidx, personi(tselect2))');
        fcorr(vidx, repidx )  = corr(data(vidx, foodi(tselect1))', data(vidx, foodi(tselect2))');
        lcorr(vidx, repidx )  = corr(data(vidx, loci(tselect1))', data(vidx, loci(tselect2))');
        ccorr(vidx, repidx )  = corr(data(vidx, controli(tselect1))', data(vidx, controli(tselect2))');

%         fprintf('rep %d    \n', repidx);
    end
end


%% 
pcorrn = atanh(pcorr); 
fcorrn = atanh(fcorr); 
lcorrn = atanh(lcorr); 
ccorrn = atanh(ccorr); 

plotdata = [mean(pcorrn, 2), mean(fcorrn, 2), mean(lcorrn, 2), mean(ccorrn, 2)]; 


%% baseline 

bslcorr = zeros(size(data, 1), nrep);

for vidx = 1: size(data, 1)
    if mod(vidx, 500) == 0
    fprintf('voxel %d   \n', vidx);
    end 
for repidx = 1:nrep
    tselect1 = randperm(800, ntrials); 
    tselect2 = randperm(800, ntrials); 
    bslcorr(vidx, repidx )  = corr(data(vidx, tselect1)', data(vidx, tselect2)');

    %         fprintf('rep %d    \n', repidx);
end
end 

bslcorrn = atanh(bslcorr); 
bslcorrn = mean(bslcorrn, 2); 
thrr = sort(bslcorrn(:));
thr = thrr(floor(length(thrr)*0.95));

%% voxel by voxel comparison 
pcorrn_bsl = pcorrn - atanh(bslcorr); 
fcorrn_bsl = fcorrn - atanh(bslcorr); 
lcorrn_bsl = lcorrn - atanh(bslcorr); 
ccorrn_bsl = ccorrn - atanh(bslcorr); 

%%% todo t test against 0? 
%%% this map is also sparse 
%% try saving map 
gmmask = spm_vol(fullfile(mridatapath, 'coreg', 'gm_mask_thr05_func.nii')); 
gmvol = spm_read_vols(gmmask); 
gmvol = logical(gmvol);

dataloc = find(gmvol); 

controlmap = zeros(size(gmvol)); 
foodmap = zeros(size(gmvol)); 
locmap = zeros(size(gmvol)); 
personmap = zeros(size(gmvol)); 

for i = 1: length(nanidx)
    currentloc = dataloc(nanidx(i)); % todo. any way to check if this is absolutely correct? 
    cc = plotdata(i, 4); 
    cp = plotdata(i, 1); 
    cf = plotdata(i, 2); 
    cl = plotdata(i, 3); 
    
    controlmap(currentloc) = cc; 
    foodmap(currentloc) = cf; 
    locmap(currentloc) = cl; 
    personmap(currentloc) = cp; 

end 


%% save map -- maps should be thresholded at 0.04; 
personvol = spm_create_vol(gmmask); 
personvol.fname = fullfile(outdir, 'rsa_person_diff.nii');
spm_write_vol(personvol,personmap); 

foodvol = spm_create_vol(gmmask); 
foodvol.fname = fullfile(outdir, 'rsa_food.nii');
spm_write_vol(foodvol,foodmap); 

locvol = spm_create_vol(gmmask); 
locvol.fname = fullfile(outdir, 'rsa_loc.nii');
spm_write_vol(locvol,locmap); 

controlvol = spm_create_vol(gmmask); 
controlvol.fname = fullfile(outdir, 'rsa_control.nii');
spm_write_vol(controlvol,controlmap); 

%%% looks very veeeerrryyyy sparse and super noisy 


%% get ROI locations  
olfroi = {'AON', 'pirF', 'pirT', 'TU'};
amgroi = {'ACo','CeA', 'MeA', 'PAC', 'LA', 'BMA', 'BLA', 'PCo'};
allroi = {'AON', 'pirF', 'pirT', 'TU','ACo','CeA', 'MeA', 'PAC', 'LA', 'BMA', 'BLA', 'PCo'};

gmmask = spm_vol(fullfile(mridatapath, 'coreg', 'gm_mask_thr05_func.nii')); 
gmvol = spm_read_vols(gmmask); 
gmloc = find(gmvol); 

roiloc_all = zeros(length(allroi), length(gmloc)); 
%%% get overlap between ROI and GM mask
for roiidx = 1:length(allroi)
    roimask = spm_vol(fullfile(mridatapath, 'coreg', [allroi{roiidx}, '_func_thr02.nii'])); 
    roivol = spm_read_vols(roimask); 
    roiloc = find(roivol); 

    gmroiloc = intersect(gmloc, roiloc); 
    roi1d = ismember(gmloc, gmroiloc); 
    
    roiloc_all(roiidx, :) = roi1d; 
end 

%%


roimask = spm_vol(fullfile(mridatapath, 'coreg', [allroi{roiidx}, '_func_thr02.nii']));
roivol = spm_read_vols(roimask);
roiloc = find(roivol);

gmroiloc = intersect(gmloc, roiloc);
roi1d = ismember(gmloc, gmroiloc);

roiloc_all(roiidx, :) = roi1d;

%% get beta per ROI 
plotdata = zeros(length(allroi), 800); 
for roiidx = 1:length(allroi)
    b = modelmd(find(roiloc_all(roiidx, :)), :); 
    b = mean(b, 1); 
    plotdata(roiidx, :) = b; 
end 
%% 

%% try some other sanity check measure 
% compare odor vs no odor beta 

