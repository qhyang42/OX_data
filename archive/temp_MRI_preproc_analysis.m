%% 
SUBJNAMES = {'240711_fMRI_OX_NWU_AS', ...
        '240723_fMRI_OX_NWU_LS', ...
        '240814_fMRI_OX_NWU_JN', ...
        '240816_fMRI_OX_NWU_RR'}; 

session_count = [4, 3, 4, 4]; % number of sessions for each subj so far. EDIT as needed. 
subjidx = 1; % enter subj name here 


%% 
subjname = SUBJNAMES{subjidx}; 
wkdir =     '/Volumes/ExtremeSSD/OX_DATA/MRI'; 
datapath = fullfile(wkdir, ['subj_', num2str(subjidx)], 'nifti'); % sn:subject's name

%% 
behaviordata = load(fullfile(wkdir, 'pilotsubj2_ses1.mat')); 
datastart = behaviordata.datastart; 
dataend = behaviordata.dataend; 
data = behaviordata.data; 

mripulse = cell(6,1);
events = cell(6,1);
resp = cell(6,1);

for i = 1:6
    mripulse{i} = data(datastart(2, i+2):dataend(2, i+2)); 
    resp{i} = data(datastart(1, i+2):dataend(1, i+2)); 
    events{i} = data(datastart(3, i+2):dataend(3, i+2)); 

end 


%% 
%%%%% MB2: run1 443/1.4, run4 440/1.4
%%%%% MB3: run3 993/0.62, run6 998/0.62
%%%%% MB4: run2 951/0.645, run5 954/0.645

nframes = [443; 951; 993; 440; 954; 998]; 
runtime = [443*1.4; 951*0.645; 993*0.62; 440*1.4; 954*0.645; 998*0.62]; 

%% 
event_onsets = zeros(20,6); 
cue_onsets = zeros(20, 6); 

for evidx = 1:6
    currentevents = events{evidx};
    eventbi = currentevents>0.1;
    event_real = eventbi;
    for i = 1:length(eventbi)-1
        if eventbi(i+1) == eventbi(i)
            event_real(i+1) = 0;
        end
    end
    evonsets = find(event_real>0);
    event_o = evonsets(5:2:end);
    event_q = evonsets(4:2:end-1);

    mrionset = find(mripulse{evidx} > 1, 1); 
    event_o = event_o - mrionset; % set MRI onset to time 0 
    event_o = event_o/1000; % divide by fs to get time.  
    event_onsets(:, evidx) = event_o; 
    
    event_q = event_q - mrionset; % set MRI onset to time 0 
    event_q = event_q/1000;
    cue_onsets(:, evidx) = event_q; 

end


%% calculate offset for each acqiusition 


mb2offset = runtime(1); 
mb4offset = runtime(2); 
mb3offset = runtime(3);

mb2event_onsets = [event_onsets(:, 1), event_onsets(:, 4)]; 
mb4event_onsets = [event_onsets(:, 2), event_onsets(:, 5)]; 
mb3event_onsets = [event_onsets(:, 3), event_onsets(:, 6)]; 

mb2event_onsets_all = [event_onsets(:, 1), event_onsets(:, 4) + mb2offset]; 
mb4event_onsets_all = [event_onsets(:, 2), event_onsets(:, 5) + mb4offset]; 
mb3event_onsets_all = [event_onsets(:, 3), event_onsets(:, 6) + mb3offset]; 

mb2event_onsets_vec = reshape(mb2event_onsets_all, [numel(mb2event_onsets_all), 1]); 
mb3event_onsets_vec = reshape(mb3event_onsets_all, [numel(mb3event_onsets_all), 1]); 
mb4event_onsets_vec = reshape(mb4event_onsets_all, [numel(mb4event_onsets_all), 1]); 

%% make cue onsets. 
mb2cue_onsets = [cue_onsets(:, 1), cue_onsets(:, 4)]; 
mb4cue_onsets = [cue_onsets(:, 2), cue_onsets(:, 5)]; 
mb3cue_onsets = [cue_onsets(:, 3), cue_onsets(:, 6)]; 

mb2cue_onsets_all = [cue_onsets(:, 1), cue_onsets(:, 4) + mb2offset]; 
mb4cue_onsets_all = [cue_onsets(:, 2), cue_onsets(:, 5) + mb4offset]; 
mb3cue_onsets_all = [cue_onsets(:, 3), cue_onsets(:, 6) + mb3offset]; 

mb2cue_onsets_vec = reshape(mb2cue_onsets_all, [numel(mb2cue_onsets_all), 1]); 
mb3cue_onsets_vec = reshape(mb3cue_onsets_all, [numel(mb3cue_onsets_all), 1]); 
mb4cue_onsets_vec = reshape(mb4cue_onsets_all, [numel(mb4cue_onsets_all), 1]); 



% event_onsets_vec = reshape(event_onsets_all, [numel(event_onsets_all), 1]); 
 
% cue_onsets_all = offsets+cue_onsets;
% cue_onsets_vec = reshape(cue_onsets_all, [numel(cue_onsets_all), 1]); 
%% 

cd '/Users/qhyang/Desktop/OX_pilot02/nifti'
datadir = 'MB3/'; 
load(fullfile(datadir, 'odor_onsets.mat')); 
% pirmask = MRIread('T1/r_MB3_T1_pir_mask_bin.nii.gz'); 
% pirmask = MRIread('T1/MB3_temp_cue_mask.nii.gz'); % mask covering significant cue response 
pirmask = MRIread('T1/MB3_mea_mask_bin.nii.gz');

pirvolloc = find(pirmask.vol >0); 
datafiles = dir(fullfile(datadir, 'r*run*.nii')); 

figure; 
hold on

plotts = []; 
for runidx = 1:2
   dataname = datafiles(runidx).name; 
   mrifile = MRIread(fullfile(datadir, dataname)); 
   mrvol  = mrifile.vol; 
   mrvolvec = reshape(mrvol, [], size(mrvol, 4)); 
   pirvol = mrvolvec(pirvolloc, :); 
   
   evs_seconds = mb3event_onsets(:, runidx); 


   evs = ceil(evs_seconds/tr);  % choose the next TR after event onset. 
   r = [-5: 15];
   r = ceil(r/tr); 
   event_epochs = zeros(length(r), length(evs)); 
   event_ts = zeros(size(pirvol, 1), length(r), length(evs)); % voxel by time by trial 
   
   for eventidx = 1: length(evs)
        currentepoch = evs(eventidx)+r; 
        event_epochs(:, eventidx) = currentepoch; 
        event_ts(:, :, eventidx) = pirvol(:, currentepoch);

   end 

   bsl = [-5, -1];  % take baseline
   bsl = ceil(bsl/tr); 
   bslloc = G_RangeLoc(r, bsl); 
   bslval = mean(event_ts(:, bslloc, :), 2);
   event_ts_bslcorr = event_ts - bslval;
   %     event_ts_bslcorr = event_ts./bslval;

   if runidx ==1
       event_ts_all = event_ts;
       event_ts_bsl_all = event_ts_bslcorr; 
   else
       event_ts_all = cat(3, event_ts_all, event_ts);
       event_ts_bsl_all = cat(3, event_ts_bsl_all, event_ts_bslcorr); 
   end


   plotts(:, runidx) = mean(squeeze(mean(event_ts_bslcorr, 1)), 2); 
%    plotts = plotts - mean(plotts); 
   plot(r*tr, plotts(:, runidx));
   
end 

% figure; 
% hold on; 
% 
% plot(r*tr, mean(plotts, 1)); 
% patchPlot(plotts, r*tr); 
xlabel('time(s)'); 
ylabel('BOLD')
% xlim([-5, 30]);


%% load cuelist. The cuelist is the same as pilotsubj1 
category = cell(20, 6); 
cueidx = zeros(20, 6); 
odor = zeros(20, 6); 

cuedir = '/Users/qhyang/Desktop/OX_PILOT01_VS_20240203_OX_PILOT01_VS/pilot_subj_1_cuelist'; 

for runidx  = 1: 6
    cue = load(fullfile(cuedir, ['cuelist_sess1_run', num2str(runidx), '.mat'])); 
    category(:, runidx) = cue.cuelist.category';
    cueidx(:, runidx) = cue.cuelist.cueidx;
    odor(:, runidx) = cue.cuelist.odor; 
end 

personi = find(strcmp(category, 'PERSON'));
foodi = find(strcmp(category, 'FOOD'));
loci = find(strcmp(category, 'LOCATION')); 
controli = find(strcmp(category, 'CONTROL')); 

categoryidx = zeros(20, 6); 
categoryidx(personi) = 1; 
categoryidx(foodi) = 2;
categoryidx(loci) = 3; 
categoryidx(controli) = 4; 

mb2catvec = categoryidx(:, [1, 4]); 
mb2catvec = mb2catvec(:); 

mb4catvec = categoryidx(:, [2, 5]); 
mb4catvec = mb4catvec(:); 

mb3catvec = categoryidx(:, [3, 6]); 
mb3catvec = mb3catvec(:); 

%% 
line_colors = [[0, 0.4470, 0.7410];...	          
          	[0.8500, 0.3250, 0.0980];...	       
          	[0.9290, 0.6940, 0.1250];...	          
          	[0.4940, 0.1840, 0.5560];...	          
          	[0.4660, 0.6740, 0.1880];...	          
          	[0.3010, 0.7450, 0.9330];...	          
          	[0.6350, 0.0780, 0.1840]]; 	


% line_colors10 = [0, 18, 25; ... 
%     0, 95, 115; ...
%     10, 147, 150; ...
%     148, 210, 189; ... 
%     233, 216, 166; ...
%     238, 155, 0; ...
%     202, 103, 2; ...
%     187, 62, 3; ...
%     174, 32, 18; ...
%     155, 34, 38;]; 
% 
% line_colors10 = line_colors10/255; 


% line_colors10 = ["#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2"]; 


line_colors10 = colormap(parula(10));



 r = [-5: 15];
r = ceil(r/tr); 

cond_type = 'category'; %| 'odor'; 

catvec = mb3catvec; 
switch cond_type
    case 'category' 
        p = event_ts_bsl_all(:, :, catvec == 1); 
        pmean = squeeze(mean(p, 1)); 
        pmean = pmean'; 
        f = event_ts_bsl_all(:, :, catvec == 2);
        fmean = squeeze(mean(f, 1)); 
        fmean = fmean'; 
        l = event_ts_bsl_all(:, :, catvec == 3);
        lmean = squeeze(mean(l, 1)); 
        lmean = lmean'; 
        c = event_ts_bsl_all(:, :, catvec == 4);
        cmean = squeeze(mean(c, 1)); 
        cmean = cmean'; 

%         figure; 
%         hold on 
% 
%         plot(r, mean(squeeze(mean(p, 1)), 2)); 
%         plot(r, mean(squeeze(mean(f, 1)), 2)); 
%         plot(r, mean(squeeze(mean(l, 1)), 2)); 
%         plot(r, mean(squeeze(mean(c, 1)), 2)); 
% 
%         xlabel('Time(s)'); 
%         ylabel('BOLD'); 
%         legend({'PERSON', 'FOOD', 'LOCATION', 'CONTROL'}); 

        figure; 
        hold on 

        [~, lobj(1)] = patchPlot(pmean, r, line_colors(1, :), '-', 1, 2); 
        [~, lobj(2)] = patchPlot(fmean, r, line_colors(2, :), '-', 1, 2); 
        [~, lobj(3)] = patchPlot(lmean, r, line_colors(4, :), '-', 1, 2);
        [~, lobj(4)] = patchPlot(cmean, r, line_colors(5, :), '-', 1, 2); 


        xlabel('Time(s)'); 
        ylabel('BOLD'); 
        legend(lobj(1:4), {'PERSON', 'FOOD', 'LOCATION', 'CONTROL'}); 





    case 'odor' 
    
        figure; 
        hold on 
        
        label = cell(10, 1); 
        for odoridx = 1:10 
            currento = event_ts_bsl_all(:, :, odorvec == odoridx);    
            omean = squeeze(mean(currento, 1)); 
            omean = omean'; 
            
            [~, lobj(odoridx)] = patchPlot(omean, r, line_colors10(odoridx, :), '-', 1, 2); 
            label{odoridx} = ['ODOR ', num2str(odoridx)]; 
            
        end 
        xlabel('Time(s)'); 
        ylabel('BOLD'); 
        legend(lobj(1:10), label); 


end 


%% try plot along first PC 

pcadata = squeeze(mean(event_ts_bsl_all, 3));
pcadata = pcadata';
[coef, score, ~, ~, explained] = pca(pcadata);


p = event_ts_bsl_all(:, :, catvec == 1);

f = event_ts_bsl_all(:, :, catvec == 2);

l = event_ts_bsl_all(:, :, catvec == 3);

c = event_ts_bsl_all(:, :, catvec == 4);

pmean = squeeze(mean(p, 3));
lmean = squeeze(mean(l, 3));
fmean = squeeze(mean(f, 3));
cmean = squeeze(mean(c, 3));
lscore = lmean'*coef;
pscore = pmean'*coef;
fscore = fmean'*coef;
cscore = cmean'*coef;

figure;
hold on
plot(r, pscore(:, 1)); 
plot(r, fscore(:, 1)); 
plot(r, lscore(:, 1)); 
plot(r, cscore(:, 1)); 
xlabel('Time(s)');
ylabel('BOLD');
legend({'PERSON', 'FOOD', 'LOCATION', 'CONTROL'});



%% make motion parameter regressors 
% 
% mp = load('rp_run1_optcom_bold.txt');
% motion_regressor= [mp, [zeros(1,6); diff(mp)], mp.^2, [zeros(1,6); diff(mp).^2]];


motion_regressor = []; 
for runidx = 1:8
    mpfile = ['rp_run', num2str(runidx), '_optcom_bold.txt']; 
    mp = load(mpfile); 
    current_reg = [mp, [zeros(1,6); diff(mp)], mp.^2, [zeros(1,6); diff(mp).^2]];
    motion_regressor = [motion_regressor; current_reg]; 
end 

%%%%% adjust spm.mat for sessions 
% spm_fmri_concatenate('SPM.mat', nframes');

% nframes = [443; 951; 993; 440; 954; 998]; 
% spm_fmri_concatenate('SPM.mat', [951; 954]');



%% behavior results 

rsp_dir =     '/Users/qhyang/Desktop/OX_pilot02/pilot_subj_2'; 

rsp = zeros(20, 6); 
int_rsp = zeros(20, 6); 
for runidx = 1:6
    rspdata = load(fullfile(rsp_dir, ['pilot_subj_2_session1_run', num2str(runidx), '_results.mat'])); 
    for t = 1:20 
        if isempty(rspdata.outMat{t}{4})
            rsp(t, runidx) = 700;
        else 
        rsp(t, runidx) = rspdata.outMat{t}{4}; 
        end 

        if isempty(rspdata.outMat{t}{5})
            int_rsp(t, runidx) = 700; 
        else 
        int_rsp(t, runidx) = rspdata.outMat{t}{5}; 
        end 


    end 
end 


%% 

run_rsp = [int_rsp(:, 1); int_rsp(:, 4)]; 
posloc = run_rsp > 800; % right side 

postrial = event_ts_bsl_all(:, :, posloc);
negtrial = event_ts_bsl_all(:, :, ~posloc); 

posmean = squeeze(mean(postrial, 1)); 
negmean = squeeze(mean(negtrial, 1)); 


r = [-5: 15];
r = ceil(r/tr); 

figure;
hold on

[~, lobj(1)] = patchPlot(posmean', r, line_colors(1, :), '-', 1, 2);
[~, lobj(2)] = patchPlot(negmean', r, line_colors(2, :), '-', 1, 2);


xlabel('Time(s)');
ylabel('BOLD');
legend(lobj(1:2), {'POS', 'NEG'});

%% 
