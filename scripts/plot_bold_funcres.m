%% plot functionally restricted BOLD time series 
%% require ROI in func space, cuelist, all event markers from labchart_data_preproc.m 

%% 
SUBJNAMES = {'240711_fMRI_OX_NWU_AS', ...
        '240723_fMRI_OX_NWU_LS', ...
        '240814_fMRI_OX_NWU_JN', ...
        '240816_fMRI_OX_NWU_RR', ...
        '240816_fMRI_OX_NWU_BN'}; 

session_count = [4, 13, 10, 9, 6]; % number of sessions for each subj so far. EDIT as needed. 
TR= 0.76; 

%% enter subjID and session ID here 

% wkdir = '/Volumes/ExtremeSSD/OX_DATA/nifti'; 

subjidx = 3; % enter subjidx here 
subjname = ['subj_', num2str(subjidx)];  
subjname_real = SUBJNAMES{subjidx}; 

% mridir = '/Volumes/ExtremeSSD/OX_DATA/MRI'; 
mridir =  '/Users/qhyang/Desktop/OX_DATA/MRI'; 

mridatapath = fullfile(mridir, ['subj_', num2str(subjidx)], 'nifti'); % sn:subject's name

% evdir = '/Volumes/ExtremeSSD/OX_DATA/labchart'; 
evdir =  '/Users/qhyang/Desktop/OX_DATA/labchart'; 

%% get odor and category variable 
[odor, category] = OX_get_odor(subjname); 

personi = find(strcmp(category, 'PERSON'));
foodi = find(strcmp(category, 'FOOD'));
loci = find(strcmp(category, 'LOCATION')); 
controli = find(strcmp(category, 'CONTROL')); 

categoryidx = zeros(length(category), 1); 
categoryidx(personi) = 1; 
categoryidx(foodi) = 2;
categoryidx(loci) = 3; 
categoryidx(controli) = 4; 


catvec = reshape(categoryidx, [numel(categoryidx), 1]); 
odorvec = reshape(odor, [numel(categoryidx), 1]); 

%% plot
% drawroi = {'mea', 'pir'}; 
% plotidx = 2; 

drawroi = {'AON', 'pirF', 'pirT', 'TU'};
drawroi_amg = {'ACo','CeA', 'MeA', 'PAC', 'LA', 'BMA', 'BLA', 'PCo'};
drawroi_all = {'AON', 'pirF', 'pirT', 'TU', 'ACo','CeA', 'MeA', 'PAC', 'LA', 'BMA', 'BLA', 'PCo'};
%% 
load(fullfile(evdir, [subjname, '_events.mat'])); 
%% ROI loop 

roits = cell(length(drawroi_all), 1); 

for plotidx = 1: length(drawroi_all)
    plotmask = MRIread(fullfile(mridir, subjname, 'nifti', 'coreg', [drawroi_all{plotidx}, '_func_thr02.nii']));
    plotvolloc = find(plotmask.vol >0);
    plotloc = plotvolloc;
    funcfiles = func_list(subjname_real, mridatapath, sessi, sessf);
    
    % read all volumes 
    plotvol_all = []; 
    for runidx = 1:length(funcfiles )
        funcname = funcfiles{runidx};
        mrifile = MRIread(funcname);
        mrvol  = mrifile.vol;
        mrvolvec = reshape(mrvol, [], size(mrvol, 4));
        plotvol = mrvolvec(plotloc, :);
        plotvol_mean = mean(plotvol, 2);
        plotvol = double(plotvol);
        plotvol_per = plotvol./plotvol_mean*100;
        plotvol_all = [plotvol_all, plotvol_per]; 

    end 

    evs = ceil(event_onsets_vec/TR);
    r = [-5: 15];
    event_epochs = zeros(length(r), length(evs));
    event_ts = zeros(size(plotvol, 1), length(r), length(evs)); % voxel by time by trial

    for eventidx = 1: length(evs)
        currentepoch = evs(eventidx)+r;
        event_epochs(:, eventidx) = currentepoch;
        event_ts(:, :, eventidx) = plotvol_all(:, currentepoch);

    end

    roits{plotidx} = event_ts; 
end 


%% save roits 
save(fullfile('MRI', subjname, 'event_ts_percentage.mat'), 'roits', 'drawroi_all', 'odorvec', 'catvec');

%% plot percent change 
line_colors = [[0, 0.4470, 0.7410];...
    [0.8500, 0.3250, 0.0980];...
    [0.9290, 0.6940, 0.1250];...
    [0.4940, 0.1840, 0.5560];...
    [0.4660, 0.6740, 0.1880];...
    [0.3010, 0.7450, 0.9330];...
    [0.6350, 0.0780, 0.1840]];

r = [-5: 15];
figure; 
for plotidx = 1: length(drawroi_all)
    p = roits{plotidx}(:, :, catvec == 1);
    pmean = squeeze(mean(p, 1));
    pmean = pmean';
    f = roits{plotidx}(:, :, catvec == 2);
    fmean = squeeze(mean(f, 1));
    fmean = fmean';
    l = roits{plotidx}(:, :, catvec == 3);
    lmean = squeeze(mean(l, 1));
    lmean = lmean';
    c = roits{plotidx}(:, :, catvec == 4);
    cmean = squeeze(mean(c, 1));
    cmean = cmean';

  
    subplot(6, 2, plotidx);
    hold on

    [~, lobj(1)] = patchPlot(pmean, r, line_colors(1, :), '-', 1, 2);
    [~, lobj(2)] = patchPlot(fmean, r, line_colors(2, :), '-', 1, 2);
    [~, lobj(3)] = patchPlot(lmean, r, line_colors(4, :), '-', 1, 2);
    [~, lobj(4)] = patchPlot(cmean, r, line_colors(5, :), '-', 1, 2);


    xlabel('Time(TR)');
    ylabel('BOLD');
%     ylim([99, 101]);
    lgd = legend(lobj(1:4), {'PERSON', 'FOOD', 'LOCATION', 'CONTROL'});
    lgd.Location = 'best';

    title(['subj ', num2str(subjidx),', ', drawroi_all{plotidx}, ', average BOLD percentage change per category']);
end

%% restriction 
roivpick = []; 
npick = []; 
for plotidx = 1: length(drawroi_all)
    event_ts = roits{plotidx}; 

    %%% 
    event_pk = max(event_ts, [], 2); 
    event_sd = std(event_ts, [], 2); 
    
    pkgmm_model = fitgmdist(event_pk(:), 2); 
    %%% 
    peakthr = min(pkgmm_model.mu); 
    vpick = squeeze(mean(event_pk, 3))>peakthr; 
    
    roivpick{plotidx} = vpick; 
    npick(plotidx) = sum(vpick); 


end 


%%
line_colors = [[0, 0.4470, 0.7410];...
    [0.8500, 0.3250, 0.0980];...
    [0.9290, 0.6940, 0.1250];...
    [0.4940, 0.1840, 0.5560];...
    [0.4660, 0.6740, 0.1880];...
    [0.3010, 0.7450, 0.9330];...
    [0.6350, 0.0780, 0.1840]];

r = [-5: 15];
figure; 
for plotidx = 1: length(drawroi_all)
    p = roits{plotidx}(roivpick{plotidx}, :, catvec == 1);
    pmean = squeeze(mean(p, 1));
    pmean = pmean';
    f = roits{plotidx}(roivpick{plotidx}, :, catvec == 2);
    fmean = squeeze(mean(f, 1));
    fmean = fmean';
    l = roits{plotidx}(roivpick{plotidx}, :, catvec == 3);
    lmean = squeeze(mean(l, 1));
    lmean = lmean';
    c = roits{plotidx}(roivpick{plotidx}, :, catvec == 4);
    cmean = squeeze(mean(c, 1));
    cmean = cmean';

  
    subplot(6, 2, plotidx);
    hold on

    [~, lobj(1)] = patchPlot(pmean, r, line_colors(1, :), '-', 1, 2);
    [~, lobj(2)] = patchPlot(fmean, r, line_colors(2, :), '-', 1, 2);
    [~, lobj(3)] = patchPlot(lmean, r, line_colors(4, :), '-', 1, 2);
    [~, lobj(4)] = patchPlot(cmean, r, line_colors(5, :), '-', 1, 2);


    xlabel('Time(TR)');
    ylabel('BOLD');
    ylim([99.5, 101]);
    lgd = legend(lobj(1:4), {'PERSON', 'FOOD', 'LOCATION', 'CONTROL'});
    lgd.Location = 'best';

    title(['subj ', num2str(subjidx),', ', drawroi_all{plotidx}, ', average BOLD percentage change per category']);
end

%%

SetPrintProp(gcf, 0.4, 0.8); 

%%%% minimal overlapping between cue contrast mask and olfactory regions 


%%% 
% print(fullfile(mridir, subjname, 'BOLD_all_preCueBsl'), '-dpdf', '-fillpage'); 
%%     
 

%         evs_seconds = event_onsets(:, runidx);
% 
%         evs = ceil(evs_seconds/TR);  % choose the next TR after event onset.
%         r = [-5: 15];
%         event_epochs = zeros(length(r), length(evs));
%         event_ts = zeros(size(plotvol, 1), length(r), length(evs)); % voxel by time by trial
% 
%         for eventidx = 1: length(evs)
%             currentepoch = evs(eventidx)+r;
%             event_epochs(:, eventidx) = currentepoch;
%             event_ts(:, :, eventidx) = plotvol(:, currentepoch);
% 
%         end
% 
%         bsl = [-3, -1];  % take -3~-1 TR as baseline
%         bslloc = G_RangeLoc(r, bsl);
%         bslval = mean(event_ts(:, bslloc, :), 2);
%         event_ts_bslcorr = event_ts - bslval;
%         %     event_ts_bslcorr = event_ts./bslval;
% 
%         if runidx ==1
%             event_ts_all = event_ts;
%             event_ts_bsl_all = event_ts_bslcorr;
%         else
%             event_ts_all = cat(3, event_ts_all, event_ts);
%             event_ts_bsl_all = cat(3, event_ts_bsl_all, event_ts_bslcorr);
%         end


        %    plotts = mean(squeeze(mean(event_ts_bslcorr, 1)), 2);
        %    plotts = plotts - mean(plotts);
        %    plot(r, plotts);

%     end
% 
% end 



%%
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
