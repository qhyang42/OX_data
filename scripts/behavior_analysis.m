%% Descriptive behavior correlation boxplots from behavior.mat
% Retain the original 100 sampled odor profiles and mixed-context baseline.
% Correlations share trials/profiles; no inferential p-values are computed.

%% Paths and reproducible sampling (both measures, all subjects 2-6)
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root, 'utils', 'OX_utilities'));
wkdir = fullfile(project_root, 'behavior');
output_dir = fullfile(project_root, 'results', 'behavior');
if isfolder(output_dir)
    output_dir = fullfile(output_dir, char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS')));
end
mkdir(output_dir);
seed = 20260909;
rng(seed, 'twister');
summary = table();
trial_table = table();
condition_summary = table();
measurelabel = {'pleasantness', 'intensity'};

for subject_id = 2:6
    subjname = sprintf('subj_%d', subject_id);
    input_file = fullfile(wkdir, subjname, 'behavior.mat');
    behavior = load(input_file);
    metadata = OX_load_trial_metadata(subject_id);
    assert(isequal(double(behavior.odor(:)), metadata.odor));
    assert(isequal(upper(strtrim(string(behavior.category(:)))), metadata.context));
    odor = metadata.odor;
    category = cellstr(metadata.context);
    valence_all = double(behavior.valence_all(:));
    intensity_all = double(behavior.intensity_all(:));
    assert(numel(valence_all) == 800 && numel(intensity_all) == 800);
    assert(all(accumarray(metadata.run_id, 1) == 10));
    assert(all(arrayfun(@(c) sum(metadata.context == c), ...
        ["CONTROL","FOOD","PERSON","LOCATION"]) == 200));
    assert(~any(isinf([valence_all; intensity_all])), 'Infinite rating in %s.', subjname);
    subject_trials = addvars(metadata, repmat(subject_id,800,1), ...
        valence_all, intensity_all, 'NewVariableNames', ...
        {'subject_id','pleasantness','intensity'});
    trial_table = [trial_table; subject_trials];
for measureidx = 1:2
    rng_state = rng;
%% 
%%%% index context labels 
catidx = zeros(size(odor)); 
for i = 1: length(odor) 
    this_cat = category{i};
    switch this_cat
        case 'CONTROL'
            x = 1;
        case 'FOOD'
            x = 2;
        case 'PERSON'
            x = 3;
        case 'LOCATION'
            x = 4;
    end 
    
    catidx(i) = x;

end 
catlabels = {'CONTROL', 'FOOD', 'PERSON', 'LOCATION'};

%%%% group trials according to odor 
data = cell(length(unique(odor)), max(catidx));

for i = 1:length(odor)
    
    thisodor = odor(i); 
    thiscat = catidx(i);
    thisv = valence_all(i); 
    thisi = intensity_all(i); 
    appi = [thisv, thisi]; 
    data{thisodor, thiscat} = [data{thisodor, thiscat}; appi ]; 
end 

%% Validate all conditions; report missing ratings without imputation.
for n = 1:20
    for c = 1:4
        ratings = data{n,c}(:,measureidx);
        assert(~isempty(ratings), 'Missing odor/context condition.');
        valid = ratings(~isnan(ratings));
        assert(~isempty(valid), 'No ratings for %s odor %d context %d.', subjname,n,c);
        condition_summary = [condition_summary; table(subject_id, ...
            string(measurelabel{measureidx}), n, string(catlabels{c}), ...
            numel(ratings), numel(valid), sum(isnan(ratings)), mean(valid), ...
            std(valid), median(valid), 'VariableNames', ...
            {'subject_id','measure','odor','context','n_trials','n_valid', ...
             'n_missing','mean','sd','median'})];
    end
end

%% calculate inter-odor trajectroy similarity for each category  
odoridx = 1:20; 
nrep = 100; 
v = zeros(nrep, length(odoridx), max(catidx)); 
for repidx = 1:nrep
    for c = 1: max(catidx)

        for n = 1: length(odoridx)
            ratings = data{n, c};
            ratings = ratings(:, measureidx);
            ratings = ratings(~isnan(ratings));
            v(repidx, n, c) = ratings(randperm(length(ratings), 1));
        end
    end

end

%%% 
corr_c = zeros(nrep*(nrep-1), max(catidx)); % all the off diagonal correlation coefficient  
for c = 1:max(catidx)
    v_c = squeeze(v(:, :, c)); 
    r = corrcoef(v_c'); 
    r = r(~eye(nrep));
%     r = tanh(r); 
    corr_c(:, c) = r; 
end 

%% cross condition correlations 
rcross = zeros(size(corr_c, 1), 1); 
for repidx = 1:size(corr_c, 1) 
    c = randperm(max(catidx), 2); % pick 2 random context  
    i = randperm(size(v, 1), 2); % pick 1 random dot set 
    thisr = corrcoef(squeeze(v(i(1), :, c(1))), squeeze(v(i(2),:, c(2)))); 
    thisr = thisr(1,2); % off diagonal 
    rcross(repidx, 1) = thisr; 
end 
%% calculate null distribution for the correlation coefficient 
nperm = 100; 
vperm = zeros(nperm, length(odoridx)); 
for repidx = 1: nperm 
    for n = 1: length(odoridx)
        c = randperm(max(catidx), 1); 
        ratings = data{n, c}; 
        ratings = ratings(:, measureidx);
            ratings = ratings(~isnan(ratings));
        vperm(repidx, n) = ratings(randperm(length(ratings), 1)); 
    end 
end 

rperm = corrcoef(vperm'); 
rperm = rperm(~eye(nperm));
rperm = sort(rperm); 
%% add mean rperm line to the box plot 
bsl = median(atanh(rperm)); 


fig = figure('Visible', 'off', 'Position', [100 100 1100 650]);
boxh = boxplot([atanh(corr_c), atanh(rcross)]); 
hold on 
ph = plot([0:6], bsl*ones(1,7), 'k', 'LineStyle', '-.'); 
ylabel('Intertrial correlation (Fisher z)');
xlabel('context'); 
xticklabels([catlabels, {'cross condition'}]); 
title([subjname, ', ', measurelabel{measureidx}], 'Interpreter', 'none');
set(boxh, 'LineWidth' , 1.5); 
set(ph, 'LineWidth', 1.5); 
% ylim([-0.5, 1]); 

saveas(fig, fullfile(output_dir, [subjname, '_', measurelabel{measureidx}, '_baseline.png']));
savefig(fig, fullfile(output_dir, [subjname, '_', measurelabel{measureidx}, '_baseline.fig']));
close(fig);

%% add rperm as box plot 
bsl = median(atanh(rperm)); 
rperm_box = atanh(rperm); 


fig = figure('Visible', 'off', 'Position', [100 100 1100 650]);
boxh = boxplot([atanh(corr_c), atanh(rcross), rperm_box]); 
hold on 
ph = plot([0:6], bsl*ones(1,7), 'k', 'LineStyle', '-.'); 
ylabel('Intertrial correlation (Fisher z)');
xlabel('context'); 
xticklabels([catlabels, {'cross condition', 'permuted'}]); 
title([subjname, ', ', measurelabel{measureidx}], 'Interpreter', 'none');
set(boxh, 'LineWidth' , 1.5); 
set(ph, 'LineWidth', 1.5); 

saveas(fig, fullfile(output_dir, [subjname, '_', measurelabel{measureidx}, '_resampled.png']));
savefig(fig, fullfile(output_dir, [subjname, '_', measurelabel{measureidx}, '_resampled.fig']));
close(fig);

% These are descriptive resampling distributions, not independent samples.
z_values = [atanh(corr_c), atanh(rcross), rperm_box];
assert(all(isfinite(z_values(:))), 'Nonfinite sampled correlations in %s.', subjname);
box_labels = [string(catlabels), "cross condition", "permuted"];
for c = 1:6
    q = prctile(z_values(:,c), [25 50 75]);
    summary = [summary; table(subject_id, string(measurelabel{measureidx}), ...
        box_labels(c), size(z_values,1), mean(z_values(:,c)), std(z_values(:,c)), ...
        q(1), q(2), q(3), bsl, 'VariableNames', ...
        {'subject_id','measure','box','n_sampled_correlations','mean_z','sd_z', ...
         'q25_z','median_z','q75_z','mixed_context_baseline_median_z'})];
end
save(fullfile(output_dir, [subjname, '_', measurelabel{measureidx}, '_stats.mat']), ...
    'subject_id','measureidx','measurelabel','input_file','metadata','seed', ...
    'rng_state','nrep','nperm','v','vperm','corr_c','rcross','rperm', ...
    'z_values','box_labels','bsl');
fprintf('%s %s: 800 trials, %d missing ratings.\n', subjname, ...
    measurelabel{measureidx}, sum(isnan(subject_trials{:,measurelabel{measureidx}})));
end
end
writetable(summary, fullfile(output_dir, 'boxplot_descriptive_stats.csv'));
writetable(condition_summary, fullfile(output_dir, 'rating_stats_by_subject_odor_context.csv'));
writetable(trial_table, fullfile(output_dir, 'trial_ratings.csv'));
save(fullfile(output_dir, 'behavior_summary.mat'), ...
    'summary','condition_summary','trial_table','seed');
fprintf('Saved behavior figures and statistics to %s\n', output_dir);
