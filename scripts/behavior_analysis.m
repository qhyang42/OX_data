%% stats on behavior scatter plot 

%%%% plot a trace for a set of dots within context. 
%%%% plot traces for 1000 sets of random dots 
%%%% get a p value for context modulation of valence

%%% use behavior.mat for each participant 

%% 
% wkdir =  '/Volumes/ExtremeSSD/OX_DATA/behavior'; 
wkdir =     '/Users/qhyang/Desktop/OX_DATA/behavior'; 

% subjname = 'subj_4';
SUBJNAMES = {'subj_1', 'subj_2', 'subj_3', 'subj_4', 'subj_5', 'subj_6'}; 
sidx = 6; 
measureidx = 1; % 1 for valence, 2 for intensity

% for sidx = 1:4
% for measureidx = 1:2
    subjname = SUBJNAMES{sidx}; 
    load(fullfile( wkdir, subjname, 'behavior.mat')); 
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

%% deal with empty values -- this is a temporary measure for not enough data. 
%%%% using control value from the same odor for the missing part 

%%% for subj_3 
% data(2, 3) = data(2, 1); 
% data(12, 4) = data(12, 1); 
% odoridx = 1:20; 

%%%% pick the odors that have been used in all conditions
% %%%% for subj1 2 and 4
% if strcmp(subjname, 'subj_3')
%     %%% for subj_3
%     data(2, 3) = data(2, 1);
%     data(12, 4) = data(12, 1);
%     odoridx = 1:20;
% 
% else
%     emptycells = cellfun(@isempty, data);
%     rowswempty = any(emptycells, 2);
%     data_clean = data(~rowswempty, :);
%     odorlabels = odorlabels(~rowswempty);
%     odoridx = 1:20;
%     odoridx = odoridx(~rowswempty);
% 
%     data = data_clean;
% 
% end 

%% calculate inter-odor trajectroy similarity for each category  
odoridx = 1:20; 
nrep = 100; 
v = zeros(nrep, length(odoridx), max(catidx)); 
for repidx = 1:nrep
    for c = 1: max(catidx)

        for n = 1: length(odoridx)
            ratings = data{n, c};
            ratings = ratings(:, measureidx); 
            v(repidx, n, c) = ratings(randperm(length(ratings), 1));
        end
    end

end

%%% 
corr_c = zeros(nrep*(nrep-1), max(catidx)); % all the off diagonal correlation coefficient  
for c = 1:max(catidx)
    v_c = squeeze(v(:, :, c)); 
    r = corrcoef(v_c'); 
    r = r(r~=1); 
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
        vperm(repidx, n) = ratings(randperm(length(ratings), 1)); 
    end 
end 

rperm = corrcoef(vperm'); 
rperm = rperm(rperm ~= 1); 
rperm = sort(rperm); 
%% add mean rperm line to the box plot 
bsl = median(atanh(rperm)); 
subjlabel = {'AS', 'LS', 'JN', 'RR', 'BN', 'VS'};
measurelabel = {'valence', 'intensity'}; 
figure; 
boxh = boxplot([atanh(corr_c), atanh(rcross)]); 
hold on 
ph = plot([0:6], bsl*ones(1,7), 'k', 'LineStyle', '-.'); 
ylabel('intertrial correlation'); 
xlabel('context'); 
xticklabels([catlabels, {'cross condition'}]); 
% title([subjlabel{sidx}, ', ', measurelabel{measureidx}]);
title('subj_6'); 
set(boxh, 'LineWidth' , 1.5); 
set(ph, 'LineWidth', 1.5); 
% ylim([-0.5, 1]); 
% end 
% end 

%% add rperm as box plot 
bsl = median(atanh(rperm)); 
rperm_box = atanh(rperm); 
subjlabel = {'AS', 'LS', 'JN', 'RR', 'BN', 'VS'};
measurelabel = {'valence', 'intensity'}; 
figure; 
boxh = boxplot([atanh(corr_c), atanh(rcross), rperm_box]); 
hold on 
ph = plot([0:6], bsl*ones(1,7), 'k', 'LineStyle', '-.'); 
ylabel('intertrial correlation'); 
xlabel('context'); 
xticklabels([catlabels, {'cross condition', 'permuted'}]); 
title([subjlabel{sidx}, ', ', measurelabel{measureidx}]);
set(boxh, 'LineWidth' , 1.5); 
set(ph, 'LineWidth', 1.5); 
 