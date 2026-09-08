%% Paths and acquisition layout
wkdir = '/Users/qhyang/Desktop/OX_DATA/behavior';
cuelistdir = '/Users/qhyang/Desktop/OX_DATA/cuelist';
SUBJNAMES = {'subj_2', 'subj_3', 'subj_4', 'subj_5', 'subj_6'};
SESSION_RUN_COUNTS = { ...
    [10, 10, 10, 5, 5, 3, 5, 5, 5, 6, 5, 6, 5], ... % subj_2
    [10, 8, 8, 8, 8, 6, 8, 8, 7, 9], ...             % subj_3
    [8, 7, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5], ... % subj_4
    [5, 5, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 6, 5, 5], ... % subj_5
    [8, 8, 8, 8, 8, 6, 8, 9, 9, 8]};                 % subj_6
odorlabels = {'anise', 'brownie', 'orange', '1-nonanol', 'chicken', 'rose', 'dish soap', 'coffee', 'coconut', 'benzaldehyde', ...
    'Air', 'Ethylbenzene', 'banana', 'Whiskey', 'Diethylpyrazine', 'Garlic', 'Methyl salicylate', 'Tea tree', 'Parmesan', 'Salsa'}; 
c = 756; % this is the screen center and default x position 

%% Extract and save behavior data for subjects 2-6
for sidx = 1:length(SUBJNAMES)
    subjname = SUBJNAMES{sidx};
    run_counts = SESSION_RUN_COUNTS{sidx};
    valence_all = [];
    intensity_all = [];
    odor = [];
    category = {};
    nruns = 0;

    for sesidx = 1:length(run_counts)
        for runidx = 1:run_counts(sesidx)
            result_file = fullfile(wkdir, subjname, ...
                [subjname, '_session', num2str(sesidx), '_run', num2str(runidx), '_results.mat']);
            run_data = load(result_file);
            nruns = nruns + 1;

            valence = [];
            intensity = [];

            if isfield(run_data, 'datalabel')
                % read cell 2 and 3 as valence and intensity rating
                for n = 1:length(run_data.outMat)
                    if isempty(run_data.outMat{n}{2})
                        valence(n) = NaN;
                    else
                        valence(n) = run_data.outMat{n}{2};
                    end
                    if isempty(run_data.outMat{n}{3})
                        intensity(n) = NaN;
                    else
                        intensity(n) = run_data.outMat{n}{3};
                    end
                end

            else
                % read 4 and 5
                for n = 1:length(run_data.outMat)
                    if isempty(run_data.outMat{n}{4})
                        valence(n) = NaN;
                    else
                        valence(n) = run_data.outMat{n}{4};
                    end
                    if isempty(run_data.outMat{n}{5})
                        intensity(n) = NaN;
                    else
                        intensity(n) = run_data.outMat{n}{5};
                    end
                end
            end

            valence_all = [valence_all, valence];
            intensity_all = [intensity_all, intensity];

            cue_file = fullfile(cuelistdir, subjname, ...
                ['session', num2str(sesidx)], ...
                ['cuelist_sess', num2str(sesidx), '_run', num2str(runidx), '.mat']);
            cue_data = load(cue_file);
            odor = [odor; cue_data.cuelist.odor];
            category = [category, cue_data.cuelist.category];

            fprintf('%s: loaded session %d run %d.\n', subjname, sesidx, runidx);
        end
    end

    valence_all = valence_all - c;
    intensity_all = intensity_all - c; % center results

    save(fullfile(wkdir, subjname, 'behavior.mat'), ...
        'intensity_all', 'valence_all', 'odor', 'category', 'odorlabels');
    fprintf('%s: saved %d runs and %d trials.\n', subjname, nruns, length(valence_all));
end

%% Exploratory plotting (run this block separately after extraction)
plot_subjname = 'subj_6';
load(fullfile(wkdir, plot_subjname, 'behavior.mat'));

figure;
hold on
plot(valence_all);
plot(intensity_all);
legend({'valence', 'intensity'});
xlabel('n trial');
ylabel('rating');

%% Exploratory plot: results by odor

valence_c = []; 
valence_f = []; 
valence_l = []; 
valence_p = []; 

intensity_c = []; 
intensity_f = []; 
intensity_p = []; 
intensity_l = []; 

odor_c = []; 
odor_f = []; 
odor_p = []; 
odor_l = []; 

for n = 1: length(odor)
    this_cat = category{n}; 
    switch this_cat 
        case 'CONTROL'
            valence_c = [valence_c, valence_all(n)]; 
            intensity_c = [intensity_c, intensity_all(n)]; 
            odor_c = [odor_c, odor(n)]; 
        case 'FOOD'
            valence_f = [valence_f, valence_all(n)]; 
            intensity_f = [intensity_f, intensity_all(n)]; 
            odor_f = [odor_f, odor(n)]; 

        case 'PERSON'
            valence_p = [valence_p, valence_all(n)]; 
            intensity_p = [intensity_p, intensity_all(n)]; 
            odor_p = [odor_p, odor(n)]; 
        case 'LOCATION'
            valence_l = [valence_l, valence_all(n)]; 
            intensity_l = [intensity_l, intensity_all(n)]; 
            odor_l = [odor_l, odor(n)]; 
    end 
end 

%% plot 

figure; 
hold on 
plot(1*ones(length(valence_c), 1), valence_c, '.'); 
plot(2*ones(length(valence_f), 1), valence_f, '.'); 
plot(3*ones(length(valence_p), 1), valence_p, '.'); 
plot(4*ones(length(valence_l), 1), valence_l, '.'); 

%% plot rating for each odor in each category 


plotcolors = [[0, 0.4470, 0.7410];...	          
          	[0.8500, 0.3250, 0.0980];...	       
          	[0.9290, 0.6940, 0.1250];...	          
          	[0.4940, 0.1840, 0.5560];...	          
          	[0.4660, 0.6740, 0.1880];...	          
          	[0.3010, 0.7450, 0.9330];...	          
          	[0.6350, 0.0780, 0.1840]]; 	

figure; 
hold on;
for n = 1:length(odor)
    % get odor 
    this_odor = odor(n); 
    j = 0.2*(rand(1)-1); 
    thisval = valence_all(n);
    this_cat = category{n};
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

    scatter(this_odor+j, thisval, 50, plotcolors(x, :), "filled"); 

end 
scatter (max(odor)+1, -100, 50, plotcolors(1, :), "filled"); 
text(max(odor)+1.5, -100, 'CONTROL'); 

scatter (max(odor)+1, 0, 50, plotcolors(2, :), "filled"); 
text(max(odor)+1.5, 0, 'FOOD'); 

scatter (max(odor)+1, 100, 50, plotcolors(3, :), "filled"); 
text(max(odor)+1.5, 100, 'PERSON'); 

scatter (max(odor)+1, 200, 50, plotcolors(4, :), "filled"); 
text(max(odor)+1.5, 200, 'LOCATION'); 

%%% line zero 
plot([0:max(odor)], zeros(1, length([0:max(odor)])), 'k'); 

xlim([0, max(odor)+3]); 
xticks([1:max(odor)]); 
xlabel('odor'); 
ylabel('valence rating'); 
xticklabels(odorlabels); 

title(plot_subjname);
%% NOTES 
% quality check for subject 1
% valence and intensity ratings are correlated. shouldn't be ? 

