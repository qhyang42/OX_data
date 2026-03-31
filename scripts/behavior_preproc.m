%%% subj_1: AS 
%%% subj_2: LS 
%%% subj_3: JN 


%% 
wkdir =     '/Users/qhyang/Desktop/OX_DATA/behavior'; 
subjname = 'subj_6';
cuelistdir = '/Users/qhyang/Desktop/OX_DATA/cuelist'; 
% sesidx = 1; 
% runidx = 1; 
odorlabels = {'anise', 'brownie', 'orange', '1-nonanol', 'chicken', 'rose', 'dish soap', 'coffee', 'coconut', 'benzaldehyde', ...
    'Air', 'Ethylbenzene', 'banana', 'Whiskey', 'Diethylpyrazine', 'Garlic', 'Methyl salicylate', 'Tea tree', 'Parmesan', 'Salsa'}; 
c = 756; % this is the screen center and default x position 

%% read behavior data  
valence_all = []; 
intensity_all = []; 
nruns = 0; 

for sesidx = 1: 16
    for runidx = 1:10
        try
            load(fullfile(wkdir, subjname, [subjname, '_session', num2str(sesidx), '_run', num2str(runidx), '_results.mat']));
            fprintf('load session %d run %d.   \n', sesidx, runidx); 
            nruns = nruns+1; 
        catch
            continue
        end
        varlist = who();
        valence = [];
        intensity = [];

        if sum(ismember(varlist, 'datalabel'))>0
            % read cell 2 and 3 as valence and intensity rating
            for n = 1: length(outMat)
                if isempty(outMat{n}{2})
                    valence(n) = c;
                else
                    valence(n) = outMat{n}{2};
                end
                if isempty(outMat{n}{3})
                    intensity(n) = c;
                else
                    intensity(n) = outMat{n}{3};
                end
            end

        else
            % read 4 and 5
            for n = 1: length(outMat)
                if isempty(outMat{n}{4})
                    valence(n) = c;
                else
                    valence(n) = outMat{n}{4};
                end
                if isempty(outMat{n}{5})
                    intensity(n) = c;
                else
                    intensity(n) = outMat{n}{5};
                end
            end


        end

        valence_all = [valence_all, valence]; 
        intensity_all = [intensity_all, intensity]; 

    end
end


valence_all = valence_all - c; 
intensity_all = intensity_all - c; % center results  

%% plot 
figure; 
hold on 
plot(valence_all); 
plot(intensity_all); 
legend({'valence', 'intensity'}); 
xlabel('n trial'); 
ylabel('rating'); 

%% read cuelist 
odor = []; 
category = {}; 
for sesidx = 1: 16
    for runidx = 1:10
        try
            load(fullfile(cuelistdir, subjname, ['session', num2str(sesidx)], ['cuelist_sess', num2str(sesidx), '_run', num2str(runidx), '.mat']));
            fprintf('load session %d run %d.   \n', sesidx, runidx); 

        catch
            continue
        end

    odor = [odor; cuelist.odor]; 
    category = [category, cuelist.category]; 
    end 
end     

% odor = odor(1:160); 
% category = category(1:160);


%% plot results by odor 

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

title(subjname); 

%% 
save(fullfile(wkdir, subjname, 'behavior.mat'), 'intensity_all', 'valence_all', 'odor', 'category', 'odorlabels');
%% NOTES 
% quality check for subject 1
% valence and intensity ratings are correlated. shouldn't be ? 



