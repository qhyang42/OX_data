%% Group consistency and same-odor CONTROL displacement (descriptive).
% Extend behavior_analysis.m using its saved resamples; validate against inputs.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'utils','OX_utilities'));
source = fullfile(root,'results','behavior');
out = fullfile(source,'group_consistency_and_trial_absolute_displacement_20260915');
if isfolder(out), out = [out '_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]; end
mkdir(out);
for d = ["consistency/group","displacement/individual","displacement/group","tables"]
    mkdir(fullfile(out,d));
end
subjects = 2:6; measures = ["pleasantness","intensity"];
contexts = ["CONTROL","FOOD","PERSON","LOCATION"];
labels = [contexts,"Cross-context","Mixed-context"];
colors = [0.45 .49 .53; .90 .47 .20; .31 .55 .76; .43 .65 .43; .65 .48 .70; .65 .65 .65];
trials = table(); conditions = table(); displacements = table();
subject_summary = table(); group_summary = table(); consistency = table();
means = zeros(5,2,20,4); counts = means; zmedian = zeros(5,2,6);
for s = 1:5
    id = subjects(s); input = fullfile(root,'behavior',sprintf('subj_%d',id),'behavior.mat');
    b = load(input);
    if s == 1, odorlabels = string(b.odorlabels); else, assert(isequal(string(b.odorlabels),odorlabels)); end
    meta = OX_load_trial_metadata(id);
    assert(isequal(double(b.odor(:)),meta.odor));
    assert(isequal(upper(strtrim(string(b.category(:)))),meta.context));
    assert(all(accumarray(meta.run_id,1)==10));
    ratings = [double(b.valence_all(:)),double(b.intensity_all(:))];
    assert(isequal(size(ratings),[800 2]) && ~any(isinf(ratings(:))));
    trials = [trials; addvars(meta,repmat(id,800,1),ratings(:,1),ratings(:,2), ...
        'NewVariableNames',{'subject_id','pleasantness','intensity'})];
    for m = 1:2
        saved = load(fullfile(source,sprintf('subj_%d_%s_stats.mat',id,measures(m))));
        assert(isequaln(saved.metadata,meta) && all(isfinite(saved.z_values(:))));
        assert(isequal(string(saved.box_labels(1:4)),contexts));
        zmedian(s,m,:) = median(saved.z_values,1);
        for c = 1:6
            consistency = [consistency; table(id,measures(m),labels(c),median(saved.z_values(:,c)), ...
                mean(saved.z_values(:,c)),'VariableNames',{'subject_id','measure','category','median_z','mean_z'})];
        end
        for o = 1:20
            for c = 1:4
                x = ratings(meta.odor==o & meta.context==contexts(c),m);
                assert(~isempty(x) && any(isfinite(x)));
                means(s,m,o,c) = mean(x,'omitnan'); counts(s,m,o,c) = sum(isfinite(x));
                conditions = [conditions; table(id,measures(m),o,contexts(c),numel(x),sum(isfinite(x)), ...
                    sum(isnan(x)),mean(x,'omitnan'),'VariableNames', ...
                    {'subject_id','measure','odor','context','n_trials','n_valid','n_missing','mean_rating'})];
            end
        end
    end
    assert(all(arrayfun(@(c)sum(meta.context==c),contexts)==200));
end
% Ensure saved consistency inputs still match current raw ratings, including NaNs.
prior = load(fullfile(source,'behavior_summary.mat'),'trial_table');
assert(isequaln(prior.trial_table,trials),'Saved consistency input changed: rerun behavior_analysis first.');
% Each semantic trial is referenced to its subject/odor CONTROL mean.
assert(height(trials)==4000 && height(conditions)==800);
participant_means = zeros(5,2,3);
for s = 1:5
    for m = 1:2
        rows = trials(trials.subject_id==subjects(s) & trials.context~="CONTROL",:);
        rating = rows{:,measures(m)};
        baseline = reshape(means(s,m,rows.odor,1),[],1);
        baseline_n = reshape(counts(s,m,rows.odor,1),[],1);
        absolute_displacement = abs(rating-baseline);
        rows = addvars(rows,repmat(measures(m),height(rows),1),rating,baseline,baseline_n, ...
            absolute_displacement,'NewVariableNames', ...
            {'measure','rating','control_mean','control_n_valid','absolute_displacement'});
        displacements = [displacements; rows];
        for c = 1:3
            x = absolute_displacement(rows.context==contexts(c+1));
            assert(numel(x)==200 && any(isfinite(x)));
            participant_means(s,m,c)=mean(x,'omitnan');
            subject_summary = [subject_summary; table(subjects(s),measures(m),contexts(c+1), ...
                numel(x),sum(isfinite(x)),sum(isnan(x)),mean(x,'omitnan'),median(x,'omitnan'),std(x,'omitnan'), ...
                'VariableNames',{'subject_id','measure','context','n_trials','n_valid','n_missing', ...
                'mean_absolute_displacement','median_absolute_displacement','sd_absolute_displacement'})];
        end
    end
end
group_ylim = [0 ceil(max(participant_means(:))*1.15/10)*10];
consistency_ylim = [floor(min(zmedian(:))*10)/10 ceil(max(zmedian(:))*10)/10];
for m = 1:2
    y = squeeze(zmedian(:,m,:));
    f = figure('Visible','off','Color','w','Position',[100 100 1150 630]);
    ax = axes(f); drawboxes(ax,y,labels,colors,true);
    ylim(ax,consistency_ylim);
    ylabel(ax,'Odor-profile consistency (Fisher z)');
    title(ax,upperFirst(measures(m)) + " | group consistency");
    subtitle(ax,'Dots: participant medians of saved resamples; diamonds: equal-weight mean (n = 5)');
    yline(ax,mean(y(:,6)),'--','Mean mixed-context reference','LabelHorizontalAlignment','left');
    saveplot(f,fullfile(out,'consistency','group',measures(m)+'_group_consistency'));
    for c = 1:6
        group_summary = [group_summary; table(measures(m),labels(c),mean(y(:,c)),std(y(:,c)), ...
            'VariableNames',{'measure','category','mean_participant_median_z','sd_participant_median_z'})];
    end
    for s = 1:5
        y = nan(200,3);
        for c = 1:3
            mask = displacements.subject_id==subjects(s) & displacements.measure==measures(m) & ...
                displacements.context==contexts(c+1);
            y(:,c)=displacements.absolute_displacement(mask);
        end
        stem = string(fullfile(out,'displacement','individual',sprintf('subj_%d_%s',subjects(s),measures(m))));
        displacementplot(y,contexts(2:4),colors(2:4,:),false, ...
            sprintf('Subject %d | %s',subjects(s),measures(m)),stem,[0 ceil(max(y(:))*1.1/10)*10]);
    end
    y = squeeze(participant_means(:,m,:));
    stem = fullfile(out,'displacement','group',measures(m)+'_group');
    displacementplot(y,contexts(2:4),colors(2:4,:),true, ...
        upperFirst(measures(m))+' | group absolute displacement',stem,group_ylim);
end
writetable(trials,fullfile(out,'tables','trial_ratings.csv'));
writetable(conditions,fullfile(out,'tables','condition_ratings.csv'));
writetable(consistency,fullfile(out,'tables','consistency_participant_summary.csv'));
writetable(group_summary,fullfile(out,'tables','consistency_group_summary.csv'));
writetable(displacements,fullfile(out,'tables','absolute_displacement_by_trial.csv'));
writetable(subject_summary,fullfile(out,'tables','displacement_participant_summary.csv'));
group_displacement = groupsummary(subject_summary,{'measure','context'},{'mean','std'}, ...
    {'mean_absolute_displacement'});
writetable(group_displacement,fullfile(out,'tables','displacement_group_summary.csv'));
save(fullfile(out,'behavior_results.mat'),'trials','conditions','consistency','group_summary', ...
    'displacements','subject_summary','group_displacement','means','counts','participant_means','group_ylim','consistency_ylim','zmedian','subjects','contexts','measures','source');
copyfile(fullfile(root,'scripts','behavior_group_and_displacement_methods.md'),fullfile(out,'README.md'));
copyfile(mfilename('fullpath')+".m",fullfile(out,'behavior_group_and_displacement_source.m'));
fprintf('OUTPUT: %s\n',out);
disp(group_displacement);

function drawboxes(ax,y,labels,colors,participant)
hold(ax,'on');
for c = 1:size(y,2)
    boxchart(ax,repmat(c,size(y,1),1),y(:,c),'BoxFaceColor',colors(c,:), ...
        'BoxFaceAlpha',.23,'MarkerStyle','none','BoxWidth',.5);
    offsets = linspace(-.14,.14,size(y,1))';
    if participant
        for s = 1:size(y,1)
            scatter(ax,c+offsets(s),y(s,c),55,colors(c,:),'filled','MarkerEdgeColor','w');
            text(ax,c+offsets(s)+.025,y(s,c),sprintf(' %d',s+1),'FontSize',9,'Color',[.25 .25 .25]);
        end
    else
        scatter(ax,c+offsets,y(:,c),30,colors(c,:),'filled','MarkerFaceAlpha',.75);
    end
    scatter(ax,c,mean(y(:,c),'omitnan'),85,'k','d','filled');
end
set(ax,'XTick',1:numel(labels),'XTickLabel',labels,'FontSize',12,'LineWidth',1,'Box','off');
xlim(ax,[.5 numel(labels)+.5]); grid(ax,'on'); ax.XGrid='off'; ax.GridAlpha=.12;
end
function displacementplot(y,labels,colors,participant,heading,stem,limits)
f=figure('Visible','off','Color','w','Position',[100 100 1100 650]);
ax=axes(f);drawboxes(ax,y,labels,colors,participant);
ylabel(ax,'|Trial rating - same-odor CONTROL mean| (rating units)');
title(ax,heading);ylim(ax,limits);
if participant
 subtitle(ax,'Dots: participant means across valid trials (labels = subject ID); diamonds: group mean; n = 5');
else
 subtitle(ax,'Dots: individual valid trials; diamonds: context mean; missing ratings omitted');
end
saveplot(f,stem);
end
function saveplot(f,stem)
stem=string(stem);
exportgraphics(f,stem+'.png','Resolution',180);exportgraphics(f,stem+'.pdf','ContentType','vector');
savefig(f,stem+'.fig');close(f);
end
function s=upperFirst(s)
s=char(s);s(1)=upper(s(1));s=string(s);
end
