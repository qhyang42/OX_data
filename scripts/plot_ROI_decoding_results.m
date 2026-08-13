%% ROI decoding results 

%% load data 
wkdir =  '/Users/qhyang/Desktop/OX_DATA/MRI/'; 
subjname = {'subj_1','subj_2', 'subj_3', 'subj_4', 'subj_5', 'subj_6' }; 
context_result = [];
odor_result = []; 
for subjidx = 2:6
    context_result{subjidx} = load(fullfile(wkdir, subjname{subjidx}, 'nifti', 'single_trial_by_category',...
        'roi_decoding_context_template_loro', ['context_template_subj', num2str(subjidx), '_loro_results.mat'])); 

  odor_result{subjidx} = load(fullfile(wkdir, subjname{subjidx}, 'nifti', 'single_trial_by_category',...
        'roi_decoding_odor_template_loro', ['odor_template_subj', num2str(subjidx), '_loro_results.mat'])); 
end 

%% plot dots per subject 

figure;
hold on
for roiidx = 1:66
    x = roiidx + 0.1*randn([5, 1]);
    y = []; 
    for subjidx = 2:6
        y = [y, odor_result{subjidx}.results.summary.accuracy(roiidx)]; 
    end 
    scatter(x, y); 

end 
xticks([1:1:66]); 
xticklabels(odor_result{2}.results.roi_metadata.roi_name); 
plot([1:66], 0.05*ones([1, 66]), 'k-'); 
title('odor'); 

figure; 
hold on
for roiidx = 1:66
    x = roiidx + 0.1*randn([5, 1]);
    y = []; 
    for subjidx = 2:6
        y = [y, context_result{subjidx}.results.summary.accuracy(roiidx)]; 
    end 
    scatter(x, y); 

end 
xticks([1:1:66]); 
xticklabels(context_result{2}.results.roi_metadata.roi_name); 
plot([1:66], 0.25*ones([1, 66]), 'k-'); 
title('context'); 

figure; 
hold on
for roiidx = 1:66
    x = roiidx + 0.1*randn([5, 1]);
    y = []; 
    for subjidx = 2:6
        y = [y, context_result{subjidx}.results.summary.n_gm_overlap(roiidx)]; 
    end 
    scatter(x, y); 

end 
xticks([1:1:66]); 
xticklabels(context_result{2}.results.roi_metadata.roi_name); 
plot([1:66], 0.25*ones([1, 66]), 'k-'); 
title('Voxel count'); 

%% bilateral ROIs only 
bilateralidx = 1:3:66; 

figure;
hold on
i = 0; 
for roiidx = bilateralidx
    i = i+1; 
    x = i + 0.1*randn([5, 1]);
    y = []; 
    for subjidx = 2:6
        y = [y, odor_result{subjidx}.results.summary.accuracy(roiidx)]; 
    end 
    scatter(x, y); 

end 
xticks([1:1:length(bilateralidx)]); 
xticklabels(odor_result{2}.results.roi_metadata.roi_name(bilateralidx)); 
plot([1:length(bilateralidx)], 0.05*ones([1, length(bilateralidx)]), 'k-'); 
title('odor'); 

figure; 
hold on
i = 0; 
for roiidx = bilateralidx
    i = i+1; 
    x = i + 0.1*randn([5, 1]);
    y = []; 
    for subjidx = 2:6
        y = [y, context_result{subjidx}.results.summary.accuracy(roiidx)]; 
    end 
    scatter(x, y); 

end 
xticks([1:1:length(bilateralidx)]); 
xticklabels(context_result{2}.results.roi_metadata.roi_name(bilateralidx)); 
plot([1:length(bilateralidx)], 0.25*ones([1, length(bilateralidx)]), 'k-'); 
title('context'); 

% figure; 
% hold on
% i = 0; 
% for roiidx = bilateralidx
%     i = i+1; 
%     x = i + 0.1*randn([5, 1]);
%     y = []; 
%     for subjidx = 2:6
%         y = [y, context_result{subjidx}.results.summary.n_gm_overlap(roiidx)]; 
%     end 
%     scatter(x, y); 
% 
% end 
% xticks([1:1:length(bilateralidx)]); 
% xticklabels(context_result{2}.results.roi_metadata.roi_name(bilateralidx)); 
% plot([1:length(bilateralidx)], 0.25*ones([1, length(bilateralidx)]), 'k-'); 
% title('v count'); 

%% Correct recall for each context in each bilateral ROI
subjectidx = 2:6;
context_names = string(context_result{subjectidx(1)}.results.classes.values);
nContexts = numel(context_names);

% Recall is the diagonal of the row-normalized confusion matrix:
% P(predicted context = c | true context = c).
context_recall = nan(numel(subjectidx), nContexts, numel(bilateralidx));
for s = 1:numel(subjectidx)
    confusion = context_result{subjectidx(s)}.results.confusion_matrices_normalized;
    for r = 1:numel(bilateralidx)
        context_recall(s, :, r) = diag(squeeze(confusion(bilateralidx(r), :, :)));
    end
end

figure;
tl = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
for r = 1:numel(bilateralidx)
    nexttile;
    recall = context_recall(:, :, r);
    hold on;
    scatter(repmat(1:nContexts, numel(subjectidx), 1), recall, 18, ...
        'filled', 'MarkerFaceAlpha', 0.45);
    plot(1:nContexts, mean(recall, 1, 'omitnan'), '-ok', ...
        'LineWidth', 1.2, 'MarkerFaceColor', 'k');
    yline(1/nContexts, '--');
    xticks(1:nContexts);
    xticklabels(context_names);
    xtickangle(30);
    % ylim([0 1]);
    title(context_result{subjectidx(1)}.results.roi_metadata.roi_name(bilateralidx(r)), ...
        'Interpreter', 'none');
end
xlabel(tl, 'Context');
ylabel(tl, 'Correct recall');
sgtitle('Context recall by bilateral ROI');



%% 
