%% ROI decoding results 

%% load data 
wkdir =  '/Users/qhyang/Desktop/OX_DATA/MRI/'; 
subjname = {'subj_1','subj_2', 'subj_3', 'subj_4', 'subj_5', 'subj_6' }; 
context_result = [];
odor_result = []; 
for subjidx = 2:6
    context_result_sniff{subjidx} = load(fullfile(wkdir, subjname{subjidx}, 'nifti', 'sniff_single_trial_by_category_physio',...
        'sniff_roi_decoding_primary_context_template_loro_physio', ['context_template_subj', num2str(subjidx), '_loro_results.mat']));

  % odor_result{subjidx} = load(fullfile(wkdir, subjname{subjidx}, 'nifti', 'sniff_single_trial_by_category_physio',...
  %       'roi_decoding_odor_template_loro_physio', ['odor_template_subj', num2str(subjidx), '_loro_results.mat']));

  context_result_countdown{subjidx} = load(fullfile(wkdir, subjname{subjidx}, 'nifti', 'countdown_single_trial_by_category_physio',...
        'countdown_roi_decoding_primary_context_template_loro_physio', ['context_template_subj', num2str(subjidx), '_loro_results.mat']));

end 

%% plot dots per subject 

figure;
hold on
for roiidx = 1:12
    x = roiidx + 0.1*randn([5, 1]);
    y = []; 
    for subjidx = 2:6
        y = [y, context_result_sniff{subjidx}.results.summary.accuracy(roiidx)]; 
    end 
    scatter(x, y); 

end 
xticks([1:1:12]); 
xticklabels(context_result_sniff{2}.results.roi_metadata.roi_name); 
plot([1:12], 0.25*ones([1, 12]), 'k-'); 
title('sniff'); 

figure; 
hold on
for roiidx = 1:12
    x = roiidx + 0.1*randn([5, 1]);
    y = []; 
    for subjidx = 2:6
        y = [y, context_result_countdown{subjidx}.results.summary.accuracy(roiidx)]; 
    end 
    scatter(x, y); 

end 
xticks([1:1:12]); 
xticklabels(context_result_countdown{2}.results.roi_metadata.roi_name); 
plot([1:12], 0.25*ones([1, 12]), 'k-'); 
title('countdown'); 

% figure; 
% hold on
% for roiidx = 1:66
%     x = roiidx + 0.1*randn([5, 1]);
%     y = []; 
%     for subjidx = 2:6
%         y = [y, context_result{subjidx}.results.summary.n_gm_overlap(roiidx)]; 
%     end 
%     scatter(x, y); 
% 
% end 
% xticks([1:1:66]); 
% xticklabels(context_result{2}.results.roi_metadata.roi_name); 
% plot([1:66], 0.25*ones([1, 66]), 'k-'); 
% title('Voxel count'); 


%% Correct recall for each context in each bilateral ROI
bilateralidx = 1:12; 
subjectidx = 2:6;
context_names = string(context_result_sniff{subjectidx(1)}.results.classes.values);
nContexts = numel(context_names);

% Recall is the diagonal of the row-normalized confusion matrix:
% P(predicted context = c | true context = c).
context_recall = nan(numel(subjectidx), nContexts, numel(bilateralidx));
for s = 1:numel(subjectidx)
    confusion = context_result_sniff{subjectidx(s)}.results.confusion_matrices_normalized;
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
    title(context_result_sniff{subjectidx(1)}.results.roi_metadata.roi_name(bilateralidx(r)), ...
        'Interpreter', 'none');
end
xlabel(tl, 'Context');
ylabel(tl, 'Correct recall');
sgtitle('Context recall by bilateral ROI');

%% Per-subject context decoding accuracy in bilateral ROIs


context_result = context_result_sniff; 

subjectidx = 2:6;
nSubjects = numel(subjectidx);
nBilateral = numel(bilateralidx);
bilateral_accuracy = nan(nSubjects, nBilateral);
bilateral_significant = false(nSubjects, nBilateral);

for s = 1:nSubjects
    summary = context_result{subjectidx(s)}.results.summary;
    bilateral_accuracy(s, :) = summary.accuracy(bilateralidx);
    bilateral_significant(s, :) = summary.sig_fdr05(bilateralidx); % FDR q < 0.05
end

% Fixed offsets keep each subject visible at every ROI.
x = repmat(1:nBilateral, nSubjects, 1) + ...
    repmat(linspace(-0.18, 0.18, nSubjects)', 1, nBilateral);

figure;
hold on;
h_nonsig = scatter(x(~bilateral_significant), ...
    bilateral_accuracy(~bilateral_significant), 32, [0.65 0.65 0.65], 'filled');
h_sig = scatter(x(bilateral_significant), ...
    bilateral_accuracy(bilateral_significant), 36, [0 0.45 0.74], 'filled');
yline(context_result{subjectidx(1)}.results.classes.chance, '--k', 'Chance');
xticks(1:nBilateral);
xticklabels(context_result{subjectidx(1)}.results.roi_metadata.roi_name(bilateralidx));
xtickangle(45);
xlim([0.5 nBilateral + 0.5]);
ylabel('Decoding accuracy');
title('Per-subject context decoding accuracy');
legend([h_sig h_nonsig], {'FDR q < 0.05', 'Not significant'}, ...
    'Location', 'best');


%% Per-subject odor decoding accuracy in bilateral ROIs
subjectidx = 2:6;
nSubjects = numel(subjectidx);
nBilateral = numel(bilateralidx);
bilateral_accuracy = nan(nSubjects, nBilateral);
bilateral_significant = false(nSubjects, nBilateral);

for s = 1:nSubjects
    summary = odor_result{subjectidx(s)}.results.summary;
    bilateral_accuracy(s, :) = summary.accuracy(bilateralidx);
    bilateral_significant(s, :) = summary.sig_fdr05(bilateralidx); % FDR q < 0.05
end

% Fixed offsets keep each subject visible at every ROI.
x = repmat(1:nBilateral, nSubjects, 1) + ...
    repmat(linspace(-0.18, 0.18, nSubjects)', 1, nBilateral);

figure;
hold on;
h_nonsig = scatter(x(~bilateral_significant), ...
    bilateral_accuracy(~bilateral_significant), 32, [0.65 0.65 0.65], 'filled');
h_sig = scatter(x(bilateral_significant), ...
    bilateral_accuracy(bilateral_significant), 36, [0 0.45 0.74], 'filled');
yline(odor_result{subjectidx(1)}.results.classes.chance, '--k', 'Chance');
xticks(1:nBilateral);
xticklabels(odor_result{subjectidx(1)}.results.roi_metadata.roi_name(bilateralidx));
xtickangle(45);
xlim([0.5 nBilateral + 0.5]);
ylabel('Decoding accuracy');
title('Per-subject odor decoding accuracy: bilateral ROIs');
legend([h_sig h_nonsig], {'FDR q < 0.05', 'Not significant'}, ...
    'Location', 'best');



%% Change in context decoding: odor-aligned minus countdown-aligned
odor_aligned_context = cell(size(subjname));
countdown_aligned_context = cell(size(subjname));

for subjidx = subjectidx
    result_name = ['context_template_subj', num2str(subjidx), '_loro_results.mat'];
    odor_aligned_context{subjidx} = load(fullfile(wkdir, subjname{subjidx}, ...
        'nifti', 'sniff_single_trial_by_category_physio', ...
        'roi_decoding_context_template_loro_physio', result_name));
    countdown_aligned_context{subjidx} = load(fullfile(wkdir, subjname{subjidx}, ...
        'nifti', 'countdown_single_trial_by_category_physio', ...
        'roi_decoding_context_template_loro_physio', result_name));

    assert(isequal(odor_aligned_context{subjidx}.results.roi_metadata.roi_name, ...
        countdown_aligned_context{subjidx}.results.roi_metadata.roi_name), ...
        'ROI ordering differs between alignments for subject %d.', subjidx);
end

accuracy_change = nan(nSubjects, nBilateral);
for s = 1:nSubjects
    subjidx = subjectidx(s);
    odor_accuracy = odor_aligned_context{subjidx}.results.summary.accuracy(bilateralidx);
    countdown_accuracy = countdown_aligned_context{subjidx}.results.summary.accuracy(bilateralidx);
    accuracy_change(s, :) = odor_accuracy - countdown_accuracy;
end

figure;
hold on;
subject_offsets = linspace(-0.18, 0.18, nSubjects);
x_subject = repmat(1:nBilateral, nSubjects, 1) + subject_offsets';
h_subject = scatter(x_subject(:), accuracy_change(:), 28, ...
    [0.78 0.78 0.78], 'filled');
h_mean = scatter(1:nBilateral, mean(accuracy_change, 1, 'omitnan'), ...
    46, [0 0.45 0.74], 'filled');
yline(0, '--k');
xticks(1:nBilateral);
xticklabels(odor_aligned_context{subjectidx(1)}.results.roi_metadata.roi_name(bilateralidx));
xtickangle(45);
xlim([0.5 nBilateral + 0.5]);
ylabel('\Delta accuracy (odor-aligned - countdown-aligned)');
title('Change in context decoding accuracy by bilateral ROI');
legend([h_subject h_mean], {'Subjects', 'Mean'}, 'Location', 'best');

%% 
