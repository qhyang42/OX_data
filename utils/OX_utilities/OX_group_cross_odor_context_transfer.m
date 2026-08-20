function group_results = OX_group_cross_odor_context_transfer(subject_ids, varargin)
%OX_GROUP_CROSS_ODOR_CONTEXT_TRANSFER Group inference for context transfer.
%
%   group_results = OX_group_cross_odor_context_transfer(subject_ids, ...)
%
% Overall evidence and accuracy use separate selected-ROI max-statistic
% families. Context-specific evidence and recall use separate selected
% ROI-by-context families. All tests are exact one-sided subject sign flips.

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'subject_ids', @(x) isnumeric(x) && isvector(x) && ...
    all(isfinite(x)) && all(x == round(x)));
addParameter(p, 'MRIRoot', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'ROISelection', 'old', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'OutputDir', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'MakePlots', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'SaveOutputs', true, @(x) islogical(x) && isscalar(x));
parse(p, subject_ids, varargin{:});
opts = p.Results;
roi_selection = OX_normalize_decoding_roi_selection(opts.ROISelection);
opts.ROISelection = roi_selection;
run_inference_self_tests();

subject_ids = subject_ids(:)';
n_subjects = numel(subject_ids);
assert(n_subjects >= 2, 'At least two subjects are required.');
mri_root = resolve_mri_root(opts.MRIRoot);
if strlength(string(opts.OutputDir)) == 0
    if strcmp(roi_selection, 'old')
        output_name = 'roi_cross_odor_context_template_loro';
    else
        output_name = sprintf('roi_%s_cross_odor_context_template_loro', roi_selection);
    end
    output_dir = fullfile(mri_root, 'group', output_name);
else
    output_dir = char(string(opts.OutputDir));
end

subject_results = cell(n_subjects, 1);
for subject_idx = 1:n_subjects
    subject_id = subject_ids(subject_idx);
    if strcmp(roi_selection, 'old')
        subject_output_name = 'roi_cross_odor_context_template_loro';
    else
        subject_output_name = sprintf('roi_%s_cross_odor_context_template_loro', roi_selection);
    end
    filename = fullfile(mri_root, sprintf('subj_%d', subject_id), 'nifti', ...
        'single_trial_by_category', subject_output_name, ...
        sprintf('cross_odor_context_template_subj%d_loro_results.mat', subject_id));
    assert(isfile(filename), 'Missing subject result: %s', filename);
    loaded = load(filename, 'results');
    subject_results{subject_idx} = loaded.results;
end

reference = subject_results{1};
assert_result_selection(reference, roi_selection, subject_ids(1));
roi_metadata = reference.roi_metadata;
roi_names = string(roi_metadata.roi_name);
n_rois = height(roi_metadata);
odor_values = reference.odor_values(:);
n_odors = numel(odor_values);
contexts = string(reference.analysis.semantic_context_order(:));
assert(isequal(contexts, ["PERSON"; "FOOD"; "LOCATION"]), ...
    'Unexpected semantic context ordering.');

overall_evidence = nan(n_subjects, n_rois);
overall_accuracy = nan(n_subjects, n_rois);
overall_evidence_null = nan(n_subjects, n_rois);
overall_accuracy_null = nan(n_subjects, n_rois);
by_odor_evidence = nan(n_subjects, n_rois, n_odors);
by_odor_accuracy = nan(n_subjects, n_rois, n_odors);
by_odor_evidence_null = nan(n_subjects, n_rois, n_odors);
by_odor_accuracy_null = nan(n_subjects, n_rois, n_odors);
by_context_evidence = nan(n_subjects, n_rois, 3);
by_context_recall = nan(n_subjects, n_rois, 3);
by_context_evidence_null = nan(n_subjects, n_rois, 3);
by_context_recall_null = nan(n_subjects, n_rois, 3);
confusion = nan(n_subjects, n_rois, 3, 3);
confusion_normalized = nan(n_subjects, n_rois, 3, 3);

for subject_idx = 1:n_subjects
    item = subject_results{subject_idx};
    assert_result_selection(item, roi_selection, subject_ids(subject_idx));
    assert(isequal(string(item.roi_metadata.roi_name), roi_names), ...
        'ROI ordering differs for subject %d.', subject_ids(subject_idx));
    assert(isequal(item.odor_values(:), odor_values), ...
        'Odor ordering differs for subject %d.', subject_ids(subject_idx));
    assert(isequal(string(item.analysis.semantic_context_order(:)), contexts), ...
        'Context ordering differs for subject %d.', subject_ids(subject_idx));
    assert(~isempty(item.null.overall_evidence), ...
        'Subject %d has no permutation null.', subject_ids(subject_idx));
    overall_evidence(subject_idx, :) = item.overall.evidence;
    overall_accuracy(subject_idx, :) = item.overall.accuracy;
    overall_evidence_null(subject_idx, :) = item.null.mean_overall_evidence;
    overall_accuracy_null(subject_idx, :) = item.null.mean_overall_accuracy;
    by_odor_evidence(subject_idx, :, :) = item.by_odor.evidence;
    by_odor_accuracy(subject_idx, :, :) = item.by_odor.accuracy;
    by_odor_evidence_null(subject_idx, :, :) = item.null.mean_odor_evidence;
    by_odor_accuracy_null(subject_idx, :, :) = item.null.mean_odor_accuracy;
    by_context_evidence(subject_idx, :, :) = item.by_context.evidence;
    by_context_recall(subject_idx, :, :) = item.by_context.recall;
    by_context_evidence_null(subject_idx, :, :) = item.null.mean_context_evidence;
    by_context_recall_null(subject_idx, :, :) = item.null.mean_context_recall;
    confusion(subject_idx, :, :, :) = item.confusion_matrices;
    confusion_normalized(subject_idx, :, :, :) = item.confusion_matrices_normalized;
end
assert(all(isfinite([overall_evidence, overall_accuracy]), 'all'), ...
    'All bilateral ROIs must have finite results in every subject.');

overall_evidence_centered = overall_evidence - overall_evidence_null;
overall_accuracy_centered = overall_accuracy - overall_accuracy_null;
by_odor_evidence_centered = by_odor_evidence - by_odor_evidence_null;
by_odor_accuracy_centered = by_odor_accuracy - by_odor_accuracy_null;
by_context_evidence_centered = by_context_evidence - by_context_evidence_null;
by_context_recall_centered = by_context_recall - by_context_recall_null;

overall_evidence_inference = exact_sign_flip(overall_evidence_centered);
overall_accuracy_inference = exact_sign_flip(overall_accuracy_centered);
context_evidence_inference = exact_sign_flip( ...
    reshape(by_context_evidence_centered, n_subjects, []));
context_recall_inference = exact_sign_flip( ...
    reshape(by_context_recall_centered, n_subjects, []));

overall_table = make_overall_table(roi_metadata, n_subjects, ...
    overall_evidence, overall_evidence_null, overall_evidence_inference, ...
    overall_accuracy, overall_accuracy_null, overall_accuracy_inference);
context_table = make_context_table(roi_metadata, contexts, n_subjects, ...
    by_context_evidence, by_context_evidence_null, context_evidence_inference, ...
    by_context_recall, by_context_recall_null, context_recall_inference);
odor_table = make_odor_table(roi_metadata, odor_values, n_subjects, ...
    by_odor_evidence, by_odor_evidence_null, by_odor_accuracy, by_odor_accuracy_null);
subject_table = make_subject_table(subject_ids, roi_names, overall_evidence, ...
    overall_evidence_null, overall_accuracy, overall_accuracy_null, ...
    by_context_evidence, by_context_recall, contexts);
subject_odor_table = make_subject_odor_table(subject_ids, roi_names, odor_values, ...
    by_odor_evidence, by_odor_evidence_null, by_odor_accuracy, by_odor_accuracy_null);

group_results = struct();
group_results.analysis = struct( ...
    'name', 'cross_odor_generalization_of_semantic_context', ...
    'roi_selection', roi_selection, ...
    'subjects', subject_ids, 'n_subjects', n_subjects, ...
    'semantic_context_order', {cellstr(contexts)}, ...
    'odor_values', odor_values, ...
    'primary_alternative', 'cross-odor context evidence > empirical null', ...
    'secondary_alternative', 'cross-odor three-way accuracy > empirical null', ...
    'inference', 'exact one-sided subject sign flips', ...
    'n_exact_sign_flips', 2^n_subjects, ...
    'minimum_attainable_p', 1/(2^n_subjects), ...
    'multiple_comparison_control', 'max-statistic FWE; BH-FDR also reported');
group_results.roi_metadata = roi_metadata;
group_results.subject_level = struct( ...
    'overall_evidence', overall_evidence, ...
    'overall_evidence_null_centered', overall_evidence_centered, ...
    'overall_accuracy', overall_accuracy, ...
    'overall_accuracy_null_centered', overall_accuracy_centered, ...
    'by_odor_evidence', by_odor_evidence, ...
    'by_odor_evidence_null_centered', by_odor_evidence_centered, ...
    'by_odor_accuracy', by_odor_accuracy, ...
    'by_odor_accuracy_null_centered', by_odor_accuracy_centered, ...
    'by_context_evidence', by_context_evidence, ...
    'by_context_evidence_null_centered', by_context_evidence_centered, ...
    'by_context_recall', by_context_recall, ...
    'by_context_recall_null_centered', by_context_recall_centered, ...
    'confusion_matrices', confusion, ...
    'confusion_matrices_normalized', confusion_normalized);
group_results.inference = struct( ...
    'overall_evidence', overall_evidence_inference, ...
    'overall_accuracy', overall_accuracy_inference, ...
    'context_evidence', context_evidence_inference, ...
    'context_recall', context_recall_inference);
group_results.overall_table = overall_table;
group_results.context_table = context_table;
group_results.odor_table = odor_table;
group_results.subject_table = subject_table;
group_results.subject_odor_table = subject_odor_table;
group_results.output_dir = output_dir;
group_results.matlab_version = version;

if opts.SaveOutputs || opts.MakePlots
    if ~isfolder(output_dir); mkdir(output_dir); end
end
if opts.SaveOutputs
    save(fullfile(output_dir, 'group_cross_odor_context_template_results.mat'), ...
        'group_results', '-v7.3');
    writetable(overall_table, fullfile(output_dir, 'group_overall_cross_odor_context.tsv'), ...
        'FileType', 'text', 'Delimiter', '\t');
    writetable(context_table, fullfile(output_dir, 'group_cross_odor_by_context.tsv'), ...
        'FileType', 'text', 'Delimiter', '\t');
    writetable(odor_table, fullfile(output_dir, 'group_cross_odor_by_target_odor.tsv'), ...
        'FileType', 'text', 'Delimiter', '\t');
    writetable(subject_table, fullfile(output_dir, 'subject_level_cross_odor_summary.tsv'), ...
        'FileType', 'text', 'Delimiter', '\t');
    writetable(subject_odor_table, fullfile(output_dir, 'subject_target_odor_summary.tsv'), ...
        'FileType', 'text', 'Delimiter', '\t');
end
if opts.MakePlots
    create_group_plots(output_dir, roi_names, contexts, odor_values, ...
        overall_evidence, overall_evidence_null, overall_accuracy, overall_accuracy_null, ...
        by_odor_evidence, by_odor_accuracy, by_context_evidence, by_context_recall, ...
        confusion_normalized, reference.chance_accuracy);
end
fprintf(['Cross-odor group inference complete: %d subjects, %d sign flips; ' ...
    'minimum one-sided p = %.5f.\n'], n_subjects, 2^n_subjects, 1/(2^n_subjects));
end

function assert_result_selection(result, expected, subject_id)
if isfield(result, 'analysis') && isfield(result.analysis, 'roi_selection')
    actual = OX_normalize_decoding_roi_selection(result.analysis.roi_selection);
else
    actual = 'old';
end
assert(strcmp(actual, expected), ...
    'Subject %d result uses ROISelection=%s; expected %s.', ...
    subject_id, actual, expected);
end

function root = resolve_mri_root(requested)
if strlength(string(requested)) > 0
    root = char(string(requested)); assert(isfolder(root), 'MRIRoot does not exist: %s', root); return;
end
candidates = {'/Users/qhyang/Desktop/OX_DATA/MRI', '/Volumes/ExtremeSSD/OX_DATA/MRI'};
root = '';
for idx = 1:numel(candidates)
    if isfolder(candidates{idx}); root = candidates{idx}; break; end
end
assert(~isempty(root), 'Could not auto-detect MRIRoot.');
end

function inference = exact_sign_flip(effects)
[n_subjects, n_tests] = size(effects);
assert(all(isfinite(effects), 'all'), 'Sign-flip effects must be finite.');
n_signs = 2^n_subjects;
sign_matrix = ones(n_signs, n_subjects);
for subject_idx = 1:n_subjects
    sign_matrix(:, subject_idx) = 2*bitget((0:n_signs-1)', subject_idx)-1;
end
null_statistics = (sign_matrix*effects)/n_subjects;
observed = mean(effects, 1);
p_uncorrected = sum(null_statistics >= observed, 1)/n_signs;
max_null = max(null_statistics, [], 2);
p_fwe = sum(max_null >= observed, 1)/n_signs;
q_fdr = bh_fdr(p_uncorrected(:))';
inference = struct('observed_null_centered_mean', observed, ...
    'p_uncorrected', p_uncorrected, 'q_fdr', q_fdr, ...
    'p_fwe_maxstat', p_fwe, 'significant_fwe05', p_fwe <= .05, ...
    'null_statistics', null_statistics, 'max_null', max_null, ...
    'n_exact_sign_flips', n_signs, 'alternative', 'greater');
assert(numel(observed) == n_tests, 'Sign-flip output size mismatch.');
end

function run_inference_self_tests()
x = exact_sign_flip(ones(5, 2));
assert(all(x.p_uncorrected == 1/32) && all(x.p_fwe_maxstat >= x.p_uncorrected), ...
    'Exact sign-flip self-test failed.');
assert(max(abs(bh_fdr([.01;.04;.03])-[.03;.04;.04])) < 1e-12, ...
    'BH-FDR self-test failed.');
end

function q = bh_fdr(p_values)
p_values = p_values(:); [sorted, order] = sort(p_values); n = numel(sorted);
adjusted = sorted.*n./(1:n)'; adjusted = flipud(cummin(flipud(adjusted)));
q = nan(n,1); q(order) = min(adjusted,1);
end

function output = make_overall_table(metadata, n_subjects, evidence, evidence_null, ...
        evidence_inference, accuracy, accuracy_null, accuracy_inference)
output = metadata; output.n_subjects = repmat(n_subjects, height(metadata), 1);
output.mean_evidence = mean(evidence,1)';
output.mean_evidence_null = mean(evidence_null,1)';
output.mean_evidence_null_centered = evidence_inference.observed_null_centered_mean';
output.p_evidence = evidence_inference.p_uncorrected';
output.q_evidence_fdr = evidence_inference.q_fdr';
output.p_evidence_fwe = evidence_inference.p_fwe_maxstat';
output.sig_evidence_fwe05 = evidence_inference.significant_fwe05';
output.mean_accuracy = mean(accuracy,1)';
output.mean_accuracy_null = mean(accuracy_null,1)';
output.mean_accuracy_null_centered = accuracy_inference.observed_null_centered_mean';
output.p_accuracy = accuracy_inference.p_uncorrected';
output.q_accuracy_fdr = accuracy_inference.q_fdr';
output.p_accuracy_fwe = accuracy_inference.p_fwe_maxstat';
output.sig_accuracy_fwe05 = accuracy_inference.significant_fwe05';
end

function output = make_context_table(metadata, contexts, n_subjects, evidence, evidence_null, ...
        evidence_inference, recall, recall_null, recall_inference)
n_rois = height(metadata);
output = metadata(repmat((1:n_rois)',3,1),:);
output.context = repelem(contexts(:),n_rois);
output.n_subjects = repmat(n_subjects,height(output),1);
output.mean_evidence = reshape(squeeze(mean(evidence,1)),[],1);
output.mean_evidence_null = reshape(squeeze(mean(evidence_null,1)),[],1);
output.mean_evidence_null_centered = evidence_inference.observed_null_centered_mean';
output.p_evidence = evidence_inference.p_uncorrected';
output.q_evidence_fdr = evidence_inference.q_fdr';
output.p_evidence_fwe = evidence_inference.p_fwe_maxstat';
output.sig_evidence_fwe05 = evidence_inference.significant_fwe05';
output.mean_recall = reshape(squeeze(mean(recall,1)),[],1);
output.mean_recall_null = reshape(squeeze(mean(recall_null,1)),[],1);
output.mean_recall_null_centered = recall_inference.observed_null_centered_mean';
output.p_recall = recall_inference.p_uncorrected';
output.q_recall_fdr = recall_inference.q_fdr';
output.p_recall_fwe = recall_inference.p_fwe_maxstat';
output.sig_recall_fwe05 = recall_inference.significant_fwe05';
end

function output = make_odor_table(metadata, odors, n_subjects, evidence, evidence_null, accuracy, accuracy_null)
n_rois = height(metadata); n_odors = numel(odors);
output = metadata(repmat((1:n_rois)',n_odors,1),:);
output.target_odor = repelem(odors(:),n_rois);
output.n_subjects = repmat(n_subjects,height(output),1);
output.mean_evidence = reshape(squeeze(mean(evidence,1)),[],1);
output.mean_evidence_null = reshape(squeeze(mean(evidence_null,1)),[],1);
output.mean_evidence_null_centered = output.mean_evidence-output.mean_evidence_null;
output.mean_accuracy = reshape(squeeze(mean(accuracy,1)),[],1);
output.mean_accuracy_null = reshape(squeeze(mean(accuracy_null,1)),[],1);
output.mean_accuracy_null_centered = output.mean_accuracy-output.mean_accuracy_null;
end

function output = make_subject_table(subject_ids, roi_names, evidence, evidence_null, ...
        accuracy, accuracy_null, context_evidence, context_recall, contexts)
n_subjects = numel(subject_ids); n_rois = numel(roi_names);
subject_id = repelem(subject_ids(:),n_rois); roi_name = repmat(roi_names(:),n_subjects,1);
output = table(subject_id,roi_name,reshape(evidence',[],1),reshape(evidence_null',[],1), ...
    reshape((evidence-evidence_null)',[],1),reshape(accuracy',[],1), ...
    reshape(accuracy_null',[],1),reshape((accuracy-accuracy_null)',[],1), ...
    'VariableNames',{'subject_id','roi_name','evidence','evidence_null', ...
    'evidence_null_centered','accuracy','accuracy_null','accuracy_null_centered'});
for context_idx=1:3
    label=lower(char(contexts(context_idx)));
    output.([label '_evidence'])=reshape(squeeze(context_evidence(:,:,context_idx))',[],1);
    output.([label '_recall'])=reshape(squeeze(context_recall(:,:,context_idx))',[],1);
end
end

function output = make_subject_odor_table(subject_ids, roi_names, odors, evidence, evidence_null, accuracy, accuracy_null)
n_subjects=numel(subject_ids); n_rois=numel(roi_names); n_odors=numel(odors);
subject_id=repelem(subject_ids(:),n_rois*n_odors);
roi_name=repmat(repelem(roi_names(:),n_odors),n_subjects,1);
target_odor=repmat(repmat(odors(:),n_rois,1),n_subjects,1);
% Arrange each subject as ROI-major, then odor within ROI.
evidence_long=reshape(permute(evidence,[3 2 1]),[],1);
evidence_null_long=reshape(permute(evidence_null,[3 2 1]),[],1);
accuracy_long=reshape(permute(accuracy,[3 2 1]),[],1);
accuracy_null_long=reshape(permute(accuracy_null,[3 2 1]),[],1);
output=table(subject_id,roi_name,target_odor,evidence_long,evidence_null_long, ...
    evidence_long-evidence_null_long,accuracy_long,accuracy_null_long, ...
    accuracy_long-accuracy_null_long,'VariableNames',{'subject_id','roi_name','target_odor', ...
    'evidence','evidence_null','evidence_null_centered','accuracy','accuracy_null', ...
    'accuracy_null_centered'});
end

function create_group_plots(output_dir, roi_names, contexts, odors, evidence, evidence_null, ...
        accuracy, accuracy_null, odor_evidence, odor_accuracy, context_evidence, context_recall, ...
        confusion, chance)
plot_paired(evidence_null,evidence,roi_names,'Evidence','Cross-odor context evidence', ...
    fullfile(output_dir,'overall_cross_odor_evidence'),NaN);
plot_paired(accuracy_null,accuracy,roi_names,'Accuracy','Cross-odor context accuracy', ...
    fullfile(output_dir,'overall_cross_odor_accuracy'),chance);
plot_heatmap(squeeze(mean(odor_evidence,1)),roi_names,string(odors), ...
    'Target odor','Cross-odor context evidence by target odor', ...
    fullfile(output_dir,'target_odor_evidence'));
plot_heatmap(squeeze(mean(odor_accuracy,1)),roi_names,string(odors), ...
    'Target odor','Cross-odor context accuracy by target odor', ...
    fullfile(output_dir,'target_odor_accuracy'));
plot_context(context_evidence,roi_names,contexts,'Evidence', ...
    'Cross-odor evidence by true context',fullfile(output_dir,'context_evidence'),NaN);
plot_context(context_recall,roi_names,contexts,'Recall', ...
    'Cross-odor recall by true context',fullfile(output_dir,'context_recall'),chance);
plot_confusions(squeeze(mean(confusion,1)),roi_names,contexts, ...
    fullfile(output_dir,'context_confusion_matrices'));
end

function plot_paired(null_values,observed,roi_names,y_label,figure_title,stem,reference)
n_rois=numel(roi_names); fig=figure('Visible','off','Color','w','Position',[50 50 1500 1100]);
layout=tiledlayout(5,5,'TileSpacing','compact','Padding','compact');
for r=1:n_rois
    nexttile; hold on; plot([1 2],[null_values(:,r),observed(:,r)]','-','Color',[.75 .75 .75]);
    scatter(ones(size(observed,1),1),null_values(:,r),16,'filled');
    scatter(2*ones(size(observed,1),1),observed(:,r),16,'filled');
    plot([1 2],[mean(null_values(:,r)),mean(observed(:,r))],'-ok','LineWidth',1.2,'MarkerFaceColor','k');
    if isfinite(reference); yline(reference,'--k'); end
    xlim([.7 2.3]); xticks([1 2]); xticklabels({'Null','Observed'});
    title(roi_names(r),'Interpreter','none','FontSize',8);
end
ylabel(layout,y_label); title(layout,figure_title); export_plot(fig,stem);
end

function plot_heatmap(values,roi_names,column_labels,x_label,figure_title,stem)
fig=figure('Visible','off','Color','w','Position',[50 50 1500 850]);
imagesc(values); colorbar; yticks(1:numel(roi_names)); yticklabels(roi_names);
xticks(1:numel(column_labels)); xticklabels(column_labels); xtickangle(45);
xlabel(x_label); ylabel('ROI'); title(figure_title); set(gca,'TickLabelInterpreter','none');
export_plot(fig,stem);
end

function plot_context(values,roi_names,contexts,y_label,figure_title,stem,reference)
n_rois=numel(roi_names); fig=figure('Visible','off','Color','w','Position',[50 50 1500 1100]);
layout=tiledlayout(5,5,'TileSpacing','compact','Padding','compact');
for r=1:n_rois
    nexttile; hold on; x=squeeze(values(:,r,:)); plot(1:3,x','-','Color',[.75 .75 .75]);
    scatter(repmat(1:3,size(x,1),1),x,15,'filled');
    plot(1:3,mean(x,1),'-ok','LineWidth',1.2,'MarkerFaceColor','k');
    if isfinite(reference); yline(reference,'--k'); end
    xticks(1:3); xticklabels(contexts); xtickangle(35); title(roi_names(r),'Interpreter','none','FontSize',8);
end
xlabel(layout,'True context'); ylabel(layout,y_label); title(layout,figure_title); export_plot(fig,stem);
end

function plot_confusions(values,roi_names,contexts,stem)
n_rois=numel(roi_names); fig=figure('Visible','off','Color','w','Position',[50 50 1500 1100]);
layout=tiledlayout(5,5,'TileSpacing','compact','Padding','compact');
for r=1:n_rois
    nexttile; imagesc(squeeze(values(r,:,:)),[0 1]); axis image;
    xticks(1:3); yticks(1:3); xticklabels(contexts); yticklabels(contexts); xtickangle(35);
    title(roi_names(r),'Interpreter','none','FontSize',8);
end
xlabel(layout,'Predicted context'); ylabel(layout,'True context');
title(layout,'Normalized cross-odor context confusion matrices'); colorbar; export_plot(fig,stem);
end

function export_plot(fig,stem)
exportgraphics(fig,[stem '.png'],'Resolution',200);
exportgraphics(fig,[stem '.pdf'],'ContentType','vector'); close(fig);
end
