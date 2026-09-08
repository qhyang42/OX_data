function outputs = run_exploratory_rsa_first_batch(varargin)
%RUN_EXPLORATORY_RSA_FIRST_BATCH Broad descriptive ROI RSA analyses.
%   OUTPUTS = RUN_EXPLORATORY_RSA_FIRST_BATCH() runs the six exploratory
%   analyses requested for subjects 2--6 and the five core olfactory ROIs.
%   Results are written below RDMs/exploratory_RSA_results.

parser = inputParser;
parser.addParameter('MakePlots', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('OutputDir', '', @(x) ischar(x) || isstring(x));
parser.parse(varargin{:});
make_plots = parser.Results.MakePlots;

script_path = mfilename('fullpath');
assert(~isempty(script_path), 'Run this function from its saved file.');
project_root = fileparts(fileparts(script_path));
rdm_root = fullfile(project_root, 'RDMs');
output_root = char(parser.Results.OutputDir);
if isempty(output_root)
    output_root = fullfile(rdm_root, 'exploratory_RSA_results');
end
table_root = fullfile(output_root, 'tables');
figure_root = fullfile(output_root, 'figures');
ensure_directory(output_root);
ensure_directory(table_root);
ensure_directory(figure_root);
subject_root = fullfile(output_root, 'subject_results');
ensure_directory(subject_root);
for analysis_index = 1:6
    ensure_directory(fullfile(figure_root, sprintf('analysis%d', analysis_index)));
end

subject_ids = 2:6;
roi_names = ["AON", "PirF", "PirT", "olfAMG", "olfOFC"];
metric_names = ["simple", "crossnobis"];
version_names = ["linear", "rank"];
context_names = ["PERSON", "FOOD", "LOCATION", "CONTROL"];
predictor_names = ["odor", "context", "pleasantness", "intensity"];
perceptual_names = ["pleasantness", "intensity"];
context_pairs = nchoosek(1:numel(context_names), 2);
upper80 = triu(true(80), 1);
upper20 = triu(true(20), 1);

categorical_file = fullfile(rdm_root, 'categorical_RDMs.mat');
assert(isfile(categorical_file), 'Missing categorical RDM file: %s', categorical_file);
categorical = load(categorical_file);
assert_common_order(categorical.condition_metadata, categorical.context_order, ...
    categorical.odor_ids, context_names);
validate_input_rdm(categorical.D_odor, 1, 'categorical odor identity');
validate_input_rdm(categorical.D_context, 1, 'categorical context identity');
assert(all(ismember(categorical.D_odor(:), [0 1])) && ...
    all(ismember(categorical.D_context(:), [0 1])), ...
    'Categorical RDMs must be binary same-identity indicators.');
condition_order = categorical.condition_metadata;
writetable(condition_order, fullfile(output_root, 'condition_ordering.csv'));
upper_indices = find(upper80);
[upper_i, upper_j] = ind2sub([80, 80], upper_indices);
upper_triangle_index = table((1:numel(upper_indices))', upper_i, upper_j, ...
    string(condition_order.condition_context(upper_i)), ...
    string(condition_order.condition_context(upper_j)), ...
    condition_order.condition_odor(upper_i), condition_order.condition_odor(upper_j), ...
    'VariableNames', {'observation','condition_i','condition_j','context_i', ...
    'context_j','odor_i','odor_j'});
writetable(upper_triangle_index, fullfile(output_root, 'rdm_upper_triangle_index.csv'));
context_pair_table = table((1:6)', context_names(context_pairs(:,1))', ...
    context_names(context_pairs(:,2))', ...
    (context_names(context_pairs(:,1)) + "_vs_" + context_names(context_pairs(:,2)))', ...
    'VariableNames', {'pair_index','context_1','context_2','context_pair'});
writetable(context_pair_table, fullfile(output_root, 'context_pairs.csv'));
save(fullfile(output_root, 'condition_ordering.mat'), 'condition_order', ...
    'context_names', 'subject_ids', 'roi_names', 'metric_names');

dummy_fit4 = empty_fit(4, 0);
dummy_fit2 = empty_fit(2, 0);
omnibus_rows = repmat(omnibus_row(0, "", "", "", dummy_fit4), 0, 1);
predictor_corr_rows = repmat(predictor_corr_row(0, "", "", "", "", "", NaN, 0), 0, 1);
agreement_rows = repmat(agreement_row(0, "", NaN, NaN, 0, NaN, NaN, NaN, NaN, NaN), 0, 1);
rdm_point_rows = repmat(rdm_point_row(0, "", 0, 0, 0, NaN, NaN), 0, 1);
geometry_rows = repmat(geometry_row(0, "", "", "", "", "", NaN, NaN, 0), 0, 1);
behavior_geometry_rows = repmat(geometry_row(0, "", "", "", "", "", NaN, NaN, 0), 0, 1);
context_rsa_rows = repmat(context_rsa_row(0, "", "", "", "", dummy_fit2), 0, 1);
warping_rows = repmat(warping_row(0, "", "", "", "", dummy_fit2), 0, 1);
displacement_rows = repmat(displacement_row(0, "", "", "", dummy_fit2), 0, 1);
displacement_coef_rows = repmat(struct('subject_id', 0, 'subject', "", 'roi', "", ...
    'neural_metric', "", 'model', "", 'term', "", 'coefficient', NaN), 0, 1);
displacement_observation_rows = repmat(struct('subject_id', 0, 'subject', "", ...
    'roi', "", 'neural_metric', "", 'observation', 0, 'odor_id', 0, ...
    'context_pair', "", 'neural_distance', NaN, ...
    'abs_delta_pleasantness', NaN, 'abs_delta_intensity', NaN), 0, 1);
input_qc_rows = repmat(struct('subject_id', 0, 'subject', "", ...
    'missing_pleasantness_trials', 0, 'missing_intensity_trials', 0, ...
    'min_trials_per_condition', 0, 'max_trials_per_condition', 0, ...
    'n_conditions_not_equal_10', 0, 'status', ""), 0, 1);
qc_rows = repmat(qc_row(0, "", "", "", dummy_fit2), 0, 1);

fprintf('[EXPLORATORY RSA] Output: %s\n', output_root);
for subject_id = subject_ids
    subject_name = "subj_" + string(subject_id);
    fprintf('[EXPLORATORY RSA] Loading %s\n', subject_name);
    behavior_file = fullfile(rdm_root, subject_name + "_behavioral_RDMs.mat");
    neural_file = fullfile(rdm_root, 'neural', subject_name + "_neural_RDMs.mat");
    assert(isfile(behavior_file), 'Missing behavioral RDM file: %s', behavior_file);
    assert(isfile(neural_file), 'Missing neural RDM file: %s', neural_file);
    behavior = load(behavior_file);
    neural_loaded = load(neural_file, 'results');
    neural = neural_loaded.results;
    assert_common_order(behavior.condition_metadata, behavior.context_order, ...
        behavior.odor_ids, context_names);
    validate_input_rdm(behavior.D_pleasantness, 0, subject_name + " pleasantness");
    validate_input_rdm(behavior.D_intensity, 0, subject_name + " intensity");
    assert(isequal(condition_order.condition_index, behavior.condition_metadata.condition_index) && ...
        isequal(string(condition_order.condition_context), string(behavior.condition_metadata.condition_context)) && ...
        isequal(condition_order.condition_odor, behavior.condition_metadata.condition_odor), ...
        'Behavioral condition order mismatch for %s.', subject_name);
    assert(isequal(condition_order.condition_index, neural.condition_metadata.condition_index) && ...
        isequal(string(condition_order.condition_context), string(neural.condition_metadata.condition_context)) && ...
        isequal(condition_order.condition_odor, neural.condition_metadata.condition_odor), ...
        'Neural condition order mismatch for %s.', subject_name);
    missing_p = sum(behavior.total_trials_by_condition - behavior.valid_pleasantness_trials);
    missing_i = sum(behavior.total_trials_by_condition - behavior.valid_intensity_trials);
    count_deviation = sum(behavior.total_trials_by_condition ~= 10);
    input_status = "ok";
    warnings = strings(0,1);
    if missing_p > 0, warnings(end+1) = "missing_pleasantness_trials"; end %#ok<AGROW>
    if missing_i > 0, warnings(end+1) = "missing_intensity_trials"; end %#ok<AGROW>
    if count_deviation > 0, warnings(end+1) = "condition_count_imbalance"; end %#ok<AGROW>
    if ~isempty(warnings), input_status = strjoin(warnings, '|'); end
    input_qc_rows(end+1) = struct('subject_id', subject_id, 'subject', subject_name, ... %#ok<AGROW>
        'missing_pleasantness_trials', missing_p, 'missing_intensity_trials', missing_i, ...
        'min_trials_per_condition', min(behavior.total_trials_by_condition), ...
        'max_trials_per_condition', max(behavior.total_trials_by_condition), ...
        'n_conditions_not_equal_10', count_deviation, 'status', input_status);

    behavior_vectors = [behavior.D_pleasantness(upper80), behavior.D_intensity(upper80)];
    categorical_vectors = [categorical.D_odor(upper80), categorical.D_context(upper80)];
    omnibus_X = [categorical_vectors, behavior_vectors];

    % Behavioral context-geometry is independent of ROI and neural metric.
    behavior_rdms = {behavior.D_pleasantness, behavior.D_intensity};
    for measure_index = 1:2
        blocks = extract_context_blocks(behavior_rdms{measure_index});
        for pair_index = 1:size(context_pairs, 1)
            c1 = context_pairs(pair_index, 1);
            c2 = context_pairs(pair_index, 2);
            x = blocks{c1}(upper20);
            y = blocks{c2}(upper20);
            [r_pearson, n_valid_p] = pair_correlation(x, y, false);
            [r_spearman, n_valid_s] = pair_correlation(x, y, true);
            behavior_geometry_rows(end+1) = geometry_row(subject_id, "", "", ... %#ok<AGROW>
                perceptual_names(measure_index), context_names(c1), context_names(c2), ...
                r_pearson, r_spearman, min(n_valid_p, n_valid_s));
            qc_rows(end+1) = correlation_qc_row(subject_id, "", ... %#ok<AGROW>
                perceptual_names(measure_index), "3_behavior_geometry_" + ...
                context_names(c1) + "_vs_" + context_names(c2), ...
                numel(x), min(n_valid_p,n_valid_s), r_pearson, r_spearman);
        end
    end

    for roi_index = 1:numel(roi_names)
        roi_name = roi_names(roi_index);
        saved_roi_names = string({neural.roi_results.roi_name});
        match = find(strcmpi(saved_roi_names, roi_name));
        assert(isscalar(match), 'Expected one %s ROI in %s.', roi_name, subject_name);
        roi = neural.roi_results(match);
        D_by_metric = {double(roi.D_neural_simple), double(roi.D_neural_crossnobis)};
        for metric_index = 1:2
            metric_name = metric_names(metric_index);
            D_neural = D_by_metric{metric_index};
            validate_rdm(D_neural, subject_name, roi_name, metric_name);
            y80 = D_neural(upper80);

            % Analysis 1: omnibus multiple-regression RSA.
            for version_index = 1:2
                version_name = version_names(version_index);
                fit = fit_standardized_ols(y80, omnibus_X, predictor_names, version_name == "rank");
                omnibus_rows(end+1) = omnibus_row(subject_id, roi_name, metric_name, ... %#ok<AGROW>
                    version_name, fit);
                qc_rows(end+1) = qc_row(subject_id, roi_name, metric_name, ... %#ok<AGROW>
                    "1_omnibus_" + version_name, fit);
                for p1 = 1:4
                    for p2 = 1:4
                        predictor_corr_rows(end+1) = predictor_corr_row(subject_id, ... %#ok<AGROW>
                            roi_name, metric_name, version_name, predictor_names(p1), ...
                            predictor_names(p2), fit.predictor_correlation(p1, p2), fit.n_valid);
                    end
                end
            end

            % Analysis 3: context-specific odor geometry.
            neural_blocks = extract_context_blocks(D_neural);
            for pair_index = 1:size(context_pairs, 1)
                c1 = context_pairs(pair_index, 1);
                c2 = context_pairs(pair_index, 2);
                x = neural_blocks{c1}(upper20);
                y = neural_blocks{c2}(upper20);
                [r_pearson, n_valid_p] = pair_correlation(x, y, false);
                [r_spearman, n_valid_s] = pair_correlation(x, y, true);
                geometry_rows(end+1) = geometry_row(subject_id, roi_name, metric_name, ... %#ok<AGROW>
                    "neural", context_names(c1), context_names(c2), ...
                    r_pearson, r_spearman, min(n_valid_p, n_valid_s));
                qc_rows(end+1) = correlation_qc_row(subject_id, roi_name, metric_name, ... %#ok<AGROW>
                    "3_neural_geometry_" + context_names(c1) + "_vs_" + context_names(c2), ...
                    numel(x), min(n_valid_p,n_valid_s), r_pearson, r_spearman);
            end

            % Analysis 4: context-specific perceptual RSA.
            pleasant_blocks = extract_context_blocks(behavior.D_pleasantness);
            intensity_blocks = extract_context_blocks(behavior.D_intensity);
            for context_index = 1:4
                context_X = [pleasant_blocks{context_index}(upper20), ...
                    intensity_blocks{context_index}(upper20)];
                context_y = neural_blocks{context_index}(upper20);
                for version_index = 1:2
                    version_name = version_names(version_index);
                    fit = fit_standardized_ols(context_y, context_X, ...
                        perceptual_names, version_name == "rank");
                    context_rsa_rows(end+1) = context_rsa_row(subject_id, ... %#ok<AGROW>
                        roi_name, metric_name, context_names(context_index), version_name, fit);
                    qc_rows(end+1) = qc_row(subject_id, roi_name, metric_name, ... %#ok<AGROW>
                        "4_context_RSA_" + context_names(context_index) + "_" + version_name, fit);
                end
            end

            % Analysis 5: context-dependent geometry warping.
            for pair_index = 1:size(context_pairs, 1)
                c1 = context_pairs(pair_index, 1);
                c2 = context_pairs(pair_index, 2);
                pair_name = context_names(c1) + "_vs_" + context_names(c2);
                delta_y = neural_blocks{c1}(upper20) - neural_blocks{c2}(upper20);
                delta_X = [pleasant_blocks{c1}(upper20) - pleasant_blocks{c2}(upper20), ...
                    intensity_blocks{c1}(upper20) - intensity_blocks{c2}(upper20)];
                for version_index = 1:2
                    version_name = version_names(version_index);
                    fit = fit_standardized_ols(delta_y, delta_X, perceptual_names, ...
                        version_name == "rank");
                    warping_rows(end+1) = warping_row(subject_id, roi_name, ... %#ok<AGROW>
                        metric_name, pair_name, version_name, fit);
                    qc_rows(end+1) = qc_row(subject_id, roi_name, metric_name, ... %#ok<AGROW>
                        "5_warping_" + pair_name + "_" + version_name, fit);
                end
            end

            % Analysis 6: same-odor, cross-context displacement.
            [disp_y, disp_X, disp_pair, disp_odor] = displacement_data( ...
                D_neural, behavior.D_pleasantness, behavior.D_intensity, context_names);
            [fit_basic, coef_basic] = fit_factor_ols(disp_y, disp_X, ...
                perceptual_names, disp_pair, strings(size(disp_pair)), false);
            [fit_odor, coef_odor] = fit_factor_ols(disp_y, disp_X, ...
                perceptual_names, disp_pair, "odor" + string(disp_odor), true);
            displacement_rows(end+1) = displacement_row(subject_id, roi_name, ... %#ok<AGROW>
                metric_name, "context_pair", fit_basic);
            displacement_rows(end+1) = displacement_row(subject_id, roi_name, ... %#ok<AGROW>
                metric_name, "context_pair_plus_odor", fit_odor);
            displacement_coef_rows = append_coefficient_rows(displacement_coef_rows, ...
                subject_id, roi_name, metric_name, "context_pair", coef_basic);
            displacement_coef_rows = append_coefficient_rows(displacement_coef_rows, ...
                subject_id, roi_name, metric_name, "context_pair_plus_odor", coef_odor);
            qc_rows(end+1) = qc_row(subject_id, roi_name, metric_name, ... %#ok<AGROW>
                "6_displacement_context_pair", fit_basic);
            qc_rows(end+1) = qc_row(subject_id, roi_name, metric_name, ... %#ok<AGROW>
                "6_displacement_context_pair_plus_odor", fit_odor);
            for observation_index = 1:numel(disp_y)
                displacement_observation_rows(end+1) = struct( ... %#ok<AGROW>
                    'subject_id', subject_id, 'subject', subject_name, ...
                    'roi', roi_name, 'neural_metric', metric_name, ...
                    'observation', observation_index, 'odor_id', disp_odor(observation_index), ...
                    'context_pair', disp_pair(observation_index), ...
                    'neural_distance', disp_y(observation_index), ...
                    'abs_delta_pleasantness', disp_X(observation_index,1), ...
                    'abs_delta_intensity', disp_X(observation_index,2));
            end
        end

        % Analysis 2: neural-distance and omnibus-beta concordance.
        simple_vector = D_by_metric{1}(upper80);
        crossnobis_vector = D_by_metric{2}(upper80);
        [pearson_r, n_valid_p] = pair_correlation(simple_vector, crossnobis_vector, false);
        [spearman_r, n_valid_s] = pair_correlation(simple_vector, crossnobis_vector, true);
        valid = isfinite(simple_vector) & isfinite(crossnobis_vector);
        cross_valid = crossnobis_vector(valid);
        agreement_rows(end+1) = agreement_row(subject_id, roi_name, pearson_r, ... %#ok<AGROW>
            spearman_r, min(n_valid_p, n_valid_s), mean(cross_valid < 0), ...
            mean(cross_valid), median(cross_valid), mean(simple_vector(valid)), ...
            median(simple_vector(valid)));
        qc_rows(end+1) = correlation_qc_row(subject_id, roi_name, "simple_vs_crossnobis", ... %#ok<AGROW>
            "2_metric_concordance", numel(simple_vector), min(n_valid_p,n_valid_s), ...
            pearson_r, spearman_r);
        indices = find(upper80);
        [condition_i, condition_j] = ind2sub([80, 80], indices);
        for observation_index = 1:numel(indices)
            rdm_point_rows(end+1) = rdm_point_row(subject_id, roi_name, ... %#ok<AGROW>
                observation_index, condition_i(observation_index), condition_j(observation_index), ...
                simple_vector(observation_index), crossnobis_vector(observation_index));
        end
    end
end

% Convert, augment, and save all detailed tables.
omnibus_table = struct2table(omnibus_rows);
predictor_correlation_table = struct2table(predictor_corr_rows);
rdm_agreement_table = struct2table(agreement_rows);
rdm_points_table = struct2table(rdm_point_rows);
context_geometry_table = struct2table(geometry_rows);
behavior_geometry_table = struct2table(behavior_geometry_rows);
context_geometry_matrices_table = build_geometry_matrix_table( ...
    context_geometry_table, behavior_geometry_table, context_names);
context_perceptual_rsa_table = struct2table(context_rsa_rows);
warping_table = struct2table(warping_rows);
warping_summary_table = summarize_warping(warping_table);
displacement_table = struct2table(displacement_rows);
displacement_coefficients_table = struct2table(displacement_coef_rows);
displacement_observations_table = struct2table(displacement_observation_rows);
input_qc_table = struct2table(input_qc_rows);
qc_table = struct2table(qc_rows);
summary_table = build_summary_table(omnibus_table, rdm_agreement_table, ...
    context_geometry_table, context_perceptual_rsa_table, warping_table, displacement_table);

write_table(omnibus_table, table_root, 'analysis1_omnibus_RSA');
write_table(predictor_correlation_table, table_root, 'analysis1_predictor_correlations');
write_table(rdm_agreement_table, table_root, 'analysis2_RDM_agreement');
write_table(rdm_points_table, table_root, 'analysis2_RDM_points');
write_table(context_geometry_table, table_root, 'analysis3_neural_context_geometry');
write_table(behavior_geometry_table, table_root, 'analysis3_behavior_context_geometry');
write_table(context_geometry_matrices_table, table_root, 'analysis3_context_geometry_matrices');
write_table(context_perceptual_rsa_table, table_root, 'analysis4_context_perceptual_RSA');
write_table(warping_table, table_root, 'analysis5_geometry_warping');
write_table(warping_summary_table, table_root, 'analysis5_geometry_warping_summary');
write_table(displacement_table, table_root, 'analysis6_same_odor_displacement');
write_table(displacement_coefficients_table, table_root, 'analysis6_all_model_coefficients');
write_table(displacement_observations_table, table_root, 'analysis6_same_odor_observations');
write_table(input_qc_table, table_root, 'QC_input_status');
write_table(qc_table, table_root, 'QC_model_status');
write_table(summary_table, table_root, 'all_analyses_summary');

for subject_id = subject_ids
    subject_name = "subj_" + string(subject_id);
    subject_results = struct();
    subject_results.omnibus = omnibus_table(omnibus_table.subject_id == subject_id, :);
    subject_results.predictor_correlations = predictor_correlation_table(predictor_correlation_table.subject_id == subject_id, :);
    subject_results.rdm_agreement = rdm_agreement_table(rdm_agreement_table.subject_id == subject_id, :);
    subject_results.rdm_points = rdm_points_table(rdm_points_table.subject_id == subject_id, :);
    subject_results.context_geometry = context_geometry_table(context_geometry_table.subject_id == subject_id, :);
    subject_results.behavior_geometry = behavior_geometry_table(behavior_geometry_table.subject_id == subject_id, :);
    subject_results.context_geometry_matrices = context_geometry_matrices_table( ...
        context_geometry_matrices_table.subject_id == subject_id | ...
        context_geometry_matrices_table.subject_id == 0, :);
    subject_results.context_perceptual_rsa = context_perceptual_rsa_table(context_perceptual_rsa_table.subject_id == subject_id, :);
    subject_results.warping = warping_table(warping_table.subject_id == subject_id, :);
    subject_results.warping_summary = warping_summary_table(warping_summary_table.subject_id == subject_id, :);
    subject_results.displacement = displacement_table(displacement_table.subject_id == subject_id, :);
    subject_results.displacement_coefficients = displacement_coefficients_table(displacement_coefficients_table.subject_id == subject_id, :);
    subject_results.displacement_observations = displacement_observations_table(displacement_observations_table.subject_id == subject_id, :);
    subject_results.summary = summary_table(summary_table.subject_id == subject_id, :);
    subject_results.model_qc = qc_table(qc_table.subject_id == subject_id, :);
    subject_results.input_qc = input_qc_table(input_qc_table.subject_id == subject_id, :);
    save(fullfile(subject_root, subject_name + "_exploratory_RSA_results.mat"), ...
        'subject_results', '-v7.3');
end

save(fullfile(output_root, 'exploratory_RSA_all_results.mat'), ...
    'omnibus_table', 'predictor_correlation_table', 'rdm_agreement_table', ...
    'rdm_points_table', 'context_geometry_table', 'behavior_geometry_table', ...
    'context_geometry_matrices_table', ...
    'context_perceptual_rsa_table', 'warping_table', 'warping_summary_table', ...
    'displacement_table', 'displacement_coefficients_table', 'qc_table', ...
    'displacement_observations_table', 'input_qc_table', 'summary_table', ...
    'condition_order', 'upper_triangle_index', 'context_pair_table', 'subject_ids', 'roi_names', ...
    'metric_names', 'context_names', '-v7.3');

if make_plots
    fprintf('[EXPLORATORY RSA] Creating figures\n');
    make_all_plots(figure_root, omnibus_table, predictor_correlation_table, ...
        rdm_agreement_table, rdm_points_table, context_geometry_table, ...
        behavior_geometry_table, context_perceptual_rsa_table, ...
        warping_table, displacement_table, displacement_observations_table, ...
        subject_ids, roi_names, ...
        metric_names, context_names);
end

write_readme(output_root, categorical_file, nnz(qc_table.status ~= "ok"), ...
    height(qc_table), nnz(input_qc_table.status ~= "ok"), make_plots);
run_exploratory_rsa_analysis6_pair_slopes('OutputDir', output_root, 'MakePlots', make_plots);
outputs = struct('output_dir', output_root, 'summary_table', summary_table, ...
    'qc_table', qc_table, 'omnibus_table', omnibus_table);
fprintf('[EXPLORATORY RSA] Complete: %d summary rows, %d QC warnings.\n', ...
    height(summary_table), nnz(qc_table.status ~= "ok"));
end

function assert_common_order(metadata, context_order, odor_ids, expected_contexts)
assert(height(metadata) == 80, 'Expected 80 conditions.');
assert(isequal(string(context_order(:)), expected_contexts(:)), 'Unexpected context order.');
assert(isequal(double(odor_ids(:)), (1:20)'), 'Unexpected odor order.');
assert(isequal(double(metadata.condition_index(:)), (1:80)'), 'Unexpected condition indices.');
assert(isequal(string(metadata.condition_context(:)), repelem(expected_contexts(:), 20)), ...
    'Conditions are not context-major.');
assert(isequal(double(metadata.condition_odor(:)), repmat((1:20)', 4, 1)), ...
    'Conditions are not odor-minor within context.');
end

function validate_rdm(D, subject_name, roi_name, metric_name)
assert(isequal(size(D), [80, 80]), '%s %s %s RDM is not 80-by-80.', ...
    subject_name, roi_name, metric_name);
assert(all(isfinite(D), 'all'), '%s %s %s RDM contains nonfinite values.', ...
    subject_name, roi_name, metric_name);
assert(max(abs(D - D'), [], 'all', 'omitnan') < 1e-8, ...
    '%s %s %s RDM is not symmetric.', subject_name, roi_name, metric_name);
assert(all(abs(diag(D)) < 1e-8 | isnan(diag(D))), ...
    '%s %s %s RDM diagonal is nonzero.', subject_name, roi_name, metric_name);
end

function validate_input_rdm(D, expected_diagonal, label)
assert(isequal(size(D), [80, 80]), '%s RDM is not 80-by-80.', label);
assert(isequal(isfinite(D), isfinite(D')), '%s RDM has asymmetric nonfinite entries.', label);
finite_difference = abs(D - D');
finite_difference(~isfinite(finite_difference)) = 0;
assert(max(finite_difference, [], 'all') < 1e-8, '%s RDM is not symmetric.', label);
assert(all(isfinite(D), 'all'), '%s RDM contains nonfinite values.', label);
assert(all(abs(diag(D) - expected_diagonal) < 1e-8), ...
    '%s RDM has an unexpected diagonal.', label);
end

function blocks = extract_context_blocks(D)
blocks = cell(4, 1);
for context_index = 1:4
    indices = (context_index - 1) * 20 + (1:20);
    blocks{context_index} = D(indices, indices);
end
end

function fit = fit_standardized_ols(y, X, predictor_names, rank_transform)
y = double(y(:));
X = double(X);
valid = isfinite(y) & all(isfinite(X), 2);
fit = empty_fit(numel(predictor_names), sum(valid));
fit.n_total = numel(y);
fit.n_valid = sum(valid);
fit.n_dropped = fit.n_total - fit.n_valid;
if fit.n_valid <= size(X, 2) + 1
    fit.status = "too_few_valid_rows";
    return
end
y = y(valid);
X = X(valid, :);
if rank_transform
    y = tied_rank_local(y);
    for column = 1:size(X, 2)
        X(:, column) = tied_rank_local(X(:, column));
    end
end
[y, y_ok] = zscore_local(y);
for column = 1:size(X, 2)
    [X(:, column), x_ok(column)] = zscore_local(X(:, column)); %#ok<AGROW>
end
if ~y_ok
    fit.status = "constant_outcome";
    return
elseif ~all(x_ok)
    fit.status = "constant_predictor:" + strjoin(predictor_names(~x_ok), '|');
    return
end
design = [ones(size(X, 1), 1), X];
fit.design_rank = rank(design);
if fit.design_rank < size(design, 2)
    fit.status = "rank_deficient_design";
    return
end
b = design \ y;
y_hat = design * b;
sse = sum((y - y_hat).^2);
sst = sum((y - mean(y)).^2);
fit.beta = b(2:end)';
fit.r2 = 1 - sse / sst;
fit.adjusted_r2 = 1 - (1 - fit.r2) * (fit.n_valid - 1) / ...
    (fit.n_valid - size(X, 2) - 1);
fit.predictor_correlation = corrcoef(X);
fit.delta_r2 = nan(1, size(X, 2));
fit.partial_r2 = nan(1, size(X, 2));
for predictor_index = 1:size(X, 2)
    reduced_X = X;
    reduced_X(:, predictor_index) = [];
    reduced_design = [ones(size(reduced_X, 1), 1), reduced_X];
    reduced_residual = y - reduced_design * (reduced_design \ y);
    reduced_sse = sum(reduced_residual.^2);
    reduced_r2 = 1 - reduced_sse / sst;
    fit.delta_r2(predictor_index) = max(0, fit.r2 - reduced_r2);
    if reduced_sse > 0
        fit.partial_r2(predictor_index) = max(0, (reduced_sse - sse) / reduced_sse);
    end
end
fit.status = "ok";
end

function fit = empty_fit(n_predictors, n_valid)
fit = struct('beta', nan(1, n_predictors), 'r2', NaN, 'adjusted_r2', NaN, ...
    'delta_r2', nan(1, n_predictors), 'partial_r2', nan(1, n_predictors), ...
    'predictor_correlation', nan(n_predictors), 'n_total', NaN, ...
    'n_valid', n_valid, 'n_dropped', NaN, 'design_rank', NaN, 'status', "not_fit");
end

function [fit, coefficient_table] = fit_factor_ols(y, continuous_X, continuous_names, ...
    pair_labels, odor_labels, include_odor)
y = double(y(:));
continuous_X = double(continuous_X);
pair_labels = string(pair_labels(:));
odor_labels = string(odor_labels(:));
valid = isfinite(y) & all(isfinite(continuous_X), 2) & pair_labels ~= "";
if include_odor
    valid = valid & odor_labels ~= "";
end
fit = empty_fit(numel(continuous_names), sum(valid));
fit.n_total = numel(y);
fit.n_valid = sum(valid);
fit.n_dropped = fit.n_total - fit.n_valid;
y = y(valid);
continuous_X = continuous_X(valid, :);
pair_labels = pair_labels(valid);
odor_labels = odor_labels(valid);
[y, y_ok] = zscore_local(y);
x_ok = false(1, size(continuous_X, 2));
for column = 1:size(continuous_X, 2)
    [continuous_X(:, column), x_ok(column)] = zscore_local(continuous_X(:, column));
end
[pair_dummy, pair_terms] = treatment_dummy(pair_labels, "context_pair");
factor_X = pair_dummy;
factor_terms = pair_terms;
if include_odor
    [odor_dummy, odor_terms] = treatment_dummy(odor_labels, "odor");
    factor_X = [factor_X, odor_dummy];
    factor_terms = [factor_terms, odor_terms];
end
term_names = ["intercept", continuous_names, factor_terms];
coefficient_table = table(term_names(:), nan(numel(term_names), 1), ...
    'VariableNames', {'term', 'coefficient'});
if ~y_ok
    fit.status = "constant_outcome";
    return
elseif ~all(x_ok)
    fit.status = "constant_predictor:" + strjoin(continuous_names(~x_ok), '|');
    return
end
design = [ones(size(y)), continuous_X, factor_X];
fit.design_rank = rank(design);
if fit.n_valid <= size(design, 2) || fit.design_rank < size(design, 2)
    fit.status = "rank_deficient_design";
    return
end
b = design \ y;
coefficient_table.coefficient = b;
fit.beta = b(2:1+numel(continuous_names))';
residual = y - design * b;
sse = sum(residual.^2);
sst = sum((y - mean(y)).^2);
fit.r2 = 1 - sse / sst;
fit.adjusted_r2 = 1 - (1 - fit.r2) * (fit.n_valid - 1) / ...
    (fit.n_valid - size(design, 2));
fit.status = "ok";
end

function [dummy, term_names] = treatment_dummy(labels, prefix)
levels = unique(labels, 'stable');
dummy = zeros(numel(labels), max(0, numel(levels) - 1));
term_names = strings(1, size(dummy, 2));
for level_index = 2:numel(levels)
    dummy(:, level_index - 1) = labels == levels(level_index);
    term_names(level_index - 1) = prefix + "_" + levels(level_index);
end
end

function [z, ok] = zscore_local(x)
x = double(x(:));
sigma = std(x, 0);
ok = isfinite(sigma) && sigma > sqrt(eps);
if ok
    z = (x - mean(x)) ./ sigma;
else
    z = nan(size(x));
end
end

function ranks = tied_rank_local(x)
[sorted_x, order] = sort(x(:));
ranks_sorted = zeros(size(sorted_x));
start_index = 1;
while start_index <= numel(sorted_x)
    end_index = start_index;
    while end_index < numel(sorted_x) && sorted_x(end_index + 1) == sorted_x(start_index)
        end_index = end_index + 1;
    end
    ranks_sorted(start_index:end_index) = mean(start_index:end_index);
    start_index = end_index + 1;
end
ranks = zeros(size(ranks_sorted));
ranks(order) = ranks_sorted;
end

function [r, n_valid] = pair_correlation(x, y, rank_transform)
x = double(x(:));
y = double(y(:));
valid = isfinite(x) & isfinite(y);
n_valid = sum(valid);
if n_valid < 3 || std(x(valid)) <= sqrt(eps) || std(y(valid)) <= sqrt(eps)
    r = NaN;
    return
end
x = x(valid);
y = y(valid);
if rank_transform
    x = tied_rank_local(x);
    y = tied_rank_local(y);
end
C = corrcoef(x, y);
r = C(1, 2);
end

function [y, X, pair_labels, odor_ids] = displacement_data(D, Dp, Di, context_names)
context_pairs = nchoosek(1:4, 2);
y = nan(120, 1);
X = nan(120, 2);
pair_labels = strings(120, 1);
odor_ids = nan(120, 1);
row = 0;
for pair_index = 1:size(context_pairs, 1)
    c1 = context_pairs(pair_index, 1);
    c2 = context_pairs(pair_index, 2);
    pair_name = context_names(c1) + "_vs_" + context_names(c2);
    for odor_id = 1:20
        row = row + 1;
        i = (c1 - 1) * 20 + odor_id;
        j = (c2 - 1) * 20 + odor_id;
        y(row) = D(i, j);
        X(row, :) = [Dp(i, j), Di(i, j)];
        pair_labels(row) = pair_name;
        odor_ids(row) = odor_id;
    end
end
assert(row == 120, 'Same-odor displacement extraction did not yield 120 rows.');
end

function row = omnibus_row(subject_id, roi, metric, version, fit)
row = base_model_row(subject_id, roi, metric, version, fit);
row.beta_odor = fit.beta(1); row.beta_context = fit.beta(2);
row.beta_pleasantness = fit.beta(3); row.beta_intensity = fit.beta(4);
row.delta_r2_odor = fit.delta_r2(1); row.delta_r2_context = fit.delta_r2(2);
row.delta_r2_pleasantness = fit.delta_r2(3); row.delta_r2_intensity = fit.delta_r2(4);
row.partial_r2_odor = fit.partial_r2(1); row.partial_r2_context = fit.partial_r2(2);
row.partial_r2_pleasantness = fit.partial_r2(3); row.partial_r2_intensity = fit.partial_r2(4);
end

function row = context_rsa_row(subject_id, roi, metric, context, version, fit)
row = base_model_row(subject_id, roi, metric, version, fit);
row.context = context;
row.beta_pleasantness = fit.beta(1); row.beta_intensity = fit.beta(2);
row.delta_r2_pleasantness = fit.delta_r2(1); row.delta_r2_intensity = fit.delta_r2(2);
row.partial_r2_pleasantness = fit.partial_r2(1); row.partial_r2_intensity = fit.partial_r2(2);
end

function row = warping_row(subject_id, roi, metric, pair_name, version, fit)
row = base_model_row(subject_id, roi, metric, version, fit);
row.context_pair = pair_name;
row.beta_pleasantness = fit.beta(1); row.beta_intensity = fit.beta(2);
row.delta_r2_pleasantness = fit.delta_r2(1); row.delta_r2_intensity = fit.delta_r2(2);
row.partial_r2_pleasantness = fit.partial_r2(1); row.partial_r2_intensity = fit.partial_r2(2);
end

function row = displacement_row(subject_id, roi, metric, model, fit)
row = base_model_row(subject_id, roi, metric, "linear", fit);
row.model = model;
row.beta_pleasantness = fit.beta(1); row.beta_intensity = fit.beta(2);
end

function row = base_model_row(subject_id, roi, metric, version, fit)
row = struct('subject_id', subject_id, 'subject', "subj_" + string(subject_id), ...
    'roi', roi, 'neural_metric', metric, 'version', version, ...
    'n_total', fit.n_total, 'n_valid', fit.n_valid, 'n_dropped', fit.n_dropped, ...
    'design_rank', fit.design_rank, 'r2', fit.r2, ...
    'adjusted_r2', fit.adjusted_r2, 'status', fit.status);
end

function row = predictor_corr_row(subject_id, roi, metric, version, p1, p2, r, n)
row = struct('subject_id', subject_id, 'subject', "subj_" + string(subject_id), ...
    'roi', roi, 'neural_metric', metric, 'version', version, ...
    'predictor_1', p1, 'predictor_2', p2, 'correlation', r, 'n_valid', n);
end

function row = agreement_row(subject_id, roi, pearson_r, spearman_r, n_valid, ...
    prop_negative, cross_mean, cross_median, simple_mean, simple_median)
row = struct('subject_id', subject_id, 'subject', "subj_" + string(subject_id), ...
    'roi', roi, 'pearson_r', pearson_r, 'spearman_r', spearman_r, ...
    'n_valid', n_valid, 'proportion_crossnobis_negative', prop_negative, ...
    'mean_crossnobis', cross_mean, 'median_crossnobis', cross_median, ...
    'mean_simple', simple_mean, 'median_simple', simple_median);
end

function row = rdm_point_row(subject_id, roi, observation, i, j, simple, crossnobis)
row = struct('subject_id', subject_id, 'subject', "subj_" + string(subject_id), ...
    'roi', roi, 'observation', observation, 'condition_i', i, 'condition_j', j, ...
    'simple_distance', simple, 'crossnobis_distance', crossnobis);
end

function row = geometry_row(subject_id, roi, metric, source, c1, c2, pearson_r, spearman_r, n)
row = struct('subject_id', subject_id, 'subject', "subj_" + string(subject_id), ...
    'roi', roi, 'neural_metric', metric, 'source', source, ...
    'context_1', c1, 'context_2', c2, 'pearson_r', pearson_r, ...
    'spearman_r', spearman_r, 'n_valid', n);
end

function row = qc_row(subject_id, roi, metric, analysis, fit)
row = struct('subject_id', subject_id, 'subject', "subj_" + string(subject_id), ...
    'roi', roi, 'neural_metric', metric, 'analysis', analysis, ...
    'n_total', fit.n_total, 'n_valid', fit.n_valid, 'n_dropped', fit.n_dropped, ...
    'design_rank', fit.design_rank, 'status', fit.status);
end

function row = correlation_qc_row(subject_id, roi, metric, analysis, n_total, n_valid, r1, r2)
status = "ok";
if n_valid < 3
    status = "too_few_valid_rows";
elseif ~isfinite(r1) || ~isfinite(r2)
    status = "constant_or_invalid_correlation";
end
row = struct('subject_id', subject_id, 'subject', "subj_" + string(subject_id), ...
    'roi', roi, 'neural_metric', metric, 'analysis', analysis, ...
    'n_total', n_total, 'n_valid', n_valid, 'n_dropped', n_total - n_valid, ...
    'design_rank', NaN, 'status', status);
end

function rows = append_coefficient_rows(rows, subject_id, roi, metric, model, coefficients)
for index = 1:height(coefficients)
    new_row = struct('subject_id', subject_id, 'subject', "subj_" + string(subject_id), ...
        'roi', roi, 'neural_metric', metric, 'model', model, ...
        'term', string(coefficients.term(index)), ...
        'coefficient', coefficients.coefficient(index));
    rows(end+1) = new_row; %#ok<AGROW>
end
end

function table_out = summarize_warping(warping_table)
keys = unique(warping_table(:, {'subject_id','subject','roi','neural_metric','version'}), 'rows');
template = struct('subject_id', 0, 'subject', "", 'roi', "", 'neural_metric', "", ...
    'version', "", 'mean_beta_pleasantness', NaN, 'median_beta_pleasantness', NaN, ...
    'mean_beta_intensity', NaN, 'median_beta_intensity', NaN, 'mean_r2', NaN, ...
    'n_context_pairs', 0, 'n_valid_fits', 0);
rows = repmat(template, 0, 1);
for index = 1:height(keys)
    mask = warping_table.subject_id == keys.subject_id(index) & ...
        warping_table.roi == keys.roi(index) & ...
        warping_table.neural_metric == keys.neural_metric(index) & ...
        warping_table.version == keys.version(index);
    rows(end+1) = struct('subject_id', keys.subject_id(index), ... %#ok<AGROW>
        'subject', keys.subject(index), 'roi', keys.roi(index), ...
        'neural_metric', keys.neural_metric(index), 'version', keys.version(index), ...
        'mean_beta_pleasantness', mean(warping_table.beta_pleasantness(mask), 'omitnan'), ...
        'median_beta_pleasantness', median(warping_table.beta_pleasantness(mask), 'omitnan'), ...
        'mean_beta_intensity', mean(warping_table.beta_intensity(mask), 'omitnan'), ...
        'median_beta_intensity', median(warping_table.beta_intensity(mask), 'omitnan'), ...
        'mean_r2', mean(warping_table.r2(mask), 'omitnan'), ...
        'n_context_pairs', nnz(mask), ...
        'n_valid_fits', nnz(warping_table.status(mask) == "ok"));
end
table_out = struct2table(rows);
end

function matrix_table = build_geometry_matrix_table(neural_table, behavior_table, contexts)
template = struct('subject_id', 0, 'subject', "", 'roi', "", ...
    'neural_metric', "", 'source', "", 'correlation', "", ...
    'context_row', "", 'context_column', "", 'row_index', 0, ...
    'column_index', 0, 'similarity', NaN, 'is_across_subject_average', false);
rows = repmat(template, 0, 1);
neural_keys = unique(neural_table(:, {'subject_id','subject','roi','neural_metric','source'}), 'rows');
for key_index = 1:height(neural_keys)
    T = neural_table(neural_table.subject_id == neural_keys.subject_id(key_index) & ...
        neural_table.roi == neural_keys.roi(key_index) & ...
        neural_table.neural_metric == neural_keys.neural_metric(key_index), :);
    for correlation = ["pearson", "spearman"]
        M = geometry_matrix(T, correlation + "_r", contexts);
        rows = append_matrix_rows(rows, neural_keys.subject_id(key_index), ...
            neural_keys.subject(key_index), neural_keys.roi(key_index), ...
            neural_keys.neural_metric(key_index), "neural", correlation, M, contexts, false);
    end
end
behavior_keys = unique(behavior_table(:, {'subject_id','subject','source'}), 'rows');
for key_index = 1:height(behavior_keys)
    T = behavior_table(behavior_table.subject_id == behavior_keys.subject_id(key_index) & ...
        behavior_table.source == behavior_keys.source(key_index), :);
    for correlation = ["pearson", "spearman"]
        M = geometry_matrix(T, correlation + "_r", contexts);
        rows = append_matrix_rows(rows, behavior_keys.subject_id(key_index), ...
            behavior_keys.subject(key_index), "", "", behavior_keys.source(key_index), ...
            correlation, M, contexts, false);
    end
end
individual = struct2table(rows);
group_keys = unique(individual(:, {'roi','neural_metric','source','correlation'}), 'rows');
for key_index = 1:height(group_keys)
    T = individual(individual.roi == group_keys.roi(key_index) & ...
        individual.neural_metric == group_keys.neural_metric(key_index) & ...
        individual.source == group_keys.source(key_index) & ...
        individual.correlation == group_keys.correlation(key_index), :);
    M = nan(4);
    for row_index = 1:4
        for column_index = 1:4
            mask = T.row_index == row_index & T.column_index == column_index;
            M(row_index,column_index) = fisher_mean(T.similarity(mask));
        end
    end
    rows = append_matrix_rows(rows, 0, "across_subject_fisher_mean", ...
        group_keys.roi(key_index), group_keys.neural_metric(key_index), ...
        group_keys.source(key_index), group_keys.correlation(key_index), M, contexts, true);
end
matrix_table = struct2table(rows);
end

function rows = append_matrix_rows(rows, subject_id, subject, roi, metric, source, correlation, M, contexts, is_group)
for row_index = 1:4
    for column_index = 1:4
        rows(end+1) = struct('subject_id', subject_id, 'subject', subject, ... %#ok<AGROW>
            'roi', roi, 'neural_metric', metric, 'source', source, ...
            'correlation', correlation, 'context_row', contexts(row_index), ...
            'context_column', contexts(column_index), 'row_index', row_index, ...
            'column_index', column_index, 'similarity', M(row_index,column_index), ...
            'is_across_subject_average', is_group);
    end
end
end

function value = fisher_mean(values)
values = double(values(:));
values = values(isfinite(values));
if isempty(values)
    value = NaN;
    return
end
bound = 1 - 1e-12;
values = min(max(values, -bound), bound);
value = tanh(mean(atanh(values)));
end

function summary = build_summary_table(omnibus, agreement, geometry, context_rsa, warping, displacement)
subject_ids = unique(omnibus.subject_id, 'stable')';
rois = unique(omnibus.roi, 'stable')';
metrics = unique(omnibus.neural_metric, 'stable')';
analyses = ["1_omnibus", "2_metric_concordance", "3_context_geometry", ...
    "4_context_perceptual_RSA", "5_geometry_warping", "6_same_odor_displacement"];
rows = repmat(summary_row_template(), 0, 1);
for subject_id = subject_ids
    for roi = rois
        for metric = metrics
            for analysis = analyses
                row = summary_row_template();
                row.subject_id = subject_id; row.subject = "subj_" + string(subject_id);
                row.roi = roi; row.neural_metric = metric; row.analysis = analysis;
                switch analysis
                    case "1_omnibus"
                        T = omnibus(omnibus.subject_id == subject_id & omnibus.roi == roi & ...
                            omnibus.neural_metric == metric, :);
                        L = T(T.version == "linear", :); R = T(T.version == "rank", :);
                        row.linear_beta_odor = scalar_or_nan(L.beta_odor);
                        row.linear_beta_context = scalar_or_nan(L.beta_context);
                        row.linear_beta_pleasantness = scalar_or_nan(L.beta_pleasantness);
                        row.linear_beta_intensity = scalar_or_nan(L.beta_intensity);
                        row.linear_r2 = scalar_or_nan(L.r2); row.linear_adjusted_r2 = scalar_or_nan(L.adjusted_r2);
                        row.rank_beta_odor = scalar_or_nan(R.beta_odor);
                        row.rank_beta_context = scalar_or_nan(R.beta_context);
                        row.rank_beta_pleasantness = scalar_or_nan(R.beta_pleasantness);
                        row.rank_beta_intensity = scalar_or_nan(R.beta_intensity);
                        row.rank_r2 = scalar_or_nan(R.r2); row.rank_adjusted_r2 = scalar_or_nan(R.adjusted_r2);
                        row.n_models_expected = 2; row.n_models_valid = nnz(T.status == "ok");
                        row.min_n_valid = min(T.n_valid); row.status = join_status(T.status);
                    case "2_metric_concordance"
                        T = agreement(agreement.subject_id == subject_id & agreement.roi == roi, :);
                        row.pearson_r = scalar_or_nan(T.pearson_r); row.spearman_r = scalar_or_nan(T.spearman_r);
                        row.proportion_crossnobis_negative = scalar_or_nan(T.proportion_crossnobis_negative);
                        row.mean_crossnobis = scalar_or_nan(T.mean_crossnobis);
                        row.median_crossnobis = scalar_or_nan(T.median_crossnobis);
                        row.n_models_expected = 1;
                        row.n_models_valid = nnz(isfinite(T.pearson_r) & isfinite(T.spearman_r));
                        row.min_n_valid = scalar_or_nan(T.n_valid); row.status = status_from_count(row.n_models_valid, 1);
                    case "3_context_geometry"
                        T = geometry(geometry.subject_id == subject_id & geometry.roi == roi & ...
                            geometry.neural_metric == metric, :);
                        row.n_models_expected = 6; row.n_models_valid = nnz(isfinite(T.pearson_r) & isfinite(T.spearman_r));
                        if row.n_models_valid == row.n_models_expected
                            row.mean_pearson_r = fisher_mean(T.pearson_r);
                            row.mean_spearman_r = fisher_mean(T.spearman_r);
                        end
                        row.min_n_valid = min(T.n_valid); row.status = status_from_count(row.n_models_valid, 6);
                    case "4_context_perceptual_RSA"
                        T = context_rsa(context_rsa.subject_id == subject_id & context_rsa.roi == roi & ...
                            context_rsa.neural_metric == metric, :);
                        for version = ["linear", "rank"]
                            V = T(T.version == version, :);
                            if height(V) == 4 && all(V.status == "ok")
                                row.(version + "_mean_beta_pleasantness") = mean(V.beta_pleasantness);
                                row.(version + "_mean_beta_intensity") = mean(V.beta_intensity);
                                row.(version + "_r2") = mean(V.r2);
                                row.(version + "_adjusted_r2") = mean(V.adjusted_r2);
                            end
                        end
                        row.n_models_expected = 8; row.n_models_valid = nnz(T.status == "ok");
                        row.min_n_valid = min(T.n_valid); row.status = join_status(T.status);
                    case "5_geometry_warping"
                        T = warping(warping.subject_id == subject_id & warping.roi == roi & ...
                            warping.neural_metric == metric, :);
                        for version = ["linear", "rank"]
                            V = T(T.version == version, :);
                            if height(V) == 6 && all(V.status == "ok")
                                row.(version + "_mean_beta_pleasantness") = mean(V.beta_pleasantness);
                                row.(version + "_median_beta_pleasantness") = median(V.beta_pleasantness);
                                row.(version + "_mean_beta_intensity") = mean(V.beta_intensity);
                                row.(version + "_median_beta_intensity") = median(V.beta_intensity);
                                row.(version + "_r2") = mean(V.r2);
                                row.(version + "_adjusted_r2") = mean(V.adjusted_r2);
                            end
                        end
                        row.n_models_expected = 12; row.n_models_valid = nnz(T.status == "ok");
                        row.min_n_valid = min(T.n_valid); row.status = join_status(T.status);
                    case "6_same_odor_displacement"
                        T = displacement(displacement.subject_id == subject_id & displacement.roi == roi & ...
                            displacement.neural_metric == metric, :);
                        B = T(T.model == "context_pair", :); O = T(T.model == "context_pair_plus_odor", :);
                        row.linear_beta_pleasantness = scalar_or_nan(B.beta_pleasantness);
                        row.linear_beta_intensity = scalar_or_nan(B.beta_intensity);
                        row.linear_r2 = scalar_or_nan(B.r2); row.linear_adjusted_r2 = scalar_or_nan(B.adjusted_r2);
                        row.odor_fe_beta_pleasantness = scalar_or_nan(O.beta_pleasantness);
                        row.odor_fe_beta_intensity = scalar_or_nan(O.beta_intensity);
                        row.odor_fe_r2 = scalar_or_nan(O.r2); row.odor_fe_adjusted_r2 = scalar_or_nan(O.adjusted_r2);
                        row.n_models_expected = 2; row.n_models_valid = nnz(T.status == "ok");
                        row.min_n_valid = min(T.n_valid); row.status = join_status(T.status);
                end
                rows(end+1) = row; %#ok<AGROW>
            end
        end
    end
end
summary = struct2table(rows);
assert(height(summary) == numel(subject_ids) * numel(rois) * numel(metrics) * 6, ...
    'Master summary does not have one row per subject x ROI x metric x analysis.');
end

function row = summary_row_template()
row = struct('subject_id', 0, 'subject', "", 'roi', "", 'neural_metric', "", ...
    'analysis', "", 'linear_beta_odor', NaN, 'linear_beta_context', NaN, ...
    'linear_beta_pleasantness', NaN, 'linear_beta_intensity', NaN, ...
    'linear_mean_beta_pleasantness', NaN, 'linear_median_beta_pleasantness', NaN, ...
    'linear_mean_beta_intensity', NaN, 'linear_median_beta_intensity', NaN, ...
    'linear_r2', NaN, 'linear_adjusted_r2', NaN, ...
    'rank_beta_odor', NaN, 'rank_beta_context', NaN, ...
    'rank_beta_pleasantness', NaN, 'rank_beta_intensity', NaN, ...
    'rank_mean_beta_pleasantness', NaN, 'rank_median_beta_pleasantness', NaN, ...
    'rank_mean_beta_intensity', NaN, 'rank_median_beta_intensity', NaN, ...
    'rank_r2', NaN, 'rank_adjusted_r2', NaN, 'pearson_r', NaN, ...
    'spearman_r', NaN, 'mean_pearson_r', NaN, 'mean_spearman_r', NaN, ...
    'proportion_crossnobis_negative', NaN, 'mean_crossnobis', NaN, ...
    'median_crossnobis', NaN, 'odor_fe_beta_pleasantness', NaN, ...
    'odor_fe_beta_intensity', NaN, 'odor_fe_r2', NaN, ...
    'odor_fe_adjusted_r2', NaN, 'n_models_expected', 0, ...
    'n_models_valid', 0, 'min_n_valid', NaN, 'status', "not_evaluated");
end

function value = scalar_or_nan(values)
if isempty(values), value = NaN; else, value = values(1); end
end

function status = status_from_count(actual, expected)
if actual == expected, status = "ok"; else, status = "incomplete:" + actual + "_of_" + expected; end
end

function status = join_status(values)
bad = unique(values(values ~= "ok"));
if isempty(bad), status = "ok"; else, status = strjoin(bad, '|'); end
end

function write_table(T, output_root, stem)
writetable(T, fullfile(output_root, [stem '.csv']));
save(fullfile(output_root, [stem '.mat']), 'T');
end

function make_all_plots(root, omnibus, predictor_corr, agreement, points, geometry, ...
    behavior_geometry, context_rsa, warping, displacement, displacement_observations, ...
    subject_ids, rois, metrics, contexts)
colors = lines(numel(subject_ids));
% Analysis 1: separate neural metric and model version.
for metric = metrics
    for version = ["linear", "rank"]
        subset = omnibus(omnibus.neural_metric == metric & omnibus.version == version, :);
        f = new_figure(sprintf('Omnibus beta: %s %s', metric, version), [1200 760]);
        beta_vars = {'beta_odor','beta_context','beta_pleasantness','beta_intensity'};
        beta_labels = {'Odor identity (same=1)','Context identity (same=1)','Pleasantness','Intensity'};
        t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        for p = 1:4
            ax = nexttile(t); hold(ax, 'on');
            plot_subject_points(ax, subset, beta_vars{p}, rois, subject_ids, colors);
            yline(ax, 0, ':k'); title(ax, beta_labels{p}); ylabel(ax, 'Standardized beta');
        end
        title(t, sprintf('Omnibus multiple-regression RSA | %s | %s', metric, version));
        save_figure(f, fullfile(root, 'analysis1', sprintf('omnibus_betas_%s_%s', metric, version)));
    end
end
f = new_figure('Predictor correlations', [1050 480]);
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for v = 1:2
    version = ["linear", "rank"]; version = version(v);
    subset = predictor_corr(predictor_corr.version == version, :);
    names = ["odor","context","pleasantness","intensity"];
    matrix = nan(4);
    for i = 1:4, for j = 1:4
        mask = subset.predictor_1 == names(i) & subset.predictor_2 == names(j);
        matrix(i,j) = mean(subset.correlation(mask), 'omitnan');
    end, end
    ax = nexttile(t); imagesc(ax, matrix, [-1 1]); axis(ax, 'square'); colorbar(ax);
    ax.XTick = 1:4; ax.YTick = 1:4; ax.XTickLabel = names; ax.YTickLabel = names;
    xtickangle(ax, 30); title(ax, version + " predictors (subject mean)");
end
title(t, 'Predictor correlation matrices');
save_figure(f, fullfile(root, 'analysis1', 'predictor_correlation_heatmap'));

% Analysis 2.
f = new_figure('RDM agreement', [1150 480]);
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:2
    ax = nexttile(t); hold(ax, 'on');
    variable = ["pearson_r", "spearman_r"]; variable = variable(k);
    plot_subject_points(ax, agreement, char(variable), rois, subject_ids, colors);
    ylabel(ax, variable + " correlation"); title(ax, upper(extractBefore(variable, '_')));
end
title(t, 'Upper-triangle agreement: simple 1-r vs crossnobis');
save_figure(f, fullfile(root, 'analysis2', 'RDM_agreement_by_ROI'));
f = new_figure('RDM scatter by subject', [1500 1200]);
t = tiledlayout(numel(rois), numel(subject_ids), 'TileSpacing', 'compact', 'Padding', 'compact');
for r = 1:numel(rois), for s = 1:numel(subject_ids)
    ax = nexttile(t); mask = points.roi == rois(r) & points.subject_id == subject_ids(s);
    scatter(ax, points.simple_distance(mask), points.crossnobis_distance(mask), 5, '.', ...
        'MarkerEdgeAlpha', .18); yline(ax, 0, ':k');
    if r == 1, title(ax, "subj_" + subject_ids(s)); end
    if s == 1, ylabel(ax, rois(r)); end
    if r == numel(rois), xlabel(ax, '1-r'); end
end, end
title(t, 'Neural RDM scatterplots (all 3,160 distances per panel)');
save_figure(f, fullfile(root, 'analysis2', 'simple_vs_crossnobis_scatter_subject_ROI'));
f = new_figure('Distance distributions', [1000 1050]);
t = tiledlayout(numel(rois), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for r = 1:numel(rois)
    mask = points.roi == rois(r);
    ax = nexttile(t); histogram(ax, points.simple_distance(mask), 45, 'Normalization', 'probability');
    ylabel(ax, rois(r)); if r == 1, title(ax, 'Simple 1-r'); end
    ax = nexttile(t); histogram(ax, points.crossnobis_distance(mask), 45, 'Normalization', 'probability');
    xline(ax, 0, ':k'); if r == 1, title(ax, 'Crossnobis'); end
end
title(t, 'Distance distributions pooled across subjects');
save_figure(f, fullfile(root, 'analysis2', 'distance_distributions_by_ROI'));
for version = ["linear", "rank"]
    f = new_figure('Beta agreement', [1150 760]);
    t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    beta_vars = {'beta_odor','beta_context','beta_pleasantness','beta_intensity'};
    for p = 1:4
        ax = nexttile(t); hold(ax, 'on');
        for r = 1:numel(rois)
            a = omnibus(omnibus.roi == rois(r) & omnibus.version == version, :);
            xs = nan(numel(subject_ids),1); ys = xs;
            for s = 1:numel(subject_ids)
                xs(s) = a{a.subject_id == subject_ids(s) & a.neural_metric == "simple", beta_vars{p}};
                ys(s) = a{a.subject_id == subject_ids(s) & a.neural_metric == "crossnobis", beta_vars{p}};
            end
            scatter(ax, xs, ys, 36, repmat(r, size(xs)), 'filled', 'DisplayName', rois(r));
        end
        add_identity_line(ax); xlabel(ax, 'Simple beta'); ylabel(ax, 'Crossnobis beta');
        title(ax, erase(beta_vars{p}, 'beta_'));
        if p == 1
            colormap(ax, lines(numel(rois))); clim(ax, [1 numel(rois)]);
            legend(ax, 'Location', 'best', 'Interpreter', 'none');
        end
    end
    title(t, "Omnibus beta agreement | " + version);
    save_figure(f, fullfile(root, 'analysis2', "omnibus_beta_agreement_" + version));
end

% Analysis 3 heatmaps: six panels (five subjects plus mean).
for metric = metrics
    for correlation = ["pearson", "spearman"]
        field = correlation + "_r";
        for r = 1:numel(rois)
            subset = geometry(geometry.roi == rois(r) & geometry.neural_metric == metric, :);
            plot_geometry_panels(subset, field, contexts, subject_ids, ...
                "Neural " + rois(r) + " | " + metric + " | " + correlation, ...
                fullfile(root, 'analysis3', "neural_geometry_" + rois(r) + "_" + metric + "_" + correlation));
        end
    end
end
for correlation = ["pearson", "spearman"]
    field = correlation + "_r";
    f = new_figure('ROI geometry comparison', [1100 520]);
    t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:2
        current_metric = metrics(k);
        ax = nexttile(t); hold(ax, 'on');
        subset = geometry(geometry.neural_metric == current_metric, :);
        plot_geometry_subject_means(ax, subset, char(field), rois, subject_ids, colors);
        title(ax, current_metric); ylabel(ax, 'Mean cross-context similarity');
    end
    title(t, "Context-geometry preservation | " + correlation);
    save_figure(f, fullfile(root, 'analysis3', "ROI_mean_cross_context_" + correlation));
end
for measure = ["pleasantness", "intensity"]
    for correlation = ["pearson", "spearman"]
        subset = behavior_geometry(behavior_geometry.source == measure, :);
        plot_geometry_panels(subset, correlation + "_r", contexts, subject_ids, ...
            upper(extractBefore(measure, 2)) + extractAfter(measure, 1) + " geometry | " + correlation, ...
            fullfile(root, 'analysis3', measure + "_geometry_" + correlation));
    end
end

% Analysis 4 trajectories, with both neural metrics in each panel.
for version = ["linear", "rank"]
    for beta = ["pleasantness", "intensity"]
        f = new_figure('Context-specific RSA', [1250 800]);
        t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        for r = 1:numel(rois)
            ax = nexttile(t); hold(ax, 'on');
            for m = 1:2, for s = 1:numel(subject_ids)
                mask = context_rsa.roi == rois(r) & context_rsa.version == version & ...
                    context_rsa.neural_metric == metrics(m) & context_rsa.subject_id == subject_ids(s);
                values = order_context_values(context_rsa(mask,:), "beta_" + beta, contexts);
                style = '-'; if m == 2, style = '--'; end
                plot(ax, 1:4, values, style, 'Color', [colors(s,:) .55], 'LineWidth', 1);
            end, end
            yline(ax, 0, ':k'); ax.XTick = 1:4; ax.XTickLabel = contexts; xtickangle(ax, 25);
            title(ax, rois(r)); ylabel(ax, 'Standardized beta');
        end
        key_ax = nexttile(t, 6); hold(key_ax, 'on');
        h1 = plot(key_ax, nan, nan, '-k', 'LineWidth', 1.5);
        h2 = plot(key_ax, nan, nan, '--k', 'LineWidth', 1.5);
        legend(key_ax, [h1 h2], {'simple','crossnobis'}, 'Location', 'northwest'); axis(key_ax, 'off');
        title(t, beta + " beta by context | " + version + " | individual trajectories");
        save_figure(f, fullfile(root, 'analysis4', "context_" + beta + "_beta_" + version));
    end
end

% Analysis 5.
pair_order = unique(warping.context_pair, 'stable');
for version = ["linear", "rank"]
    for metric = metrics
        f = new_figure('Warping by pair', [1250 800]);
        t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        for r = 1:numel(rois)
            ax = nexttile(t); hold(ax, 'on');
            for s = 1:numel(subject_ids)
                mask = warping.roi == rois(r) & warping.version == version & ...
                    warping.neural_metric == metric & warping.subject_id == subject_ids(s);
                p = order_pair_values(warping(mask,:), 'beta_pleasantness', pair_order);
                q = order_pair_values(warping(mask,:), 'beta_intensity', pair_order);
                plot(ax, 1:6, p, '-', 'Color', [colors(s,:) .55]);
                plot(ax, 1:6, q, '--', 'Color', [colors(s,:) .55]);
            end
            yline(ax, 0, ':k'); ax.XTick = 1:6; ax.XTickLabel = pair_order; xtickangle(ax, 35);
            ax.TickLabelInterpreter = 'none';
            title(ax, rois(r)); ylabel(ax, 'Standardized beta');
        end
        key_ax = nexttile(t, 6); hold(key_ax, 'on');
        h1 = plot(key_ax, nan, nan, '-k', 'LineWidth', 1.5);
        h2 = plot(key_ax, nan, nan, '--k', 'LineWidth', 1.5);
        legend(key_ax, [h1 h2], {'pleasantness','intensity'}, 'Location', 'northwest'); axis(key_ax, 'off');
        title(t, "Geometry warping by context pair | " + metric + " | " + version);
        save_figure(f, fullfile(root, 'analysis5', "warping_by_pair_" + metric + "_" + version));
    end
    f = new_figure('Mean warping by ROI', [1150 480]);
    t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for b = 1:2
        ax = nexttile(t); hold(ax, 'on'); beta_var = ["beta_pleasantness","beta_intensity"]; beta_var=beta_var(b);
        plot_warping_means(ax, warping(warping.version == version,:), char(beta_var), rois, metrics, subject_ids);
        if b == 1, add_metric_legend(ax); end
        yline(ax,0,':k'); title(ax, erase(beta_var,'beta_')); ylabel(ax,'Mean beta across six pairs');
    end
    title(t, "Mean geometry warping by ROI | " + version);
    save_figure(f, fullfile(root, 'analysis5', "warping_mean_by_ROI_" + version));
    f = new_figure('Pleasantness vs intensity warping', [1150 480]);
    t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for m = 1:2
        ax = nexttile(t); hold(ax, 'on'); subset=warping(warping.version==version & warping.neural_metric==metrics(m),:);
        scatter(ax, subset.beta_pleasantness, subset.beta_intensity, 18, subset.subject_id, 'filled', 'MarkerFaceAlpha',.55);
        xline(ax,0,':k'); yline(ax,0,':k'); add_identity_line(ax); axis(ax,'square');
        xlabel(ax,'Pleasantness beta'); ylabel(ax,'Intensity beta'); title(ax,metrics(m));
    end
    title(t, "Pleasantness-driven vs intensity-driven warping | " + version);
    save_figure(f, fullfile(root, 'analysis5', "warping_pleasantness_vs_intensity_" + version));
end

% Analysis 6 raw scatterplots and model coefficients.
pair_levels = unique(displacement_observations.context_pair, 'stable');
pair_colors = lines(numel(pair_levels));
for metric = metrics
    for predictor = ["pleasantness", "intensity"]
        xvar = "abs_delta_" + predictor;
        f = new_figure('Same-odor displacement scatter', [1250 800]);
        t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
        for r=1:numel(rois)
            ax=nexttile(t); hold(ax,'on');
            for p=1:numel(pair_levels)
                mask=displacement_observations.roi==rois(r) & ...
                    displacement_observations.neural_metric==metric & ...
                    displacement_observations.context_pair==pair_levels(p);
                scatter(ax,displacement_observations{mask,char(xvar)}, ...
                    displacement_observations.neural_distance(mask),12, ...
                    pair_colors(p,:),'filled','MarkerFaceAlpha',.38, ...
                    'DisplayName',pair_levels(p));
            end
            yline(ax,0,':k'); xlabel(ax,"Absolute change in " + predictor);
            ylabel(ax,'Neural distance'); title(ax,rois(r));
        end
        key_ax = nexttile(t,6); hold(key_ax,'on'); handles=gobjects(numel(pair_levels),1);
        for p=1:numel(pair_levels)
            handles(p)=scatter(key_ax,nan,nan,25,pair_colors(p,:),'filled');
        end
        legend(key_ax,handles,pair_levels,'Interpreter','none','Location','northwest'); axis(key_ax,'off');
        title(t,"Same-odor neural displacement vs absolute " + predictor + " change | " + metric);
        save_figure(f,fullfile(root,'analysis6',"displacement_vs_"+predictor+"_"+metric));
    end
end
for model = ["context_pair", "context_pair_plus_odor"]
    f = new_figure('Same-odor displacement coefficients', [1150 480]);
    t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    for b=1:2
        ax=nexttile(t); hold(ax,'on'); beta_var=["beta_pleasantness","beta_intensity"]; beta_var=beta_var(b);
        subset=displacement(displacement.model==model,:);
        plot_metric_subject_points(ax,subset,char(beta_var),rois,metrics,subject_ids);
        if b == 1, add_metric_legend(ax); end
        yline(ax,0,':k'); title(ax,erase(beta_var,'beta_')); ylabel(ax,'Standardized coefficient');
    end
    title(t,"Same-odor displacement coefficients | " + replace(model,'_',' '));
    save_figure(f,fullfile(root,'analysis6',"displacement_coefficients_"+model));
end
end

function plot_subject_points(ax, T, variable, rois, subject_ids, colors)
for r = 1:numel(rois)
    values = nan(numel(subject_ids), 1);
    for s = 1:numel(subject_ids)
        mask = T.roi == rois(r) & T.subject_id == subject_ids(s);
        values(s) = mean(T{mask, variable}, 'omitnan');
    end
    scatter(ax, r + linspace(-.13,.13,numel(subject_ids))', values, 35, colors, 'filled');
    plot(ax, [r-.2 r+.2], repmat(mean(values,'omitnan'),1,2), 'k-', 'LineWidth', 2);
end
ax.XTick = 1:numel(rois); ax.XTickLabel = rois; xtickangle(ax, 25); xlim(ax,[.5 numel(rois)+.5]);
end

function plot_geometry_subject_means(ax,T,variable,rois,subject_ids,colors)
for r=1:numel(rois)
    vals=nan(numel(subject_ids),1);
    for s=1:numel(subject_ids)
        mask=T.roi==rois(r)&T.subject_id==subject_ids(s); vals(s)=fisher_mean(T{mask,variable});
    end
    scatter(ax,r+linspace(-.13,.13,numel(vals))',vals,35,colors,'filled');
    plot(ax,[r-.2 r+.2],repmat(mean(vals,'omitnan'),1,2),'k-','LineWidth',2);
end
ax.XTick=1:numel(rois); ax.XTickLabel=rois; xtickangle(ax,25); xlim(ax,[.5 numel(rois)+.5]); yline(ax,0,':k');
end

function plot_geometry_panels(T, variable, contexts, subject_ids, plot_title, output_stem)
f = new_figure(plot_title, [1200 720]); t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
matrices = nan(4,4,numel(subject_ids));
for s=1:numel(subject_ids)
    matrices(:,:,s)=geometry_matrix(T(T.subject_id==subject_ids(s),:),variable,contexts);
    ax=nexttile(t); imagesc(ax,matrices(:,:,s),[-1 1]); axis(ax,'square');
    ax.XTick=1:4; ax.YTick=1:4; ax.XTickLabel=contexts; ax.YTickLabel=contexts; xtickangle(ax,25);
    title(ax,"subj_"+subject_ids(s));
end
group_matrix = nan(4);
for row_index=1:4, for column_index=1:4
    group_matrix(row_index,column_index)=fisher_mean(squeeze(matrices(row_index,column_index,:)));
end, end
ax=nexttile(t); imagesc(ax,group_matrix,[-1 1]); axis(ax,'square'); colorbar(ax);
ax.XTick=1:4; ax.YTick=1:4; ax.XTickLabel=contexts; ax.YTickLabel=contexts; xtickangle(ax,25); title(ax,'Across-subject mean');
title(t,plot_title); save_figure(f,output_stem);
end

function M=geometry_matrix(T,variable,contexts)
M=eye(4);
for i=1:height(T)
    a=find(contexts==T.context_1(i)); b=find(contexts==T.context_2(i));
    if isscalar(a)&&isscalar(b), M(a,b)=T{i,char(variable)}; M(b,a)=M(a,b); end
end
end

function values=order_context_values(T,variable,contexts)
values=nan(1,numel(contexts));
for k=1:numel(contexts), mask=T.context==contexts(k); if any(mask), values(k)=T{find(mask,1),char(variable)}; end, end
end

function values=order_pair_values(T,variable,pairs)
values=nan(1,numel(pairs));
for k=1:numel(pairs), mask=T.context_pair==pairs(k); if any(mask), values(k)=T{find(mask,1),variable}; end, end
end

function plot_warping_means(ax,T,variable,rois,metrics,subject_ids)
offset=[-.13 .13]; symbols={'o','s'};
for r=1:numel(rois), for m=1:2
    vals=nan(numel(subject_ids),1);
    for s=1:numel(subject_ids)
        mask=T.roi==rois(r)&T.neural_metric==metrics(m)&T.subject_id==subject_ids(s);
        vals(s)=mean(T{mask,variable},'omitnan');
    end
    scatter(ax,repmat(r+offset(m),size(vals)),vals,34,symbols{m},'filled','DisplayName',metrics(m));
end,end
ax.XTick=1:numel(rois); ax.XTickLabel=rois; xtickangle(ax,25); xlim(ax,[.5 numel(rois)+.5]);
end

function plot_metric_subject_points(ax,T,variable,rois,metrics,subject_ids)
offset=[-.13 .13]; symbols={'o','s'};
for r=1:numel(rois), for m=1:2
    vals=nan(numel(subject_ids),1);
    for s=1:numel(subject_ids)
        mask=T.roi==rois(r)&T.neural_metric==metrics(m)&T.subject_id==subject_ids(s);
        vals(s)=mean(T{mask,variable},'omitnan');
    end
    scatter(ax,repmat(r+offset(m),size(vals)),vals,34,symbols{m},'filled');
end,end
ax.XTick=1:numel(rois); ax.XTickLabel=rois; xtickangle(ax,25); xlim(ax,[.5 numel(rois)+.5]);
end

function add_metric_legend(ax)
h1 = scatter(ax, nan, nan, 34, 'o', 'filled', 'MarkerFaceColor', [.25 .25 .25]);
h2 = scatter(ax, nan, nan, 34, 's', 'filled', 'MarkerFaceColor', [.25 .25 .25]);
legend_handle = legend(ax, [h1 h2], {'simple','crossnobis'}, ...
    'Location', 'best', 'Interpreter', 'none');
legend_handle.AutoUpdate = 'off';
end

function add_identity_line(ax)
limits=[min([ax.XLim ax.YLim]) max([ax.XLim ax.YLim])]; plot(ax,limits,limits,':k'); xlim(ax,limits); ylim(ax,limits);
end

function f=new_figure(name,dimensions)
f=figure('Visible','off','Color','w','Name',name,'Position',[50 50 dimensions]);
end

function save_figure(f,stem)
exportgraphics(f,[char(stem) '.png'],'Resolution',180); close(f);
end

function ensure_directory(path)
if ~isfolder(path), mkdir(path); end
end

function write_readme(output_root, categorical_file, n_warnings, n_qc, n_input_warnings, make_plots)
readme = fullfile(output_root, 'README.md');
fid = fopen(readme, 'w'); assert(fid >= 0, 'Could not create %s.', readme);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '# RSA first-batch exploratory analysis\n\n');
fprintf(fid, 'Generated by `scripts/run_exploratory_rsa_first_batch.m`. These are descriptive/exploratory outputs; no group inferential tests were run.\n\n');
fprintf(fid, '## Condition ordering and inputs\n\n');
fprintf(fid, '- Subjects: `subj_2` through `subj_6`.\n');
fprintf(fid, '- ROIs: AON, PirF, PirT, olfAMG, olfOFC.\n');
fprintf(fid, '- Neural metrics: simple correlation distance (`1-r`) and crossnobis. Negative crossnobis distances are retained.\n');
fprintf(fid, '- Order: PERSON odors 1-20, FOOD 1-20, LOCATION 1-20, CONTROL 1-20. See `condition_ordering.csv`.\n');
fprintf(fid, '- Neural input: `RDMs/neural/subj_N_neural_RDMs.mat`; behavioral input: `RDMs/subj_N_behavioral_RDMs.mat`; categorical input: `%s`.\n', categorical_file);
fprintf(fid, '- The saved categorical predictors are similarity-coded: 1=same odor/context, 0=different. They are preserved exactly, so positive categorical beta means larger neural distance for same-identity pairs after adjustment.\n\n');
fprintf(fid, '## What the analyses test\n\n');
fprintf(fid, '1. Omnibus RSA asks whether odor identity, context identity, pleasantness distance, and intensity distance explain unique neural-distance variance.\n');
fprintf(fid, '2. Metric concordance checks whether simple and crossnobis RDMs and their omnibus betas agree.\n');
fprintf(fid, '3. Context geometry correlates the six pairs of within-context 20-by-20 odor geometries for neural and perceptual RDMs.\n');
fprintf(fid, '4. Context-specific perceptual RSA estimates pleasantness and intensity effects within each context.\n');
fprintf(fid, '5. Geometry warping asks whether context-pair changes in perceptual odor distances track changes in neural odor distances.\n');
fprintf(fid, '6. Same-odor displacement asks whether the neural distance between the same odor across contexts scales with absolute rating change, controlling context pair, with a secondary odor-fixed-effect model.\n\n');
fprintf(fid, '## Estimation and files\n\n');
fprintf(fid, '- Analysis 6 follow-up: simple-distance scatterplots include separate pooled, univariate linear fits for each context pair. Raw slopes and standardized betas are saved separately in `tables/analysis6_context_pair_slopes.csv`, `.mat`, and a compact `.md` table. These descriptive fits do not adjust for subject, odor, or the other rating. Rerun independently with `run_exploratory_rsa_analysis6_pair_slopes`.\n');
fprintf(fid, '- Linear models z-score continuous outcome and predictors. Rank models use average tied ranks followed by z-scoring. All fits include an intercept.\n');
fprintf(fid, '- Delta-R2 is full-model R2 minus the corresponding reduced-model R2 on identical valid rows. Partial R2 is `(SSE_reduced-SSE_full)/SSE_reduced`.\n');
fprintf(fid, '- Correlation summaries and across-subject context-geometry matrices use Fisher-z averaging before conversion back to r.\n');
fprintf(fid, '- `tables/` contains tidy CSV plus matching MAT tables. `analysis3_context_geometry_matrices` contains every saved 4-by-4 subject matrix and Fisher-z across-subject average in tidy cell form. `exploratory_RSA_all_results.mat` collects all tables. `all_analyses_summary` is the compact subject-by-ROI-by-metric summary.\n');
fprintf(fid, '- `subject_results/` contains one MAT bundle per subject. `rdm_upper_triangle_index.csv` and `context_pairs.csv` make every vector/pair mapping explicit.\n');
n_figures = 54 * double(make_plots);
fprintf(fid, '- `figures/analysis1` through `figures/analysis6` contain the requested descriptive plots (plots generated: %d).\n\n', n_figures);
fprintf(fid, '## QC\n\n');
fprintf(fid, '- Every model records total, valid, and dropped RDM entries plus design rank and status in `tables/QC_model_status`.\n');
fprintf(fid, '- Nonfinite rows are removed jointly within a model and explicitly counted; constant predictors/outcomes and rank-deficient fits return NA estimates with a non-`ok` status.\n');
fprintf(fid, '- This run recorded %d non-`ok` statuses across %d QC rows. Inspect the QC table before interpretation.\n', n_warnings, n_qc);
fprintf(fid, '- Input-level QC recorded %d subjects with at least one source-data warning; see `tables/QC_input_status.csv` for rating missingness and condition-count imbalance.\n', n_input_warnings);
end
