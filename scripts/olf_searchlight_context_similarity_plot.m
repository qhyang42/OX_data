%% Olfactory context-similarity searchlight: flattened ROI topology plot
% This exploratory script loads one participant's searchlight result and
% draws a deliberately non-anatomical 2-D ROI topology map.
%
% What the plot preserves:
%   * each ROI's superior-inferior-collapsed boundary;
%   * separate left- and right-hemisphere topology;
%   * anterior-to-posterior ROI ordering;
%   * each hemisphere's native projected AP rows without repetition;
%   * the effect value at significant searchlight centers.
%
% What the plot does NOT preserve:
%   * superior-inferior position;
%   * anatomical anterior-posterior distances between ROIs;
%   * an accurate location in a brain.
%
% Run the sections in order. The settings in Section 1 are intended to be
% changed interactively.

%% 1. Editable plotting settings
subject_id = 6;                 % completed participants: 2:6
context_name = "PERSON";       % "PERSON", "FOOD", or "LOCATION"

% Which statistic supplies marker colors. The significance mask is always
% the joint max-statistic FWE p<=.05 result unless ShowOnlySignificant=false.
% Available result-map suffixes:
%   "delta_z", "null_z", "p_uncorrected", or "p_fwe_maxstat"
map_suffix = "delta_z";
show_only_significant = false;   % false shows every valid searchlight center

% BoundaryMode controls the solid ROI outline:
%   "anatomical" = complete input ROI mask (recommended)
%   "tested"     = ROI intersected with the final functional restriction
boundary_mode = "anatomical";

% More than one center can land at the same left-right/AP location after z
% is removed. Choose how those collisions are summarized:
%   "max_abs" = retain the value with the largest absolute magnitude
%   "mean"    = average the values
projection_reducer = "max_abs";

% Leave empty to infer limits from displayed values. Set this manually to
% use identical scales across figures, for example [0, 0.5].
color_limits = [];

save_figure = false;
output_resolution = 220;

%%% 2. Configure the project and resolve participant paths
% Keep this path explicit. This is safer than deriving it from
% mfilename('fullpath') when individual sections are run interactively.
project_root = '~/Desktop/OX_DATA/';

% SPM does not expand "~", so turn the explicit interactive path into an
% absolute path before constructing any NIfTI filename.
if startsWith(project_root, ['~' filesep])
    project_root = fullfile(getenv('HOME'), project_root(3:end));
end
assert(isfolder(project_root), 'Project root does not exist: %s', project_root);

% Configure SPM and the other project dependencies only when SPM is not
% already available in the current interactive MATLAB session.
if exist('spm_vol', 'file') ~= 2
    addpath(fullfile(project_root, 'environment'));
    setup_ox;
end

context_name = upper(string(context_name));
assert(ismember(context_name, ["PERSON", "FOOD", "LOCATION"]), ...
    'context_name must be PERSON, FOOD, or LOCATION.');
assert(ismember(string(map_suffix), ["delta_z", "null_z", ...
    "p_uncorrected", "p_fwe_maxstat"]), ...
    ['map_suffix must be delta_z, null_z, p_uncorrected, ' ...
     'or p_fwe_maxstat.']);
assert(ismember(string(boundary_mode), ["anatomical", "tested"]), ...
    'boundary_mode must be anatomical or tested.');
assert(ismember(string(projection_reducer), ["max_abs", "mean"]), ...
    'projection_reducer must be max_abs or mean.');

subject_name = sprintf('subj_%d', subject_id);
nifti_dir = fullfile(project_root, 'MRI', subject_name, 'nifti');
result_dir = fullfile(nifti_dir, 'olf_context_similarity_searchlight');
prefix = sprintf('olf_context_similarity_subj%d', subject_id);
context_stem = lower(char(context_name));

effect_file = fullfile(result_dir, sprintf('%s_%s_%s.nii', ...
    prefix, context_stem, map_suffix));
significance_file = fullfile(result_dir, sprintf( ...
    '%s_%s_sig_fwe_p05.nii', prefix, context_stem));
valid_center_file = fullfile(result_dir, ...
    sprintf('%s_valid_center_mask.nii', prefix));
final_mask_file = fullfile(result_dir, ...
    sprintf('%s_final_restriction_mask.nii', prefix));

required_files = {effect_file, significance_file, valid_center_file, ...
    final_mask_file};
for file_idx = 1:numel(required_files)
    assert(isfile(required_files{file_idx}), ...
        'Missing searchlight output: %s', required_files{file_idx});
end

%%% 3. Load the result, significance, and analysis-support maps
effect_header = spm_vol(effect_file);
effect_volume = spm_read_vols(effect_header);
significance_header = spm_vol(significance_file);
significant_centers = spm_read_vols(significance_header) > 0;
valid_header = spm_vol(valid_center_file);
valid_centers = spm_read_vols(valid_header) > 0;
final_header = spm_vol(final_mask_file);
final_mask = spm_read_vols(final_header) > 0;

assert_same_geometry(significance_header, effect_header, ...
    significance_file, effect_file);
assert_same_geometry(valid_header, effect_header, ...
    valid_center_file, effect_file);
assert_same_geometry(final_header, effect_header, ...
    final_mask_file, effect_file);

if show_only_significant
    displayed_centers = significant_centers & valid_centers & ...
        isfinite(effect_volume);
    center_description = 'joint max-stat FWE p<=.05 centers';
else
    displayed_centers = valid_centers & isfinite(effect_volume);
    center_description = 'all valid centers';
end

fprintf('\nTopology plot: %s | %s | %s\n', ...
    subject_name, context_name, center_description);
fprintf('Displayed centers before ROI assignment: %d\n', ...
    nnz(displayed_centers));

%%% 4. Load the seven ROI masks used by the searchlight analysis
roi_names = ["TU"; "AON"; "PirF"; "PirT"; "olfAMG"; "EC"; "HIPP"];
roi_groups = ["secondary"; "primary"; "primary"; "primary"; ...
              "primary"; "primary"; "primary"];
roi_basenames = ["olf_TU_bilateral_func_thr02.nii"; ...
                 "AON_bilateral_func_thr02.nii"; ...
                 "PirF_bilateral_func_thr02.nii"; ...
                 "PirT_bilateral_func_thr02.nii"; ...
                 "olfAMG_bilateral_func_thr02.nii"; ...
                 "EC_bilateral_func_thr02.nii"; ...
                 "HIPP_bilateral_func_thr02.nii"];

n_rois = numel(roi_names);
roi_masks = false([effect_header.dim, n_rois]);
roi_files = strings(n_rois, 1);
roi_ap_centroid_mm = nan(n_rois, 1);
roi_si_centroid_mm = nan(n_rois, 1);
roi_reference_k = nan(n_rois, 1);

for roi_idx = 1:n_rois
    roi_files(roi_idx) = fullfile(nifti_dir, 'coreg', 'roi_decoding', ...
        roi_groups(roi_idx), roi_basenames(roi_idx));
    assert(isfile(roi_files(roi_idx)), ...
        'Missing ROI mask: %s', roi_files(roi_idx));
    roi_header = spm_vol(char(roi_files(roi_idx)));
    assert_same_geometry(roi_header, effect_header, ...
        roi_files(roi_idx), effect_file);
    roi_mask = spm_read_vols(roi_header) > 0;
    roi_masks(:, :, :, roi_idx) = roi_mask;

    [voxel_i, voxel_j, voxel_k] = ind2sub(effect_header.dim, find(roi_mask));
    world = effect_header.mat * [voxel_i'; voxel_j'; voxel_k'; ...
        ones(1, numel(voxel_i))];
    roi_ap_centroid_mm(roi_idx) = mean(world(2, :));
    roi_si_centroid_mm(roi_idx) = mean(world(3, :));
    roi_reference_k(roi_idx) = round(mean(voxel_k));
end

roi_table = table(roi_names, roi_groups, roi_files, ...
    roi_ap_centroid_mm, roi_si_centroid_mm, ...
    'VariableNames', {'ROI', 'MaskGroup', 'MaskFile', ...
    'APCentroidMM', 'SICentroidMM'});
roi_table = sortrows(roi_table, 'APCentroidMM', 'descend');
disp(roi_table(:, {'ROI', 'APCentroidMM', 'SICentroidMM'}));

%%% 5. Project each ROI and assemble one m-by-n grid per hemisphere
% World x<0 defines the left hemisphere and world x>0 the right. The
% projected ROI blocks are stacked anterior-to-posterior with no blank rows
% between blocks. Each output contains a numeric map grid, a logical data
% cell grid, and one logical contour grid per ROI.
hemisphere_names = ["Left", "Right"];
n_hemispheres = numel(hemisphere_names);
world_x_volume = make_world_x_volume(effect_header);
hemisphere_masks = cat(4, world_x_volume < 0, world_x_volume > 0);

projected_anatomical = cell(n_rois, n_hemispheres);
projected_tested = cell(n_rois, n_hemispheres);
projected_effect = cell(n_rois, n_hemispheres);

for roi_idx = 1:n_rois
    for hemisphere_idx = 1:n_hemispheres
        hemisphere_roi = roi_masks(:, :, :, roi_idx) & ...
            hemisphere_masks(:, :, :, hemisphere_idx);
        projected_anatomical{roi_idx, hemisphere_idx} = ...
            any(hemisphere_roi, 3);
        projected_tested{roi_idx, hemisphere_idx} = ...
            any(hemisphere_roi & final_mask, 3);

        effect_in_roi = effect_volume;
        effect_in_roi(~(hemisphere_roi & displayed_centers)) = NaN;
        projected_effect{roi_idx, hemisphere_idx} = ...
            collapse_effect_along_z(effect_in_roi, projection_reducer);
    end
end

[~, anterior_to_posterior_order] = sort(roi_ap_centroid_mm, 'descend');
[world_x_grid, world_ap_grid] = make_world_coordinate_grids( ...
    effect_header, round((effect_header.dim(3) + 1) / 2));

if boundary_mode == "tested"
    projected_boundaries = projected_tested;
else
    projected_boundaries = projected_anatomical;
end

topology_grids = cell(n_hemispheres, 1);
for hemisphere_idx = 1:n_hemispheres
    topology_grids{hemisphere_idx} = build_topology_grid( ...
        projected_effect(:, hemisphere_idx), ...
        projected_anatomical(:, hemisphere_idx), ...
        projected_boundaries(:, hemisphere_idx), ...
        anterior_to_posterior_order, ...
        world_x_grid, world_ap_grid);
end

% Convenient workspace aliases for interactive inspection/modification.
left_map_grid = topology_grids{1}.map_grid;
left_data_cell_grid = topology_grids{1}.data_cell_grid;
right_map_grid = topology_grids{2}.map_grid;
right_data_cell_grid = topology_grids{2}.data_cell_grid;
left_native_roi_heights = topology_grids{1}.roi_native_heights;
right_native_roi_heights = topology_grids{2}.roi_native_heights;

%% 6. Display each composite grid directly with imagesc, then contour
figure_handle = figure('Color', 'w', 'Position', [80, 80, 1200, 820]);
layout_handle = tiledlayout(figure_handle, 1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
axes_handles = gobjects(n_hemispheres, 1);
roi_colors = lines(n_rois);

all_displayed_values = effect_volume(displayed_centers & ...
    any(roi_masks, 4));
all_displayed_values = all_displayed_values(isfinite(all_displayed_values));
if isempty(color_limits) && ~isempty(all_displayed_values)
    shared_color_limits = [min(all_displayed_values), ...
        max(all_displayed_values)];
    if shared_color_limits(1) == shared_color_limits(2)
        padding = max(0.01, abs(shared_color_limits(1)) * 0.05);
        shared_color_limits = shared_color_limits + [-padding, padding];
    end
elseif ~isempty(color_limits)
    assert(isnumeric(color_limits) && numel(color_limits) == 2 && ...
        color_limits(1) < color_limits(2), ...
        'color_limits must be empty or [minimum maximum].');
    shared_color_limits = color_limits;
else
    shared_color_limits = [];
end

for hemisphere_idx = 1:n_hemispheres
    topology = topology_grids{hemisphere_idx};
    axes_handle = nexttile(layout_handle, hemisphere_idx);
    axes_handles(hemisphere_idx) = axes_handle;

    % AlphaData is the simplest transparent/white non-data-cell method.
    % map_grid is numeric everywhere; only true data_cell_grid cells show.
    imagesc(axes_handle, topology.map_grid, ...
        'AlphaData', double(topology.data_cell_grid));
    axes_handle.Color = 'w';
    hold(axes_handle, 'on');
    box(axes_handle, 'off');

    % Draw contours only after the one imagesc call, keeping outlines on top.
    for roi_idx = 1:n_rois
        roi_contour = topology.contour_grids(:, :, roi_idx);
        if any(roi_contour, 'all')
            contour(axes_handle, double(roi_contour), [0.5, 0.5], ...
                'Color', roi_colors(roi_idx, :), 'LineWidth', 1.8);
        end
    end

    yticks(axes_handle, topology.roi_row_centers);
    yticklabels(axes_handle, roi_names(anterior_to_posterior_order));
    xticks(axes_handle, topology.x_tick_positions);
    xticklabels(axes_handle, compose('%.0f', topology.x_tick_world_mm));
    xlim(axes_handle, [1, size(topology.map_grid, 2)]);
    ylim(axes_handle, [1, size(topology.map_grid, 1)]);
    axes_handle.YDir = 'reverse';
    axes_handle.TickDir = 'out';
    axes_handle.FontSize = 11;
    if hemisphere_idx == 2
        axes_handle.YAxisLocation = 'right';
    end
    xlabel(axes_handle, 'World x coordinate (mm)');
    title(axes_handle, sprintf('%s hemisphere', ...
        hemisphere_names(hemisphere_idx)), 'FontSize', 13);
    colormap(axes_handle, turbo(256));
    if ~isempty(shared_color_limits)
        clim(axes_handle, shared_color_limits);
    else
        text(axes_handle, 0.5, 0.02, ...
            'No centers survive the selected display threshold.', ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'Color', [0.35, 0.35, 0.35], 'FontAngle', 'italic');
    end
    text(axes_handle, 0.98, 0.98, 'Anterior', 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'FontWeight', 'bold');
    text(axes_handle, 0.98, 0.02, 'Posterior', 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'FontWeight', 'bold');
end

title(layout_handle, sprintf('%s | %s context similarity | %s', ...
    strrep(subject_name, '_', ' '), context_name, center_description), ...
    'Interpreter', 'none', 'FontSize', 14);
subtitle(layout_handle, [ ...
    'Independent native-row topology grids; no AP-row repetition; ' ...
    'transparent cells contain no data; contours are overlaid last'], ...
    'FontSize', 10, 'FontWeight', 'normal');

if ~isempty(shared_color_limits)
    colorbar_handle = colorbar(axes_handles(2));
    colorbar_handle.Layout.Tile = 'east';
    colorbar_handle.Label.String = strrep(map_suffix, '_', ' ');
end

%% 7. Optional export
output_file = fullfile(result_dir, sprintf( ...
    '%s_%s_%s_topology.png', prefix, context_stem, map_suffix));
if save_figure
    exportgraphics(figure_handle, output_file, ...
        'Resolution', output_resolution);
    fprintf('Saved topology figure: %s\n', output_file);
else
    fprintf('Figure not saved. Set save_figure=true and rerun Section 7.\n');
end

%% Local helper functions
function [world_x_grid, world_ap_grid] = make_world_coordinate_grids( ...
        header, reference_k)
% Project an i-j grid at one representative k coordinate through the NIfTI
% affine. World x/y prevents accidental radiological-axis flips.
[voxel_i, voxel_j] = ndgrid(1:header.dim(1), 1:header.dim(2));
voxel_k = repmat(reference_k, size(voxel_i));
world = header.mat * [voxel_i(:)'; voxel_j(:)'; voxel_k(:)'; ...
    ones(1, numel(voxel_i))];
world_x_grid = reshape(world(1, :), size(voxel_i));
world_ap_grid = reshape(world(2, :), size(voxel_i));
end

function world_x_volume = make_world_x_volume(header)
% Evaluate world x at every voxel so hemisphere assignment remains valid
% even if the functional image affine is oblique.
[voxel_i, voxel_j, voxel_k] = ndgrid(1:header.dim(1), ...
    1:header.dim(2), 1:header.dim(3));
world_x = header.mat(1, 1) * voxel_i + ...
    header.mat(1, 2) * voxel_j + ...
    header.mat(1, 3) * voxel_k + header.mat(1, 4);
world_x_volume = reshape(world_x, header.dim);
end

function topology = build_topology_grid(effect_cells, anatomical_cells, ...
        boundary_cells, roi_order, ...
        world_x_grid, world_ap_grid)
% Put every ROI into one numeric m-by-n grid. One transparent padding cell
% around the grid lets contour close masks that touch an outside edge.
n_rois = numel(roi_order);
oriented_effect = cell(n_rois, 1);
oriented_anatomical = cell(n_rois, 1);
oriented_boundary = cell(n_rois, 1);
native_row_ranges = cell(n_rois, 1);
roi_native_heights = zeros(n_rois, 1);
anatomical_union = [];

for roi_idx = 1:n_rois
    [x_world, ~, oriented_effect{roi_idx}] = ...
        orient_projected_matrix(world_x_grid, world_ap_grid, ...
        effect_cells{roi_idx});
    [~, ~, oriented_anatomical{roi_idx}] = ...
        orient_projected_matrix(world_x_grid, world_ap_grid, ...
        anatomical_cells{roi_idx});
    [~, ~, oriented_boundary{roi_idx}] = ...
        orient_projected_matrix(world_x_grid, world_ap_grid, ...
        boundary_cells{roi_idx});
    if isempty(anatomical_union)
        anatomical_union = oriented_anatomical{roi_idx};
    else
        anatomical_union = anatomical_union | oriented_anatomical{roi_idx};
    end
    occupied_rows = find(any(oriented_anatomical{roi_idx}, 2));
    assert(~isempty(occupied_rows), ...
        'A projected ROI has no anatomical rows in this hemisphere.');
    native_row_ranges{roi_idx} = occupied_rows(1):occupied_rows(end);
    roi_native_heights(roi_idx) = numel(native_row_ranges{roi_idx});
end

occupied_columns = find(any(anatomical_union, 1));
assert(~isempty(occupied_columns), 'No projected ROI columns were found.');
column_range = occupied_columns(1):occupied_columns(end);
n_columns = numel(column_range);
n_rows = sum(roi_native_heights(roi_order));

% Padding is deliberately part of the m-by-n grid and remains non-data.
map_grid = zeros(n_rows + 2, n_columns + 2);
data_cell_grid = false(size(map_grid));
contour_grids = false(n_rows + 2, n_columns + 2, n_rois);
roi_row_centers = nan(n_rois, 1);
row_cursor = 2;

for order_idx = 1:n_rois
    roi_idx = roi_order(order_idx);
    source_rows = native_row_ranges{roi_idx};
    target_rows = row_cursor:(row_cursor + numel(source_rows) - 1);
    target_columns = 2:(n_columns + 1);

    effect_block = oriented_effect{roi_idx}( ...
        source_rows, column_range);
    block_has_data = isfinite(effect_block);
    effect_block(~block_has_data) = 0;
    map_grid(target_rows, target_columns) = effect_block;
    data_cell_grid(target_rows, target_columns) = block_has_data;
    contour_grids(target_rows, target_columns, roi_idx) = ...
        oriented_boundary{roi_idx}(source_rows, column_range);
    roi_row_centers(order_idx) = mean(target_rows);
    row_cursor = target_rows(end) + 1;
end

inner_x_world = x_world(column_range);
n_x_ticks = min(4, n_columns);
inner_tick_index = unique(round(linspace(1, n_columns, n_x_ticks)));
topology.map_grid = map_grid;
topology.data_cell_grid = data_cell_grid;
topology.contour_grids = contour_grids;
topology.roi_row_centers = roi_row_centers;
topology.roi_native_heights = roi_native_heights;
topology.x_tick_positions = inner_tick_index + 1;
topology.x_tick_world_mm = inner_x_world(inner_tick_index);
end

function [x_world, ap_world, image_matrix] = orient_projected_matrix( ...
        world_x_grid, world_ap_grid, projected_matrix)
% Convert i-by-j voxel arrays into image rows-by-columns. Columns increase
% in world x; rows proceed from anterior to posterior for direct imagesc.
x_world = mean(world_x_grid, 2);
ap_world = mean(world_ap_grid, 1)';
image_matrix = projected_matrix';

if x_world(1) > x_world(end)
    x_world = flip(x_world);
    image_matrix = flip(image_matrix, 2);
end
if ap_world(1) < ap_world(end)
    ap_world = flip(ap_world);
    image_matrix = flip(image_matrix, 1);
end
end

function assert_same_geometry(image_header, reference_header, ...
        image_file, reference_file)
assert(isequal(image_header.dim, reference_header.dim), ...
    'Dimension mismatch between %s and %s.', image_file, reference_file);
assert(max(abs(image_header.mat(:) - reference_header.mat(:))) < 1e-4, ...
    'Affine mismatch between %s and %s.', image_file, reference_file);
end

function collapsed = collapse_effect_along_z(effect_volume, reducer)
% Reduce centers that collide after removing the superior-inferior axis.
% max_abs keeps the strongest signed value; mean gives an average summary.
if reducer == "mean"
    collapsed = mean(effect_volume, 3, 'omitmissing');
    collapsed(all(~isfinite(effect_volume), 3)) = NaN;
    return;
end

absolute_effect = abs(effect_volume);
absolute_effect(~isfinite(absolute_effect)) = -Inf;
[~, strongest_k] = max(absolute_effect, [], 3);
[voxel_i, voxel_j] = ndgrid(1:size(effect_volume, 1), ...
    1:size(effect_volume, 2));
linear_index = sub2ind(size(effect_volume), voxel_i, voxel_j, strongest_k);
collapsed = effect_volume(linear_index);
collapsed(all(~isfinite(effect_volume), 3)) = NaN;
end
