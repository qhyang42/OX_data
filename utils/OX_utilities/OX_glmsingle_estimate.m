function results = OX_glmsingle_estimate(subject_id, event_type, event_duration)
%OX_GLMSINGLE_ESTIMATE Estimate category-wise single-trial responses.
%   RESULTS = OX_GLMSINGLE_ESTIMATE(SUBJECT_ID, EVENT_TYPE) runs GLMsingle
%   for cue, countdown, or sniff events with zero-duration events.
%
%   RESULTS = OX_GLMSINGLE_ESTIMATE(SUBJECT_ID, EVENT_TYPE, EVENT_DURATION)
%   uses EVENT_DURATION in seconds.

if nargin < 3 || isempty(event_duration)
    event_duration = 0;
end

% info = setup_ox;
% project_root = info.projectRoot;
project_root = '/Users/qhyang/Desktop/OX_DATA'; 
TR = 0.76;

subject = sprintf('subj_%d', subject_id);
mri_data_path = fullfile(project_root, 'MRI', subject, 'nifti');

events = load(fullfile(project_root, 'labchart', 'extracted_events', ...
    sprintf('subj%d_events_bm.mat', subject_id)));

event_type = lower(string(event_type));
switch event_type
    case "cue"
        event_onsets = events.cue_onsets_vec;
    case "countdown"
        pulse_events = load(fullfile(project_root, 'labchart', ...
            'extracted_events', sprintf('subj_%d_events.mat', subject_id)), ...
            'event_onsets_vec');
        event_onsets = pulse_events.event_onsets_vec - 3;
    case "sniff"
        event_onsets = events.event_onsets_vec;
    otherwise
        error('event_type must be cue, countdown, or sniff.');
end
event_TRs = ceil(event_onsets/TR);

[~, category] = OX_get_odor(subject);
design = zeros(sum(events.nframes), 4);
design(event_TRs(strcmp(category, 'PERSON')), 1) = 1;
design(event_TRs(strcmp(category, 'FOOD')), 2) = 1;
design(event_TRs(strcmp(category, 'LOCATION')), 3) = 1;
design(event_TRs(strcmp(category, 'CONTROL')), 4) = 1;

runs = OX_discover_functional_runs(subject_id, ...
    'ProjectRoot', project_root, 'ExpectedRuns', numel(events.nframes));
functional_files = {runs.functional_file}';

mask_header = spm_vol(fullfile(mri_data_path, 'coreg', ...
    'gm_mask_thr05_func.nii'));
mask = logical(spm_read_vols(mask_header));

volume_headers = spm_vol(functional_files);
[~, XYZmm] = spm_read_vols(volume_headers{1});
XYZvx = round(volume_headers{1}(1).mat \ ...
    [XYZmm; ones(1, size(XYZmm, 2))]);
all_volume_headers = cell2mat(volume_headers);
masked_voxels = XYZvx(:, mask(:));
voxel_activity = single(spm_get_data(all_volume_headers, masked_voxels)');
voxel_activity(any(isnan(voxel_activity), 2), :) = [];

number_of_runs = numel(runs);
run_designs = cell(1, number_of_runs);
run_data = cell(1, number_of_runs);
run_confounds = cell(1, number_of_runs);
confound_dir = fullfile(mri_data_path, 'glmsingle_confounds');

for run_index = 1:number_of_runs
    first_frame = sum(events.nframes(1:run_index - 1)) + 1;
    last_frame = sum(events.nframes(1:run_index));
    frame_indices = first_frame:last_frame;

    run_designs{run_index} = design(frame_indices, :);
    run_data{run_index} = voxel_activity(:, frame_indices);

    confound_file = fullfile(confound_dir, sprintf( ...
        '%s_session%02d_run%02d_confounds.mat', subject, ...
        runs(run_index).session, runs(run_index).run));
    confounds = load(confound_file, 'confounds');
    run_confounds{run_index} = confounds.confounds;
end

options.extraregressors = run_confounds;
output_dir = fullfile(mri_data_path, sprintf( ...
    '%s_single_trial_by_category_physio', char(event_type)));
results = GLMestimatesingletrial(run_designs, run_data, event_duration, ...
    TR, output_dir, options);
end
