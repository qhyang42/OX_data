%% Putamen/A1 control semantic-context split-half similarity: GM only.
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root, 'environment'));
setup_ox;
subject_ids = 2:6;
mri_root = fullfile(project_root, 'MRI');
output_dir = fullfile(mri_root, 'group', 'context_split_half_similarity_physio_control');
assert(~isfolder(output_dir), 'Output already exists; preserve completed results.');
roi_names = ["control_GM_putamen", "control_A1"];
OX_test_context_split_half_similarity();
% Keep process workers and their files inside the project workspace.
cluster = parcluster('Processes');
cluster.JobStorageLocation = fullfile(project_root, 'tmp', 'control_split_half_pool');
if ~isfolder(cluster.JobStorageLocation), mkdir(cluster.JobStorageLocation); end
if isempty(gcp('nocreate')), parpool(cluster, 4); end
for subject_id = subject_ids
    roi_dir = fullfile(mri_root, sprintf('subj_%d',subject_id), ...
        'nifti','coreg','roi_decoding','control');
    % 'old' invokes the resolver's explicit flat-directory override only;
    % all mask paths point to control/, never to archived old_rois/.
    results = OX_roi_context_split_half_similarity(subject_id, ...
        'MRIRoot',mri_root,'ROISelection','old','ROIDir',roi_dir, ...
        'ROINames',roi_names,'UseFunctionalRestriction',false, ...
        'NumSplits',200,'RandomSeed',1,'NumPermutations',5000, ...
        'PermutationSeed',1001,'UseParallel',true, ...
        'StratifyBySession',true,'MinVoxels',10, ...
        'OutputDir',output_dir,'SaveOutputs',true);
    assert(all(string(results.summary.status)=="ok"));
    assert(all(results.roi_metadata.n_features_used>=10));
    assert(~results.preprocessing.use_functional_restriction);
end
group_results = OX_group_context_split_half_similarity(subject_ids, ...
    'MRIRoot',mri_root,'OutputDir',output_dir,'MakePlots',true,'SaveOutputs',true);
