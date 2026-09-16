function export_rsa_tsnr_inputs()
% Export exact saved RSA features and canonical run mapping for tSNR QC.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'environment')); setup_ox;
out=fullfile(root,'results','ROI_tSNR');
assert(~isfile(fullfile(out,'inputs.json')),'Preserve existing exported inputs.');
mkdir(out);
subjects={};
for s=2:6
    runs=OX_discover_functional_runs(s,'ProjectRoot',root);
    meta=OX_load_trial_metadata(s);
    events=load(fullfile(root,'labchart','extracted_events',sprintf('subj%d_events_bm.mat',s)),'nframes');
    assert(numel(events.nframes)==80);
    for k=1:80
        m=meta(meta.run_id==k,:);
        assert(height(m)==10 && all(m.session_id==runs(k).session) && all(m.run_in_session==runs(k).run));
        runs(k).expected_frames=events.nframes(k);
    end
    gmfile=fullfile(root,'MRI',sprintf('subj_%d',s),'nifti','coreg','gm_mask_thr05_func.nii');
    gm=find(spm_read_vols(spm_vol(gmfile))>0);
    rois={};
    for family={'neural_TU','neural','neural_control'}
        if s==5 && strcmp(family{1},'neural_TU'), continue; end
        source=fullfile(root,'RDMs',family{1},sprintf('subj_%d_neural_RDMs.mat',s));
        N=load(source,'results'); R=N.results;
        assert(isequal(R.trial_metadata,meta));
        for j=1:numel(R.roi_results)
            r=R.roi_results(j); idx=r.model_feature_indices;
            assert(numel(idx)==r.n_voxels && r.n_voxels>=10);
            rois{end+1}=struct('name',r.roi_name,'mask_file',r.mask_file, ...
                'source',source,'n_voxels',r.n_voxels, ...
                'linear_indices_zero_based',gm(idx)-1, ...
                'functional_mask',R.preprocessing.functional_mask); %#ok<AGROW>
        end
    end
    subjects{end+1}=struct('subject',s,'gm_file',gmfile,'runs',runs,'rois',[rois{:}]); %#ok<AGROW>
end
fid=fopen(fullfile(out,'inputs.json'),'w'); assert(fid>=0);
fprintf(fid,'%s',jsonencode([subjects{:}],PrettyPrint=true)); fclose(fid);
end
