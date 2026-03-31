function motion_regressor = make_motion_regressors (subjidx, nsess_i, nsess_f)
%%%% make motion regressors for OX project data 

SUBJNAMES = {'240711_fMRI_OX_NWU_AS', ...
        '240723_fMRI_OX_NWU_LS', ...
        '240814_fMRI_OX_NWU_JN', ...
        '240816_fMRI_OX_NWU_RR', ...
        '241018_fMRI_OX_NWU_BN', ...
        '250117_fMRI_OX_NWU_VS'}; 
%% 
subjname = SUBJNAMES{subjidx}; 
% wkdir =     '/Volumes/ExtremeSSD/OX_DATA/MRI'; 
wkdir =     '/Users/qhyang/Desktop/OX_DATA/MRI'; 
datapath = fullfile(wkdir, ['subj_', num2str(subjidx)], 'nifti'); % sn:subject's name

% rcntr = 0; % run counter
motion_regressor = [];

for sess = nsess_i:nsess_f
    % sess_2_run_1 is counted as 5 if sess_1 had 4 runs
    path_ = fullfile(datapath, 'func/');
    n = dir(fullfile(path_, [subjname, '_', num2str(sess), '_*.nii']));
    nfiles = length(n);

    for i=1:nfiles
%         rcntr = rcntr+1;
        thisfile = dir(fullfile(path_, ['rp_', subjname, '_', num2str(sess), '_Run', num2str(i), '_*.txt'])); % read motion parameter files 
        fname = fullfile(path_, thisfile.name); 
        mp = load(fname ); 
        current_reg = [mp, [zeros(1,6); diff(mp)], mp.^2, [zeros(1,6); diff(mp).^2]];
        motion_regressor = [motion_regressor; current_reg]; 
    end 

end 
end 
