%% 
% realignment, coreg and smoothing for OX project 
% preproc is done per subject 
% realign wb image -- realign func -- coreg between anat and wb --
% smoothing of func (vivek used fwhm = 2 here. I have 2mm iso voxel size. I
% should probably increase the window size?) 

%% BEFORE PREPROCESSING! 
% make sure files names are continuous and no repeat. 
%% 
SUBJNAMES = {'240711_fMRI_OX_NWU_AS', ...
        '240723_fMRI_OX_NWU_LS', ...
        '240814_fMRI_OX_NWU_JN', ...
        '240816_fMRI_OX_NWU_RR', ...
        '241018_fMRI_OX_NWU_BN', ...
        '250117_fMRI_OX_NWU_VS'}; 

session_count = [4, 13, 10, 15, 16, 10]; % number of sessions for each subj so far. EDIT as needed. 
subjidx = 5; % enter subj name here 

%% 
subjname = SUBJNAMES{subjidx}; 
wkdir = '/Users/qhyang/Desktop/OX_DATA/MRI'; 
datapath = fullfile(wkdir, ['subj_', num2str(subjidx)], 'nifti'); % sn:subject's name

matlabbatch = [];

%% REALIGNMENT OF WHOLE BRAIN
filename_wb = wb_realign_list(subjname, datapath, 1, session_count(subjidx)); 

%% 
matlabbatch{1}.spm.spatial.realign.estwrite.data = filename_wb;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

%% Functional scan realignment

filename_func = func_realign_list(subjname, datapath, 1, session_count(subjidx)); 

%%%% exception -- invalid run files are now located in trash_run folder in
%%%% each subj's func folder. 

%% 
% Corregistration needs files from different runs in one big array.
fimages = cat(1,filename_func{:});
% fimages = filename_func;

matlabbatch{2}.spm.spatial.realign.estwrite.data = filename_func;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 0; % register to 1st 
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [0 1];
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

%% Anatomical realignment and avraging 
%%% align and average across all anatomical images 
filename_anat = anat_realign_list(subjname, datapath, 1, session_count(subjidx)); 

%% 
matlabbatch{3}.spm.spatial.realign.estwrite.data = filename_anat;
matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{3}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{3}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{3}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{3}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{3}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

% this should give us a mean anat template 

%% 
spm_jobman('run', matlabbatch(1:3));    % SPM jobman will run the batch

%% COREGISTRATION- func to anat 
% Anatomical and whole brain coreg

% define path for anatomical image
anatpath = fullfile(datapath, 'anat');
[~, n] = fileparts(filename_anat{1}); 
anatfile = ['mean', n, '.nii'];
fname = fullfile(anatpath, anatfile); % mean anatomical image
wbpath = fullfile(datapath, 'wb');
[~, wbn] = fileparts(filename_wb{1});
wbfile = ['mean', wbn, '.nii'];
wbfname = fullfile(wbpath, wbfile); % grab the whole brain file 

matlabbatch{4}.spm.spatial.coreg.estimate.ref = {fname};
matlabbatch{4}.spm.spatial.coreg.estimate.source = {wbfname};
matlabbatch{4}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.tol = ...
    [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

% Whole brain and functional scan coreg
% fname = fullfile(wbpath, sprintf('mean%s', wbfile));
funcpath = fullfile(datapath, 'func'); 
[~, funcfile] = fileparts(fimages{1}); 
matlabbatch{5}.spm.spatial.coreg.estimate.ref = {wbfname};
fname = fullfile(funcpath, sprintf('mean%s.nii', funcfile));
matlabbatch{5}.spm.spatial.coreg.estimate.source = {fname};
matlabbatch{5}.spm.spatial.coreg.estimate.other = fimages;
matlabbatch{5}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{5}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{5}.spm.spatial.coreg.estimate.eoptions.tol =...
    [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{5}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];


%% COREGISTRATION- anat to func
% Anatomical and whole brain coreg

% define path for anatomical image
anatpath = fullfile(datapath, 'anat');
[~, n] = fileparts(filename_anat{1}); 
anatfile = ['mean', n, '.nii'];
fname = fullfile(anatpath, anatfile); % mean anatomical image
wbpath = fullfile(datapath, 'wb');
[~, wbn] = fileparts(filename_wb{1});
wbfile = ['mean', wbn, '.nii'];
wbfname = fullfile(wbpath, wbfile); % grab the whol e brain file 

matlabbatch{4}.spm.spatial.coreg.estimate.ref = {wbfname};
matlabbatch{4}.spm.spatial.coreg.estimate.source = {fname};
matlabbatch{4}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.tol = ...
    [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

% Whole brain and functional scan coreg
% fname = fullfile(wbpath, sprintf('mean%s', wbfile));
funcpath = fullfile(datapath, 'func'); 
[~, funcfile] = fileparts(fimages{1}); 
matlabbatch{5}.spm.spatial.coreg.estimate.source = {wbfname};
fname = fullfile(funcpath, sprintf('mean%s.nii', funcfile));
matlabbatch{5}.spm.spatial.coreg.estimate.ref = {fname};
% matlabbatch{5}.spm.spatial.coreg.estimate.other = fimages;
matlabbatch{5}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{5}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{5}.spm.spatial.coreg.estimate.eoptions.tol =...
    [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{5}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];


%% Normalize
fname = fullfile(anatpath, anatfile);
matlabbatch{6}.spm.spatial.normalise.estwrite.subj.vol = {fname};
matlabbatch{6}.spm.spatial.normalise.estwrite.subj.resample =  {fname};
matlabbatch{6}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{6}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{6}.spm.spatial.normalise.estwrite.eoptions.tpm = {fullfile(spm('dir'),'tpm','TPM.nii')};
matlabbatch{6}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
matlabbatch{6}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{6}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch{6}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch{6}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70; 78 76 85];
matlabbatch{6}.spm.spatial.normalise.estwrite.woptions.vox = [1 1 1];
matlabbatch{6}.spm.spatial.normalise.estwrite.woptions.interp = 4;

%     % Normalise write % we don't normalize all the functional images at
%     this point. all analysis are done in the native space. 
% fname = fullfile(anatpath, sprintf('y_%s', anatfile));
% matlabbatch{7}.spm.spatial.normalise.write.subj.def = {fname};
% matlabbatch{7}.spm.spatial.normalise.write.subj.resample = fimages;
% matlabbatch{7}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
% matlabbatch{7}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
% matlabbatch{7}.spm.spatial.normalise.write.woptions.interp = 4;

%% REALIGN-RESLICE FUNCTIONAL IMAGES
matlabbatch{8}.spm.spatial.realign.write.data = fimages;
matlabbatch{8}.spm.spatial.realign.write.roptions.which = [2 1];
matlabbatch{8}.spm.spatial.realign.write.roptions.interp = 4;
matlabbatch{8}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
matlabbatch{8}.spm.spatial.realign.write.roptions.mask = 1;
matlabbatch{8}.spm.spatial.realign.write.roptions.prefix = 'r';

%%%% realign all functional images 
% Assuming 'filename' contains all the volumes for each run
% and 'fimages' contains the first volume of each run, realigned

%%%% NOTE this is trash 

% for r = 1:length(fimages)
%     first_frame = fimages{r};  % Realigned first frame of run r
% 
%     % Get all the volumes for this run (from the original filename variable)
%     run_volumes = filename{r};
% 
%     % Create a batch to realign all volumes to the first frame
%     matlabbatch{r}.spm.spatial.realign.estwrite.data = {run_volumes};  % All volumes for run r
%     matlabbatch{r}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
%     matlabbatch{r}.spm.spatial.realign.estwrite.eoptions.sep = 4;
%     matlabbatch{r}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
%     matlabbatch{r}.spm.spatial.realign.estwrite.eoptions.rtm = 0; % Set to 0 to realign to the first volume
%     matlabbatch{r}.spm.spatial.realign.estwrite.eoptions.interp = 2;
%     matlabbatch{r}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
%     matlabbatch{r}.spm.spatial.realign.estwrite.eoptions.weight = '';
%     matlabbatch{r}.spm.spatial.realign.estwrite.roptions.which = [2 1]; % Save realigned images
%     matlabbatch{r}.spm.spatial.realign.estwrite.roptions.interp = 4;
%     matlabbatch{r}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
%     matlabbatch{r}.spm.spatial.realign.estwrite.roptions.mask = 1;
%     matlabbatch{r}.spm.spatial.realign.estwrite.roptions.prefix = 'r2';  % Prefix for realigned images
% end



%% Spatial smoothing
filename_sm = sm_list(subjname, datapath, 1, session_count(subjidx)); 
%% 
matlabbatch{9}.spm.spatial.smooth.data = filename_sm;
matlabbatch{9}.spm.spatial.smooth.fwhm = [2.5 2.5 2.5]; % orig voxel size is 2mm iso
matlabbatch{9}.spm.spatial.smooth.dtype = 0;
matlabbatch{9}.spm.spatial.smooth.im = 0;
matlabbatch{9}.spm.spatial.smooth.prefix = 's';


% %% test smoothing on machine 
% matlabbatch.spm.spatial.smooth.data = {'~/Desktop/OX_DATA/testsm.nii'};
% matlabbatch.spm.spatial.smooth.fwhm = [5 5 5]; % orig voxel size is 2mm iso
% matlabbatch.spm.spatial.smooth.dtype = 0;
% matlabbatch.spm.spatial.smooth.im = 0;
% matlabbatch.spm.spatial.smooth.prefix = 's';

%% 
spm_jobman('run', matlabbatch(8:9));    % SPM jobman will run the batch


%% functions 

%% 
function filename = wb_realign_list(subjname, datapath, nsess_i, nsess_f)
%%% make file list for whole brain realignment 
% nsess = nsess_f - nsess_i +1; 
filename = [];
ss=0;

for sess = nsess_i:nsess_f

    % path for whole brain images.
    path_ = fullfile(datapath, 'wb'); 
    % *.nii will give all the files.
    n = dir(fullfile(path_, [subjname, '_', num2str(sess), '_*.nii']));
    if ~isempty(n)
        ss=ss+1;
        fname = fullfile(path_, n.name);
        % SPM job file has format (----.nii,1)
        filename{ss}{1,1} = sprintf('%s,1', fname);
    end
end
end 

%% 
function filename = func_realign_list(subjname, datapath, nsess_i, nsess_f)

%%%%% NOTE. This is now working for 4D nifti volumes. 

% rnctr is run counter
rcntr = 0;
filename = [];

for sess = nsess_i:nsess_f

    path_ = fullfile(datapath, 'func/');
    n = dir(fullfile(path_, [subjname, '_', num2str(sess), '_*.nii']));
    nfiles = length(n); 

    %     if rcntr==1
    %         funcpath = path_;
    %         funcfile = n(1).name;
    %     end

    for i=1:length(n)
        rcntr = rcntr+1; % sess_2_run_1 is counted as 5 if sess_1 had 4 runs
        thisfile = dir(fullfile(path_, [subjname, '_', num2str(sess), '_Run', num2str(i), '_*.nii']));
        fname = fullfile(path_, thisfile.name);
        % Different runs in different cells, add ",1" for spm
%         filename{rcntr}{i,1} = sprintf('%s,1', fname);
        filename{rcntr, 1}{1,1} = sprintf('%s', fname);
    end

end


end

%% 
function filename = anat_realign_list(subjname, datapath, nsess_i, nsess_f)
%%% make file list for whole brain realignment 

filename = [];
ss=0;

for sess = nsess_i:nsess_f
    
    % path for whole brain images.
    path_ = fullfile(datapath, 'anat'); 
    
    % *.nii will give all the files.
    n = dir(fullfile(path_, [subjname, '_', num2str(sess), '_*.nii']));
    if ~isempty(n)
    ss=ss+1;
    fname = fullfile(path_, n.name);
    % SPM job file has format (----.nii,1)
    filename{ss}{1,1} = sprintf('%s,1', fname);
    end
end
end 
%% 
function filename = sm_list(subjname, datapath, nsess_i, nsess_f)
%%%% make list for functional files that needs smoothing.
%%%% NOTE: SPM can handle 4D nifti now 
ii=0;
filename = [];
for sess = nsess_i:nsess_f
    path_ = fullfile(datapath, 'func');
    n = dir(fullfile(path_, [subjname, '_', num2str(sess), '_*.nii']));

    %             for i=1:length(n)
    %                 ii=ii+1;
    %                 fname = fullfile(path_, sprintf('r%s',n(i).name));
    %                 filename_sm{ii,1} = sprintf('%s,1', fname);
    %             end

    %             for i=1:length(n)
    %                 ii = ii+1;
    %                 thisfile = dir(fullfile(path_, [subjname, '_', num2str(sess), '_Run', num2str(i), '_*.nii']));
    %                 fname = fullfile(path_, ['r', thisfile.name]);
    %                 % Different runs in different cells, add ",1" for spm
    %                 filename{ii, 1} = sprintf('%s,1', fname);
    %             end


    for i=1:length(n)

        ii = ii+1;
        thisfile = dir(fullfile(path_, [subjname, '_', num2str(sess), '_Run', num2str(i), '_*.nii']));
        fname_r = fullfile(path_, ['r', thisfile.name]);

        fname = fullfile(path_, [thisfile.name]);
        filename{ii, 1} = sprintf('%s', fname_r);

        %         nifti_header = niftiinfo(fname);
        %         nvols = nifti_header.ImageSize(4);

        %         filename{ii} = [];
        %         for vol = 1:nvols
        %             % Different runs in different cells, add ",1" for spm
        %             filename{ii}{vol,1} = sprintf('%s,%d', fname_r, vol);
        %         end
        %     end
    end
end
end


% function filename = sm_list(datapath, fimages)
% %%% make cell array for smoothing. All images in a single cell array
% filename = [];
% path_ = fullfile(datapath, 'func');
% for n = 1: length(fimages)
%     thisimg = fimages{n};
%     [~, fname] = fileparts(thisimg);
%     fname_s = fullfile(path_, ['r', fname, '.nii']);
%     filename{n, 1} = fname_s;
% 
% end
% end



%% 
function filename = func_realign_list_all(subjname, datapath, nsess_i, nsess_f)
%%% DO NOT USE. 
% %%% change this. SPM can handle 4D nifti now. 
rcntr = 0; % run counter
filename = [];

for sess = nsess_i:nsess_f
    % sess_2_run_1 is counted as 5 if sess_1 had 4 runs
    path_ = fullfile(datapath, 'func/');
    n = dir(fullfile(path_, [subjname, '_', num2str(sess), '*.nii']));
    nfiles = length(n);

    for i=1:nfiles

        rcntr = rcntr+1;
        thisfile = dir(fullfile(path_, [subjname, '_', num2str(sess), '_Run', num2str(i), '_*.nii']));
        fname = fullfile(path_, thisfile.name);

        nifti_header = niftiinfo(fname);
        nvols = nifti_header.ImageSize(4);

        filename{rcntr} = [];
        for vol = 1:nvols
            % Different runs in different cells, add ",1" for spm
            filename{rcntr}{vol,1} = sprintf('%s,%d', fname, vol);
        end
    end

end
end
%% 
% function filename = all_frame_realign_list(subjname, datapath, nsess_i, nsess_f, r_fimages)
% %%%% get all frames for realignment 
% 
%     % Initialize output cell array
%     filename = {};
% 
%     % Loop over each session
%     for sess = nsess_i:nsess_f
%         % Path to functional data
%         path_ = fullfile(datapath, 'func/');
%         
%         % Get list of all runs (NIfTI files) for this session
%         run_files = dir(fullfile(path_, [subjname, '_', num2str(sess), '_Run*.nii']));
%         nruns = length(run_files);
% 
%         if nruns == 0
%             warning('No runs found for session %d.', sess);
%             continue;
%         end
% 
%         % Sort the run files numerically based on the run number
%         run_names = {run_files.name};  % Extract the filenames
%         run_numbers = regexp(run_names, '_Run(\d+)', 'tokens');  % Extract the run numbers using regular expression
%         run_numbers = cellfun(@(x) str2double(x{1}), run_numbers);  % Convert to numeric
% 
%         % Sort the run_files based on the extracted run numbers
%         [~, sort_idx] = sort(run_numbers);
%         run_files = run_files(sort_idx);
% 
%         % Initialize the cell array for this session
%         session_files = {};
% 
%         % Loop over each run in the session (now correctly sorted)
%         for run = 1:nruns
%             % Get the NIfTI file for this run
%             run_fname = fullfile(path_, run_files(run).name);
% 
%             % Use niftiinfo to quickly retrieve header info
%             nifti_header = niftiinfo(run_fname);
%             nvols = nifti_header.ImageSize(4);  % Number of volumes (4th dimension)
% 
%             % Initialize cell array for the volumes in this run
%             vol_files = cell(nvols, 1);
% 
%             % Add the realigned first frame to the beginning of the volume list
%             vol_files{1} = r_fimages{sess - nsess_i + 1}{run};  % Reference the realigned first frame from fimages
% 
%             % Loop over each subsequent volume in the run, starting from the 2nd volume
%             for vol = 2:nvols
%                 % Store each volume with the correct format
%                 vol_files{vol} = sprintf('%s,%d', run_fname, vol);
%             end
% 
%             % Add this run's volumes to the session's run list
%             session_files{run} = vol_files;
%         end
% 
%         % Add this session's runs to the output filename cell array
%         filename{sess - nsess_i + 1} = session_files;
%     end
% 
% end
% 
% 
% 
% 
% end 