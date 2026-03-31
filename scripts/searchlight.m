%% searchlight MVPA using GLMsingle estimated single trial beta 

%% 
SUBJNAMES = {'240711_fMRI_OX_NWU_AS', ...
        '240723_fMRI_OX_NWU_LS', ...
        '240814_fMRI_OX_NWU_JN', ...
        '240816_fMRI_OX_NWU_RR', ...
        '241018_fMRI_OX_NWU_BN', ...
        '250117_fMRI_OX_NWU_VS'}; 

session_count = [4, 13, 10, 15, 16, 10]; % number of sessions for each subj so far. EDIT as needed. 
TR= 0.76; 

subjidx = 5; % enter subjidx here 
subjname = ['subj_', num2str(subjidx)];  
subjname_real = SUBJNAMES{subjidx}; 

% mridir = '/Volumes/ExtremeSSD/OX_DATA/MRI'; 
mridir = '/Users/qhyang/Desktop/OX_DATA/MRI'; 

mridatapath = fullfile(mridir, ['subj_', num2str(subjidx)], 'nifti'); % sn:subject's name

% evdir = '/Volumes/ExtremeSSD/OX_DATA/labchart'; 
evdir = '/Users/qhyang/Desktop/OX_DATA/labchart'; 

nruns = 80; 
stimdur = 2;

outdir = fullfile(mridatapath, 'single_trial_by_category');

load(fullfile(outdir, 'TYPED_FITHRF_GLMDENOISE_RR.mat')); 
[odor, category] = OX_get_odor(subjname); 

personi = find(strcmp(category, 'PERSON'));
foodi = find(strcmp(category, 'FOOD'));
loci = find(strcmp(category, 'LOCATION')); 
controli = find(strcmp(category, 'CONTROL')); 

categoryidx = zeros(size(category)); 
categoryidx(personi) = 1; 
categoryidx(foodi) = 2; 
categoryidx(loci) = 3; 
categoryidx(controli) = 4; 


modelmd = squeeze(modelmd);

gmmask = spm_vol(fullfile(mridatapath, 'coreg', 'gm_mask_thr05_func.nii')); 
gm_mask = spm_read_vols(gmmask); 

rmask = R2>0.5; 
%% building serachlight neighborhoods 

gm_inds = find(gm_mask);                      % 50277 x 1
idx_vol = nan(size(gm_mask));                 % 104 x 96 x 30
idx_vol(gm_inds) = 1:length(gm_inds);         % maps to modelmd rows

modelmd_inds = find(rmask > 0);               % indices into modelmd (not into 3D)
gm_voxel_inds = gm_inds(modelmd_inds);        % linear indices into 3D brain
[xc, yc, zc] = ind2sub(size(gm_mask), gm_voxel_inds);

searchlight_radius = 2;
neighborhoods = cell(length(modelmd_inds), 1);

for vi = 1:length(modelmd_inds)
    cx = xc(vi); cy = yc(vi); cz = zc(vi);

    [xg, yg, zg] = ndgrid(cx-searchlight_radius:cx+searchlight_radius, ...
                          cy-searchlight_radius:cy+searchlight_radius, ...
                          cz-searchlight_radius:cz+searchlight_radius);

    % Keep within bounds
    valid = xg > 0 & xg <= size(gm_mask,1) & ...
            yg > 0 & yg <= size(gm_mask,2) & ...
            zg > 0 & zg <= size(gm_mask,3);

    xg = xg(valid); yg = yg(valid); zg = zg(valid);

    % Optional: spherical
    dist = sqrt((xg - cx).^2 + (yg - cy).^2 + (zg - cz).^2);
    keep = dist <= searchlight_radius;

    xg = xg(keep); yg = yg(keep); zg = zg(keep);

    % Convert to modelmd indices using idx_vol
    lin_idx = sub2ind(size(gm_mask), xg, yg, zg);
    model_idx = idx_vol(lin_idx);
    neighborhoods{vi} = model_idx(~isnan(model_idx));
end

%% category 
% Parallel pool
if isempty(gcp('nocreate'))
    parpool; % opens default pool
end

nVoxels_total = length(neighborhoods); 
nClasses = numel(unique(categoryidx));

% Allocate results
voxel_acc = nan(nVoxels_total, 1);
acc_per_class = nan(nVoxels_total, nClasses);

fprintf('Running searchlight decoding...\n');
parfor vi = 1:nVoxels_total
    voxel_inds = neighborhoods{vi};
    if numel(voxel_inds) < 10
        continue;
    end

    X = modelmd(voxel_inds, :)'; % trials x features
    y = categoryidx';               % trials x 1
%     y = odor; 
    X = zscore(X);

    preds = zeros(size(y));
    for t = 1:length(y)
        Xtrain = X; Xtrain(t,:) = [];
        ytrain = y; ytrain(t) = [];

        class_labels = unique(ytrain);
        template = arrayfun(@(c) mean(Xtrain(ytrain==c,:),1), class_labels, 'UniformOutput', false);
        template = cat(1, template{:});
        test_corr = corr(template', X(t,:)');
        [~, max_idx] = max(test_corr);
        preds(t) = class_labels(max_idx);
    end

    voxel_acc(vi) = mean(preds == y);

    true_labels = y;

    for c = 1:nClasses
        class_id = c;
        is_this_class = true_labels == class_id;
        acc_per_class(vi, c) = mean(preds(is_this_class) == true_labels(is_this_class));
    end

end

acc_per_class(acc_per_class == 1) = 0; 

gm_inds = find(gm_mask);                     % all valid gray matter voxels
final_mask_inds = gm_inds(rmask > 0);        % linear indices into gm_mask for selected centers
searchlight_accuracy = zeros(size(gm_mask));   % 104 x 96 x 30
searchlight_accuracy_person = zeros(size(gm_mask));   % 104 x 96 x 30
searchlight_accuracy_food = zeros(size(gm_mask));   % 104 x 96 x 30
searchlight_accuracy_loc = zeros(size(gm_mask));   % 104 x 96 x 30

for vi = 1:length(final_mask_inds)
    searchlight_accuracy(final_mask_inds(vi)) = voxel_acc(vi);
    searchlight_accuracy_person(final_mask_inds(vi)) = acc_per_class(vi, 1);
    searchlight_accuracy_food(final_mask_inds(vi)) = acc_per_class(vi, 2);
    searchlight_accuracy_loc(final_mask_inds(vi)) = acc_per_class(vi, 3);

end

fprintf('Done! \n');




%% Save

% Copy header and modify it for new data
out_nii = gmmask;
out_nii.fname = fullfile(outdir, 'category_searchlight_accuracy.nii');   % output filename
out_nii.dt = [16, 0];                         % set data type to float32

% Write volume using SPM
spm_write_vol(out_nii, single(searchlight_accuracy));

%% save class accuracy 
out_nii.fname = fullfile(outdir, 'person_searchlight_accuracy.nii');   % output filename
% Write volume using SPM
spm_write_vol(out_nii, single(searchlight_accuracy_person));

out_nii.fname = fullfile(outdir, 'food_searchlight_accuracy.nii');   % output filename
% Write volume using SPM
spm_write_vol(out_nii, single(searchlight_accuracy_food));

out_nii.fname = fullfile(outdir, 'loc_searchlight_accuracy.nii');   % output filename
% Write volume using SPM
spm_write_vol(out_nii, single(searchlight_accuracy_loc));



%% 
save('category_searchlight_accuracy.mat', 'voxel_acc', 'acc_per_class'); 
%% permuted null dist for all category 
nPerms = 100;
nSample = 1000;
nVoxelsTotal = length(neighborhoods);
nTrials = size(modelmd, 2);

null_accuracy = nan(nSample, nPerms);      % results per sampled voxel
voxel_indices = nan(nSample, nPerms);      % which voxels were sampled

for p = 1:nPerms
    fprintf('Permutation %d/%d\n', p, nPerms);

    y_perm = categoryidx(randperm(nTrials));
%     y_perm = odor(randperm(nTrials));
    acc_this_perm = nan(nSample, 1);

    sampled_voxels = randperm(nVoxelsTotal, nSample);
    % Extract sampled neighborhoods directly using parfor index
    this_neighborhoods = neighborhoods(sampled_voxels);  % 1 x nSample

    parfor i = 1:nSample
        voxel_inds = this_neighborhoods{i};
        if numel(voxel_inds) < 5
            continue;
        end

        X = modelmd(voxel_inds, :)';
        X = zscore(X);
        preds = zeros(nTrials, 1);

        for t = 1:nTrials
            Xtrain = X; Xtrain(t,:) = [];
            ytrain = y_perm; ytrain(t) = [];

            class_labels = unique(ytrain);
            template = arrayfun(@(c) mean(Xtrain(ytrain == c,:), 1), ...
                                class_labels, 'UniformOutput', false);
            template = cat(1, template{:});
            test_corr = corr(template', X(t,:)');
            [~, max_idx] = max(test_corr);
            preds(t) = class_labels(max_idx);
        end

        % Only update acc_this_perm(i) — always defined
        acc_this_perm(i) = mean(preds' == y_perm);
    end

    null_accuracy(:, p) = acc_this_perm;
end


%% FDR correction for single null dist

nVoxels = size(voxel_acc, 1);
nPerms = size(null_accuracy, 2);
nSample = size(null_accuracy, 1);

p_uncorrected = nan(nVoxels, 1);
p_fdr_corrected = nan(nVoxels, 1);

null_accuracy(null_accuracy == 1) = nan; 


% Flatten null distribution across permutations
null_dist = null_accuracy(:); 

% Compute uncorrected p-values for each voxel
pvals = nan(nVoxels, 1);
for v = 1:nVoxels
    real_acc = voxel_acc(v);
    pvals(v) = mean(null_dist >= real_acc);
end

% FDR correction
[~, pvals_fdr] = fdr(pvals);

% Store
p_uncorrected(:, 1) = pvals;
p_fdr_corrected(:, 1) = pvals_fdr;

% Get all voxels with FDR-corrected p < 0.05
sig_voxels = p_fdr_corrected(:) < 0.05;

if any(sig_voxels)
    % Accuracy threshold is the lowest accuracy among those significant
    acc_threshold_fdr05 = min(voxel_acc(sig_voxels));
else
    acc_threshold_fdr05 = NaN;  % No significant voxels
end

%% permuted null for each category 
nSample = 1000;
nPerms = 100;

nClasses = numel(unique(categoryidx));
nTrials = size(modelmd, 2);
nVoxelsTotal = length(neighborhoods);

null_accuracy = nan(nSample, nClasses, nPerms);
voxel_indices = nan(nSample, nPerms);

for p = 1:nPerms
    fprintf('Permutation %d / %d\n', p, nPerms);
    
    % Shuffle labels
    y_perm = categoryidx(randperm(nTrials));
    
    % Sample a subset of voxels
    sampled_voxels = randperm(nVoxelsTotal, nSample);
    voxel_indices(:, p) = sampled_voxels;
    this_neighborhoods = neighborhoods(sampled_voxels);

    % Preallocate this permutation's accuracy
    acc_this_perm = nan(nSample, nClasses);

    parfor i = 1:nSample
        voxel_inds = this_neighborhoods{i};
        if numel(voxel_inds) < 10
            continue;
        end

        % Extract data: [trials x features]
        X = modelmd(voxel_inds, :)';
        X = zscore(X);
        preds = zeros(nTrials, 1);

        for t = 1:nTrials
            Xtrain = X; Xtrain(t,:) = [];
            ytrain = y_perm; ytrain(t) = [];

            class_labels = unique(ytrain);
            template = arrayfun(@(c) mean(Xtrain(ytrain == c, :), 1), ...
                                class_labels, 'UniformOutput', false);
            template = cat(1, template{:});
            test_corr = corr(template', X(t,:)');
            [~, max_idx] = max(test_corr);
            preds(t) = class_labels(max_idx);
        end

        % Per-class decoding accuracy
        for c = 1:nClasses
            true_class_trials = (y_perm == c);
            acc_this_perm(i, c) = mean(preds(true_class_trials)' == y_perm(true_class_trials));
        end
    end

    null_accuracy(:, :, p) = acc_this_perm;
end

%% FDR correction 

nVoxels = size(acc_per_class, 1);
nClasses = size(acc_per_class, 2);
nPerms = size(null_accuracy, 3);
nSample = size(null_accuracy, 1);

p_uncorrected = nan(nVoxels, nClasses);
p_fdr_corrected = nan(nVoxels, nClasses);

null_accuracy(null_accuracy == 1) = nan; 

for c = 1:nClasses
    % Flatten null distribution across permutations
    null_dist = reshape(null_accuracy(:, c, :), nSample * nPerms, 1);

    % Compute uncorrected p-values for each voxel
    pvals = nan(nVoxels, 1);
    for v = 1:nVoxels
        real_acc = acc_per_class(v, c);
        pvals(v) = mean(null_dist >= real_acc);
    end

    % FDR correction
    [~, pvals_fdr] = fdr(pvals);

    % Store
    p_uncorrected(:, c) = pvals;
    p_fdr_corrected(:, c) = pvals_fdr;
end

% Initialize threshold matrix
acc_threshold_fdr05 = nan(1, nClasses);

for c = 1:nClasses
    % Get all voxels with FDR-corrected p < 0.05
    sig_voxels = p_fdr_corrected(:, c) < 0.05;

    if any(sig_voxels)
        % Accuracy threshold is the lowest accuracy among those significant
        acc_threshold_fdr05(c) = min(acc_per_class(sig_voxels, c));
    else
        acc_threshold_fdr05(c) = NaN;  % No significant voxels
    end
end



%% odor 

% Parallel pool
if isempty(gcp('nocreate'))
    parpool; % opens default pool
end

nVoxels_total = length(neighborhoods); 
nClasses = numel(unique(categoryidx));

% Allocate results
voxel_acc = nan(nVoxels_total, 1);
acc_per_class = nan(nVoxels_total, nClasses);

fprintf('Running searchlight decoding...\n');
parfor vi = 1:nVoxels_total
    voxel_inds = neighborhoods{vi};
    if numel(voxel_inds) < 10
        continue;
    end

    X = modelmd(voxel_inds, :)'; % trials x features
    y = odor; 
    X = zscore(X);

    preds = zeros(size(y));
    for t = 1:length(y)
        Xtrain = X; Xtrain(t,:) = []; % leave one out 
        ytrain = y; ytrain(t) = [];

        class_labels = unique(ytrain);
        template = arrayfun(@(c) mean(Xtrain(ytrain==c,:),1), class_labels, 'UniformOutput', false);
        template = cat(1, template{:});
        test_corr = corr(template', X(t,:)'); % by reversing this, we'll find odors that are farthest way from each other in representation distance 
        [~, max_idx] = max(test_corr);
        preds(t) = class_labels(max_idx);
    end

    voxel_acc(vi) = mean(preds == y);

    true_labels = y;

    for c = 1:nClasses
        class_id = c;
        is_this_class = true_labels == class_id;
        acc_per_class(vi, c) = mean(preds(is_this_class) == true_labels(is_this_class));
    end

end

acc_per_class(acc_per_class == 1) = 0; 

gm_inds = find(gm_mask);                     % all valid gray matter voxels
final_mask_inds = gm_inds(rmask > 0);        % linear indices into gm_mask for selected centers
searchlight_accuracy = zeros(size(gm_mask));   % 104 x 96 x 30
searchlight_accuracy_person = zeros(size(gm_mask));   % 104 x 96 x 30
searchlight_accuracy_food = zeros(size(gm_mask));   % 104 x 96 x 30
searchlight_accuracy_loc = zeros(size(gm_mask));   % 104 x 96 x 30

for vi = 1:length(final_mask_inds)
    searchlight_accuracy(final_mask_inds(vi)) = voxel_acc(vi);
    searchlight_accuracy_person(final_mask_inds(vi)) = acc_per_class(vi, 1);
    searchlight_accuracy_food(final_mask_inds(vi)) = acc_per_class(vi, 2);
    searchlight_accuracy_loc(final_mask_inds(vi)) = acc_per_class(vi, 3);

end

fprintf('Done! \n');



%% save map 

out_nii = gmmask;
out_nii.fname = fullfile(outdir, 'odor_searchlight_accuracy.nii');   % output filename
out_nii.dt = [16, 0];                         % set data type to float32

% Write volume using SPM
spm_write_vol(out_nii, single(searchlight_accuracy));

out_nii.fname = fullfile(outdir, 'person_odor_searchlight_accuracy.nii');   % output filename
% Write volume using SPM
spm_write_vol(out_nii, single(searchlight_accuracy_person));

out_nii.fname = fullfile(outdir, 'food_odor_searchlight_accuracy.nii');   % output filename
% Write volume using SPM
spm_write_vol(out_nii, single(searchlight_accuracy_food));

out_nii.fname = fullfile(outdir, 'loc_odor_searchlight_accuracy.nii');   % output filename
% Write volume using SPM
spm_write_vol(out_nii, single(searchlight_accuracy_loc));


%% save accuracy file 
save(fullfile(outdir, 'odor_searchlight_accuracy.mat'), 'voxel_acc', 'acc_per_class'); 


%% single null dist for all odors 
nPerms = 100;
nSample = 1000;
nVoxelsTotal = length(neighborhoods);
nTrials = size(modelmd, 2);

null_accuracy = nan(nSample, nPerms);      % results per sampled voxel
voxel_indices = nan(nSample, nPerms);      % which voxels were sampled

for p = 1:nPerms
    fprintf('Permutation %d/%d\n', p, nPerms);

%     y_perm = categoryidx(randperm(nTrials));
    y_perm = odor(randperm(nTrials));
    acc_this_perm = nan(nSample, 1);

    sampled_voxels = randperm(nVoxelsTotal, nSample);
    % Extract sampled neighborhoods directly using parfor index
    this_neighborhoods = neighborhoods(sampled_voxels);  % 1 x nSample

    parfor i = 1:nSample
        voxel_inds = this_neighborhoods{i};
        if numel(voxel_inds) < 5
            continue;
        end

        X = modelmd(voxel_inds, :)';
        X = zscore(X);
        preds = zeros(nTrials, 1);

        for t = 1:nTrials
            Xtrain = X; Xtrain(t,:) = [];
            ytrain = y_perm; ytrain(t) = [];

            class_labels = unique(ytrain);
            template = arrayfun(@(c) mean(Xtrain(ytrain == c,:), 1), ...
                                class_labels, 'UniformOutput', false);
            template = cat(1, template{:});
            test_corr = corr(template', X(t,:)');
            [~, max_idx] = max(test_corr);
            preds(t) = class_labels(max_idx);
        end

        % Only update acc_this_perm(i) — always defined
        acc_this_perm(i) = mean(preds == y_perm);
    end

    null_accuracy(:, p) = acc_this_perm;
end

%% FDR correction for odor null dist 

nVoxels = size(voxel_acc, 1);
nPerms = size(null_accuracy, 2);
nSample = size(null_accuracy, 1);

p_uncorrected = nan(nVoxels, 1);
p_fdr_corrected = nan(nVoxels, 1);

null_accuracy(null_accuracy == 1) = nan; 


% Flatten null distribution across permutations
null_dist = null_accuracy(:); 

% Compute uncorrected p-values for each voxel
pvals = nan(nVoxels, 1);
for v = 1:nVoxels
    real_acc = voxel_acc(v);
    pvals(v) = mean(null_dist >= real_acc);
end

% FDR correction
[~, pvals_fdr] = fdr(pvals);

% Store
p_uncorrected(:, 1) = pvals;
p_fdr_corrected(:, 1) = pvals_fdr;

% Get all voxels with FDR-corrected p < 0.05
sig_voxels = p_fdr_corrected(:) < 0.05;

if any(sig_voxels)
    % Accuracy threshold is the lowest accuracy among those significant
    acc_threshold_fdr05 = min(voxel_acc(sig_voxels));
else
    acc_threshold_fdr05 = NaN;  % No significant voxels
end


%% permuted null for each odors in each category 
nSample = 1000;
nPerms = 100;

nClasses = numel(unique(categoryidx));
nTrials = size(modelmd, 2);
nVoxelsTotal = length(neighborhoods);

null_accuracy = nan(nSample, nClasses, nPerms);
voxel_indices = nan(nSample, nPerms);

for p = 1:nPerms
    fprintf('Permutation %d / %d\n', p, nPerms);
    
    % Shuffle labels
    y_perm = odor(randperm(nTrials));
    label_perm = categoryidx(randperm(nTrials));
    
    % Sample a subset of voxels
    sampled_voxels = randperm(nVoxelsTotal, nSample);
    voxel_indices(:, p) = sampled_voxels;
    this_neighborhoods = neighborhoods(sampled_voxels);

    % Preallocate this permutation's accuracy
    acc_this_perm = nan(nSample, nClasses);

    parfor i = 1:nSample
        voxel_inds = this_neighborhoods{i};
        if numel(voxel_inds) < 10
            continue;
        end

        % Extract data: [trials x features]
        X = modelmd(voxel_inds, :)';
        X = zscore(X);
        preds = zeros(nTrials, 1);

        for t = 1:nTrials
            Xtrain = X; Xtrain(t,:) = [];
            ytrain = y_perm; ytrain(t) = [];

            class_labels = unique(ytrain);
            template = arrayfun(@(c) mean(Xtrain(ytrain == c, :), 1), ...
                                class_labels, 'UniformOutput', false);
            template = cat(1, template{:});
            test_corr = corr(template', X(t,:)');
            [~, max_idx] = max(test_corr);
            preds(t) = class_labels(max_idx);
        end

        % Per-class decoding accuracy
        for c = 1:nClasses
            true_class_trials = (label_perm == c);
            acc_this_perm(i, c) = mean(preds(true_class_trials) == y_perm(true_class_trials));
        end
    end

    null_accuracy(:, :, p) = acc_this_perm;
end

%% FDR correction 

nVoxels = size(acc_per_class, 1);
nClasses = size(acc_per_class, 2);
nPerms = size(null_accuracy, 3);
nSample = size(null_accuracy, 1);

p_uncorrected = nan(nVoxels, nClasses);
p_fdr_corrected = nan(nVoxels, nClasses);

null_accuracy(null_accuracy == 1) = nan; 

for c = 1:nClasses
    % Flatten null distribution across permutations
    null_dist = reshape(null_accuracy(:, c, :), nSample * nPerms, 1);

    % Compute uncorrected p-values for each voxel
    pvals = nan(nVoxels, 1);
    for v = 1:nVoxels
        real_acc = acc_per_class(v, c);
        pvals(v) = mean(null_dist >= real_acc);
    end

    % FDR correction
    [~, pvals_fdr] = fdr(pvals);

    % Store
    p_uncorrected(:, c) = pvals;
    p_fdr_corrected(:, c) = pvals_fdr;
end

% Initialize threshold matrix
acc_threshold_fdr05 = nan(1, nClasses);

for c = 1:nClasses
    % Get all voxels with FDR-corrected p < 0.05
    sig_voxels = p_fdr_corrected(:, c) < 0.05;

    if any(sig_voxels)
        % Accuracy threshold is the lowest accuracy among those significant
        acc_threshold_fdr05(c) = min(acc_per_class(sig_voxels, c));
    else
        acc_threshold_fdr05(c) = NaN;  % No significant voxels
    end
end


%% sanity check for shape of searchlight
% % Assumes:
% % - neighborhoods: cell array of voxel indices into modelmd (precomputed)
% % - final_mask_inds: linear indices into gm_mask
% % - gm_mask: 3D volume, size [104 x 96 x 30]
% 
% % Get modelmd indices in the first neighborhood
% first_nb = neighborhoods{15000};  % indices into modelmd
% fprintf('Neighborhood size: %d voxels\n', numel(first_nb));
% 
% % Convert modelmd indices to 3D coordinates via gm_mask
% gm_inds = find(gmvol);  % 50277 x 1
% voxel_lin_inds = gm_inds(first_nb);  % linear indices into 3D volume
% [xv, yv, zv] = ind2sub(size(gmvol), voxel_lin_inds);
% 
% % Create an empty binary volume and mark neighborhood voxels
% nb_img = false(size(gmvol));
% nb_img(voxel_lin_inds) = true;
% 
% % Plot as orthogonal slices
% figure;
% ortho_slice = @(v) imagesc(squeeze(v(:,:,round(size(v,3)/2)))');  % simple orthoslice
% subplot(1,3,1); imagesc(squeeze(nb_img(:, :, round(mean(zv))))'); axis equal tight; title('Axial');
% subplot(1,3,2); imagesc(squeeze(nb_img(:, round(mean(yv)), :))'); axis equal tight; title('Coronal');
% subplot(1,3,3); imagesc(squeeze(nb_img(round(mean(xv)), :, :))'); axis equal tight; title('Sagittal');

