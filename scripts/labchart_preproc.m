%% lab chart data preprocessing 



SUBJNAMES = {'240711_fMRI_OX_NWU_AS', ...
        '240723_fMRI_OX_NWU_LS', ...
        '240814_fMRI_OX_NWU_JN', ...
        '240816_fMRI_OX_NWU_RR', ...
        '241018_fMRI_OX_NWU_BN', ...
        '250117_fMRI_OX_NWU_VS'}; 

session_count = [4, 13, 10, 15, 16, 10]; % number of sessions for each subj so far. EDIT as needed. 

%% enter subjID and session ID here 

% wkdir = '/Volumes/ExtremeSSD/OX_DATA/labchart'; 
wkdir = '/Users/qhyang/Desktop/OX_DATA/labchart'; 

subjidx = 6; % enter subjidx here 

subjname = ['subj_', num2str(subjidx)];  
subjname_real = SUBJNAMES{subjidx}; 

%%% we can reject all labchart data that are less than 300000 points and more than 420000. This
%%% is between 5 and 7 min. 
%%% then check channel 2 and reject all runs with incorrect number of
%%% pulses. 

% mridir = '/Volumes/ExtremeSSD/OX_DATA/MRI'; 
mridir = '/Users/qhyang/Desktop/OX_DATA/MRI'; 
mridatapath = fullfile(mridir, ['subj_', num2str(subjidx)], 'nifti'); % sn:subject's name

%% 

eventdata = {}; 
nrun = []; 
for sessionID = 1: session_count(subjidx) 
    
    try 
        load(fullfile(wkdir, 'raw_data', [subjname, '_ses', num2str(sessionID)])); 
    catch
        continue 
    end 
    
    var = who;

    if sum(ismember(var, 'data')) == 1
      tvec = dataend(1,:)-datastart(1,:); 
      tidx = ((tvec > 350000) + (tvec < 420000))>1; 

      
      %%% exceptions 
      if subjidx == 3 && sessionID == 1
          tidx(4) =0; 
      end 
      if subjidx == 3 && sessionID == 7
          tidx(3) = 1; 
      end 
      
      if subjidx == 4 && sessionID == 3
          tidx(2) =0; 
      end 
    
      if subjidx == 6 && sessionID == 6
         tidx(7) = 0; 
      end 

      %%% exception: subj4ses12 run1(block2) has missing scanner TTL at the
      %%% beginning due to BNC cable being plugged in after scanner start
      %%% time (visually confirmed with a baseline shift). Pad the scanner
      %%% TTL so it matches nframes at the beginning.
      if subjidx == 4 && sessionID == 12
          blockidx = 2;
          mri_data = data(datastart(4, blockidx):dataend(4, blockidx));
          mri_binary = mri_data > 1;
          mri_rises = find(diff([false, mri_binary]) == 1);
          mri_falls = find(diff([mri_binary, false]) == -1);

          mri_file = dir(fullfile(mridatapath, 'func', ...
              [subjname_real, '_', num2str(sessionID), '_Run1_*.nii']));
          if numel(mri_file) ~= 1
              error('Expected one subject 4 session 12 run 1 NIfTI, found %d.', ...
                  numel(mri_file));
          end
          mri_info = niftiinfo(fullfile(mri_file.folder, mri_file.name));
          expected_nframes = mri_info.ImageSize(4);
          n_missing_ttls = expected_nframes - numel(mri_rises);

          if n_missing_ttls < 0
              error('LabChart block 2 has more scanner TTLs than fMRI frames.');
          elseif n_missing_ttls > 0
              tr_samples = round(median(diff(mri_rises)));
              pulse_width = round(median(mri_falls - mri_rises + 1));
              pulse_template = mri_data(mri_rises(1):mri_rises(1) + pulse_width - 1);
              padded_rises = mri_rises(1) - (n_missing_ttls:-1:1) * tr_samples;

              if padded_rises(1) < 1
                  error('Not enough leading samples to pad the missing scanner TTLs.');
              end
              for onset = padded_rises
                  mri_data(onset:onset + pulse_width - 1) = pulse_template;
              end
              data(datastart(4, blockidx):dataend(4, blockidx)) = mri_data;
          end

          padded_count = sum(diff([false, mri_data > 1]) == 1);
          if padded_count ~= expected_nframes
              error('Padded scanner TTL count (%d) does not match nframes (%d).', ...
                  padded_count, expected_nframes);
          end
      end


      data_chunks = cell(sum(tidx), 1); 
      
      datastart = datastart(:, tidx); 
      dataend = dataend(:, tidx); 
      for i = 1: size(datastart, 2)
          data_chunks{i, 1} = [data(datastart(1, i):dataend(1, i)); ...
              data(datastart(2, i):dataend(2, i)); ...
              data(datastart(4, i):dataend(4, i))]; 
      end

    else 
        var_resp = who('C1B*'); 
        var_daq = who('C2B*'); 
        var_mri = who('C4B*'); 
        [~, idx] = sort(cellfun(@(x) sscanf(x, 'C1B%d'), var_resp));
        var_resp = var_resp(idx);
        [~, idx] = sort(cellfun(@(x) sscanf(x, 'C2B%d'), var_daq));
        var_daq = var_daq(idx);
        [~, idx] = sort(cellfun(@(x) sscanf(x, 'C4B%d'), var_mri));
        var_mri = var_mri(idx);
        tvec = []; 
        for n = 1: length(var_resp)
            tvec(n ) = length(eval(var_resp{n}));
        end 
        tidx = ((tvec > 300000) + (tvec < 420000))>1;
        data_chunks = cell(sum(tidx), 1);
        
        var_resp = var_resp(tidx); 
        var_daq = var_daq(tidx);
        var_mri = var_mri(tidx); 
        %%% there might be an exception in subj4 session2 

        for i = 1: length(var_resp)
            data_chunks{i, 1} = [eval(var_resp{i}), ...
                eval(var_daq{i}), ...
                eval(var_mri{i})]';
        end 
      

    end  

%%% concatenate all of data chunks across runs and save a run counter 
    eventdata = [eventdata; data_chunks]; 
    nrun = [nrun, sum(tidx)]; 
    
    clearvars -except eventdata nrun subjname wkdir subjname_real session_count mridatapath subjidx
end 

%% double check. all channel 2 should have 23 pulses 

for runidx = 1: length(eventdata)
    current_events = eventdata{runidx}(2, :);
    figure; 
    plot(current_events); 
end


%% manully check nrun and make sure they match with n of scans!!! 
% 
% exceptions: 
%%% subj3 session 1 has repeated run1s. use second one (lab chart data
%%% start from segment 6) 

%% save raw labchart data  
% save("labchart/subj4_events_raw.mat", "eventdata"); 


%% crude sniff onset with SNIFF pulse 

nruns = length(eventdata); 
mripulse = cell(nruns,1);
events = cell(nruns,1);
resp = cell(nruns,1);

for i = 1: nruns

    mripulse{i} = eventdata{i}(3, :);
    events{i} = eventdata{i}(2, :);
    resp{i} = eventdata{i}(1, :);
end 

event_onsets = zeros(10,nruns); 
cue_onsets = zeros(10, nruns); 

for evidx = 1:nruns
    currentevents = events{evidx};
    eventbi = currentevents>0.1;
    event_real = eventbi;
    for i = 1:length(eventbi)-1
        if eventbi(i+1) == eventbi(i)
            event_real(i+1) = 0;
        end
    end
    evonsets = find(event_real>0);
    event_o = evonsets(5:2:end);
    event_q = evonsets(4:2:end-1);

    mrionset = find(mripulse{evidx} > 1, 1); 
    event_o = event_o - mrionset; % set MRI onset to time 0 
    event_o = event_o/1000;
    event_onsets(:, evidx) = event_o; 
    
    event_q = event_q - mrionset; % set MRI onset to time 0 
    event_q = event_q/1000;
    cue_onsets(:, evidx) = event_q; 

end


%% check if all cue onsets are >0
goodtrials = ones(numel(cue_onsets), 1); 
goodtrials(cue_onsets<0) = 0; 
goodtrials = logical(goodtrials); 
%% get nframes for each session from MRI data 
nframes = getnframes(subjname_real, mridatapath, 1, session_count(subjidx)); 
%% !!!!!! NOTE for subject 5 BN. last two sessions may have some runs where scanner is started too late and cue onset is outside of scanner TTL. 
%% 
offsets = cumsum(nframes); 
offsets = offsets';
offsets = [0, offsets(1:nruns-1)];

%%% convert offsets from TRs to seconds 
TR = 0.76; % TR in seconds 0
offsets = offsets*TR; 

event_onsets_all = offsets+event_onsets;
event_onsets_vec = reshape(event_onsets_all, [numel(event_onsets_all), 1]); 

cue_onsets_all = offsets+cue_onsets;
cue_onsets_vec = reshape(cue_onsets_all, [numel(cue_onsets_all), 1]); 


%% save event onsets 
sessi = 1; 
sessf = session_count(subjidx); 
save('labchart/extracted_events/subj_6_events.mat', 'event_onsets', 'cue_onsets', "event_onsets_vec", "event_onsets_all", "cue_onsets_vec", "cue_onsets_all", "sessi", "sessf", 'nframes', "offsets", "goodtrials"); 

%% functions 
function nframes = getnframes(subjname, datapath, nsess_i, nsess_f)

rcntr = 0; % run counter
nframes = [];

for sess = nsess_i:nsess_f
    % sess_2_run_1 is counted as 5 if sess_1 had 4 runs
    path_ = fullfile(datapath, 'func/');
    n = dir(fullfile(path_, [subjname, '_', num2str(sess), '_*.nii']));
    nfiles = length(n);

    for i=1:nfiles

        rcntr = rcntr+1;
        thisfile = dir(fullfile(path_, [subjname, '_', num2str(sess), '_Run', num2str(i), '_*.nii']));
        fname = fullfile(path_, thisfile.name);

        nifti_header = niftiinfo(fname);
        nvols = nifti_header.ImageSize(4);

        nframes(rcntr, 1) = nvols;
    end

end
end
